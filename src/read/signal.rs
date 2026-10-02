//! Per-base signal: bedGraph, `samtools depth`, and a bare column of values.
//!
//! Three formats that give a value to every base, read as half-open
//! `(start, end, value)` spans for a coverage track and as windows for a
//! window track, which takes
//! bedGraph alone. bedGraph is 0-based and half-open and passes straight
//! through, every base of an interval taking the interval's value. `samtools
//! depth` is 1-based, one position to a line, so one comes off on the way in.
//! A bare column names no coordinate at all and starts at the left edge of the
//! region, so the same file over a different window is a different figure.
//!
//! # Which of the three a file is, and where the guess breaks
//!
//! Nothing in these files says which format they are, so the column count is
//! the evidence: four is bedGraph, three is depth, one is a bare value. It is
//! read off the first line carrying data and the rest of the file has to keep
//! to it, since a file that changes width halfway is one whose positions
//! cannot be trusted.
//!
//! `samtools depth` over several files writes a depth column per file, so two
//! of them come to four columns and read as a bedGraph, each record becoming a
//! run of bases at the height of the second sample: a figure that is wrong and
//! does not look wrong. Overlap tells them apart, since a bedGraph never puts
//! two intervals over the same base and depth records read as intervals cover
//! each other nearly everywhere. An overlap is therefore refused, with a
//! message naming `--format depth` for the file it really was and
//! `--format bedgraph` for a bedGraph that was merely out of order. Three
//! columns has no such tell: a BED3 carries no value column whose absence
//! could be noticed, so `chr1 100 200` reads as depth and becomes one position
//! at 0-based 99 carrying a depth of 200.
//!
//! # What stops the read, and what goes past without a word
//!
//! A row on another sequence, or one lying entirely outside the window, is
//! skipped: handing over a whole genome file to draw one locus is the ordinary
//! way to use this. A row that does not parse is the other case. The file is
//! then not what it was taken for, and a reader that stepped over the row would
//! draw a figure with data missing from it and say nothing, so it stops on the
//! line and names the field that would not read.
//!
//! # Every sequence at once
//!
//! [`genome_spans`] and [`genome_windows`] read the same files with no window,
//! each row kept on the sequence it names, for a figure laid across a whole
//! genome. A bare column of values has no sequence to be laid on there, and is
//! refused. The command line draws these when a coverage or a window file is
//! named with no place.

use crate::{Region, Window};

use super::Format;
use super::{columns, lines, number, ReadError};

/// Reads a value per position, as 0-based half-open `(start, end, value)`
/// spans.
///
/// Spans rather than one entry per base, because a bedGraph row is one row
/// however many bases it covers. Expanded here, a file tiling a whole
/// chromosome became one entry per base of it and a kilobyte of input asked
/// for six gigabytes; the span form costs the size of the file.
/// [`CoverageTrack::from_spans`](crate::CoverageTrack::from_spans) lays them
/// over the region once.
///
/// A position no span covers stays at zero, which is what a depth of zero
/// means and what a bedGraph leaves out. That is a decision rather than an
/// accident, and it is why this reader hands back what the file stated rather
/// than what it did not.
///
/// The format is told by the shape of a line unless `format` says otherwise:
/// four columns is bedGraph (`chrom start end value`, 0-based half-open, and
/// every base of the interval gets the value), three is `samtools depth`
/// (`chrom pos depth`, 1-based), and one is a bare column of values starting
/// at the left edge of the region.
///
/// Rows on another sequence than `region.seq()` are skipped, so a whole genome
/// file can be handed over and only the window comes back. So are rows outside
/// the region, which is what keeps a genome-wide file from being widened into
/// memory.
pub fn spans(
    text: &str,
    region: &Region,
    format: Option<Format>,
) -> Result<Vec<(u64, u64, f64)>, ReadError> {
    let mut found = Vec::new();
    fold_spans(text, region, format, |start, end, value| {
        found.push((start, end, value))
    })?;
    Ok(found)
}

/// The same reading as [`spans`], handing each span to `each` as the line it
/// came from is read instead of collecting them, and answering how many there
/// were.
///
/// `samtools depth` writes a line per base, so collecting first costs
/// twenty-four bytes per base of the window before a single one is used. Over
/// ten million bases that was 231 MB standing beside the 152 MB of file text
/// it was read from and the 76 MB of the track being built out of it, and the
/// caller here paints each span as it arrives and never holds the middle one.
pub(crate) fn fold_spans(
    text: &str,
    region: &Region,
    format: Option<Format>,
    mut each: impl FnMut(u64, u64, f64),
) -> Result<usize, ReadError> {
    fold(text, Rows::In(region), format, |_, start, end, value| {
        each(start, end, value)
    })
}

/// A stretch of bases and the value over every one of them, as
/// `(start, end, value)`, 0-based and half-open.
pub type Span = (u64, u64, f64);

/// Every value of a bedGraph or a `samtools depth` file, on every sequence it
/// names, as a signal drawn across a whole genome reads it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GenomeSignal {
    /// Each sequence the file names, in the order it first names them, with
    /// its spans, 0-based and half-open on that sequence, in the order of the
    /// file.
    pub sequences: Vec<(String, Vec<Span>)>,
}

/// Reads every value of a bedGraph or a `samtools depth` file, on every
/// sequence it names, for a signal drawn across a whole genome rather than
/// over one window of it.
///
/// The file is read as [`spans`] reads it, shape, `--format` and all, and
/// only no row is left out for being somewhere else. Two intervals of a
/// bedGraph overlap only on one sequence: the first row of the second
/// sequence starts back at nought, and is no sign of anything.
///
/// # Errors
///
/// What [`spans`] refuses, and a bare column of values, which names no
/// sequence to lay its values on.
///
/// ```
/// use karyon::read::signal::genome_spans;
///
/// let text = "chr2\t0\t100\t5\nchr1\t0\t50\t3\nchr1\t50\t80\t4\n";
/// let read = genome_spans(text, None)?;
/// assert_eq!(read.sequences[0], ("chr2".to_string(), vec![(0, 100, 5.0)]));
/// assert_eq!(read.sequences[1].1, [(0, 50, 3.0), (50, 80, 4.0)]);
/// # Ok::<(), karyon::read::ReadError>(())
/// ```
pub fn genome_spans(text: &str, format: Option<Format>) -> Result<GenomeSignal, ReadError> {
    let mut found = Sequences::default();
    fold_genome_spans(text, format, |sequence, start, end, value| {
        found.of(sequence).push((start, end, value))
    })?;
    Ok(GenomeSignal {
        sequences: found.sequences,
    })
}

/// The same reading as [`genome_spans`], handing each span to `each` with
/// the sequence its row names as the line is read, and answering how many
/// there were, so a caller laying a whole genome's depth out never holds the
/// list, as [`fold_spans`] does over a window.
pub(crate) fn fold_genome_spans(
    text: &str,
    format: Option<Format>,
    each: impl FnMut(&str, u64, u64, f64),
) -> Result<usize, ReadError> {
    fold(text, Rows::All, format, each)
}

/// Which rows of a file are read: the ones over one window of one sequence,
/// or every one, on whatever sequence it names.
#[derive(Debug, Clone, Copy)]
enum Rows<'a> {
    In(&'a Region),
    All,
}

/// Each sequence a file names, in the order it first names them, with what
/// its rows hold, found by name without a search for each row: a file's rows
/// come a sequence at a time, so the last one is asked about again far more
/// often than any other.
struct Sequences<T> {
    sequences: Vec<(String, Vec<T>)>,
    index: std::collections::HashMap<String, usize>,
    last: Option<usize>,
}

impl<T> Default for Sequences<T> {
    fn default() -> Self {
        Sequences {
            sequences: Vec::new(),
            index: std::collections::HashMap::new(),
            last: None,
        }
    }
}

impl<T> Sequences<T> {
    /// What the rows on `name` hold so far.
    fn of(&mut self, name: &str) -> &mut Vec<T> {
        let at = match self.last {
            Some(at) if self.sequences[at].0 == name => at,
            _ => match self.index.get(name) {
                Some(at) => *at,
                None => {
                    self.index.insert(name.to_string(), self.sequences.len());
                    self.sequences.push((name.to_string(), Vec::new()));
                    self.sequences.len() - 1
                }
            },
        };
        self.last = Some(at);
        &mut self.sequences[at].1
    }
}

/// How far the last interval on each sequence reached, which is what says
/// whether the next one overlaps it. The sequence the rows are on now is
/// kept apart from the others, since a file's rows come a sequence at a time
/// and a lookup by name for each row of a ten million row file is time spent
/// on nothing.
#[derive(Default)]
struct Reached<'t> {
    now: Option<(&'t str, u64)>,
    earlier: std::collections::HashMap<&'t str, u64>,
}

impl<'t> Reached<'t> {
    /// The end of the last interval on `sequence`, where there was one, with
    /// `end` kept in its place.
    fn swap(&mut self, sequence: &'t str, end: u64) -> Option<u64> {
        match &mut self.now {
            Some((name, reach)) if *name == sequence => Some(std::mem::replace(reach, end)),
            _ => {
                if let Some((name, reach)) = self.now.take() {
                    self.earlier.insert(name, reach);
                }
                let before = self.earlier.remove(sequence);
                self.now = Some((sequence, end));
                before
            }
        }
    }
}

/// [`fold_spans`] and [`fold_genome_spans`], reading the rows `rows` asks
/// for and handing each to `each` with the sequence it names.
fn fold(
    text: &str,
    rows: Rows<'_>,
    format: Option<Format>,
    mut each: impl FnMut(&str, u64, u64, f64),
) -> Result<usize, ReadError> {
    let asked = match format {
        Some(Format::BedGraph) => Some(Shape::BedGraph),
        Some(Format::Depth) => Some(Shape::Depth),
        Some(Format::Values) => Some(Shape::Values),
        // BED and GFF3 name intervals with a strand and a name, not a value
        // per base. Reading one here would either fail on a column that holds
        // a gene name or, worse, take a score for a depth, so it is refused.
        Some(Format::Bed | Format::Gff3) => {
            return Err(ReadError::whole(
                "a coverage track reads a value per base, so --format takes bedgraph, depth or values here",
            ))
        }
        None => None,
    };

    // The shape is decided once, on the first line that carries data, and the
    // rest of the file has to keep to it. `samtools depth` and bedGraph both
    // start with a sequence name, so the column count is the only thing that
    // tells them apart, and a file that changes count halfway is a file whose
    // positions cannot be trusted rather than one to guess at line by line.
    let mut shape = asked;
    let mut read = 0usize;
    // Where the next bare value lands, that shape carrying no position of its
    // own. Across a whole genome it has nowhere to land, and is refused.
    let mut next = match rows {
        Rows::In(region) => region.start(),
        Rows::All => 0,
    };
    // How far the last bedGraph interval on each sequence reached, which is
    // what says whether the intervals overlap and so whether this is a
    // bedGraph at all.
    let mut reached = Reached::default();

    for (at, line) in lines(text) {
        let fields = columns(line);
        let this = match shape {
            // `samtools depth` over more than one file writes one depth column
            // per file, so an asked-for depth takes three columns or more and
            // reads the first sample.
            Some(Shape::Depth) if asked.is_some() && fields.len() >= 3 => Shape::Depth,
            Some(known) if fields.len() != known.width() => {
                let reason = if asked.is_some() {
                    format!(
                        "--format {} reads {}, and this line has {}",
                        known.word(),
                        count(known.width()),
                        count(fields.len())
                    )
                } else {
                    format!(
                        "the file began as {} with {}, and this line has {}",
                        known.name(),
                        count(known.width()),
                        count(fields.len())
                    )
                };
                return Err(ReadError::at(at, reason));
            }
            Some(known) => known,
            None => Shape::sniff(fields.len()).ok_or_else(|| {
                ReadError::at(
                    at,
                    format!(
                        "expected 4 columns for bedGraph, 3 for samtools depth or 1 for a bare value, and this line has {}",
                        count(fields.len())
                    ),
                )
            })?,
        };
        shape = Some(this);

        match this {
            Shape::BedGraph => {
                if let Rows::In(region) = rows {
                    if fields[0] != region.seq() {
                        continue;
                    }
                }
                let start: u64 = number(fields[1], "start", at)?;
                let end: u64 = number(fields[2], "end", at)?;
                if end < start {
                    // The overlap check below is the one that names the way
                    // out, and on real `samtools depth` output over two files
                    // it is never reached: the second column is a position and
                    // the third is a depth, so the very first record has an end
                    // smaller than its start and stops here instead. A reader
                    // who got the bare sentence had nothing to act on, so the
                    // guidance is on both refusals rather than on the rarer
                    // one.
                    return Err(if asked.is_none() {
                        ReadError::at(
                            at,
                            "end is before start, so this is not a bedGraph. \
                             samtools depth over more than one file also writes four \
                             columns, and its second column is a position rather than \
                             an end: pass --format depth to read it as that, or \
                             --format bedgraph to insist",
                        )
                    } else {
                        ReadError::at(at, "end is before start")
                    });
                }
                // Two bedGraph intervals never overlap. `samtools depth` over
                // two files has four columns too, and reading one as a bedGraph
                // turns each depth record into a run of bases at the height of
                // the second sample, which is a plausible looking figure of
                // nothing. Overlap is what tells them apart, so it is refused
                // here rather than guessed at.
                // Checked a sequence at a time: the first row of the next
                // sequence starts back at nought, and overlaps nothing.
                if asked.is_none() {
                    if let Some(previous) = reached.swap(fields[0], end) {
                        if start < previous {
                            return Err(ReadError::at(
                                at,
                                "these intervals overlap, so this is not a bedGraph. \
                                 samtools depth over more than one file also writes four \
                                 columns: pass --format depth to read it as that, or \
                                 --format bedgraph to insist",
                            ));
                        }
                    }
                }
                let value: f64 = number(fields[3], "value", at)?;
                // The span, not one entry per base of it. A row is a row
                // however wide it is: expanded here, a bedGraph tiling a whole
                // chromosome became one pair per base of it, and a kilobyte of
                // input asked for six gigabytes.
                let (from, to) = match rows {
                    Rows::In(region) => (start.max(region.start()), end.min(region.end())),
                    Rows::All => (start, end),
                };
                if to > from {
                    each(fields[0], from, to, value);
                    read += 1;
                }
            }
            Shape::Depth => {
                if let Rows::In(region) = rows {
                    if fields[0] != region.seq() {
                        continue;
                    }
                }
                let pos: u64 = number(fields[1], "position", at)?;
                if pos == 0 {
                    return Err(ReadError::at(
                        at,
                        "samtools depth is 1-based, and 0 is not a position",
                    ));
                }
                let depth: f64 = number(fields[2], "depth", at)?;
                // 1-based inclusive to 0-based.
                let pos = pos - 1;
                if rows.at(pos) {
                    each(fields[0], pos, pos + 1, depth);
                    read += 1;
                }
            }
            Shape::Values => {
                let Rows::In(region) = rows else {
                    return Err(ReadError::at(
                        at,
                        "a bare column of values names no sequence, and across a whole \
                         genome each value is laid on the sequence its row names: write the \
                         place first, as chr1, or give the values their sequences as bedGraph",
                    ));
                };
                let value: f64 = number(fields[0], "value", at)?;
                let pos = next;
                next += 1;
                if region.contains(pos) {
                    each(region.seq(), pos, pos + 1, value);
                    read += 1;
                }
            }
        }
    }

    Ok(read)
}

impl Rows<'_> {
    /// Whether a position on a sequence that is read at all is kept.
    fn at(self, pos: u64) -> bool {
        match self {
            Rows::In(region) => region.contains(pos),
            Rows::All => true,
        }
    }
}

/// Reads intervals with a value each, for a window track.
///
/// bedGraph, 0-based half-open, kept as intervals rather than flattened to one
/// value per base, since a window track draws the window and not the base.
pub fn windows(text: &str, region: &Region) -> Result<Vec<Window>, ReadError> {
    let mut found = Vec::new();
    window_rows(text, Rows::In(region), |_, window| found.push(window))?;
    Ok(found)
}

/// Every window of a bedGraph, on every sequence it names, for windows drawn
/// across a whole genome: each sequence in the order the file first names it,
/// with its windows, 0-based and half-open on that sequence.
///
/// # Errors
///
/// What [`windows`] refuses.
///
/// ```
/// use karyon::read::signal::genome_windows;
///
/// let text = "1\t0\t500000\t0.25\n2\t0\t500000\t-0.5\n";
/// let read = genome_windows(text)?;
/// assert_eq!(read.len(), 2);
/// assert_eq!(read[1].0, "2");
/// assert_eq!(read[1].1[0].value, -0.5);
/// # Ok::<(), karyon::read::ReadError>(())
/// ```
pub fn genome_windows(text: &str) -> Result<Vec<(String, Vec<Window>)>, ReadError> {
    let mut found = Sequences::default();
    window_rows(text, Rows::All, |sequence, window| {
        found.of(sequence).push(window)
    })?;
    Ok(found.sequences)
}

/// [`windows`] and [`genome_windows`], handing each window `rows` asks for to
/// `each` with the sequence its row names.
fn window_rows(
    text: &str,
    rows: Rows<'_>,
    mut each: impl FnMut(&str, Window),
) -> Result<(), ReadError> {
    for (at, line) in lines(text) {
        let fields = columns(line);
        if fields.len() < 4 {
            return Err(ReadError::at(
                at,
                format!(
                    "a window file is bedGraph, chrom start end value, and this line has {}",
                    count(fields.len())
                ),
            ));
        }
        if let Rows::In(region) = rows {
            if fields[0] != region.seq() {
                continue;
            }
        }
        let start: u64 = number(fields[1], "start", at)?;
        let end: u64 = number(fields[2], "end", at)?;
        if end < start {
            return Err(ReadError::at(at, "end is before start"));
        }
        let value: f64 = number(fields[3], "value", at)?;
        // A window that does not reach the region is not drawn and does not
        // count towards the value axis either, so it is left behind here
        // rather than handed over. The bounds of the ones that do reach it are
        // passed through whole: the track clips them to the region itself, and
        // a window cut short would draw as a window of another size.
        if let Rows::In(region) = rows {
            if end <= region.start() || start >= region.end() {
                continue;
            }
        }
        each(fields[0], Window::new(start, end, value));
    }
    Ok(())
}

/// Which of the three per-base shapes a file is written in.
#[derive(Debug, Clone, Copy)]
enum Shape {
    /// `chrom start end value`, 0-based half-open, every base of the interval
    /// taking the value.
    BedGraph,
    /// `chrom pos depth`, 1-based, one position per line.
    Depth,
    /// One value per line, starting at the left edge of the region.
    Values,
}

impl Shape {
    /// The shape a line of this many columns is in, since the count is the
    /// only difference between them. A four column depth file does not exist,
    /// so four columns is bedGraph and three is depth.
    fn sniff(width: usize) -> Option<Shape> {
        // Only the three exact counts are guessed at. A wider file is read as
        // depth when `--format depth` says so, and never by guessing, because
        // four columns is far more often a bedGraph.
        Some(match width {
            4 => Shape::BedGraph,
            3 => Shape::Depth,
            1 => Shape::Values,
            _ => return None,
        })
    }

    /// How many columns a line of this shape has.
    fn width(self) -> usize {
        match self {
            Shape::BedGraph => 4,
            Shape::Depth => 3,
            Shape::Values => 1,
        }
    }

    /// The shape in prose, for the message about a file that changed shape.
    fn name(self) -> &'static str {
        match self {
            Shape::BedGraph => "bedGraph",
            Shape::Depth => "samtools depth",
            Shape::Values => "a bare column of values",
        }
    }

    /// The word `--format` takes for this shape.
    fn word(self) -> &'static str {
        match self {
            Shape::BedGraph => "bedgraph",
            Shape::Depth => "depth",
            Shape::Values => "values",
        }
    }
}

/// A column count with its noun, so that a message says `1 column` and not
/// `1 columns`.
fn count(columns: usize) -> String {
    if columns == 1 {
        "1 column".to_string()
    } else {
        format!("{columns} columns")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    const DEPTH: &str = "\
# samtools depth -a -r NC_000962.3:761100-761104 aln.bam
NC_000962.3\t761100\t12
NC_000962.3\t761101\t14
NC_000962.3\t761102\t0
";

    const BEDGRAPH: &str = "\
track type=bedGraph name=coverage
chr2L\t100\t103\t5
chr2L\t103\t105\t9
";

    #[test]
    fn samtools_depth_is_one_based_so_position_761100_lands_on_761099() {
        let spans = spans(DEPTH, &region("NC_000962.3:761100-761104"), None).unwrap();
        assert_eq!(
            spans,
            vec![
                (761_099, 761_100, 12.0),
                (761_100, 761_101, 14.0),
                (761_101, 761_102, 0.0)
            ]
        );
    }

    #[test]
    fn a_bedgraph_interval_is_zero_based_half_open_and_covers_every_base_in_it() {
        let spans = spans(BEDGRAPH, &region("chr2L:1-200"), None).unwrap();
        // 100 to 103 is three bases and does not include 103, which the next
        // interval starts on. The span says so without being expanded.
        assert_eq!(spans, vec![(100, 103, 5.0), (103, 105, 9.0)]);
    }

    #[test]
    fn a_bare_column_of_values_starts_at_the_left_edge_of_the_region() {
        let text = "0.5\n0.25\n0.75\n";
        let spans = spans(text, &region("Chr4:501-600"), None).unwrap();
        assert_eq!(
            spans,
            vec![(500, 501, 0.5), (501, 502, 0.25), (502, 503, 0.75)]
        );
    }

    #[test]
    fn values_past_the_right_edge_of_the_region_are_left_out() {
        let spans = spans("1\n2\n3\n4\n", &region("Chr4:1-2"), None).unwrap();
        assert_eq!(spans, vec![(0, 1, 1.0), (1, 2, 2.0)]);
    }

    #[test]
    fn a_genome_wide_interval_is_clipped_to_the_region() {
        // One row across a whole chromosome is one span the width of the
        // window, not a million entries. This is the thing that made a forty
        // byte file cost six gigabytes.
        let spans = spans("chrX\t0\t1000000\t3\n", &region("chrX:11-15"), None).unwrap();
        assert_eq!(spans, vec![(10, 15, 3.0)]);
    }

    #[test]
    fn rows_on_another_sequence_or_outside_the_region_are_not_data() {
        let text = "\
chr7\t100\t200\t5
NC_045512.2\t50\t60\t1
NC_045512.2\t100\t101\t7
";
        let spans = spans(text, &region("NC_045512.2:101-110"), None).unwrap();
        assert_eq!(spans, vec![(100, 101, 7.0)]);
    }

    #[test]
    fn a_file_that_changes_shape_names_the_line_it_changed_on() {
        let text = "chrM\t10\t20\t4\nchrM\t21\t4\n";
        let error = spans(text, &region("chrM:1-100"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("bedGraph"), "{error}");
    }

    #[test]
    fn a_line_that_is_none_of_the_three_shapes_says_what_the_three_are() {
        let error = spans("chrM\t10\t20\t4\t+\n", &region("chrM:1-100"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("bedGraph"), "{error}");
        assert!(error.to_string().contains("5 columns"), "{error}");
    }

    #[test]
    fn a_malformed_number_names_its_line_counting_the_ones_that_were_skipped() {
        let text = "#depth\n\nSL2.40ch01\t10\t7\nSL2.40ch01\t11\tNA\n";
        let error = spans(text, &region("SL2.40ch01:1-100"), None).unwrap_err();
        assert_eq!(error.line, 4);
        assert_eq!(error.to_string(), "line 4: depth is not a number: \"NA\"");
    }

    #[test]
    fn a_position_of_zero_in_a_one_based_file_is_an_error_and_not_an_underflow() {
        let error = spans("MT\t0\t5\n", &region("MT:1-100"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("1-based"), "{error}");
    }

    #[test]
    fn an_interval_that_ends_before_it_starts_is_an_error_rather_than_an_empty_range() {
        let error = spans("chr3\t200\t100\t1\n", &region("chr3:1-500"), None).unwrap_err();
        assert_eq!(error.line, 1);
        let error = windows("chr3\t200\t100\t1\n", &region("chr3:1-500")).unwrap_err();
        assert_eq!(error.line, 1);
    }

    #[test]
    fn the_format_flag_is_taken_over_the_shape_of_the_line() {
        // Three columns sniff as depth, so asking for bedGraph has to be the
        // thing that decides, and has to say so when the line cannot be one.
        let error = spans(
            "chr7\t100\t7\n",
            &region("chr7:1-200"),
            Some(Format::BedGraph),
        )
        .unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("--format bedgraph"), "{error}");
    }

    #[test]
    fn samtools_depth_over_several_files_is_read_by_asking_for_it() {
        // `samtools depth a.bam b.bam` writes one depth column per file, and
        // the first sample is the one drawn.
        let text = "amplicon\t1\t3000\t2900\namplicon\t2\t3010\t2880\n";
        let spans = spans(text, &region("amplicon:1-5"), Some(Format::Depth)).unwrap();
        assert_eq!(spans, vec![(0, 1, 3000.0), (1, 2, 3010.0)]);
    }

    #[test]
    fn a_depth_file_of_several_samples_is_not_guessed_at_as_a_bedgraph() {
        // Four columns, so it sniffs as bedGraph, and read that way each depth
        // record becomes a run of three thousand bases at the height of the
        // second sample. The intervals overlap, which no bedGraph does, and
        // that is what says the guess was wrong.
        let text = "amplicon\t1\t3000\t2900\namplicon\t2\t3010\t2880\n";
        let error = spans(text, &region("amplicon:1-5"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("--format depth"), "{error}");
    }

    #[test]
    fn real_depth_output_stops_on_the_first_record_and_still_says_how_to_read_it() {
        // The test above is the case where a depth is larger than the position
        // it sits at, which lets two records overlap. Real output is the other
        // way round: a position of ten thousand and a depth of thirty, so the
        // third column is smaller than the second and the very first record
        // ends before it starts. That refusal has to carry the same way out,
        // because it is the one a reader will actually meet.
        let text = "chr1\t10000\t30\t41\nchr1\t10001\t31\t40\n";
        let error = spans(text, &region("chr1:1-20000"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("--format depth"), "{error}");
    }

    #[test]
    fn insisting_on_bedgraph_gets_the_plain_refusal_and_no_guess() {
        // A reader who named the format is not guessing, so neither is this.
        let text = "chr1\t10000\t30\t41\n";
        let error = spans(text, &region("chr1:1-20000"), Some(Format::BedGraph)).unwrap_err();
        assert!(error.to_string().contains("end is before start"), "{error}");
        assert!(!error.to_string().contains("--format depth"), "{error}");
    }

    #[test]
    fn a_bedgraph_that_is_merely_out_of_order_says_how_to_insist() {
        let text = "chr9\t100\t200\t1\nchr9\t50\t100\t2\n";
        assert!(spans(text, &region("chr9:1-300"), None).is_err());
        let got = spans(text, &region("chr9:1-300"), Some(Format::BedGraph)).unwrap();
        assert_eq!(got, vec![(100, 200, 1.0), (50, 100, 2.0)]);
    }

    #[test]
    fn a_column_of_values_can_be_asked_for_by_name() {
        let got = spans("7\n8\n", &region("ChrUn:3-10"), Some(Format::Values)).unwrap();
        assert_eq!(got, vec![(2, 3, 7.0), (3, 4, 8.0)]);
    }

    #[test]
    fn an_interval_format_is_refused_rather_than_read_as_a_signal() {
        let error = spans(
            "chr1\t10\t20\tgene\n",
            &region("chr1:1-100"),
            Some(Format::Bed),
        )
        .unwrap_err();
        assert_eq!(error.line, 0);
        assert!(
            error.to_string().contains("bedgraph, depth or values"),
            "{error}"
        );
    }

    #[test]
    fn an_empty_file_is_no_pairs_rather_than_an_error() {
        assert!(spans("# nothing here\n", &region("chr1:1-100"), None)
            .unwrap()
            .is_empty());
    }

    const PNPS: &str = "\
track type=bedGraph name=pnps
NC_045512.2\t0\t1000\t1.4
NC_045512.2\t1000\t2000\t0.6
AF086833.2\t0\t1000\t9.9
";

    #[test]
    fn a_window_keeps_the_zero_based_half_open_bounds_the_file_gave_it() {
        let found = windows(PNPS, &region("NC_045512.2:1-2000")).unwrap();
        assert_eq!(
            found,
            vec![Window::new(0, 1000, 1.4), Window::new(1000, 2000, 0.6)]
        );
        // The second window starts where the first ends, and 1000 belongs to
        // the second one.
        assert_eq!(found[1].start, 1000);
    }

    #[test]
    fn a_window_that_does_not_reach_the_region_is_left_out() {
        let text = "chr3\t0\t100\t1\nchr3\t100\t200\t2\nchr3\t900\t1000\t3\n";
        // 101-200 is 100 to 200 with a zero base, so the window ending at 100
        // stops one base short of it.
        let found = windows(text, &region("chr3:101-200")).unwrap();
        assert_eq!(found, vec![Window::new(100, 200, 2.0)]);
    }

    #[test]
    fn a_window_hanging_over_the_edge_of_the_region_keeps_its_own_bounds() {
        let found = windows(
            "Pf3D7_01_v3\t0\t5000\t0.2\n",
            &region("Pf3D7_01_v3:101-200"),
        )
        .unwrap();
        assert_eq!(found, vec![Window::new(0, 5000, 0.2)]);
    }

    #[test]
    fn a_window_file_that_is_not_bedgraph_names_the_line() {
        let error = windows("chr3\t0\t100\n", &region("chr3:1-200")).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(
            error.to_string().contains("chrom start end value"),
            "{error}"
        );
    }

    #[test]
    fn a_window_with_a_value_that_is_not_a_number_names_the_line() {
        let text = "chr3\t0\t100\t1\nchr3\t100\t200\t.\n";
        let error = windows(text, &region("chr3:1-200")).unwrap_err();
        assert_eq!(error.to_string(), "line 2: value is not a number: \".\"");
    }

    /// Every row, on every sequence, each sequence in the order the file
    /// first names it and its spans in the file's order, none cut to a window.
    #[test]
    fn genome_spans_reads_every_sequence_in_file_order() {
        let text = "\
track type=bedGraph
chr2\t0\t100\t5
chr10\t50\t60\t1.5
chr2\t100\t300000000\t7
chr10\t60\t60\t9
";
        let read = genome_spans(text, None).unwrap().sequences;
        assert_eq!(
            read,
            vec![
                (
                    "chr2".to_string(),
                    vec![(0, 100, 5.0), (100, 300_000_000, 7.0)]
                ),
                // An interval of no bases is no span, as over a window.
                ("chr10".to_string(), vec![(50, 60, 1.5)]),
            ]
        );
        // samtools depth, one position a line from one.
        let depth = genome_spans("X\t1\t4\nY\t10\t2\nX\t2\t5\n", None)
            .unwrap()
            .sequences;
        assert_eq!(
            depth,
            vec![
                ("X".to_string(), vec![(0, 1, 4.0), (1, 2, 5.0)]),
                ("Y".to_string(), vec![(9, 10, 2.0)]),
            ]
        );
        let mut count = 0;
        assert_eq!(
            fold_genome_spans(text, None, |_, _, _, _| count += 1).unwrap(),
            3
        );
        assert_eq!(count, 3);
    }

    /// The first row of the next sequence starts back at nought, which is
    /// no overlap; a row that goes back on its own sequence is one, wherever
    /// the rows of other sequences fell between.
    #[test]
    fn overlap_is_checked_within_a_sequence_not_across_them() {
        let sorted = "chr1\t0\t100\t1\nchr1\t100\t200\t2\nchr2\t0\t100\t3\n";
        assert_eq!(genome_spans(sorted, None).unwrap().sequences.len(), 2);
        let back = "chr1\t0\t100\t1\nchr2\t0\t100\t3\nchr1\t50\t150\t2\n";
        let error = genome_spans(back, None).unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.reason.contains("these intervals overlap"), "{error}");
        // Over a window, the rows of another sequence never counted.
        let spans = spans(back, &region("chr2:1-100"), None).unwrap();
        assert_eq!(spans, vec![(0, 100, 3.0)]);
        // And returning to a sequence after another is no overlap by itself.
        let apart = "chr1\t0\t100\t1\nchr2\t0\t100\t3\nchr1\t100\t150\t2\n";
        let read = genome_spans(apart, None).unwrap().sequences;
        assert_eq!(read[0].1, [(0, 100, 1.0), (100, 150, 2.0)]);
    }

    /// A bare value has no sequence to be laid on, and across a genome it
    /// is refused rather than laid on the first.
    #[test]
    fn a_bare_column_has_no_place_on_a_genome() {
        let error = genome_spans("0.5\n0.25\n", None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("names no sequence"), "{error}");
        let asked = genome_spans("0.5\n", Some(Format::Values)).unwrap_err();
        assert!(asked.reason.contains("names no sequence"), "{asked}");
    }

    #[test]
    fn genome_windows_keep_each_sequence_s_windows_whole() {
        let text = "1\t0\t500000\t0.25\n2\t0\t500000\t-0.5\n1\t500000\t900000\t1\n";
        let read = genome_windows(text).unwrap();
        assert_eq!(
            read,
            vec![
                (
                    "1".to_string(),
                    vec![
                        Window::new(0, 500_000, 0.25),
                        Window::new(500_000, 900_000, 1.0)
                    ]
                ),
                ("2".to_string(), vec![Window::new(0, 500_000, -0.5)]),
            ]
        );
        let error = genome_windows("1\t0\t100\n").unwrap_err();
        assert!(error.reason.contains("chrom start end value"), "{error}");
    }
}
