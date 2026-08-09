//! Per-base signal: bedGraph, `samtools depth`, and a bare column of values.
//!
//! Three formats that give a value to every base, read as `(position, value)`
//! pairs for a coverage track and as spans for a window track, which takes
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

use crate::{Region, Window};

use super::Format;
use super::{columns, lines, number, ReadError};

/// Reads a value per position, as 0-based `(position, value)` pairs.
///
/// Pairs rather than a dense array, because
/// [`CoverageTrack::from_pairs`](crate::CoverageTrack::from_pairs) already
/// lays them over the region and a sparse file should not have to be widened
/// on the way in. A position the file says nothing about stays at zero, which
/// is what a depth of zero means and what a bedGraph leaves out.
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
pub fn dense(
    text: &str,
    region: &Region,
    format: Option<Format>,
) -> Result<Vec<(u64, f64)>, ReadError> {
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
    let mut pairs: Vec<(u64, f64)> = Vec::new();
    // Where the next bare value lands, that shape carrying no position of its own.
    let mut next = region.start();
    // How far the last bedGraph interval reached, which is what says whether
    // the intervals overlap and so whether this is a bedGraph at all.
    let mut reached: Option<u64> = None;

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
                if fields[0] != region.seq() {
                    continue;
                }
                let start: u64 = number(fields[1], "start", at)?;
                let end: u64 = number(fields[2], "end", at)?;
                if end < start {
                    return Err(ReadError::at(at, "end is before start"));
                }
                // Two bedGraph intervals never overlap. `samtools depth` over
                // two files has four columns too, and reading one as a bedGraph
                // turns each depth record into a run of bases at the height of
                // the second sample, which is a plausible looking figure of
                // nothing. Overlap is what tells them apart, so it is refused
                // here rather than guessed at.
                if asked.is_none() {
                    if let Some(previous) = reached {
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
                    reached = Some(end);
                }
                let value: f64 = number(fields[3], "value", at)?;
                // An interval covers every base in it, so a genome-wide file
                // would otherwise become one pair per base of the genome. The
                // clip is what keeps that from happening.
                let from = start.max(region.start());
                let to = end.min(region.end());
                for pos in from..to {
                    pairs.push((pos, value));
                }
            }
            Shape::Depth => {
                if fields[0] != region.seq() {
                    continue;
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
                if region.contains(pos) {
                    pairs.push((pos, depth));
                }
            }
            Shape::Values => {
                let value: f64 = number(fields[0], "value", at)?;
                let pos = next;
                next += 1;
                if region.contains(pos) {
                    pairs.push((pos, value));
                }
            }
        }
    }

    Ok(pairs)
}

/// Reads intervals with a value each, for a window track.
///
/// bedGraph, 0-based half-open, kept as intervals rather than flattened to one
/// value per base, since a window track draws the window and not the base.
pub fn windows(text: &str, region: &Region) -> Result<Vec<Window>, ReadError> {
    let mut found = Vec::new();
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
        if fields[0] != region.seq() {
            continue;
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
        if end <= region.start() || start >= region.end() {
            continue;
        }
        found.push(Window::new(start, end, value));
    }
    Ok(found)
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
        let pairs = dense(DEPTH, &region("NC_000962.3:761100-761104"), None).unwrap();
        assert_eq!(
            pairs,
            vec![(761_099, 12.0), (761_100, 14.0), (761_101, 0.0)]
        );
    }

    #[test]
    fn a_bedgraph_interval_is_zero_based_half_open_and_covers_every_base_in_it() {
        let pairs = dense(BEDGRAPH, &region("chr2L:1-200"), None).unwrap();
        // 100 to 103 is three bases and does not include 103, which the next
        // interval starts on.
        assert_eq!(
            pairs,
            vec![(100, 5.0), (101, 5.0), (102, 5.0), (103, 9.0), (104, 9.0),]
        );
    }

    #[test]
    fn a_bare_column_of_values_starts_at_the_left_edge_of_the_region() {
        let text = "0.5\n0.25\n0.75\n";
        let pairs = dense(text, &region("Chr4:501-600"), None).unwrap();
        assert_eq!(pairs, vec![(500, 0.5), (501, 0.25), (502, 0.75)]);
    }

    #[test]
    fn values_past_the_right_edge_of_the_region_are_left_out() {
        let pairs = dense("1\n2\n3\n4\n", &region("Chr4:1-2"), None).unwrap();
        assert_eq!(pairs, vec![(0, 1.0), (1, 2.0)]);
    }

    #[test]
    fn a_genome_wide_interval_is_clipped_to_the_region() {
        let pairs = dense("chrX\t0\t1000000\t3\n", &region("chrX:11-15"), None).unwrap();
        assert_eq!(
            pairs,
            vec![(10, 3.0), (11, 3.0), (12, 3.0), (13, 3.0), (14, 3.0)]
        );
    }

    #[test]
    fn rows_on_another_sequence_or_outside_the_region_are_not_data() {
        let text = "\
chr7\t100\t200\t5
NC_045512.2\t50\t60\t1
NC_045512.2\t100\t101\t7
";
        let pairs = dense(text, &region("NC_045512.2:101-110"), None).unwrap();
        assert_eq!(pairs, vec![(100, 7.0)]);
    }

    #[test]
    fn a_file_that_changes_shape_names_the_line_it_changed_on() {
        let text = "chrM\t10\t20\t4\nchrM\t21\t4\n";
        let error = dense(text, &region("chrM:1-100"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("bedGraph"), "{error}");
    }

    #[test]
    fn a_line_that_is_none_of_the_three_shapes_says_what_the_three_are() {
        let error = dense("chrM\t10\t20\t4\t+\n", &region("chrM:1-100"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("bedGraph"), "{error}");
        assert!(error.to_string().contains("5 columns"), "{error}");
    }

    #[test]
    fn a_malformed_number_names_its_line_counting_the_ones_that_were_skipped() {
        let text = "#depth\n\nSL2.40ch01\t10\t7\nSL2.40ch01\t11\tNA\n";
        let error = dense(text, &region("SL2.40ch01:1-100"), None).unwrap_err();
        assert_eq!(error.line, 4);
        assert_eq!(error.to_string(), "line 4: depth is not a number: \"NA\"");
    }

    #[test]
    fn a_position_of_zero_in_a_one_based_file_is_an_error_and_not_an_underflow() {
        let error = dense("MT\t0\t5\n", &region("MT:1-100"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("1-based"), "{error}");
    }

    #[test]
    fn an_interval_that_ends_before_it_starts_is_an_error_rather_than_an_empty_range() {
        let error = dense("chr3\t200\t100\t1\n", &region("chr3:1-500"), None).unwrap_err();
        assert_eq!(error.line, 1);
        let error = windows("chr3\t200\t100\t1\n", &region("chr3:1-500")).unwrap_err();
        assert_eq!(error.line, 1);
    }

    #[test]
    fn the_format_flag_is_taken_over_the_shape_of_the_line() {
        // Three columns sniff as depth, so asking for bedGraph has to be the
        // thing that decides, and has to say so when the line cannot be one.
        let error = dense(
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
        let pairs = dense(text, &region("amplicon:1-5"), Some(Format::Depth)).unwrap();
        assert_eq!(pairs, vec![(0, 3000.0), (1, 3010.0)]);
    }

    #[test]
    fn a_depth_file_of_several_samples_is_not_guessed_at_as_a_bedgraph() {
        // Four columns, so it sniffs as bedGraph, and read that way each depth
        // record becomes a run of three thousand bases at the height of the
        // second sample. The intervals overlap, which no bedGraph does, and
        // that is what says the guess was wrong.
        let text = "amplicon\t1\t3000\t2900\namplicon\t2\t3010\t2880\n";
        let error = dense(text, &region("amplicon:1-5"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("--format depth"), "{error}");
    }

    #[test]
    fn a_bedgraph_that_is_merely_out_of_order_says_how_to_insist() {
        let text = "chr9\t100\t200\t1\nchr9\t50\t100\t2\n";
        assert!(dense(text, &region("chr9:1-300"), None).is_err());
        let pairs = dense(text, &region("chr9:1-300"), Some(Format::BedGraph)).unwrap();
        assert_eq!(pairs.len(), 150);
    }

    #[test]
    fn a_column_of_values_can_be_asked_for_by_name() {
        let pairs = dense("7\n8\n", &region("ChrUn:3-10"), Some(Format::Values)).unwrap();
        assert_eq!(pairs, vec![(2, 7.0), (3, 8.0)]);
    }

    #[test]
    fn an_interval_format_is_refused_rather_than_read_as_a_signal() {
        let error = dense(
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
        assert!(dense("# nothing here\n", &region("chr1:1-100"), None)
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
}
