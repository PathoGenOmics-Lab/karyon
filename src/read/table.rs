//! A value per sample per site, as a table.
//!
//! One row per sample and one column per site: an allele fraction, a genotype
//! written as a number, a per-sample depth. The coordinate convention is the
//! header's, and it counts from one, the way a VCF writes positions and the way
//! every tool that makes such a table prints them, so a site comes out at
//! `position - 1`. Columns outside the region are dropped here rather than
//! carried to the track, since a genome-wide table is mostly not the window on
//! display.
//!
//! # The header is the whole coordinate system
//!
//! No data row says where its values sit. A value's position is the column it
//! is in and nothing else, so the header has to be read before a row means
//! anything, and a row that does not line up with it column for column holds
//! values that belong to no site in particular. That stops the read, naming the
//! sample and both counts.
//!
//! The first field of the header is the corner of the table, and it counts as
//! one exactly when it does not read as a position. Nothing has to be said on
//! the command line about which of the two kinds of table this is.
//!
//! # A blank cell is a statement, not a zero
//!
//! Tables of this shape come with holes in them, a sample that was never typed
//! at a site being ordinary rather than an accident, so empty, `.` and `NA` are
//! all read as no value at all. Filling a hole with a zero would invent a
//! measurement, and the figure would then carry a claim nobody made. Anything
//! else that is not a number stops the read on its line and names the column.

use crate::{MatrixRow, Region};

use super::{columns, lines, number, ReadError};

/// Reads a matrix whose header names the sites and whose first column names
/// the samples.
///
/// The header is the positions, 1-based, with the first field either empty or
/// a word such as `sample`. Every row after it is a name followed by one value
/// per site. A value that is empty, `.` or `NA` is missing rather than zero.
///
/// Returns the site positions, 0-based, and one row per sample.
pub fn matrix(text: &str, region: &Region) -> Result<(Vec<u64>, Vec<MatrixRow>), ReadError> {
    // A VCF's header is comments to this reader, so its first record was read
    // as the header of sites and refused for an `ID` where a position goes,
    // which says nothing of what the file is or what reads it.
    if is_vcf(text) {
        return Err(ReadError::whole(
            "this is a VCF, and --matrix reads a table of sites across a header; \
             --genotypes draws a VCF's samples, a row each",
        ));
    }
    let mut rest = lines(text);
    let Some((header_line, header)) = rest.next() else {
        // A file with no header has no sites and no samples. The caller is the
        // one that decides an empty track is worth an error message.
        return Ok((Vec::new(), Vec::new()));
    };

    let head = columns(header);
    let first = corner(&head);
    let mut sites = Vec::with_capacity(head.len().saturating_sub(first));
    for (index, field) in head.iter().enumerate().skip(first) {
        let column = index + 1;
        let position: u64 = number(
            field.trim(),
            &format!("the site position in column {column}"),
            header_line,
        )?;
        if position == 0 {
            return Err(ReadError::at(
                header_line,
                format!("the site position in column {column} is 0, and the header is 1-based"),
            ));
        }
        // 1-based inclusive to the 0-based coordinate the rest of the crate uses.
        sites.push(position - 1);
    }

    // Which columns survive the region, kept as indices so a row can be cut
    // down to the same columns as the header without comparing anything twice.
    let kept: Vec<usize> = (0..sites.len())
        .filter(|index| region.contains(sites[*index]))
        .collect();
    let inside: Vec<u64> = kept.iter().map(|index| sites[*index]).collect();

    let mut rows = Vec::new();
    for (line, row) in rest {
        let fields = columns(row);
        let name = fields.first().copied().unwrap_or_default().trim();
        let values = fields.get(1..).unwrap_or_default();
        if values.len() != sites.len() {
            return Err(ReadError::at(
                line,
                format!(
                    "sample {:?} has {} values and the header names {} sites",
                    name,
                    values.len(),
                    sites.len()
                ),
            ));
        }
        let mut cells = Vec::with_capacity(inside.len());
        for index in &kept {
            cells.push(cell(values[*index], index + 2, line)?);
        }
        rows.push(MatrixRow::new(name, cells));
    }

    Ok((inside, rows))
}

/// Whether the text opens with a VCF's header: its `##fileformat=VCF` line,
/// or the `#CHROM POS` line that names its columns.
fn is_vcf(text: &str) -> bool {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    text.lines()
        .map(|line| line.trim_end_matches('\r'))
        .take_while(|line| line.starts_with('#') || line.trim().is_empty())
        .any(|line| {
            line.starts_with("##fileformat=VCF")
                || columns(line).get(..2) == Some(&["#CHROM", "POS"][..])
        })
}

/// The windows of a table, each 0-based and half-open, and a row of values
/// per sample, one value per window.
pub type WindowMatrix = (Vec<(u64, u64)>, Vec<MatrixRow>);

/// The windows of one sequence, each 0-based and half-open, with the value of
/// each sample there.
pub type WindowValues = Vec<((u64, u64), Vec<f64>)>;

/// Every window of a table of windows, on every sequence it names, as a
/// heatmap drawn across a whole genome reads it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GenomeWindows {
    /// The samples, in the order of the header, or of a long table's rows.
    pub samples: Vec<String>,
    /// Each sequence the table names, in the order it first names them, with
    /// its windows, each 0-based and half-open on that sequence and in the
    /// order of the file, and the value of each sample there, in the order
    /// of `samples`: missing where the cell is empty, `.` or `NA`, and where
    /// a long table gives the sample no row in the window.
    pub sequences: Vec<(String, WindowValues)>,
}

/// Which rows of a table are read: the ones over one window of one sequence,
/// or every one, on whatever sequence it names.
#[derive(Debug, Clone, Copy)]
enum Rows<'a> {
    In(&'a Region),
    All,
}

impl Rows<'_> {
    /// Whether a window on `sequence` from `start` to `end` is read.
    fn keep(self, sequence: &str, start: u64, end: u64) -> bool {
        match self {
            Rows::In(region) => {
                sequence == region.seq() && start < region.end() && end > region.start()
            }
            Rows::All => true,
        }
    }
}

/// Reads a value per sample per window, as `bedtools unionbedg` writes it.
///
/// Each row is a window, a sequence, a start and an end, 0-based and
/// half-open as BED is, followed by one value per sample. The header names
/// the samples: `chrom start end S1 S2` as `unionbedg -header` writes it, or
/// `#'chr' 'start' 'end' 'S1.bam'` as deepTools does, its quotes and its hash
/// taken off. A file with no header names its samples by their column. A
/// value that is empty, `.` or `NA` is missing rather than zero.
///
/// Returns the windows that reach into the region, and one row per sample
/// with a value per window, in the order of the file.
///
/// A table in the long form, a window and a sample to a row, as `chrom start
/// end sample value`, is read too, as [`long_windows`] reads it: told by a
/// fourth column that names a sample rather than holding a value.
pub fn windows(text: &str, region: &Region) -> Result<WindowMatrix, ReadError> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    if is_long(text) {
        return long_windows(text, region);
    }
    let mut spans = Vec::new();
    let mut values: Vec<Vec<f64>> = Vec::new();
    let names = wide(text, Rows::In(region), |_, span, cells| {
        spans.push(span);
        if values.is_empty() {
            values = vec![Vec::new(); cells.len()];
        }
        for (sample, value) in cells.into_iter().enumerate() {
            values[sample].push(value);
        }
    })?;
    let rows = names
        .into_iter()
        .zip(values.into_iter().chain(std::iter::repeat(Vec::new())))
        .map(|(name, cells)| MatrixRow::new(name, cells))
        .collect();
    Ok((spans, rows))
}

/// Reads every window of a table of windows, on every sequence it names,
/// for a heatmap laid across a whole genome: the table [`windows`] reads,
/// wide or long, with no window left out for being somewhere else.
///
/// # Errors
///
/// What [`windows`] refuses.
///
/// ```
/// use karyon::read::table::genome_windows;
///
/// let text = "chrom\tstart\tend\tS1\tS2\n2\t0\t500\t30\t.\n1\t0\t500\t28\t31\n";
/// let read = genome_windows(text)?;
/// assert_eq!(read.samples, ["S1", "S2"]);
/// assert_eq!(read.sequences[1].0, "1");
/// assert_eq!(read.sequences[1].1[0], ((0, 500), vec![28.0, 31.0]));
/// assert!(read.sequences[0].1[0].1[1].is_nan());
/// # Ok::<(), karyon::read::ReadError>(())
/// ```
pub fn genome_windows(text: &str) -> Result<GenomeWindows, ReadError> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    if is_long(text) {
        let (samples, sequences) = long(text, Rows::All)?;
        let sequences = sequences
            .into_iter()
            .map(|(name, cells)| {
                let (spans, values) = matrixed(&samples, cells)?;
                let windows = spans
                    .into_iter()
                    .enumerate()
                    .map(|(column, span)| (span, values.iter().map(|row| row[column]).collect()))
                    .collect();
                Ok((name, windows))
            })
            .collect::<Result<_, ReadError>>()?;
        return Ok(GenomeWindows { samples, sequences });
    }
    let mut sequences: Vec<(String, WindowValues)> = Vec::new();
    let mut index: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let samples = wide(text, Rows::All, |sequence, span, cells| {
        // A file's rows come a sequence at a time, so the last one is asked
        // about again far more often than any other.
        let at = match sequences.last() {
            Some((last, _)) if last == sequence => sequences.len() - 1,
            _ => *index.entry(sequence.to_string()).or_insert_with(|| {
                sequences.push((sequence.to_string(), Vec::new()));
                sequences.len() - 1
            }),
        };
        sequences[at].1.push((span, cells));
    })?;
    Ok(GenomeWindows { samples, sequences })
}

/// [`windows`] and [`genome_windows`] of a wide table, handing each window
/// `rows` asks for to `each`, with the sequence its row names and the value
/// of each sample there, and returning the names of the samples.
fn wide(
    text: &str,
    rows: Rows<'_>,
    mut each: impl FnMut(&str, (u64, u64), Vec<f64>),
) -> Result<Vec<String>, ReadError> {
    let mut names: Option<Vec<String>> = None;
    // How many values the first window holds, which a table with no header
    // names its samples by.
    let mut width = 0;
    let mut kept = false;
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        let raw = raw.trim_end_matches('\r');
        if raw.trim().is_empty() {
            continue;
        }
        let fields = columns(raw.trim_start_matches('#'));
        let is_data = fields.len() >= 3
            && !raw.starts_with('#')
            && fields[1].trim().parse::<u64>().is_ok()
            && fields[2].trim().parse::<u64>().is_ok();
        if !is_data {
            // The first line that is not a window names the samples, where
            // it names at least one; a comment above it is only a comment.
            if names.is_none() && !kept && fields.len() >= 4 {
                names = Some(
                    fields[3..]
                        .iter()
                        .map(|name| {
                            name.trim()
                                .trim_matches(|c| c == '\'' || c == '"')
                                .to_string()
                        })
                        .collect(),
                );
                continue;
            }
            if raw.starts_with('#') || raw.starts_with("track ") || raw.starts_with("browser ") {
                continue;
            }
            return Err(ReadError::at(
                line,
                "a window is a sequence, a start and an end, then a value per sample",
            ));
        }
        let samples = fields.len() - 3;
        let wanted = names.as_ref().map_or(width, Vec::len);
        if samples == 0 || (wanted > 0 && samples != wanted) {
            return Err(ReadError::at(
                line,
                format!(
                    "this window has {samples} values and {} names {wanted} samples",
                    if names.is_some() {
                        "the header"
                    } else {
                        "the first window"
                    }
                ),
            ));
        }
        if width == 0 {
            width = samples;
        }
        let sequence = fields[0].trim();
        if let Rows::In(region) = rows {
            if sequence != region.seq() {
                continue;
            }
        }
        let start: u64 = number(fields[1].trim(), "the start", line)?;
        let end: u64 = number(fields[2].trim(), "the end", line)?;
        if end <= start || !rows.keep(sequence, start, end) {
            continue;
        }
        let cells = fields[3..]
            .iter()
            .enumerate()
            .map(|(sample, field)| cell(field, sample + 4, line))
            .collect::<Result<Vec<f64>, ReadError>>()?;
        kept = true;
        each(sequence, (start, end), cells);
    }
    Ok(names.unwrap_or_else(|| {
        (0..width)
            .map(|column| format!("column {}", column + 4))
            .collect()
    }))
}

/// The names a column of samples goes by in the header of a long table.
const SAMPLE: &[&str] = &[
    "sample",
    "samples",
    "name",
    "id",
    "sample_id",
    "sample_name",
];

/// Whether a table of windows is in the long form, a window and a sample to a
/// row: its header calls the fourth column a sample, or its first window
/// holds a word there, where a wide table holds a value.
///
/// A long table names its samples on its rows, in the order the file first
/// gives them, so the rows over one window do not say every sample there is
/// or their order, and the command line reads one whole where it reads a
/// wide one through an index.
pub(crate) fn is_long(text: &str) -> bool {
    for raw in text.lines() {
        let raw = raw.trim_end_matches('\r');
        if raw.trim().is_empty() {
            continue;
        }
        let fields = columns(raw.trim_start_matches('#'));
        let is_data = fields.len() >= 3
            && !raw.starts_with('#')
            && fields[1].trim().parse::<u64>().is_ok()
            && fields[2].trim().parse::<u64>().is_ok();
        if !is_data {
            if fields.len() == 5 {
                let named = fields[3].trim().trim_matches(|c| c == '\'' || c == '"');
                if SAMPLE.iter().any(|name| named.eq_ignore_ascii_case(name)) {
                    return true;
                }
            }
            continue;
        }
        let fourth = fields.get(3).map(|field| field.trim()).unwrap_or_default();
        return fields.len() == 5
            && !fourth.is_empty()
            && fourth != "."
            && fourth != "NA"
            && fourth.parse::<f64>().is_err();
    }
    false
}

/// Reads a value per sample per window from a table in the long form: a
/// sequence, a start and an end, 0-based and half-open as BED is, then the
/// sample and its value, one sample of one window to a row.
///
/// The windows are every span a row names, in order along the sequence, and
/// the samples every name, in the order the file first gives them. A sample
/// the file gives no row for in a window is missing there rather than nought,
/// and two rows for one sample in one window, or two windows that overlap
/// without being the same window, stop the read: a heatmap has one cell for
/// each sample in each window, and neither has one answer.
pub fn long_windows(text: &str, region: &Region) -> Result<WindowMatrix, ReadError> {
    let (samples, mut sequences) = long(text, Rows::In(region))?;
    let cells = sequences.pop().map(|(_, cells)| cells).unwrap_or_default();
    let (spans, values) = matrixed(&samples, cells)?;
    let rows = samples
        .into_iter()
        .zip(values)
        .map(|(name, cells)| MatrixRow::new(name, cells))
        .collect();
    Ok((spans, rows))
}

/// One row of a long table: its window, which of the samples it is, its
/// value and the line it is on.
type LongCell = ((u64, u64), usize, f64, usize);

/// The samples a long table names and each sequence's rows of it.
type LongRows = (Vec<String>, Vec<(String, Vec<LongCell>)>);

/// One sequence's windows of a long table and each sample's value in each.
type Columns = (Vec<(u64, u64)>, Vec<Vec<f64>>);

/// The rows of a long table `rows` asks for, each sequence in the order the
/// file first names it with its rows, and the samples they name, in the
/// order the file first gives them.
fn long(text: &str, rows: Rows<'_>) -> Result<LongRows, ReadError> {
    let mut samples: Vec<String> = Vec::new();
    let mut sequences: Vec<(String, Vec<LongCell>)> = Vec::new();
    let mut index: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for (line, raw) in lines(text) {
        let fields = columns(raw.trim_start_matches('#'));
        let is_data = fields.len() >= 3
            && !raw.starts_with('#')
            && fields[1].trim().parse::<u64>().is_ok()
            && fields[2].trim().parse::<u64>().is_ok();
        if !is_data {
            continue;
        }
        if fields.len() != 5 {
            return Err(ReadError::at(
                line,
                format!(
                    "{} columns: a long table is a sequence, a start, an end, a sample and \
                     its value",
                    fields.len()
                ),
            ));
        }
        let sequence = fields[0].trim();
        if let Rows::In(region) = rows {
            if sequence != region.seq() {
                continue;
            }
        }
        let start: u64 = number(fields[1].trim(), "the start", line)?;
        let end: u64 = number(fields[2].trim(), "the end", line)?;
        if end <= start || !rows.keep(sequence, start, end) {
            continue;
        }
        let name = fields[3].trim();
        let sample = match samples.iter().position(|seen| seen == name) {
            Some(at) => at,
            None => {
                samples.push(name.to_string());
                samples.len() - 1
            }
        };
        let row = ((start, end), sample, cell(fields[4], 5, line)?, line);
        let at = *index.entry(sequence.to_string()).or_insert_with(|| {
            sequences.push((sequence.to_string(), Vec::new()));
            sequences.len() - 1
        });
        sequences[at].1.push(row);
    }
    Ok((samples, sequences))
}

/// The windows of one sequence's rows of a long table, in order along it,
/// and each sample's value in each, missing where it has no row.
fn matrixed(samples: &[String], cells: Vec<LongCell>) -> Result<Columns, ReadError> {
    let mut spans: Vec<(u64, u64)> = cells.iter().map(|(span, ..)| *span).collect();
    spans.sort_unstable();
    spans.dedup();
    if let Some(pair) = spans.windows(2).find(|pair| pair[0].1 > pair[1].0) {
        return Err(ReadError::whole(format!(
            "two windows overlap, {}-{} and {}-{}: a heatmap has one cell for each sample \
             in each window",
            pair[0].0, pair[0].1, pair[1].0, pair[1].1
        )));
    }
    let mut values = vec![vec![f64::NAN; spans.len()]; samples.len()];
    let mut given = vec![vec![false; spans.len()]; samples.len()];
    for (span, sample, value, line) in cells {
        let column = spans.binary_search(&span).expect("every span is listed");
        if std::mem::replace(&mut given[sample][column], true) {
            return Err(ReadError::at(
                line,
                format!(
                    "a second value for {} in the window {}-{}",
                    samples[sample], span.0, span.1
                ),
            ));
        }
        values[sample][column] = value;
    }
    Ok((spans, values))
}

/// Where the positions start in the header row.
///
/// The corner of the table is either empty or a word such as `sample`, so a
/// header whose first field reads as a position has no corner at all. That
/// happens on its own: a file separated by spaces instead of tabs cannot carry
/// an empty first field, because the run of whitespace collapses into the
/// separator and the field disappears with it.
fn corner(head: &[&str]) -> usize {
    match head.first() {
        Some(field) => usize::from(field.trim().parse::<u64>().is_err()),
        None => 0,
    }
}

/// Reads one cell, where the three spellings of nothing become a missing value.
///
/// [`MatrixRow`] says missing with a value that is not finite, which is what
/// lets the track leave a hole where a sample was never typed instead of
/// drawing the zero end of the ramp, and those are different claims.
fn cell(field: &str, column: usize, line: usize) -> Result<f64, ReadError> {
    let field = field.trim();
    if field.is_empty() || field == "." || field == "NA" {
        return Ok(f64::NAN);
    }
    number(field, &format!("the value in column {column}"), line)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    #[test]
    fn windows_are_read_as_unionbedg_and_deeptools_write_them() {
        let bedtools = "chrom\tstart\tend\tS1\tS2\nchr1\t0\t100\t5\tNA\n\
                        chr1\t100\t200\t0\t7.5\nchr2\t0\t100\t1\t1\n";
        let (spans, rows) = windows(bedtools, &region("chr1:1-150")).unwrap();
        assert_eq!(spans, vec![(0, 100), (100, 200)]);
        assert_eq!(rows.len(), 2);
        assert_eq!((rows[0].name.as_str(), rows[1].name.as_str()), ("S1", "S2"));
        assert_eq!(rows[1].value(0), None, "NA is missing, not zero");
        assert_eq!(rows[0].value(1), Some(0.0));
        let deeptools = "#'chr'\t'start'\t'end'\t'a.bam'\t'b.bam'\nchr1\t0\t100\t3\t4\n";
        let (_, rows) = windows(deeptools, &region("chr1:1-100")).unwrap();
        assert_eq!(rows[1].name, "b.bam");
        // No header: the samples are named by their column.
        let (_, rows) = windows("chr1\t0\t100\t3\t4\n", &region("chr1:1-100")).unwrap();
        assert_eq!(rows[0].name, "column 4");
        let error = windows(
            "chrom\tstart\tend\tS1\tS2\nchr1\t0\t100\t3\n",
            &region("chr1:1-100"),
        )
        .unwrap_err();
        assert!(
            error.reason.contains("the header names 2"),
            "{}",
            error.reason
        );
    }

    /// A long table, a window and a sample to a row, reads as the wide table
    /// with the same values does, and a window a sample has no row in is
    /// missing there.
    #[test]
    fn a_long_table_reads_as_the_wide_one() {
        let wide = "chrom\tstart\tend\tS1\tS2\nchr1\t0\t100\t5\tNA\n\
                    chr1\t100\t200\t0\t7.5\n";
        let long = "chrom\tstart\tend\tsample\tdepth\nchr1\t100\t200\tS1\t0\n\
                    chr1\t0\t100\tS1\t5\nchr1\t100\t200\tS2\t7.5\nchr2\t0\t100\tS3\t1\n";
        let (wide_spans, wide_rows) = windows(wide, &region("chr1:1-200")).unwrap();
        let (long_spans, long_rows) = windows(long, &region("chr1:1-200")).unwrap();
        assert_eq!(long_spans, wide_spans);
        assert_eq!(long_rows.len(), 2, "S3 is on another sequence");
        for (long, wide) in long_rows.iter().zip(&wide_rows) {
            assert_eq!(long.name, wide.name);
            assert_eq!(
                (long.value(0), long.value(1)),
                (wide.value(0), wide.value(1))
            );
        }
        // Without a header, a word in the fourth column names a sample.
        let bare = "chr1\t0\t100\tS1\t3\nchr1\t0\t100\tS2\t4\n";
        let (_, rows) = windows(bare, &region("chr1:1-100")).unwrap();
        assert_eq!(rows[1].name, "S2");
        // Two answers for one cell, or windows that overlap, are refused.
        let twice = "chr1\t0\t100\tS1\t3\nchr1\t0\t100\tS1\t4\n";
        let error = windows(twice, &region("chr1:1-100")).unwrap_err();
        assert!(
            error.reason.contains("a second value for S1"),
            "{}",
            error.reason
        );
        let overlap = "chr1\t0\t100\tS1\t3\nchr1\t50\t150\tS1\t4\n";
        let error = windows(overlap, &region("chr1:1-200")).unwrap_err();
        assert!(error.reason.contains("overlap"), "{}", error.reason);
    }

    /// Every window of every sequence, wide or long, each sequence in the
    /// order the file first names it and its windows in the file's order,
    /// with each sample's value there in the header's order: the same
    /// windows and values a window over each sequence reads, and the same
    /// refusals.
    #[test]
    fn a_table_of_windows_is_read_across_every_sequence_it_names() {
        let wide = "chrom\tstart\tend\tS1\tS2\n\
                    chr2\t0\t100\t5\tNA\n\
                    chr1\t0\t100\t1\t2\n\
                    chr2\t100\t150\t6\t7\n\
                    chr1\t100\t200\t3\t.\n";
        let read = genome_windows(wide).unwrap();
        assert_eq!(read.samples, ["S1", "S2"]);
        let names: Vec<&str> = read
            .sequences
            .iter()
            .map(|(name, _)| name.as_str())
            .collect();
        assert_eq!(names, ["chr2", "chr1"]);
        let chr2: Vec<(u64, u64)> = read.sequences[0].1.iter().map(|(span, _)| *span).collect();
        assert_eq!(chr2, [(0, 100), (100, 150)]);
        assert_eq!(read.sequences[1].1[0].1, [1.0, 2.0]);
        assert!(read.sequences[0].1[0].1[1].is_nan(), "NA is missing");
        // As each sequence's own window reads it.
        let (spans, rows) = windows(wide, &region("chr1:1-200")).unwrap();
        let across: Vec<(u64, u64)> = read.sequences[1].1.iter().map(|(span, _)| *span).collect();
        assert_eq!(across, spans);
        assert_eq!(rows[0].value(1), Some(read.sequences[1].1[1].1[0]));

        // A long table: a sample with no row on a sequence is missing there.
        let long = "chrom\tstart\tend\tsample\tdepth\n\
                    chr1\t100\t200\tS1\t0\nchr1\t0\t100\tS1\t5\n\
                    chr2\t0\t100\tS2\t1\nchr1\t0\t100\tS2\t7.5\n";
        let read = genome_windows(long).unwrap();
        assert_eq!(read.samples, ["S1", "S2"]);
        assert_eq!(read.sequences[0].0, "chr1");
        assert_eq!(read.sequences[0].1[0].0, (0, 100));
        assert_eq!(read.sequences[0].1[0].1, [5.0, 7.5]);
        assert!(read.sequences[0].1[1].1[1].is_nan(), "S2 has no row there");
        assert!(
            read.sequences[1].1[0].1[0].is_nan(),
            "S1 has no row on chr2"
        );

        // Refused as a window refuses it, wherever the row is.
        let short = "chrom\tstart\tend\tS1\tS2\nchr1\t0\t100\t1\t2\nchr9\t0\t100\t3\n";
        let error = genome_windows(short).unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.reason.contains("the header names 2"), "{error}");
        let overlap = "chr1\t0\t100\tS1\t3\nchr2\t0\t100\tS1\t3\nchr2\t50\t150\tS1\t4\n";
        let error = genome_windows(overlap).unwrap_err();
        assert!(error.reason.contains("overlap"), "{error}");
    }

    #[test]
    fn a_header_position_is_one_based_and_the_site_is_not() {
        // rpoB in Mycobacterium tuberculosis, where 761155 is the codon 450
        // position everyone quotes, and the site behind it is 761154.
        let text = "sample\t761155\t761160\nERR001\t1\t0\n";
        let (sites, rows) = matrix(text, &region("NC_000962.3:761,000-761,200")).unwrap();
        assert_eq!(sites, vec![761_154, 761_159]);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].name, "ERR001");
        assert_eq!(rows[0].value(0), Some(1.0));
    }

    /// A VCF's header is comments here, so its first record was read as the
    /// header and refused for an `ID` where a position goes, which said
    /// nothing of what the file is. It names the track that reads one.
    #[test]
    fn a_vcf_is_refused_and_the_track_for_it_named() {
        let region = region("chr1:1-1000");
        let vcf = "##fileformat=VCFv4.2\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                   chr1\t100\trs1\tC\tT\t.\t.\t.\tGT\t1\n";
        let headless = &vcf[vcf.find("#CHROM").unwrap()..];
        for text in [vcf, headless] {
            let error = matrix(text, &region).unwrap_err();
            assert_eq!(error.line, 0, "the whole file, not one line");
            assert!(error.reason.contains("--genotypes"), "{}", error.reason);
        }
        // A comment above a table is still a table.
        let (sites, _) = matrix("# typed in 2026\nsample\t101\nS1\t1\n", &region).unwrap();
        assert_eq!(sites, vec![100]);
    }

    #[test]
    fn an_empty_corner_is_a_corner_and_so_is_a_word() {
        let named = "sample\t101\t201\nBY4741\t1\t0\n";
        let bare = "\t101\t201\nBY4741\t1\t0\n";
        let region = region("chrIV:1-1000");
        let (left, _) = matrix(named, &region).unwrap();
        let (right, _) = matrix(bare, &region).unwrap();
        assert_eq!(left, vec![100, 200]);
        assert_eq!(right, left);
    }

    #[test]
    fn a_space_separated_header_with_no_corner_still_names_the_sites() {
        // Spaces collapse, so an empty corner cannot survive the separator.
        // A header that starts with a position is therefore all positions.
        let text = "101 201\nCol-0 1 0\nLer-0 0 1\n";
        let (sites, rows) = matrix(text, &region("Chr5:1-1000")).unwrap();
        assert_eq!(sites, vec![100, 200]);
        assert_eq!(rows[1].name, "Ler-0");
        assert_eq!(rows[1].value(1), Some(1.0));
    }

    #[test]
    fn empty_dot_and_na_are_missing_rather_than_zero() {
        let text = "sample\t101\t201\t301\t401\nEPI_0001\t0\t.\tNA\t\n";
        let (_, rows) = matrix(text, &region("NC_045512.2:1-1000")).unwrap();
        let row = &rows[0];
        // A typed zero is a value and the three spellings of nothing are not.
        assert_eq!(row.value(0), Some(0.0));
        assert_eq!(row.value(1), None);
        assert_eq!(row.value(2), None);
        assert_eq!(row.value(3), None);
        assert!(row.values[1].is_nan());
    }

    #[test]
    fn a_row_of_the_wrong_length_names_the_line_and_both_counts() {
        let text = "sample\t101\t201\t301\nERR001\t1\t0\t1\nERR002\t1\t0\n";
        let error = matrix(text, &region("NC_000962.3:1-1000")).unwrap_err();
        assert_eq!(error.line, 3);
        assert_eq!(
            error.to_string(),
            "line 3: sample \"ERR002\" has 2 values and the header names 3 sites"
        );
    }

    #[test]
    fn a_site_outside_the_region_takes_its_column_with_it() {
        let text = "sample\t101\t9001\t201\nDSM20231\t1\t9\t0\n";
        let (sites, rows) = matrix(text, &region("NZ_CP007433.1:1-1000")).unwrap();
        assert_eq!(sites, vec![100, 200]);
        // The value that belonged to the dropped column is dropped with it, so
        // the row still lines up with the sites column for column.
        assert_eq!(rows[0].values.len(), 2);
        assert_eq!(rows[0].value(0), Some(1.0));
        assert_eq!(rows[0].value(1), Some(0.0));
    }

    #[test]
    fn comments_and_blank_lines_do_not_move_the_line_numbers() {
        let text = "# a genotype matrix\n\nsample\t101\t201\nERR001\t1\tx\n";
        let error = matrix(text, &region("NC_000962.3:1-1000")).unwrap_err();
        assert_eq!(error.line, 4);
        assert!(error.to_string().contains("the value in column 3"));
    }

    #[test]
    fn a_header_field_that_is_not_a_position_is_an_error_on_the_header_line() {
        let text = "sample\t101\tchr1:201\nERR001\t1\t0\n";
        let error = matrix(text, &region("NC_000962.3:1-1000")).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("the site position in column 3"));
    }

    #[test]
    fn a_zero_in_a_one_based_header_is_an_error_and_not_an_underflow() {
        let text = "sample\t0\t201\nERR001\t1\t0\n";
        let error = matrix(text, &region("NC_000962.3:1-1000")).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("1-based"), "{error}");
    }

    #[test]
    fn a_file_with_nothing_in_it_is_not_an_error_here() {
        let (sites, rows) = matrix("\n# only a comment\n", &region("chr1:1-1000")).unwrap();
        assert!(sites.is_empty());
        assert!(rows.is_empty());
    }

    #[test]
    fn a_name_with_a_space_in_it_survives_a_tab_separated_file() {
        let text = "sample\t101\nisolate 12\t1\n";
        let (_, rows) = matrix(text, &region("chr1:1-1000")).unwrap();
        assert_eq!(rows[0].name, "isolate 12");
    }
}
