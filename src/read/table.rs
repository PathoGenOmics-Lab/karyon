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

/// The windows of a table, each 0-based and half-open, and a row of values
/// per sample, one value per window.
pub type WindowMatrix = (Vec<(u64, u64)>, Vec<MatrixRow>);

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
pub fn windows(text: &str, region: &Region) -> Result<WindowMatrix, ReadError> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut names: Option<Vec<String>> = None;
    let mut spans = Vec::new();
    let mut values: Vec<Vec<f64>> = Vec::new();
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
            if names.is_none() && spans.is_empty() && fields.len() >= 4 {
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
        let wanted = names.as_ref().map_or(values.len(), Vec::len);
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
        if values.is_empty() {
            values = vec![Vec::new(); samples];
        }
        if fields[0].trim() != region.seq() {
            continue;
        }
        let start: u64 = number(fields[1].trim(), "the start", line)?;
        let end: u64 = number(fields[2].trim(), "the end", line)?;
        if end <= start || start >= region.end() || end <= region.start() {
            continue;
        }
        spans.push((start, end));
        for (sample, field) in fields[3..].iter().enumerate() {
            values[sample].push(cell(field, sample + 4, line)?);
        }
    }
    let names = names.unwrap_or_else(|| {
        (0..values.len())
            .map(|column| format!("column {}", column + 4))
            .collect()
    });
    let rows = names
        .into_iter()
        .zip(values.into_iter().chain(std::iter::repeat(Vec::new())))
        .map(|(name, cells)| MatrixRow::new(name, cells))
        .collect();
    Ok((spans, rows))
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
