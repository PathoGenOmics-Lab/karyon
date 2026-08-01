//! A value per sample per site, as a table.
//!
//! The header is the one place in this file where the coordinate convention
//! matters: the positions it names are 1-based, the way a VCF writes them and
//! the way every tool that makes such a table prints them, so a site comes out
//! at `position - 1`. Columns outside the region are dropped here rather than
//! carried to the track, since a genome-wide genotype table is mostly not the
//! window on display.

use karyon::{MatrixRow, Region};

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
