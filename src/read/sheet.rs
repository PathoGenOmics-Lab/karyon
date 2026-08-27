//! A sample sheet: one row per named thing, one column per fact about it.
//!
//! ```text
//! sample  lineage  host    depth
//! S001    L4       human   72.5
//! S002    L2       bovine  61.0
//! S003    L4               48.2
//! ```
//!
//! The first line is the header, the first column is the name, and every other
//! column is an attribute the header names. Nothing here is a coordinate and
//! nothing here is drawn on its own: the rows join by name to the rows of a
//! track that already has some, and become the strips beside them.
//!
//! # The first line is the header, and there is no way to say otherwise
//!
//! Every other reader in this module can tell its rows from its headers by
//! shape, because a BED row has a number in column two and a GFF3 row has nine
//! columns. A sheet has no shape at all: any row is a name and some words, and
//! so is the header. So the header is not detected, it is a rule, and a file
//! whose first line is data loses that line to the column names, which is
//! visible in the drawn figure rather than silent.
//!
//! # A blank field is not a value
//!
//! A sample with no entry for a column has no entry, and the strip beside it is
//! drawn as absent rather than as a category. `NA` is read the same way, and so
//! is anything that parses as a number and is not one: a column of depths with
//! `NaN` in it is a column with a depth missing, not a column with a sample
//! whose depth is the word NaN.
//!
//! The cost is one honest ambiguity. A categorical column whose levels really
//! are the two-letter codes for continents has a level `NA` that is read as a
//! missing value, and the only thing that separates those two files is what
//! they mean.
//!
//! # What is refused
//!
//! A row whose field count is not the header's, because a sheet is a rectangle
//! and a short row is a file that was cut rather than a sample about which less
//! is known. A repeated name, because a strip drawn from the second row of a
//! pair says nothing about which of them it came from. A repeated or empty
//! column name, for the same reason a repeated sample name is refused. And a
//! header of one column, which names things and says nothing about them.

use std::collections::{BTreeMap, BTreeSet};

use crate::tree::{AnnotationValue, Annotations};

use super::{columns, lines, ReadError};

/// The rows a sheet holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Sheet {
    /// Each name and what the sheet knows about it.
    ///
    /// Keyed by name rather than ordered by row, because the order a figure
    /// draws these in belongs to the track that has the rows, and a phylogeny
    /// attached to that track will have reordered them already.
    pub rows: BTreeMap<String, Annotations>,
    /// The attribute columns, in the order the header named them.
    pub columns: Vec<String>,
    /// Rows in the file, before anything is joined to anything.
    pub records: usize,
    /// Fields that held no value, over every row and column.
    ///
    /// The figure draws these as absent, so the number is what says whether a
    /// sheet of mostly empty strips is a sheet of unknowns or a file read with
    /// the wrong separator.
    pub blank: usize,
}

impl Sheet {
    /// The names the sheet holds, in order.
    pub fn names(&self) -> impl Iterator<Item = &str> + '_ {
        self.rows.keys().map(String::as_str)
    }

    /// How many of `rows` the sheet has an entry for.
    ///
    /// The number a caller checks before drawing: a sheet that names none of
    /// the rows it was joined to draws every strip as absent, which is a figure
    /// that looks finished and is not.
    pub fn covers<'a>(&self, rows: impl IntoIterator<Item = &'a str>) -> usize {
        rows.into_iter()
            .filter(|name| self.rows.contains_key(*name))
            .count()
    }
}

/// Reads a sample sheet out of a tab or whitespace separated table.
///
/// The first line is the header. Its first field names the column the row
/// names are in and is not itself an attribute; every other field is an
/// attribute, kept in the order it was written.
///
/// A field is a number when it parses as one, `true` or `false` when it spells
/// one, and text otherwise. An empty field, the word `NA`, and anything that
/// parses to a number that is not a number are all absent, which is a state of
/// its own and not a value.
///
/// # Errors
///
/// Returns the line that is not part of a rectangle: a header of one column, a
/// repeated or empty column name, a row whose field count is not the header's,
/// or a name the sheet has already used.
pub fn sheet(text: &str) -> Result<Sheet, ReadError> {
    let mut rows = lines(text);
    let (at, head) = rows
        .next()
        .ok_or_else(|| ReadError::whole("a sheet begins with a header and this file is empty"))?;

    let head = columns(head);
    if head.len() < 2 {
        return Err(ReadError::at(
            at,
            format!(
                "a sheet names things in column one and says something about them in the rest, and this header has {} column{}",
                head.len(),
                if head.len() == 1 { "" } else { "s" }
            ),
        ));
    }

    let mut seen = BTreeSet::new();
    for name in &head[1..] {
        let name = name.trim();
        if name.is_empty() {
            return Err(ReadError::at(
                at,
                "a column with no name cannot be asked for",
            ));
        }
        if !seen.insert(name.to_string()) {
            return Err(ReadError::at(
                at,
                format!("two columns are both called {name:?}"),
            ));
        }
    }
    let names: Vec<String> = head[1..]
        .iter()
        .map(|name| name.trim().to_string())
        .collect();

    let mut found = Sheet {
        rows: BTreeMap::new(),
        columns: names.clone(),
        records: 0,
        blank: 0,
    };

    for (at, line) in rows {
        let cols = columns(line);
        if cols.len() != head.len() {
            return Err(ReadError::at(
                at,
                format!(
                    "the header has {} columns and this row has {}",
                    head.len(),
                    cols.len()
                ),
            ));
        }
        found.records += 1;

        let name = cols[0].trim().to_string();
        if name.is_empty() {
            return Err(ReadError::at(at, "a row with no name joins to nothing"));
        }

        let mut held = Annotations::new();
        for (column, field) in names.iter().zip(&cols[1..]) {
            match value(field) {
                Some(value) => {
                    held.insert(column.clone(), value);
                }
                None => found.blank += 1,
            }
        }

        if found.rows.insert(name.clone(), held).is_some() {
            return Err(ReadError::at(at, format!("{name:?} is named twice")));
        }
    }

    Ok(found)
}

/// One field, or `None` where the field says there is nothing.
fn value(field: &str) -> Option<AnnotationValue> {
    let field = field.trim();
    if field.is_empty() || field == "NA" {
        return None;
    }
    if let Ok(number) = field.parse::<f64>() {
        // A field that parses to a value that is not a number is a field
        // saying its number is missing, which the empty field says too. Kept
        // as a number it would reach a strip and a tooltip spelling `NaN`,
        // and a reader would have to know that this crate writes an absent
        // annotation the same way to tell the two apart.
        if number.is_nan() {
            return None;
        }
        return Some(AnnotationValue::Number(number));
    }
    match field {
        "true" => Some(AnnotationValue::Boolean(true)),
        "false" => Some(AnnotationValue::Boolean(false)),
        other => Some(AnnotationValue::Text(other.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SHEET: &str = "\
sample\tlineage\thost\tdepth
S001\tL4\thuman\t72.5
S002\tL2\tbovine\t61
S003\tL4\t\t48.2
";

    #[test]
    fn the_first_line_names_the_columns_and_the_rest_are_rows() {
        let found = sheet(SHEET).expect("a sheet");
        assert_eq!(found.columns, ["lineage", "host", "depth"]);
        assert_eq!(found.records, 3);
        assert_eq!(found.rows.len(), 3);
        assert_eq!(
            found.rows["S001"]["lineage"],
            AnnotationValue::Text("L4".to_string())
        );
        assert_eq!(found.rows["S002"]["depth"], AnnotationValue::Number(61.0));
    }

    #[test]
    fn a_blank_field_leaves_the_row_without_that_column() {
        // Not an empty string and not a zero. The strip beside S003 has to be
        // drawn as absent, and the only way a drawer can know that is for the
        // value to be missing rather than to be something meaning missing.
        let found = sheet(SHEET).expect("a sheet");
        assert!(!found.rows["S003"].contains_key("host"));
        assert_eq!(found.blank, 1);
    }

    #[test]
    fn na_and_a_number_that_is_not_a_number_are_both_absent() {
        let found = sheet("id\tx\na\tNA\nb\tNaN\nc\t3\n").expect("a sheet");
        assert!(!found.rows["a"].contains_key("x"));
        assert!(!found.rows["b"].contains_key("x"));
        assert_eq!(found.rows["c"]["x"], AnnotationValue::Number(3.0));
        assert_eq!(found.blank, 2);
    }

    #[test]
    fn a_field_is_a_number_a_flag_or_text_in_that_order() {
        let found = sheet("id\tv\na\t2.5\nb\ttrue\nc\tfalse\nd\tL4\n").expect("a sheet");
        assert_eq!(found.rows["a"]["v"], AnnotationValue::Number(2.5));
        assert_eq!(found.rows["b"]["v"], AnnotationValue::Boolean(true));
        assert_eq!(found.rows["c"]["v"], AnnotationValue::Boolean(false));
        assert_eq!(
            found.rows["d"]["v"],
            AnnotationValue::Text("L4".to_string())
        );
    }

    #[test]
    fn a_row_that_is_not_the_headers_width_is_refused() {
        // A short row is a file that was cut, not a sample less is known
        // about: the blank field already says that, and says it in a column.
        let error = sheet("id\ta\tb\nS1\tx\n").unwrap_err();
        assert_eq!(error.line, 2);
        assert!(
            error.reason.contains("3 columns and this row has 2"),
            "{error}"
        );

        let error = sheet("id\ta\nS1\tx\ty\n").unwrap_err();
        assert!(
            error.reason.contains("2 columns and this row has 3"),
            "{error}"
        );
    }

    #[test]
    fn a_name_used_twice_is_refused() {
        let error = sheet("id\ta\nS1\tx\nS1\ty\n").unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.reason.contains("named twice"), "{error}");
    }

    #[test]
    fn a_header_that_names_nothing_or_names_it_twice_is_refused() {
        let error = sheet("id\nS1\n").unwrap_err();
        assert!(error.reason.contains("this header has 1 column"), "{error}");

        let error = sheet("id\ta\ta\nS1\tx\ty\n").unwrap_err();
        assert!(error.reason.contains("both called \"a\""), "{error}");

        let empty = sheet("id\ta\t\nS1\tx\ty\n").unwrap_err();
        assert!(empty.reason.contains("no name"), "{empty}");
    }

    #[test]
    fn an_empty_file_is_a_file_without_a_header() {
        let error = sheet("").unwrap_err();
        assert!(error.reason.contains("empty"), "{error}");
        let error = sheet("# only a comment\n").unwrap_err();
        assert!(error.reason.contains("empty"), "{error}");
    }

    #[test]
    fn how_many_rows_a_sheet_covers_is_what_says_a_join_found_nothing() {
        let found = sheet(SHEET).expect("a sheet");
        assert_eq!(found.covers(["S001", "S002", "nothing"]), 2);
        assert_eq!(found.covers(["ERR001", "ERR002"]), 0);
    }

    #[test]
    fn spaces_separate_a_sheet_that_has_no_tabs_in_it() {
        let found = sheet("sample lineage\nS1 L4\n").expect("a sheet");
        assert_eq!(found.columns, ["lineage"]);
        assert_eq!(
            found.rows["S1"]["lineage"],
            AnnotationValue::Text("L4".to_string())
        );
    }
}
