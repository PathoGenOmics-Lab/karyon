//! Pairs of places and a value between them: linkage, contacts, epistasis.
//!
//! Three shapes, told apart by their first line.
//!
//! PLINK's `--r2` table, `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`, its columns
//! found by those names, so `--r` and `D'` tables read too. The positions are
//! 1-based, as PLINK writes them.
//!
//! BEDPE, six columns of two stretches, `chrom1 start1 end1 chrom2 start2
//! end2`, 0-based and half-open as BED is, as `cooler dump --join` writes a
//! contact map and loop callers write their loops. The value is the first
//! number after the six, a count or a score, and a pair with none is drawn at
//! one.
//!
//! A table of your own, its columns found by what its header calls them:
//! two positions (`pos1` and `pos2`, `site_a` and `site_b`...), a value (`r2`,
//! `score`, `count`...) and, where there is one, a sequence. Positions there
//! count from one. Without a header, three columns are two positions and a
//! value.

use crate::region::Region;
use crate::track::pairs::Pair;

use super::{fields, lines, number, ReadError};

/// The names a column of first positions goes by.
const FIRST: &[&str] = &[
    "pos1",
    "pos_1",
    "pos_a",
    "posa",
    "position1",
    "site1",
    "site_1",
    "site_a",
    "bp_a",
    "i",
    "first",
    "start1",
];
/// The names a column of second positions goes by.
const SECOND: &[&str] = &[
    "pos2",
    "pos_2",
    "pos_b",
    "posb",
    "position2",
    "site2",
    "site_2",
    "site_b",
    "bp_b",
    "j",
    "second",
    "start2",
];
/// The names a column of values goes by.
const VALUE: &[&str] = &[
    "r2",
    "r^2",
    "r²",
    "r",
    "dp",
    "d'",
    "value",
    "score",
    "count",
    "weight",
    "strength",
    "mi",
    "correlation",
    "frequency",
    "contacts",
];
/// The names a column of sequences goes by.
const SEQUENCE: &[&str] = &[
    "chr_a",
    "chrom",
    "chr",
    "chromosome",
    "seqname",
    "sequence",
    "chrom1",
];

/// The column whose header is one of `names`, compared without case.
fn column(header: &[&str], names: &[&str]) -> Option<usize> {
    header.iter().position(|field| {
        names
            .iter()
            .any(|name| field.trim().eq_ignore_ascii_case(name))
    })
}

/// Whether a column of values called `name` is a correlation, which runs
/// from nought to one whatever the data: r², r, D'.
pub fn is_correlation(name: &str) -> bool {
    ["r2", "r^2", "r²", "r", "dp", "d'"]
        .iter()
        .any(|known| name.trim().eq_ignore_ascii_case(known))
}

/// The names PLINK gives the places of its pairs, as `(position, name)` with
/// the position 0-based, from the `SNP_A` and `SNP_B` columns beside `BP_A`
/// and `BP_B`. Empty for a table with no such columns, which names nothing.
///
/// A table of linkage with a lead, as `--ld-snp` writes one, names the lead
/// on every row, and a scan whose own table names no variant can still call
/// its lead by name.
pub fn names(text: &str) -> Vec<(u64, String)> {
    let mut rows = lines(text);
    let Some((_, first)) = rows.next() else {
        return Vec::new();
    };
    let head = fields(first.trim());
    let sides = [
        (column(&head, &["bp_a"]), column(&head, &["snp_a"])),
        (column(&head, &["bp_b"]), column(&head, &["snp_b"])),
    ];
    let mut out: std::collections::BTreeMap<u64, String> = std::collections::BTreeMap::new();
    for (_, row) in rows {
        let cells = fields(row.trim());
        for (at, name) in sides {
            let (Some(at), Some(name)) = (at, name) else {
                continue;
            };
            let position = cells
                .get(at)
                .and_then(|field| field.trim().parse::<u64>().ok());
            let (Some(position), Some(name)) = (position.filter(|p| *p > 0), cells.get(name))
            else {
                continue;
            };
            out.entry(position - 1)
                .or_insert_with(|| name.trim().to_string());
        }
    }
    out.into_iter().collect()
}

/// Reads the pairs with both places on the region's sequence, and what the
/// header calls their values, where it names them.
///
/// A pair with one place on another sequence is left out, as is any pair on
/// another sequence: a contact map of a whole genome is handed over to draw
/// one chromosome. A value that is empty, `.`, `NA` or `nan` is a pair with no
/// answer, kept and not drawn.
pub fn pairs(text: &str, region: &Region) -> Result<(Vec<Pair>, Option<String>), ReadError> {
    let mut rows = lines(text).peekable();
    let Some(&(first_line, first)) = rows.peek() else {
        return Ok((Vec::new(), None));
    };
    let head = fields(first.trim());
    let numeric = |field: &str| field.trim().parse::<f64>().is_ok();
    let bedpe = head.len() >= 6
        && numeric(head[1])
        && numeric(head[2])
        && numeric(head[4])
        && numeric(head[5]);
    if bedpe {
        return bedpe_pairs(rows, region).map(|pairs| (pairs, None));
    }
    let has_header = !head.iter().take(2).all(|field| numeric(field));
    let (at_first, at_second, at_value, at_sequence) = if has_header {
        rows.next();
        let (Some(a), Some(b)) = (column(&head, FIRST), column(&head, SECOND)) else {
            return Err(ReadError::at(
                first_line,
                format!(
                    "no two columns of positions: they are found by their header, as \
                     pos1 and pos2, BP_A and BP_B or site_a and site_b, and this header is {}",
                    head.join(" ")
                ),
            ));
        };
        let value = column(&head, VALUE);
        (a, b, value, column(&head, SEQUENCE))
    } else {
        (0, 1, Some(2), None)
    };
    // PLINK names the second place's sequence apart from the first's.
    let at_second_sequence = has_header
        .then(|| column(&head, &["chr_b", "chrom2"]))
        .flatten();
    let mut out = Vec::new();
    for (line, row) in rows {
        let cells = fields(row.trim());
        let field = |at: usize| cells.get(at).copied().unwrap_or_default().trim();
        let on_sequence = |at: Option<usize>| at.map_or(true, |at| field(at) == region.seq());
        if !on_sequence(at_sequence) || !on_sequence(at_second_sequence) {
            continue;
        }
        let a: u64 = number(field(at_first), "the first position", line)?;
        let b: u64 = number(field(at_second), "the second position", line)?;
        if a == 0 || b == 0 {
            return Err(ReadError::at(
                line,
                "a position of 0: positions in this table count from 1",
            ));
        }
        let value = match at_value {
            Some(at) => value(field(at), line)?,
            None => 1.0,
        };
        out.push(Pair::new(a - 1, b - 1, value));
    }
    let named = at_value
        .filter(|_| has_header)
        .map(|at| head[at].trim().to_string());
    Ok((out, named))
}

/// Reads BEDPE rows: two stretches, 0-based and half-open, and a value.
fn bedpe_pairs<'a>(
    rows: impl Iterator<Item = (usize, &'a str)>,
    region: &Region,
) -> Result<Vec<Pair>, ReadError> {
    let mut out = Vec::new();
    for (line, row) in rows {
        let fields = fields(row.trim());
        if fields.len() < 6 {
            return Err(ReadError::at(
                line,
                format!(
                    "{} columns, and a BEDPE row is two stretches of three",
                    fields.len()
                ),
            ));
        }
        if fields[0].trim() != region.seq() || fields[3].trim() != region.seq() {
            continue;
        }
        let at = |index: usize, what: &str| number::<u64>(fields[index].trim(), what, line);
        let first = (at(1, "the first start")?, at(2, "the first end")?);
        let second = (at(4, "the second start")?, at(5, "the second end")?);
        let value = fields[6..]
            .iter()
            .find_map(|field| field.trim().parse::<f64>().ok())
            .unwrap_or(1.0);
        out.push(Pair::spans(first, second, value));
    }
    Ok(out)
}

/// A value, where the spellings of nothing are a pair with no answer.
fn value(field: &str, line: usize) -> Result<f64, ReadError> {
    if field.is_empty()
        || field == "."
        || field.eq_ignore_ascii_case("na")
        || field.eq_ignore_ascii_case("nan")
    {
        return Ok(f64::NAN);
    }
    number(field, "the value", line)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    /// PLINK names both places of a pair, and a lead is named on every row.
    #[test]
    fn plink_names_the_places_of_its_pairs() {
        let ld = " CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n\
                   1 754400 rs42 1 655400 rs7 0.05\n1 754400 rs42 1 657600 rs8 0.05\n";
        let named = names(ld);
        assert_eq!(
            named,
            vec![
                (655_399, "rs7".to_string()),
                (657_599, "rs8".to_string()),
                (754_399, "rs42".to_string())
            ]
        );
        assert!(names("pos1 pos2 score\n10 20 1\n").is_empty());
    }

    #[test]
    fn plink_writes_linkage_with_its_own_names_and_counts_from_one() {
        let ld =
            " CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2 \n\
                  \x20    1         1001         rs1      1         1401         rs2     0.923 \n\
                  \x20    1         1001         rs1      2         1401         rs9     0.5 \n";
        let (read, named) = pairs(ld, &region("1:1-5,000")).unwrap();
        assert_eq!(read, vec![Pair::new(1000, 1400, 0.923)]);
        assert!(is_correlation(&named.unwrap()));
    }

    #[test]
    fn bedpe_is_two_stretches_as_bed_writes_them() {
        let contacts = "chr2\t0\t10000\tchr2\t20000\t30000\t57\n\
                        chr2\t0\t10000\tchr3\t0\t10000\t9\n\
                        chr2\t10000\t20000\tchr2\t10000\t20000\tloop\t12.5\n";
        let (read, _) = pairs(contacts, &region("chr2:1-40,000")).unwrap();
        assert_eq!(
            read,
            vec![
                Pair::spans((0, 10_000), (20_000, 30_000), 57.0),
                Pair::spans((10_000, 20_000), (10_000, 20_000), 12.5),
            ]
        );
    }

    #[test]
    fn a_table_of_your_own_is_read_by_its_header_or_as_three_columns() {
        let table = "site_a,site_b,score\n10,250,0.8\n10,40,NA\n";
        let (read, named) = pairs(table, &region("gene:1-300")).unwrap();
        assert_eq!(named.as_deref(), Some("score"));
        assert_eq!(read[0], Pair::new(9, 249, 0.8));
        assert!(read[1].value.is_nan(), "NA is no answer, not zero");
        let (bare, _) = pairs("10\t250\t0.8\n", &region("gene:1-300")).unwrap();
        assert_eq!(bare, vec![Pair::new(9, 249, 0.8)]);
        let error = pairs("left\tright\tscore\n1\t2\t3\n", &region("gene:1-300")).unwrap_err();
        assert!(error.reason.contains("pos1 and pos2"), "{}", error.reason);
    }
}
