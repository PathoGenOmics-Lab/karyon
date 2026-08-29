//! Methylation one molecule at a time, from a Bismark extractor file.
//!
//! `bismark_methylation_extractor` writes one row per cytosine per read:
//!
//! ```text
//! 1  read name        the SAM query name; both mates of a pair carry it
//! 2  + or -           the call again, and not the strand
//! 3  sequence
//! 4  position, 1-based
//! 5  call             Z z X x H h U u
//! ```
//!
//! Position counts from one, so one is taken off it. A call moved a base lands
//! on the other cytosine of the same CpG, where it looks entirely reasonable.
//!
//! The letter is both the context and the answer: `Z` is a methylated CpG and
//! `z` an unmethylated one, `X`/`x` are CHG and `H`/`h` are CHH. Column two is
//! the case of that letter written out again, so it says nothing the letter did
//! not, and in particular it is not the strand.
//!
//! # A grid is built by position and never by order
//!
//! [`BisulfiteTrack`](crate::BisulfiteTrack) is a matrix: one shared list of
//! sites, and per molecule one call for each of them, in that order.
//! [`Molecule::calls`](crate::Molecule::calls) is indexed by the site's place in
//! the list and nothing checks that a caller got it right, so a vector built by
//! pushing calls in the order the file happened to list them puts every call
//! after the first gap one column to the left. That is a methylation pattern
//! that never existed, drawn as cleanly as one that did.
//!
//! So every row here is built at the full width of the site list and written
//! into by position. A site a molecule never covered stays absent, which the
//! track draws as nothing at all rather than as the ring it draws for a
//! cytosine that was measured and found unmethylated.
//!
//! # The row is the fragment
//!
//! Both mates of a pair carry one query name, and they are two halves of one
//! molecule, so they are one row. Where they overlap and disagree about a
//! cytosine, neither call is kept: the molecule's state there is what the two
//! reads could not agree on, and either answer drawn is a coin toss shown as a
//! measurement.
//!
//! # The two cytosines of a CpG stay two sites
//!
//! A CpG is a C on each strand, one base apart, and this file has a row for
//! each. They are not folded together, because folding needs to know which
//! strand a call came from and this file does not say: column two is the call,
//! and the strand lives in the name of the file the extractor wrote, which a
//! reader handed a `&str` cannot see. Two columns a base apart are what was
//! measured, and pairing them would be a reconstruction offered as a reading.

use std::collections::{BTreeMap, BTreeSet};

use crate::{Molecule, Region};

use super::{columns, lines, number, ReadError};

/// The header line the extractor writes, unless it was told not to.
const HEADER: &str = "Bismark methylation extractor version";

/// The molecules a file holds for one context, and what it held besides.
#[derive(Debug, Clone, PartialEq)]
pub struct Calls {
    /// The sites every molecule is measured against, ascending.
    pub sites: Vec<u64>,
    /// One row per fragment, in the order the reads were first seen, with one
    /// call per site by construction.
    pub molecules: Vec<Molecule>,
    /// Rows in the file, before any filter.
    pub records: usize,
    /// Rows naming another sequence.
    pub passed_over: usize,
    /// Rows on this sequence that fall outside the window.
    pub off_region: usize,
    /// Rows calling a cytosine in another context.
    pub other_context: usize,
    /// Cells where the two mates of one fragment disagreed.
    ///
    /// Neither call is kept. A cytosine two reads of one molecule answered
    /// differently is not a cytosine that was measured, and drawing either
    /// answer shows a coin toss as a result.
    pub contradicted: usize,
}

/// Every cytosine context in the file, with how many rows each has.
///
/// Ordered, so the same file answers the same way twice. The question that
/// comes before drawing anything, since a file written with every context in it
/// is three grids and only one of them is a track.
pub fn contexts(text: &str) -> Result<BTreeMap<String, usize>, ReadError> {
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for (at, line) in body(text) {
        let cols = columns(line);
        check(&cols, at)?;
        *seen.entry(context(cols[4], at)?.to_string()).or_insert(0) += 1;
    }
    Ok(seen)
}

/// The molecules for one context, inside one window.
///
/// The position counts from one and comes back one lower. The site list is
/// every position any kept row named, ascending, and every molecule carries one
/// call for each of them.
///
/// Rows on another sequence than `region.seq()` are skipped, and so are rows in
/// another context and rows outside the window. All three are counted.
///
/// # Errors
///
/// Returns the first row that is not an extractor row, whose call letter is not
/// one of the eight, whose position is nought, or whose second column
/// contradicts the case of its own call letter. A file that parses and holds
/// nothing for this context in this window is not an error here: [`Calls`] says
/// which of the three reasons it was.
pub fn molecules(text: &str, region: &Region, wanted: &str) -> Result<Calls, ReadError> {
    let mut found = Calls {
        sites: Vec::new(),
        molecules: Vec::new(),
        records: 0,
        passed_over: 0,
        off_region: 0,
        other_context: 0,
        contradicted: 0,
    };

    // `None` in the inner map is a cell the two mates disagreed about, which is
    // not the same as a cell neither of them reached: that one has no entry.
    let mut per_read: BTreeMap<String, BTreeMap<u64, Option<bool>>> = BTreeMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut sites: BTreeSet<u64> = BTreeSet::new();

    for (at, line) in body(text) {
        let cols = columns(line);
        check(&cols, at)?;
        found.records += 1;

        if cols[2] != region.seq() {
            found.passed_over += 1;
            continue;
        }
        // Before the position joins the site list. A cytosine of another
        // context left in becomes a column every molecule of this one is
        // absent from, which draws as a site nobody could call.
        if context(cols[4], at)? != wanted {
            found.other_context += 1;
            continue;
        }

        let modified = modified(cols[4], at)?;
        // Column two is the case of the call written out again. Where the two
        // disagree the file is not the one the flag claimed it was, and the
        // cheapest thing that could be wrong with it is the column order.
        let stated = match cols[1] {
            "+" => true,
            "-" => false,
            other => {
                return Err(ReadError::at(
                    at,
                    format!("column two of an extractor row is + or -, this one is {other:?}"),
                ))
            }
        };
        if stated != modified {
            return Err(ReadError::at(
                at,
                format!(
                    "column two says {:?} and the call {:?} says the other, so the columns are \
                     not the ones this reader was told to expect",
                    cols[1], cols[4]
                ),
            ));
        }

        let stated: u64 = number(cols[3], "position", at)?;
        let Some(pos) = stated.checked_sub(1) else {
            return Err(ReadError::at(
                at,
                "the extractor counts from 1, so 0 is not a position",
            ));
        };
        // Before the site list too. A site outside the window is drawn past the
        // edge of the plot and clipped away, and still counted in the row's own
        // tooltip, so the figure says it measured what nobody can see.
        if pos < region.start() || pos >= region.end() {
            found.off_region += 1;
            continue;
        }
        sites.insert(pos);

        // Both mates of a pair carry the fragment's name, and a fragment is one
        // molecule.
        let read = cols[0]
            .strip_suffix("/1")
            .or_else(|| cols[0].strip_suffix("/2"))
            .unwrap_or(cols[0]);
        let calls = match per_read.get_mut(read) {
            Some(calls) => calls,
            None => {
                order.push(read.to_string());
                per_read.entry(read.to_string()).or_default()
            }
        };
        match calls.get(&pos) {
            None => {
                calls.insert(pos, Some(modified));
            }
            Some(Some(before)) if *before != modified => {
                found.contradicted += 1;
                calls.insert(pos, None);
            }
            // The same answer twice, or a cell already given up on.
            Some(_) => {}
        }
    }

    found.sites = sites.into_iter().collect();
    // Written by position into a row of the full width, which is the whole of
    // what keeps a call in its own column.
    for name in order {
        let calls = &per_read[&name];
        let mut row = vec![None; found.sites.len()];
        // The molecule's own calls, each finding its column, rather than every
        // column asking this molecule whether it was reached. A fragment
        // covers a few dozen sites and the window holds tens of thousands, so
        // asking the other way round was twenty thousand molecules times
        // thirty thousand sites of lookups to place a million calls, and 2.12
        // of the 2.38 seconds the figure took. Every position here was put
        // into the site list by the same line that put it here, so the search
        // finds it.
        for (site, call) in calls {
            if let Ok(column) = found.sites.binary_search(site) {
                row[column] = *call;
            }
        }
        found.molecules.push(Molecule::new(name, row));
    }

    Ok(found)
}

/// The lines that are rows, without the version line the extractor writes.
///
/// That line carries no `#` and no `@`, so nothing upstream drops it, and the
/// column splitter falls back to whitespace when there is no tab, so it comes
/// apart into exactly five fields like a row of data. It has to go by name.
fn body(text: &str) -> impl Iterator<Item = (usize, &str)> {
    lines(text).filter(|(_, line)| !line.trim_start().starts_with(HEADER))
}

/// Refuses a row that is not an extractor row.
fn check(cols: &[&str], at: usize) -> Result<(), ReadError> {
    // Five columns as standard and eight with `--yacht`, which adds where the
    // read started and ended. Neither of those is needed to place a call.
    if cols.len() != 5 && cols.len() != 8 {
        return Err(ReadError::at(
            at,
            format!(
                "a Bismark extractor row is 5 columns, or 8 with --yacht, and this one has {}",
                cols.len()
            ),
        ));
    }
    Ok(())
}

/// The cytosine context a call letter names.
fn context(call: &str, at: usize) -> Result<&'static str, ReadError> {
    Ok(match call {
        "Z" | "z" => "CpG",
        "X" | "x" => "CHG",
        "H" | "h" => "CHH",
        "U" | "u" => "Unknown",
        other => {
            return Err(ReadError::at(
                at,
                format!("a methylation call is one of ZzXxHhUu, this one is {other:?}"),
            ))
        }
    })
}

/// Whether a call letter says the cytosine was modified.
fn modified(call: &str, at: usize) -> Result<bool, ReadError> {
    // Through `context` first, so an unknown letter fails as an unknown letter
    // rather than as an unmethylated one by virtue of its case.
    context(call, at)?;
    Ok(call.chars().all(|letter| letter.is_ascii_uppercase()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two fragments over three CpGs, one CHH row, one row elsewhere, and the
    /// version line the extractor writes.
    const TEXT: &str = "\
Bismark methylation extractor version v0.25.1
r1/1\t+\tchr1\t101\tZ
r1/1\t-\tchr1\t201\tz
r1/2\t+\tchr1\t301\tZ
r2/1\t-\tchr1\t101\tz
r2/1\t+\tchr1\t301\tZ
r2/1\t-\tchr1\t150\th
r3/1\t+\tchr2\t101\tZ
";

    fn window() -> Region {
        Region::new("chr1", 0, 1_000).unwrap()
    }

    #[test]
    fn a_position_moves_back_one_because_the_extractor_counts_from_one() {
        // A call moved a base lands on the other cytosine of the same CpG,
        // which is a real site and the wrong one.
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        assert_eq!(found.sites, vec![100, 200, 300]);
    }

    #[test]
    fn a_call_lands_in_its_own_column_however_many_sites_a_read_missed() {
        // The whole reason this reader exists. r2 says nothing at site 200, so
        // a row built by pushing in file order would put its call for 300 into
        // the column for 200, and the figure would show a pattern nobody
        // measured, drawn as cleanly as one somebody did.
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        let r2 = found
            .molecules
            .iter()
            .find(|m| m.name == "r2")
            .expect("r2 is a molecule");
        assert_eq!(r2.calls.len(), found.sites.len());
        assert_eq!(r2.calls, vec![Some(false), None, Some(true)]);
    }

    #[test]
    fn a_site_a_molecule_did_not_cover_is_absent_and_not_unmethylated() {
        // The track draws nothing for absent and a ring for unmethylated, and
        // the two are different statements about the same cytosine.
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        let r2 = found.molecules.iter().find(|m| m.name == "r2").unwrap();
        assert_eq!(r2.call(1), None, "an uncovered site was given a call");
        assert_eq!(r2.call(0), Some(false), "a measured site lost its call");
    }

    #[test]
    fn both_mates_of_a_pair_are_one_molecule() {
        // They are two halves of one fragment. Two rows would draw one molecule
        // twice and halve the apparent coverage of each.
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        assert_eq!(found.molecules.len(), 2);
        let r1 = found.molecules.iter().find(|m| m.name == "r1").unwrap();
        assert_eq!(
            r1.calls,
            vec![Some(true), Some(false), Some(true)],
            "the second mate's call did not reach the fragment"
        );
    }

    #[test]
    fn two_mates_that_disagree_leave_the_cell_unmeasured() {
        // Either answer drawn is a coin toss shown as a measurement, and the
        // track has a spelling for a cell nobody could call.
        let text = "r1/1\t+\tchr1\t101\tZ\nr1/2\t-\tchr1\t101\tz\n";
        let found = molecules(text, &window(), "CpG").unwrap();
        assert_eq!(found.contradicted, 1);
        assert_eq!(found.molecules[0].calls, vec![None]);
        assert_eq!(found.sites, vec![100], "the site itself is still a site");
    }

    #[test]
    fn a_cytosine_of_another_context_is_not_a_column_of_this_grid() {
        // Left in, it becomes a site every molecule of this context is absent
        // from, which the figure draws as a cytosine nobody could call.
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        assert!(!found.sites.contains(&149), "a CHH site joined a CpG grid");
        assert_eq!(found.other_context, 1);

        // And the file says what it holds before anything is drawn.
        let held = contexts(TEXT).unwrap();
        assert_eq!(held.get("CpG"), Some(&6));
        assert_eq!(held.get("CHH"), Some(&1));
    }

    #[test]
    fn a_site_outside_the_window_is_counted_rather_than_drawn_past_the_edge() {
        // Kept, it is drawn beyond the plot and clipped away, and still counted
        // in the row's own tooltip, so the figure says it measured what nobody
        // can see.
        let narrow = molecules(TEXT, &Region::new("chr1", 0, 250).unwrap(), "CpG").unwrap();
        assert_eq!(narrow.sites, vec![100, 200]);
        assert_eq!(narrow.off_region, 2);
        assert!(narrow.molecules.iter().all(|m| m.calls.len() == 2));
    }

    #[test]
    fn rows_on_another_sequence_are_counted_rather_than_dropped_in_silence() {
        let found = molecules(TEXT, &window(), "CpG").unwrap();
        assert_eq!(found.passed_over, 1);
        assert_eq!(found.records, 7, "the version line was counted as a row");
    }

    #[test]
    fn the_version_line_is_not_a_row_of_data() {
        // It carries no comment marker, so nothing upstream drops it, and it
        // splits into five fields like a row.
        assert_eq!(
            columns("Bismark methylation extractor version v0.25.1").len(),
            5
        );
        assert!(molecules(TEXT, &window(), "CpG").is_ok());
        assert_eq!(contexts(TEXT).unwrap().values().sum::<usize>(), 7);
    }

    #[test]
    fn a_row_that_contradicts_itself_stops_the_read() {
        // Column two is the case of the call written out again, so the two
        // disagreeing means these are not the columns this reader expects, and
        // the cheapest thing that could be wrong is their order.
        let error = molecules("r1\t+\tchr1\t101\tz\n", &window(), "CpG").unwrap_err();
        assert!(error.reason.contains("the other"), "{error}");
    }

    #[test]
    fn a_line_that_is_not_an_extractor_row_stops_the_read_and_says_which() {
        let error = molecules("r1\t+\tchr1\n", &window(), "CpG").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("5 columns"), "{error}");

        let error = molecules("r1\t+\tchr1\t101\tQ\n", &window(), "CpG").unwrap_err();
        assert!(error.reason.contains("ZzXxHhUu"), "{error}");

        let error = molecules("r1\t+\tchr1\t0\tZ\n", &window(), "CpG").unwrap_err();
        assert!(error.reason.contains("counts from 1"), "{error}");
    }
}
