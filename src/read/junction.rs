//! Splice junctions, as an aligner counts them.
//!
//! ```text
//! chr1  14830  14969  2  2  1  14  3  40
//! ```
//!
//! Nine columns, tab separated, which is what STAR writes to `SJ.out.tab`:
//! sequence, first base of the intron, last base of the intron, strand, intron
//! motif, whether an annotation held it, uniquely mapping reads across it,
//! multi-mapping reads across it, and the longest spliced overhang.
//!
//! # A third coordinate convention, and it is neither of its neighbours
//!
//! The two coordinates are 1-based and inclusive **on the intron**, not on the
//! exons either side. So the first base of the intron is column two minus one
//! counted from nought, and the intron made half-open ends at column three
//! unchanged. That is minus one on the start and the identity on the end, which
//! is what GFF3 does; a reader that copied [`point`](super::point), where both
//! come down by one, would draw every arc a base short at its right foot.
//!
//! # What is not read, and why
//!
//! Multi-mapping reads are kept and are never added to the unique ones. A read
//! that mapped in four places is one read and four pieces of evidence, and
//! adding it to a count the figure draws a thickness from would make a repeat
//! look like an expressed isoform.
//!
//! Column six is nought for a junction no annotation held. That is read as
//! stated novelty, and a file whose sixth column is nought everywhere is a file
//! read against no annotation, which cannot be told apart from a file whose
//! junctions are all novel. Both come back as `Some(false)`, and the caller who
//! knows which one they have is the one who can say.

use crate::{Junction, Motif, Region, Strand};

use super::{columns, lines, number, ReadError};

/// The junctions a file holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Junctions {
    /// The junctions, in the order the file listed them.
    pub junctions: Vec<Junction>,
    /// Rows in the file, before any filter.
    pub records: usize,
    /// Rows naming another sequence.
    pub other_sequence: usize,
    /// Rows on this sequence that touch no base of the window.
    pub off_region: usize,
    /// Rows no uniquely mapping read crossed.
    ///
    /// Kept in [`Junctions::junctions`] and left for the track to hold back, so
    /// that the figure can print how many it did not draw. The count is here
    /// too, for a caller that wants the number without drawing anything.
    pub uncrossed: usize,
}

/// Reads splice junctions out of STAR's `SJ.out.tab`.
///
/// The intron's coordinates are 1-based and inclusive, so one comes off the
/// start and the end passes through unchanged.
///
/// # Errors
///
/// Returns the first row without nine columns, whose coordinates are not
/// numbers, whose start is nought in a file counting from one, or whose end is
/// before its start.
pub fn junctions(text: &str, region: &Region) -> Result<Junctions, ReadError> {
    let mut found = Junctions::default();

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 9 {
            return Err(ReadError::at(
                at,
                format!(
                    "an SJ.out.tab row has nine columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        found.records += 1;
        if cols[0] != region.seq() {
            found.other_sequence += 1;
            continue;
        }

        let first: u64 = number(cols[1], "intron start", at)?;
        let last: u64 = number(cols[2], "intron end", at)?;
        if first == 0 {
            return Err(ReadError::at(
                at,
                "SJ.out.tab counts from 1, so 0 is not an intron start",
            ));
        }
        if last < first {
            return Err(ReadError::at(at, "the intron ends before it starts"));
        }
        // Inclusive on the intron: the start comes down by one and the end,
        // made half-open, is the number as written.
        let start = first - 1;
        let end = last;
        if end <= region.start() || start >= region.end() {
            found.off_region += 1;
            continue;
        }

        let unique: u32 = number(cols[6], "uniquely mapping reads", at)?;
        if unique == 0 {
            // Counted here and kept, rather than dropped here and counted
            // twice. Deciding not to draw a junction nobody crossed belongs to
            // the track, which is also what prints how many it held back: a
            // reader that filtered first would take that number off the figure.
            found.uncrossed += 1;
        }
        let multi: u32 = number(cols[7], "multi-mapping reads", at)?;
        let overhang: u32 = number(cols[8], "overhang", at)?;

        let mut junction = Junction::new(start, end, unique)
            .multi(multi)
            .overhang(overhang)
            .strand(match cols[3] {
                "1" => Strand::Forward,
                "2" => Strand::Reverse,
                _ => Strand::Unknown,
            });
        if let Some(motif) = motif(cols[4]) {
            junction = junction.motif(motif);
        }
        if let Some(annotated) = annotated(cols[5]) {
            junction = junction.annotated(annotated);
        }
        found.junctions.push(junction);
    }

    Ok(found)
}

/// STAR's six motif codes, folded to four.
///
/// The pairs differ only in which strand the intron is on, and the strand is
/// its own column, so keeping six would be the same fact written twice.
fn motif(code: &str) -> Option<Motif> {
    Some(match code {
        "0" => Motif::Noncanonical,
        "1" | "2" => Motif::GtAg,
        "3" | "4" => Motif::GcAg,
        "5" | "6" => Motif::AtAc,
        _ => return None,
    })
}

/// Whether an annotation held the junction, or nothing where the column does
/// not say.
fn annotated(code: &str) -> Option<bool> {
    match code {
        "0" => Some(false),
        "1" => Some(true),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> Region {
        Region::new("chr1", 0, 100_000).unwrap()
    }

    #[test]
    fn the_intron_is_inclusive_so_only_the_start_comes_down_by_one() {
        // A reader that copied point.rs would take one off both and draw every
        // arc a base short at its right foot.
        let found =
            junctions("chr1\t14830\t14969\t2\t2\t1\t14\t3\t40\n", &region()).expect("junctions");
        let junction = &found.junctions[0];
        assert_eq!(junction.start(), 14_829);
        assert_eq!(junction.end(), 14_969);
        assert_eq!(junction.len(), 140);
    }

    #[test]
    fn a_junction_no_unique_read_crossed_is_counted_and_handed_on() {
        // Handed on rather than dropped, because the track is what prints how
        // many it held back, and a reader that filtered first would take that
        // number off the figure.
        let found = junctions(
            "chr1\t100\t200\t1\t1\t1\t0\t9\t30\nchr1\t300\t400\t1\t1\t1\t5\t0\t30\n",
            &region(),
        )
        .expect("junctions");
        assert_eq!(found.records, 2);
        assert_eq!(found.uncrossed, 1);
        assert_eq!(found.junctions.len(), 2);
        assert!(!found.junctions[0].is_observed());
        assert!(found.junctions[1].is_observed());
    }

    #[test]
    fn multi_mapping_reads_are_kept_apart_from_the_unique_ones() {
        // A read that mapped in four places is one read and four pieces of
        // evidence. Added together, a repeat looks like an expressed isoform.
        let found =
            junctions("chr1\t100\t200\t1\t1\t1\t7\t900\t30\n", &region()).expect("junctions");
        assert_eq!(found.junctions[0].reads(), 7);
        assert_eq!(found.junctions[0].multi_reads(), Some(900));
    }

    #[test]
    fn the_six_motif_codes_fold_to_four_and_the_strand_keeps_its_own_column() {
        for (code, expected) in [
            ("0", Motif::Noncanonical),
            ("1", Motif::GtAg),
            ("2", Motif::GtAg),
            ("3", Motif::GcAg),
            ("5", Motif::AtAc),
        ] {
            let text = format!("chr1\t100\t200\t2\t{code}\t1\t5\t0\t30\n");
            let found = junctions(&text, &region()).expect("junctions");
            assert_eq!(found.junctions[0].intron_motif(), Some(expected), "{code}");
            assert_eq!(found.junctions[0].on_strand(), Strand::Reverse);
        }
    }

    #[test]
    fn novelty_is_read_as_stated_and_never_guessed() {
        let novel = junctions("chr1\t100\t200\t1\t1\t0\t5\t0\t30\n", &region()).expect("j");
        assert_eq!(novel.junctions[0].in_annotation(), Some(false));
        let known = junctions("chr1\t100\t200\t1\t1\t1\t5\t0\t30\n", &region()).expect("j");
        assert_eq!(known.junctions[0].in_annotation(), Some(true));
        let quiet = junctions("chr1\t100\t200\t1\t1\t.\t5\t0\t30\n", &region()).expect("j");
        assert_eq!(quiet.junctions[0].in_annotation(), None);
    }

    #[test]
    fn rows_elsewhere_are_counted_rather_than_read() {
        let found = junctions(
            "chr2\t100\t200\t1\t1\t1\t5\t0\t30\n\
             chr1\t500000\t500100\t1\t1\t1\t5\t0\t30\n\
             chr1\t100\t200\t1\t1\t1\t5\t0\t30\n",
            &region(),
        )
        .expect("junctions");
        assert_eq!(found.other_sequence, 1);
        assert_eq!(found.off_region, 1);
        assert_eq!(found.junctions.len(), 1);
    }

    #[test]
    fn a_row_that_is_not_an_sj_row_is_refused_by_line() {
        let error = junctions("chr1\t100\t200\n", &region()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("nine columns"), "{error}");

        let error = junctions("chr1\t0\t200\t1\t1\t1\t5\t0\t30\n", &region()).unwrap_err();
        assert!(error.reason.contains("counts from 1"), "{error}");

        let error = junctions("chr1\t500\t200\t1\t1\t1\t5\t0\t30\n", &region()).unwrap_err();
        assert!(error.reason.contains("ends before it starts"), "{error}");
    }
}
