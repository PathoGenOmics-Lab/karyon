//! Modified bases, as `modkit pileup` counts them.
//!
//! bedMethyl is nine BED-like columns and nine of counts, one row per position
//! per strand per modification:
//!
//! ```text
//! 1  chrom                  10  valid coverage
//! 2  start, 0-based         11  percent modified
//! 3  end, exclusive         12  count modified
//! 4  modification code      13  count canonical
//! 5  score                  14  count other modification
//! 6  strand                 15  count delete
//! 7  thick start            16  count fail
//! 8  thick end              17  count diff
//! 9  colour                 18  count nocall
//! ```
//!
//! Columns two and three are 0-based and half-open, which is this crate's own
//! convention, so nothing is added or taken off. It is worth saying because the
//! two readers next door do convert, and a methylation call moved one base left
//! lands on the other strand's partner in a CpG, where it looks entirely
//! reasonable.
//!
//! # One file is several tracks
//!
//! Column four is which modification was counted, and a dual-mode run writes
//! `m` and `h` rows at the same cytosine. Stacked on one axis those become two
//! marks at one position with nothing saying which is which, and
//! [`MethylationTrack::hemimethylated`](crate::MethylationTrack::hemimethylated)
//! pairs by position, so it would call a symmetrically modified CpG a strand
//! difference that no pileup reported.
//!
//! So a read names its modification. [`codes`] says what a file holds, and
//! [`sites`] takes the one to draw. The code is compared on its first field,
//! since `modkit --motif` writes `m,CG,0` where a plain run writes `m`.
//!
//! # Nought reads is not nought per cent
//!
//! `modkit` writes a row for a position it could not call: valid coverage nought,
//! every count nought, and column eleven `0.00`. Passed through that is a mark
//! on the baseline whose tooltip says `0% modified`, which is a measurement, and
//! the position was not measured. Those rows are skipped and counted.
//!
//! The fraction is column twelve over column ten rather than column eleven,
//! which is the same quotient rendered to two decimals and is the one field a
//! pileup can write the word `NaN` into. `"NaN".parse::<f64>()` succeeds, and
//! [`MethylSite`]'s only guard is a `clamp`, which propagates one, so the mark
//! disappears from the figure and leaves a tooltip stating a hard number.

use std::collections::BTreeMap;

use crate::{MethylSite, Region, Strand};

use super::{columns, lines, number, ReadError};

/// The sites a bedMethyl holds for one modification, and what it held besides.
#[derive(Debug, Clone, PartialEq)]
pub struct Calls {
    /// The sites, in the order the file listed them.
    pub sites: Vec<MethylSite>,
    /// Rows in the file, before any filter.
    pub records: usize,
    /// Rows counting a different modification.
    pub other_code: usize,
    /// Rows naming another sequence.
    pub passed_over: usize,
    /// Rows on this sequence that fall outside the window.
    pub off_region: usize,
    /// Rows where nothing passed the caller's threshold.
    ///
    /// A position with no valid coverage was not measured, and the file says so
    /// by writing nought in every count. It is not a position measured at
    /// nought per cent, and the two are one mark apart on the page.
    pub no_coverage: usize,
}

/// Every modification code in the file, with how many rows each has.
///
/// Ordered, so the same file answers the same way twice, and the question a
/// caller has to settle before drawing anything: one file is one track only
/// when it counted one modification.
pub fn codes(text: &str) -> Result<BTreeMap<String, usize>, ReadError> {
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for (at, line) in lines(text) {
        let cols = fields(line, at)?;
        *seen.entry(code(&cols[3]).to_string()).or_insert(0) += 1;
    }
    Ok(seen)
}

/// The calls for one modification, inside one window.
///
/// Coordinates pass straight through: bedMethyl counts from nought and its end
/// is exclusive, which is what the rest of the crate counts in.
///
/// Rows on another sequence than `region.seq()` are skipped, and so are rows
/// counting a different modification, rows outside the window, and rows with no
/// valid coverage. All four are counted rather than dropped in silence.
///
/// # Errors
///
/// Returns the first row that is not bedMethyl, whose strand is neither `+` nor
/// `-`, or whose counts do not make a fraction. A file that parses and holds
/// nothing for this modification in this window is not an error here:
/// [`Calls`] says which of the four reasons it was, and the caller decides.
pub fn sites(text: &str, region: &Region, wanted: &str) -> Result<Calls, ReadError> {
    let mut found = Calls {
        sites: Vec::new(),
        records: 0,
        other_code: 0,
        passed_over: 0,
        off_region: 0,
        no_coverage: 0,
    };

    for (at, line) in lines(text) {
        let cols = fields(line, at)?;
        found.records += 1;

        if cols[0] != region.seq() {
            found.passed_over += 1;
            continue;
        }
        if code(&cols[3]) != wanted {
            found.other_code += 1;
            continue;
        }

        let start: u64 = number(&cols[1], "start", at)?;
        if start < region.start() || start >= region.end() {
            found.off_region += 1;
            continue;
        }

        // `.` is what a strand-combined pileup writes, and it is not a strand.
        // Read as one it puts every call in the forward lane, in the forward
        // colour, with a tooltip that says forward, and the band then shows one
        // strand's worth of data claiming to be both.
        let strand = match cols[5].as_str() {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            other => {
                return Err(ReadError::at(
                    at,
                    format!(
                        "the strand column is + or -, this one is {other:?}; a pileup with its \
                         strands combined has no strand to draw"
                    ),
                ))
            }
        };

        let coverage: u32 = number(&cols[9], "valid coverage", at)?;
        if coverage == 0 {
            found.no_coverage += 1;
            continue;
        }
        let modified: u32 = number(&cols[11], "count modified", at)?;
        if modified > coverage {
            return Err(ReadError::at(
                at,
                format!("{modified} reads modified out of {coverage} valid"),
            ));
        }

        found.sites.push(MethylSite::new(
            start,
            strand,
            f64::from(modified) / f64::from(coverage),
            coverage,
        ));
    }

    Ok(found)
}

/// The columns of one row, however the file spaced them.
///
/// bedMethyl is written with tabs through column ten and spaces after it by at
/// least one of the tools that produce it, and [`columns`] splits on tabs the
/// moment a line holds one, so such a row arrives as ten fields with the last
/// nine glued into the tenth. Splitting the tail again is the difference
/// between reading the counts and reading a length.
fn fields(line: &str, at: usize) -> Result<Vec<String>, ReadError> {
    let mut cols: Vec<String> = columns(line).into_iter().map(str::to_string).collect();
    if cols.len() < 18 {
        if let Some(tail) = cols.pop() {
            cols.extend(tail.split_whitespace().map(str::to_string));
        }
    }
    if cols.len() < 18 {
        return Err(ReadError::at(
            at,
            format!(
                "a bedMethyl row has 18 columns, this one has {}; `modkit pileup` writes \
                 nine counts after the nine BED ones",
                cols.len()
            ),
        ));
    }
    Ok(cols)
}

/// The modification a row counted, without the motif a `--motif` run appends.
fn code(field: &str) -> &str {
    field.split(',').next().unwrap_or(field)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two modifications at one cytosine, one uncallable position, and a row on
    /// another sequence, which is every shape this reader has to tell apart.
    const BED: &str = "\
chr1\t99\t100\tm\t40\t+\t99\t100\t255,0,0\t40\t95.00\t38\t2\t0\t0\t3\t0\t1
chr1\t99\t100\th\t40\t+\t99\t100\t255,0,0\t40\t5.00\t2\t38\t0\t0\t3\t0\t1
chr1\t150\t151\tm\t0\t-\t150\t151\t255,0,0\t0\t0.00\t0\t0\t0\t0\t12\t0\t4
chr1\t200\t201\tm\t20\t-\t200\t201\t255,0,0\t20\t50.00\t10\t10\t0\t0\t0\t0\t0
chr2\t99\t100\tm\t40\t+\t99\t100\t255,0,0\t40\t95.00\t38\t2\t0\t0\t3\t0\t1
";

    fn window() -> Region {
        Region::new("chr1", 0, 1_000).unwrap()
    }

    #[test]
    fn a_position_passes_straight_through_because_bedmethyl_counts_from_nought() {
        // The reader next door takes one off a GFF3 start. Doing it here moves
        // every call onto its partner base in a CpG, where it looks right.
        let found = sites(BED, &window(), "m").unwrap();
        assert_eq!(found.sites[0].pos, 99);
        assert_eq!(found.sites[1].pos, 200);
    }

    #[test]
    fn the_fraction_is_the_counts_and_not_the_rendered_percentage() {
        let found = sites(BED, &window(), "m").unwrap();
        // 38 of 40, which is 0.95 exactly, where column eleven says 95.00.
        assert_eq!(found.sites[0].fraction, 0.95);
        assert_eq!(found.sites[0].coverage, 40);
        assert_eq!(found.sites[1].fraction, 0.5);
    }

    #[test]
    fn a_file_says_which_modifications_it_holds_before_anything_is_drawn() {
        let held = codes(BED).unwrap();
        assert_eq!(held.get("m"), Some(&4));
        assert_eq!(held.get("h"), Some(&1));

        // And a read takes one of them, counting the rest rather than stacking
        // two modifications of one cytosine on one axis.
        let found = sites(BED, &window(), "m").unwrap();
        assert_eq!(found.sites.len(), 2);
        assert_eq!(found.other_code, 1);
    }

    #[test]
    fn a_motif_run_writes_the_code_with_its_motif_and_still_matches() {
        let text =
            "chr1\t99\t100\tm,CG,0\t40\t+\t99\t100\t0,0,0\t40\t95.00\t38\t2\t0\t0\t0\t0\t0\n";
        assert_eq!(codes(text).unwrap().get("m"), Some(&1));
        assert_eq!(sites(text, &window(), "m").unwrap().sites.len(), 1);
    }

    #[test]
    fn nought_valid_reads_is_not_nought_per_cent() {
        // The row at 150 has no valid coverage and column eleven says 0.00.
        // Read as a measurement it is a mark on the baseline saying this
        // cytosine is unmethylated, and the pileup said it could not tell.
        let found = sites(BED, &window(), "m").unwrap();
        assert_eq!(found.no_coverage, 1);
        assert!(
            found.sites.iter().all(|site| site.pos != 150),
            "a position nobody could call was drawn at nought per cent"
        );
        assert!(found.sites.iter().all(|site| site.coverage > 0));
    }

    #[test]
    fn a_fraction_that_is_not_a_number_cannot_be_built_from_counts() {
        // Column eleven is the one field a pileup writes the word into, and
        // `"NaN".parse::<f64>()` succeeds. Taking the quotient of two integer
        // counts instead makes the value unrepresentable rather than guarded.
        let text = "chr1\t99\t100\tm\t40\t+\t99\t100\t0,0,0\t40\tNaN\t38\t2\t0\t0\t0\t0\t0\n";
        let found = sites(text, &window(), "m").unwrap();
        assert_eq!(found.sites[0].fraction, 0.95);
        assert!(found.sites[0].fraction.is_finite());
    }

    #[test]
    fn more_modified_reads_than_valid_ones_stops_the_read() {
        let text = "chr1\t99\t100\tm\t40\t+\t99\t100\t0,0,0\t10\t95.00\t38\t2\t0\t0\t0\t0\t0\n";
        let error = sites(text, &window(), "m").unwrap_err();
        assert!(error.reason.contains("out of 10"), "{error}");
    }

    #[test]
    fn a_strand_combined_pileup_is_refused_rather_than_drawn_as_one_strand() {
        // `.` becomes Strand::Unknown, which is not reverse, so every call
        // lands in the forward lane in the forward colour with a tooltip that
        // says forward, and half the band is empty for a reason nobody stated.
        let text = "chr1\t99\t100\tm\t40\t.\t99\t100\t0,0,0\t40\t95.00\t38\t2\t0\t0\t0\t0\t0\n";
        let error = sites(text, &window(), "m").unwrap_err();
        assert!(error.reason.contains("strands combined"), "{error}");
    }

    #[test]
    fn rows_elsewhere_are_counted_rather_than_dropped_in_silence() {
        let found = sites(BED, &window(), "m").unwrap();
        assert_eq!(found.passed_over, 1, "another sequence went unmentioned");
        assert_eq!(found.records, 5);

        let narrow = sites(BED, &Region::new("chr1", 300, 400).unwrap(), "m").unwrap();
        assert!(narrow.sites.is_empty());
        // All three m rows on chr1, since the window filter comes before the
        // coverage one and the uncallable row was never going to be drawn.
        assert_eq!(narrow.off_region, 3);
        assert_eq!(narrow.records, 5, "the file did hold rows");
    }

    #[test]
    fn a_row_spaced_rather_than_tabbed_still_reads() {
        // Tabs through column ten and spaces after is what at least one tool
        // writes, and the column splitter takes a tab as proof of a tab file,
        // so the last nine fields arrive glued into the tenth.
        let mixed = "chr1\t99\t100\tm\t40\t+\t99\t100\t255,0,0\t40 95.00 38 2 0 0 3 0 1\n";
        let found = sites(mixed, &window(), "m").unwrap();
        assert_eq!(found.sites.len(), 1);
        assert_eq!(found.sites[0].fraction, 0.95);

        let spaced = "chr1 99 100 m 40 + 99 100 255,0,0 40 95.00 38 2 0 0 3 0 1\n";
        assert_eq!(sites(spaced, &window(), "m").unwrap().sites.len(), 1);
    }

    #[test]
    fn a_line_that_is_not_bedmethyl_stops_the_read_and_says_which() {
        let short = "chr1\t99\t100\tm\t40\t+\n";
        let error = sites(short, &window(), "m").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("18 columns"), "{error}");
    }
}
