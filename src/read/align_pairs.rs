//! PAF: one alignment per line, carrying two coordinate systems at once.
//!
//! `minimap2` writes this by default, and so do `miniasm`, `winnowmap`,
//! `wfmash` and `paftools.js` converting a nucmer delta. Twelve mandatory
//! columns, then optional `tag:type:value` fields this reader does not need:
//!
//! ```text
//! 1  query name          7  target length
//! 2  query length        8  target start, 0-based
//! 3  query start         9  target end, half-open
//! 4  query end          10  residue matches
//! 5  strand, + or -     11  alignment block length
//! 6  target name        12  mapping quality
//! ```
//!
//! Both coordinate pairs are already 0-based and half-open, which is the
//! crate's own convention, so nothing is added or taken off here. That makes
//! PAF the one format in this directory that needs no conversion, and it is
//! worth saying out loud because every other reader in here does.
//!
//! # A block belongs to a pair of sequences, and the track only knows one
//!
//! [`AlignmentBlock`] holds two spans and no names. Two
//! rows of a whole-genome PAF can therefore look identical to a track while
//! describing alignments against different chromosomes, and stacking them on
//! one vertical axis would draw a comparison nobody made.
//!
//! So a read names its target. [`blocks`] takes the query the figure is on and
//! the target it is being compared with, keeps the rows that match both, and
//! returns how many rows it passed over. That count is the point: a caller that
//! draws 40 blocks from a file holding 4,000 should be able to say so, the same
//! way [`Map`](crate::Map) says how many locations it could not place.
//!
//! [`targets`] answers the question that comes first, which is what a file even
//! contains, so a caller can choose rather than guess.

use std::collections::BTreeMap;

use crate::track::AlignmentBlock;

use super::{columns, lines, number, ReadError};

/// Every query and target pair in the file, with how many rows each has.
///
/// Ordered, so the same file always answers in the same order, and a caller
/// picking "the first one" picks the same one twice.
pub fn pairs(text: &str) -> Result<BTreeMap<(String, String), usize>, ReadError> {
    let mut seen: BTreeMap<(String, String), usize> = BTreeMap::new();
    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 12 {
            return Err(ReadError::at(
                at,
                format!(
                    "a PAF line has at least 12 columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        *seen
            .entry((cols[0].to_string(), cols[5].to_string()))
            .or_insert(0) += 1;
    }
    Ok(seen)
}

/// Target names in the file, for one query, most-aligned first.
pub fn targets(text: &str, query: &str) -> Result<Vec<(String, usize)>, ReadError> {
    let mut out: Vec<(String, usize)> = pairs(text)?
        .into_iter()
        .filter(|((q, _), _)| q == query)
        .map(|((_, t), n)| (t, n))
        .collect();
    // Most rows first, and by name where they tie, so it is a stable answer.
    out.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(out)
}

/// What one PAF says about one pair of sequences.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignments {
    /// The blocks between the two, in the order the file listed them.
    pub blocks: Vec<AlignmentBlock>,
    /// The target's whole length, from column seven.
    ///
    /// The track wants it and the alignments cannot supply it: blocks cover
    /// the part that aligned, so the longest of them says how far the
    /// alignment reached and never how long the sequence is. `None` when no
    /// row named this pair, since then nothing said.
    pub target_length: Option<u64>,
    /// Rows about some other pair of sequences.
    ///
    /// Not an error. A PAF of a whole assembly against a whole reference is
    /// the ordinary case, and one figure shows one pair, so passing over the
    /// rest is right. Saying nothing about it would not be.
    pub passed_over: usize,
}

/// The alignments between one query and one target.
pub fn blocks(text: &str, query: &str, target: &str) -> Result<Alignments, ReadError> {
    let mut blocks = Vec::new();
    let mut passed_over = 0usize;
    let mut target_length = None;

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 12 {
            return Err(ReadError::at(
                at,
                format!(
                    "a PAF line has at least 12 columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        if cols[0] != query || cols[5] != target {
            passed_over += 1;
            continue;
        }

        let query_start: u64 = number(cols[2], "query start", at)?;
        let query_end: u64 = number(cols[3], "query end", at)?;
        let target_start: u64 = number(cols[7], "target start", at)?;
        let target_end: u64 = number(cols[8], "target end", at)?;

        // Column five is the whole of what says the target runs backwards, and
        // a PAF writes the target's coordinates ascending either way, so the
        // flag has to be carried rather than inferred from the numbers.
        let reversed = match cols[4] {
            "+" => false,
            "-" => true,
            other => {
                return Err(ReadError::at(
                    at,
                    format!("the strand column is + or -, this one is {other:?}"),
                ))
            }
        };

        let mut block = AlignmentBlock::new(query_start, query_end, target_start, target_end)
            .reversed(reversed);

        // Identity is matches over block length. A block length of nought is
        // not a block, and dividing by it would put a number on the page that
        // no aligner reported, so the block keeps no identity at all instead.
        let matches: u64 = number(cols[9], "residue matches", at)?;
        let span: u64 = number(cols[10], "alignment block length", at)?;
        if span > 0 {
            block = block.identity(matches as f64 / span as f64);
        }
        // Column seven is the same number on every row of a pair, so the first
        // one is as good as any, and a file that disagrees with itself about
        // it is a file two runs were concatenated into.
        let stated: u64 = number(cols[6], "target length", at)?;
        match target_length {
            None => target_length = Some(stated),
            Some(first) if first != stated => {
                return Err(ReadError::at(
                    at,
                    format!("{target} is {stated} long here and {first} long earlier"),
                ))
            }
            Some(_) => {}
        }
        blocks.push(block);
    }

    Ok(Alignments {
        blocks,
        target_length,
        passed_over,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two queries against two targets, one of them on the reverse strand.
    const PAF: &str = "\
q1\t1000\t0\t100\t+\tt1\t2000\t500\t600\t95\t100\t60
q1\t1000\t200\t300\t-\tt1\t2000\t900\t1000\t80\t100\t60
q1\t1000\t400\t500\t+\tt2\t3000\t10\t110\t99\t100\t60
q2\t800\t0\t50\t+\tt1\t2000\t0\t50\t50\t50\t60
";

    #[test]
    fn both_coordinate_pairs_come_through_untouched() {
        // PAF is already 0-based and half-open, so this reader is the one in
        // the directory that must NOT move a coordinate.
        let Alignments { blocks, .. } = blocks(PAF, "q1", "t1").unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(
            (blocks[0].query_start, blocks[0].query_end),
            (0, 100),
            "a query span moved"
        );
        assert_eq!(
            (blocks[0].target_start, blocks[0].target_end),
            (500, 600),
            "a target span moved"
        );
    }

    #[test]
    fn the_strand_column_is_the_only_thing_that_says_backwards() {
        let Alignments { blocks, .. } = blocks(PAF, "q1", "t1").unwrap();
        assert!(!blocks[0].reversed);
        // The second row's target coordinates ascend exactly as the first's do,
        // so nothing but column five distinguishes them.
        assert!(blocks[1].target_start < blocks[1].target_end);
        assert!(blocks[1].reversed, "a reverse block was read as forward");
    }

    #[test]
    fn rows_against_another_sequence_are_counted_rather_than_drawn() {
        // The whole reason this reader takes two names. Without the filter
        // these four rows become four blocks on one axis, describing a
        // comparison that was never made.
        let Alignments {
            blocks,
            passed_over,
            ..
        } = blocks(PAF, "q1", "t1").unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(passed_over, 2, "rows for another pair went missing quietly");
    }

    #[test]
    fn identity_is_matches_over_the_block_and_nothing_when_there_is_no_block() {
        let Alignments { blocks: found, .. } = blocks(PAF, "q1", "t1").unwrap();
        assert_eq!(found[0].identity, Some(0.95));

        let zero = "q1\t10\t0\t5\t+\tt1\t10\t0\t5\t0\t0\t60\n";
        let Alignments {
            blocks: no_span, ..
        } = blocks(zero, "q1", "t1").unwrap();
        assert_eq!(
            no_span[0].identity, None,
            "a block of no length was given an identity anyway"
        );
    }

    #[test]
    fn a_file_says_what_it_holds_before_anything_is_drawn() {
        let found = targets(PAF, "q1").unwrap();
        assert_eq!(found, vec![("t1".to_string(), 2), ("t2".to_string(), 1)]);
        assert_eq!(pairs(PAF).unwrap().len(), 3);
    }

    #[test]
    fn a_line_that_is_not_paf_stops_the_read_and_says_which() {
        let short = "q1\t1000\t0\t100\t+\tt1\n";
        let error = blocks(short, "q1", "t1").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("at least 12 columns"), "{error}");

        let strand = "q1\t1000\t0\t100\t?\tt1\t2000\t0\t100\t99\t100\t60\n";
        let error = blocks(strand, "q1", "t1").unwrap_err();
        assert!(error.reason.contains("+ or -"), "{error}");
    }
}
