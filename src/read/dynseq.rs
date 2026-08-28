//! Per-base attribution scores, as a bedGraph and nothing else.
//!
//! ```text
//! chr1  1000  1001  0.42
//! chr1  1001  1002  -0.13
//! chr1  1002  1003  0.05
//! ```
//!
//! # One shape, on purpose
//!
//! [`signal::dense`](super::signal::dense) reads three shapes and tells them
//! apart by how many columns a row has. Two of the three mean nothing here. The
//! output of `samtools depth` is a read depth, which is never negative and is
//! not a contribution; a bare column of values carries no coordinates, so a
//! model scored over one stretch of a chromosome would be laid over another.
//! Reading either of them as an attribution would draw a figure, and the figure
//! would be wrong in a way nothing on it could show.
//!
//! So this reads bedGraph, four columns, 0-based and half-open, passed straight
//! through with no conversion.
//!
//! # A base no row covers is not a base scoring nought
//!
//! Rows are what the file states, and the gaps between them are what it does
//! not. They come back as gaps, and
//! [`DynseqTrack::from_pairs`](crate::DynseqTrack::from_pairs) leaves them
//! unscored rather than filling them in.

use crate::Region;

use super::{columns, lines, number, ReadError};

/// The scores a file holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Scores {
    /// Position and score, one per base the file states.
    pub pairs: Vec<(u64, f64)>,
    /// Rows in the file, before any filter.
    pub records: usize,
    /// Rows naming another sequence.
    pub other_sequence: usize,
    /// Rows on this sequence that touch no base of the window.
    pub off_region: usize,
    /// Rows whose value is a stated no score.
    pub no_score: usize,
}

/// Reads per-base attribution scores out of a bedGraph.
///
/// A row covering several bases gives that score to each of them, which is what
/// a bedGraph means. Only the part of a row inside `region` is expanded, so a
/// genome-wide file does not become a genome-wide vector.
///
/// # Errors
///
/// Returns the first row without four columns, whose coordinates are not
/// numbers, or whose end is before its start.
pub fn scores(text: &str, region: &Region) -> Result<Scores, ReadError> {
    let mut found = Scores::default();

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 4 {
            return Err(ReadError::at(
                at,
                format!(
                    "a bedGraph row has four columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        found.records += 1;
        if cols[0] != region.seq() {
            found.other_sequence += 1;
            continue;
        }

        let start: u64 = number(cols[1], "start", at)?;
        let end: u64 = number(cols[2], "end", at)?;
        if end < start {
            return Err(ReadError::at(at, "end is before start"));
        }
        if end <= region.start() || start >= region.end() {
            found.off_region += 1;
            continue;
        }

        // A value that is not a number is a base nobody scored, and it stays
        // one: emitted as a pair it would reach a glyph height.
        let Ok(score) = cols[3].parse::<f64>() else {
            found.no_score += 1;
            continue;
        };
        if !score.is_finite() {
            found.no_score += 1;
            continue;
        }

        let from = start.max(region.start());
        let to = end.min(region.end());
        for pos in from..to {
            found.pairs.push((pos, score));
        }
    }

    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> Region {
        Region::new("chr1", 1_000, 1_010).unwrap()
    }

    #[test]
    fn a_row_gives_its_score_to_every_base_it_covers() {
        let found = scores("chr1\t1000\t1003\t0.42\n", &region()).expect("scores");
        assert_eq!(found.pairs, [(1_000, 0.42), (1_001, 0.42), (1_002, 0.42)]);
    }

    #[test]
    fn only_the_part_inside_the_window_is_expanded() {
        // A genome-wide file must not become a genome-wide vector.
        let found = scores("chr1\t0\t1000000\t0.1\n", &region()).expect("scores");
        assert_eq!(found.pairs.len(), 10);
        assert_eq!(found.pairs[0].0, 1_000);
        assert_eq!(found.pairs[9].0, 1_009);
    }

    #[test]
    fn a_negative_score_survives_and_a_value_that_is_not_a_number_does_not() {
        let found = scores(
            "chr1\t1000\t1001\t-0.5\nchr1\t1001\t1002\tNaN\nchr1\t1002\t1003\tNA\n",
            &region(),
        )
        .expect("scores");
        assert_eq!(found.pairs, [(1_000, -0.5)]);
        assert_eq!(found.no_score, 2);
        assert_eq!(found.records, 3);
    }

    #[test]
    fn rows_elsewhere_are_counted_rather_than_read() {
        let found = scores(
            "chr2\t1000\t1001\t0.5\nchr1\t50\t60\t0.5\nchr1\t1000\t1001\t0.5\n",
            &region(),
        )
        .expect("scores");
        assert_eq!(found.other_sequence, 1);
        assert_eq!(found.off_region, 1);
        assert_eq!(found.pairs.len(), 1);
    }

    #[test]
    fn a_row_that_is_not_a_bedgraph_row_is_refused_by_line() {
        let error = scores("chr1\t1000\t1001\n", &region()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("four columns"), "{error}");

        let error = scores("chr1\t1005\t1001\t0.5\n", &region()).unwrap_err();
        assert!(error.reason.contains("end is before start"), "{error}");
    }
}
