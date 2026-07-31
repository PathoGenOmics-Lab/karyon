//! Genomic intervals and the coordinate convention used throughout the crate.

use std::fmt;

use crate::error::Error;

/// A genomic interval on one sequence.
///
/// # Coordinates
///
/// Internally a region is **0-based, half-open**, the BED convention: the first
/// base of a sequence is `0`, and `end` is one past the last included base. This
/// is the convention every constructor and accessor of this crate uses, unless
/// its documentation says otherwise.
///
/// Two places speak the other dialect, because that is what a user reads on
/// screen: [`Region::parse`] accepts the **1-based inclusive** locus strings
/// that samtools, IGV and the GFF/VCF formats use, and [`Display`](fmt::Display)
/// prints them back in that same form. So `Region::parse("chr1:101-200")`
/// yields `start = 100`, `end = 200`, a 100 bp region, and prints again as
/// `chr1:101-200`.
///
/// ```
/// use karyon::Region;
///
/// let region = Region::parse("NC_000962.3:761,100-761,200").unwrap();
/// assert_eq!(region.start(), 761_099); // 0-based
/// assert_eq!(region.end(), 761_200); // exclusive
/// assert_eq!(region.len(), 101);
/// assert_eq!(region.to_string(), "NC_000962.3:761100-761200");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    seq: String,
    start: u64,
    end: u64,
}

impl Region {
    /// Builds a region from 0-based half-open coordinates.
    ///
    /// # Errors
    ///
    /// Returns [`Error::EmptyRegion`] when `end <= start`. A figure needs at
    /// least one base to draw.
    pub fn new(seq: impl Into<String>, start: u64, end: u64) -> Result<Self, Error> {
        if end <= start {
            return Err(Error::EmptyRegion { start, end });
        }
        Ok(Region {
            seq: seq.into(),
            start,
            end,
        })
    }

    /// Parses a `seq:start-end` locus string in the 1-based inclusive
    /// convention used by samtools and IGV.
    ///
    /// Thousands separators (`,` and `_`) are accepted and ignored, and the
    /// sequence name may itself contain colons: the split happens at the last
    /// colon.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidLocus`] when the string has no colon, no dash in
    /// the coordinate part, non-numeric coordinates, a start below 1, or an end
    /// before the start.
    pub fn parse(locus: &str) -> Result<Self, Error> {
        let invalid = |reason: &'static str| Error::InvalidLocus {
            input: locus.to_string(),
            reason,
        };

        let (seq, span) = locus
            .rsplit_once(':')
            .ok_or_else(|| invalid("expected seq:start-end"))?;
        if seq.is_empty() {
            return Err(invalid("sequence name is empty"));
        }
        let (start, end) = span
            .split_once('-')
            .ok_or_else(|| invalid("expected start-end after the colon"))?;

        let number = |raw: &str| -> Result<u64, Error> {
            let cleaned: String = raw
                .chars()
                .filter(|c| *c != ',' && *c != '_' && !c.is_whitespace())
                .collect();
            if cleaned.is_empty() {
                return Err(invalid("missing coordinate"));
            }
            cleaned
                .parse::<u64>()
                .map_err(|_| invalid("coordinates must be non-negative integers"))
        };

        let start = number(start)?;
        let end = number(end)?;
        if start == 0 {
            return Err(invalid("1-based coordinates start at 1, not 0"));
        }
        if end < start {
            return Err(invalid("end is before start"));
        }
        // 1-based inclusive to 0-based half-open.
        Region::new(seq, start - 1, end)
    }

    /// Name of the sequence this region sits on.
    pub fn seq(&self) -> &str {
        &self.seq
    }

    /// First base of the region, 0-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// One past the last base of the region.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Length of the region in bases. Always at least 1.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Always `false`: a region cannot be constructed empty.
    ///
    /// Present only because clippy asks for it next to [`Region::len`].
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Whether a 0-based position falls inside the region.
    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos < self.end
    }

    /// First base of the region as it is shown to users, 1-based.
    pub fn display_start(&self) -> u64 {
        self.start + 1
    }

    /// Last base of the region as it is shown to users, 1-based inclusive.
    pub fn display_end(&self) -> u64 {
        self.end
    }
}

impl fmt::Display for Region {
    /// Prints the 1-based inclusive locus string, the form [`Region::parse`]
    /// reads back.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}-{}",
            self.seq,
            self.display_start(),
            self.display_end()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_rejects_empty_and_inverted_regions() {
        assert_eq!(
            Region::new("chr1", 100, 100),
            Err(Error::EmptyRegion {
                start: 100,
                end: 100
            })
        );
        assert!(Region::new("chr1", 100, 99).is_err());
        assert!(Region::new("chr1", 100, 101).is_ok());
    }

    #[test]
    fn parse_converts_from_one_based_inclusive() {
        let r = Region::parse("chr1:1-10").unwrap();
        assert_eq!(r.start(), 0);
        assert_eq!(r.end(), 10);
        assert_eq!(r.len(), 10);
    }

    #[test]
    fn parse_accepts_thousands_separators() {
        let a = Region::parse("chr1:1,000-2,000").unwrap();
        let b = Region::parse("chr1:1_000-2_000").unwrap();
        assert_eq!(a, b);
        assert_eq!(a.start(), 999);
        assert_eq!(a.end(), 2000);
    }

    #[test]
    fn parse_splits_at_the_last_colon() {
        let r = Region::parse("gi|123|ref|NC_1.1:5-8").unwrap();
        assert_eq!(r.seq(), "gi|123|ref|NC_1.1");
        assert_eq!(r.start(), 4);
    }

    #[test]
    fn parse_rejects_malformed_input() {
        for bad in [
            "chr1",
            "chr1:100",
            ":100-200",
            "chr1:0-200",
            "chr1:200-100",
            "chr1:a-b",
            "chr1:-200",
        ] {
            assert!(Region::parse(bad).is_err(), "{bad} should not parse");
        }
    }

    #[test]
    fn display_round_trips_through_parse() {
        let r = Region::parse("NC_000962.3:761100-761200").unwrap();
        assert_eq!(Region::parse(&r.to_string()).unwrap(), r);
    }

    #[test]
    fn contains_respects_the_half_open_end() {
        let r = Region::new("chr1", 10, 20).unwrap();
        assert!(r.contains(10));
        assert!(r.contains(19));
        assert!(!r.contains(20));
        assert!(!r.contains(9));
    }
}
