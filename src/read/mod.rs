//! Readers for the formats a genomics shell already produces.
//!
//! Every format here is line based text, which is what keeps the dependency
//! count at zero. Binary formats are not read at all: BAM, CRAM and BCF come
//! in through a pipe, as `samtools depth`, `samtools view` or `bcftools view`
//! already write exactly what these readers take.
//!
//! # Text in, values out, and no path anywhere
//!
//! Every function here takes a `&str` and returns the vector a track is built
//! from. None of them opens a file, which is deliberate: the caller decides
//! where the text came from, so a path, standard input, an HTTP response and a
//! string in a test are all the same to a reader, and nothing in this crate
//! needs permission to read a disk.
//!
//! ```
//! use karyon::{plot, read};
//!
//! let bed = "chr1\t1000\t1500\tgeneA\t0\t+\n";
//! let region = karyon::Region::parse("chr1:1-2,000")?;
//! let features = read::interval::features(bed, &region, None)?;
//!
//! let svg = plot("chr1:1-2,000")?.add_features(features).to_svg();
//! assert!(svg.contains("geneA"));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # The one place the convention is not uniform
//!
//! Everywhere else in the project a position is 0-based and half-open. Here it
//! is whatever the specification of the file says, and a reader that takes it
//! for the other one is out by a base and quiet about it, so every reader
//! states which convention it is reading and every one has a test in `audit`
//! that pins a known position.
//!
//! That is also what `tests/properties.rs` checks without picking the
//! position: the same interval written as BED and as GFF3 has to come back as
//! the same pair of numbers, whichever base it starts on.
//!
//! # Which formats count from one
//!
//! Half-open and 0-based, passed straight through: BED, bedGraph, cytoBand.
//!
//! Inclusive and 1-based, so one is taken off the start on the way in: GFF3,
//! VCF, SAM and the output of `samtools depth`.
//!
//! # Silence is for rows that were never ours
//!
//! A whole file is handed over and one locus is asked for, so a row naming
//! another sequence and a row outside the window are the ordinary case and go
//! past without a word. A row that will not parse is the other kind. It says
//! the file is not what the flag claimed it was, and reading on would hand
//! someone a figure with data missing from it and nothing said, so the read
//! stops there and the [`ReadError`] carries the line it stopped on.

pub mod align;
pub mod align_pairs;
#[cfg(test)]
mod audit;
pub mod bisulfite;
pub mod clade;
pub mod domain;
pub mod interval;
pub mod locus;
pub mod methyl;
pub mod point;
pub mod seq;
pub mod sheet;
pub mod signal;
pub mod split;
pub mod structural;
pub mod table;

use std::fmt;

/// How a file should be read, when the shape of it is not enough to tell.
///
/// Several of these formats are four columns of text that differ only in what
/// the columns mean, and two of them count from different places. Guessing
/// between them is how a signal ends up drawn one base to the left of where it
/// was measured, so where the shape is ambiguous the reader is told rather
/// than left to work it out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// `chrom start end value`, 0-based half-open.
    BedGraph,
    /// `chrom pos depth`, the 1-based output of `samtools depth`.
    Depth,
    /// One value per line, starting at the left edge of the region.
    Values,
    /// `chrom start end [name] [score] [strand]`, 0-based half-open.
    Bed,
    /// Nine columns, 1-based inclusive, attributes in the ninth.
    Gff3,
}

impl Format {
    /// Parses a format name, as `--format` spells it.
    ///
    /// The aliases are the ones the files themselves come labelled with: `bg`
    /// for bedGraph, and `gtf` for GFF3, which is close enough in the columns
    /// these readers use.
    pub fn parse(word: &str) -> Option<Format> {
        Some(match word {
            "bedgraph" | "bg" => Format::BedGraph,
            "depth" => Format::Depth,
            "values" => Format::Values,
            "bed" => Format::Bed,
            "gff3" | "gff" | "gtf" => Format::Gff3,
            _ => return None,
        })
    }
}

/// What went wrong while reading, and where.
#[derive(Debug)]
pub struct ReadError {
    /// The line number, counting from one, or zero for a whole file problem.
    pub line: usize,
    /// What was wrong with it.
    pub reason: String,
}

impl ReadError {
    /// An error on a numbered line.
    pub fn at(line: usize, reason: impl Into<String>) -> Self {
        ReadError {
            line,
            reason: reason.into(),
        }
    }

    /// An error about the file as a whole.
    pub fn whole(reason: impl Into<String>) -> Self {
        ReadError {
            line: 0,
            reason: reason.into(),
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(f, "{}", self.reason)
        } else {
            write!(f, "line {}: {}", self.line, self.reason)
        }
    }
}

impl std::error::Error for ReadError {}

/// The lines worth looking at, numbered from one.
///
/// Blank lines and the comment and header lines every one of these formats
/// carries are dropped: `#` for BED, GFF3, VCF and SAM headers, `@` for the
/// SAM header proper, and the `track` and `browser` lines a UCSC file starts
/// with.
///
/// A UCSC header has to carry a `key=value` to be taken as one, because a
/// sequence may be called `track`, and in a space separated file its rows are
/// otherwise indistinguishable from the header and disappear without a word.
///
/// A UTF-8 byte order mark at the very start of the text is dropped. It is a
/// mark on the file rather than the first character of the first field, and
/// tooling on Windows and spreadsheet exports write one, so leaving it in would
/// glue it to the first field of the first row and lose that row silently.
pub fn lines(text: &str) -> impl Iterator<Item = (usize, &str)> {
    // The mark carries no newline, so the numbering is unchanged by this.
    text.strip_prefix('\u{feff}')
        .unwrap_or(text)
        .lines()
        .enumerate()
        .map(|(index, line)| (index + 1, line.trim_end_matches(['\r'])))
        .filter(|(_, line)| {
            let line = line.trim();
            let ucsc =
                (line.starts_with("track ") || line.starts_with("browser ")) && line.contains('=');
            !line.is_empty() && !line.starts_with('#') && !line.starts_with('@') && !ucsc
        })
}

/// Splits on tabs, falling back to whitespace when there is not one.
///
/// Every format here is tab separated on paper, and space separated in about
/// half the files that exist.
pub fn columns(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').collect()
    } else {
        line.split_whitespace().collect()
    }
}

/// Parses a coordinate, saying which column failed rather than just failing.
pub fn number<T: std::str::FromStr>(value: &str, what: &str, line: usize) -> Result<T, ReadError> {
    value
        .parse::<T>()
        .map_err(|_| ReadError::at(line, format!("{what} is not a number: {value:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comments_headers_and_blank_lines_are_not_data() {
        let text = "\
#comment
track name=x

@HD\tVN:1.6
browser position=chr1
chr1\t10\t20\t5
";
        let kept: Vec<&str> = lines(text).map(|(_, line)| line).collect();
        assert_eq!(kept, vec!["chr1\t10\t20\t5"]);
    }

    #[test]
    fn a_sequence_called_track_keeps_its_rows() {
        // A UCSC header carries a key=value. A row of a sequence named `track`
        // in a space separated file does not, and used to disappear.
        let kept: Vec<&str> = lines("track 99 100 x\nbrowser 5 6 y\n")
            .map(|(_, line)| line)
            .collect();
        assert_eq!(kept, vec!["track 99 100 x", "browser 5 6 y"]);
    }

    #[test]
    fn a_byte_order_mark_does_not_swallow_the_first_row() {
        // A BED written by tooling that stamps a UTF-8 BOM on the file. The
        // mark used to stay glued to the first field, so `chr1` on line 1 read
        // as a different sequence and geneA was dropped without a word.
        let text = "\u{feff}chr1\t100\t200\tgeneA\nchr1\t300\t400\tgeneB\n";
        let kept: Vec<(usize, &str)> = lines(text).collect();
        assert_eq!(
            kept,
            vec![(1, "chr1\t100\t200\tgeneA"), (2, "chr1\t300\t400\tgeneB"),]
        );
        assert_eq!(columns(kept[0].1)[0], "chr1");
    }

    #[test]
    fn the_line_number_counts_the_lines_that_were_skipped() {
        let text = "#one\n\nchr1\t10\t20\t5\n";
        let (number, _) = lines(text).next().unwrap();
        assert_eq!(number, 3);
    }

    #[test]
    fn a_space_separated_file_still_reads() {
        assert_eq!(columns("chr1 10 20"), vec!["chr1", "10", "20"]);
        assert_eq!(columns("chr1\t10\t20"), vec!["chr1", "10", "20"]);
    }

    #[test]
    fn a_tab_file_with_spaces_in_a_field_keeps_the_field_whole() {
        assert_eq!(
            columns("chr1\t10\t20\tsome name"),
            vec!["chr1", "10", "20", "some name"]
        );
    }
}
