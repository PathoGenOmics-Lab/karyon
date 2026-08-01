//! Readers for the formats a genomics shell already produces.
//!
//! Every format here is line based text, which is what keeps the dependency
//! count at zero. Binary formats are not read at all: BAM, CRAM and BCF come
//! in through a pipe, as `samtools depth`, `samtools view` or `bcftools view`
//! already write exactly what these readers take.
//!
//! # Coordinates
//!
//! This is the one place in the project where the convention is not uniform,
//! and getting it wrong is silent, so every reader states which it is reading
//! and every one has a test that pins a known position.
//!
//! Half-open and 0-based, passed straight through: BED, bedGraph, cytoBand.
//!
//! Inclusive and 1-based, so one is taken off the start on the way in: GFF3,
//! VCF, SAM and the output of `samtools depth`.

pub mod align;
#[cfg(test)]
mod audit;
pub mod interval;
pub mod point;
pub mod seq;
pub mod signal;
pub mod table;

use std::fmt;

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
pub fn lines(text: &str) -> impl Iterator<Item = (usize, &str)> {
    text.lines()
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
