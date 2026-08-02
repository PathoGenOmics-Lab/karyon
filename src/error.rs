//! Error type returned by fallible constructors.
//!
//! One enum, four variants, no dependency and no source chain. Four things in
//! the crate can fail, and each of them is a value being built out of text or
//! out of names the caller chose: a region with nothing in it
//! ([`Region::new`](crate::Region::new)), a locus string that is not
//! `seq:start-end` ([`Region::parse`](crate::Region::parse)), a Newick string
//! that is not a tree ([`Tree::parse_newick`](crate::Tree::parse_newick)), and
//! a genome that names one sequence twice
//! ([`Genome::checked`](crate::Genome::checked)). Everything else that returns
//! a `Result` either forwards one of those, the way [`plot`](fn@crate::plot)
//! forwards the locus parse, or is writing a file.
//!
//! # No `?` in the middle of a figure
//!
//! Failure sits at the two ends of a program that draws: reading a locus or a
//! tree at the start, writing the SVG at the end. Nothing between them is
//! fallible, which is what lets a figure be one chain of builder calls and
//! [`Figure::to_svg`](crate::Figure::to_svg) hand back a `String`. The
//! conversion into [`std::io::Error`] below is what joins the two ends up.
//!
//! # What is not an error
//!
//! Data that cannot be drawn as asked is drawn as well as it can be. A coverage
//! array that runs out before the region does stops where it stops, a feature
//! that begins before the window keeps its real start and is clipped, and a
//! width too small to leave a plotting area, or not a number at all, is raised
//! to the smallest width that renders. Each of those could be a variant here,
//! and each would trade a figure with one clipped rectangle in it for no figure
//! at all.

use std::fmt;

/// Anything that can go wrong while building a figure.
///
/// Rendering itself never fails: once a [`Figure`](crate::Figure) exists, its
/// inputs have already been validated.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Two sequences of a genome share a name.
    DuplicateSequence {
        /// The name that appeared twice.
        name: String,
    },
    /// A region whose end is not strictly greater than its start.
    ///
    /// Regions are half-open, so a one-base region is `start..start + 1`.
    EmptyRegion {
        /// Requested start, 0-based.
        start: u64,
        /// Requested end, exclusive.
        end: u64,
    },
    /// A locus string that does not look like `seq:start-end`.
    InvalidLocus {
        /// The string that failed to parse.
        input: String,
        /// Why it failed.
        reason: &'static str,
    },
    /// A Newick string that does not describe a tree.
    InvalidNewick {
        /// Why it failed.
        reason: &'static str,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::DuplicateSequence { name } => {
                write!(f, "the genome has more than one sequence called {name}")
            }
            Error::EmptyRegion { start, end } => write!(
                f,
                "empty region: end ({end}) must be greater than start ({start})"
            ),
            Error::InvalidLocus { input, reason } => {
                write!(f, "invalid locus {input:?}: {reason}")
            }
            Error::InvalidNewick { reason } => write!(f, "invalid Newick tree: {reason}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<Error> for std::io::Error {
    /// So that a locus string and the file it renders to can share one `?`.
    ///
    /// A program that draws figures spends its errors on writing files, so its
    /// functions return [`std::io::Result`]. Without this, the one call that
    /// parses a region is the only thing in such a function that cannot use
    /// `?`, which is a poor reason to reach for `unwrap`.
    ///
    /// ```
    /// use karyon::plot;
    ///
    /// fn draw() -> std::io::Result<String> {
    ///     Ok(plot("chr1:1-1000")?.add_coverage(vec![30.0; 1000]).to_svg())
    /// }
    ///
    /// assert!(draw().unwrap().starts_with("<svg"));
    /// ```
    fn from(error: Error) -> Self {
        std::io::Error::new(std::io::ErrorKind::InvalidInput, error)
    }
}
