//! Genomic track plots, rendered to standalone SVG.
//!
//! `karyon` draws the kind of figure a genome browser draws: a stack of tracks
//! over one shared coordinate axis, so a read depth profile, the reference
//! bases, the gene models and the variant calls all line up on the same
//! position. It has no runtime dependencies, does no I/O beyond an optional
//! `save_svg`, and emits plain SVG 1.1 that opens unchanged in a browser, in
//! Inkscape and in Illustrator.
//!
//! # Coordinates
//!
//! Positions are **0-based and half-open** everywhere, the BED convention. The
//! two exceptions are the ones a reader sees: [`Region::parse`] accepts the
//! 1-based inclusive locus strings that samtools and IGV use, and tick labels
//! are printed in that same 1-based form. A VCF `POS` or a GFF `start` is
//! therefore `pos - 1` on the way in.
//!
//! # Example
//!
//! ```
//! use karyon::{
//!     AxisTrack, CoverageTrack, Feature, FeatureTrack, Figure, Region,
//!     SequenceTrack, Strand, Variant, VariantTrack,
//! };
//!
//! let region = Region::parse("NC_000962.3:761000-761200").unwrap();
//!
//! let depth: Vec<f64> = (0..201).map(|i| 55.0 - (i as f64 % 23.0)).collect();
//! let bases: Vec<u8> = b"ACGT".iter().cycle().take(201).copied().collect();
//!
//! let svg = Figure::new(region)
//!     .title("rpoB resistance determining region")
//!     .push(CoverageTrack::new(760_999, depth).label("depth"))
//!     .push(SequenceTrack::new(760_999, bases).label("reference"))
//!     .push(FeatureTrack::new(vec![
//!         Feature::new(761_040, 761_160).name("RRDR").strand(Strand::Forward),
//!     ]).label("genes"))
//!     .push(VariantTrack::new(vec![
//!         Variant::new(761_109).value(0.98).category("rpoB S450L"),
//!     ]).label("variants"))
//!     .push(AxisTrack::new())
//!     .to_svg();
//!
//! assert!(svg.starts_with("<svg"));
//! ```
//!
//! # Extending
//!
//! The six track types shipped here are implementations of one small trait,
//! [`Track`], with no privileged access to the figure. A track type the crate
//! does not have is around thirty lines: see the example on [`Track`].
//!
//! # Sequence logos
//!
//! [`LogoTrack`] draws the classic information content logo and, through
//! [`LogoScaling::EnrichmentDepletion`], the enrichment and depletion variant
//! that stacks symbols on both sides of a baseline. Only the second one can
//! show a symbol that is absent rather than merely uncommon, which is the
//! whole argument for it. Symbols are arbitrary strings, so amino acid codes
//! and k-mers plot as readily as bases.

#![forbid(unsafe_code)]

pub mod error;
pub mod figure;
pub mod region;
pub mod scale;
pub mod svg;
pub mod theme;
pub mod track;

pub use crate::error::Error;
pub use crate::figure::{Figure, Margin};
pub use crate::region::Region;
pub use crate::scale::Scale;
pub use crate::svg::{Anchor, SvgWriter};
pub use crate::theme::{BaseColors, Theme};
pub use crate::track::{
    Aggregate, AxisTrack, CoverageStyle, CoverageTrack, DrawContext, Feature, FeatureTrack,
    LogoColumn, LogoScaling, LogoStack, LogoTrack, Rect, SequenceTrack, StackOrder, Strand, Track,
    Variant, VariantStyle, VariantTrack,
};
