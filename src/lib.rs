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
//! # Circular sequences
//!
//! A plasmid, an organelle genome and most bacterial chromosomes have no ends,
//! and drawing one as a line puts an edge where the biology has none. [`Rings`] maps position to an angle instead, so
//! it is a container of its own rather than a track: annotation, composition
//! and variants go on concentric rings, and chords across the middle join the
//! two ends of a rearrangement. A [`Rings`] plot and a [`Figure`] can share one
//! [`Panels`] sheet, which is all the [`Drawing`] trait is for.
//!
//! # Whole genomes
//!
//! A figure is one region on one sequence, which is right for a locus and
//! wrong for an assembly. [`Genome`] lays several sequences end to end and
//! hands back the one region that covers them, so every track here works
//! across all of them at once: the genome-wide association figure is a
//! [`ManhattanTrack`] over a `Genome` with its
//! [`boundaries`](Genome::boundaries) handed to
//! [`bands`](ManhattanTrack::bands), and [`GenomeTrack`] draws the sequences
//! underneath so a reader can see where one ends and the next begins.
//!
//! # Protein coordinates
//!
//! A variant in a coding sequence is named by residue rather than by base:
//! BRAF V600E, TP53 R175H, rpoB S450L. [`CodonTrack`] is the
//! [`AxisTrack`] that can be pointed at with those names. It partitions a coding
//! sequence into codons, numbers them, translates them where there is room for a
//! letter, and counts from the far end on the reverse strand, which is where the
//! arithmetic usually goes wrong.
//!
//! # Extending
//!
//! The twenty-nine track types shipped here are implementations of one small trait,
//! [`Track`], with no privileged access to the figure. A track type the crate
//! does not have is around thirty lines: see the example on [`Track`].
//!
//! # Sequence logos
//!
//! [`LogoTrack`] draws the classic information content logo and six other
//! [`LogoScore`]s, five of which measure a symbol against a background and can
//! therefore hang it below the baseline. [`LogoTrack::edlogo`] is the one to
//! reach for first: log odds recentred on the column median, the plot that
//! shows a symbol which is absent rather than merely uncommon.
//!
//! Where the baseline sits is a separate choice, [`Centering`], and symbols
//! are arbitrary strings, so amino acid codes and k-mers plot as readily as
//! bases.
//!
//! How much of a logo to believe is [`LogoTrack::stabilize`], which shrinks
//! each column towards the background by the amount its sample size can
//! support. Without it, four aligned sequences that happen to agree produce the
//! same confident column as four thousand. The machinery is in [`dash`] and
//! works on any compositional data.

#![forbid(unsafe_code)]

pub mod dash;
pub mod error;
pub mod figure;
pub mod genome;
pub mod panels;
pub mod region;
pub mod rings;
pub mod scale;
pub mod svg;
pub mod theme;
pub mod track;
pub mod tree;

pub use crate::error::Error;
pub use crate::figure::{Figure, Margin};
pub use crate::genome::{Chromosome, Genome};
pub use crate::panels::Panels;
pub use crate::region::Region;
pub use crate::rings::{
    AxisRing, Drawing, FeatureRing, MarkerRing, Polar, Ring, RingContext, Rings, SignalRing,
};
pub use crate::scale::Scale;
pub use crate::svg::{Anchor, SvgWriter};
pub use crate::theme::{wash, BaseColors, Theme};
pub use crate::track::{
    strand_color, Aggregate, AlignmentBlock, Association, AxisTrack, Band, BisulfiteTrack,
    CellScale, Centering, CigarOp, CladeBlock, CladeTrack, CodonTrack, CoverageStyle,
    CoverageTrack, DotplotTrack, DrawContext, Feature, FeatureTrack, GeneShape, GenomeTrack,
    Homology, IdeogramTrack, Legend, LegendItem, LegendTrack, Locus, LocusTrack, LogoColumn,
    LogoScore, LogoStack, LogoTrack, ManhattanTrack, Marker, MatrixRow, MatrixTrack, MethylSite,
    MethylationTrack, Molecule, Move, MsaColoring, MsaDisplay, MsaSequence, MsaTrack, Orf,
    OrfTrack, PileupLayout, PileupTrack, Read, ReadColoring, Rect, Segment, SequenceTrack, SnpSite,
    SnpTrack, SplitRead, SplitReadTrack, SplitSegment, SquiggleTrack, StackOrder, Stain, Strand,
    StructuralTrack, StructuralVariant, SvKind, SyntenyTrack, TanglegramTrack, Terminator, Track,
    TranscriptionUnit, TranscriptionUnitTrack, TreeShape, TreeTrack, Variant, VariantStyle,
    VariantTrack, Window, WindowStyle, WindowTrack,
};
