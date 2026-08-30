//! Genomic track plots, rendered to standalone SVG.
//!
//! `karyon` draws the kind of figure a genome browser draws: a stack of tracks
//! over one shared coordinate axis, so a read depth profile, the reference
//! bases, the gene models and the variant calls all line up on the same
//! position. It has no runtime dependencies, does no I/O beyond an optional
//! `save_svg`, and emits plain SVG 1.1 that opens unchanged in a browser, in
//! Inkscape and in Illustrator.
//!
//! # 0-based inside, 1-based where a reader looks
//!
//! Positions are **0-based and half-open** everywhere, the BED convention. The
//! two exceptions are the ones a reader sees: [`Region::parse`] accepts the
//! 1-based inclusive locus strings that samtools and IGV use, and tick labels
//! are printed in that same 1-based form. A VCF `POS` or a GFF `start` is
//! therefore `pos - 1` on the way in.
//!
//! [`read`] does that subtraction so nobody has to remember which formats need
//! it, and `tests/properties.rs` checks that the same interval written in a
//! 0-based format and a 1-based one comes back as the same pair of numbers.
//!
//! # Files are text, and text is somebody else's problem
//!
//! [`read`] turns BED, bedGraph, GFF3, VCF, SAM, cytoBand, `samtools depth`,
//! FASTA and Newick into the vectors the tracks take. Every one of those
//! functions takes a `&str`: nothing in this crate opens a path to read one,
//! so where the text came from stays the caller's decision, and the dependency
//! count stays at zero because all nine formats are lines of text.
//!
//! # Four tracks, one region, one call each
//!
//! ```
//! use karyon::{plot, Feature, Strand, Variant};
//!
//! let depth: Vec<f64> = (0..201).map(|i| 55.0 - (i as f64 % 23.0)).collect();
//! let bases: Vec<u8> = b"ACGT".iter().cycle().take(201).copied().collect();
//!
//! let svg = plot("NC_000962.3:761000-761200")
//!     .unwrap()
//!     .title("rpoB resistance determining region")
//!     .add_coverage(depth)
//!     .label("depth")
//!     .add_sequence(bases)
//!     .label("reference")
//!     .add_features(vec![
//!         Feature::new(761_040, 761_160).name("RRDR").strand(Strand::Forward),
//!     ])
//!     .label("genes")
//!     .add_variants(vec![
//!         Variant::new(761_109).value(0.98).category("rpoB S450L"),
//!     ])
//!     .label("variants")
//!     .to_svg();
//!
//! assert!(svg.starts_with("<svg"));
//! ```
//!
//! # Neither layer draws what the other cannot
//!
//! [`plot()`] is the short way to write a figure down: one call per track, in
//! the order they stack, with the region held once and the coordinate ruler
//! filled in. [`Figure`] is what it builds, and the layer to reach for when a
//! track is built by an alternative constructor, read back before it is drawn,
//! or passed around. See the [`plot`](mod@crate::plot) module for where the
//! short form stops.
//!
//! # Drawing a circle as a line invents an edge
//!
//! A plasmid, an organelle genome and most bacterial chromosomes have no ends.
//! [`Rings`] maps position to an angle instead, so it is a container of its own
//! rather than a track: annotation, composition and variants go on concentric
//! rings, and chords across the middle join the two ends of a rearrangement. A
//! [`Rings`] plot and a [`Figure`] can share one [`Panels`] sheet, which is all
//! the [`Drawing`] trait is for.
//!
//! # Right for a locus, wrong for an assembly
//!
//! A figure is one region on one sequence. [`Genome`] lays several sequences
//! end to end and hands back the one region that covers them, so every track
//! here works across all of them at once: the genome-wide association figure
//! is a [`ManhattanTrack`] over a `Genome` with its
//! [`boundaries`](Genome::boundaries) handed to
//! [`bands`](ManhattanTrack::bands), and [`GenomeTrack`] draws the sequences
//! underneath so a reader can see where one ends and the next begins.
//!
//! # A residue is a coordinate too
//!
//! A variant in a coding sequence is named by residue rather than by base:
//! BRAF V600E, TP53 R175H, rpoB S450L. [`CodonTrack`] is the [`AxisTrack`] that
//! can be pointed at with those names. It partitions a coding sequence into
//! codons, numbers them, translates them where there is room for a letter, and
//! counts from the far end on the reverse strand, which is where the arithmetic
//! usually goes wrong.
//!
//! # A track type the crate lacks is thirty lines
//!
//! The thirty-six track types shipped here are implementations of one small
//! trait, [`Track`], with no privileged access to the figure. See the example
//! on [`Track`] for the whole of what a new one has to do.
//!
//! # Information content is not the only score
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

/// The version of this crate.
///
/// The command prints it for `--version`, and so does the WebAssembly build the
/// documentation site's playground runs. That build is a separate crate with a
/// version of its own, so reading `CARGO_PKG_VERSION` there answered for the
/// shim rather than for the library and reported `0.0.0` for a release.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod cli;
pub mod dash;
pub mod error;
pub mod figure;
pub mod genome;
pub mod map;
pub mod panels;
pub mod plot;
pub mod read;
pub mod region;
pub mod rings;
pub mod scale;
pub mod style;
pub mod svg;
pub mod theme;
pub mod track;
pub mod tree;

pub use crate::error::Error;
pub use crate::figure::{Figure, Margin};
pub use crate::genome::{Chromosome, Genome};
pub use crate::map::{
    GeoFlow, GeoLocation, GeoPosition, GeoProjection, Map, PhyloConnector, PhyloMap,
};
pub use crate::panels::Panels;
pub use crate::plot::{plot, Plot};
pub use crate::read::{Format, ReadError};
pub use crate::region::Region;
pub use crate::rings::{
    AxisRing, Drawing, FeatureRing, MarkerRing, Polar, Ring, RingContext, Rings, SignalRing,
};
pub use crate::scale::Scale;
pub use crate::style::{
    AxisFormat, Density, Emphasis, LinePattern, MarkStyle, QuantitativeAxis, ReferenceLine,
    RenderProfile, Symbol, VisualTokens,
};
pub use crate::svg::{fit_text, Anchor, SvgWriter};
pub use crate::theme::{wash, BaseColors, Theme};
pub use crate::track::{
    strand_color, Aggregate, AlignmentBlock, AncestralStateLayer, Association, AxisTrack, Band,
    BisulfiteTrack, BranchEventLayer, BranchGeometry, BranchIntervalLayer, BranchRateMixture,
    CellScale, Centering, CigarOp, CladeBlock, CladeHighlight, CladeTrack, CodonTrack, CopyNumber,
    CopyNumberSegment, CopyNumberTrack, CoverageStyle, CoverageTrack, DomainArchitecture,
    DomainFeature, DomainTrack, DotplotTrack, DrawContext, DynseqTrack, Feature, FeatureTrack,
    GeneShape, GenomeTrack, Homology, HomoplasyLayer, IdeogramTrack, Junction, JunctionTrack,
    Legend, LegendItem, LegendTrack, Locus, LocusTrack, LogoColumn, LogoScore, LogoStack,
    LogoTrack, ManhattanTrack, Marker, MatrixRow, MatrixTrack, MethylSite, MethylationTrack,
    Molecule, Motif, Move, MsaColoring, MsaDisplay, MsaSequence, MsaTrack, NodeGlyph,
    NodeGlyphStyle, NodeGlyphTarget, Orf, OrfTrack, PhylodynamicPoint, PhylodynamicScale,
    PhylodynamicTrack, PileupLayout, PileupTrack, RadialDirection, Read, ReadColoring, Rect,
    Segment, SelectionEvidence, SelectionSite, SelectionTrack, SequenceTrack, SnpSite, SnpTrack,
    SplitRead, SplitReadTrack, SplitSegment, SquiggleTrack, StackOrder, Stain, Strand,
    StructuralTrack, StructuralVariant, SupportStyle, SurveillanceMetric, SurveillanceObservation,
    SurveillanceStyle, SurveillanceTrack, SvKind, SyntenyTrack, TangleLabels, TangleTieStyle,
    TanglegramTrack, Terminator, Track, TraitColumn, TraitScale, TraitStyle, TranscriptionUnit,
    TranscriptionUnitTrack, TreeProjection, TreeShape, TreeTrack, Variant, VariantStyle,
    VariantTrack, Window, WindowStyle, WindowTrack,
};
pub use crate::tree::{
    AnnotationValue, Annotations, Clade, Mutation, Mutations, Placement, TimeDirection, Tree,
};
