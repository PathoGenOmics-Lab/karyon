//! A short way to build a figure: one call per track, in the order they stack.
//!
//! [`Figure`] is the layer underneath. It takes tracks that are already built,
//! which is what a track passed around, stored, or returned by a function
//! needs. Most figures are not like that. Most are a handful of tracks written
//! out once, in place, and there the plumbing is the bulk of the code: a
//! [`Region`] to unwrap, a `use` line naming every track type, a `push` around
//! each one, and the window start repeated on every track that takes one even
//! though the figure already holds it.
//!
//! [`plot`] holds the region and names the track types, so none of that is
//! written twice.
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
//!     .add_variants(vec![Variant::new(761_109).value(0.98)])
//!     .label("variants")
//!     .to_svg();
//!
//! assert!(svg.starts_with("<svg"));
//! ```
//!
//! # One track is pending at a time
//!
//! A plot holds the track the last `add_` built, so the call after it still
//! talks about that track. [`Plot::label`] names it, and [`Plot::adjust`] hands
//! it to a closure for anything else:
//!
//! ```
//! use karyon::{plot, Aggregate};
//!
//! let ratios = plot("chr1:1-1000")
//!     .unwrap()
//!     .add_coverage(vec![30.0; 1000])
//!     .label("depth")
//!     .adjust(|track| track.aggregate(Aggregate::Min).height(70.0));
//!
//! assert!(ratios.to_svg().contains(">depth<"));
//! ```
//!
//! The closure is handed the concrete track type, so every builder method on it
//! is in reach and a name that is not on it fails to compile. Which track it
//! reaches is a matter of where the call sits, and that part is not checked: an
//! `adjust` written one `add_` too late configures the next track instead, and
//! compiles, because `height` and `label` are on nearly every track.
//!
//! Adding the next track puts the previous one away, and so does anything that
//! renders. [`Plot::save`] hands back a plot with nothing pending, and so does
//! [`Plot::add_track`].
//!
//! # Two things happen without being asked
//!
//! Both can be undone.
//!
//! An [`AxisTrack`] goes on the bottom, since a figure without coordinates
//! along it is rarely what anyone meant. [`Plot::add_axis`] puts it somewhere
//! else instead, and [`Plot::remove_axis`] leaves it out.
//!
//! Tracks that lay an array along the axis start at the left edge of the
//! region. [`Plot::add_coverage`] is the depth of the region on display, so
//! there is no start to pass. When the array starts somewhere else, the `_at`
//! form takes the start: [`Plot::add_coverage_at`].
//!
//! What the axis counts is not always bases. An alignment is indexed by column
//! and a raw signal by sample, so a plot holding one of those wants a region in
//! that unit, as in `plot("alignment:1-320")`. The region is a coordinate
//! system, not a claim about a genome.
//!
//! # Where a plot stops and what takes over
//!
//! A plot is a stack of tracks over one axis and nothing more. [`Panels`] for a
//! sheet of figures and [`Rings`] for a circular sequence are still there, and
//! [`Plot::into_figure`] hands over the [`Figure`] it built so a plot can be the
//! start of one of those. A track that an alternative constructor built, or one
//! read back before it is drawn, needs none of that: [`Plot::add_track`] takes a
//! finished track without leaving the plot.
//!
//! [`Panels`]: crate::Panels
//! [`Rings`]: crate::Rings

use std::fmt;
use std::io;
use std::path::Path;

use crate::error::Error;
use crate::figure::{Figure, Margin};
use crate::genome::Genome;
use crate::region::Region;
use crate::style::{Density, RenderProfile};
use crate::theme::Theme;
use crate::track::{
    AlignmentBlock, Association, AxisTrack, Band, BisulfiteTrack, CladeBlock, CladeTrack,
    CodonTrack, CopyNumberSegment, CopyNumberTrack, CoverageTrack, DomainArchitecture, DomainTrack,
    DotplotTrack, Feature, FeatureTrack, GenomeTrack, IdeogramTrack, Legend, LegendTrack, Locus,
    LocusTrack, LogoColumn, LogoTrack, ManhattanTrack, MatrixRow, MatrixTrack, MethylSite,
    MethylationTrack, Molecule, MsaSequence, MsaTrack, OrfTrack, PhylodynamicPoint,
    PhylodynamicTrack, PileupTrack, Read, SelectionSite, SelectionTrack, SequenceTrack, SnpSite,
    SnpTrack, SplitRead, SplitReadTrack, SquiggleTrack, Strand, StructuralTrack, StructuralVariant,
    SurveillanceObservation, SurveillanceTrack, SyntenyTrack, TanglegramTrack, Track,
    TranscriptionUnit, TranscriptionUnitTrack, TreeTrack, Variant, VariantTrack, Window,
    WindowTrack,
};
use crate::tree::Tree;

/// Starts a plot over the region named by a locus string.
///
/// The string is the 1-based inclusive form samtools and IGV use, as parsed by
/// [`Region::parse`]. Start from a [`Region`] you already hold with
/// [`Plot::over`].
///
/// ```
/// use karyon::plot;
///
/// let svg = plot("chr7:140,753,000-140,754,000").unwrap().to_svg();
/// assert!(svg.starts_with("<svg"));
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidLocus`] for anything [`Region::parse`] will not
/// take: no colon, no dash after it, coordinates that are not numbers, a start
/// below 1, or an end before the start.
pub fn plot(locus: &str) -> Result<Plot, Error> {
    Region::parse(locus).map(Plot::over)
}

/// The empty slot at the head of a plot, with no track pending.
///
/// It is the default of [`Plot`]'s type parameter, so it stands for a fresh
/// plot and for whatever [`Plot::add_track`], [`Plot::done`] and [`Plot::save`]
/// hand back. Naming it is never necessary: `Plot` on its own already means
/// `Plot<Empty>`. It is deliberately not a [`Track`], so a figure cannot be
/// asked to draw it:
///
/// ```compile_fail
/// use karyon::{plot::Empty, Figure, Region};
///
/// Figure::new(Region::parse("chr1:1-100").unwrap()).push(Empty);
/// ```
///
/// and neither [`Plot::adjust`] nor [`Plot::label`] is on a plot holding one,
/// so naming a track before there is one is caught at compile time:
///
/// ```compile_fail
/// karyon::plot("chr1:1-100").unwrap().label("depth");
/// ```
///
/// ```compile_fail
/// karyon::plot("chr1:1-100").unwrap().adjust(|track| track);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Empty;

mod sealed {
    /// Closes [`Slot`](super::Slot) and [`Named`](super::Named) to this crate.
    pub trait Sealed {}
}

/// What a [`Plot`] can hold in its pending slot.
///
/// Implemented for every track type the crate ships and for [`Empty`], which
/// holds nothing. It is sealed: naming it in a bound is fine, implementing it
/// is not. A track from outside the crate goes in through [`Plot::add_track`].
pub trait Slot: sealed::Sealed + Sized {
    /// Hands the track to `figure`, or leaves the figure alone for [`Empty`].
    fn push_onto(self, figure: Figure) -> Figure;
}

impl sealed::Sealed for Empty {}

impl Slot for Empty {
    fn push_onto(self, figure: Figure) -> Figure {
        figure
    }
}

/// A track that can be given a name for the left gutter.
///
/// Implemented for every track type the crate ships, which is what lets
/// [`Plot::label`] name the track it just added without being told the type.
/// It is the same builder method the track already has, under a name that does
/// not collide with [`Track::label`]. [`Empty`] does not implement it, so
/// `label` and [`Plot::adjust`] are not on a plot with nothing pending.
pub trait Named: Slot {
    /// Returns the track with `label` shown in the left gutter.
    fn named(self, label: String) -> Self;
}

macro_rules! tracks {
    ($($track:ty),+ $(,)?) => {
        $(
            impl sealed::Sealed for $track {}

            impl Slot for $track {
                fn push_onto(self, figure: Figure) -> Figure {
                    figure.push(self)
                }
            }

            impl Named for $track {
                fn named(self, label: String) -> Self {
                    self.label(label)
                }
            }
        )+
    };
}

tracks![
    AxisTrack,
    BisulfiteTrack,
    CladeTrack,
    CodonTrack,
    CopyNumberTrack,
    CoverageTrack,
    DomainTrack,
    DotplotTrack,
    FeatureTrack,
    GenomeTrack,
    IdeogramTrack,
    LegendTrack,
    LocusTrack,
    LogoTrack,
    ManhattanTrack,
    MatrixTrack,
    MethylationTrack,
    MsaTrack,
    OrfTrack,
    PhylodynamicTrack,
    PileupTrack,
    SelectionTrack,
    SequenceTrack,
    SnpTrack,
    SplitReadTrack,
    SquiggleTrack,
    StructuralTrack,
    SurveillanceTrack,
    SyntenyTrack,
    TanglegramTrack,
    TranscriptionUnitTrack,
    TreeTrack,
    VariantTrack,
    WindowTrack,
];

/// A figure under construction, holding the track that was added last.
///
/// The type parameter is that track. It is what makes [`Plot::adjust`] see the
/// concrete type, and it changes with every `add_` call, so it is rarely worth
/// writing down. `Plot` on its own means `Plot<Empty>` and will not take a plot
/// mid-chain: a function that wants one is generic over `T: Slot`, or takes the
/// [`Figure`] that [`Plot::into_figure`] hands back.
///
/// See the [module documentation](self) for what a plot fills in on its own.
pub struct Plot<T = Empty> {
    figure: Figure,
    /// The track added last, not yet handed to the figure. It waits here so
    /// `label` and `adjust` can still reach it.
    pending: Option<T>,
    /// Whether an axis still needs to be appended at the end.
    wants_axis: bool,
}

impl<T> fmt::Debug for Plot<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Plot")
            .field("region", &self.figure.region().to_string())
            .field("tracks", &self.figure.track_count())
            .field("pending", &self.pending.is_some())
            .finish()
    }
}

impl Plot<Empty> {
    /// Starts a plot over `region`.
    ///
    /// [`plot`] is the same thing from a locus string.
    pub fn over(region: Region) -> Self {
        Plot {
            figure: Figure::new(region),
            pending: None,
            wants_axis: true,
        }
    }

    /// Wraps a figure that already holds tracks.
    ///
    /// The axis is not filled in: a figure built by hand is taken as finished.
    fn resume(figure: Figure) -> Self {
        Plot {
            figure,
            pending: None,
            wants_axis: false,
        }
    }

    /// Puts `track` in the pending slot. The slot is known to be free.
    fn park<N: Slot>(self, track: N) -> Plot<N> {
        Plot {
            figure: self.figure,
            pending: Some(track),
            wants_axis: self.wants_axis,
        }
    }
}

impl<T: Slot> Plot<T> {
    /// Hands the pending track to the figure, leaving the slot free.
    fn settle(self) -> Plot<Empty> {
        let Plot {
            mut figure,
            pending,
            wants_axis,
        } = self;
        if let Some(track) = pending {
            figure = track.push_onto(figure);
        }
        Plot {
            figure,
            pending: None,
            wants_axis,
        }
    }

    /// Sets the image width in pixels. Defaults to 900.
    pub fn width(mut self, width: f64) -> Self {
        self.figure = self.figure.width(width);
        self
    }

    /// Sets the title drawn above the tracks.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.figure = self.figure.title(title);
        self
    }

    /// Replaces the theme. [`Theme::dark`] is the other one that ships.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.figure = self.figure.theme(theme);
        self
    }

    /// Applies a named palette, type scale and track density together.
    pub fn profile(mut self, profile: RenderProfile) -> Self {
        self.figure = self.figure.profile(profile);
        self
    }

    /// Scales fonts, margins, label gutters, gaps and rounded corners together.
    pub fn visual_scale(mut self, factor: f64) -> Self {
        self.figure = self.figure.visual_scale(factor);
        self
    }

    /// Sets how tightly repeated data rows are packed.
    pub fn density(mut self, density: Density) -> Self {
        self.figure = self.figure.density(density);
        self
    }

    /// Replaces the whitespace around the plotting area.
    pub fn margin(mut self, margin: Margin) -> Self {
        self.figure = self.figure.margin(margin);
        self
    }

    /// Sets the width of the left gutter holding track labels.
    pub fn label_width(mut self, width: f64) -> Self {
        self.figure = self.figure.label_width(width);
        self
    }

    /// Sets the vertical gap between tracks in pixels.
    pub fn track_gap(mut self, gap: f64) -> Self {
        self.figure = self.figure.track_gap(gap);
        self
    }

    /// Leaves out the locus string printed at the top right.
    ///
    /// Worth doing when the horizontal axis is not genomic coordinates, as it
    /// is not for a squiggle or an alignment column index.
    pub fn remove_region_label(mut self) -> Self {
        self.figure = self.figure.show_region_label(false);
        self
    }

    /// Leaves out the axis this plot would otherwise append at the bottom.
    ///
    /// A tanglegram or a bare tree has no coordinate to put ticks on. It
    /// suppresses the automatic axis only: an axis put somewhere with
    /// [`Plot::add_axis`], or handed over to [`Plot::add_track`], stays where
    /// it was put.
    pub fn remove_axis(mut self) -> Self {
        self.wants_axis = false;
        self
    }

    /// The region the plot is drawn over.
    pub fn region(&self) -> &Region {
        self.figure.region()
    }

    /// Hands over the figure, with the pending track and the axis in place.
    ///
    /// This is the way out to [`Panels`](crate::Panels), to
    /// [`Figure::to_svg_with_id_prefix`], or to anything else that wants a
    /// [`Figure`].
    pub fn into_figure(self) -> Figure {
        let plot = self.settle();
        if plot.wants_axis {
            plot.figure.push(AxisTrack::new())
        } else {
            plot.figure
        }
    }

    /// Renders the figure to a standalone SVG document.
    pub fn to_svg(self) -> String {
        self.into_figure().to_svg()
    }

    /// Writes the SVG to `path`, and gives the plot back.
    ///
    /// Handing the plot back is what lets one stack be rendered twice:
    ///
    /// ```no_run
    /// use karyon::{plot, Theme};
    ///
    /// plot("plasmid:1-2000")
    ///     .unwrap()
    ///     .add_coverage(vec![30.0; 2000])
    ///     .save("example.svg")?
    ///     .theme(Theme::dark())
    ///     .save("example-dark.svg")?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    ///
    /// A save closes the stack. It puts the pending track away and fills in the
    /// axis, so saving twice does not draw two of them, and a track added after
    /// a save sits below the axis rather than above it.
    ///
    /// # Errors
    ///
    /// Returns whatever the write returned.
    pub fn save(self, path: impl AsRef<Path>) -> io::Result<Plot<Empty>> {
        let figure = self.into_figure();
        figure.save_svg(path)?;
        Ok(Plot::resume(figure))
    }

    /// Puts the pending track away without adding another.
    ///
    /// The plot's type stops naming it, which is what a stack built in a loop
    /// or behind a condition needs, since every arm has to have one type:
    ///
    /// ```
    /// use karyon::plot;
    ///
    /// let mut figure = plot("chr2:1-1000").unwrap();
    /// for (name, depth) in [("tumour", 40.0), ("normal", 22.0)] {
    ///     figure = figure.add_coverage(vec![depth; 1000]).label(name).done();
    /// }
    ///
    /// assert_eq!(figure.into_figure().track_count(), 3);
    /// ```
    pub fn done(self) -> Plot<Empty> {
        self.settle()
    }

    /// Adds a track that is already built.
    ///
    /// This is the way in for a track an alternative constructor made, one read
    /// back before it is drawn, or one from outside the crate. Name it before
    /// handing it over: nothing stays pending afterwards, so [`Plot::label`]
    /// and [`Plot::adjust`] do not reach it.
    ///
    /// It is not the way to add the coordinate ruler. An [`AxisTrack`] pushed
    /// through here is in addition to the one appended at the bottom, which
    /// draws two. [`Plot::add_axis`] is the one that knows.
    pub fn add_track(self, track: impl Track + 'static) -> Plot<Empty> {
        let mut plot = self.settle();
        plot.figure = plot.figure.push(track);
        plot
    }

    /// Adds a boxed track, mirroring [`Figure::push_boxed`].
    pub fn add_boxed(self, track: Box<dyn Track>) -> Plot<Empty> {
        let mut plot = self.settle();
        plot.figure = plot.figure.push_boxed(track);
        plot
    }

    /// Coordinate ticks and their labels.
    ///
    /// Adding one explicitly is how it goes anywhere other than the bottom; the
    /// axis is not appended a second time.
    pub fn add_axis(self) -> Plot<AxisTrack> {
        let mut plot = self.settle();
        plot.wants_axis = false;
        plot.park(AxisTrack::new())
    }

    /// Per-read methylation calls at a set of sites, one row per molecule.
    pub fn add_bisulfite(
        self,
        sites: impl Into<Vec<u64>>,
        molecules: impl Into<Vec<Molecule>>,
    ) -> Plot<BisulfiteTrack> {
        self.settle().park(BisulfiteTrack::new(sites, molecules))
    }

    /// Intervals carried by whole clades of a tree, aligned to its leaf rows.
    pub fn add_clades(self, tree: Tree, blocks: impl Into<Vec<CladeBlock>>) -> Plot<CladeTrack> {
        self.settle().park(CladeTrack::new(tree, blocks.into()))
    }

    /// A codon ruler over a coding sequence, numbering residues rather than
    /// bases.
    pub fn add_codons(self, start: u64, end: u64, strand: Strand) -> Plot<CodonTrack> {
        self.settle().park(CodonTrack::new(start, end, strand))
    }

    /// Segmented copy number, on a ladder of whole copies.
    ///
    /// `ploidy` is where balanced sits and there is no default for it: this
    /// crate does not know what it is drawing, and a rule in the wrong place
    /// swaps every gain for a loss.
    pub fn add_copy_number(
        self,
        segments: impl Into<Vec<CopyNumberSegment>>,
        ploidy: f64,
    ) -> Plot<CopyNumberTrack> {
        self.settle()
            .park(CopyNumberTrack::at_ploidy(segments, ploidy))
    }

    /// Read depth, one value per base, starting at the left edge of the region.
    pub fn add_coverage(self, values: impl Into<Vec<f64>>) -> Plot<CoverageTrack> {
        let start = self.figure.region().start();
        self.settle().park(CoverageTrack::new(start, values))
    }

    /// Read depth whose first value is at `start` rather than at the left edge.
    pub fn add_coverage_at(self, start: u64, values: impl Into<Vec<f64>>) -> Plot<CoverageTrack> {
        self.settle().park(CoverageTrack::new(start, values))
    }

    /// Alignment blocks as a dotplot, query against target.
    pub fn add_dotplot(self, blocks: impl Into<Vec<AlignmentBlock>>) -> Plot<DotplotTrack> {
        self.settle().park(DotplotTrack::new(blocks))
    }

    /// Domains, motifs or other interval architectures, one sequence per row.
    pub fn add_domains(self, rows: impl Into<Vec<DomainArchitecture>>) -> Plot<DomainTrack> {
        self.settle().park(DomainTrack::new(rows))
    }

    /// Gene models and other intervals, packed into as few rows as they need.
    pub fn add_features(self, features: impl Into<Vec<Feature>>) -> Plot<FeatureTrack> {
        self.settle().park(FeatureTrack::new(features))
    }

    /// Every sequence of a genome laid end to end, for a whole genome view.
    pub fn add_genome(self, genome: Genome) -> Plot<GenomeTrack> {
        self.settle().park(GenomeTrack::new(genome))
    }

    /// A chromosome drawn as stained bands, with the region marked on it.
    pub fn add_ideogram(self, length: u64, bands: impl Into<Vec<Band>>) -> Plot<IdeogramTrack> {
        self.settle().park(IdeogramTrack::new(length, bands))
    }

    /// A key for the colours used elsewhere in the figure.
    pub fn add_legend(self, legend: Legend) -> Plot<LegendTrack> {
        self.settle().park(LegendTrack::new(legend))
    }

    /// One gene cluster per row, with homology ribbons between them.
    pub fn add_loci(self, loci: impl Into<Vec<Locus>>) -> Plot<LocusTrack> {
        self.settle().park(LocusTrack::new(loci))
    }

    /// A sequence logo, one stack per column, starting at the left edge of the
    /// region.
    pub fn add_logo(self, columns: impl Into<Vec<LogoColumn>>) -> Plot<LogoTrack> {
        let start = self.figure.region().start();
        self.settle().park(LogoTrack::new(start, columns))
    }

    /// A sequence logo whose first column is at `start`.
    pub fn add_logo_at(self, start: u64, columns: impl Into<Vec<LogoColumn>>) -> Plot<LogoTrack> {
        self.settle().park(LogoTrack::new(start, columns))
    }

    /// Association statistics, one point per position.
    pub fn add_manhattan(self, points: impl Into<Vec<Association>>) -> Plot<ManhattanTrack> {
        self.settle().park(ManhattanTrack::new(points))
    }

    /// A value per sample per site, as a grid of coloured cells.
    pub fn add_matrix(
        self,
        sites: impl Into<Vec<u64>>,
        rows: impl Into<Vec<MatrixRow>>,
    ) -> Plot<MatrixTrack> {
        self.settle().park(MatrixTrack::new(sites, rows))
    }

    /// Methylation calls, one row per strand.
    pub fn add_methylation(self, sites: impl Into<Vec<MethylSite>>) -> Plot<MethylationTrack> {
        self.settle().park(MethylationTrack::new(sites))
    }

    /// A multiple sequence alignment, one row per sequence.
    pub fn add_msa(self, sequences: impl Into<Vec<MsaSequence>>) -> Plot<MsaTrack> {
        self.settle().park(MsaTrack::new(sequences))
    }

    /// Open reading frames found in the sequence, three rows per strand,
    /// starting at the left edge of the region.
    pub fn add_orfs(self, seq: impl Into<Vec<u8>>) -> Plot<OrfTrack> {
        let start = self.figure.region().start();
        self.settle().park(OrfTrack::new(start, seq))
    }

    /// Open reading frames in a sequence whose first base is at `start`.
    pub fn add_orfs_at(self, start: u64, seq: impl Into<Vec<u8>>) -> Plot<OrfTrack> {
        self.settle().park(OrfTrack::new(start, seq))
    }

    /// Aligned reads, stacked into rows that do not overlap.
    pub fn add_pileup(self, reads: impl Into<Vec<Read>>) -> Plot<PileupTrack> {
        self.settle().park(PileupTrack::new(reads))
    }

    /// A fitted population, reproductive-number or growth trajectory through
    /// time, with optional uncertainty intervals.
    pub fn add_phylodynamics(
        self,
        points: impl Into<Vec<PhylodynamicPoint>>,
    ) -> Plot<PhylodynamicTrack> {
        self.settle().park(PhylodynamicTrack::new(points))
    }

    /// Site-wise molecular-selection estimates and their independent evidence.
    pub fn add_selection(self, sites: impl Into<Vec<SelectionSite>>) -> Plot<SelectionTrack> {
        self.settle().park(SelectionTrack::new(sites))
    }

    /// Reference bases, starting at the left edge of the region.
    pub fn add_sequence(self, seq: impl Into<Vec<u8>>) -> Plot<SequenceTrack> {
        let start = self.figure.region().start();
        self.settle().park(SequenceTrack::new(start, seq))
    }

    /// Reference bases whose first base is at `start`.
    pub fn add_sequence_at(self, start: u64, seq: impl Into<Vec<u8>>) -> Plot<SequenceTrack> {
        self.settle().park(SequenceTrack::new(start, seq))
    }

    /// Alleles at variable sites, one row per sample.
    pub fn add_snps(
        self,
        names: impl IntoIterator<Item = impl Into<String>>,
        sites: impl Into<Vec<SnpSite>>,
    ) -> Plot<SnpTrack> {
        let names: Vec<String> = names.into_iter().map(Into::into).collect();
        self.settle().park(SnpTrack::new(names, sites))
    }

    /// Reads whose alignment is split across the region, with the joins drawn.
    pub fn add_split_reads(self, reads: impl Into<Vec<SplitRead>>) -> Plot<SplitReadTrack> {
        self.settle().park(SplitReadTrack::new(reads.into()))
    }

    /// Raw signal against sample index, starting at the left edge of the
    /// region.
    pub fn add_squiggle(self, signal: impl Into<Vec<f64>>) -> Plot<SquiggleTrack> {
        let start = self.figure.region().start() as usize;
        self.settle().park(SquiggleTrack::new(start, signal))
    }

    /// Raw signal whose first sample is at index `start`.
    pub fn add_squiggle_at(self, start: usize, signal: impl Into<Vec<f64>>) -> Plot<SquiggleTrack> {
        self.settle().park(SquiggleTrack::new(start, signal))
    }

    /// Structural variants, each drawn as the rearrangement it is.
    pub fn add_structural(
        self,
        variants: impl Into<Vec<StructuralVariant>>,
    ) -> Plot<StructuralTrack> {
        self.settle().park(StructuralTrack::new(variants))
    }

    /// Observed lineage counts and denominators through time.
    pub fn add_surveillance(
        self,
        observations: impl Into<Vec<SurveillanceObservation>>,
    ) -> Plot<SurveillanceTrack> {
        self.settle().park(SurveillanceTrack::new(observations))
    }

    /// Alignment blocks as ribbons between two sequences.
    pub fn add_synteny(self, blocks: impl Into<Vec<AlignmentBlock>>) -> Plot<SyntenyTrack> {
        self.settle().park(SyntenyTrack::new(blocks))
    }

    /// Two trees face to face, with their shared leaves joined.
    pub fn add_tanglegram(self, left: Tree, right: Tree) -> Plot<TanglegramTrack> {
        self.settle().park(TanglegramTrack::new(left, right))
    }

    /// Transcription units, from start site to terminator.
    pub fn add_transcription_units(
        self,
        units: impl Into<Vec<TranscriptionUnit>>,
    ) -> Plot<TranscriptionUnitTrack> {
        self.settle()
            .park(TranscriptionUnitTrack::new(units.into()))
    }

    /// A phylogeny, with its leaves as the rows of the tracks beside it.
    pub fn add_tree(self, tree: Tree) -> Plot<TreeTrack> {
        self.settle().park(TreeTrack::new(tree))
    }

    /// Variant calls, one mark per position.
    pub fn add_variants(self, variants: impl Into<Vec<Variant>>) -> Plot<VariantTrack> {
        self.settle().park(VariantTrack::new(variants))
    }

    /// A statistic computed in windows, as a line, bars or steps.
    pub fn add_windows(self, windows: impl Into<Vec<Window>>) -> Plot<WindowTrack> {
        self.settle().park(WindowTrack::new(windows))
    }
}

impl<T: Named> Plot<T> {
    /// Applies `f` to the track added last.
    ///
    /// The closure is handed the track's own type, so every builder method on
    /// it is available:
    ///
    /// ```
    /// use karyon::{plot, Aggregate, CoverageStyle};
    ///
    /// let svg = plot("chr1:1-1000")
    ///     .unwrap()
    ///     .add_coverage(vec![30.0; 1000])
    ///     .adjust(|track| {
    ///         track
    ///             .aggregate(Aggregate::Min)
    ///             .style(CoverageStyle::Area)
    ///             .height(70.0)
    ///     })
    ///     .to_svg();
    ///
    /// assert!(svg.contains("<svg"));
    /// ```
    ///
    /// Calling it before any track has been added does not compile, because
    /// [`Empty`] is not a track for the closure to take.
    pub fn adjust(mut self, f: impl FnOnce(T) -> T) -> Self {
        self.pending = self.pending.map(f);
        self
    }

    /// Names the track added last, for the left gutter.
    ///
    /// ```
    /// use karyon::plot;
    ///
    /// let svg = plot("chr2:1-1000")
    ///     .unwrap()
    ///     .add_coverage(vec![30.0; 1000])
    ///     .label("depth")
    ///     .to_svg();
    ///
    /// assert!(svg.contains("depth"));
    /// ```
    pub fn label(self, label: impl Into<String>) -> Self {
        let label = label.into();
        self.adjust(|track| track.named(label))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::track::Aggregate;

    fn window() -> Plot {
        plot("NC_000962.3:761001-762000").unwrap()
    }

    #[test]
    fn a_locus_string_sets_the_region() {
        let plot = window();
        assert_eq!(plot.region().seq(), "NC_000962.3");
        assert_eq!(plot.region().start(), 761_000);
        assert_eq!(plot.region().end(), 762_000);
    }

    #[test]
    fn a_bad_locus_string_is_an_error() {
        assert!(plot("not a locus").is_err());
        assert!(plot("NC_000962.3:900-800").is_err());
    }

    #[test]
    fn the_axis_is_added_at_the_bottom_on_its_own() {
        let figure = window().add_coverage(vec![30.0; 1000]).into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn an_explicit_axis_is_not_doubled() {
        let figure = window()
            .add_axis()
            .add_coverage(vec![30.0; 1000])
            .into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn a_removed_axis_stays_removed() {
        let figure = window()
            .add_coverage(vec![30.0; 1000])
            .remove_axis()
            .into_figure();
        assert_eq!(figure.track_count(), 1);
    }

    #[test]
    fn removing_the_axis_before_adding_tracks_also_works() {
        let figure = window()
            .remove_axis()
            .add_coverage(vec![30.0; 1000])
            .into_figure();
        assert_eq!(figure.track_count(), 1);
    }

    #[test]
    fn an_empty_plot_is_still_a_figure() {
        let figure = window().into_figure();
        assert_eq!(figure.track_count(), 1);
        assert!(figure.to_svg().starts_with("<svg"));
    }

    #[test]
    fn every_added_track_reaches_the_figure() {
        let figure = window()
            .add_coverage(vec![30.0; 1000])
            .add_sequence(vec![b'A'; 1000])
            .add_features(vec![Feature::new(761_100, 761_400)])
            .add_variants(vec![Variant::new(761_200)])
            .into_figure();
        assert_eq!(figure.track_count(), 5);
    }

    #[test]
    fn the_label_lands_on_the_track_it_follows() {
        let svg = window()
            .add_coverage(vec![30.0; 1000])
            .label("depth")
            .add_sequence(vec![b'A'; 1000])
            .label("reference")
            .to_svg();
        assert!(svg.contains(">depth<"));
        assert!(svg.contains(">reference<"));
    }

    #[test]
    fn adjust_reaches_the_track_that_was_added_last() {
        // Min and Max disagree wherever the window is wider than the pixels, so
        // an aggregate that did not land shows up as a different drawing.
        let depth: Vec<f64> = (0..1000)
            .map(|i| if i % 7 == 0 { 2.0 } else { 60.0 })
            .collect();
        let plain = window().add_coverage(depth.clone()).to_svg();
        let floored = window()
            .add_coverage(depth)
            .adjust(|track| track.aggregate(Aggregate::Min))
            .to_svg();
        assert_ne!(plain, floored);
    }

    #[test]
    fn adjust_before_a_track_is_added_is_harmless() {
        // Not reachable through `Plot<Empty>`, whose closure has no track to
        // take, but the slot being empty must not panic on the way through.
        let figure = window()
            .add_coverage(vec![30.0; 1000])
            .settle()
            .into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn tracks_that_take_an_array_start_at_the_left_edge() {
        let from_region = window().add_sequence(vec![b'A'; 1000]).to_svg();
        let spelled_out = window().add_sequence_at(761_000, vec![b'A'; 1000]).to_svg();
        assert_eq!(from_region, spelled_out);
    }

    #[test]
    fn the_at_form_can_start_somewhere_else() {
        let left = window().add_coverage(vec![30.0; 500]).to_svg();
        let right = window().add_coverage_at(761_500, vec![30.0; 500]).to_svg();
        assert_ne!(left, right);
    }

    #[test]
    fn a_plot_renders_the_same_figure_as_the_tracks_built_by_hand() {
        let depth: Vec<f64> = (0..1000).map(|i| 40.0 - (i as f64 / 50.0)).collect();
        let genes = vec![Feature::new(761_100, 761_400)
            .name("gene")
            .strand(Strand::Forward)];

        let by_hand = Figure::new(Region::parse("NC_000962.3:761001-762000").unwrap())
            .title("a locus")
            .push(CoverageTrack::new(761_000, depth.clone()).label("depth"))
            .push(FeatureTrack::new(genes.clone()).label("genes"))
            .push(AxisTrack::new())
            .to_svg();

        let by_plot = window()
            .title("a locus")
            .add_coverage(depth)
            .label("depth")
            .add_features(genes)
            .label("genes")
            .to_svg();

        assert_eq!(by_hand, by_plot);
    }

    #[test]
    fn figure_settings_survive_the_chain() {
        let svg = window()
            .width(1200.0)
            .title("wide")
            .theme(Theme::dark())
            .add_coverage(vec![30.0; 1000])
            .to_svg();
        assert!(svg.contains("width=\"1200\""));
        assert!(svg.contains(">wide<"));
    }

    #[test]
    fn the_region_label_can_be_left_out() {
        let with = window().add_coverage(vec![30.0; 1000]).to_svg();
        let without = window()
            .add_coverage(vec![30.0; 1000])
            .remove_region_label()
            .to_svg();
        let locus = Region::parse("NC_000962.3:761001-762000")
            .unwrap()
            .to_string();
        let drawn = |svg: &str| -> bool {
            svg.split("<text").skip(1).any(|piece| {
                piece
                    .split("</text>")
                    .next()
                    .is_some_and(|t| t.contains(&locus))
            })
        };
        assert!(drawn(&with));
        // Only the drawn label goes; the document keeps its accessible name.
        assert!(!drawn(&without));
        assert!(without.contains("<title"));
    }

    #[test]
    fn a_track_from_outside_can_still_be_pushed() {
        let figure = window()
            .add_track(CoverageTrack::new(761_000, vec![30.0; 1000]).label("by hand"))
            .into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn add_track_after_add_does_not_drop_the_earlier_track() {
        // `add_track` settles the pending track first, so the two coexist.
        let figure = window()
            .add_coverage(vec![30.0; 1000])
            .label("depth")
            .add_track(AxisTrack::new())
            .into_figure();
        assert_eq!(figure.track_count(), 3);
    }

    // ------------------------------------------------------------------
    // Audit probes. Written to settle claims about `plot.rs`, not to be
    // kept as they are: the ones that fail are the findings.
    // ------------------------------------------------------------------

    use crate::track::{CigarOp, SplitSegment, Stain, SvKind};

    fn newick() -> Tree {
        Tree::parse_newick("((a:0.1,b:0.2):0.3,c:0.4);").unwrap()
    }

    /// Every `add_` in the module, in one chain. Nothing may fall out.
    #[test]
    fn no_add_method_drops_the_track_before_it() {
        let figure = window()
            .add_bisulfite(
                vec![761_010u64],
                vec![Molecule::new("m1", vec![Some(true)])],
            )
            .add_clades(
                newick(),
                vec![CladeBlock::new(761_010, 761_020, ["a", "b"])],
            )
            .add_codons(761_000, 761_090, Strand::Forward)
            .add_coverage(vec![30.0; 1000])
            .add_coverage_at(761_100, vec![30.0; 100])
            .add_dotplot(vec![AlignmentBlock::new(
                761_000, 761_100, 761_000, 761_100,
            )])
            .add_domains(vec![DomainArchitecture::new("s1", 1_000)])
            .add_features(vec![Feature::new(761_100, 761_400)])
            .add_genome(Genome::new([("chr1", 1_000u64)]))
            .add_ideogram(1_000_000, vec![Band::new(0, 1_000, Stain::Gneg)])
            .add_legend(Legend::new().key("a", "#000000"))
            .add_loci(vec![Locus::new("l1", vec![Feature::new(761_000, 761_100)])])
            .add_logo(vec![LogoColumn::acgt(1.0, 0.0, 0.0, 0.0)])
            .add_logo_at(761_500, vec![LogoColumn::acgt(0.0, 1.0, 0.0, 0.0)])
            .add_manhattan(vec![Association::new(761_100, 3.0)])
            .add_matrix(vec![761_100u64], vec![MatrixRow::new("s1", vec![0.5])])
            .add_methylation(vec![MethylSite::new(761_100, Strand::Forward, 0.5, 10)])
            .add_msa(vec![MsaSequence::new("s1", b"ACGT".to_vec())])
            .add_orfs(vec![b'A'; 1000])
            .add_orfs_at(761_100, vec![b'A'; 100])
            .add_phylodynamics(vec![PhylodynamicPoint::new(761_100, 240.0)])
            .add_pileup(vec![Read::new(761_100, vec![CigarOp::Match(100)])])
            .add_selection(vec![SelectionSite::new(761_100).rates(0.4, 0.2)])
            .add_sequence(vec![b'A'; 1000])
            .add_sequence_at(761_100, vec![b'C'; 100])
            .add_snps(
                vec!["s1".to_string()],
                vec![SnpSite::new(761_100, b'A', vec![b'G'])],
            )
            .add_split_reads(vec![SplitRead::new(vec![SplitSegment::new(
                761_000,
                761_100,
                Strand::Forward,
            )])])
            .add_squiggle(vec![1.0; 100])
            .add_squiggle_at(10, vec![1.0; 100])
            .add_structural(vec![StructuralVariant::new(
                761_100,
                761_200,
                SvKind::Deletion,
            )])
            .add_surveillance(vec![SurveillanceObservation::new(761_100, "L1", 42, 100)])
            .add_synteny(vec![AlignmentBlock::new(
                761_000, 761_100, 761_000, 761_100,
            )])
            .add_tanglegram(newick(), newick())
            .add_transcription_units(vec![TranscriptionUnit::new(
                761_000,
                761_500,
                Strand::Forward,
            )])
            .add_tree(newick())
            .add_variants(vec![Variant::new(761_200)])
            .add_windows(vec![Window::new(761_000, 761_100, 1.0)])
            .add_track(AxisTrack::new())
            .add_axis()
            .into_figure();
        assert_eq!(figure.track_count(), 39);
    }

    #[test]
    fn every_way_a_locus_string_can_be_wrong_is_an_invalid_locus() {
        // `Region::parse` rejects a start of 0 and an end before the start
        // before `Region::new` is ever reached, so `EmptyRegion` cannot come
        // out of here and the `# Errors` section does not promise it.
        for bad in [
            "not a locus",
            "NC_000962.3",
            "NC_000962.3:900",
            "NC_000962.3:900-800",
            "NC_000962.3:0-100",
            "NC_000962.3:a-b",
        ] {
            let err = plot(bad).unwrap_err();
            assert!(
                matches!(err, Error::InvalidLocus { .. }),
                "plot({bad:?}) gave {err:?}"
            );
        }
    }

    #[test]
    fn remove_axis_suppresses_the_automatic_axis_and_only_that() {
        // An axis put somewhere on purpose is not something to take back, and
        // the doc says so. Suppressing the automatic one works from either
        // side of the call that would have added it.
        let put_there = window()
            .add_axis()
            .add_coverage(vec![30.0; 1000])
            .remove_axis()
            .into_figure();
        assert_eq!(put_there.track_count(), 2);

        let before = window().remove_axis().add_coverage(vec![30.0; 1000]);
        let after = window().add_coverage(vec![30.0; 1000]).remove_axis();
        assert_eq!(before.into_figure().track_count(), 1);
        assert_eq!(after.into_figure().track_count(), 1);
    }

    #[test]
    fn add_axis_is_the_one_that_knows_about_the_automatic_axis() {
        // `add_track` sees an `impl Track` and cannot tell a ruler from
        // anything else, so a ruler handed to it is a second ruler. That is
        // what `add_track` documents, and `add_axis` is the way to mean one.
        let handed_over = window()
            .add_coverage(vec![30.0; 1000])
            .add_track(AxisTrack::new())
            .into_figure();
        let asked_for = window()
            .add_coverage(vec![30.0; 1000])
            .add_axis()
            .into_figure();
        assert_eq!(handed_over.track_count(), 3);
        assert_eq!(asked_for.track_count(), 2);
    }

    #[test]
    fn a_track_added_after_a_save_sits_below_the_axis() {
        let path = std::env::temp_dir().join("karyon_plot_after_save.svg");
        let carried_on = window()
            .add_coverage(vec![30.0; 1000])
            .save(&path)
            .unwrap()
            .add_sequence(vec![b'A'; 1000])
            .into_figure()
            .to_svg();
        let _ = std::fs::remove_file(&path);

        let region = || Region::parse("NC_000962.3:761001-762000").unwrap();
        let in_order = Figure::new(region())
            .push(CoverageTrack::new(761_000, vec![30.0; 1000]))
            .push(SequenceTrack::new(761_000, vec![b'A'; 1000]))
            .push(AxisTrack::new())
            .to_svg();
        let axis_in_the_middle = Figure::new(region())
            .push(CoverageTrack::new(761_000, vec![30.0; 1000]))
            .push(AxisTrack::new())
            .push(SequenceTrack::new(761_000, vec![b'A'; 1000]))
            .to_svg();
        // A save closes the stack: the axis is already in it, so a later
        // track goes under the ruler rather than over it, which is what the
        // doc on `save` says.
        assert_eq!(carried_on, axis_in_the_middle);
        assert_ne!(carried_on, in_order);
    }

    #[test]
    fn saving_twice_does_not_draw_two_axes() {
        let first = std::env::temp_dir().join("karyon_plot_save_once.svg");
        let second = std::env::temp_dir().join("karyon_plot_save_twice.svg");
        let figure = window()
            .add_coverage(vec![30.0; 1000])
            .save(&first)
            .unwrap()
            .save(&second)
            .unwrap()
            .into_figure();
        let _ = std::fs::remove_file(&first);
        let _ = std::fs::remove_file(&second);
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn a_stack_can_be_built_in_a_loop() {
        // Every arm of a loop has to have one type, so the pending track has
        // to be put away before the next turn. That is what `done` is for.
        let mut figure = window();
        for name in ["tumour", "normal", "germline"] {
            figure = figure.add_coverage(vec![30.0; 1000]).label(name).done();
        }
        let svg = figure.into_figure().to_svg();
        for name in ["tumour", "normal", "germline"] {
            assert!(svg.contains(&format!(">{name}<")));
        }
    }

    #[test]
    fn a_boxed_track_goes_in_like_any_other() {
        let track: Box<dyn Track> = Box::new(CoverageTrack::new(761_000, vec![30.0; 1000]));
        let figure = window().add_boxed(track).into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    #[test]
    fn sample_names_do_not_have_to_be_owned_strings() {
        let sites = vec![SnpSite::new(761_100, b'A', vec![b'A', b'G'])];
        let figure = window().add_snps(["S01", "S02"], sites).into_figure();
        assert_eq!(figure.track_count(), 2);
    }

    /// Every `add_` whose arguments could be transposed, against the
    /// constructor it stands for.
    #[test]
    fn paired_arguments_reach_the_constructor_in_order() {
        fn bare(track: impl Track + 'static) -> String {
            Figure::new(Region::parse("NC_000962.3:761001-762000").unwrap())
                .push(track)
                .to_svg()
        }

        let sites = vec![761_100u64, 761_200, 761_300];
        let molecules = vec![
            Molecule::new("m1", vec![Some(true), Some(false), None]),
            Molecule::new("m2", vec![None, Some(true), Some(true)]),
        ];
        assert_eq!(
            window()
                .add_bisulfite(sites.clone(), molecules.clone())
                .remove_axis()
                .to_svg(),
            bare(BisulfiteTrack::new(sites.clone(), molecules)),
        );

        let rows = vec![
            MatrixRow::new("s1", vec![0.1, 0.9, 0.5]),
            MatrixRow::new("s2", vec![0.8, 0.2, 0.4]),
        ];
        assert_eq!(
            window()
                .add_matrix(sites.clone(), rows.clone())
                .remove_axis()
                .to_svg(),
            bare(MatrixTrack::new(sites.clone(), rows)),
        );

        let names = vec!["S01".to_string(), "S02".to_string()];
        let snps = vec![
            SnpSite::new(761_100, b'A', vec![b'A', b'G']),
            SnpSite::new(761_300, b'C', vec![b'T', b'C']),
        ];
        assert_eq!(
            window()
                .add_snps(names.clone(), snps.clone())
                .remove_axis()
                .to_svg(),
            bare(SnpTrack::new(names, snps)),
        );

        let bands = vec![
            Band::new(0, 500_000, Stain::Gneg),
            Band::new(500_000, 1_000_000, Stain::Gpos75),
        ];
        assert_eq!(
            window()
                .add_ideogram(1_000_000, bands.clone())
                .remove_axis()
                .to_svg(),
            bare(IdeogramTrack::new(1_000_000, bands)),
        );

        // Both arguments are a `Tree`, so a transposition would compile.
        let left = Tree::parse_newick("((a:0.1,b:0.1):0.2,c:0.3);").unwrap();
        let right = Tree::parse_newick("(a:0.4,(b:0.1,c:0.1):0.2);").unwrap();
        let tangled = window()
            .add_tanglegram(left.clone(), right.clone())
            .remove_axis()
            .to_svg();
        assert_eq!(
            tangled,
            bare(TanglegramTrack::new(left.clone(), right.clone()))
        );
        assert_ne!(tangled, bare(TanglegramTrack::new(right, left)));

        // The codon ruler's start and end are both `u64`, so a transposition
        // would compile. It would also be invisible: `CodonTrack::new` puts a
        // reversed pair back in order, so the two spellings draw the same
        // ruler and only the pass-through can be checked.
        let codons = window()
            .add_codons(761_100, 761_400, Strand::Forward)
            .remove_axis()
            .to_svg();
        assert_eq!(
            codons,
            bare(CodonTrack::new(761_100, 761_400, Strand::Forward))
        );
        assert_eq!(
            bare(CodonTrack::new(761_400, 761_100, Strand::Forward)),
            bare(CodonTrack::new(761_100, 761_400, Strand::Forward))
        );
    }

    /// A squiggle's axis is sample index; the plot casts the region start.
    #[test]
    fn a_squiggle_starts_at_the_left_edge_of_the_region() {
        let over_samples = Plot::over(Region::new("read", 0, 2_000).unwrap());
        let implicit = over_samples.add_squiggle(vec![50.0; 2_000]).to_svg();
        let spelled_out = Plot::over(Region::new("read", 0, 2_000).unwrap())
            .add_squiggle_at(0, vec![50.0; 2_000])
            .to_svg();
        assert_eq!(implicit, spelled_out);

        let shifted = Plot::over(Region::new("read", 500, 2_500).unwrap())
            .add_squiggle(vec![50.0; 2_000])
            .to_svg();
        let shifted_by_hand = Plot::over(Region::new("read", 500, 2_500).unwrap())
            .add_squiggle_at(500, vec![50.0; 2_000])
            .to_svg();
        assert_eq!(shifted, shifted_by_hand);
    }
}
