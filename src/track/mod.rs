//! The [`Track`] trait and the track types shipped with the crate.
//!
//! A track owns one horizontal band of the figure. It knows how tall it wants
//! to be and how to draw itself, and it never decides where it sits: the
//! [`Figure`](crate::Figure) stacks the bands and hands each track the
//! rectangle it may paint in, already clipped.
//!
//! # What makes a track a track
//!
//! [`Track::draw`] has to put its marks where [`DrawContext::scale`] says, and
//! that is the entry requirement. A band that ignores the scale and spreads its
//! data evenly across the width it was given is a bar chart that was handed
//! genomic data: right on its own, wrong the moment it is stacked, because a
//! reader is entitled to run a finger down the figure and read every band it
//! crosses at one position.
//!
//! One track breaks that deliberately. [`IdeogramTrack`] draws the whole
//! sequence rather than the region on display, since that is the only way it
//! can answer where the region is: a legitimate exception is not a track that
//! forgot the scale but one whose subject is elsewhere.
//!
//! # Height is asked for, position is not
//!
//! A track answers [`Track::height`] and nothing else about its geometry. The
//! figure asks once, lays the bands out in push order and clips each one, so a
//! track that draws past its edges loses the overflow rather than landing in
//! its neighbour.
//!
//! [`Track::y_axis_width`] is the one request a track cannot answer alone, so
//! the figure settles it once for every band rather than band by band.

pub mod axis;
pub mod bisulfite;
pub mod clade;
pub mod codon;
pub mod copy_number;
pub mod coverage;
pub mod domain;
pub mod dynseq;
pub mod feature;
pub mod genome;
pub mod ideogram;
pub mod junction;
pub mod legend;
pub mod locus;
pub mod logo;
pub mod manhattan;
pub mod matrix;
pub mod methylation;
pub mod msa;
pub mod orf;
pub mod pairs;
pub mod phylodynamics;
pub mod pileup;
pub mod selection;
pub mod sequence;
pub mod snp;
pub mod split;
pub mod squiggle;
pub mod structural;
pub mod surveillance;
pub mod synteny;
pub mod tanglegram;
pub mod traits;
pub mod transcript;
pub mod tree;
pub mod variant;
pub mod window;

pub use axis::AxisTrack;
pub use bisulfite::{BisulfiteTrack, Molecule};
pub use clade::{CladeBlock, CladeTrack};
pub use codon::CodonTrack;
pub use copy_number::{CopyNumber, CopyNumberSegment, CopyNumberTrack};
pub use coverage::{Aggregate, CoverageStyle, CoverageTrack};
pub use domain::{DomainArchitecture, DomainFeature, DomainTrack};
pub use dynseq::DynseqTrack;
pub use feature::{strand_color, Feature, FeatureTrack, Strand};
pub use genome::GenomeTrack;
pub use ideogram::{Band, IdeogramTrack, Stain};
pub use junction::{Junction, JunctionTrack, Motif};
pub use legend::{Legend, LegendItem, LegendTrack, Marker};
pub use locus::{GeneShape, Homology, Locus, LocusTrack};
pub use logo::{Centering, LogoColumn, LogoScore, LogoStack, LogoTrack, StackOrder};
pub use manhattan::{Association, ManhattanTrack};
pub use matrix::{CellScale, MatrixRow, MatrixTrack};
pub use methylation::{MethylSite, MethylationTrack};
pub use msa::{MsaColoring, MsaDisplay, MsaSequence, MsaTrack, ResidueClass};
pub use orf::{Orf, OrfTrack};
pub use pairs::{Pair, PairStyle, PairTrack};
pub use phylodynamics::{PhylodynamicPoint, PhylodynamicScale, PhylodynamicTrack};
pub use pileup::{CigarOp, PileupLayout, PileupTrack, Read, ReadColoring, Segment};
pub use selection::{SelectionEvidence, SelectionSite, SelectionTrack};
pub use sequence::SequenceTrack;
pub use snp::{SnpSite, SnpTrack};
pub use split::{SplitRead, SplitReadTrack, SplitSegment};
pub use squiggle::{Move, SquiggleTrack};
pub use structural::{StructuralTrack, StructuralVariant, SvKind};
pub use surveillance::{
    SurveillanceMetric, SurveillanceObservation, SurveillanceStyle, SurveillanceTrack,
};
pub use synteny::{AlignmentBlock, DotplotTrack, SyntenyTrack};
pub use tanglegram::{TangleLabels, TangleTieStyle, TanglegramTrack};
pub use transcript::{Terminator, TranscriptionUnit, TranscriptionUnitTrack};
pub use tree::{
    AncestralStateLayer, BranchEventLayer, BranchGeometry, BranchIntervalLayer, BranchRateMixture,
    CladeHighlight, HomoplasyLayer, NodeGlyph, NodeGlyphStyle, NodeGlyphTarget, RadialDirection,
    SupportStyle, TraitColumn, TraitScale, TraitStyle, TreeProjection, TreeShape, TreeTrack,
};
pub use variant::{Variant, VariantStyle, VariantTrack};
pub use window::{Window, WindowStyle, WindowTrack};

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::SvgWriter;
use crate::theme::Theme;

/// A line of small print at the bottom right of a band, over a wash of the
/// page so the marks under it do not read through: what the track left out,
/// said where it was left out.
pub(crate) fn band_note(ctx: &mut DrawContext<'_>, text: &str) {
    let band = ctx.band;
    let size = ctx.theme.font_size - 1.0;
    let width = crate::svg::text_width(text, size) + 6.0;
    ctx.svg.rect_opacity(
        band.right() - width,
        band.bottom() - size - 3.0,
        width,
        size + 3.0,
        &ctx.theme.background,
        0.8,
    );
    ctx.svg.text(
        band.right() - 3.0,
        band.bottom() - 2.0,
        text,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::End,
    );
}

/// What a panel says of the tips of its tree that have no row.
pub(crate) fn tips_without_rows(count: usize) -> String {
    match count {
        1 => "1 tip of the tree has no row".to_string(),
        n => format!("{n} tips of the tree have no row"),
    }
}

/// An axis-aligned rectangle in output pixels, y growing downwards.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    /// Left edge.
    pub x: f64,
    /// Top edge.
    pub y: f64,
    /// Width.
    pub w: f64,
    /// Height.
    pub h: f64,
}

impl Rect {
    /// Right edge.
    pub fn right(&self) -> f64 {
        self.x + self.w
    }

    /// Bottom edge, the largest y inside the rectangle.
    pub fn bottom(&self) -> f64 {
        self.y + self.h
    }

    /// Vertical centre.
    pub fn mid_y(&self) -> f64 {
        self.y + self.h / 2.0
    }
}

/// The path of an arc springing from `baseline` at both feet and reaching
/// `apex` in the middle.
///
/// One quadratic, whose control point is put twice as far from the baseline as
/// the apex asked for, since a quadratic reaches only half way to its control.
/// Two tracks draw arcs and they draw the same arc: a structural variant
/// joining its breakpoints and a splice junction spanning its intron differ in
/// what they mean and in how high they go, not in their shape, and two copies
/// of this formula would be two shapes the first time one of them was tuned.
pub(crate) fn arc_path(x0: f64, x1: f64, baseline: f64, apex: f64) -> String {
    let control = baseline - (baseline - apex) * 2.0;
    format!(
        "M{} {}Q{} {} {} {}",
        crate::svg::num(x0),
        crate::svg::num(baseline),
        crate::svg::num((x0 + x1) / 2.0),
        crate::svg::num(control),
        crate::svg::num(x1),
        crate::svg::num(baseline)
    )
}

/// The runs of `points` a line may be drawn through, split wherever a value is
/// missing.
///
/// A line is a claim about every place it crosses, not only about its
/// vertices. Joined across a `None`, it runs straight from the last value
/// before the hole to the first one after, and draws a value at every point in
/// between where the data holds none. So a line stops at a missing value and
/// starts again after it, and the tracks that draw one over values that can be
/// missing split it here rather than each deciding for itself where a line may
/// go.
pub(crate) fn unbroken<T>(points: impl IntoIterator<Item = Option<T>>) -> Vec<Vec<T>> {
    let mut runs = Vec::new();
    let mut run = Vec::new();
    for point in points {
        match point {
            Some(point) => run.push(point),
            None if !run.is_empty() => runs.push(std::mem::take(&mut run)),
            None => {}
        }
    }
    if !run.is_empty() {
        runs.push(run);
    }
    runs
}

/// Everything a track needs in order to draw one band.
pub struct DrawContext<'a> {
    /// Where to write the SVG elements.
    pub svg: &'a mut SvgWriter,
    /// Shared horizontal mapping, identical for every track in the figure.
    pub scale: &'a Scale,
    /// Shared colours and fonts.
    pub theme: &'a Theme,
    /// The band this track may paint in. Output is clipped to it.
    pub band: Rect,
    /// The strip immediately left of the band, as wide as this track asked for
    /// in [`Track::y_axis_width`].
    ///
    /// Somewhere to put a value axis. It is zero width unless the track asked,
    /// and the clip covers it along with the band, so a track that wants no
    /// axis cannot accidentally draw outside itself.
    pub axis: Rect,
    /// The strip immediately right of the band, as wide as this track asked
    /// for in [`Track::right_axis_width`], for a second value axis. Zero
    /// width unless the track asked, and inside the clip when it did.
    pub right_axis: Rect,
    /// The region on display.
    pub region: &'a Region,
    /// Scale applied to track-internal marks and repeated row geometry.
    ///
    /// Use [`DrawContext::px`] for fixed pixel measurements owned by a custom
    /// track. Genomic x coordinates still come exclusively from `scale`.
    pub visual_scale: f64,
}

impl DrawContext<'_> {
    /// Scales one track-internal pixel measurement for the active profile and
    /// density without changing genomic coordinates.
    pub fn px(&self, value: f64) -> f64 {
        value * self.visual_scale
    }
}

/// One horizontal band of a figure.
///
/// Implement this to add a track type the crate does not ship. The contract is
/// small: report a height for a given scale, then draw inside the band you are
/// given.
///
/// ```
/// use karyon::{DrawContext, Figure, Rect, Region, Scale, Track};
///
/// /// A band of flat colour, the smallest useful track.
/// struct Rule(f64);
///
/// impl Track for Rule {
///     fn height(&self, _scale: &Scale) -> f64 {
///         self.0
///     }
///
///     fn draw(&self, ctx: &mut DrawContext<'_>) {
///         let band = ctx.band;
///         ctx.svg.rect(band.x, band.mid_y(), band.w, 1.0, &ctx.theme.rule);
///     }
/// }
///
/// let svg = Figure::new(Region::parse("chr1:1-1000").unwrap())
///     .push(Rule(20.0))
///     .to_svg();
/// assert!(svg.contains("<rect"));
/// ```
pub trait Track {
    /// Height of the band in pixels.
    ///
    /// The scale is passed in because some tracks change height with zoom: a
    /// [`FeatureTrack`] needs more rows when features overlap on screen.
    fn height(&self, scale: &Scale) -> f64;

    /// Text for the left gutter, or `None` for an unlabelled track.
    ///
    /// The figure reserves the gutter only when at least one track wants it.
    fn label(&self) -> Option<&str> {
        None
    }

    /// What this track is, as a sentence says it: `a ruler`, `a phylogeny`,
    /// `aligned reads`.
    ///
    /// The figure's description names each track by its label, and a track
    /// with none by this. It counted the ruler and the key among the tracks
    /// and named neither, so four tracks were three names long.
    fn noun(&self) -> &str {
        "a track"
    }

    /// What the value axis measures, written under the name in the gutter,
    /// such as `-log10 p` or `AF`, or `None` for an axis that needs no title.
    ///
    /// A word set after the top tick as a unit reads as part of that number,
    /// and `10 -log10 p` reads as ten minus something. Under the track's name
    /// it is the axis's title, and the ticks keep their numbers to themselves.
    /// A symbol that is a unit, the `x` of `50x`, stays on the tick.
    fn axis_title(&self) -> Option<&str> {
        None
    }

    /// The key a reader needs for this track as drawn over `region` at
    /// `px_per_bp` pixels a base, where colour alone carries something there,
    /// for a caller who builds one: [`Figure::key`](crate::Figure::key) asks
    /// every track at the zoom it draws at.
    ///
    /// A base drawn as a block of colour is a base whose letter is not on the
    /// page, and its colour names it only to a reader who knows the palette.
    fn key(&self, _region: &Region, _px_per_bp: f64, _theme: &Theme) -> Option<legend::Legend> {
        None
    }

    /// How much room this track wants for a value axis, in pixels.
    ///
    /// A track that returns more than zero gets [`DrawContext::axis`], a strip
    /// between the labels and the plotting area, and may draw its own ticks
    /// there. The figure reserves the widest request across every track, so the
    /// plotting areas still line up.
    ///
    /// Only quantitative tracks want this. A pileup or a sequence has no value
    /// to put a number on, and the default of zero says so.
    fn y_axis_width(&self, _theme: &Theme) -> f64 {
        0.0
    }

    /// How much room this track wants right of its band, for a second value
    /// axis, in pixels.
    ///
    /// A track that draws two quantities over one band, as a scan with the
    /// recombination rate laid over it, reads the second off a scale on the
    /// right, and gets [`DrawContext::right_axis`] for it. The figure takes
    /// the widest request across every track from every plotting area, as it
    /// does on the left, so the bands still end at one x. The default of zero
    /// asks for nothing and moves nothing.
    fn right_axis_width(&self, _theme: &Theme) -> f64 {
        0.0
    }

    /// Whether this track is drawn against the shared coordinate axis.
    ///
    /// Nearly every track is, which is what makes a figure a figure: the ruler
    /// along the bottom measures the window, and each track lays its data on
    /// that same window. So the default is yes, and a track that is not has to
    /// say so.
    ///
    /// A phylogeny is not. Its x is a branch length or a depth in the tree, and
    /// the window it is handed exists only because a figure needs one, so the
    /// ruler under it was measuring nothing: `tree:1-1` drew a rule the width
    /// of the figure with a single tick labelled `1`. A plot appends its axis
    /// only where something is measured against it.
    fn on_coordinates(&self) -> bool {
        true
    }

    /// Whether this track shows where the window is, so that a figure holding
    /// it prints the locus at its top right.
    ///
    /// Everything on the coordinates does, and that is the default. An
    /// ideogram is not on them, since it draws the whole chromosome, and it
    /// marks the window on it, so it says yes too. A phylogeny says neither:
    /// a locus above one names a window the tree is not drawn in.
    fn shows_region(&self) -> bool {
        self.on_coordinates()
    }

    /// Draws the track inside `ctx.band`.
    fn draw(&self, ctx: &mut DrawContext<'_>);
}

/// The vertices of every `<polyline>` in a document, one list per element, in
/// the order they were drawn.
#[cfg(test)]
pub(crate) fn polylines(svg: &str) -> Vec<Vec<(f64, f64)>> {
    svg.split("<polyline points=\"")
        .skip(1)
        .filter_map(|rest| rest.split('"').next())
        .map(|list| {
            list.split(' ')
                .filter_map(|pair| {
                    let (x, y) = pair.split_once(',')?;
                    Some((x.parse().ok()?, y.parse().ok()?))
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_line_is_split_at_every_missing_value_and_nowhere_else() {
        let runs = unbroken([None, Some(1), Some(2), None, None, Some(3), None]);
        assert_eq!(runs, vec![vec![1, 2], vec![3]]);
        assert_eq!(unbroken([Some(1), Some(2)]), vec![vec![1, 2]]);
        assert!(unbroken::<u8>([None, None]).is_empty());
    }

    #[test]
    fn rect_edges_are_consistent() {
        let r = Rect {
            x: 10.0,
            y: 20.0,
            w: 100.0,
            h: 40.0,
        };
        assert_eq!(r.right(), 110.0);
        assert_eq!(r.bottom(), 60.0);
        assert_eq!(r.mid_y(), 40.0);
    }
}
