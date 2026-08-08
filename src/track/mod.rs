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
pub mod coverage;
pub mod feature;
pub mod genome;
pub mod ideogram;
pub mod legend;
pub mod locus;
pub mod logo;
pub mod manhattan;
pub mod matrix;
pub mod methylation;
pub mod msa;
pub mod orf;
pub mod pileup;
pub mod sequence;
pub mod snp;
pub mod split;
pub mod squiggle;
pub mod structural;
pub mod synteny;
pub mod tanglegram;
pub mod transcript;
pub mod tree;
pub mod variant;
pub mod window;

pub use axis::AxisTrack;
pub use bisulfite::{BisulfiteTrack, Molecule};
pub use clade::{CladeBlock, CladeTrack};
pub use codon::CodonTrack;
pub use coverage::{Aggregate, CoverageStyle, CoverageTrack};
pub use feature::{strand_color, Feature, FeatureTrack, Strand};
pub use genome::GenomeTrack;
pub use ideogram::{Band, IdeogramTrack, Stain};
pub use legend::{Legend, LegendItem, LegendTrack, Marker};
pub use locus::{GeneShape, Homology, Locus, LocusTrack};
pub use logo::{Centering, LogoColumn, LogoScore, LogoStack, LogoTrack, StackOrder};
pub use manhattan::{Association, ManhattanTrack};
pub use matrix::{CellScale, MatrixRow, MatrixTrack};
pub use methylation::{MethylSite, MethylationTrack};
pub use msa::{MsaColoring, MsaDisplay, MsaSequence, MsaTrack, ResidueClass};
pub use orf::{Orf, OrfTrack};
pub use pileup::{CigarOp, PileupLayout, PileupTrack, Read, ReadColoring, Segment};
pub use sequence::SequenceTrack;
pub use snp::{SnpSite, SnpTrack};
pub use split::{SplitRead, SplitReadTrack, SplitSegment};
pub use squiggle::{Move, SquiggleTrack};
pub use structural::{StructuralTrack, StructuralVariant, SvKind};
pub use synteny::{AlignmentBlock, DotplotTrack, SyntenyTrack};
pub use tanglegram::TanglegramTrack;
pub use transcript::{Terminator, TranscriptionUnit, TranscriptionUnitTrack};
pub use tree::{RadialDirection, TraitColumn, TraitScale, TreeProjection, TreeShape, TreeTrack};
pub use variant::{Variant, VariantStyle, VariantTrack};
pub use window::{Window, WindowStyle, WindowTrack};

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::SvgWriter;
use crate::theme::Theme;

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

    /// Draws the track inside `ctx.band`.
    fn draw(&self, ctx: &mut DrawContext<'_>);
}

#[cfg(test)]
mod tests {
    use super::*;

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
