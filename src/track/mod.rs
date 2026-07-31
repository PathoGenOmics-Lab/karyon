//! The [`Track`] trait and the track types shipped with the crate.
//!
//! A track owns one horizontal band of the figure. It knows how tall it wants
//! to be and how to draw itself, and it never decides where it sits: the
//! [`Figure`](crate::Figure) stacks the bands and hands each track the
//! rectangle it may paint in, already clipped.

pub mod axis;
pub mod coverage;
pub mod feature;
pub mod ideogram;
pub mod logo;
pub mod pileup;
pub mod sequence;
pub mod variant;

pub use axis::AxisTrack;
pub use coverage::{Aggregate, CoverageStyle, CoverageTrack};
pub use feature::{Feature, FeatureTrack, Strand};
pub use ideogram::{Band, IdeogramTrack, Stain};
pub use logo::{Centering, LogoColumn, LogoScore, LogoStack, LogoTrack, StackOrder};
pub use pileup::{CigarOp, PileupLayout, PileupTrack, Read, ReadColoring, Segment};
pub use sequence::SequenceTrack;
pub use variant::{Variant, VariantStyle, VariantTrack};

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
    /// The region on display.
    pub region: &'a Region,
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
