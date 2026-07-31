//! Mapping between genomic coordinates and horizontal pixels.

use crate::region::Region;

/// Linear map from base positions to x coordinates in the output image.
///
/// Every track in a figure shares one scale, which is what keeps their x axes
/// aligned. A base at 0-based position `p` occupies the pixel span
/// `[x(p), x(p + 1))`, so [`Scale::x`] returns the *left edge* of a base, not
/// its centre. Use [`Scale::x_center`] when you want the middle of a base, as a
/// variant marker does.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scale {
    start: f64,
    span: f64,
    x0: f64,
    width: f64,
}

impl Scale {
    /// Builds the scale that maps `region` onto `width` pixels starting at `x0`.
    pub fn new(region: &Region, x0: f64, width: f64) -> Self {
        Scale {
            start: region.start() as f64,
            span: region.len() as f64,
            x0,
            width,
        }
    }

    /// Left edge of the base at 0-based position `pos`.
    ///
    /// Positions outside the region map outside the plotting area; tracks are
    /// clipped to their band, so nothing spills into a neighbour.
    pub fn x(&self, pos: u64) -> f64 {
        self.x_at(pos as f64)
    }

    /// Centre of the base at 0-based position `pos`.
    pub fn x_center(&self, pos: u64) -> f64 {
        self.x_at(pos as f64 + 0.5)
    }

    /// Left edge of a fractional position, for sub-base precision.
    pub fn x_at(&self, pos: f64) -> f64 {
        self.x0 + (pos - self.start) / self.span * self.width
    }

    /// Inverse of [`Scale::x_at`]: the fractional base position under a pixel.
    pub fn pos_at_x(&self, x: f64) -> f64 {
        self.start + (x - self.x0) / self.width * self.span
    }

    /// Bases covered by one pixel. Above 1 the data must be binned to be drawn.
    pub fn bp_per_px(&self) -> f64 {
        self.span / self.width
    }

    /// Pixels covered by one base. Above roughly 8 there is room for a letter.
    pub fn px_per_bp(&self) -> f64 {
        self.width / self.span
    }

    /// The region on display, 0-based half-open.
    ///
    /// A track is handed the scale before it is handed a region, and some of
    /// them have to know what is on screen in order to say how tall they are:
    /// a pileup packs only the reads in view, so its height follows the view.
    pub fn bounds(&self) -> (u64, u64) {
        let start = self.start.max(0.0) as u64;
        (start, start + self.span.max(0.0) as u64)
    }

    /// Left edge of the plotting area.
    pub fn x0(&self) -> f64 {
        self.x0
    }

    /// Width of the plotting area in pixels.
    pub fn width(&self) -> f64 {
        self.width
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scale() -> Scale {
        let region = Region::new("chr1", 100, 200).unwrap();
        Scale::new(&region, 50.0, 500.0)
    }

    #[test]
    fn region_edges_map_to_plot_edges() {
        let s = scale();
        assert!((s.x(100) - 50.0).abs() < 1e-9);
        assert!((s.x(200) - 550.0).abs() < 1e-9);
    }

    #[test]
    fn a_base_spans_one_pixel_column_width() {
        let s = scale();
        assert!((s.x(101) - s.x(100) - 5.0).abs() < 1e-9);
        assert!((s.x_center(100) - 52.5).abs() < 1e-9);
    }

    #[test]
    fn pos_at_x_inverts_x_at() {
        let s = scale();
        for pos in [100.0, 137.5, 199.0] {
            let round_trip = s.pos_at_x(s.x_at(pos));
            assert!((round_trip - pos).abs() < 1e-9, "{pos}");
        }
    }

    #[test]
    fn the_region_can_be_read_back_off_the_scale() {
        assert_eq!(scale().bounds(), (100, 200));
        let whole = Region::new("NC_000962.3", 0, 4_411_532).unwrap();
        assert_eq!(Scale::new(&whole, 0.0, 900.0).bounds(), (0, 4_411_532));
    }

    #[test]
    fn resolution_helpers_agree() {
        let s = scale();
        assert!((s.bp_per_px() - 0.2).abs() < 1e-9);
        assert!((s.px_per_bp() - 5.0).abs() < 1e-9);
    }
}
