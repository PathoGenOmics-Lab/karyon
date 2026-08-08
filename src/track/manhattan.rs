//! Association statistics along the sequence: the Manhattan plot.
//!
//! One point per test, height by significance, and a line where significance
//! starts. The plot is named for the skyline that appears when a real signal
//! stacks a run of neighbouring markers into a tower, and it is read as a
//! texture rather than as a set of markers.
//!
//! # A tower is only a tower against a crowd
//!
//! Which settles the drawing. Points are small by default, and a point above
//! the threshold is given a ring rather than a larger disc: the ring keeps it
//! findable where the crowd is densest without swallowing the neighbours it has
//! to be seen against.
//!
//! # Where significance starts is not decided here
//!
//! [`ManhattanTrack::threshold`] takes a number from the caller,
//! [`ManhattanTrack::significant`] returns nothing at all until it has one, and
//! no point is coloured as a hit without it. What is being plotted is the
//! caller's too: [`Association::from_p_value`] converts a p-value into the
//! conventional `-log10(p)`, but a test statistic or a Bayes factor plots the
//! same way and only [`ManhattanTrack::unit`] tells the axis which of them it
//! is showing.

use crate::scale::Scale;
use crate::style::{Emphasis, LinePattern, QuantitativeAxis, Symbol};
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// One tested position and what the test said about it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Association {
    /// Position, 0-based.
    pub pos: u64,
    /// The plotted statistic, conventionally `-log10(p)`.
    pub value: f64,
}

impl Association {
    /// A point at a 0-based position.
    pub fn new(pos: u64, value: f64) -> Self {
        Association { pos, value }
    }

    /// A point from a p-value, converted to `-log10(p)`.
    ///
    /// A p-value of zero has no logarithm, so it is floored at the smallest
    /// positive double rather than becoming infinity. Anything at that floor is
    /// already far above any threshold you would draw.
    pub fn from_p_value(pos: u64, p: f64) -> Self {
        let p = if p.is_finite() {
            p.max(f64::MIN_POSITIVE)
        } else {
            1.0
        };
        Association {
            pos,
            value: -p.log10(),
        }
    }
}

/// Points scored by significance, with a threshold line.
///
/// The x axis is genomic, so this draws one chromosome or one region of one.
/// A genome-wide plot that lays every chromosome end to end is a different
/// coordinate system, and this crate does not pretend otherwise: give it a
/// region spanning the sequence you are testing.
///
/// ```
/// use karyon::{Association, Figure, ManhattanTrack, Region};
///
/// let hits = vec![
///     Association::from_p_value(1_000, 0.4),
///     Association::from_p_value(2_000, 3e-9),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 0, 3_000).unwrap())
///     .push(ManhattanTrack::new(hits).genome_wide_threshold().label("GWAS"))
///     .to_svg();
/// assert!(svg.contains("<circle"));
/// ```
#[derive(Debug, Clone)]
pub struct ManhattanTrack {
    points: Vec<Association>,
    label: Option<String>,
    height: f64,
    radius: f64,
    max: Option<f64>,
    threshold: Option<f64>,
    color: Option<String>,
    significant_color: Option<String>,
    unit: String,
    show_scale: bool,
    bands: Vec<u64>,
    axis: QuantitativeAxis,
}

impl ManhattanTrack {
    /// A track holding `points`.
    pub fn new(points: impl Into<Vec<Association>>) -> Self {
        ManhattanTrack {
            points: points.into(),
            label: None,
            height: 90.0,
            radius: 2.2,
            max: None,
            threshold: None,
            color: None,
            significant_color: None,
            unit: String::new(),
            show_scale: true,
            bands: Vec::new(),
            axis: QuantitativeAxis::new(),
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(8.0);
        self
    }

    /// Alternates the point colour at each of these positions.
    ///
    /// For a figure whose axis is several sequences laid end to end: hand it
    /// [`Genome::boundaries`](crate::Genome::boundaries) and the chromosomes
    /// come out in alternating shades, which is the only thing telling a
    /// reader of a genome-wide plot where one ends and the next begins. A
    /// point past the last boundary keeps the last shade.
    pub fn bands(mut self, boundaries: impl Into<Vec<u64>>) -> Self {
        self.bands = boundaries.into();
        self.bands.sort_unstable();
        self
    }

    /// Which band a position falls in, counting from zero.
    pub fn band_of(&self, position: u64) -> usize {
        self.bands
            .iter()
            .rposition(|start| *start <= position)
            .unwrap_or(0)
    }

    /// Sets the radius of a point.
    ///
    /// Points are small on purpose: a Manhattan plot is read as a texture with
    /// towers in it, and a marker big enough to admire individually turns the
    /// texture into a smear.
    pub fn radius(mut self, radius: f64) -> Self {
        self.radius = radius.max(0.3);
        self
    }

    /// Pins the top of the y axis.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self.axis.max = Some(max);
        self
    }

    /// Uses a shared quantitative-axis contract for range, ticks, units and
    /// reference lines.
    pub fn axis(mut self, axis: QuantitativeAxis) -> Self {
        self.max = axis.max;
        self.unit = axis.unit.clone();
        self.axis = axis;
        self
    }

    /// Draws a horizontal line at `value` and colours what rises above it.
    pub fn threshold(mut self, value: f64) -> Self {
        self.threshold = Some(value);
        self
    }

    /// The usual genome-wide threshold, `-log10(5e-8)`, about 7.3.
    ///
    /// It is a Bonferroni correction for a million independent tests, which is
    /// the convention in human GWAS and frequently the wrong number everywhere
    /// else: what it should be follows from how many independent tests were
    /// really run, and a shorter genome or stronger linkage leaves far fewer
    /// than a million. Set your own if you know it.
    pub fn genome_wide_threshold(self) -> Self {
        self.threshold(-(5e-8f64).log10())
    }

    /// Sets the colour of points below the threshold.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets the colour of points above it.
    pub fn significant_color(mut self, color: impl Into<String>) -> Self {
        self.significant_color = Some(color.into());
        self
    }

    /// Sets a unit suffix for the axis labels, such as `" -log10 p"`.
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self.axis.unit = self.unit.clone();
        self
    }

    /// Draws or hides the value axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// The points in the track.
    pub fn points(&self) -> &[Association] {
        &self.points
    }

    /// Points at or above the threshold, in input order.
    ///
    /// Empty when no threshold was set, since without one no point is
    /// significant and calling some of them so would be an invention.
    pub fn significant(&self) -> Vec<Association> {
        let Some(threshold) = self.threshold else {
            return Vec::new();
        };
        self.points
            .iter()
            .filter(|p| p.value.is_finite() && p.value >= threshold)
            .copied()
            .collect()
    }

    /// Top of the axis: the pinned maximum, or the tallest point with a little
    /// headroom, or the threshold when every point is under it.
    fn ceiling(&self) -> f64 {
        if let Some(max) = self.axis.max.or(self.max) {
            return max.max(1e-9);
        }
        let tallest = self
            .points
            .iter()
            .map(|p| p.value)
            .filter(|v| v.is_finite())
            .fold(0.0f64, f64::max);
        let with_threshold = self.threshold.map_or(tallest, |t| tallest.max(t));
        (with_threshold * 1.08).max(1e-9)
    }
}

impl Track for ManhattanTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.points.is_empty() {
            return 0.0;
        }
        let (floor, ceiling) = self.axis.resolve(0.0, self.ceiling());
        let labels = [self.axis.label(floor), self.axis.label(ceiling)];
        labels
            .iter()
            .map(|label| text_width(label, theme.font_size - 1.0))
            .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let baseline = band.bottom();
        ctx.svg.line(
            band.x,
            baseline - 0.5,
            band.right(),
            baseline - 0.5,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );

        if self.points.is_empty() {
            return;
        }
        let (floor, ceiling) = self.axis.resolve(0.0, self.ceiling());
        let y_of =
            |value: f64| baseline - ((value - floor) / (ceiling - floor)).clamp(0.0, 1.0) * band.h;

        for value in self.axis.values(floor, ceiling).into_iter().skip(1) {
            ctx.svg.line(
                band.x,
                y_of(value),
                band.right(),
                y_of(value),
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
        }

        let plain = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.muted.clone());
        let significant = self
            .significant_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());

        if let Some(threshold) = self.threshold {
            if threshold >= floor && threshold <= ceiling {
                let y = y_of(threshold);
                ctx.svg.line_pattern(
                    band.x,
                    y,
                    band.right(),
                    y,
                    &significant,
                    ctx.theme.tokens.stroke,
                    LinePattern::Dashed,
                );
            }
        }

        for reference in &self.axis.references {
            if !reference.value.is_finite() || reference.value < floor || reference.value > ceiling
            {
                continue;
            }
            let style = ctx.theme.mark_style(reference.emphasis);
            let color = if reference.emphasis == Emphasis::Alert {
                &significant
            } else {
                &ctx.theme.muted
            };
            ctx.svg.line_pattern(
                band.x,
                y_of(reference.value),
                band.right(),
                y_of(reference.value),
                color,
                style.stroke_width,
                reference.pattern,
            );
        }

        for point in &self.points {
            if !ctx.region.contains(point.pos) || !point.value.is_finite() {
                continue;
            }
            let above = self.threshold.is_some_and(|t| point.value >= t);
            // Every other sequence a shade lighter, which is what separates one
            // chromosome from the next when the axis is all of them.
            let banded = if self.bands.len() > 1 && self.band_of(point.pos) % 2 == 1 {
                mix(&plain, ctx.theme.surface(), 0.42)
            } else {
                plain.clone()
            };
            let color = if above { &significant } else { &banded };
            let x = ctx.scale.x_center(point.pos);
            let y = y_of(point.value);
            let radius = self.radius * ctx.visual_scale;
            let symbol = if above {
                Symbol::Diamond
            } else {
                ctx.theme.symbol(self.band_of(point.pos))
            };
            if above {
                // A hit is worth a ring, so it stays a point where the texture
                // around it is densest.
                ctx.svg.symbol_ringed(
                    x,
                    y,
                    radius + ctx.theme.tokens.hairline,
                    symbol,
                    color,
                    ctx.theme.surface(),
                    ctx.theme.tokens.hairline,
                );
            } else {
                ctx.svg.symbol(x, y, radius, symbol, color);
            }
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            for value in self.axis.values(floor, ceiling) {
                ctx.svg.text(
                    right,
                    (y_of(value) + size * 0.35)
                        .max(band.y + size)
                        .min(baseline - size * 0.22),
                    &self.axis.label(value),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
            for reference in &self.axis.references {
                let Some(label) = reference.label.as_deref() else {
                    continue;
                };
                if reference.value >= floor && reference.value <= ceiling {
                    ctx.svg.text(
                        band.x + ctx.theme.tokens.label_gap,
                        (y_of(reference.value) - ctx.theme.tokens.row_gap).max(band.y + size),
                        label,
                        &ctx.theme.foreground,
                        size,
                        Anchor::Start,
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 0, 10_000).unwrap()
    }

    #[test]
    fn a_p_value_becomes_minus_log10_of_itself() {
        assert!((Association::from_p_value(0, 0.1).value - 1.0).abs() < 1e-12);
        assert!((Association::from_p_value(0, 1e-8).value - 8.0).abs() < 1e-12);
        assert_eq!(Association::from_p_value(0, 1.0).value, 0.0);
    }

    #[test]
    fn a_p_value_of_zero_stays_finite() {
        let point = Association::from_p_value(0, 0.0);
        assert!(point.value.is_finite());
        assert!(
            point.value > 300.0,
            "and is still enormous: {}",
            point.value
        );
    }

    #[test]
    fn the_genome_wide_threshold_is_the_conventional_one() {
        let track = ManhattanTrack::new(Vec::new()).genome_wide_threshold();
        let threshold = track.threshold.unwrap();
        assert!((threshold - 7.30103).abs() < 1e-5, "got {threshold}");
    }

    #[test]
    fn banding_alternates_at_each_boundary() {
        let track = ManhattanTrack::new(Vec::new()).bands(vec![0u64, 1_000, 1_600]);
        assert_eq!(track.band_of(0), 0);
        assert_eq!(track.band_of(999), 0);
        assert_eq!(track.band_of(1_000), 1);
        assert_eq!(track.band_of(1_599), 1);
        assert_eq!(track.band_of(9_000), 2, "past the last boundary");
        // Boundaries arrive in whatever order they are given.
        assert_eq!(
            ManhattanTrack::new(Vec::new())
                .bands(vec![1_600u64, 0, 1_000])
                .band_of(1_200),
            1
        );
    }

    #[test]
    fn banded_sequences_come_out_in_two_shades() {
        let points = vec![Association::new(10, 3.0), Association::new(1_200, 3.0)];
        let plain = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points.clone()))
            .to_svg();
        let banded = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).bands(vec![0u64, 1_000]))
            .to_svg();
        assert_ne!(plain, banded, "nothing alternated");
        let base = Theme::light().muted;
        assert!(banded.contains(&mix(&base, Theme::light().surface(), 0.42)));
    }

    #[test]
    fn one_boundary_is_not_a_banding() {
        // A single sequence has nothing to alternate with.
        let points = vec![Association::new(10, 3.0)];
        let one = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points.clone()).bands(vec![0u64]))
            .to_svg();
        let none = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points))
            .to_svg();
        assert_eq!(one, none);
    }

    #[test]
    fn nothing_is_significant_without_a_threshold() {
        let points = vec![Association::new(1, 50.0), Association::new(2, 0.1)];
        let track = ManhattanTrack::new(points.clone());
        assert!(track.significant().is_empty());
        assert_eq!(
            ManhattanTrack::new(points)
                .threshold(1.0)
                .significant()
                .len(),
            1
        );
    }

    #[test]
    fn the_ceiling_leaves_headroom_but_a_pin_does_not() {
        let points = vec![Association::new(1, 10.0)];
        let fitted = ManhattanTrack::new(points.clone());
        let pinned = ManhattanTrack::new(points).max(10.0);
        assert!(fitted.ceiling() > 10.0);
        assert_eq!(pinned.ceiling(), 10.0);
    }

    #[test]
    fn a_threshold_above_every_point_still_fits_on_the_axis() {
        let track = ManhattanTrack::new(vec![Association::new(1, 2.0)]).threshold(7.3);
        assert!(track.ceiling() >= 7.3, "the line would be off the top");
    }

    #[test]
    fn an_empty_track_asks_for_no_axis_and_draws_a_baseline() {
        let track = ManhattanTrack::new(Vec::new());
        assert_eq!(track.y_axis_width(&Theme::light()), 0.0);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("<line"));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn hits_and_misses_are_drawn_in_different_colours() {
        let points = vec![Association::new(1_000, 1.0), Association::new(2_000, 9.0)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                ManhattanTrack::new(points)
                    .threshold(5.0)
                    .color("#111111")
                    .significant_color("#ff0000"),
            )
            .to_svg();
        assert!(svg.contains("#111111"), "the quiet point");
        assert!(svg.contains("#ff0000"), "the hit and the threshold line");
    }

    #[test]
    fn non_finite_values_are_skipped_rather_than_drawn_at_zero() {
        let points = vec![
            Association::new(1_000, f64::NAN),
            Association::new(2_000, 4.0),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).show_scale(false))
            .to_svg();
        assert_eq!(svg.matches("<circle").count(), 1);
    }

    #[test]
    fn points_outside_the_region_are_not_drawn() {
        let points = vec![Association::new(50_000, 4.0)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).show_scale(false))
            .to_svg();
        assert!(!svg.contains("<circle"));
    }
}
