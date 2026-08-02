//! Per-base quantitative signal: read depth, GC content, mappability.
//!
//! Dense values from a start position are what [`CoverageTrack::new`] takes;
//! scattered points come in through [`CoverageTrack::from_pairs`]. The profile
//! is drawn upwards from the floor of the band, so the quantity has to be one
//! whose zero is the bottom in fact and not by convention: a signed or centred
//! statistic belongs in [`WindowTrack`](crate::WindowTrack).
//!
//! # What a pixel column throws away
//!
//! A region wider than the figure has more bases than columns, so each column
//! is a summary of the bases beneath it. Which summary is [`Aggregate`], and
//! the default is [`Aggregate::Max`], so a single base spike survives being
//! drawn at a megabase.
//!
//! The cost is the mirror image, and it is quiet: a maximum cannot show a hole.
//! Across a 5 Mb view one pixel spans thousands of bases, and a 300 bp deletion
//! inside an otherwise well covered gene contributes nothing to the maximum of
//! its column, so the profile is drawn flat over it. [`Aggregate::Min`] takes
//! the floor of each column instead, and the dropout appears. A column with
//! nothing under it at all is skipped rather than drawn at zero, so a gap in
//! the input stays a gap in the profile.
//!
//! # The ceiling moves unless it is pinned
//!
//! Left alone, the top of the band is the largest value on screen, so the
//! profile fills the band whatever the depth is. [`CoverageTrack::max`] fixes
//! the top and is taken literally; an automatic ceiling is lifted by six per
//! cent instead, so the tallest point reads as a peak rather than as something
//! that ran out of room. [`CoverageTrack::log_scale`] is the other way to make
//! a wide range fit into one band.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Track};

/// How a coverage track is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageStyle {
    /// Filled area under the profile. The default, and the one that survives
    /// being printed small.
    Area,
    /// Outline only, for overlaying several samples.
    Line,
    /// One bar per pixel column, closest to what a genome browser shows.
    Bars,
}

/// How values are combined when one pixel covers many bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aggregate {
    /// Highest value in the column. Keeps narrow spikes visible, which is
    /// usually what you want when hunting for duplications or dropouts.
    Max,
    /// Mean of the column. Truer to the overall level, hides single-base
    /// spikes.
    Mean,
    /// Lowest value in the column. The one to use when looking for dropouts
    /// that a mean would smooth away.
    Min,
}

/// A quantitative signal sampled once per base.
///
/// Values are stored densely from a start position, which is the shape
/// `samtools depth` output arrives in. When a pixel covers more than one base
/// the column is reduced with [`Aggregate`]; the SVG therefore has at most one
/// point per pixel however wide the region is.
///
/// ```
/// use karyon::{CoverageTrack, Figure, Region};
///
/// let depth: Vec<f64> = (0..1000).map(|i| 30.0 + (i % 17) as f64).collect();
/// let svg = Figure::new(Region::parse("chr1:1-1000").unwrap())
///     .push(CoverageTrack::new(0, depth).label("depth"))
///     .to_svg();
/// assert!(svg.contains("depth"));
/// ```
#[derive(Debug, Clone)]
pub struct CoverageTrack {
    start: u64,
    values: Vec<f64>,
    label: Option<String>,
    height: f64,
    color: Option<String>,
    style: CoverageStyle,
    aggregate: Aggregate,
    max: Option<f64>,
    log_scale: bool,
    fill_opacity: Option<f64>,
    show_max: bool,
}

impl CoverageTrack {
    /// A track whose `values[i]` describes base `start + i`, 0-based.
    ///
    /// Values need not cover the whole region: anything outside is not
    /// drawn, and non-finite values are treated as missing.
    pub fn new(start: u64, values: impl Into<Vec<f64>>) -> Self {
        CoverageTrack {
            start,
            values: values.into(),
            label: None,
            height: 60.0,
            color: None,
            style: CoverageStyle::Area,
            aggregate: Aggregate::Max,
            max: None,
            log_scale: false,
            fill_opacity: None,
            show_max: true,
        }
    }

    /// A track built from sparse `(position, value)` pairs.
    ///
    /// Positions are 0-based. The dense buffer spans `region`, so memory is
    /// proportional to the region on display, not to the genome. Pairs outside
    /// the region are ignored, and positions not listed stay at zero.
    pub fn from_pairs(region: &Region, pairs: impl IntoIterator<Item = (u64, f64)>) -> Self {
        let mut values = vec![0.0; region.len() as usize];
        for (pos, value) in pairs {
            if region.contains(pos) {
                values[(pos - region.start()) as usize] = value;
            }
        }
        CoverageTrack::new(region.start(), values)
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(1.0);
        self
    }

    /// Overrides the colour, which defaults to the theme accent.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Chooses area, line or bars.
    pub fn style(mut self, style: CoverageStyle) -> Self {
        self.style = style;
        self
    }

    /// Chooses how a pixel column summarises the bases under it.
    pub fn aggregate(mut self, aggregate: Aggregate) -> Self {
        self.aggregate = aggregate;
        self
    }

    /// Pins the top of the y axis, instead of taking the visible maximum.
    ///
    /// Pin it when comparing samples side by side, or the eye will read two
    /// different scales as the same one.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self
    }

    /// Plots `log10(1 + value)` so that a 1000x depth range fits in one band.
    ///
    /// The axis label still reports the value in the original units.
    pub fn log_scale(mut self, enabled: bool) -> Self {
        self.log_scale = enabled;
        self
    }

    /// Sets the fill opacity, between 0 and 1.
    ///
    /// Left alone, each style picks what suits it: an area is a wash under a
    /// drawn line, so it fills at a fifth and lets whatever is behind it show
    /// through, while bars are the mark itself and stay solid.
    pub fn fill_opacity(mut self, opacity: f64) -> Self {
        self.fill_opacity = Some(opacity.clamp(0.0, 1.0));
        self
    }

    /// Shows or hides the small maximum-value annotation.
    pub fn show_max(mut self, show: bool) -> Self {
        self.show_max = show;
        self
    }

    /// Largest finite value inside `region`, or `None` when nothing overlaps.
    pub fn visible_max(&self, region: &Region) -> Option<f64> {
        let (lo, hi) = self.visible_slice(region)?;
        self.values[lo..hi]
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            })
    }

    /// Index range of `self.values` overlapping `region`.
    fn visible_slice(&self, region: &Region) -> Option<(usize, usize)> {
        let end = self.start + self.values.len() as u64;
        if region.end() <= self.start || region.start() >= end {
            return None;
        }
        let lo = region.start().max(self.start) - self.start;
        let hi = region.end().min(end) - self.start;
        if hi <= lo {
            return None;
        }
        Some((lo as usize, hi as usize))
    }

    /// Reduces the bases in `[lo, hi)` to the single value a pixel shows.
    fn sample(&self, lo: f64, hi: f64) -> Option<f64> {
        let origin = self.start as f64;
        let first = (lo - origin).floor().max(0.0);
        let last = (hi - origin).ceil().min(self.values.len() as f64);
        if last <= first {
            return None;
        }
        let slice = &self.values[first as usize..last as usize];
        let mut count = 0usize;
        let mut acc = match self.aggregate {
            Aggregate::Max => f64::NEG_INFINITY,
            Aggregate::Min => f64::INFINITY,
            Aggregate::Mean => 0.0,
        };
        for value in slice.iter().copied().filter(|v| v.is_finite()) {
            count += 1;
            acc = match self.aggregate {
                Aggregate::Max => acc.max(value),
                Aggregate::Min => acc.min(value),
                Aggregate::Mean => acc + value,
            };
        }
        if count == 0 {
            return None;
        }
        Some(match self.aggregate {
            Aggregate::Mean => acc / count as f64,
            _ => acc,
        })
    }

    fn transform(&self, value: f64) -> f64 {
        if self.log_scale {
            (1.0 + value.max(0.0)).log10()
        } else {
            value
        }
    }
}

impl Track for CoverageTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_max {
            return 0.0;
        }
        // Room for the widest label this track could print, which is the
        // ceiling rather than the zero underneath it.
        let widest = self
            .max
            .map(format_value)
            .unwrap_or_else(|| "999.9k".to_string());
        text_width(&widest, theme.font_size - 1.0) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());

        let baseline = band.bottom();
        // Half a stroke inside the band. A one pixel rule centred on the clip
        // edge keeps only half its ink and renders at half the weight of the
        // ceiling rule drawn with identical parameters, which AxisTrack
        // already worked around with the same half pixel.
        ctx.svg.line(
            band.x,
            baseline - 0.5,
            band.right(),
            baseline - 0.5,
            &ctx.theme.rule,
            1.0,
        );

        let ceiling = self
            .max
            .or_else(|| self.visible_max(ctx.region))
            .filter(|m| m.is_finite() && *m > 0.0);
        let Some(ceiling) = ceiling else {
            return;
        };
        // A little headroom, so the tallest point is a peak rather than
        // something that ran out of band. A pinned maximum is taken literally,
        // because that is the whole reason for pinning one.
        let headroom = if self.max.is_some() { 1.0 } else { 1.06 };
        let top = self.transform(ceiling) * headroom;
        if top <= 0.0 {
            return;
        }

        let columns = band.w.max(1.0).ceil() as usize;
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(columns);
        for column in 0..columns {
            let x = band.x + column as f64;
            let lo = ctx.scale.pos_at_x(x);
            let hi = ctx.scale.pos_at_x(x + 1.0);
            let Some(value) = self.sample(lo, hi) else {
                continue;
            };
            let scaled = (self.transform(value) / top).clamp(0.0, 1.0);
            let y = baseline - scaled * band.h;
            match self.style {
                CoverageStyle::Bars => {
                    ctx.svg.rect_opacity(
                        x,
                        y,
                        1.0,
                        baseline - y,
                        &color,
                        self.fill_opacity.unwrap_or(1.0),
                    );
                }
                _ => points.push((x, y)),
            }
        }

        match self.style {
            CoverageStyle::Bars => {}
            CoverageStyle::Line => ctx.svg.polyline(&points, &color, 2.0),
            CoverageStyle::Area => {
                if points.len() >= 2 {
                    let mut d = String::with_capacity(points.len() * 14);
                    d.push('M');
                    d.push_str(&num(points[0].0));
                    d.push(' ');
                    d.push_str(&num(baseline));
                    for (x, y) in &points {
                        d.push_str(" L");
                        d.push_str(&num(*x));
                        d.push(' ');
                        d.push_str(&num(*y));
                    }
                    let last = points[points.len() - 1].0;
                    d.push_str(" L");
                    d.push_str(&num(last));
                    d.push(' ');
                    d.push_str(&num(baseline));
                    d.push('Z');
                    // A wash under a drawn line, rather than a saturated block.
                    // The line is what carries the shape; the fill only says
                    // which side of it is under the curve.
                    ctx.svg.path(&d, &color, self.fill_opacity.unwrap_or(0.18));
                    ctx.svg.polyline(&points, &color, 2.0);
                }
            }
        }

        if self.show_max {
            self.draw_axis(ctx, ceiling, top);
        }
    }
}

impl CoverageTrack {
    /// Draws the value axis in the strip the figure reserved for it.
    ///
    /// Two ticks, zero and the ceiling, plus a hairline across the plot at the
    /// top of the scale. Two is enough: a coverage track is read for its shape
    /// and its order of magnitude, and a ladder of six gridlines would be more
    /// ink than the profile it is measuring.
    fn draw_axis(&self, ctx: &mut DrawContext<'_>, ceiling: f64, top: f64) {
        let band = ctx.band;
        let size = ctx.theme.font_size - 1.0;
        let baseline = band.bottom();
        // Where the ceiling lands once the headroom is accounted for.
        let ceiling_y = baseline - (self.transform(ceiling) / top) * band.h;

        ctx.svg.line(
            band.x,
            ceiling_y,
            band.right(),
            ceiling_y,
            &ctx.theme.rule,
            1.0,
        );

        if ctx.axis.w <= 0.0 {
            return;
        }
        let right = ctx.axis.right() - 4.0;
        let suffix = if self.log_scale { " log" } else { "" };
        ctx.svg.text(
            right,
            (ceiling_y + size * 0.35).max(band.y + size * 0.78),
            &format!("{}{}", format_value(ceiling), suffix),
            &ctx.theme.muted,
            size,
            Anchor::End,
        );
        ctx.svg
            .text(right, baseline, "0", &ctx.theme.muted, size, Anchor::End);
    }
}

/// Compact label for the top of the y axis.
///
/// A depth of 72.04 is labelled `72`, not `72.0`: the extra digit is noise from
/// however the values were computed, and it reads as precision the data does
/// not have.
fn format_value(value: f64) -> String {
    if value >= 1_000_000.0 {
        format!("{}M", trim(value / 1e6))
    } else if value >= 1_000.0 {
        format!("{}k", trim(value / 1e3))
    } else {
        trim(value)
    }
}

/// One decimal at most, and none at all when it would be a zero.
fn trim(value: f64) -> String {
    let rounded = (value * 10.0).round() / 10.0;
    if rounded == rounded.trunc() {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded:.1}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> Region {
        Region::new("chr1", 0, 100).unwrap()
    }

    #[test]
    fn sample_takes_the_maximum_of_a_binned_column() {
        let track = CoverageTrack::new(0, vec![1.0, 9.0, 2.0, 3.0]);
        assert_eq!(track.sample(0.0, 4.0), Some(9.0));
    }

    #[test]
    fn sample_can_average_or_take_the_minimum() {
        let values = vec![1.0, 9.0, 2.0, 4.0];
        let mean = CoverageTrack::new(0, values.clone()).aggregate(Aggregate::Mean);
        assert_eq!(mean.sample(0.0, 4.0), Some(4.0));
        let min = CoverageTrack::new(0, values).aggregate(Aggregate::Min);
        assert_eq!(min.sample(0.0, 4.0), Some(1.0));
    }

    #[test]
    fn sample_returns_one_base_when_zoomed_past_single_base_resolution() {
        let track = CoverageTrack::new(10, vec![5.0, 6.0, 7.0]);
        assert_eq!(track.sample(11.2, 11.6), Some(6.0));
    }

    #[test]
    fn sample_outside_the_data_is_missing_rather_than_zero() {
        let track = CoverageTrack::new(10, vec![5.0, 6.0]);
        assert_eq!(track.sample(0.0, 5.0), None);
        assert_eq!(track.sample(50.0, 60.0), None);
    }

    #[test]
    fn non_finite_values_are_ignored_not_propagated() {
        let track = CoverageTrack::new(0, vec![f64::NAN, 4.0, f64::INFINITY]);
        assert_eq!(track.sample(0.0, 3.0), Some(4.0));
        assert_eq!(track.visible_max(&region()), Some(4.0));
    }

    #[test]
    fn a_fully_missing_column_has_no_maximum() {
        let track = CoverageTrack::new(0, vec![f64::NAN, f64::NAN]);
        assert_eq!(track.sample(0.0, 2.0), None);
        assert_eq!(track.visible_max(&region()), None);
    }

    #[test]
    fn visible_max_only_looks_inside_the_region() {
        let track = CoverageTrack::new(0, vec![1.0, 2.0, 100.0, 3.0]);
        let window = Region::new("chr1", 0, 2).unwrap();
        assert_eq!(track.visible_max(&window), Some(2.0));
    }

    #[test]
    fn data_disjoint_from_the_region_has_no_visible_slice() {
        let track = CoverageTrack::new(1000, vec![1.0, 2.0]);
        assert_eq!(track.visible_slice(&region()), None);
        assert_eq!(track.visible_max(&region()), None);
    }

    #[test]
    fn from_pairs_places_values_at_their_positions() {
        let window = Region::new("chr1", 100, 105).unwrap();
        let track = CoverageTrack::from_pairs(&window, [(101, 7.0), (104, 3.0), (900, 99.0)]);
        assert_eq!(track.values, vec![0.0, 7.0, 0.0, 0.0, 3.0]);
        assert_eq!(track.visible_max(&window), Some(7.0));
    }

    #[test]
    fn log_scale_compresses_the_range() {
        let track = CoverageTrack::new(0, vec![0.0]).log_scale(true);
        assert!((track.transform(9.0) - 1.0).abs() < 1e-12);
        assert_eq!(track.transform(0.0), 0.0);
    }

    #[test]
    fn the_axis_labels_both_ends_of_the_scale() {
        use crate::figure::Figure;
        let depth: Vec<f64> = (0..500).map(|i| (i % 87) as f64).collect();
        let svg = Figure::new(Region::parse("chr1:1-500").unwrap())
            .show_region_label(false)
            .push(CoverageTrack::new(0, depth).label("depth"))
            .to_svg();
        assert!(svg.contains(">86</text>"), "the ceiling should be labelled");
        assert!(svg.contains(">0</text>"), "and so should the floor");
    }

    #[test]
    fn a_pinned_maximum_is_what_the_axis_says() {
        use crate::figure::Figure;
        let svg = Figure::new(Region::parse("chr1:1-500").unwrap())
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 500]).max(250.0).label("d"))
            .to_svg();
        assert!(svg.contains(">250</text>"));
        assert!(!svg.contains(">10</text>"));
    }

    #[test]
    fn value_labels_shorten_large_numbers() {
        assert_eq!(format_value(42.0), "42");
        assert_eq!(format_value(1500.0), "1.5k");
        assert_eq!(format_value(2_500_000.0), "2.5M");
        assert_eq!(format_value(0.5), "0.5");
    }

    #[test]
    fn value_labels_drop_a_decimal_that_rounds_to_zero() {
        assert_eq!(format_value(72.04), "72");
        assert_eq!(format_value(72.44), "72.4");
        assert_eq!(format_value(2_000_400.0), "2M");
    }
}
