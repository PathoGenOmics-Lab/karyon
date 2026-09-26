//! Time-varying phylodynamic estimates with uncertainty kept visible.
//!
//! A skyline, effective population size, reproductive number or lineage
//! growth trajectory is an estimate with a time coordinate and often an
//! interval. [`PhylodynamicTrack`] draws the interval as a quiet ribbon, the
//! estimate as a line and any epidemiological reference as an independent
//! guide. It renders upstream inference; it does not fit a coalescent model.

use crate::scale::Scale;
use crate::style::{legible_ticks, LinePattern, QuantitativeAxis};
use crate::svg::{num, text_rounded, text_width_strong, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{unbroken, DrawContext, Track};

/// Vertical transformation for a phylodynamic estimate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PhylodynamicScale {
    /// Values are positioned without transformation.
    #[default]
    Linear,
    /// Positive values are positioned on a base-ten logarithmic axis.
    Log10,
}

/// One time point in a phylodynamic trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct PhylodynamicPoint {
    /// Time coordinate, interpreted through the figure's shared x scale.
    pub time: u64,
    /// Point estimate.
    pub estimate: f64,
    lower: Option<f64>,
    upper: Option<f64>,
    label: Option<String>,
}

impl PhylodynamicPoint {
    /// Creates a point estimate at integer time coordinate `time`.
    pub fn new(time: u64, estimate: f64) -> Self {
        PhylodynamicPoint {
            time,
            estimate,
            lower: None,
            upper: None,
            label: None,
        }
    }

    /// Adds a lower and upper uncertainty bound.
    pub fn interval(mut self, lower: f64, upper: f64) -> Self {
        if lower.is_finite() && upper.is_finite() && upper >= lower {
            self.lower = Some(lower);
            self.upper = Some(upper);
        }
        self
    }

    /// Adds a human-readable time, model or event label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Returns the supplied uncertainty interval, when valid.
    pub fn bounds(&self) -> Option<(f64, f64)> {
        Some((self.lower?, self.upper?))
    }
}

/// A skyline-like trajectory with an optional credible/confidence ribbon.
#[derive(Debug, Clone)]
pub struct PhylodynamicTrack {
    points: Vec<PhylodynamicPoint>,
    label: Option<String>,
    height: f64,
    scale: PhylodynamicScale,
    unit: String,
    color: Option<String>,
    reference: Option<(f64, String)>,
    show_points: bool,
    show_interval: bool,
}

impl PhylodynamicTrack {
    /// Creates a trajectory from `points`, sorted by time when drawn.
    pub fn new(points: impl Into<Vec<PhylodynamicPoint>>) -> Self {
        PhylodynamicTrack {
            points: points.into(),
            label: None,
            height: 126.0,
            scale: PhylodynamicScale::Linear,
            unit: String::new(),
            color: None,
            reference: None,
            show_points: true,
            show_interval: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the total track height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        if height.is_finite() {
            self.height = height.max(44.0);
        }
        self
    }

    /// Chooses a linear or base-ten logarithmic vertical scale.
    pub fn scale(mut self, scale: PhylodynamicScale) -> Self {
        self.scale = scale;
        self
    }

    /// Adds a unit such as `Ne`, `R`, or `lineages / year` to tooltips and axis.
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self
    }

    /// Overrides the theme accent used for the estimate.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Adds an independent reference guide such as `R = 1`.
    pub fn reference(mut self, value: f64, label: impl Into<String>) -> Self {
        if value.is_finite() {
            self.reference = Some((value, label.into()));
        }
        self
    }

    /// Shows or hides point markers without removing exact tooltips.
    pub fn show_points(mut self, show: bool) -> Self {
        self.show_points = show;
        self
    }

    /// Shows or hides uncertainty ribbons.
    pub fn show_interval(mut self, show: bool) -> Self {
        self.show_interval = show;
        self
    }

    /// Returns the source points in their supplied order.
    pub fn points(&self) -> &[PhylodynamicPoint] {
        &self.points
    }

    fn transformed(&self, value: f64) -> Option<f64> {
        if !value.is_finite() {
            return None;
        }
        match self.scale {
            PhylodynamicScale::Linear => Some(value),
            PhylodynamicScale::Log10 if value > 0.0 => Some(value.log10()),
            PhylodynamicScale::Log10 => None,
        }
    }
}

/// Tick labels as plain numbers, with the unit once on the highest.
///
/// Plain rather than compact: an effective population size is read off this
/// axis as a number, so it is written as one, 10000 and not 10k, and a value
/// too large to write out goes in exponent form, which is still a number.
/// Each label carries the decimals its own size needs, since a log axis runs
/// from 0.01 to 100 with no one step between them.
fn plain_labels(ticks: &[f64], unit: &str) -> Vec<String> {
    let top = ticks
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let step = ticks
        .windows(2)
        .map(|pair| (pair[1] - pair[0]).abs())
        .filter(|gap| *gap > 0.0)
        .fold(f64::INFINITY, f64::min);
    ticks
        .iter()
        .map(|&value| {
            let own = if value != 0.0 && value.abs() < 1.0 {
                -(value.abs().log10().floor())
            } else {
                0.0
            };
            let from_step = if step.is_finite() && step < 1.0 {
                -(step.log10().floor())
            } else {
                0.0
            };
            let decimals = own.max(from_step).clamp(0.0, 12.0) as u32;
            let number = text_rounded(value, decimals);
            if value == top && !unit.is_empty() {
                format!("{number} {unit}")
            } else {
                number
            }
        })
        .collect()
}

/// Powers of ten inside a range given in log10 units, with the 2 and 5
/// between them when the range spans too few decades to hold two powers.
fn log_ticks(low: f64, high: f64) -> Vec<f64> {
    if !(low.is_finite() && high.is_finite()) || high <= low {
        return Vec::new();
    }
    let first = low.floor() as i32;
    let last = high.ceil() as i32;
    let inside = |value: f64| {
        let exponent = value.log10();
        exponent >= low - 1e-9 && exponent <= high + 1e-9 && value.is_finite()
    };
    // Every power when there are a handful, every second or fifth or tenth
    // when a diverged estimate stretches the axis over hundreds of decades.
    let decades = (last - first).max(1);
    let stride = [1, 2, 5, 10, 20, 50, 100]
        .into_iter()
        .find(|stride| decades / stride <= 6)
        .unwrap_or(100);
    let powers: Vec<f64> = (first..=last)
        .filter(|k| k.rem_euclid(stride) == 0)
        .filter_map(|k| power_of_ten(1, k))
        .filter(|v| inside(*v))
        .collect();
    if powers.len() >= 2 {
        return powers;
    }
    let multiples: Vec<f64> = (first..=last)
        .flat_map(|k| [1, 2, 5].map(|m| power_of_ten(m, k)))
        .flatten()
        .filter(|v| inside(*v))
        .collect();
    if multiples.len() >= 2 {
        return multiples;
    }
    // Less than a factor of two and a half between the ends, with no round
    // value inside: the ends themselves, so the panel still has a scale.
    [low, high]
        .iter()
        .map(|exponent| 10f64.powf(*exponent))
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect()
}

/// `mantissa` times ten to the `exponent`, as the nearest double to the
/// decimal rather than a product that drifts: `10f64.powi(300)` is
/// 1.0000000000000006e300, which is not a label anyone should have to read.
fn power_of_ten(mantissa: u8, exponent: i32) -> Option<f64> {
    format!("{mantissa}e{exponent}").parse().ok()
}

impl Track for PhylodynamicTrack {
    fn noun(&self) -> &str {
        "a phylodynamic trajectory"
    }

    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, _theme: &Theme) -> f64 {
        58.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        // Each survivor is kept beside the value it was transformed to, rather
        // than transformed again at every use. The second read used to fall
        // back to the bottom of the domain when it came out empty, which is a
        // value standing in for the absence of one: an estimate nothing could
        // be computed from would have been drawn as the smallest estimate in
        // the panel, and the tooltip beside it would have said otherwise.
        let mut points: Vec<(&PhylodynamicPoint, f64)> = self
            .points
            .iter()
            .filter(|point| ctx.region.contains(point.time))
            .filter_map(|point| Some((point, self.transformed(point.estimate)?)))
            .collect();
        points.sort_by_key(|(point, _)| point.time);
        if points.is_empty() {
            return;
        }

        let header = ctx.px(19.0).min(ctx.band.h * 0.24);
        let plot_top = ctx.band.y + header + ctx.px(2.0);
        let plot_bottom = ctx.band.bottom() - ctx.px(7.0);
        let plot_height = (plot_bottom - plot_top).max(2.0);
        let mut values = Vec::new();
        for (point, estimate) in &points {
            values.push(*estimate);
            if self.show_interval {
                if let Some((lower, upper)) = point.bounds() {
                    if let Some(value) = self.transformed(lower) {
                        values.push(value);
                    }
                    if let Some(value) = self.transformed(upper) {
                        values.push(value);
                    }
                }
            }
        }
        if let Some((reference, _)) = &self.reference {
            if let Some(value) = self.transformed(*reference) {
                values.push(value);
            }
        }
        let lowest = values.iter().copied().fold(f64::MAX, f64::min);
        let highest = values.iter().copied().fold(f64::MIN, f64::max);
        let padding = if highest <= lowest {
            highest.abs().max(1.0) * 0.12
        } else {
            (highest - lowest) * 0.08
        };
        // Both ends are widened through the same guard, because both used to
        // leave the range the same way: an estimate of `f64::MAX` plus any
        // breathing room at all is infinity, and so is `-f64::MAX` minus it.
        let minimum = widened(lowest, -padding);
        let maximum = widened(highest, padding);
        // The distance between the two ends is not a number in every panel
        // where both ends are. An estimate of `1e308` above one of `-1e308` is
        // what a diverged optimiser writes, two ordinary finite numbers, and
        // their difference overflows on its own. Every coordinate in this
        // track is a fraction of that distance, so the overflow does not stay
        // where it happened: it arrives at the writer as `inf / inf`, which is
        // not a number, in every circle, every ribbon vertex and every axis
        // label. Halving both ends before subtracting brings them back inside
        // the range and leaves the fraction exactly the fraction it was, since
        // halving a double is exact at any magnitude that can overflow. The
        // factor is one wherever the subtraction already worked, so no panel
        // that reads correctly today moves by a thousandth of a pixel.
        let shrink = if (maximum - minimum).is_finite() {
            1.0
        } else {
            0.5
        };
        let span = maximum * shrink - minimum * shrink;
        let y_of =
            |value: f64| plot_bottom - (value * shrink - minimum * shrink) / span * plot_height;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());

        let ribbon = self.show_interval && points.iter().any(|(point, _)| point.bounds().is_some());
        self.draw_header(ctx, &color, ribbon);
        // Round values in the units the estimate is in, not four equal cuts
        // of the transformed range: on a log axis those came out as 77.988 and
        // 2167.255, which are nobody's idea of a population size.
        let ticks: Vec<f64> = match self.scale {
            PhylodynamicScale::Linear => QuantitativeAxis::new().values(minimum, maximum),
            PhylodynamicScale::Log10 => log_ticks(minimum, maximum),
        };
        let size = ctx.theme.font_size - 1.0;
        let placed: Vec<f64> = ticks
            .iter()
            .map(|&original| match self.scale {
                PhylodynamicScale::Linear => original,
                PhylodynamicScale::Log10 => original.log10(),
            })
            .collect();
        let shown = legible_ticks(&placed, y_of, size);
        let labels = plain_labels(&ticks, &self.unit);
        for ((&original, &transformed), label) in ticks.iter().zip(&placed).zip(&labels) {
            if !shown.contains(&transformed) {
                continue;
            }
            // A rule is a promise that this height means this value, so a rule
            // whose value cannot be written is a promise this track cannot
            // keep. The log axis is where it happens: raising ten to the top of
            // a domain that reaches `f64::MAX` overflows, and the rule used to
            // be drawn anyway, labelled `inf`. Nothing above the largest number
            // there is can be named, and an axis that names it anyway is a
            // figure that will be read and believed. The bottom of the domain
            // always survives this, so the panel is never left with no scale.
            if !original.is_finite() {
                continue;
            }
            let y = y_of(transformed);
            ctx.svg.line(
                ctx.band.x,
                y,
                ctx.band.right(),
                y,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
            ctx.svg.text(
                ctx.axis.right() - ctx.px(4.0),
                y + size * 0.35,
                label,
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }

        if let Some((reference, label)) = &self.reference {
            if let Some(value) = self.transformed(*reference) {
                let y = y_of(value);
                ctx.svg.begin_titled(&format!(
                    "reference {} {}",
                    text_rounded(*reference, 5),
                    label
                ));
                ctx.svg.line_pattern(
                    ctx.band.x,
                    y,
                    ctx.band.right(),
                    y,
                    ctx.theme.color(1),
                    ctx.theme.tokens.hairline.max(1.0),
                    LinePattern::Dashed,
                );
                ctx.svg.end_group();
            }
        }

        if self.show_interval {
            let mut segments = Vec::<Vec<(f64, f64, f64)>>::new();
            let mut segment = Vec::new();
            for (point, _) in &points {
                let transformed_bounds = point.bounds().and_then(|(lower, upper)| {
                    Some((self.transformed(lower)?, self.transformed(upper)?))
                });
                if let Some((lower, upper)) = transformed_bounds {
                    segment.push((ctx.scale.x_center(point.time), y_of(lower), y_of(upper)));
                } else if !segment.is_empty() {
                    segments.push(std::mem::take(&mut segment));
                }
            }
            if !segment.is_empty() {
                segments.push(segment);
            }
            for segment in segments.into_iter().filter(|segment| segment.len() >= 2) {
                let mut path = format!("M {} {}", num(segment[0].0), num(segment[0].2));
                for point in segment.iter().skip(1) {
                    path.push_str(&format!(" L {} {}", num(point.0), num(point.2)));
                }
                for point in segment.iter().rev() {
                    path.push_str(&format!(" L {} {}", num(point.0), num(point.1)));
                }
                path.push_str(" Z");
                ctx.svg.begin_titled("uncertainty interval");
                ctx.svg.path(&path, &color, 0.16);
                ctx.svg.end_group();
            }
        }

        // The line stops at a time whose estimate has no place on this axis,
        // because it is not a finite number or, on a log axis, not above
        // nought, and starts again after it. Drawn through, it ran straight
        // from the estimate before to the one after, as though the trajectory
        // had passed between them at a time where nothing put it, which the
        // ribbon already refused to draw. A time that holds another estimate,
        // one that can be placed, is not a gap.
        let mut steps: Vec<(u64, Option<(f64, f64)>)> = points
            .iter()
            .map(|(point, estimate)| {
                let at = (ctx.scale.x_center(point.time), y_of(*estimate));
                (point.time, Some(at))
            })
            .collect();
        steps.extend(
            self.points
                .iter()
                .filter(|point| ctx.region.contains(point.time))
                .filter(|point| self.transformed(point.estimate).is_none())
                .filter(|point| {
                    points
                        .binary_search_by_key(&point.time, |(placed, _)| placed.time)
                        .is_err()
                })
                .map(|point| (point.time, None)),
        );
        steps.sort_by_key(|(time, _)| *time);
        for line in unbroken(steps.into_iter().map(|(_, step)| step)) {
            ctx.svg
                .polyline(&line, &color, ctx.theme.tokens.stroke.max(1.8));
        }
        for (point, estimate) in points {
            let x = ctx.scale.x_center(point.time);
            let y = y_of(estimate);
            ctx.svg.begin_titled(&point_title(point, &self.unit));
            if self.show_points {
                ctx.svg.circle_ringed(
                    x,
                    y,
                    ctx.theme.tokens.marker_radius * 0.72,
                    &color,
                    ctx.theme.surface(),
                    ctx.theme.tokens.hairline.max(0.8),
                );
            } else {
                ctx.svg.circle(x, y, 1.2, &color);
            }
            ctx.svg.end_group();
        }

        // Last, so no point is drawn over the words, and at the left, over the
        // line, as a scan writes its threshold. At the right end they sat on
        // the latest estimate, which is the one a reader looks at first.
        if let Some((reference, label)) = &self.reference {
            if let Some(value) = self.transformed(*reference).filter(|_| !label.is_empty()) {
                let y = y_of(value);
                let size = ctx.theme.font_size * 0.74;
                let gap = ctx.theme.tokens.row_gap;
                let baseline = if y - gap - size * 0.8 >= plot_top {
                    y - gap
                } else {
                    y + gap + size * 0.8
                };
                ctx.svg.text(
                    ctx.band.x + ctx.theme.tokens.label_gap,
                    baseline,
                    label,
                    ctx.theme.color(1),
                    size,
                    Anchor::Start,
                );
            }
        }
    }
}

impl PhylodynamicTrack {
    /// A key to what is drawn: the line is the estimate, the ribbon where
    /// there is one its interval, and the scale where it is not the plain one.
    /// It used to read `trajectory · linear`, which named the default and
    /// left the ribbon to be guessed at.
    fn draw_header(&self, ctx: &mut DrawContext<'_>, color: &str, ribbon: bool) {
        let size = ctx.theme.font_size * 0.76;
        let top = ctx.band.y + ctx.px(2.0);
        let height = ctx.px(15.0);
        let middle = top + height / 2.0;
        let baseline = middle + size * 0.34;
        let swatch = ctx.px(17.0);
        let gap = ctx.px(6.0);
        let log = self.scale == PhylodynamicScale::Log10;
        let mut width = ctx.px(8.0) + swatch + gap + text_width_strong("estimate", size);
        if ribbon {
            width += ctx.px(12.0) + swatch + gap + text_width_strong("interval", size);
        }
        if log {
            width += ctx.px(12.0) + text_width_strong("log scale", size);
        }
        width += ctx.px(8.0);
        let chip = mix(ctx.theme.surface(), &ctx.theme.rule, 0.32);
        ctx.svg
            .rect_rounded(ctx.band.x + 2.0, top, width, height, height / 2.0, &chip);
        let mut x = ctx.band.x + 2.0 + ctx.px(8.0);
        ctx.svg.line(x, middle, x + swatch, middle, color, 2.0);
        x += swatch + gap;
        ctx.svg.text_bold(
            x,
            baseline,
            "estimate",
            &ctx.theme.muted,
            size,
            Anchor::Start,
        );
        x += text_width_strong("estimate", size);
        if ribbon {
            x += ctx.px(12.0);
            let tall = ctx.px(8.0);
            ctx.svg
                .rect_opacity(x, middle - tall / 2.0, swatch, tall, color, 0.16);
            x += swatch + gap;
            ctx.svg.text_bold(
                x,
                baseline,
                "interval",
                &ctx.theme.muted,
                size,
                Anchor::Start,
            );
            x += text_width_strong("interval", size);
        }
        if log {
            x += ctx.px(12.0);
            ctx.svg.text_bold(
                x,
                baseline,
                "log scale",
                &ctx.theme.muted,
                size,
                Anchor::Start,
            );
        }
    }
}

/// Moves one end of the value domain out by `padding`, staying in range.
///
/// The padding is breathing room around the data rather than part of it, so an
/// end that cannot move that far stays where the data left it. Moving it anyway
/// is how a panel of real estimates came to be drawn as an empty box: an end of
/// `f64::MAX` widened by a fraction of itself is infinity, every height in the
/// track is measured against that end, and a height that is not a number is
/// dropped by the writer rather than guessed at. The trajectory disappeared and
/// its tooltips stayed, promising values with nothing under them.
fn widened(bound: f64, padding: f64) -> f64 {
    let widened = bound + padding;
    if widened.is_finite() {
        widened
    } else {
        bound
    }
}

fn point_title(point: &PhylodynamicPoint, unit: &str) -> String {
    let suffix = if unit.is_empty() {
        String::new()
    } else {
        format!(" {unit}")
    };
    let mut parts = vec![
        // As the ruler under it writes it, as a surveillance panel does.
        format!("time {}", point.time.saturating_add(1)),
        format!("estimate {}{}", text_rounded(point.estimate, 6), suffix),
    ];
    if let Some((lower, upper)) = point.bounds() {
        parts.push(format!(
            "interval {}{} to {}{}",
            text_rounded(lower, 6),
            suffix,
            text_rounded(upper, 6),
            suffix
        ));
    }
    if let Some(label) = &point.label {
        parts.push(label.clone());
    }
    parts.join(" | ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Figure, Region};

    /// A tooltip calls a time what the ruler under it calls it: time 0 is
    /// the first, which a ruler numbers 1. It said 0, so week 12 was week 11
    /// on hover.
    #[test]
    fn a_time_is_written_as_the_ruler_writes_it() {
        let svg = Figure::new(Region::new("week", 0, 3).unwrap())
            .push(PhylodynamicTrack::new(vec![
                PhylodynamicPoint::new(0, 1.2),
                PhylodynamicPoint::new(2, 0.9),
            ]))
            .to_svg();
        assert!(svg.contains("<title>time 1 | estimate 1.2"), "{svg}");
        assert!(svg.contains("<title>time 3 | estimate 0.9"), "{svg}");
    }

    #[test]
    fn intervals_references_and_exact_values_are_independent() {
        let track = PhylodynamicTrack::new(vec![
            PhylodynamicPoint::new(1, 0.8).interval(0.5, 1.1),
            PhylodynamicPoint::new(2, 1.4)
                .interval(1.0, 1.9)
                .label("growth phase"),
        ])
        .unit("R")
        .reference(1.0, "epidemic threshold");
        let svg = Figure::new(Region::new("week", 0, 4).unwrap())
            .push(track)
            .to_svg();
        assert!(svg.contains("estimate 1.4 R"), "{svg}");
        assert!(svg.contains("interval 1 R to 1.9 R"), "{svg}");
        assert!(svg.contains("epidemic threshold"), "{svg}");
        assert!(svg.contains("fill-opacity=\"0.16\""), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }

    #[test]
    fn logarithmic_mode_omits_nonpositive_points_without_nonfinite_output() {
        let svg = Figure::new(Region::new("time", 0, 5).unwrap())
            .push(
                PhylodynamicTrack::new(vec![
                    PhylodynamicPoint::new(1, 0.0),
                    PhylodynamicPoint::new(2, 10.0).interval(2.0, 30.0),
                ])
                .scale(PhylodynamicScale::Log10),
            )
            .to_svg();
        assert!(!svg.contains("estimate 0"), "{svg}");
        assert!(svg.contains("estimate 10"), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
        assert!(!svg.contains("Infinity"), "{svg}");
    }

    #[test]
    fn a_missing_interval_breaks_the_ribbon_instead_of_being_interpolated() {
        let svg = Figure::new(Region::new("time", 0, 5).unwrap())
            .push(PhylodynamicTrack::new(vec![
                PhylodynamicPoint::new(1, 2.0).interval(1.0, 3.0),
                PhylodynamicPoint::new(2, 4.0),
                PhylodynamicPoint::new(3, 3.0).interval(2.0, 5.0),
            ]))
            .to_svg();
        assert!(!svg.contains("uncertainty interval"), "{svg}");
        assert!(svg.contains("estimate 4"), "{svg}");
    }

    #[test]
    fn a_missing_estimate_breaks_the_line_instead_of_being_interpolated() {
        // The third time has no estimate, and the line used to run straight
        // from the second to the fourth, drawing one at the third that nothing
        // estimated. The ribbon beside it already stopped there.
        let svg = Figure::new(Region::new("time", 0, 7).unwrap())
            .show_region_label(false)
            .push(PhylodynamicTrack::new(vec![
                PhylodynamicPoint::new(1, 2.0),
                PhylodynamicPoint::new(2, 3.0),
                PhylodynamicPoint::new(3, f64::NAN),
                PhylodynamicPoint::new(4, 5.0),
                PhylodynamicPoint::new(5, 4.0),
            ]))
            .to_svg();
        let lines = crate::track::polylines(&svg);
        assert_eq!(lines.len(), 2, "one line either side of the gap: {lines:?}");
        assert!(lines.iter().all(|line| line.len() == 2), "{lines:?}");
    }

    #[test]
    fn an_estimate_a_log_axis_cannot_place_breaks_the_line_too() {
        // Nought has no logarithm, so it has no height on this axis. Joined
        // across, the trajectory read as though it had stayed between its
        // neighbours when the estimate says it fell to nothing.
        let svg = Figure::new(Region::new("time", 0, 5).unwrap())
            .show_region_label(false)
            .push(
                PhylodynamicTrack::new(vec![
                    PhylodynamicPoint::new(1, 10.0),
                    PhylodynamicPoint::new(2, 0.0),
                    PhylodynamicPoint::new(3, 100.0),
                ])
                .scale(PhylodynamicScale::Log10),
            )
            .to_svg();
        let lines = crate::track::polylines(&svg);
        assert!(lines.is_empty(), "joined across the gap: {lines:?}");
        assert_eq!(
            svg.matches("<circle").count(),
            4,
            "both placed points drawn"
        );
    }

    #[test]
    fn a_time_holding_an_estimate_is_not_a_gap_because_another_row_there_is_missing() {
        // Two rows at one time, one of them empty. The time has an estimate,
        // so the line goes through it whichever order the rows came in.
        for order in [[f64::NAN, 3.0], [3.0, f64::NAN]] {
            let svg = Figure::new(Region::new("time", 0, 5).unwrap())
                .show_region_label(false)
                .push(PhylodynamicTrack::new(vec![
                    PhylodynamicPoint::new(1, 2.0),
                    PhylodynamicPoint::new(2, order[0]),
                    PhylodynamicPoint::new(2, order[1]),
                    PhylodynamicPoint::new(3, 4.0),
                ]))
                .to_svg();
            let lines = crate::track::polylines(&svg);
            assert_eq!(lines.len(), 1, "{order:?}: {lines:?}");
            assert_eq!(lines[0].len(), 3, "{order:?}: {lines:?}");
        }
    }

    /// The height of every mark in the document, in the order it was drawn.
    fn marker_heights(svg: &str) -> Vec<f64> {
        svg.split("cy=\"")
            .skip(1)
            .filter_map(|piece| piece.split('"').next()?.parse().ok())
            .collect()
    }

    /// Every vertex of the uncertainty ribbon, in the order it is visited.
    fn ribbon_vertices(svg: &str) -> Vec<(f64, f64)> {
        let d = svg
            .split("uncertainty interval</title>")
            .nth(1)
            .and_then(|rest| rest.split("d=\"").nth(1))
            .and_then(|rest| rest.split('"').next())
            .expect("the ribbon is drawn");
        let numbers: Vec<f64> = d
            .split(' ')
            .filter_map(|piece| piece.parse::<f64>().ok())
            .collect();
        numbers
            .chunks_exact(2)
            .map(|two| (two[0], two[1]))
            .collect()
    }

    #[test]
    fn padding_that_would_leave_the_range_leaves_the_domain_where_the_data_left_it() {
        // The breathing room around the data is a twelfth of it, and a twelfth
        // more than `f64::MAX` is infinity. The domain then had an end that was
        // not a number, every height in the panel is measured against that end,
        // and both estimates came out at the same fraction of an infinite
        // distance: two values as far apart as doubles go, drawn on one line at
        // the floor of the track, under an axis with no rules left on it.
        let figure = Figure::new(Region::new("time", 0, 4).unwrap())
            .show_region_label(false)
            .push(PhylodynamicTrack::new(vec![
                PhylodynamicPoint::new(1, 0.0),
                PhylodynamicPoint::new(3, f64::MAX),
            ]));
        let svg = figure.to_svg();
        let heights = marker_heights(&svg);
        assert_eq!(heights.len(), 4, "two ringed markers, two circles each");
        assert!(
            heights[0] > heights[3] + 1.0,
            "nought is drawn at {} and f64::MAX at {}",
            heights[0],
            heights[3]
        );
        // A span no round step fits is labelled at its two ends, and the top
        // one names the largest number there is rather than an overflow.
        assert_eq!(svg.matches(r#"text-anchor="end""#).count(), 2, "{svg}");
        assert!(svg.contains(">1.7976931348623157e308<"), "{svg}");
    }

    #[test]
    fn a_domain_too_wide_to_subtract_still_draws_every_point_and_its_ribbon() {
        // `1e308` above `-1e308` is two ordinary finite estimates whose distance
        // is not a finite number, which is what a diverged optimiser writes. The
        // circles and the axis went quiet, because the writer drops a coordinate
        // that is not a number rather than guessing one, and the ribbon did not:
        // a path has no such guard, so the vertex was written as the nought that
        // stands for an unwritable number and the interval was painted from the
        // trajectory to the top of the sheet, over the track above it.
        let figure = Figure::new(Region::new("time", 0, 4).unwrap())
            .show_region_label(false)
            .push(PhylodynamicTrack::new(vec![
                PhylodynamicPoint::new(1, -1e308).interval(-1e308, 0.0),
                PhylodynamicPoint::new(2, 5.0).interval(0.0, 1e308),
                PhylodynamicPoint::new(3, 1e308),
            ]));
        let (_, height) = figure.dimensions();
        let svg = figure.to_svg();
        assert_eq!(svg.matches("<circle").count(), 6, "{svg}");
        assert_eq!(svg.matches(r#"text-anchor="end""#).count(), 2, "{svg}");
        for (x, y) in ribbon_vertices(&svg) {
            assert!(
                x > 0.0 && y > 0.0 && y <= height,
                "a ribbon vertex at {x},{y} on a figure {height} tall"
            );
        }
        assert!(!svg.contains("NaN"), "{svg}");
    }

    #[test]
    fn a_log_axis_draws_no_rule_it_cannot_name() {
        // Ten raised to the top of a log domain that reaches `f64::MAX` is
        // infinity, and the rule at that height used to be drawn and labelled
        // with it. A missing rule costs the reader a gridline; a rule reading
        // `inf` costs them the figure, because it is a number they will believe.
        // The rules below it are powers of ten, every one of them nameable.
        let svg = Figure::new(Region::new("time", 0, 4).unwrap())
            .show_region_label(false)
            .push(
                PhylodynamicTrack::new(vec![
                    PhylodynamicPoint::new(1, 1.0),
                    PhylodynamicPoint::new(3, f64::MAX),
                ])
                .scale(PhylodynamicScale::Log10),
            )
            .to_svg();
        assert!(!svg.contains("inf"), "{svg}");
        let rules = svg.matches(r#"text-anchor="end""#).count();
        assert!((2..=6).contains(&rules), "{rules} rules");
        assert!(
            svg.contains(">1e300<"),
            "a power of ten near the top: {svg}"
        );
        assert_eq!(svg.matches("<circle").count(), 4, "{svg}");
    }
}
