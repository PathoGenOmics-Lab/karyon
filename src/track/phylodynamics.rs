//! Time-varying phylodynamic estimates with uncertainty kept visible.
//!
//! A skyline, effective population size, reproductive number or lineage
//! growth trajectory is an estimate with a time coordinate and often an
//! interval. [`PhylodynamicTrack`] draws the interval as a quiet ribbon, the
//! estimate as a line and any epidemiological reference as an independent
//! guide. It renders upstream inference; it does not fit a coalescent model.

use crate::scale::Scale;
use crate::style::LinePattern;
use crate::svg::{num, text_rounded, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

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

impl Track for PhylodynamicTrack {
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
        let mut points: Vec<&PhylodynamicPoint> = self
            .points
            .iter()
            .filter(|point| ctx.region.contains(point.time))
            .filter(|point| self.transformed(point.estimate).is_some())
            .collect();
        points.sort_by_key(|point| point.time);
        if points.is_empty() {
            return;
        }

        let header = ctx.px(19.0).min(ctx.band.h * 0.24);
        let plot_top = ctx.band.y + header + ctx.px(2.0);
        let plot_bottom = ctx.band.bottom() - ctx.px(7.0);
        let plot_height = (plot_bottom - plot_top).max(2.0);
        let mut values = Vec::new();
        for point in &points {
            if let Some(value) = self.transformed(point.estimate) {
                values.push(value);
            }
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
        let mut minimum = values.iter().copied().fold(f64::MAX, f64::min);
        let mut maximum = values.iter().copied().fold(f64::MIN, f64::max);
        if maximum <= minimum {
            let padding = maximum.abs().max(1.0) * 0.12;
            minimum -= padding;
            maximum += padding;
        } else {
            let padding = (maximum - minimum) * 0.08;
            minimum -= padding;
            maximum += padding;
        }
        let y_of = |value: f64| plot_bottom - (value - minimum) / (maximum - minimum) * plot_height;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());

        self.draw_header(ctx, &color);
        for index in 0..=3 {
            let fraction = index as f64 / 3.0;
            let transformed = minimum + fraction * (maximum - minimum);
            let y = y_of(transformed);
            ctx.svg.line(
                ctx.band.x,
                y,
                ctx.band.right(),
                y,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
            let original = match self.scale {
                PhylodynamicScale::Linear => transformed,
                PhylodynamicScale::Log10 => 10f64.powf(transformed),
            };
            let label = if self.unit.is_empty() {
                text_rounded(original, 3)
            } else {
                format!("{} {}", text_rounded(original, 3), self.unit)
            };
            ctx.svg.text(
                ctx.axis.right() - ctx.px(4.0),
                y + ctx.theme.font_size * 0.28,
                &label,
                &ctx.theme.muted,
                ctx.theme.font_size * 0.72,
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
                ctx.svg.text(
                    ctx.band.right() - ctx.px(4.0),
                    y - ctx.px(3.0),
                    label,
                    ctx.theme.color(1),
                    ctx.theme.font_size * 0.74,
                    Anchor::End,
                );
                ctx.svg.end_group();
            }
        }

        if self.show_interval {
            let mut segments = Vec::<Vec<(f64, f64, f64)>>::new();
            let mut segment = Vec::new();
            for point in &points {
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

        let line: Vec<(f64, f64)> = points
            .iter()
            .filter_map(|point| {
                Some((
                    ctx.scale.x_center(point.time),
                    y_of(self.transformed(point.estimate)?),
                ))
            })
            .collect();
        ctx.svg
            .polyline(&line, &color, ctx.theme.tokens.stroke.max(1.8));
        for point in points {
            let x = ctx.scale.x_center(point.time);
            let y = y_of(self.transformed(point.estimate).unwrap_or(minimum));
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
    }
}

impl PhylodynamicTrack {
    fn draw_header(&self, ctx: &mut DrawContext<'_>, color: &str) {
        let size = ctx.theme.font_size * 0.76;
        let top = ctx.band.y + ctx.px(2.0);
        let height = ctx.px(15.0);
        let chip = mix(ctx.theme.surface(), &ctx.theme.rule, 0.32);
        let scale = match self.scale {
            PhylodynamicScale::Linear => "linear",
            PhylodynamicScale::Log10 => "log10",
        };
        ctx.svg
            .rect_rounded(ctx.band.x + 2.0, top, 112.0, height, height / 2.0, &chip);
        ctx.svg.line(
            ctx.band.x + 10.0,
            top + height / 2.0,
            ctx.band.x + 27.0,
            top + height / 2.0,
            color,
            2.0,
        );
        ctx.svg.text_bold(
            ctx.band.x + 33.0,
            top + height / 2.0 + size * 0.34,
            &format!("trajectory · {scale}"),
            &ctx.theme.muted,
            size,
            Anchor::Start,
        );
    }
}

fn point_title(point: &PhylodynamicPoint, unit: &str) -> String {
    let suffix = if unit.is_empty() {
        String::new()
    } else {
        format!(" {unit}")
    };
    let mut parts = vec![
        format!("time {}", point.time),
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
}
