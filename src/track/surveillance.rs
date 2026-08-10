//! Temporal lineage surveillance as frequencies, counts and growth alerts.
//!
//! Surveillance asks which lineage is replacing which, and whether a signal is
//! large enough to inspect. [`SurveillanceTrack`] can draw stacked composition
//! or comparable trajectories while retaining observed counts, denominators
//! and alert reasons in exact SVG tooltips. It does not smooth or extrapolate
//! frequencies; upstream estimates can be supplied as counts and totals.

use std::collections::BTreeSet;

use crate::scale::Scale;
use crate::style::Symbol;
use crate::svg::{num, text_rounded, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// Quantity drawn for each lineage observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SurveillanceMetric {
    /// Count divided by the supplied total at that time.
    #[default]
    Frequency,
    /// Raw observed count.
    Count,
}

/// Arrangement of several lineage trajectories.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SurveillanceStyle {
    /// Lineage frequencies or counts accumulated into filled bands.
    #[default]
    Stacked,
    /// One independent line per lineage on a common axis.
    Lines,
}

/// One lineage observation at one integer time coordinate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurveillanceObservation {
    /// Time coordinate interpreted through the figure's x scale.
    pub time: u64,
    /// Stable lineage, clade, genotype or mutation label.
    pub lineage: String,
    /// Observed number assigned to this lineage.
    pub count: u64,
    /// Total observations used as the frequency denominator.
    pub total: u64,
    label: Option<String>,
}

impl SurveillanceObservation {
    /// Creates one observed lineage count and its denominator.
    pub fn new(time: u64, lineage: impl Into<String>, count: u64, total: u64) -> Self {
        SurveillanceObservation {
            time,
            lineage: lineage.into(),
            count,
            total,
            label: None,
        }
    }

    /// Adds a location, assay, sampling frame or alert note.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Returns the observed frequency when the denominator is positive.
    pub fn frequency(&self) -> Option<f64> {
        (self.total > 0).then_some(self.count as f64 / self.total as f64)
    }
}

/// A multi-lineage temporal surveillance panel.
#[derive(Debug, Clone)]
pub struct SurveillanceTrack {
    observations: Vec<SurveillanceObservation>,
    label: Option<String>,
    height: f64,
    metric: SurveillanceMetric,
    style: SurveillanceStyle,
    minimum_total: u64,
    frequency_alert: Option<f64>,
    growth_alert: Option<f64>,
    show_points: bool,
}

impl SurveillanceTrack {
    /// Creates a surveillance panel from observed lineage counts.
    pub fn new(observations: impl Into<Vec<SurveillanceObservation>>) -> Self {
        SurveillanceTrack {
            observations: observations.into(),
            label: None,
            height: 138.0,
            metric: SurveillanceMetric::Frequency,
            style: SurveillanceStyle::Stacked,
            minimum_total: 1,
            frequency_alert: None,
            growth_alert: None,
            show_points: true,
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
            self.height = height.max(48.0);
        }
        self
    }

    /// Chooses observed frequency or raw count.
    pub fn metric(mut self, metric: SurveillanceMetric) -> Self {
        self.metric = metric;
        self
    }

    /// Chooses stacked composition or independent trajectories.
    pub fn style(mut self, style: SurveillanceStyle) -> Self {
        self.style = style;
        self
    }

    /// Omits observations whose denominator is below this sampling floor.
    pub fn minimum_total(mut self, minimum: u64) -> Self {
        self.minimum_total = minimum;
        self
    }

    /// Flags observations at or above this frequency.
    pub fn frequency_alert(mut self, frequency: f64) -> Self {
        if frequency.is_finite() {
            self.frequency_alert = Some(frequency.clamp(0.0, 1.0));
        }
        self
    }

    /// Flags increases at or above this frequency change per observed step.
    pub fn growth_alert(mut self, change: f64) -> Self {
        if change.is_finite() && change >= 0.0 {
            self.growth_alert = Some(change);
        }
        self
    }

    /// Shows or hides source-observation markers.
    pub fn show_points(mut self, show: bool) -> Self {
        self.show_points = show;
        self
    }

    /// Returns source observations in supplied order.
    pub fn observations(&self) -> &[SurveillanceObservation] {
        &self.observations
    }

    fn valid_observations(&self, ctx: &DrawContext<'_>) -> Vec<&SurveillanceObservation> {
        self.observations
            .iter()
            .filter(|observation| ctx.region.contains(observation.time))
            .filter(|observation| observation.total >= self.minimum_total)
            .filter(|observation| match self.metric {
                SurveillanceMetric::Frequency => observation.total > 0,
                SurveillanceMetric::Count => true,
            })
            .collect()
    }

    fn value(&self, observation: &SurveillanceObservation) -> Option<f64> {
        match self.metric {
            SurveillanceMetric::Frequency => observation.frequency(),
            SurveillanceMetric::Count => Some(observation.count as f64),
        }
    }
}

impl Track for SurveillanceTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, _theme: &Theme) -> f64 {
        45.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let observations = self.valid_observations(ctx);
        if observations.is_empty() {
            return;
        }
        let mut lineages = Vec::<String>::new();
        for observation in &observations {
            if !lineages.contains(&observation.lineage) {
                lineages.push(observation.lineage.clone());
            }
        }
        let times: Vec<u64> = observations
            .iter()
            .map(|observation| observation.time)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if times.is_empty() || lineages.is_empty() {
            return;
        }

        let matrix: Vec<Vec<Option<f64>>> = lineages
            .iter()
            .map(|lineage| {
                times
                    .iter()
                    .map(|time| {
                        let matching: Vec<&SurveillanceObservation> = observations
                            .iter()
                            .filter(|observation| {
                                observation.time == *time && observation.lineage == *lineage
                            })
                            .copied()
                            .collect();
                        (matching.len() == 1)
                            .then(|| self.value(matching[0]))
                            .flatten()
                    })
                    .collect()
            })
            .collect();
        let complete_at: Vec<bool> = (0..times.len())
            .map(|time| matrix.iter().all(|series| series[time].is_some()))
            .collect();
        let maximum = match self.metric {
            SurveillanceMetric::Frequency => match self.style {
                SurveillanceStyle::Stacked => (0..times.len())
                    .filter(|time| complete_at[*time])
                    .map(|time| matrix.iter().filter_map(|series| series[time]).sum::<f64>())
                    .fold(1.0f64, f64::max),
                SurveillanceStyle::Lines => 1.0,
            },
            SurveillanceMetric::Count => match self.style {
                SurveillanceStyle::Stacked => (0..times.len())
                    .filter(|time| complete_at[*time])
                    .map(|time| matrix.iter().filter_map(|series| series[time]).sum::<f64>())
                    .fold(1.0f64, f64::max),
                SurveillanceStyle::Lines => matrix
                    .iter()
                    .flatten()
                    .flatten()
                    .copied()
                    .fold(1.0f64, f64::max),
            },
        };
        let header = ctx.px(21.0).min(ctx.band.h * 0.25);
        let top = ctx.band.y + header;
        let bottom = ctx.band.bottom() - ctx.px(6.0);
        let height = (bottom - top).max(2.0);
        let y_of = |value: f64| bottom - value / maximum * height;

        self.draw_legend(ctx, &lineages);
        for index in 0..=2 {
            let value = maximum * index as f64 / 2.0;
            let y = y_of(value);
            ctx.svg.line(
                ctx.band.x,
                y,
                ctx.band.right(),
                y,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
            let label = match self.metric {
                SurveillanceMetric::Frequency => {
                    format!("{}%", text_rounded(value * 100.0, 1))
                }
                SurveillanceMetric::Count => text_rounded(value, 1),
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

        match self.style {
            SurveillanceStyle::Stacked => {
                let mut lower = vec![0.0; times.len()];
                let segments = contiguous_true_ranges(&complete_at);
                for (lineage_index, series) in matrix.iter().enumerate() {
                    for range in &segments {
                        if range.len() < 2 {
                            continue;
                        }
                        let first = range.start;
                        let first_upper = lower[first] + series[first].unwrap_or_default();
                        let mut path = format!(
                            "M {} {}",
                            num(ctx.scale.x_center(times[first])),
                            num(y_of(first_upper))
                        );
                        for time_index in (first + 1)..range.end {
                            let upper = lower[time_index] + series[time_index].unwrap_or_default();
                            path.push_str(&format!(
                                " L {} {}",
                                num(ctx.scale.x_center(times[time_index])),
                                num(y_of(upper))
                            ));
                        }
                        for time_index in range.clone().rev() {
                            path.push_str(&format!(
                                " L {} {}",
                                num(ctx.scale.x_center(times[time_index])),
                                num(y_of(lower[time_index]))
                            ));
                        }
                        path.push_str(" Z");
                        ctx.svg.begin_titled(&format!(
                            "{} trajectory; {} complete time points",
                            lineages[lineage_index],
                            range.len()
                        ));
                        ctx.svg.path(&path, ctx.theme.color(lineage_index), 0.78);
                        ctx.svg.end_group();
                    }
                    for time_index in 0..times.len() {
                        if complete_at[time_index] {
                            lower[time_index] += series[time_index].unwrap_or_default();
                        }
                    }
                }
            }
            SurveillanceStyle::Lines => {
                for (lineage_index, series) in matrix.iter().enumerate() {
                    let mut segment = Vec::new();
                    for (time, value) in times.iter().zip(series) {
                        if let Some(value) = value {
                            segment.push((ctx.scale.x_center(*time), y_of(*value)));
                        } else {
                            draw_surveillance_line(ctx, lineage_index, &segment);
                            segment.clear();
                        }
                    }
                    draw_surveillance_line(ctx, lineage_index, &segment);
                }
            }
        }

        for (lineage_index, lineage) in lineages.iter().enumerate() {
            let mut history: Vec<&SurveillanceObservation> = observations
                .iter()
                .filter(|observation| observation.lineage == *lineage)
                .copied()
                .collect();
            history.sort_by_key(|observation| observation.time);
            for (history_index, observation) in history.iter().enumerate() {
                let value = self.value(observation).unwrap_or(0.0);
                let time_index = times.binary_search(&observation.time).unwrap_or_default();
                let unique = matrix[lineage_index][time_index].is_some();
                let complete = complete_at[time_index];
                if !unique || (self.style == SurveillanceStyle::Stacked && !complete) {
                    let reason = if !unique {
                        "duplicate lineage/time observations are not joined"
                    } else {
                        "incomplete lineage composition is not stacked"
                    };
                    ctx.svg.begin_titled(&format!(
                        "{reason} | {}",
                        observation_title(observation, None, None, None)
                    ));
                    ctx.svg.symbol_ringed(
                        ctx.scale.x_center(observation.time),
                        top + ctx.px(5.0 + lineage_index as f64 * 1.7),
                        ctx.px(2.7),
                        Symbol::Diamond,
                        &ctx.theme.muted,
                        ctx.theme.surface(),
                        ctx.px(0.8),
                    );
                    ctx.svg.end_group();
                    continue;
                }
                let plotted = match self.style {
                    SurveillanceStyle::Lines => value,
                    SurveillanceStyle::Stacked => matrix
                        .iter()
                        .take(lineage_index + 1)
                        .filter_map(|series| series[time_index])
                        .sum(),
                };
                let previous_frequency = history_index
                    .checked_sub(1)
                    .and_then(|index| history[index].frequency());
                let frequency = observation.frequency();
                let growth = frequency
                    .zip(previous_frequency)
                    .map(|(now, before)| now - before);
                let frequency_alert = self
                    .frequency_alert
                    .zip(frequency)
                    .is_some_and(|(threshold, value)| value >= threshold);
                let growth_alert = self
                    .growth_alert
                    .zip(growth)
                    .is_some_and(|(threshold, value)| value >= threshold);
                let title = observation_title(
                    observation,
                    growth,
                    frequency_alert.then_some(self.frequency_alert.unwrap_or_default()),
                    growth_alert.then_some(self.growth_alert.unwrap_or_default()),
                );
                let x = ctx.scale.x_center(observation.time);
                let y = y_of(plotted);
                ctx.svg.begin_titled(&title);
                if frequency_alert || growth_alert {
                    ctx.svg.symbol_ringed(
                        x,
                        y - ctx.px(5.0),
                        ctx.px(4.0),
                        Symbol::Triangle,
                        ctx.theme.color(1),
                        ctx.theme.surface(),
                        ctx.px(1.0),
                    );
                } else if self.show_points {
                    ctx.svg.circle_ringed(
                        x,
                        y,
                        ctx.theme.tokens.marker_radius * 0.55,
                        ctx.theme.color(lineage_index),
                        ctx.theme.surface(),
                        ctx.theme.tokens.hairline,
                    );
                } else {
                    ctx.svg.circle(x, y, 1.0, ctx.theme.color(lineage_index));
                }
                ctx.svg.end_group();
            }
        }
    }
}

fn contiguous_true_ranges(values: &[bool]) -> Vec<std::ops::Range<usize>> {
    let mut ranges = Vec::new();
    let mut start = None;
    for (index, value) in values.iter().copied().chain([false]).enumerate() {
        match (start, value) {
            (None, true) => start = Some(index),
            (Some(from), false) => {
                ranges.push(from..index);
                start = None;
            }
            _ => {}
        }
    }
    ranges
}

fn draw_surveillance_line(ctx: &mut DrawContext<'_>, lineage_index: usize, points: &[(f64, f64)]) {
    if points.len() >= 2 {
        ctx.svg.polyline(
            points,
            ctx.theme.color(lineage_index),
            ctx.theme.tokens.stroke.max(1.7),
        );
    }
}

impl SurveillanceTrack {
    fn draw_legend(&self, ctx: &mut DrawContext<'_>, lineages: &[String]) {
        let size = ctx.theme.font_size * 0.74;
        let top = ctx.band.y + ctx.px(2.0);
        let height = ctx.px(16.0);
        let chip = mix(ctx.theme.surface(), &ctx.theme.rule, 0.32);
        let mut x = ctx.band.x + ctx.px(2.0);
        for (index, lineage) in lineages.iter().enumerate() {
            let available = (ctx.band.right() - x).min(105.0);
            if available < 20.0 {
                break;
            }
            let visible = fit_lineage(lineage, (available - 21.0).max(0.0), size);
            let width = (text_width(&visible, size) + 25.0).min(available);
            ctx.svg
                .rect_rounded(x, top, width, height, height / 2.0, &chip);
            ctx.svg
                .circle(x + 9.0, top + height / 2.0, 3.0, ctx.theme.color(index));
            ctx.svg.text_bold(
                x + 16.0,
                top + height / 2.0 + size * 0.34,
                &visible,
                &ctx.theme.muted,
                size,
                Anchor::Start,
            );
            x += width + ctx.px(4.0);
        }
    }
}

fn fit_lineage(text: &str, available: f64, size: f64) -> String {
    if text_width(text, size) <= available {
        return text.to_string();
    }
    if available <= text_width("…", size) {
        return String::new();
    }
    let mut visible = String::new();
    for character in text.chars() {
        let candidate = format!("{visible}{character}…");
        if text_width(&candidate, size) > available {
            break;
        }
        visible.push(character);
    }
    visible.push('…');
    visible
}

fn observation_title(
    observation: &SurveillanceObservation,
    growth: Option<f64>,
    frequency_threshold: Option<f64>,
    growth_threshold: Option<f64>,
) -> String {
    let mut parts = vec![
        observation.lineage.clone(),
        format!("time {}", observation.time),
        format!("count {} of {}", observation.count, observation.total),
    ];
    if let Some(frequency) = observation.frequency() {
        parts.push(format!("frequency {}", text_rounded(frequency, 5)));
    }
    if let Some(growth) = growth {
        let sign = if growth > 0.0 { "+" } else { "" };
        parts.push(format!("change {sign}{}", text_rounded(growth, 5)));
    }
    if let Some(threshold) = frequency_threshold {
        parts.push(format!("frequency alert >= {}", text_rounded(threshold, 5)));
    }
    if let Some(threshold) = growth_threshold {
        parts.push(format!("growth alert >= {}", text_rounded(threshold, 5)));
    }
    if let Some(label) = &observation.label {
        parts.push(label.clone());
    }
    parts.join(" | ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Figure, Region};

    #[test]
    fn stacked_surveillance_keeps_counts_denominators_and_alert_reasons() {
        let observations = vec![
            SurveillanceObservation::new(1, "A", 90, 100),
            SurveillanceObservation::new(1, "B", 10, 100),
            SurveillanceObservation::new(2, "A", 60, 100),
            SurveillanceObservation::new(2, "B", 40, 100).label("regional cluster"),
        ];
        let svg = Figure::new(Region::new("week", 0, 4).unwrap())
            .push(
                SurveillanceTrack::new(observations)
                    .frequency_alert(0.35)
                    .growth_alert(0.20),
            )
            .to_svg();
        assert!(svg.contains("B | time 2 | count 40 of 100"), "{svg}");
        assert!(svg.contains("frequency alert &gt;= 0.35"), "{svg}");
        assert!(svg.contains("growth alert &gt;= 0.2"), "{svg}");
        assert!(svg.contains("regional cluster"), "{svg}");
        assert!(svg.contains("fill-opacity=\"0.78\""), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }

    #[test]
    fn sampling_floor_and_zero_denominators_are_explicit_omissions() {
        let svg = Figure::new(Region::new("week", 0, 5).unwrap())
            .push(
                SurveillanceTrack::new(vec![
                    SurveillanceObservation::new(1, "zero", 4, 0),
                    SurveillanceObservation::new(2, "thin", 2, 4),
                    SurveillanceObservation::new(3, "kept", 20, 40),
                ])
                .minimum_total(10),
            )
            .to_svg();
        assert!(!svg.contains("zero | time"), "{svg}");
        assert!(!svg.contains("thin | time"), "{svg}");
        assert!(svg.contains("kept | time 3"), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }

    #[test]
    fn missing_or_duplicate_observations_are_not_invented_as_zero() {
        let svg = Figure::new(Region::new("week", 0, 5).unwrap())
            .push(SurveillanceTrack::new(vec![
                SurveillanceObservation::new(1, "A", 40, 100),
                SurveillanceObservation::new(1, "B", 60, 100),
                SurveillanceObservation::new(2, "A", 50, 100),
                SurveillanceObservation::new(3, "A", 55, 100),
                SurveillanceObservation::new(3, "A", 56, 100),
                SurveillanceObservation::new(3, "B", 44, 100),
            ]))
            .to_svg();
        assert!(
            svg.contains("incomplete lineage composition is not stacked"),
            "{svg}"
        );
        assert!(
            svg.contains("duplicate lineage/time observations are not joined"),
            "{svg}"
        );
        assert!(!svg.contains("fill-opacity=\"0.78\""), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }
}
