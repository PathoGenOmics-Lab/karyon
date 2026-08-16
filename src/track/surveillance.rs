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

    /// Whether the two numbers can be read as a frequency at all.
    ///
    /// A part cannot exceed its whole, so a count above its total is not a
    /// frequency, it is a table whose numerator and denominator came from
    /// different places. That happens: two files joined on a date, a
    /// denominator counted over a week and a numerator over a day.
    ///
    /// The coordinates on a [`Map`](crate::Map) answer the same question the
    /// same way, and for the same reason. Nothing here repairs the numbers,
    /// because a repaired number is a claim nobody made.
    pub fn is_drawable(&self) -> bool {
        self.total > 0 && self.count <= self.total
    }

    /// Returns the observed frequency, when the two numbers are one.
    ///
    /// `None` rather than a number above one, which the frequency axis has no
    /// room for: it runs from nought to one, so a frequency of two and a half
    /// used to be drawn a hundred pixels above the top of the figure, where
    /// the track's own clip hid it and nothing counted it.
    pub fn frequency(&self) -> Option<f64> {
        self.is_drawable()
            .then_some(self.count as f64 / self.total as f64)
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

    /// How many observations this track cannot place, and so does not draw.
    ///
    /// The count depends on the metric, unlike the equivalent on a
    /// [`Map`](crate::Map), because the two numbers are read differently:
    /// against a frequency axis a count above its total is impossible, and
    /// counted as a count it is merely a count. Only the first is undrawable.
    ///
    /// This does not include an observation outside the drawn region, or one
    /// held back by [`minimum_total`](SurveillanceTrack::minimum_total). Those
    /// are choices about what to show. This is data that cannot be shown.
    pub fn undrawable_observation_count(&self) -> usize {
        match self.metric {
            SurveillanceMetric::Frequency => self
                .observations
                .iter()
                .filter(|observation| !observation.is_drawable())
                .count(),
            SurveillanceMetric::Count => 0,
        }
    }

    fn valid_observations(&self, ctx: &DrawContext<'_>) -> Vec<&SurveillanceObservation> {
        self.observations
            .iter()
            .filter(|observation| ctx.region.contains(observation.time))
            .filter(|observation| observation.total >= self.minimum_total)
            .filter(|observation| match self.metric {
                SurveillanceMetric::Frequency => observation.is_drawable(),
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
        let header = ctx.px(21.0).min(ctx.band.h * 0.25);
        let top = ctx.band.y + header;

        // An observation that cannot be read as a frequency never reaches the
        // trajectories below, because `valid_observations` drops it. Dropping it
        // quietly is the part that was wrong: it used to be drawn at a
        // frequency above one, which on an axis that stops at one put it above
        // the top of the figure, where this track's own clip hid it. Nobody
        // counts the dots on a plot against the rows of their table.
        //
        // It is named here, before the panel decides whether it has anything to
        // plot, because the panel used to return on an empty trajectory list
        // without a word. A table whose every row is unreadable is the table
        // most in need of being told so, and it was the one table that was told
        // nothing at all: the count was available to a caller who thought to
        // ask, and the document said the panel simply had no data.
        //
        // The muted diamond is the same one the unjoinable and unstackable rows
        // get, at the time it was observed, carrying its own numbers. That is
        // this file's existing way of saying "here is something I did not
        // draw", and a second way of saying it would be one more thing to keep
        // right.
        if self.metric == SurveillanceMetric::Frequency {
            for observation in self
                .observations
                .iter()
                .filter(|observation| ctx.region.contains(observation.time))
                .filter(|observation| !observation.is_drawable())
            {
                let reason = if observation.total == 0 {
                    "no denominator, so there is no frequency"
                } else {
                    "a count above its total is not a frequency"
                };
                name_undrawn(ctx, reason, observation, top, 0);
            }
        }

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

        // The single observation of a lineage at a time, where there is one.
        // Two rows for one lineage at one time are two answers to one question
        // and the panel joins neither, so this is `None` there. The row is kept
        // rather than only its value because the growth baseline further down
        // needs to know which rows the panel actually joined.
        let joined: Vec<Vec<Option<&SurveillanceObservation>>> = lineages
            .iter()
            .map(|lineage| {
                times
                    .iter()
                    .map(|time| {
                        let mut matching = observations.iter().filter(|observation| {
                            observation.time == *time && observation.lineage == *lineage
                        });
                        match (matching.next(), matching.next()) {
                            (Some(only), None) => Some(*only),
                            _ => None,
                        }
                    })
                    .collect()
            })
            .collect();
        let matrix: Vec<Vec<Option<f64>>> = joined
            .iter()
            .map(|series| {
                series
                    .iter()
                    .map(|observation| observation.and_then(|row| self.value(row)))
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
            for observation in &history {
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
                    name_undrawn(ctx, reason, observation, top, lineage_index);
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
                // The baseline is the last value this lineage has, not the last
                // row it has. Where two rows share a lineage and a time the
                // panel joins neither and says so, yet the change used to be
                // measured against whichever of the two the file happened to
                // list second: the same four rows sorted the other way reported
                // a rise of three quarters instead of a twentieth, and raised a
                // growth alert for it. A row the panel would not join is not a
                // number that came before.
                let previous_frequency = (0..time_index)
                    .rev()
                    .find_map(|earlier| joined[lineage_index][earlier])
                    .and_then(|previous| previous.frequency());
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

/// Radius of the muted diamond that names a row the panel did not draw.
const NOTICE_RADIUS: f64 = 2.7;

/// Width of that diamond's ring, which is drawn outside the radius.
const NOTICE_RING: f64 = 0.8;

/// Names one observation the panel did not draw, and the reason it gives.
///
/// Every refusal goes through here so that the mark, its ring and the row it
/// is kept on cannot drift apart from one another.
///
/// `step` separates several refusals at one time, one row per lineage, and it
/// stops at the last row that still holds the whole mark. It used to step on
/// without a stop: twenty lineages on the shortest band put diamonds below the
/// band, and thirty put every one of them there, where the track's own clip
/// removed them. A notice that is itself undrawn is exactly the loss the
/// notice exists to prevent, so it stops for the same reason the legend beside
/// it stops when it runs out of width.
fn name_undrawn(
    ctx: &mut DrawContext<'_>,
    reason: &str,
    observation: &SurveillanceObservation,
    top: f64,
    step: usize,
) {
    let radius = ctx.px(NOTICE_RADIUS);
    let ring = ctx.px(NOTICE_RING);
    let row = (top + ctx.px(5.0 + step as f64 * 1.7)).min(ctx.band.bottom() - radius - ring);
    ctx.svg.begin_titled(&format!(
        "{reason} | {}",
        observation_title(observation, None, None, None)
    ));
    ctx.svg.symbol_ringed(
        ctx.scale.x_center(observation.time),
        row,
        radius,
        Symbol::Diamond,
        &ctx.theme.muted,
        ctx.theme.surface(),
        ring,
    );
    ctx.svg.end_group();
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
    fn a_count_above_its_total_is_counted_rather_than_drawn_off_the_page() {
        // A part cannot exceed its whole, and a table where the numerator and
        // the denominator came from different files can still say it does. The
        // frequency axis runs from nought to one, so 250 out of 100 used to be
        // placed a hundred pixels above the top of the whole figure, where this
        // track's own clip hid it and nothing counted it.
        let track = SurveillanceTrack::new(vec![
            SurveillanceObservation::new(10, "A", 5, 100),
            SurveillanceObservation::new(20, "A", 250, 100),
            SurveillanceObservation::new(25, "A", 7, 0),
            SurveillanceObservation::new(30, "A", 40, 100),
        ])
        .style(SurveillanceStyle::Lines);

        assert_eq!(track.undrawable_observation_count(), 2);
        assert_eq!(
            SurveillanceObservation::new(0, "A", 250, 100).frequency(),
            None
        );

        let figure = Figure::new(Region::new("time", 0, 40).unwrap()).push(track);
        let (_, height) = figure.dimensions();
        let svg = figure.to_svg();

        // Nothing is placed outside the figure it belongs to.
        for piece in svg.split("cy=\"").skip(1) {
            let y: f64 = piece
                .split('"')
                .next()
                .and_then(|value| value.parse().ok())
                .expect("every cy is a number");
            assert!(
                (0.0..=height).contains(&y),
                "a mark sits at {y} on a figure {height} tall"
            );
        }

        // And what it did not draw, it says.
        assert!(
            svg.contains("a count above its total is not a frequency"),
            "{svg}"
        );
        assert!(
            svg.contains("no denominator, so there is no frequency"),
            "{svg}"
        );
    }

    #[test]
    fn a_panel_that_can_draw_none_of_its_rows_still_names_every_one_of_them() {
        // The table most in need of being told it could not be read is the one
        // where every row is unreadable, and that was the one table told
        // nothing: with no trajectory left to plot, the panel returned before it
        // reached its own notices, and the document came out looking exactly
        // like a panel that had been handed no data at all. The count was there
        // for a caller who thought to ask, which is not the same as saying it.
        let track = SurveillanceTrack::new(vec![
            SurveillanceObservation::new(1, "A", 250, 100),
            SurveillanceObservation::new(2, "B", 4, 0),
        ]);
        assert_eq!(track.undrawable_observation_count(), 2);

        let svg = Figure::new(Region::new("week", 0, 4).unwrap())
            .push(track)
            .to_svg();
        assert!(
            svg.contains(
                "a count above its total is not a frequency | A | time 1 | count 250 of 100"
            ),
            "{svg}"
        );
        assert!(
            svg.contains("no denominator, so there is no frequency | B | time 2 | count 4 of 0"),
            "{svg}"
        );
    }

    #[test]
    fn the_diamond_that_names_an_undrawn_row_stays_inside_the_band() {
        // The notices step down one row per lineage so that several at one time
        // do not land on each other, and the step used to run on without a stop.
        // Twenty lineages on the shortest band put the last of them below the
        // band, where the figure's own clip removes them: a notice about
        // something that was not drawn, itself not drawn.
        let observations: Vec<SurveillanceObservation> = (0..20u64)
            .map(|lineage| SurveillanceObservation::new(lineage + 1, format!("L{lineage}"), 1, 2))
            .collect();
        let svg = Figure::new(Region::new("week", 0, 40).unwrap())
            .push(SurveillanceTrack::new(observations).height(48.0))
            .to_svg();
        assert!(
            svg.contains("incomplete lineage composition is not stacked"),
            "{svg}"
        );

        // Measured against the clip the document itself carries, since that is
        // what a renderer crops to, rather than a rectangle worked out here.
        let clip = svg
            .split("<clipPath")
            .nth(1)
            .and_then(|inside| inside.split("/>").next())
            .expect("one clip for one track");
        let read = |name: &str| -> f64 {
            clip.split(&format!("{name}=\""))
                .nth(1)
                .and_then(|piece| piece.split('"').next())
                .and_then(|value| value.parse().ok())
                .expect("a number on the clipping rectangle")
        };
        let top = read("y");
        let bottom = top + read("height");
        for piece in svg.split("<polygon points=\"").skip(1) {
            let list = piece.split('"').next().expect("a point list");
            for point in list.split_whitespace() {
                let row: f64 = point
                    .split(',')
                    .nth(1)
                    .and_then(|value| value.parse().ok())
                    .expect("every vertex is a pair of numbers");
                // The slack is the rounding the writer applies to a coordinate,
                // not room for a mark to sit outside its band.
                assert!(
                    row >= top - 0.001 && row <= bottom + 0.001,
                    "a diamond vertex sits at {row} on a band running {top} to {bottom}"
                );
            }
        }
    }

    #[test]
    fn a_change_is_not_measured_against_a_row_the_panel_will_not_join() {
        // Two rows for one lineage at one time are two answers to one question,
        // and the panel says beside both that it joins neither. The change on
        // the next observation was still measured against whichever of the two
        // the file listed second, so the same four rows sorted the other way
        // reported a different rise and raised a growth alert for it. A time
        // the panel has no single value for is a gap, and a gap is measured
        // across, exactly as it already was for a lineage simply missing a row.
        let panel = |second: u64, third: u64| {
            Figure::new(Region::new("week", 0, 5).unwrap())
                .push(
                    SurveillanceTrack::new(vec![
                        SurveillanceObservation::new(1, "A", 90, 100),
                        SurveillanceObservation::new(2, "A", second, 100),
                        SurveillanceObservation::new(2, "A", third, 100),
                        SurveillanceObservation::new(3, "A", 95, 100),
                    ])
                    .style(SurveillanceStyle::Lines)
                    .growth_alert(0.5),
                )
                .to_svg()
        };
        let one = panel(10, 20);
        let other = panel(20, 10);
        // The two rows keep the order they were given, since each is named with
        // its own numbers. Only the observation after them is read here.
        let change = |svg: &str| -> String {
            svg.split("<title>")
                .find(|piece| piece.starts_with("A | time 3"))
                .and_then(|piece| piece.split("</title>").next())
                .expect("the last observation is named")
                .to_string()
        };
        assert_eq!(
            change(&one),
            change(&other),
            "swapping two rows the panel will not join changed the change after them"
        );
        assert_eq!(
            change(&one),
            "A | time 3 | count 95 of 100 | frequency 0.95 | change +0.05"
        );
        assert!(!one.contains("growth alert"), "{one}");
        assert!(!other.contains("growth alert"), "{other}");
    }

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
        // The two omissions are not the same kind, and the difference is the
        // point of the name. A denominator of nought is data that cannot be
        // read as a frequency, so it is drawn as a muted diamond carrying its
        // own numbers and its reason: it used to vanish, which is what this
        // test asserted before the numbers it could not use were counted.
        //
        // A total below the sampling floor is the caller's own filter. Nothing
        // is wrong with those numbers and nothing is hidden from anyone who
        // set the floor, so that one stays out of the document.
        assert!(
            svg.contains("no denominator, so there is no frequency | zero | time 1"),
            "{svg}"
        );
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
