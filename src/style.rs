//! Reusable visual contracts shared by figures, panels, rings and tracks.
//!
//! A theme answers *which* colours and typefaces a drawing uses. This module
//! answers *how strongly* and *how densely* it uses them: the same named
//! profile, emphasis level and quantitative axis mean the same thing in every
//! representation.

/// A complete, named output treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RenderProfile {
    /// Dense multi-panel figures and narrow document columns.
    Compact,
    /// Balanced defaults for papers and vector editing.
    #[default]
    Manuscript,
    /// Larger type and marks for projected slides.
    Presentation,
    /// Slightly larger interactive figures for documentation and browsers.
    Web,
    /// Manuscript geometry on a selected dark palette.
    Dark,
}

impl RenderProfile {
    /// Scale applied to typography, marks and surrounding whitespace.
    pub fn visual_scale(self) -> f64 {
        match self {
            RenderProfile::Compact => 0.88,
            RenderProfile::Manuscript | RenderProfile::Dark => 1.0,
            RenderProfile::Presentation => 1.35,
            RenderProfile::Web => 1.10,
        }
    }

    /// Density of the data-bearing bands.
    pub fn density(self) -> Density {
        match self {
            RenderProfile::Compact => Density::Compact,
            RenderProfile::Presentation => Density::Spacious,
            _ => Density::Balanced,
        }
    }

    /// Whether the profile starts from the dark palette.
    pub fn is_dark(self) -> bool {
        matches!(self, RenderProfile::Dark)
    }
}

/// How tightly repeated rows and marks are packed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Density {
    /// More data in the same vertical space.
    Compact,
    /// The default balance between scanning and packing.
    #[default]
    Balanced,
    /// More separation for presentations and touch-sized figures.
    Spacious,
}

impl Density {
    /// Scale applied to track-internal geometry.
    pub fn scale(self) -> f64 {
        match self {
            Density::Compact => 0.84,
            Density::Balanced => 1.0,
            Density::Spacious => 1.14,
        }
    }
}

/// Semantic prominence, independent of a concrete colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Emphasis {
    /// Context that should recede without disappearing.
    Muted,
    /// Ordinary data marks.
    #[default]
    Normal,
    /// The main series or locus the figure is about.
    Primary,
    /// A threshold crossing, warning or selected event.
    Alert,
}

/// A line pattern that remains distinguishable without colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LinePattern {
    /// Continuous stroke.
    #[default]
    Solid,
    /// Repeating dashes.
    Dashed,
    /// Short separated dots.
    Dotted,
}

impl LinePattern {
    /// SVG `stroke-dasharray`, or `None` for a solid line.
    pub fn dasharray(self) -> Option<&'static str> {
        match self {
            LinePattern::Solid => None,
            LinePattern::Dashed => Some("6 4"),
            LinePattern::Dotted => Some("1.5 3"),
        }
    }
}

/// Point shape, cycled alongside colour for categorical data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Symbol {
    /// Circular point.
    #[default]
    Circle,
    /// Square point.
    Square,
    /// Diamond point.
    Diamond,
    /// Upward triangle.
    Triangle,
}

/// Measured visual sizes used throughout a drawing.
#[derive(Debug, Clone, PartialEq)]
pub struct VisualTokens {
    /// Hairline guides and row separators.
    pub hairline: f64,
    /// Ordinary data stroke.
    pub stroke: f64,
    /// Selected or primary data stroke.
    pub strong_stroke: f64,
    /// Default point radius.
    pub marker_radius: f64,
    /// Length of an axis tick.
    pub tick_length: f64,
    /// Air between a label and the mark or axis it names.
    pub label_gap: f64,
    /// Default vertical gap between repeated rows.
    pub row_gap: f64,
    /// Default height of a gene or interval mark.
    pub feature_height: f64,
    /// Width of an arrow head.
    pub arrow_size: f64,
    /// Size of a legend swatch.
    pub legend_swatch: f64,
    /// Horizontal gap between legend entries.
    pub legend_gap: f64,
    /// Default opacity of an area placed under a line.
    pub area_opacity: f64,
}

impl VisualTokens {
    /// Returns the same tokens scaled as one system.
    pub fn scaled(mut self, factor: f64) -> Self {
        let factor = if factor.is_finite() {
            factor.max(0.25)
        } else {
            1.0
        };
        self.hairline *= factor;
        self.stroke *= factor;
        self.strong_stroke *= factor;
        self.marker_radius *= factor;
        self.tick_length *= factor;
        self.label_gap *= factor;
        self.row_gap *= factor;
        self.feature_height *= factor;
        self.arrow_size *= factor;
        self.legend_swatch *= factor;
        self.legend_gap *= factor;
        self
    }
}

impl Default for VisualTokens {
    fn default() -> Self {
        VisualTokens {
            hairline: 0.8,
            stroke: 1.2,
            strong_stroke: 2.0,
            marker_radius: 3.0,
            tick_length: 5.0,
            label_gap: 10.0,
            row_gap: 3.0,
            feature_height: 14.0,
            arrow_size: 8.0,
            legend_swatch: 10.0,
            legend_gap: 16.0,
            area_opacity: 0.18,
        }
    }
}

/// Concrete mark treatment derived from an [`Emphasis`] level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MarkStyle {
    /// Stroke width in pixels.
    pub stroke_width: f64,
    /// Point radius in pixels.
    pub marker_radius: f64,
    /// Opacity between zero and one.
    pub opacity: f64,
    /// Pattern used by line marks.
    pub pattern: LinePattern,
}

/// Formatting of quantitative tick labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AxisFormat {
    /// Compact values with automatic SI-like `k` and `M` suffixes.
    #[default]
    Auto,
    /// Fixed number of decimal places.
    Fixed(u8),
    /// Fraction shown as a percentage.
    Percent(u8),
}

impl AxisFormat {
    /// Formats one finite value for an axis.
    pub fn format(self, value: f64) -> String {
        if !value.is_finite() {
            return "0".to_string();
        }
        match self {
            AxisFormat::Auto => compact_number(value),
            AxisFormat::Fixed(places) => format!("{value:.places$}", places = places as usize),
            AxisFormat::Percent(places) => {
                format!("{:.places$}%", value * 100.0, places = places as usize)
            }
        }
    }
}

fn compact_number(value: f64) -> String {
    let (number, suffix) = if value.abs() >= 1_000_000.0 {
        (value / 1_000_000.0, "M")
    } else if value.abs() >= 1_000.0 {
        (value / 1_000.0, "k")
    } else {
        (value, "")
    };
    // Past the millions the suffixes run out, and a number too large to hold
    // every integer is written in exponent form rather than cast to one: the
    // cast saturates, so the largest double used to come out as the largest
    // i64 with an M after it.
    if number.abs() >= 1e15 {
        return format!("{value:e}");
    }
    let rounded = (number * 100.0).round() / 100.0;
    let text = if rounded == rounded.trunc() {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded:.2}").trim_end_matches('0').to_string()
    };
    format!("{text}{suffix}")
}

/// A labelled line drawn across a quantitative track.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceLine {
    /// Value at which the line is drawn.
    pub value: f64,
    /// Optional label placed beside its tick.
    pub label: Option<String>,
    /// Semantic prominence of the line.
    pub emphasis: Emphasis,
    /// Pattern that distinguishes it without colour.
    pub pattern: LinePattern,
}

impl ReferenceLine {
    /// An unlabelled reference line.
    pub fn new(value: f64) -> Self {
        ReferenceLine {
            value,
            label: None,
            emphasis: Emphasis::Muted,
            pattern: LinePattern::Dashed,
        }
    }

    /// Sets the label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets its prominence.
    pub fn emphasis(mut self, emphasis: Emphasis) -> Self {
        self.emphasis = emphasis;
        self
    }

    /// Sets the stroke pattern.
    pub fn pattern(mut self, pattern: LinePattern) -> Self {
        self.pattern = pattern;
        self
    }
}

/// Shared configuration for quantitative axes.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantitativeAxis {
    /// Optional lower bound.
    pub min: Option<f64>,
    /// Optional upper bound.
    pub max: Option<f64>,
    /// Approximate number of labelled ticks, including the ends.
    pub ticks: usize,
    /// Suffix placed after every tick.
    pub unit: String,
    /// Tick formatting.
    pub format: AxisFormat,
    /// Additional labelled or patterned lines.
    pub references: Vec<ReferenceLine>,
}

impl QuantitativeAxis {
    /// An automatic axis with three ticks.
    pub fn new() -> Self {
        Self::default()
    }

    /// Pins both ends of the range for comparison across plots.
    pub fn range(mut self, min: f64, max: f64) -> Self {
        if min.is_finite() && max.is_finite() {
            self.min = Some(min.min(max));
            self.max = Some(min.max(max));
        }
        self
    }

    /// Sets the approximate number of labelled ticks.
    pub fn ticks(mut self, ticks: usize) -> Self {
        self.ticks = ticks.clamp(2, 8);
        self
    }

    /// Sets the unit suffix.
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self
    }

    /// Sets tick formatting.
    pub fn format(mut self, format: AxisFormat) -> Self {
        self.format = format;
        self
    }

    /// Adds a reference line.
    pub fn reference(mut self, reference: ReferenceLine) -> Self {
        self.references.push(reference);
        self
    }

    /// Formats a value together with the unit.
    pub fn label(&self, value: f64) -> String {
        format!("{}{}", self.format.format(value), self.unit)
    }

    /// The labels for a column of ticks, with the unit written once.
    ///
    /// A unit repeated on every tick is the same word stacked three times
    /// beside a band a few lines tall, and it is most of why a narrow value
    /// axis reads as clutter. It goes on the highest tick, which is the one
    /// the eye reaches first and the one a reader asks "of what" about.
    pub fn labels(&self, values: &[f64]) -> Vec<String> {
        let top = values
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        values
            .iter()
            .map(|&value| {
                if value == top {
                    self.label(value)
                } else {
                    self.format.format(value)
                }
            })
            .collect()
    }

    /// Widens the ends nobody pinned out to round values.
    ///
    /// A ceiling of 71.46 labelled as 71.46 is a number nobody asked for: it
    /// is whatever the tallest sample happened to be, and it makes the middle
    /// tick 35.73. Rounding the free ends outwards to the step the ticks will
    /// use puts every label on a value a reader would have chosen, and the
    /// small overshoot is the headroom a profile wants anyway. A pinned end is
    /// taken literally, because pinning one is how two panels are made to
    /// agree.
    ///
    /// The rounding may use a finer step than the labels will, up to three
    /// more intervals than [`QuantitativeAxis::ticks`], because the step that
    /// suits three labels can overshoot badly: 10.2 rounded in steps of five
    /// is 15, a third of the band spent above the data, where steps of two
    /// stop at 12.
    pub fn nice(&self, min: f64, max: f64) -> (f64, f64) {
        if !(min.is_finite() && max.is_finite()) || max <= min {
            return (min, max);
        }
        let pin_lo = self.min.is_some();
        let pin_hi = self.max.is_some();
        if pin_lo && pin_hi {
            return (min, max);
        }
        let most = self.ticks.clamp(2, 8) as f64 + 3.0;
        let mut best: Option<(f64, f64, f64)> = None;
        for step in candidate_steps(max - min, most) {
            let lo = if pin_lo { min } else { floor_to(min, step) };
            let hi = if pin_hi { max } else { ceil_to(max, step) };
            if (hi - lo) / step > most + 1e-9 {
                continue;
            }
            // The narrowest range wins, so the data fills as much of the band
            // as a round ceiling allows; between equals, the finer step.
            if best.map_or(true, |(span, _, _)| hi - lo < span - 1e-9 * span.abs()) {
                best = Some((hi - lo, lo, hi));
            }
        }
        let Some((_, mut lo, mut hi)) = best else {
            return (min, max);
        };
        // Data that already sits on a round value would touch the edge of the
        // band, so the band is taken a little past the last tick instead. The
        // tick stays where it is; only the room above it grows. Zero is left
        // alone, because a count that starts at zero starts on the baseline.
        let room = (hi - lo) * 0.04;
        if !pin_hi && hi - max < (hi - lo) * 0.03 {
            hi += room;
        }
        if !pin_lo && lo != 0.0 && min - lo < (hi - lo) * 0.03 {
            lo -= room;
        }
        (lo, hi)
    }

    /// Round tick values inside `min..=max`.
    ///
    /// Multiples of a step of 1, 2, 2.5 or 5 times a power of ten, the finest
    /// one that keeps the count at most one more than [`QuantitativeAxis::ticks`].
    /// A range too narrow to hold two of them gets its own two ends instead.
    pub fn values(&self, min: f64, max: f64) -> Vec<f64> {
        if !(min.is_finite() && max.is_finite()) {
            return Vec::new();
        }
        if max <= min {
            return vec![min];
        }
        let most = self.ticks.clamp(2, 8) as f64 + 1.0;
        // The step whose ticks reach closest to both ends, so a range that was
        // rounded by `nice` gets labels on its own ends rather than a finer
        // step that happens to stop short of the ceiling.
        let mut best: Option<(f64, f64, f64, f64)> = None;
        for step in candidate_steps(max - min, most - 1.0) {
            let first = ((min / step) - 1e-9).ceil();
            let last = ((max / step) + 1e-9).floor();
            let count = last - first + 1.0;
            if count > most || count < 2.0 {
                continue;
            }
            let short = (first * step - min) + (max - last * step);
            // A quarter step only where it lands on both ends, as it does on
            // a range `nice` rounded to it. Anywhere else it is how a span of
            // years comes to be labelled 2022.5.
            if is_quarter_step(step) && short > 1e-9 * (max - min) {
                continue;
            }
            // Between two that reach the ends equally well, the one whose
            // count is nearer the number asked for, and between those the
            // coarser: whole years over half years, where both would do.
            let wanted = self.ticks.clamp(2, 8) as f64;
            let tolerance = 1e-9 * (max - min);
            let better = best.map_or(true, |(gap, held, held_count, _)| {
                if short < gap - tolerance {
                    return true;
                }
                if short > gap + tolerance {
                    return false;
                }
                let (near, held_near) = ((count - wanted).abs(), (held_count - wanted).abs());
                near < held_near || (near == held_near && step > held)
            });
            if better {
                best = Some((short, step, count, first));
            }
        }
        let Some((_, step, count, first)) = best else {
            return vec![min, max];
        };
        let places = decimals_of(step);
        let first = first as i64;
        (first..first + count as i64)
            .map(|k| tidy(k as f64 * step, places))
            .collect()
    }

    /// The width the widest of this axis's labels needs at `font_size`.
    pub(crate) fn label_room(&self, min: f64, max: f64, font_size: f64) -> f64 {
        let values = self.values(min, max);
        self.labels(&values)
            .iter()
            .map(|label| crate::svg::text_width(label, font_size))
            .fold(0.0f64, f64::max)
    }

    /// Resolves optional pinned ends against a data range and keeps a visible
    /// span even when all values are identical.
    pub fn resolve(&self, data_min: f64, data_max: f64) -> (f64, f64) {
        let mut min = self.min.filter(|v| v.is_finite()).unwrap_or(data_min);
        let mut max = self.max.filter(|v| v.is_finite()).unwrap_or(data_max);
        if min > max {
            std::mem::swap(&mut min, &mut max);
        }
        if !min.is_finite() || !max.is_finite() {
            return (0.0, 1.0);
        }
        if (max - min).abs() <= f64::EPSILON {
            let pad = min.abs().max(1.0) * 0.5;
            return (min - pad, max + pad);
        }
        (min, max)
    }
}

/// Round values across a span of time, each with its label.
///
/// Written as plain numbers with as many decimals as the step needs, never
/// with a `k`: a year is 2024 and not 2.02k. The unit goes on the latest.
/// The two ends of a dated tree are its root and its latest tip, and they are
/// rarely round, so labelled as they came a tree read "2021.85", "2022.965"
/// and "2024.08", three numbers nobody would put on an axis.
pub(crate) fn time_ticks(min: f64, max: f64, unit: Option<&str>) -> Vec<(f64, String)> {
    let ticks = QuantitativeAxis::new().ticks(4).values(min, max);
    let step = ticks
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .fold(f64::INFINITY, f64::min);
    let decimals = if step.is_finite() && step < 1.0 {
        (-(step.log10().floor())).clamp(0.0, 6.0) as u32 + u32::from(is_quarter_step(step))
    } else {
        0
    };
    let last = ticks.len().saturating_sub(1);
    ticks
        .iter()
        .enumerate()
        .map(|(index, &value)| {
            let number = crate::svg::text_rounded(value, decimals);
            let label = match unit {
                Some(unit) if index == last => format!("{number} {unit}"),
                _ => number,
            };
            (value, label)
        })
        .collect()
}

/// Steps of 1, 2, 2.5 and 5 times a power of ten, smallest first, starting
/// from the finest that could cover `span` in `intervals` steps.
fn candidate_steps(span: f64, intervals: f64) -> Vec<f64> {
    if !(span.is_finite() && span > 0.0 && intervals > 0.0) {
        return vec![1.0];
    }
    let magnitude = 10f64.powf((span / intervals).log10().floor());
    let mut steps = Vec::with_capacity(12);
    for scale in [1.0, 10.0, 100.0] {
        for multiple in [1.0, 2.0, 2.5, 5.0] {
            steps.push(multiple * magnitude * scale);
        }
    }
    steps
}

/// Whether `step` is 2.5 times a power of ten.
fn is_quarter_step(step: f64) -> bool {
    let magnitude = 10f64.powf(step.log10().floor());
    ((step / magnitude) - 2.5).abs() < 1e-9
}

fn floor_to(value: f64, step: f64) -> f64 {
    tidy(((value / step) + 1e-9).floor() * step, decimals_of(step))
}

fn ceil_to(value: f64, step: f64) -> f64 {
    tidy(((value / step) - 1e-9).ceil() * step, decimals_of(step))
}

/// Decimal places a multiple of `step` can need: one more than the step's
/// own order of magnitude, which is what a 2.5 needs.
fn decimals_of(step: f64) -> i32 {
    (-(step.log10().floor()) + 1.0).clamp(0.0, 12.0) as i32
}

/// A multiple of a step, with the float noise of the multiplication removed
/// and without the sign a zero can pick up on the way.
fn tidy(value: f64, places: i32) -> f64 {
    let factor = 10f64.powi(places);
    let rounded = (value * factor).round() / factor;
    if rounded == 0.0 {
        0.0
    } else {
        rounded
    }
}

/// The ticks from `values` that can be labelled without two labels touching.
///
/// A value axis beside a band twenty pixels tall cannot print three labels,
/// and printing them anyway stacks them into one smudge, which is what the
/// variant tracks were doing. Labels closer than a line and a half of
/// `font_size` are legible and still read as a ladder of numbers rather than
/// a scale.
///
/// The ticks are thinned by a regular stride, every second or fifth one,
/// counted from nought when nought is among them, so what survives is still a
/// round sequence: dropping labels one at a time from the bottom up left a
/// scale reading -1, -0.2 and 1. When not even two fit, the highest is kept on
/// its own, since a band with one number beside it still says how tall it is.
pub(crate) fn legible_ticks(values: &[f64], y_of: impl Fn(f64) -> f64, font_size: f64) -> Vec<f64> {
    let gap = font_size * 1.6;
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    sorted.sort_by(|a, b| a.total_cmp(b));
    sorted.dedup();
    let Some(&top) = sorted.last() else {
        return Vec::new();
    };
    if sorted.len() == 1 {
        return sorted;
    }
    let anchor = sorted.iter().position(|v| *v == 0.0).unwrap_or(0);
    // The first stride that clears and still reaches the highest tick wins,
    // since that one carries the unit. Failing that, the two ends, with nought
    // between them when it is inside: a short band read as "0 to 0.75" says
    // more than one read as "0 to 0.5" that stops short of its own top.
    let clears = |kept: &[f64]| {
        kept.windows(2)
            .all(|pair| (y_of(pair[0]) - y_of(pair[1])).abs() >= gap)
    };
    for stride in [1, 2, 5, 10, 20, 50, 100] {
        let kept: Vec<f64> = sorted
            .iter()
            .enumerate()
            .filter(|(index, _)| (*index as i64 - anchor as i64).rem_euclid(stride) == 0)
            .map(|(_, value)| *value)
            .collect();
        if kept.len() < 2 {
            break;
        }
        if kept.last() == Some(&top) && clears(&kept) {
            return kept;
        }
    }
    let bottom = sorted[0];
    if anchor > 0 && anchor + 1 < sorted.len() && clears(&[bottom, 0.0, top]) {
        return vec![bottom, 0.0, top];
    }
    if clears(&[bottom, top]) {
        return vec![bottom, top];
    }
    vec![top]
}

impl Default for QuantitativeAxis {
    fn default() -> Self {
        QuantitativeAxis {
            min: None,
            max: None,
            ticks: 3,
            unit: String::new(),
            format: AxisFormat::Auto,
            references: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profiles_name_distinct_scales_and_densities() {
        assert!(RenderProfile::Compact.visual_scale() < 1.0);
        assert!(RenderProfile::Presentation.visual_scale() > 1.0);
        assert_eq!(RenderProfile::Compact.density(), Density::Compact);
        assert!(RenderProfile::Dark.is_dark());
    }

    #[test]
    fn axis_formatting_is_compact_and_explicit() {
        assert_eq!(AxisFormat::Auto.format(1_250.0), "1.25k");
        assert_eq!(AxisFormat::Auto.format(f64::MAX), "1.7976931348623157e308");
        assert_eq!(AxisFormat::Auto.format(-1e300), "-1e300");
        assert_eq!(AxisFormat::Fixed(2).format(1.0), "1.00");
        assert_eq!(AxisFormat::Percent(1).format(0.125), "12.5%");
    }

    #[test]
    fn a_free_ceiling_is_rounded_up_to_a_value_someone_would_choose() {
        let axis = QuantitativeAxis::new();
        assert_eq!(axis.nice(0.0, 71.46), (0.0, 75.0));
        assert_eq!(axis.values(0.0, 75.0), vec![0.0, 25.0, 50.0, 75.0]);
        // On a round value already: the band goes a little past the tick.
        assert_eq!(axis.nice(0.0, 0.98), (0.0, 1.04));
        assert_eq!(axis.values(0.0, 1.04), vec![0.0, 0.5, 1.0]);
        assert_eq!(axis.nice(0.0, 449.0), (0.0, 500.0));
        // Within three per cent of the round ends, so the band goes past them.
        assert_eq!(axis.nice(-3.84, 3.84), (-4.32, 4.32));
        assert_eq!(axis.nice(0.0, 10.19), (0.0, 12.0));
        assert_eq!(axis.values(0.0, 12.0), vec![0.0, 5.0, 10.0]);
        assert_eq!(axis.values(0.0, 0.3), vec![0.0, 0.1, 0.2, 0.3]);
        // Whole years stay whole: no quarter step that misses both ends.
        assert_eq!(axis.values(2021.0, 2025.0), vec![2022.0, 2024.0]);
        assert_eq!(axis.values(2021.85, 2024.08), vec![2022.0, 2023.0, 2024.0]);
    }

    #[test]
    fn a_pinned_end_is_left_where_it_was_put() {
        let axis = QuantitativeAxis::new().range(0.0, 72.0);
        assert_eq!(axis.nice(0.0, 72.0), (0.0, 72.0));
        let top_only = QuantitativeAxis {
            max: Some(1.0),
            ..QuantitativeAxis::new()
        };
        assert_eq!(top_only.nice(-0.3, 1.0), (-0.5, 1.0));
        assert_eq!(axis.values(0.0, 72.0), vec![0.0, 20.0, 40.0, 60.0]);
    }

    #[test]
    fn dates_are_labelled_as_plain_round_numbers() {
        let labels: Vec<String> = time_ticks(2021.85, 2024.08, Some("year"))
            .into_iter()
            .map(|(_, label)| label)
            .collect();
        assert_eq!(labels, vec!["2022", "2023", "2024 year"]);
        let halves: Vec<String> = time_ticks(2021.9, 2024.1, None)
            .into_iter()
            .map(|(_, label)| label)
            .collect();
        assert_eq!(halves, vec!["2022", "2023", "2024"]);
        let fine: Vec<String> = time_ticks(0.0, 0.3, None)
            .into_iter()
            .map(|(_, label)| label)
            .collect();
        assert_eq!(fine, vec!["0", "0.1", "0.2", "0.3"]);
    }

    #[test]
    fn the_unit_is_written_once_on_the_highest_tick() {
        let axis = QuantitativeAxis::new().unit(" -log10 p");
        assert_eq!(
            axis.labels(&[0.0, 5.0, 10.0]),
            vec!["0", "5", "10 -log10 p"]
        );
    }

    #[test]
    fn labels_that_would_touch_are_thinned_but_the_top_one_stays() {
        let y_of = |v: f64| 100.0 - v * 20.0;
        // Ten pixels between neighbours, twenty between the ends.
        assert_eq!(legible_ticks(&[0.0, 0.5, 1.0], y_of, 7.0), vec![0.0, 1.0]);
        assert_eq!(legible_ticks(&[0.0, 0.5, 1.0], y_of, 14.0), vec![1.0]);
        assert_eq!(
            legible_ticks(&[0.0, 0.5, 1.0], y_of, 5.0),
            vec![0.0, 0.5, 1.0]
        );
        // Thinned by a stride counted from nought, not one at a time from the
        // bottom, so a centred scale stays a round one.
        let tenths: Vec<f64> = (-5..=5).map(|k| k as f64 * 0.2).collect();
        let tall = |v: f64| 100.0 - v * 40.0;
        let kept = legible_ticks(&tenths, tall, 7.0);
        assert!(kept.contains(&0.0), "{kept:?}");
        assert_eq!(kept, vec![-1.0, 0.0, 1.0]);
    }

    #[test]
    fn an_axis_range_is_ordered() {
        let axis = QuantitativeAxis::new().range(10.0, -2.0);
        assert_eq!(axis.min, Some(-2.0));
        assert_eq!(axis.max, Some(10.0));
        assert_eq!(axis.values(-2.0, 10.0), vec![0.0, 5.0, 10.0]);
        assert_eq!(QuantitativeAxis::new().resolve(3.0, 3.0), (1.5, 4.5));
    }
}
