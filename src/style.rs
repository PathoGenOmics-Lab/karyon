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

    /// Evenly spaced tick values, including both ends.
    pub fn values(&self, min: f64, max: f64) -> Vec<f64> {
        let ticks = self.ticks.clamp(2, 8);
        (0..ticks)
            .map(|index| min + (max - min) * index as f64 / (ticks - 1) as f64)
            .collect()
    }
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
        assert_eq!(AxisFormat::Fixed(2).format(1.0), "1.00");
        assert_eq!(AxisFormat::Percent(1).format(0.125), "12.5%");
    }

    #[test]
    fn an_axis_range_is_ordered() {
        let axis = QuantitativeAxis::new().range(10.0, -2.0);
        assert_eq!(axis.min, Some(-2.0));
        assert_eq!(axis.max, Some(10.0));
        assert_eq!(axis.values(-2.0, 10.0), vec![-2.0, 4.0, 10.0]);
        assert_eq!(QuantitativeAxis::new().resolve(3.0, 3.0), (1.5, 4.5));
    }
}
