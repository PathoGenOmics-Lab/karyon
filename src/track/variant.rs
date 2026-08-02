//! Point events along the sequence: SNPs, indels, insertion sites, peaks.
//!
//! A [`Variant`] is a position, an optional value and an optional category,
//! which is the least a mark needs in order to be placed, sized and coloured.
//! Nothing about the track is particular to sequence variation: any event that
//! happens at one base and may carry a magnitude plots the same way.
//!
//! # Lollipops until they smear, then ticks
//!
//! A lollipop is a stem with a ringed head, and the height of the stem is the
//! value against the top of the track. Density is what decides between the two
//! styles, so [`VariantStyle`] is worth choosing again whenever the number of
//! marks changes by an order of magnitude. The ring is what keeps two marks a
//! base apart reading as two: it is drawn in the background colour, so a head
//! in front of another head still has an edge.
//!
//! # Only lollipops carry a tooltip
//!
//! A lollipop is one mark per call that a pointer can land on, so each one
//! names itself: its position, its category and the number it carried.
//!
//! A tick does not, and the reason is the same one that chose ticks in the
//! first place. Reaching for [`VariantStyle::Tick`] is the caller saying the
//! marks are now too dense to be told apart, which is the statement a binned
//! track makes when it draws one point per pixel column rather than one per
//! datum. A mark a reader cannot isolate is not a mark worth naming, and
//! naming every tick in a genome-wide panel adds about two thirds again to the
//! size of the file for labels nobody can reach.
//!
//! # The axis appears only when the stems are measuring something
//!
//! Ticks ignore their values by design, and a variant with no value gets a full
//! height stem, so in either case a number at the top of the band would be one
//! invented for the picture. The axis is therefore drawn only when several
//! conditions hold at once, and the caller asking for it is only one of them.
//! Pin the ceiling with [`VariantTrack::max`] whenever two panels are meant to
//! carry the same quantity.
//!
//! # A category's colour is its place in the queue
//!
//! Categories take palette slots in order of first appearance, so the colouring
//! is a property of the order the list arrives in rather than of the names in
//! it. Appending leaves the marks already drawn alone; sorting the same
//! variants by position or by frequency hands out different colours. Two
//! figures that have to agree on what red means must be given their variants in
//! one order.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::Theme;
use crate::track::axis::group_thousands;
use crate::track::{DrawContext, Track};

/// How a variant track is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantStyle {
    /// A stem with a head, the height of the stem showing the value. Reads
    /// well up to a few hundred variants.
    Lollipop,
    /// A plain vertical tick, ignoring the value. The right choice once the
    /// variants are dense enough that heads would overlap into a smear.
    Tick,
}

/// A point event at a single base.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    /// Position, 0-based. A VCF line at `POS` becomes `POS - 1`.
    pub pos: u64,
    /// Height of the stem, in the units of the track maximum. `None` draws a
    /// full-height stem, which is what you want when the variant has no
    /// quantity attached.
    pub value: Option<f64>,
    /// Group this variant belongs to, which decides its colour and its legend
    /// entry. Typically the consequence, the lineage or the sample.
    pub category: Option<String>,
}

impl Variant {
    /// A variant at a 0-based position.
    pub fn new(pos: u64) -> Self {
        Variant {
            pos,
            value: None,
            category: None,
        }
    }

    /// Sets the value driving the stem height, for example an allele frequency.
    pub fn value(mut self, value: f64) -> Self {
        self.value = Some(value);
        self
    }

    /// Sets the category driving the colour.
    pub fn category(mut self, category: impl Into<String>) -> Self {
        self.category = Some(category.into());
        self
    }
}

/// Variants drawn as lollipops or ticks, coloured by category.
///
/// Categories are coloured in order of first appearance, so the same input
/// always produces the same figure. That determinism is the point: a figure
/// that recolours itself when a sample is added is not one you can put in a
/// paper.
///
/// ```
/// use karyon::{Figure, Region, Variant, VariantTrack};
///
/// let variants = vec![
///     Variant::new(1000).value(0.9).category("missense"),
///     Variant::new(1500).value(0.3).category("synonymous"),
/// ];
/// let svg = Figure::new(Region::parse("chr1:1-3000").unwrap())
///     .push(VariantTrack::new(variants).label("variants"))
///     .to_svg();
/// assert!(svg.contains("missense"));
/// ```
#[derive(Debug, Clone)]
pub struct VariantTrack {
    variants: Vec<Variant>,
    label: Option<String>,
    height: f64,
    style: VariantStyle,
    radius: f64,
    max: Option<f64>,
    show_legend: bool,
    show_scale: bool,
    color: Option<String>,
}

impl VariantTrack {
    /// A track holding `variants`.
    pub fn new(variants: impl Into<Vec<Variant>>) -> Self {
        VariantTrack {
            variants: variants.into(),
            label: None,
            height: 55.0,
            style: VariantStyle::Lollipop,
            radius: 4.0,
            max: None,
            show_legend: true,
            show_scale: true,
            color: None,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(2.0);
        self
    }

    /// Chooses lollipops or ticks.
    pub fn style(mut self, style: VariantStyle) -> Self {
        self.style = style;
        self
    }

    /// Sets the radius of the lollipop head.
    pub fn radius(mut self, radius: f64) -> Self {
        self.radius = radius.max(0.5);
        self
    }

    /// Pins the value that reaches the full band height.
    ///
    /// Without this the track uses the largest value on screen, which makes
    /// two panels of the same figure silently disagree about scale.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self
    }

    /// Draws or hides the category legend.
    pub fn show_legend(mut self, show: bool) -> Self {
        self.show_legend = show;
        self
    }

    /// Draws or hides the value axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Whether there is a scale worth drawing.
    ///
    /// Three things have to be true. The stems have to mean something, which
    /// rules out [`VariantStyle::Tick`], where every mark is full height by
    /// design. Some variant has to carry a value, or the ceiling of one is an
    /// invention rather than a measurement. And the caller has to want it.
    fn has_scale(&self) -> bool {
        self.show_scale
            && self.style == VariantStyle::Lollipop
            && (self.max.is_some() || self.value_ceiling().is_some())
    }

    /// Sets the colour used for variants without a category.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// The variants in the track.
    pub fn variants(&self) -> &[Variant] {
        &self.variants
    }

    /// Distinct categories, in order of first appearance.
    pub fn categories(&self) -> Vec<&str> {
        let mut seen: Vec<&str> = Vec::new();
        for variant in &self.variants {
            if let Some(category) = variant.category.as_deref() {
                if !seen.contains(&category) {
                    seen.push(category);
                }
            }
        }
        seen
    }

    /// Largest finite value carried by any variant, if any.
    fn value_ceiling(&self) -> Option<f64> {
        self.variants
            .iter()
            .filter_map(|v| v.value)
            .filter(|v| v.is_finite() && *v > 0.0)
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            })
    }
}

impl Track for VariantTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.has_scale() {
            return 0.0;
        }
        let ceiling = self.max.or_else(|| self.value_ceiling()).unwrap_or(1.0);
        text_width(&format_value(ceiling), theme.font_size - 1.0) + 8.0
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
            1.0,
        );

        let categories = self.categories();
        let default_color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());
        let ceiling = self.max.or_else(|| self.value_ceiling()).unwrap_or(1.0);

        // Leave room for the legend so a tall stem does not run through it.
        let legend_room = if self.show_legend && !categories.is_empty() {
            ctx.theme.font_size + 4.0
        } else {
            0.0
        };
        let stem_room = (band.h - legend_room - self.radius).max(2.0);

        if self.has_scale() {
            // The height a variant at the ceiling reaches, which is where the
            // top of the scale is and nowhere else.
            let ceiling_y = baseline - stem_room;
            ctx.svg.line(
                band.x,
                ceiling_y,
                band.right(),
                ceiling_y,
                &ctx.theme.rule,
                1.0,
            );
            if ctx.axis.w > 0.0 {
                let size = ctx.theme.font_size - 1.0;
                let right = ctx.axis.right() - 4.0;
                ctx.svg.text(
                    right,
                    ceiling_y + size * 0.35,
                    &format_value(ceiling),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
                ctx.svg
                    .text(right, baseline, "0", &ctx.theme.muted, size, Anchor::End);
            }
        }

        for variant in &self.variants {
            if !ctx.region.contains(variant.pos) {
                continue;
            }
            let x = ctx.scale.x_center(variant.pos);
            let color = match variant.category.as_deref() {
                Some(category) => {
                    let index = categories.iter().position(|c| *c == category).unwrap_or(0);
                    ctx.theme.color(index).to_string()
                }
                None => default_color.clone(),
            };

            match self.style {
                VariantStyle::Tick => {
                    // Not named, and the style is the reason. Choosing ticks is
                    // the caller saying the marks are too dense to be told
                    // apart, which is the same statement a binned track makes
                    // when it draws one point per pixel column: a mark a
                    // pointer cannot land on alone is not one worth naming.
                    // Naming them anyway costs about two thirds of the file.
                    ctx.svg
                        .line(x, baseline, x, baseline - stem_room, &color, 1.0);
                }
                VariantStyle::Lollipop => {
                    let fraction = match variant.value {
                        Some(value) if value.is_finite() && ceiling > 0.0 => {
                            (value / ceiling).clamp(0.0, 1.0)
                        }
                        _ => 1.0,
                    };
                    let top = baseline - fraction * stem_room;
                    // A lollipop is a stem and a head, and the group is what
                    // makes the two of them one thing. A tooltip on half of a
                    // mark is worse than none.
                    ctx.svg.begin_titled(&tooltip(variant));
                    ctx.svg.line(x, baseline, x, top, &color, 1.4);
                    // The ring is what keeps two variants a base apart reading
                    // as two variants instead of one blob.
                    ctx.svg
                        .circle_ringed(x, top, self.radius, &color, &ctx.theme.background, 1.5);
                    ctx.svg.end_group();
                }
            }
        }

        if self.show_legend && !categories.is_empty() {
            let font = ctx.theme.font_size - 1.0;
            let gap = 12.0;
            let swatch = 4.0;
            let widths: Vec<f64> = categories
                .iter()
                .map(|c| swatch + 4.0 + text_width(c, font))
                .collect();
            let total: f64 = widths.iter().sum::<f64>() + gap * (categories.len() - 1) as f64;
            let mut x = (band.right() - total).max(band.x);
            let y = band.y + font;
            for (index, category) in categories.iter().enumerate() {
                let color = ctx.theme.color(index).to_string();
                ctx.svg
                    .circle(x + swatch / 2.0, y - font * 0.35, swatch / 2.0, &color);
                ctx.svg.text(
                    x + swatch + 4.0,
                    y,
                    category,
                    &ctx.theme.muted,
                    font,
                    Anchor::Start,
                );
                x += widths[index] + gap;
            }
        }
    }
}

/// What a reader hovering one mark is told.
///
/// The noun first, then the position in the 1-based form the ruler under it
/// prints, then what kind of event it is and what number it carried. The last
/// two are optional in the data and neither is invented when it is absent: a
/// mark with nothing but a position says nothing but its noun and a position.
///
/// The noun is there because every tooltip in the crate is `what it is, where
/// it is`, and a mark that opened on a bare `31,218,401` was the one glyph on
/// a sheet answering in a different grammar from the pileup above it.
///
/// The number is named rather than printed bare, because `0.55` on its own is
/// not a statement. What it is named is `value` and not anything more specific,
/// since the track takes any quantity at all: the same band draws allele
/// fractions, peak heights and read counts, and only the caller knows which.
/// The gutter label is where that is said, and it is said once rather than on
/// every mark.
///
/// Only lollipops are named. See the module documentation for why a tick is
/// not.
fn tooltip(variant: &Variant) -> String {
    let mut text = format!(
        "variant, {}",
        group_thousands(variant.pos.saturating_add(1))
    );
    if let Some(category) = variant.category.as_deref() {
        if !category.is_empty() {
            text.push_str(", ");
            text.push_str(category);
        }
    }
    if let Some(value) = variant.value {
        if value.is_finite() {
            text.push_str(", value ");
            text.push_str(&exact_value(value));
        }
    }
    text
}

/// The value written out in full, for the tooltip and nowhere else.
///
/// [`format_value`] is an axis-label abbreviator and belongs on the axis, where
/// the room is a few characters wide and `1.2k` is the only thing that fits. A
/// tooltip has no width constraint at all and is the one place a reader goes
/// for the exact figure, so a read depth of 1,234 that came back as `1.2k` had
/// been rounded by a layout decision made somewhere else on the page.
///
/// The integer part is grouped like every other number in a tooltip, and the
/// fraction is written to a fixed two places rather than to however many the
/// value happens to need. Fixed, because the track takes allele fractions,
/// peak heights and read counts without being told which: a rule that dropped
/// the decimals from whole numbers would put `1` under a pointer and `0.98`
/// under the one beside it, and a reader comparing two marks would be reading
/// two formats. `1.00` and `0.98` line up; `1,234.00` is the price of that,
/// and it is a price paid in one trailing zero.
fn exact_value(value: f64) -> String {
    if !value.is_finite() {
        return "0.00".to_string();
    }
    let sign = if value < 0.0 { "-" } else { "" };
    let magnitude = value.abs();
    // Past the point where an f64 has hundredths to report, hand it to the
    // formatter rather than pretending the scaling below is exact.
    if magnitude >= 1e15 {
        return format!("{value:.2}");
    }
    let hundredths = (magnitude * 100.0).round() as u64;
    format!(
        "{sign}{}.{:02}",
        group_thousands(hundredths / 100),
        hundredths % 100
    )
}

/// Label for a point on the value axis.
///
/// Variant values are usually allele frequencies, which live between zero and
/// one and need decimals to say anything, but the track takes any quantity at
/// all, so large numbers still have to come out short.
fn format_value(value: f64) -> String {
    if !value.is_finite() {
        return "0".to_string();
    }
    if value >= 1_000_000.0 {
        return format!("{}M", trim(value / 1e6, 1));
    }
    if value >= 1_000.0 {
        return format!("{}k", trim(value / 1e3, 1));
    }
    if value >= 1.0 {
        return trim(value, 1);
    }
    trim(value, 2)
}

/// Rounds to `decimals` places and drops the point when nothing is left of it.
fn trim(value: f64, decimals: usize) -> String {
    let factor = 10f64.powi(decimals as i32);
    let rounded = (value * factor).round() / factor;
    if rounded == rounded.trunc() {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    #[test]
    fn categories_keep_their_first_appearance_order() {
        let track = VariantTrack::new(vec![
            Variant::new(1).category("indel"),
            Variant::new(2).category("snp"),
            Variant::new(3).category("indel"),
            Variant::new(4),
        ]);
        assert_eq!(track.categories(), vec!["indel", "snp"]);
    }

    #[test]
    fn colours_do_not_shift_when_a_variant_is_appended() {
        let base = vec![Variant::new(1).category("a"), Variant::new(2).category("b")];
        let mut extended = base.clone();
        extended.push(Variant::new(3).category("a"));
        assert_eq!(
            VariantTrack::new(base).categories(),
            VariantTrack::new(extended).categories()
        );
    }

    #[test]
    fn the_ceiling_ignores_missing_and_non_finite_values() {
        let track = VariantTrack::new(vec![
            Variant::new(1),
            Variant::new(2).value(f64::NAN),
            Variant::new(3).value(0.4),
            Variant::new(4).value(0.9),
        ]);
        assert_eq!(track.value_ceiling(), Some(0.9));
    }

    #[test]
    fn a_track_with_no_values_has_no_ceiling() {
        let track = VariantTrack::new(vec![Variant::new(1), Variant::new(2)]);
        assert_eq!(track.value_ceiling(), None);
    }

    #[test]
    fn a_scale_needs_values_to_measure() {
        let theme = Theme::light();
        let with_values = VariantTrack::new(vec![Variant::new(1).value(0.8)]);
        let without = VariantTrack::new(vec![Variant::new(1), Variant::new(2)]);
        assert!(with_values.y_axis_width(&theme) > 0.0);
        assert_eq!(
            without.y_axis_width(&theme),
            0.0,
            "full height stems have no scale to put a number on"
        );
    }

    #[test]
    fn ticks_have_no_scale_because_they_ignore_the_value() {
        let theme = Theme::light();
        let variants = vec![Variant::new(1).value(0.8)];
        let ticks = VariantTrack::new(variants.clone()).style(VariantStyle::Tick);
        let lollipops = VariantTrack::new(variants);
        assert_eq!(ticks.y_axis_width(&theme), 0.0);
        assert!(lollipops.y_axis_width(&theme) > 0.0);
    }

    #[test]
    fn the_scale_can_be_turned_off() {
        let theme = Theme::light();
        let track = VariantTrack::new(vec![Variant::new(1).value(0.8)]).show_scale(false);
        assert_eq!(track.y_axis_width(&theme), 0.0);
    }

    #[test]
    fn a_pinned_maximum_gives_a_scale_even_with_no_values() {
        let theme = Theme::light();
        let track = VariantTrack::new(vec![Variant::new(1)]).max(1.0);
        assert!(track.y_axis_width(&theme) > 0.0);
    }

    #[test]
    fn the_axis_labels_the_ceiling_and_the_floor() {
        let region = Region::parse("chr1:1-3000").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(
                VariantTrack::new(vec![
                    Variant::new(1000).value(0.64),
                    Variant::new(2000).value(0.21),
                ])
                .show_legend(false)
                .label("af"),
            )
            .to_svg();
        assert!(
            svg.contains(">0.64</text>"),
            "the tallest stem sets the top"
        );
        assert!(svg.contains(">0</text>"));
    }

    #[test]
    fn allele_frequencies_keep_their_decimals_and_counts_do_not() {
        assert_eq!(format_value(0.64), "0.64");
        assert_eq!(format_value(0.5), "0.5");
        assert_eq!(format_value(1.0), "1");
        assert_eq!(format_value(0.999), "1");
        assert_eq!(format_value(42.0), "42");
        assert_eq!(format_value(1500.0), "1.5k");
        assert_eq!(format_value(2_000_000.0), "2M");
        assert_eq!(format_value(f64::NAN), "0");
    }

    #[test]
    fn a_tooltip_prints_the_figure_and_the_axis_abbreviates_it() {
        // The axis has a few characters of room and the tooltip has none of
        // that constraint, so the two formatters are two functions. A read
        // depth of 1,234 rounded to `1.2k` had been rounded by a layout
        // decision taken somewhere else on the page.
        assert_eq!(format_value(1234.0), "1.2k");
        assert_eq!(exact_value(1234.0), "1,234.00");
        assert_eq!(exact_value(1_234_567.0), "1,234,567.00");
        // Two places always, so a column of allele fractions lines up.
        assert_eq!(exact_value(1.0), "1.00");
        assert_eq!(exact_value(0.98), "0.98");
        assert_eq!(exact_value(0.5), "0.50");
        assert_eq!(exact_value(0.0), "0.00");
        assert_eq!(exact_value(-0.25), "-0.25");
        assert_eq!(exact_value(f64::NAN), "0.00");
    }

    #[test]
    fn variants_outside_the_region_are_not_drawn() {
        let region = Region::parse("chr1:1-100").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(5_000)]).show_legend(false))
            .to_svg();
        assert!(!svg.contains("<circle"));
    }

    /// Every group a track opens has to be closed by exactly one `end_group`,
    /// or the rest of the figure ends up nested inside the last mark drawn.
    fn groups_balance(svg: &str) -> bool {
        svg.matches("<g ").count() + svg.matches("<g>").count() == svg.matches("</g>").count()
    }

    #[test]
    fn a_mark_says_where_it_is_what_it_is_and_what_it_carried() {
        let svg = Figure::new(Region::parse("chr1:761001-761300").unwrap())
            .show_region_label(false)
            .push(
                VariantTrack::new(vec![
                    Variant::new(761_154).value(1.00).category("missense"),
                    Variant::new(761_155).value(0.55).category("synonymous"),
                ])
                .show_legend(false),
            )
            .to_svg();
        // 1-based and thousands separated, the form the ruler prints, and the
        // number named rather than left bare.
        assert!(
            svg.contains("<title>variant, 761,155, missense, value 1.00</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>variant, 761,156, synonymous, value 0.55</title>"));
        assert!(groups_balance(&svg));
    }

    #[test]
    fn a_mark_with_only_a_position_claims_only_a_position() {
        let svg = Figure::new(Region::parse("chr1:761001-761300").unwrap())
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(761_200)]))
            .to_svg();
        assert!(svg.contains("<title>variant, 761,201</title>"), "{svg}");
    }

    #[test]
    fn ticks_are_not_named_and_lollipops_are() {
        // Choosing ticks is choosing not to resolve the marks one at a time,
        // and a mark a pointer cannot land on alone is not one worth naming.
        // Naming them anyway costs about two thirds of a genome-wide file.
        let variants = vec![Variant::new(761_154).value(0.55).category("missense")];
        let figure = |style| {
            Figure::new(Region::parse("chr1:761001-761300").unwrap())
                .show_region_label(false)
                .push(
                    VariantTrack::new(variants.clone())
                        .style(style)
                        .show_legend(false),
                )
                .to_svg()
        };

        let ticks = figure(VariantStyle::Tick);
        assert!(!ticks.contains("<title>"), "{ticks}");
        assert!(groups_balance(&ticks));

        let lollipops = figure(VariantStyle::Lollipop);
        assert!(
            lollipops.contains("<title>variant, 761,155, missense, value 0.55</title>"),
            "{lollipops}"
        );
        assert!(groups_balance(&lollipops));
    }

    #[test]
    fn a_mark_off_screen_is_not_named_either() {
        let svg = Figure::new(Region::parse("chr1:1-100").unwrap())
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(5_000)]).show_legend(false))
            .to_svg();
        assert!(!svg.contains("<title>"), "{svg}");
    }

    #[test]
    fn a_non_finite_value_is_left_out_rather_than_printed() {
        let svg = Figure::new(Region::parse("chr1:1-100").unwrap())
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(10).value(f64::NAN)]))
            .to_svg();
        assert!(svg.contains("<title>variant, 11</title>"), "{svg}");
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn ticks_draw_no_heads() {
        let region = Region::parse("chr1:1-100").unwrap();
        let variants = vec![Variant::new(10), Variant::new(20)];
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(VariantTrack::new(variants).style(VariantStyle::Tick))
            .to_svg();
        assert!(!svg.contains("<circle"));
        assert!(svg.contains("<line"));
    }
}
