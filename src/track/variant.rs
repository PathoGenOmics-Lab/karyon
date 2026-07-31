//! Point events along the sequence: SNPs, indels, insertion sites, peaks.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
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
            radius: 3.0,
            max: None,
            show_legend: true,
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

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let baseline = band.bottom();
        ctx.svg.line(
            band.x,
            baseline,
            band.right(),
            baseline,
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
                    ctx.svg.line(x, baseline, x, top, &color, 1.0);
                    ctx.svg.circle(x, top, self.radius, &color);
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
    fn variants_outside_the_region_are_not_drawn() {
        let region = Region::parse("chr1:1-100").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(5_000)]).show_legend(false))
            .to_svg();
        assert!(!svg.contains("<circle"));
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
