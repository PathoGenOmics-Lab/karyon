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
//! in front of another head still has an edge. Two heads closer than a head
//! is wide are one blot all the same, and are moved apart along the band, on
//! stems that rise at their calls and lean over to them; a head that touches
//! none stays where its call is.
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
//! one order, or one order of categories through
//! [`VariantTrack::category_order`], which fixes each category's slot whatever
//! the window holds. The command line hands it the categories of every call
//! of the file ranked from the most damaging consequence to the least, so a
//! category's colour does not hang on which call happens to come first in
//! the window, nor move in a zoom that leaves out a more damaging one, and
//! is one colour across the panels of a figure of several places.
//!
//! [`VariantTrack::uniform_color`] paints every call one colour instead, as a
//! figure whose calls are to match a colour already used elsewhere asks. The
//! categories keep their shapes, on the marks and in the key, so a missense
//! call and a synonymous one can still be told apart in one colour. There are
//! four shapes, and a fifth category takes the first one's again, so in one
//! colour no more than four can be told apart; the command line refuses
//! more.

use std::collections::BTreeSet;

use crate::scale::Scale;
use crate::style::{legible_ticks, Emphasis, QuantitativeAxis, Symbol};
use crate::svg::Anchor;
use crate::theme::Theme;
use crate::track::axis::group_thousands;
use crate::track::genotype::spelled;
use crate::track::{DrawContext, Legend, Track};

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
    /// The reference allele and the alternate allele this call is of, as a
    /// VCF's `REF` and `ALT` spell them, which its tooltip writes as `C>T`
    /// and, zoomed to a few pixels a base, the alternate beside its head.
    pub alleles: Option<(String, String)>,
}

impl Variant {
    /// A variant at a 0-based position.
    pub fn new(pos: u64) -> Self {
        Variant {
            pos,
            value: None,
            category: None,
            alleles: None,
        }
    }

    /// Sets the reference allele and the alternate one, as `REF` and `ALT`
    /// spell them.
    pub fn alleles(mut self, reference: impl Into<String>, alternate: impl Into<String>) -> Self {
        self.alleles = Some((reference.into(), alternate.into()));
        self
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
    uniform: Option<String>,
    axis: QuantitativeAxis,
    title: Option<String>,
    order: Vec<String>,
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
            uniform: None,
            axis: QuantitativeAxis::new(),
            title: None,
            order: Vec::new(),
        }
    }

    /// Deals the palette's slots to the categories `order` names, each the
    /// slot of its place in it, whether or not a variant of the track carries
    /// it. Categories it leaves out take the slots after it, in order of
    /// first appearance.
    ///
    /// Two tracks cut from one list at different windows then paint a
    /// category alike, which first appearance alone cannot promise: the
    /// window decides which category comes first. Name every category the
    /// whole list holds, and the zoom no longer moves a colour.
    pub fn category_order<S: Into<String>>(mut self, order: impl IntoIterator<Item = S>) -> Self {
        self.order.clear();
        for category in order {
            let category = category.into();
            if !self.order.contains(&category) {
                self.order.push(category);
            }
        }
        self
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
        self.axis.max = Some(max);
        self
    }

    /// Uses a shared quantitative-axis contract for range, ticks, units and
    /// reference lines.
    pub fn axis(mut self, axis: QuantitativeAxis) -> Self {
        self.max = axis.max;
        self.axis = axis;
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

    /// What the stems measure, written under the track's name while the
    /// value axis is drawn, such as `AF` for the allele fraction a VCF gives.
    pub fn axis_title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
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

    /// Paints every variant `color`, whatever its category, in place of a
    /// colour each category.
    ///
    /// Each category keeps the shape its place in the queue gives it, on the
    /// marks and in the key, which is then drawn in the one colour, so a
    /// missense call is still told from a synonymous one in a figure whose
    /// calls all had to be dark red. There are four shapes, and the fifth
    /// place in the queue takes the first one's again, so a track of more
    /// than four categories in one colour draws two of them alike. A tick
    /// has no shape, so ticks in one colour carry no category a reader could
    /// see, and are drawn with no key: one of shapes would name marks the
    /// band does not draw.
    pub fn uniform_color(mut self, color: impl Into<String>) -> Self {
        self.uniform = Some(color.into());
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
    /// The floor and ceiling the stems are measured against.
    ///
    /// A free ceiling is rounded to a value worth labelling, so a panel whose
    /// tallest call is 0.98 is read against 1 and not against 0.98.
    fn value_range(&self) -> (f64, f64) {
        let pinned = self.axis.max.or(self.max);
        let data_ceiling = pinned.or_else(|| self.value_ceiling()).unwrap_or(1.0);
        let (floor, ceiling) = self.axis.resolve(0.0, data_ceiling);
        if pinned.is_some() {
            (floor, ceiling)
        } else {
            self.axis.nice(floor, ceiling)
        }
    }

    fn value_ceiling(&self) -> Option<f64> {
        self.variants
            .iter()
            .filter_map(|v| v.value)
            .filter(|v| v.is_finite() && *v > 0.0)
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            })
    }

    /// The categories the track holds, each with its palette slot, in the
    /// order of their slots: the place [`VariantTrack::category_order`]
    /// gave it, else after every place it gave, by first appearance.
    fn slots(&self) -> Vec<(&str, usize)> {
        let mut unlisted = self.order.len();
        let mut slots: Vec<(&str, usize)> = self
            .categories()
            .into_iter()
            .map(|category| {
                let slot = self.order.iter().position(|known| known == category);
                let slot = slot.unwrap_or_else(|| {
                    unlisted += 1;
                    unlisted - 1
                });
                (category, slot)
            })
            .collect();
        slots.sort_by_key(|&(_, slot)| slot);
        slots
    }

    fn legend(&self, theme: &Theme) -> Legend {
        self.slots()
            .into_iter()
            .fold(Legend::new(), |legend, (category, slot)| {
                legend.symbol(category, self.ink(theme, slot), theme.symbol(slot))
            })
    }

    /// Whether the key is drawn: where it is asked for, and not over ticks in
    /// one colour, which nothing tells apart.
    fn keyed(&self) -> bool {
        self.show_legend && !(self.uniform.is_some() && self.style == VariantStyle::Tick)
    }

    /// The colour of the category in palette slot `slot`: its own, or the
    /// one every call is painted where [`VariantTrack::uniform_color`] says.
    fn ink<'a>(&'a self, theme: &'a Theme, slot: usize) -> &'a str {
        self.uniform.as_deref().unwrap_or_else(|| theme.color(slot))
    }
}

impl Track for VariantTrack {
    fn noun(&self) -> &str {
        "variant calls"
    }

    fn height(&self, scale: &Scale) -> f64 {
        if self.keyed() {
            let theme = Theme::default();
            self.height
                .max(self.legend(&theme).height(scale.width(), &theme) + self.radius + 2.0)
        } else {
            self.height
        }
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn axis_title(&self) -> Option<&str> {
        self.has_scale().then_some(self.title.as_deref()).flatten()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.has_scale() {
            return 0.0;
        }
        let (floor, ceiling) = self.value_range();
        self.axis.label_room(floor, ceiling, theme.font_size - 1.0) + 8.0
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

        let slots = self.slots();
        let default_color = self
            .uniform
            .clone()
            .or_else(|| self.color.clone())
            .unwrap_or_else(|| ctx.theme.accent.clone());
        // A colour and a symbol are a property of the category, not of the
        // variant, and both were being worked out inside the loop: a linear
        // scan of the category list and a fresh `String` for every mark. On a
        // panel of two hundred thousand calls with a handful of consequences
        // that is two hundred thousand scans and two hundred thousand
        // allocations for a handful of distinct answers.
        let palette: Vec<(String, Symbol)> = slots
            .iter()
            .map(|&(_, slot)| {
                (
                    self.ink(ctx.theme, slot).to_string(),
                    ctx.theme.symbol(slot),
                )
            })
            .collect();
        let slot_of = |category: Option<&str>| -> Option<usize> {
            category.and_then(|name| slots.iter().position(|&(c, _)| c == name))
        };
        let (floor, ceiling) = self.value_range();

        let legend = self.legend(ctx.theme);

        // Leave room for the legend so a tall stem does not run through it.
        let legend_room = if self.keyed() && !slots.is_empty() {
            legend.height(band.w, ctx.theme)
        } else {
            0.0
        };
        let radius = self.radius * ctx.visual_scale;
        let stem_room = (band.h - legend_room - radius).max(ctx.px(2.0));
        let y_of = |value: f64| {
            baseline - ((value - floor) / (ceiling - floor)).clamp(0.0, 1.0) * stem_room
        };

        let size = ctx.theme.font_size - 1.0;
        let ticks = self.axis.values(floor, ceiling);
        let shown = legible_ticks(&ticks, y_of, size);
        if self.has_scale() {
            // The height a variant at the ceiling reaches, which is where the
            // top of the scale is and nowhere else.
            for &value in &shown {
                let y = y_of(value);
                ctx.svg.line(
                    band.x,
                    y,
                    band.right(),
                    y,
                    &ctx.theme.rule,
                    ctx.theme.tokens.hairline,
                );
            }
            for reference in &self.axis.references {
                if !reference.value.is_finite()
                    || reference.value < floor
                    || reference.value > ceiling
                {
                    continue;
                }
                let y = y_of(reference.value);
                let style = ctx.theme.mark_style(reference.emphasis);
                let ink = if reference.emphasis == Emphasis::Alert {
                    ctx.theme.color(1)
                } else {
                    &ctx.theme.muted
                };
                ctx.svg.line_pattern(
                    band.x,
                    y,
                    band.right(),
                    y,
                    ink,
                    style.stroke_width,
                    reference.pattern,
                );
            }
            if ctx.axis.w > 0.0 {
                let right = ctx.axis.right() - 4.0;
                for (value, label) in ticks.iter().zip(self.axis.labels(&ticks)) {
                    if !shown.contains(value) {
                        continue;
                    }
                    let y = y_of(*value);
                    ctx.svg.text(
                        right,
                        (y + size * 0.35).min(baseline),
                        &label,
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

        // Ticks are the caller saying the marks are too dense to be told apart,
        // and this is where that statement is made good. Two ticks on the same
        // pixel column in the same colour are the same ink: a tick carries no
        // height, no symbol and no tooltip, so the second one adds a line to
        // the document and nothing to the picture. Measured on two hundred
        // thousand calls over two megabases, they landed on eight hundred and
        // sixty-seven columns, so two hundred and thirty of every two hundred
        // and thirty-one were drawing over one already there.
        //
        // Lollipops are not folded, and must not be: their height is the value
        // and their head is a thing a pointer lands on, so two at one column
        // are two different claims.
        let mut drawn: BTreeSet<(i64, usize)> = BTreeSet::new();

        // Where each lollipop's head is drawn: at its call, but for heads that
        // would cover each other, which are moved apart along the band.
        let heads: Vec<(f64, f64)> = self
            .variants
            .iter()
            .filter(|variant| ctx.region.contains(variant.pos))
            .map(|variant| {
                let top = variant
                    .value
                    .filter(|value| value.is_finite())
                    .map_or(baseline - stem_room, y_of);
                (ctx.scale.x_center(variant.pos), top)
            })
            .collect();
        let moved = if self.style == VariantStyle::Lollipop {
            apart(&heads, radius, band.x, band.right())
        } else {
            Vec::new()
        };
        // Zoomed until the reference's bases are letters, a head with room
        // on its right has its alternate allele written there.
        let lettered = ctx.scale.px_per_bp() >= crate::track::sequence::LETTER_PX;
        let letter_size = ctx.theme.font_size - 1.0;
        let mut at = 0;

        for variant in &self.variants {
            if !ctx.region.contains(variant.pos) {
                continue;
            }
            let x = ctx.scale.x_center(variant.pos);
            let slot = slot_of(variant.category.as_deref());
            let color: &str = match slot {
                Some(index) => palette
                    .get(index)
                    .map_or(default_color.as_str(), |(color, _)| color.as_str()),
                None => default_color.as_str(),
            };

            match self.style {
                VariantStyle::Tick => {
                    // The column, not the coordinate: the mark is a hairline and
                    // its own width is what decides whether two of them are one.
                    if !drawn.insert((x.round() as i64, slot.unwrap_or(usize::MAX))) {
                        continue;
                    }
                    // Not named, and the style is the reason. Choosing ticks is
                    // the caller saying the marks are too dense to be told
                    // apart, which is the same statement a binned track makes
                    // when it draws one point per pixel column: a mark a
                    // pointer cannot land on alone is not one worth naming.
                    // Naming them anyway costs about two thirds of the file.
                    ctx.svg.line(
                        x,
                        baseline,
                        x,
                        baseline - stem_room,
                        color,
                        ctx.theme.tokens.stroke,
                    );
                }
                VariantStyle::Lollipop => {
                    let top = variant
                        .value
                        .filter(|value| value.is_finite())
                        .map_or(baseline - stem_room, y_of);
                    let head = moved.get(at).copied().unwrap_or(x);
                    // A lollipop is a stem and a head, and the group is what
                    // makes the two of them one thing. A tooltip on half of a
                    // mark is worse than none.
                    ctx.svg
                        .begin_titled(&tooltip(variant, self.title.as_deref()));
                    if (head - x).abs() < 1e-9 {
                        ctx.svg
                            .line(x, baseline, x, top, color, ctx.theme.tokens.stroke);
                    } else {
                        // A head moved off its call stands on a stem that
                        // rises at the call and leans over to it at the
                        // last, so the stem still says where the call is.
                        let knee = (top + 2.0 * radius).min(baseline);
                        ctx.svg.polyline(
                            &[(x, baseline), (x, knee), (head, top)],
                            color,
                            ctx.theme.tokens.stroke,
                        );
                    }
                    // The ring is what keeps two variants a base apart reading
                    // as two variants instead of one blob.
                    let symbol = slot
                        .and_then(|index| palette.get(index))
                        .map_or(Symbol::Circle, |(_, symbol)| *symbol);
                    ctx.svg.symbol_ringed(
                        head,
                        top,
                        radius,
                        symbol,
                        color,
                        ctx.theme.surface(),
                        ctx.theme.tokens.hairline,
                    );
                    let written = variant
                        .alleles
                        .as_ref()
                        .map(|(_, alternate)| alternate.as_str())
                        .filter(|alternate| lettered && written_beside(alternate));
                    if let Some(alternate) = written {
                        let width = crate::svg::text_width(alternate, letter_size);
                        let beside = Beside {
                            radius,
                            size: letter_size,
                            width,
                            right: band.right(),
                        };
                        if beside.has_room(at, &heads, &moved) {
                            // A base in the colour the reference's letter of
                            // it is drawn in, and more than one in the ink.
                            let ink = match alternate.as_bytes() {
                                [base] => ctx.theme.bases.of(*base),
                                _ => &ctx.theme.foreground,
                            };
                            ctx.svg.text(
                                head + radius + 1.5,
                                top + letter_size * 0.35,
                                alternate,
                                ink,
                                letter_size,
                                Anchor::Start,
                            );
                        }
                    }
                    ctx.svg.end_group();
                    at += 1;
                }
            }
        }

        if self.keyed() && !slots.is_empty() {
            legend.draw(ctx.svg, band.x, band.y, band.w, ctx.theme);
        }
    }
}

/// Whether an alternate allele is written beside its head: a few bases, three
/// at most, of the letters a sequence is spelled in. A longer one, or a
/// symbolic allele such as `<DEL>`, is named in the tooltip alone.
fn written_beside(alternate: &str) -> bool {
    (1..=3).contains(&alternate.len()) && alternate.bytes().all(|b| b.is_ascii_alphabetic())
}

/// The room a head's alternate allele needs on its right.
struct Beside {
    radius: f64,
    size: f64,
    width: f64,
    right: f64,
}

impl Beside {
    /// Whether the allele beside head `at` of `heads`, their calls' centres,
    /// drawn at `drawn`, stays in the band and clear of every other head and
    /// stem: a mark whose head is lower than the letters passes under them,
    /// and any other in their way is crossed by them.
    fn has_room(&self, at: usize, heads: &[(f64, f64)], drawn: &[f64]) -> bool {
        let (x, top) = (drawn[at], heads[at].1);
        let (from, to) = (x + self.radius + 1.0, x + self.radius + 2.0 + self.width);
        if to > self.right {
            return false;
        }
        let bottom = top + self.size * 0.4;
        !heads
            .iter()
            .zip(drawn)
            .enumerate()
            .any(|(other, ((stem, high), head))| {
                let across = (head + self.radius > from && head - self.radius < to)
                    || (*stem >= from - 1.0 && *stem <= to + 1.0);
                other != at && across && high - self.radius < bottom
            })
    }
}

/// Where each head of `heads`, their centres in the order they are drawn,
/// is drawn across a band from `left` to `right`, so no two of them cover
/// each other.
///
/// Two heads whose centres are closer than a head is wide were one mark: a
/// call at 761,110 of 0.50 and one four bases on of 0.53 drew two heads of
/// four pixels a pixel apart, and the figure said nothing. Heads that touch,
/// and those that touch them, are laid in a row a head and a pixel apart,
/// in the order of their calls, centred where their calls are on average,
/// and kept inside the band; a head that touches none stays where it was, so
/// a figure with no heads that touch is drawn as it always was. A row moved
/// onto another head takes that one in and is laid out again. Heads too many
/// for the band to hold in a row are left where they are, which is the
/// density ticks are for, and so is a band so crowded that a head has
/// thousands of others within a head's width of it: two hundred thousand
/// calls over four megabases were twenty times as slow to draw, looking at
/// every pair of heads, for rows no band could hold.
fn apart(heads: &[(f64, f64)], radius: f64, left: f64, right: f64) -> Vec<f64> {
    let mut drawn: Vec<f64> = heads.iter().map(|(x, _)| *x).collect();
    let width = 2.0 * radius;
    // The pairs of heads within a head's width of each other across the
    // band, which the rounds below look at each of.
    const MOST_PAIRS: usize = 1_000_000;
    let mut across = drawn.clone();
    across.sort_by(f64::total_cmp);
    let mut pairs = 0usize;
    let mut reach = 0;
    for (at, x) in across.iter().enumerate() {
        reach = reach.max(at);
        while reach + 1 < across.len() && across[reach + 1] - x < width {
            reach += 1;
        }
        pairs += reach - at;
        if pairs > MOST_PAIRS {
            return drawn;
        }
    }
    let mut group: Vec<usize> = (0..heads.len()).collect();
    fn root(group: &mut [usize], mut at: usize) -> usize {
        while group[at] != at {
            group[at] = group[group[at]];
            at = group[at];
        }
        at
    }
    let pitch = width + 1.0;
    // A row moved onto a head takes it in, so the groups only grow, and a
    // few rounds settle every figure a lollipop is drawn for.
    for _ in 0..16 {
        let mut order: Vec<usize> = (0..heads.len()).collect();
        order.sort_by(|a, b| drawn[*a].total_cmp(&drawn[*b]));
        let mut joined = false;
        for (from, &a) in order.iter().enumerate() {
            for &b in &order[from + 1..] {
                let dx = drawn[b] - drawn[a];
                if dx >= width {
                    break;
                }
                let dy = heads[b].1 - heads[a].1;
                if dx.hypot(dy) < width {
                    let (ra, rb) = (root(&mut group, a), root(&mut group, b));
                    if ra != rb {
                        group[ra] = rb;
                        joined = true;
                    }
                }
            }
        }
        if !joined {
            break;
        }
        let mut members: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();
        for at in 0..heads.len() {
            let top = root(&mut group, at);
            members.entry(top).or_default().push(at);
        }
        for mut row in members.into_values().filter(|row| row.len() > 1) {
            let span = (row.len() - 1) as f64 * pitch;
            if span > right - left - width {
                continue;
            }
            row.sort_by(|a, b| {
                heads[*a]
                    .0
                    .total_cmp(&heads[*b].0)
                    .then(heads[*a].1.total_cmp(&heads[*b].1))
                    .then(a.cmp(b))
            });
            let middle = row.iter().map(|at| heads[*at].0).sum::<f64>() / row.len() as f64;
            let first = (middle - span / 2.0)
                .max(left + radius)
                .min(right - radius - span);
            for (index, at) in row.iter().enumerate() {
                drawn[*at] = first + index as f64 * pitch;
            }
        }
    }
    drawn
}

/// What a reader hovering one mark is told.
///
/// The noun first, then the position in the 1-based form the ruler under it
/// prints and the change, `C>T`, where the call carries its alleles, then
/// what kind of event it is and what number it carried. All but the position
/// are optional in the data and none is invented when it is absent: a mark
/// with nothing but a position says nothing but its noun and a position.
///
/// The noun is there because every tooltip in the crate is `what it is, where
/// it is`, and a mark that opened on a bare `31,218,401` was the one glyph on
/// a sheet answering in a different grammar from the pileup above it.
///
/// The number is named rather than printed bare, because `0.55` on its own is
/// not a statement. What it is named is what the caller titled the axis, `AF`
/// for the allele fractions a VCF gives, and `value` where the caller said
/// nothing, since the track takes any quantity at all: the same band draws
/// allele fractions, peak heights and read counts, and only the caller knows
/// which.
///
/// Only lollipops are named. See the module documentation for why a tick is
/// not.
fn tooltip(variant: &Variant, measure: Option<&str>) -> String {
    let mut text = format!(
        "variant, {}",
        group_thousands(variant.pos.saturating_add(1))
    );
    // The change, `C>T`, after the place it is at, each allele spelled to
    // twelve characters and counted past them, as a genotype's tooltip
    // spells it: an indel's REF can run to thousands of bases.
    if let Some((reference, alternate)) = &variant.alleles {
        text.push(' ');
        text.push_str(&spelled(reference));
        text.push('>');
        text.push_str(&spelled(alternate));
    }
    if let Some(category) = variant.category.as_deref() {
        if !category.is_empty() {
            text.push_str(", ");
            text.push_str(category);
        }
    }
    // Named as the axis names it, `AF 0.48`, where the axis has a title.
    if let Some(value) = variant.value {
        if value.is_finite() {
            text.push_str(", ");
            text.push_str(measure.unwrap_or("value"));
            text.push(' ');
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
        // The one caller filters these out. Saying `0.00` would claim a
        // measurement of zero, so if this is ever reached it says what it got.
        return format!("{value}");
    }
    let sign = if value < 0.0 { "-" } else { "" };
    let magnitude = value.abs();
    // An f64 holds consecutive integers up to 2^53, so it holds consecutive
    // hundredths up to a hundredth of that, and above it the two decimal
    // places below are the formatter's invention rather than the value's:
    // 999,999,999,999,999 comes back as `.04`. Past that the magnitude is all
    // there is to say, and it says it in seven characters rather than in the
    // three hundred and ten `f64::MAX` expands to under a pointer.
    const EXACT_HUNDREDTHS: f64 = (1u64 << 53) as f64 / 100.0;
    if magnitude >= EXACT_HUNDREDTHS {
        return format!("{value:.3e}");
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
#[cfg(test)]
fn format_value(value: f64) -> String {
    if !value.is_finite() {
        return "0".to_string();
    }
    if value >= 1_000_000.0 {
        return format!("{}M", crate::svg::text_rounded(value / 1e6, 1));
    }
    if value >= 1_000.0 {
        return format!("{}k", crate::svg::text_rounded(value / 1e3, 1));
    }
    if value >= 1.0 {
        return crate::svg::text_rounded(value, 1);
    }
    crate::svg::text_rounded(value, 2)
}

/// Rounds to `decimals` places and drops the point when nothing is left of it.
#[cfg(test)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    /// The title names what the stems measure, and there is nothing to
    /// title where no stem measures anything.
    #[test]
    fn the_axis_is_titled_only_while_it_is_drawn() {
        let valued = vec![Variant::new(10).value(0.4), Variant::new(20).value(0.9)];
        let titled = VariantTrack::new(valued.clone()).axis_title("AF");
        assert_eq!(Track::axis_title(&titled), Some("AF"));
        let svg = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .push(titled.clone().label("calls"))
            .to_svg();
        assert!(svg.contains(">AF</text>"), "{svg}");
        // And a stem's tooltip names its value the same way.
        assert!(svg.contains("<title>variant, 11, AF 0.40</title>"), "{svg}");
        let plain = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .push(VariantTrack::new(valued.clone()))
            .to_svg();
        assert!(
            plain.contains("<title>variant, 11, value 0.40</title>"),
            "{plain}"
        );
        assert_eq!(Track::axis_title(&titled.clone().show_scale(false)), None);
        assert_eq!(Track::axis_title(&titled.style(VariantStyle::Tick)), None);
        let unvalued = VariantTrack::new(vec![Variant::new(10)]).axis_title("AF");
        assert_eq!(
            Track::axis_title(&unvalued),
            None,
            "no value, no scale, no title"
        );
    }

    /// One colour for every call takes the palette off the marks and the key,
    /// and leaves each category its shape in both.
    #[test]
    fn a_uniform_colour_paints_every_call_and_keeps_each_category_s_shape() {
        let theme = Theme::default();
        let calls = vec![
            Variant::new(10).value(0.5).category("missense_variant"),
            Variant::new(40).value(0.9).category("synonymous_variant"),
            Variant::new(70).value(0.3),
        ];
        let svg = |track: VariantTrack| {
            Figure::new(Region::new("chr1", 0, 100).unwrap())
                .push(track)
                .to_svg()
        };
        let painted = svg(VariantTrack::new(calls.clone()).uniform_color("#8b0000"));
        for slot in 0..2 {
            assert!(
                !painted.contains(theme.color(slot)),
                "slot {slot} is still painted: {painted}"
            );
        }
        assert!(!painted.contains(&theme.accent), "{painted}");
        // Three stems, three heads and two keys, all of them the one colour.
        assert!(painted.matches("#8b0000").count() >= 8, "{painted}");
        assert!(painted.contains(">missense_variant</text>"), "{painted}");
        assert!(painted.contains(">synonymous_variant</text>"), "{painted}");
        // The shapes are the ones the categories take in colour: the second
        // category is a square in both figures.
        let coloured = svg(VariantTrack::new(calls));
        let squares = |svg: &str| svg.matches("<rect").count();
        assert_eq!(squares(&painted), squares(&coloured), "{painted}");
        assert!(coloured.contains(theme.color(1)), "{coloured}");
        // Ticks in one colour have no shape to key, and no key.
        let ticks = svg(VariantTrack::new(vec![
            Variant::new(10).category("missense_variant"),
            Variant::new(40).category("synonymous_variant"),
        ])
        .style(VariantStyle::Tick)
        .uniform_color("#8b0000"));
        assert!(!ticks.contains("missense_variant"), "{ticks}");
    }

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
    fn a_category_order_keeps_a_slot_whatever_the_window_holds() {
        let theme = Theme::default();
        // The window of a zoom: only the second category of the file is in
        // it, and it keeps the second colour; first appearance gave it the
        // first.
        let zoom = VariantTrack::new(vec![Variant::new(5).category("missense")]);
        let ordered = zoom.clone().category_order(["synonymous", "missense"]);
        let svg = |track: VariantTrack| {
            Figure::new(Region::new("chr1", 0, 10).unwrap())
                .push(track)
                .to_svg()
        };
        assert!(svg(zoom).contains(&format!("missense, colour {}", theme.color(0))));
        let drawn = svg(ordered);
        assert!(
            drawn.contains(&format!("missense, colour {}", theme.color(1))),
            "{drawn}"
        );
        assert!(
            !drawn.contains("synonymous"),
            "only what is drawn is keyed: {drawn}"
        );

        // A category the order leaves out comes after it, and the key runs
        // in slot order.
        let track = VariantTrack::new(vec![
            Variant::new(1).category("odd"),
            Variant::new(2).category("b"),
        ])
        .category_order(["a", "b"]);
        let slots = track.slots();
        assert_eq!(slots, [("b", 1), ("odd", 2)]);
    }

    #[test]
    fn a_categorical_legend_wraps_and_copies_each_point_shape() {
        let variants: Vec<Variant> = (0..12)
            .map(|index| {
                Variant::new(index)
                    .value(1.0)
                    .category(format!("class_{index}"))
            })
            .collect();
        let track = VariantTrack::new(variants);
        let wide = Scale::new(&Region::new("chr1", 0, 20).unwrap(), 0.0, 900.0);
        let narrow = Scale::new(&Region::new("chr1", 0, 20).unwrap(), 0.0, 150.0);
        assert!(Track::height(&track, &narrow) > Track::height(&track, &wide));
        let svg = crate::Figure::new(Region::new("chr1", 0, 20).unwrap())
            .width(300.0)
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("<circle"), "{svg}");
        assert!(svg.contains("<polygon"), "{svg}");
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
    fn the_axis_labels_a_round_ceiling_and_the_floor() {
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
            svg.contains(">0.75</text>"),
            "the tallest stem is rounded up to a value worth labelling"
        );
        assert!(!svg.contains(">0.64</text>"));
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
        // A number too large to have hundredths is given as a magnitude, not
        // as the three hundred and ten digits the formatter would expand it to.
        assert_eq!(exact_value(f64::MAX), "1.798e308");
        assert!(exact_value(f64::MAX).len() < 20);
        assert_eq!(exact_value(-f64::MAX), "-1.798e308");
        // The last magnitude whose hundredths an f64 still holds, and the
        // first one whose it does not.
        assert_eq!(exact_value(9e13), "90,000,000,000,000.00");
        assert_eq!(exact_value(1e14), "1.000e14");
        // Two places always, so a column of allele fractions lines up.
        assert_eq!(exact_value(1.0), "1.00");
        assert_eq!(exact_value(0.98), "0.98");
        assert_eq!(exact_value(0.5), "0.50");
        assert_eq!(exact_value(0.0), "0.00");
        assert_eq!(exact_value(-0.25), "-0.25");
        // The caller filters a non-finite value out, so this is unreachable
        // from a figure. If it is ever reached it says what it got, since
        // "0.00" would claim a measurement of zero that nobody made.
        assert_eq!(exact_value(f64::NAN), "NaN");
        assert_eq!(exact_value(f64::INFINITY), "inf");
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

    /// The heads of a document's lollipops, each with the title of its mark
    /// and where its stem stands at the foot of the band.
    fn heads(svg: &str) -> Vec<(String, f64, f64, f64)> {
        svg.split("<g><title>variant, ")
            .skip(1)
            .map(|mark| {
                let title = mark[..mark.find("</title>").unwrap()].to_string();
                let number = |after: &str, key: &str| -> f64 {
                    let rest = &after[after.find(&format!("{key}=\"")).unwrap() + key.len() + 2..];
                    rest[..rest.find('"').unwrap()].parse().unwrap()
                };
                let stem = match mark.find("points=\"") {
                    Some(at) if at < mark.find("<circle").unwrap() => {
                        let points = &mark[at + 8..];
                        points[..points.find(',').unwrap()].parse().unwrap()
                    }
                    _ => number(mark, "x1"),
                };
                let circle = &mark[mark.find("<circle").unwrap()..];
                (title, stem, number(circle, "cx"), number(circle, "cy"))
            })
            .collect()
    }

    /// Two heads closer than a head is wide are moved apart along the band,
    /// each on a stem that rises at its call and leans over to it, and each
    /// named at its call: a call at 761,110 of 0.50 and one four bases on of
    /// 0.53 were two heads a pixel apart, one mark. A head that touches
    /// neither stays where it is, on a straight stem.
    #[test]
    fn heads_that_would_cover_each_other_are_moved_apart() {
        let svg = Figure::new(Region::parse("NC_1:759,456-763,676").unwrap())
            .show_region_label(false)
            .push(
                VariantTrack::new(vec![
                    Variant::new(761_109).value(0.50).category("missense"),
                    Variant::new(761_113).value(0.53).category("missense"),
                    Variant::new(762_000).value(0.52).category("missense"),
                ])
                .show_legend(false),
            )
            .to_svg();
        let drawn = heads(&svg);
        assert_eq!(drawn.len(), 3, "{svg}");
        let (first, second, alone) = (&drawn[0], &drawn[1], &drawn[2]);
        assert!(first.0.starts_with("761,110,"), "{first:?}");
        assert!(second.0.starts_with("761,114,"), "{second:?}");
        // Apart by a head and a pixel, in the order of the calls, about
        // where the two calls are.
        // Written to a thousandth of a pixel.
        assert!((second.2 - first.2 - 9.0).abs() < 2e-3, "{drawn:?}");
        assert!(first.2 < first.1 && second.2 > second.1, "{drawn:?}");
        assert!(((first.2 + second.2) / 2.0 - (first.1 + second.1) / 2.0).abs() < 2e-3);
        assert!(svg.matches("<polyline").count() == 2, "{svg}");
        // The third is where its call is.
        assert_eq!(alone.1, alone.2, "{alone:?}");
        assert!(groups_balance(&svg));
    }

    /// Heads that touch none stay exactly where they were, on a straight
    /// stem: a figure with no heads that collide is drawn as it always was.
    #[test]
    fn heads_that_touch_none_are_where_their_calls_are() {
        let calls = vec![
            Variant::new(761_109).value(0.10),
            Variant::new(761_113).value(0.90),
            Variant::new(761_600).value(0.50),
        ];
        let svg = Figure::new(Region::parse("NC_1:759,456-763,676").unwrap())
            .show_region_label(false)
            .push(VariantTrack::new(calls).show_legend(false))
            .to_svg();
        assert!(!svg.contains("<polyline"), "{svg}");
        for (title, stem, x, _) in heads(&svg) {
            assert_eq!(stem, x, "{title}");
        }
    }

    /// Heads laid in a row take in a head the row reaches, and are laid out
    /// again with it; a row kept inside the band; and heads too many for the
    /// band to hold in a row left where they are.
    #[test]
    fn a_row_of_heads_takes_in_what_it_reaches_and_stays_in_the_band() {
        // Three at one height, the third far enough from the first two to
        // touch only the row they are moved into.
        let moved = apart(
            &[(100.0, 50.0), (101.0, 50.0), (110.0, 50.0)],
            4.0,
            0.0,
            400.0,
        );
        let middle = (100.0 + 101.0 + 110.0) / 3.0;
        for (at, expected) in moved.iter().zip([middle - 9.0, middle, middle + 9.0]) {
            assert!((at - expected).abs() < 1e-9, "{moved:?}");
        }
        // At the left edge the row starts a head's half inside it.
        let edge = apart(&[(1.0, 50.0), (2.0, 50.0)], 4.0, 0.0, 400.0);
        assert_eq!(edge, [4.0, 13.0]);
        // Far apart in height, they do not touch.
        let tall = apart(&[(100.0, 10.0), (100.5, 90.0)], 4.0, 0.0, 400.0);
        assert_eq!(tall, [100.0, 100.5]);
        // Fifty in a band of a hundred pixels cannot be a row.
        let crowd: Vec<(f64, f64)> = (0..50).map(|i| (50.0 + i as f64 * 0.1, 50.0)).collect();
        let left = apart(&crowd, 4.0, 0.0, 100.0);
        assert_eq!(left, crowd.iter().map(|(x, _)| *x).collect::<Vec<_>>());
        // Nor can twenty thousand over half a band of a thousand, and a band
        // that crowded is left as it is without looking at each of its pairs,
        // the two touching heads far from the crowd with it.
        let mut throng: Vec<(f64, f64)> = (0..20_000)
            .map(|i| (i as f64 * 0.025, 50.0 + (i % 7) as f64))
            .collect();
        throng.extend([(900.0, 50.0), (901.0, 50.0)]);
        let left = apart(&throng, 4.0, 0.0, 1_000.0);
        assert_eq!(left, throng.iter().map(|(x, _)| *x).collect::<Vec<_>>());
    }

    /// A mark says the change it is, after where it is: `761,139 C>T`. An
    /// allele longer than twelve bases is its first twelve and how many more,
    /// as a genotype's tooltip spells it.
    #[test]
    fn a_mark_says_its_change_after_where_it_is() {
        let long = "ACGTACGTACGTACGTACGT";
        let svg = Figure::new(Region::parse("chr1:761,001-761,300").unwrap())
            .show_region_label(false)
            .push(
                VariantTrack::new(vec![
                    Variant::new(761_138)
                        .value(0.62)
                        .category("missense_variant")
                        .alleles("C", "T"),
                    Variant::new(761_200).alleles(long, "A"),
                ])
                .axis_title("AF")
                .show_legend(false),
            )
            .to_svg();
        assert!(
            svg.contains("<title>variant, 761,139 C&gt;T, missense_variant, AF 0.62</title>"),
            "{svg}"
        );
        assert!(
            svg.contains("<title>variant, 761,201 ACGTACGTACGT+8&gt;A</title>"),
            "{svg}"
        );
    }

    /// Zoomed until a base is a letter wide, a head with room on its right
    /// has its alternate allele written there, a base in the colour of its
    /// letter; not over a wider window, not where another mark stands in the
    /// way, and not for an allele longer than three bases or one with no
    /// bases to it.
    #[test]
    fn at_a_letter_a_base_a_head_with_room_writes_its_alternate() {
        let theme = Theme::light();
        let calls = vec![
            Variant::new(761_109).value(0.35).alleles("A", "T"),
            Variant::new(761_138).value(0.62).alleles("C", "G"),
            // A taller call a base on stands where the letter would go.
            Variant::new(761_139).value(0.90).alleles("C", "A"),
            Variant::new(761_150).value(0.20).alleles("C", "<DEL>"),
            Variant::new(761_160).value(0.20).alleles("C", "CAAAA"),
        ];
        let drawn = |window: &str| {
            Figure::new(Region::parse(window).unwrap())
                .show_region_label(false)
                .push(VariantTrack::new(calls.clone()).show_legend(false))
                .to_svg()
        };
        let letters = |svg: &str| -> Vec<(String, String)> {
            svg.split("<text ")
                .skip(1)
                .filter_map(|text| {
                    let fill = text.split("fill=\"").nth(1)?.split('"').next()?;
                    let body = text.split_once('>')?.1.split("</text>").next()?;
                    Some((body.to_string(), fill.to_string()))
                })
                .filter(|(body, _)| {
                    body.len() <= 5 && body.chars().all(|c| c.is_ascii_alphabetic())
                })
                .collect()
        };
        let near = drawn("chr1:761,101-761,170");
        assert_eq!(
            letters(&near),
            [
                ("T".to_string(), theme.bases.of(b'T').to_string()),
                ("A".to_string(), theme.bases.of(b'A').to_string()),
            ],
            "{near}"
        );
        let far = drawn("chr1:760,001-763,000");
        assert!(letters(&far).is_empty(), "{far}");
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
    fn ticks_fold_onto_their_pixel_column_and_lollipops_do_not() {
        // A tick carries no height, no symbol and no tooltip, so two of them on
        // one pixel column in one colour are the same ink and the second draws
        // nothing. Measured before this: two hundred thousand calls over two
        // megabases put 230.7 ticks on every column, and the file was sixteen
        // megabytes for eight hundred and sixty-seven columns of picture.
        //
        // A lollipop is not the same: its height is its value and its head is a
        // thing a pointer lands on, so two at one column are two claims and
        // both are drawn.
        let dense: Vec<Variant> = (0..2_000u64)
            .map(|i| Variant::new(1_000 + i).value((i % 10) as f64 / 10.0))
            .collect();
        let figure = |style| {
            Figure::new(Region::parse("chr1:1-100000").unwrap())
                .show_region_label(false)
                .push(
                    VariantTrack::new(dense.clone())
                        .style(style)
                        .show_legend(false),
                )
                .to_svg()
        };

        let ticks = figure(VariantStyle::Tick);
        let drawn = ticks.matches("stroke-width=\"1.2\"").count();
        assert!(
            drawn < 40,
            "two thousand ticks over a few pixels drew {drawn} marks"
        );
        assert!(groups_balance(&ticks));

        let lollipops = figure(VariantStyle::Lollipop);
        assert!(
            lollipops.matches("<circle").count() >= 2_000,
            "lollipops must not fold: {} heads",
            lollipops.matches("<circle").count()
        );
    }

    #[test]
    fn two_categories_on_one_column_are_two_ticks() {
        // The fold is per colour, because two ticks of different colours at one
        // column are two different statements and the picture has to keep both.
        let svg = Figure::new(Region::parse("chr1:1-100000").unwrap())
            .show_region_label(false)
            .push(
                VariantTrack::new(vec![
                    Variant::new(1_000).category("missense"),
                    Variant::new(1_001).category("synonymous"),
                    Variant::new(1_002).category("missense"),
                ])
                .style(VariantStyle::Tick)
                .show_legend(false),
            )
            .to_svg();
        // The baseline is a line too, so the marks are counted by their own
        // stroke rather than by the element name.
        let ticks: Vec<&str> = svg.matches("stroke-width=\"1.2\"").collect();
        assert_eq!(ticks.len(), 2, "one tick per colour, and no more: {svg}");
        let theme = crate::Theme::light();
        assert!(
            svg.contains(&format!("stroke=\"{}\"", theme.color(0))),
            "first colour gone: {svg}"
        );
        assert!(
            svg.contains(&format!("stroke=\"{}\"", theme.color(1))),
            "second colour gone: {svg}"
        );
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
