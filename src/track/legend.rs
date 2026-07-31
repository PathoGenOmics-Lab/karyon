//! A key to the colours, as a band of its own.
//!
//! # Why this is a track
//!
//! A legend is a horizontal strip of the figure that carries no coordinates,
//! which is exactly what a [`Track`] is. Making it one means it stacks, clips
//! and lays itself out like everything else, and it goes where the caller puts
//! it rather than where a track decided to squeeze it.
//!
//! It also means the figure reserves room for it. A legend drawn into a corner
//! of a plot that already has data in that corner is a legend with a line
//! through it, which is the failure mode this exists to avoid.
//!
//! # Entries wrap
//!
//! Keys are laid across the band and wrap onto another row when they run out
//! of width, so [`Legend::height`] depends on how wide the figure is. Nothing
//! is ever dropped for want of room: a key that is not drawn is worse than a
//! legend that is two rows tall.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor, SvgWriter};
use crate::theme::{mix, wash, Theme};
use crate::track::{DrawContext, Track};

/// How one key is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Marker {
    /// A filled square, for anything drawn as an area.
    Box,
    /// A filled circle, for anything drawn as a point.
    Dot,
    /// A short stroke, for anything drawn as a line.
    Line,
    /// An empty square with an edge, for anything marked by being outlined
    /// rather than by being filled.
    Outline,
    /// A square washed in the colour and edged in it, which is how the crate
    /// draws a large filled shape. A solid swatch beside a washed arrow is a
    /// key that does not look like the thing it explains.
    Area,
}

/// One thing a legend explains.
#[derive(Debug, Clone, PartialEq)]
pub enum LegendItem {
    /// A colour and what it means.
    Key {
        /// What the colour stands for.
        label: String,
        /// The colour itself.
        color: String,
        /// How to draw the swatch.
        marker: Marker,
    },
    /// A continuous scale between two colours.
    ///
    /// For a quantity carried by shade rather than by hue: alignment identity,
    /// a distance, a depth. Drawn as a strip of steps rather than an SVG
    /// gradient, which keeps the output to elements every renderer agrees on.
    Ramp {
        /// What the scale measures.
        label: String,
        /// Colour at the low end.
        from: String,
        /// Colour at the high end.
        to: String,
        /// Number under the low end.
        low: String,
        /// Number under the high end.
        high: String,
    },
}

/// A key to what the colours in a figure mean.
///
/// ```
/// use karyon::{AxisTrack, Figure, Legend, LegendTrack, Region};
///
/// let legend = Legend::new()
///     .key("resistant", "#d55e00")
///     .key("susceptible", "#0072b2")
///     .ramp("identity", "#eeeeee", "#333333", "70%", "100%");
///
/// let svg = Figure::new(Region::new("chr1", 0, 1_000).unwrap())
///     .push(AxisTrack::new())
///     .push(LegendTrack::new(legend))
///     .to_svg();
/// assert!(svg.contains("resistant"));
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Legend {
    items: Vec<LegendItem>,
    swatch: f64,
    gap: f64,
}

impl Legend {
    /// An empty legend.
    pub fn new() -> Self {
        Legend {
            items: Vec::new(),
            swatch: 10.0,
            gap: 16.0,
        }
    }

    /// Adds a colour and what it means, drawn as a filled square.
    pub fn key(self, label: impl Into<String>, color: impl Into<String>) -> Self {
        self.marked(label, color, Marker::Box)
    }

    /// Adds a key drawn as a point, for a series plotted as points.
    pub fn dot(self, label: impl Into<String>, color: impl Into<String>) -> Self {
        self.marked(label, color, Marker::Dot)
    }

    /// Adds a key drawn as a stroke, for a series plotted as a line.
    pub fn line(self, label: impl Into<String>, color: impl Into<String>) -> Self {
        self.marked(label, color, Marker::Line)
    }

    /// Adds a key drawn the way a large filled shape is drawn: washed in the
    /// colour and edged in it.
    ///
    /// For gene arrows, matrix cells, ribbons: anything the figure fills. A
    /// solid swatch beside a washed shape is a key that does not look like the
    /// thing it explains.
    pub fn area(self, label: impl Into<String>, color: impl Into<String>) -> Self {
        self.marked(label, color, Marker::Area)
    }

    /// Adds a key drawn as an empty square with an edge.
    ///
    /// For anything the figure marks by outlining rather than by filling. A
    /// solid swatch for an outlined mark is a legend that does not match what
    /// it explains.
    pub fn outline(self, label: impl Into<String>, color: impl Into<String>) -> Self {
        self.marked(label, color, Marker::Outline)
    }

    /// Adds a key with the marker named.
    pub fn marked(
        mut self,
        label: impl Into<String>,
        color: impl Into<String>,
        marker: Marker,
    ) -> Self {
        self.items.push(LegendItem::Key {
            label: label.into(),
            color: color.into(),
            marker,
        });
        self
    }

    /// Adds a continuous scale between two colours.
    pub fn ramp(
        mut self,
        label: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        low: impl Into<String>,
        high: impl Into<String>,
    ) -> Self {
        self.items.push(LegendItem::Ramp {
            label: label.into(),
            from: from.into(),
            to: to.into(),
            low: low.into(),
            high: high.into(),
        });
        self
    }

    /// Sets the size of one swatch in pixels.
    pub fn swatch(mut self, size: f64) -> Self {
        self.swatch = size.max(2.0);
        self
    }

    /// Sets the space between one key and the next.
    pub fn gap(mut self, gap: f64) -> Self {
        self.gap = gap.max(0.0);
        self
    }

    /// The items.
    pub fn items(&self) -> &[LegendItem] {
        &self.items
    }

    /// How many keys the legend holds.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the legend explains nothing.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// How wide one item needs to be.
    fn item_width(&self, item: &LegendItem, font: f64) -> f64 {
        match item {
            LegendItem::Key { label, .. } => self.swatch + 5.0 + text_width(label, font),
            LegendItem::Ramp {
                label, low, high, ..
            } => {
                text_width(label, font)
                    + 5.0
                    + text_width(low, font)
                    + 3.0
                    + RAMP_WIDTH
                    + 3.0
                    + text_width(high, font)
            }
        }
    }

    /// The rows the items fall into at this width, each as `(x, item index)`.
    ///
    /// Wrapping rather than truncating: a key that is not drawn is worse than a
    /// legend that is two rows tall, because the reader cannot tell that
    /// anything is missing.
    pub fn rows(&self, width: f64, font: f64) -> Vec<Vec<(f64, usize)>> {
        let mut rows: Vec<Vec<(f64, usize)>> = Vec::new();
        let mut current: Vec<(f64, usize)> = Vec::new();
        let mut x = 0.0f64;
        for (index, item) in self.items.iter().enumerate() {
            let needed = self.item_width(item, font);
            if !current.is_empty() && x + needed > width {
                rows.push(std::mem::take(&mut current));
                x = 0.0;
            }
            current.push((x, index));
            x += needed + self.gap;
        }
        if !current.is_empty() {
            rows.push(current);
        }
        rows
    }

    /// How tall the legend is at this width.
    pub fn height(&self, width: f64, theme: &Theme) -> f64 {
        if self.items.is_empty() {
            return 0.0;
        }
        let font = theme.font_size - 1.0;
        let row = self.swatch.max(font) + 6.0;
        self.rows(width, font).len() as f64 * row
    }

    /// Draws the legend into a rectangle whose top left is `(x, y)`.
    pub fn draw(&self, svg: &mut SvgWriter, x: f64, y: f64, width: f64, theme: &Theme) {
        if self.items.is_empty() {
            return;
        }
        let font = theme.font_size - 1.0;
        let row_height = self.swatch.max(font) + 6.0;

        for (row, entries) in self.rows(width, font).iter().enumerate() {
            let top = y + row as f64 * row_height;
            let middle = top + row_height / 2.0;
            let baseline = middle + font * 0.35;
            for (offset, index) in entries {
                let at = x + offset;
                match &self.items[*index] {
                    LegendItem::Key {
                        label,
                        color,
                        marker,
                    } => {
                        match marker {
                            Marker::Box => svg.rect_rounded(
                                at,
                                middle - self.swatch / 2.0,
                                self.swatch,
                                self.swatch,
                                theme.corner_radius,
                                color,
                            ),
                            Marker::Dot => {
                                svg.circle(at + self.swatch / 2.0, middle, self.swatch / 2.5, color)
                            }
                            Marker::Line => svg.line(
                                at,
                                middle,
                                at + self.swatch,
                                middle,
                                color,
                                (self.swatch / 4.0).max(1.5),
                            ),
                            Marker::Area => {
                                svg.rect_rounded(
                                    at,
                                    middle - self.swatch / 2.0,
                                    self.swatch,
                                    self.swatch,
                                    theme.corner_radius,
                                    &wash(color, theme),
                                );
                                svg.rect_outline(
                                    at + 0.55,
                                    middle - self.swatch / 2.0 + 0.55,
                                    self.swatch - 1.1,
                                    self.swatch - 1.1,
                                    color,
                                    1.1,
                                );
                            }
                            Marker::Outline => svg.rect_outline(
                                at + 0.75,
                                middle - self.swatch / 2.0 + 0.75,
                                self.swatch - 1.5,
                                self.swatch - 1.5,
                                color,
                                1.5,
                            ),
                        }
                        svg.text(
                            at + self.swatch + 5.0,
                            baseline,
                            label,
                            &theme.muted,
                            font,
                            Anchor::Start,
                        );
                    }
                    LegendItem::Ramp {
                        label,
                        from,
                        to,
                        low,
                        high,
                    } => {
                        let mut cursor = at;
                        svg.text(cursor, baseline, label, &theme.muted, font, Anchor::Start);
                        cursor += text_width(label, font) + 5.0;
                        svg.text(cursor, baseline, low, &theme.muted, font, Anchor::Start);
                        cursor += text_width(low, font) + 3.0;
                        // Steps rather than a gradient, so the output stays to
                        // elements every renderer draws the same way.
                        let step = RAMP_WIDTH / RAMP_STEPS as f64;
                        for index in 0..RAMP_STEPS {
                            let t = index as f64 / (RAMP_STEPS - 1) as f64;
                            svg.rect(
                                cursor + index as f64 * step,
                                middle - self.swatch / 2.0,
                                step + 0.5,
                                self.swatch,
                                &mix(from, to, t),
                            );
                        }
                        cursor += RAMP_WIDTH + 3.0;
                        svg.text(cursor, baseline, high, &theme.muted, font, Anchor::Start);
                    }
                }
            }
        }
    }
}

/// How wide the coloured strip of a ramp is drawn.
const RAMP_WIDTH: f64 = 56.0;
/// How many steps the strip is made of.
const RAMP_STEPS: usize = 24;

/// A [`Legend`] as a band of the figure.
#[derive(Debug, Clone)]
pub struct LegendTrack {
    legend: Legend,
    label: Option<String>,
}

impl LegendTrack {
    /// A band drawing `legend`.
    pub fn new(legend: Legend) -> Self {
        LegendTrack {
            legend,
            label: None,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// The legend.
    pub fn legend(&self) -> &Legend {
        &self.legend
    }
}

impl Track for LegendTrack {
    fn height(&self, scale: &Scale) -> f64 {
        // The theme reaches the height only through the font size, and the
        // default is what the figure will use unless the caller changed it, as
        // in the feature track.
        self.legend
            .height(scale.width(), &Theme::default())
            .max(1.0)
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        self.legend.draw(ctx.svg, band.x, band.y, band.w, ctx.theme);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 0, 1_000).unwrap()
    }

    fn legend() -> Legend {
        Legend::new()
            .key("resistant", "#d55e00")
            .dot("susceptible", "#0072b2")
            .line("mean depth", "#009e73")
    }

    #[test]
    fn an_empty_legend_takes_no_room_and_draws_nothing() {
        let theme = Theme::light();
        assert_eq!(Legend::new().height(500.0, &theme), 0.0);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LegendTrack::new(Legend::new()))
            .to_svg();
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn keys_wrap_rather_than_being_dropped() {
        let theme = Theme::light();
        let mut wide = Legend::new();
        for index in 0..12 {
            wide = wide.key(format!("category number {index}"), "#0072b2");
        }
        assert_eq!(wide.rows(2_000.0, 10.0).len(), 1, "room for one row");
        let narrow = wide.rows(300.0, 10.0);
        assert!(narrow.len() > 3, "wrapped onto {} rows", narrow.len());
        // Every key survives the wrapping.
        assert_eq!(narrow.iter().map(|row| row.len()).sum::<usize>(), 12);
        assert!(wide.height(300.0, &theme) > wide.height(2_000.0, &theme));
    }

    #[test]
    fn one_key_too_wide_for_the_band_still_gets_its_own_row() {
        let long = Legend::new().key("a".repeat(400), "#0072b2");
        assert_eq!(long.rows(100.0, 10.0).len(), 1);
    }

    #[test]
    fn an_area_key_is_washed_and_edged_like_the_shape_it_explains() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LegendTrack::new(
                Legend::new().area("gene family", "#0072b2"),
            ))
            .to_svg();
        let theme = Theme::light();
        assert!(svg.contains(&wash("#0072b2", &theme)), "no wash");
        assert!(svg.contains("stroke=\"#0072b2\""), "no edge");
    }

    #[test]
    fn an_outlined_key_is_drawn_as_an_outline() {
        // A solid swatch for a mark the figure draws as an outline is a legend
        // that does not match what it explains.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LegendTrack::new(
                Legend::new().outline("unmatched", "#111111"),
            ))
            .to_svg();
        assert!(svg.contains("fill=\"none\" stroke=\"#111111\""), "{svg}");
    }

    #[test]
    fn every_marker_shape_is_drawn() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LegendTrack::new(legend()))
            .to_svg();
        assert!(svg.contains("#d55e00"), "the box");
        assert!(svg.contains("<circle"), "the dot");
        assert!(svg.contains("<line"), "the stroke");
        for word in ["resistant", "susceptible", "mean depth"] {
            assert!(svg.contains(word), "missing {word}");
        }
    }

    #[test]
    fn a_ramp_is_drawn_as_steps_between_its_two_ends() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LegendTrack::new(
                Legend::new().ramp("identity", "#ffffff", "#000000", "70%", "100%"),
            ))
            .to_svg();
        assert!(svg.contains(">identity</text>"));
        assert!(svg.contains(">70%</text>") && svg.contains(">100%</text>"));
        // A real ramp: both ends reached, and many shades between them.
        assert!(svg.contains("#000000"), "the dark end");
        let mut shades: Vec<&str> = svg
            .match_indices("<rect ")
            .filter_map(|(at, _)| {
                let rest = &svg[at..];
                let head = &rest[..rest.find("/>")?];
                let fill = head.find("fill=\"")? + 6;
                Some(&head[fill..fill + 7])
            })
            .filter(|shade| shade.starts_with('#'))
            .collect();
        shades.sort_unstable();
        shades.dedup();
        assert!(
            shades.len() > 15,
            "only {} shades: {shades:?}",
            shades.len()
        );
    }

    #[test]
    fn the_band_is_as_tall_as_the_legend_needs() {
        let scale = Scale::new(&region(), 0.0, 900.0);
        let narrow = Scale::new(&region(), 0.0, 120.0);
        let track = LegendTrack::new(legend());
        assert!(track.height(&narrow) > track.height(&scale));
        assert!(track.height(&scale) > 0.0);
    }
}
