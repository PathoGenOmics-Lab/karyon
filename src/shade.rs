//! A stretch of the axis shaded across every band laid on it.
//!
//! A genome browser calls it a region of interest: the deletion a sample lost,
//! the gene a scan peaks over, the weeks a lineage rose in, marked once and
//! read down the whole stack. It belongs to the figure and not to a track,
//! because it crosses bands and a track only ever sees its own, clipped. So it
//! is handed to [`Figure::shade`](crate::Figure::shade), in the figure's own
//! coordinates as every track's are, and the figure draws it.
//!
//! # Behind the data, edged over it
//!
//! A wash drawn over the tracks tints every mark beneath it, so a heatmap cell
//! or a lineage colour under a shade stops matching its key. The wash goes
//! behind the tracks instead, where it shows through everything a track leaves
//! empty. A track that fills its band, a heatmap above all, leaves nothing for
//! it to show through but the gaps between rows, so the two ends of the stretch
//! are drawn again over the tracks as dashed hairlines, which mark it without
//! changing a single colour a reader looks up.
//!
//! # One colour for every shade
//!
//! A shade says look here, not which category. A palette colour dealt to each
//! one would read as a key the figure does not have, and would collide with the
//! categories the tracks already paint, so every shade is the foreground ink at
//! 8% of its strength unless [`Shade::color`] says otherwise.

use crate::region::Region;
use crate::scale::Scale;
use crate::style::LinePattern;
use crate::svg::{fit_text_by, text_width, Anchor, SvgWriter};
use crate::theme::Theme;
use crate::track::feature::span_label;

/// How strong the wash behind a shaded stretch is.
///
/// Measured against the heaviest backdrop it sits under: at 8% the foreground
/// ink reads as a band on the page in both themes and leaves a coverage area
/// or a feature drawn across it as dark as it is anywhere else.
const SHADE_OPACITY: f64 = 0.08;

/// The narrowest a shade is drawn. One base of a ten kilobase window is a
/// tenth of a pixel, and a stretch nobody can see is not marked.
const SHADE_MIN_PX: f64 = 2.0;

/// Under this a shade has one edge, down its middle. A one-base shade at 4.2
/// kilobases over 720 pixels drew two dashed edges two pixels apart, which
/// read as one smudged line.
const ONE_EDGE_UNDER_PX: f64 = 6.0;

/// A stretch of the figure's axis, shaded across every band laid on it.
///
/// ```
/// use karyon::{plot, Shade};
///
/// let svg = plot("chr1:1-10,000")
///     .unwrap()
///     .add_coverage(vec![30.0; 10_000])
///     .label("depth")
///     .shade(Shade::new(4_000, 5_000).name("deletion"))
///     .to_svg();
///
/// assert!(svg.contains("<title>deletion, 4,001 to 5,000</title>"));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Shade {
    start: u64,
    end: u64,
    name: Option<String>,
    color: Option<String>,
    /// What the tooltip and the alt text call the stretch, where the figure's
    /// coordinates are not the ones a reader knows it by.
    described: Option<String>,
}

impl Shade {
    /// The stretch `start..end`, 0-based and half-open, as every track's
    /// coordinates are.
    ///
    /// A reversed pair is swapped, and an empty one is widened to one base, as
    /// [`CladeBlock::new`](crate::CladeBlock::new) and
    /// [`IdeogramTrack::highlight`](crate::IdeogramTrack::highlight) take
    /// theirs: a shade with nothing in it would mark nothing, and one running
    /// backwards is the same stretch written the other way.
    ///
    /// ```
    /// use karyon::Shade;
    ///
    /// let shade = Shade::new(10, 5);
    /// assert_eq!((shade.start(), shade.end()), (5, 10));
    /// let shade = Shade::new(7, 7);
    /// assert_eq!((shade.start(), shade.end()), (7, 8));
    /// ```
    pub fn new(start: u64, end: u64) -> Shade {
        let (start, end) = if end < start {
            (end, start)
        } else {
            (start, end)
        };
        Shade {
            start,
            end: end.max(start.saturating_add(1)),
            name: None,
            color: None,
            described: None,
        }
    }

    /// Names the stretch, at the head of its column.
    ///
    /// The figure gives the names a row of their own above the first band laid
    /// on the coordinates, and only while a named shade is in view. A name too
    /// long for the room before the next named shade is cut short there, and
    /// one with no room at all is left to the tooltip and the alt text, which
    /// always say it whole. An empty name is no name.
    pub fn name(mut self, name: impl Into<String>) -> Shade {
        let name = name.into();
        self.name = (!name.is_empty()).then_some(name);
        self
    }

    /// Sets the colour of the wash and its edges, which are otherwise the
    /// figure's foreground ink.
    ///
    /// One colour for every shade is the default, because a shade says look
    /// here and not which category, and a palette colour would read as a key
    /// the figure does not have. Give each one its own where the shades do
    /// stand for categories, and say what they are in the caption.
    pub fn color(mut self, color: impl Into<String>) -> Shade {
        self.color = Some(color.into());
        self
    }

    /// The first base shaded, 0-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// One past the last base shaded.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// What the tooltip and the alt text call the stretch in place of its
    /// coordinates, for an axis a reader does not count along: a scan across
    /// a whole genome lays its sequences end to end, and `7:1,001-2,000` is
    /// the place, where the shared axis would say 1,520,301,001.
    pub(crate) fn described(mut self, text: String) -> Shade {
        self.described = Some(text);
        self
    }

    /// Whether any of the stretch is inside `region`.
    pub(crate) fn touches(&self, region: &Region) -> bool {
        self.end > region.start() && self.start < region.end()
    }

    /// What a reader is told of it: its name, and where it is.
    pub(crate) fn said(&self) -> String {
        let span = self
            .described
            .clone()
            .unwrap_or_else(|| span_label(self.start, self.end));
        match &self.name {
            Some(name) => format!("{name}, {span}"),
            None => span,
        }
    }
}

/// One shade as it is drawn: its column in pixels, inside the plotting area,
/// and the edges of it that fall inside the window.
pub(crate) struct Column<'a> {
    shade: &'a Shade,
    /// The left of the column, widened and held inside the plotting area.
    left: f64,
    /// The right of it.
    right: f64,
    /// Where its dashed edges go: both ends where each is in the window, or
    /// one line down its middle where the column is too narrow for two.
    edges: Vec<f64>,
}

/// The shades of `shades` that touch the window, as columns over `scale`,
/// held between `plot_left` and `plot_right` so that a stretch starting
/// before the window never reaches into the label gutter or an axis strip.
pub(crate) fn columns<'a>(
    shades: &'a [Shade],
    region: &Region,
    scale: &Scale,
    plot_left: f64,
    plot_right: f64,
) -> Vec<Column<'a>> {
    let mut columns: Vec<Column<'a>> = shades
        .iter()
        .filter(|shade| shade.touches(region))
        .map(|shade| {
            let mut left = scale.x(shade.start.max(region.start()));
            let mut right = scale.x(shade.end.min(region.end()));
            let narrow = right - left < ONE_EDGE_UNDER_PX;
            let middle = (left + right) / 2.0;
            if right - left < SHADE_MIN_PX {
                left = middle - SHADE_MIN_PX / 2.0;
                right = middle + SHADE_MIN_PX / 2.0;
            }
            let left = left.max(plot_left);
            let right = right.min(plot_right);
            // An edge marks where the stretch stops, so an end beyond the
            // window has none: the stretch goes on past what is drawn.
            let edges = if narrow {
                vec![middle.clamp(plot_left, plot_right)]
            } else {
                let mut edges = Vec::with_capacity(2);
                if shade.start > region.start() {
                    edges.push(scale.x(shade.start));
                }
                if shade.end < region.end() {
                    edges.push(scale.x(shade.end));
                }
                edges
            };
            Column {
                shade,
                left,
                right,
                edges,
            }
        })
        .filter(|column| column.right > column.left)
        .collect();
    columns.sort_by(|a, b| a.left.total_cmp(&b.left));
    columns
}

/// The wash of every column down every run of bands, behind the tracks.
///
/// One group a shade, titled with its name and its span, so a pointer resting
/// anywhere in the stretch where no track has drawn is told what it is, and
/// the data drawn over it keeps answering for itself.
pub(crate) fn wash(
    svg: &mut SvgWriter,
    columns: &[Column<'_>],
    runs: &[(f64, f64)],
    theme: &Theme,
) {
    for column in columns {
        let color = column.shade.color.as_deref().unwrap_or(&theme.foreground);
        svg.begin_titled(&column.shade.said());
        for (top, bottom) in runs {
            svg.rect_opacity(
                column.left,
                *top,
                column.right - column.left,
                bottom - top,
                color,
                SHADE_OPACITY,
            );
        }
        svg.end_group();
    }
}

/// The edges of every column down every run, over the tracks: dashed
/// hairlines, which mark the stretch across a band that fills itself and
/// change no colour a reader looks up. Inert, so the marks under them keep
/// their tooltips.
///
/// `first_top` is where the first run's edges begin, under the row of names
/// rather than through it.
pub(crate) fn edges(
    svg: &mut SvgWriter,
    columns: &[Column<'_>],
    runs: &[(f64, f64)],
    first_top: f64,
    theme: &Theme,
) {
    if columns.iter().all(|column| column.edges.is_empty()) {
        return;
    }
    svg.begin_titled_inert("");
    for column in columns {
        let color = column.shade.color.as_deref().unwrap_or(&theme.muted);
        for (index, (top, bottom)) in runs.iter().enumerate() {
            let top = if index == 0 { first_top } else { *top };
            for x in &column.edges {
                svg.line_pattern(
                    *x,
                    top,
                    *x,
                    *bottom,
                    color,
                    theme.tokens.hairline,
                    LinePattern::Dashed,
                );
            }
        }
    }
    svg.end_group();
}

/// The names of the named columns, in the row the figure keeps for them
/// between `top` and `top + height`, each from the left of its column to the
/// left of the next named one, or to `plot_right`.
pub(crate) fn names(
    svg: &mut SvgWriter,
    columns: &[Column<'_>],
    top: f64,
    height: f64,
    plot_right: f64,
    theme: &Theme,
) {
    let named: Vec<(&Column<'_>, &str)> = columns
        .iter()
        .filter_map(|column| column.shade.name.as_deref().map(|name| (column, name)))
        .collect();
    let size = theme.font_size - 1.0;
    let gap = 3.0;
    for (index, (column, name)) in named.iter().enumerate() {
        let from = column.left + gap;
        let to = named
            .get(index + 1)
            .map_or(plot_right, |(next, _)| next.left)
            - gap * 2.0;
        let visible = fit_text_by(name, to - from, |text| text_width(text, size));
        // An ellipsis alone names nothing, and the tooltip names it whole.
        if visible.is_empty() || visible == "\u{2026}" {
            continue;
        }
        svg.text(
            from,
            top + height / 2.0 + size * 0.35,
            &visible,
            &theme.muted,
            size,
            Anchor::Start,
        );
    }
}

/// Whether any of `shades` has a name and touches the window, which is what
/// the row of names is laid out for.
pub(crate) fn any_named(shades: &[Shade], region: &Region) -> bool {
    shades
        .iter()
        .any(|shade| shade.name.is_some() && shade.touches(region))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reversed_or_empty_shade_is_one_base_at_least() {
        let reversed = Shade::new(10, 5);
        assert_eq!((reversed.start(), reversed.end()), (5, 10));
        let empty = Shade::new(7, 7);
        assert_eq!((empty.start(), empty.end()), (7, 8));
        // At the top of the range the base saturates rather than wrapping.
        let last = Shade::new(u64::MAX, u64::MAX);
        assert_eq!((last.start(), last.end()), (u64::MAX, u64::MAX));
    }

    #[test]
    fn an_empty_name_is_no_name_and_a_shade_says_its_span_counted_from_one() {
        assert_eq!(Shade::new(1_000, 2_000).name("").said(), "1,001 to 2,000");
        assert_eq!(
            Shade::new(1_000, 2_000).name("deletion").said(),
            "deletion, 1,001 to 2,000"
        );
        assert_eq!(
            Shade::new(1_000, 2_000)
                .name("deletion")
                .described("7:1,001-2,000".to_string())
                .said(),
            "deletion, 7:1,001-2,000"
        );
    }

    #[test]
    fn a_column_stays_inside_the_plotting_area_and_has_an_edge_only_where_it_ends_in_view() {
        let region = Region::new("chr1", 1_000, 2_000).unwrap();
        let scale = Scale::new(&region, 100.0, 1_000.0);
        let shades = [Shade::new(500, 1_500)];
        let drawn = columns(&shades, &region, &scale, 100.0, 1_100.0);
        assert_eq!(drawn.len(), 1);
        assert_eq!((drawn[0].left, drawn[0].right), (100.0, 600.0));
        assert_eq!(drawn[0].edges, vec![600.0], "no edge at the window's start");
        // Outside the window is nothing at all.
        let away = [Shade::new(5_000, 6_000)];
        assert!(columns(&away, &region, &scale, 100.0, 1_100.0).is_empty());
    }
}
