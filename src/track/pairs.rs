//! Two places at once: linkage between variants, contacts between the bins of
//! a chromosome, epistasis between sites, a loop from an enhancer to its gene.
//!
//! A [`Pair`] is two stretches of the sequence and a value measured between
//! them. [`PairTrack`] draws a set of them one of two ways, and which way is a
//! question about how many there are rather than about what they mean.
//!
//! # A triangle, where most pairs were measured
//!
//! Linkage disequilibrium between every two variants of a window, or contacts
//! between every two bins, is a matrix, and half of it says everything since
//! the pair of `a` and `b` is the pair of `b` and `a`. Turned a quarter and
//! hung under the axis, that half is a triangle: each pair is a cell below the
//! point half way between its two places, as deep as they are far apart, so a
//! block of variants inherited together is a dark triangle under the stretch
//! it covers. Each variant is given the stretch from half way to the one
//! before it to half way to the one after, so the cells tile the triangle
//! with no gaps at any spacing, and each still sits under its own place.
//!
//! # Arcs, where a few were
//!
//! A handful of pairs spread along a gene, a set of epistatic sites or the
//! loops a caller found, would be a few cells lost in an empty triangle. Each
//! is an arc from one place to the other instead, as tall as they are far
//! apart, the strongest drawn last so none is hidden under a weaker one.
//!
//! In both, the colour is the value, on a ramp from nought to the largest
//! value drawn, which [`Track::key`] hands to the figure's key.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::text_rounded;
use crate::theme::{mix, Theme};
use crate::track::{arc_path, DrawContext, Track};

/// How far off the page the bottom of the ramp sits, as a matrix has it.
const ZERO_TINT: f64 = 0.1;

/// Above this many cells a triangle carries no tooltip on each: a thousand
/// variants are half a million cells, and a figure that grew by a title on
/// every one of them would be too big to open.
const TITLED_CELLS: usize = 2_500;

/// Two stretches of the sequence and a value measured between them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pair {
    /// The first stretch, 0-based and half-open.
    pub first: (u64, u64),
    /// The second.
    pub second: (u64, u64),
    /// What was measured between them: an r², a contact count, a score. A
    /// value that is not a number is a pair with no answer, and is not drawn.
    pub value: f64,
}

impl Pair {
    /// Two single bases, at 0-based `first` and `second`, as two variants in
    /// linkage are.
    pub fn new(first: u64, second: u64, value: f64) -> Self {
        Pair {
            first: (first, first.saturating_add(1)),
            second: (second, second.saturating_add(1)),
            value,
        }
    }

    /// Two stretches, each 0-based and half-open, as two bins of a contact
    /// map or the two anchors of a loop are.
    pub fn spans(first: (u64, u64), second: (u64, u64), value: f64) -> Self {
        Pair {
            first,
            second,
            value,
        }
    }

    /// The two stretches, the one that starts first first.
    fn ordered(&self) -> ((u64, u64), (u64, u64)) {
        if self.second.0 < self.first.0 {
            (self.second, self.first)
        } else {
            (self.first, self.second)
        }
    }

    /// Whether this pair is a single base paired with another, rather than
    /// two stretches.
    fn is_points(&self) -> bool {
        self.first.1 == self.first.0 + 1 && self.second.1 == self.second.0 + 1
    }
}

/// How a [`PairTrack`] draws its pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PairStyle {
    /// Every pair a cell of a triangle hung under the axis, for pairs
    /// measured between most of their places: linkage, contacts.
    #[default]
    Triangle,
    /// Every pair an arc from one place to the other, for a few pairs far
    /// apart: epistasis, co-occurrence, loops.
    Arcs,
}

impl PairStyle {
    /// The style a set of pairs is best drawn in: a triangle where most
    /// places were measured against the next place along, and arcs where
    /// few were.
    ///
    /// Next to each other is where a triangle is full. Linkage within a
    /// window, as PLINK writes it by default, and every contact map measure
    /// each place against its neighbours, however far out they stop; a set of
    /// epistatic sites or of loops joins places far apart and seldom the
    /// next one. Counting all the pairs the places could make instead called
    /// a hundred variants in PLINK's window of ten a sparse set.
    pub fn for_pairs(pairs: &[Pair]) -> PairStyle {
        let mut places: Vec<(u64, u64)> = pairs
            .iter()
            .flat_map(|pair| [pair.first, pair.second])
            .collect();
        places.sort_unstable();
        places.dedup();
        if places.len() < 3 {
            return PairStyle::Arcs;
        }
        let measured: std::collections::BTreeSet<((u64, u64), (u64, u64))> =
            pairs.iter().map(Pair::ordered).collect();
        let next = places
            .windows(2)
            .filter(|two| measured.contains(&(two[0], two[1])))
            .count();
        if next * 2 >= places.len() - 1 {
            PairStyle::Triangle
        } else {
            PairStyle::Arcs
        }
    }
}

/// Pairs of places along the sequence, as a triangle or as arcs.
///
/// ```
/// use karyon::{Figure, Pair, PairTrack, Region};
///
/// // Three variants, the first two in strong linkage.
/// let pairs = vec![
///     Pair::new(1_000, 1_400, 0.92),
///     Pair::new(1_000, 2_600, 0.10),
///     Pair::new(1_400, 2_600, 0.15),
/// ];
/// let svg = Figure::new(Region::new("chr1", 800, 3_000).unwrap())
///     .push(PairTrack::new(pairs).label("r²"))
///     .to_svg();
/// assert!(svg.contains("1,001 and 1,401: 0.92"));
/// ```
#[derive(Debug, Clone)]
pub struct PairTrack {
    pairs: Vec<Pair>,
    style: PairStyle,
    label: Option<String>,
    height: Option<f64>,
    ceiling: Option<f64>,
    floor: Option<f64>,
    hue: Option<String>,
    log: bool,
}

impl PairTrack {
    /// A track of `pairs`, drawn as a triangle until told otherwise.
    pub fn new(pairs: impl Into<Vec<Pair>>) -> Self {
        PairTrack {
            pairs: pairs.into(),
            style: PairStyle::Triangle,
            label: None,
            height: None,
            ceiling: None,
            floor: None,
            hue: None,
            log: false,
        }
    }

    /// Colours on a log scale, as a contact map is read: contacts fall by
    /// orders of magnitude with distance, and on a straight ramp the
    /// diagonal is the only thing with a colour.
    pub fn log_scale(mut self, log: bool) -> Self {
        self.log = log;
        self
    }

    /// Chooses a triangle or arcs.
    pub fn style(mut self, style: PairStyle) -> Self {
        self.style = style;
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels. A triangle is as deep as its farthest
    /// pair is wide until told, up to 240 pixels, and arcs are 90 tall.
    pub fn height(mut self, pixels: f64) -> Self {
        if pixels.is_finite() && pixels > 0.0 {
            self.height = Some(pixels.max(16.0));
        }
        self
    }

    /// The value the ramp saturates at, instead of the largest one drawn: 1
    /// for an r², so a window whose strongest linkage is 0.4 does not look
    /// complete.
    pub fn ceiling(mut self, value: f64) -> Self {
        if value.is_finite() && value > 0.0 {
            self.ceiling = Some(value);
        }
        self
    }

    /// Draws only the pairs whose value is at least this.
    pub fn threshold(mut self, value: f64) -> Self {
        if value.is_finite() {
            self.floor = Some(value);
        }
        self
    }

    /// The hue at the top of the ramp, the theme's accent by default.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.hue = Some(color.into());
        self
    }

    /// The pairs, as they were given.
    pub fn pairs(&self) -> &[Pair] {
        &self.pairs
    }

    /// The pairs that are drawn: a value, above the threshold, and both
    /// places in reach of the window from `low` to `high`.
    fn drawn(&self, (low, high): (u64, u64)) -> Vec<&Pair> {
        let touches = |(start, end): (u64, u64)| end > start && start < high && end > low;
        self.pairs
            .iter()
            .filter(|pair| pair.value.is_finite())
            .filter(|pair| self.floor.map_or(true, |floor| pair.value >= floor))
            .filter(|pair| touches(pair.first) && touches(pair.second))
            .collect()
    }

    /// The value the ramp tops out at.
    fn top(&self, drawn: &[&Pair]) -> Option<f64> {
        self.ceiling.or_else(|| {
            drawn
                .iter()
                .map(|pair| pair.value)
                .fold(None, |top: Option<f64>, value| {
                    Some(top.map_or(value, |t| t.max(value)))
                })
                .filter(|top| *top > 0.0)
        })
    }

    fn color_of(&self, value: f64, top: f64, theme: &Theme) -> String {
        let hue = self.hue.clone().unwrap_or_else(|| theme.accent.clone());
        let fraction = if self.log {
            value.max(0.0).ln_1p() / top.ln_1p()
        } else {
            value / top
        }
        .clamp(0.0, 1.0);
        mix(
            theme.surface(),
            &hue,
            ZERO_TINT + (1.0 - ZERO_TINT) * fraction,
        )
    }

    /// The stretch each place owns in a triangle: a stretch as given, and a
    /// single base from half way to the place before it to half way to the
    /// one after, so cells tile the triangle whatever the spacing.
    fn cells(drawn: &[&Pair]) -> Vec<((u64, u64), (f64, f64))> {
        let mut places: Vec<(u64, u64)> = drawn
            .iter()
            .flat_map(|pair| [pair.first, pair.second])
            .collect();
        places.sort_unstable();
        places.dedup();
        if !drawn.iter().all(|pair| pair.is_points()) {
            return places
                .into_iter()
                .map(|place| (place, (place.0 as f64, place.1 as f64)))
                .collect();
        }
        let centres: Vec<f64> = places.iter().map(|place| place.0 as f64 + 0.5).collect();
        (0..places.len())
            .map(|index| {
                let here = centres[index];
                let before = index
                    .checked_sub(1)
                    .map_or(here - 0.5, |previous| (centres[previous] + here) / 2.0);
                let after = centres
                    .get(index + 1)
                    .map_or(here + 0.5, |next| (here + next) / 2.0);
                // The ends reach as far past the last place as the half
                // step before it, so the outermost cells are not slivers.
                let before = if index == 0 {
                    here - (after - here)
                } else {
                    before
                };
                let after = if index + 1 == places.len() {
                    here + (here - before)
                } else {
                    after
                };
                (places[index], (before, after))
            })
            .collect()
    }

    /// How deep a triangle of these pairs is at this scale, before it is fitted
    /// to a band: half the width its farthest pair spans.
    fn natural_depth(&self, scale: &Scale) -> f64 {
        let drawn = self.drawn(scale.bounds());
        let cells = Self::cells(&drawn);
        let extent = |place: (u64, u64)| {
            cells
                .iter()
                .find(|(known, _)| *known == place)
                .map_or((place.0 as f64, place.1 as f64), |(_, span)| *span)
        };
        drawn
            .iter()
            .map(|pair| {
                let (a, b) = pair.ordered();
                (scale.x_at(extent(b).1) - scale.x_at(extent(a).0)) / 2.0
            })
            .fold(0.0, f64::max)
    }

    fn draw_triangle(&self, ctx: &mut DrawContext<'_>, drawn: &[&Pair], top: f64) {
        let cells = Self::cells(drawn);
        let extent = |place: (u64, u64)| {
            cells
                .iter()
                .find(|(known, _)| *known == place)
                .map_or((place.0 as f64, place.1 as f64), |(_, span)| *span)
        };
        let depth = self.natural_depth(ctx.scale).max(1.0);
        let room = (ctx.band.h - ctx.px(4.0)).max(1.0);
        // Squeezed to the band where the farthest pair is deeper than it,
        // and never stretched: a cell's shape then says the pair is short.
        let squeeze = (room / depth).min(1.0);
        let baseline = ctx.band.y + ctx.px(2.0);
        let at = |u: f64, v: f64| {
            let (xu, xv) = (ctx.scale.x_at(u), ctx.scale.x_at(v));
            ((xu + xv) / 2.0, baseline + (xv - xu).abs() / 2.0 * squeeze)
        };
        let titled = drawn.len() <= TITLED_CELLS;
        // Weakest first, so the colour of a strong pair is the one on top
        // where two cells share an edge.
        let mut order: Vec<&&Pair> = drawn.iter().collect();
        order.sort_by(|a, b| a.value.total_cmp(&b.value));
        for pair in order {
            let (a, b) = pair.ordered();
            let ((sa, ea), (sb, eb)) = (extent(a), extent(b));
            let corners = if a == b {
                vec![at(sa, sa), at(ea, ea), at(sa, ea)]
            } else {
                vec![at(ea, sb), at(ea, eb), at(sa, eb), at(sa, sb)]
            };
            if titled {
                ctx.svg.begin_titled(&pair_title(pair));
            }
            ctx.svg
                .polygon(&corners, &self.color_of(pair.value, top, ctx.theme));
            if titled {
                ctx.svg.end_group();
            }
        }
    }

    fn draw_arcs(&self, ctx: &mut DrawContext<'_>, drawn: &[&Pair], top: f64) {
        let baseline = ctx.band.bottom() - ctx.px(1.0);
        let room = (ctx.band.h - ctx.px(6.0)).max(1.0);
        let middle = |(start, end): (u64, u64)| ctx.scale.x_at((start + end) as f64 / 2.0);
        let widest = drawn
            .iter()
            .map(|pair| (middle(pair.second) - middle(pair.first)).abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let mut order: Vec<&&Pair> = drawn.iter().collect();
        order.sort_by(|a, b| a.value.total_cmp(&b.value));
        for pair in order {
            let (a, b) = pair.ordered();
            let (x0, x1) = (middle(a), middle(b));
            // As tall as the pair is wide against the widest, so two arcs
            // of one length are one height and a near pair is a low one.
            let apex = baseline - room * ((x1 - x0) / widest).clamp(0.08, 1.0);
            let fraction = if self.log {
                pair.value.max(0.0).ln_1p() / top.ln_1p()
            } else {
                pair.value / top
            }
            .clamp(0.0, 1.0);
            let width = ctx.px(0.8 + 2.4 * fraction);
            ctx.svg.begin_titled(&pair_title(pair));
            ctx.svg.path_stroked(
                &arc_path(x0, x1, baseline, apex),
                &self.color_of(pair.value, top, ctx.theme),
                width,
            );
            ctx.svg.end_group();
        }
    }
}

/// What a reader hovering a pair is told: both places, counted from one as
/// the ruler counts them, and the value.
fn pair_title(pair: &Pair) -> String {
    let place = |(start, end): (u64, u64)| {
        if end == start + 1 {
            crate::track::axis::group_thousands(start + 1)
        } else {
            format!(
                "{}-{}",
                crate::track::axis::group_thousands(start + 1),
                crate::track::axis::group_thousands(end)
            )
        }
    };
    let (a, b) = pair.ordered();
    format!(
        "{} and {}: {}",
        place(a),
        place(b),
        text_rounded(pair.value, 4)
    )
}

impl Track for PairTrack {
    fn noun(&self) -> &str {
        match self.style {
            PairStyle::Triangle => "a triangle of pairs",
            PairStyle::Arcs => "arcs between pairs",
        }
    }

    fn height(&self, scale: &Scale) -> f64 {
        if let Some(height) = self.height {
            return height;
        }
        match self.style {
            PairStyle::Triangle => {
                let depth = self.natural_depth(scale);
                (depth + 4.0).clamp(24.0, 240.0).ceil()
            }
            PairStyle::Arcs => 90.0,
        }
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn key(
        &self,
        region: &Region,
        _px_per_bp: f64,
        theme: &Theme,
    ) -> Option<crate::track::legend::Legend> {
        let drawn = self.drawn((region.start(), region.end()));
        let top = self.top(&drawn)?;
        let hue = self.hue.clone().unwrap_or_else(|| theme.accent.clone());
        let named = self.label.clone().unwrap_or_else(|| "value".to_string());
        Some(crate::track::legend::Legend::new().ramp(
            if self.log {
                format!("{named}, log scale")
            } else {
                named
            },
            mix(theme.surface(), &hue, ZERO_TINT),
            hue,
            "0",
            text_rounded(top, 2),
        ))
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let drawn = self.drawn((ctx.region.start(), ctx.region.end()));
        let Some(top) = self.top(&drawn) else {
            return;
        };
        match self.style {
            PairStyle::Triangle => self.draw_triangle(ctx, &drawn, top),
            PairStyle::Arcs => self.draw_arcs(ctx, &drawn, top),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;

    fn region() -> Region {
        Region::new("chr1", 0, 1_000).unwrap()
    }

    #[test]
    fn a_dense_set_is_a_triangle_and_a_sparse_one_is_arcs() {
        let dense: Vec<Pair> = (0..6u64)
            .flat_map(|a| ((a + 1)..6).map(move |b| Pair::new(a * 100, b * 100, 0.5)))
            .collect();
        assert_eq!(PairStyle::for_pairs(&dense), PairStyle::Triangle);
        let sparse = vec![
            Pair::new(10, 900, 1.0),
            Pair::new(200, 600, 0.3),
            Pair::new(50, 70, 0.1),
            Pair::new(300, 800, 0.2),
        ];
        assert_eq!(PairStyle::for_pairs(&sparse), PairStyle::Arcs);
    }

    #[test]
    fn each_pair_is_drawn_and_named_by_its_places_from_one() {
        let pairs = vec![
            Pair::new(100, 400, 0.9),
            Pair::new(100, 800, 0.1),
            Pair::new(400, 800, 0.3),
        ];
        for style in [PairStyle::Triangle, PairStyle::Arcs] {
            let svg = Figure::new(region())
                .push(PairTrack::new(pairs.clone()).style(style))
                .to_svg();
            assert!(svg.contains("101 and 401: 0.9"), "{style:?}: {svg}");
            assert!(svg.contains("401 and 801: 0.3"), "{style:?}");
        }
    }

    #[test]
    fn a_pair_below_the_threshold_or_without_a_value_is_not_drawn() {
        let pairs = vec![
            Pair::new(100, 400, 0.9),
            Pair::new(100, 800, 0.1),
            Pair::new(400, 800, f64::NAN),
        ];
        let svg = Figure::new(region())
            .push(PairTrack::new(pairs).threshold(0.5))
            .to_svg();
        assert!(svg.contains("101 and 401"));
        assert!(!svg.contains("101 and 801"), "under the threshold");
        assert!(!svg.contains("NaN"), "no value is no cell");
    }

    #[test]
    fn the_cells_of_a_triangle_meet_whatever_the_spacing() {
        let drawn = [
            Pair::new(100, 130, 1.0),
            Pair::new(130, 700, 1.0),
            Pair::new(100, 700, 1.0),
        ];
        let refs: Vec<&Pair> = drawn.iter().collect();
        let cells = PairTrack::cells(&refs);
        for window in cells.windows(2) {
            assert_eq!(
                window[0].1 .1, window[1].1 .0,
                "a gap or an overlap between two cells"
            );
        }
    }

    #[test]
    fn the_key_is_a_ramp_to_the_strongest_pair_or_the_ceiling() {
        let pairs = vec![Pair::new(100, 400, 0.4)];
        let track = PairTrack::new(pairs.clone()).label("r²");
        let theme = Theme::light();
        let key = track.key(&region(), 1.0, &theme).unwrap();
        assert!(format!("{key:?}").contains("0.4"));
        let pinned = PairTrack::new(pairs).label("r²").ceiling(1.0);
        assert!(format!("{:?}", pinned.key(&region(), 1.0, &theme).unwrap()).contains("\"1\""));
    }
}
