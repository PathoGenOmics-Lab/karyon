//! How a pangenome grows as genomes are added to it.
//!
//! # The idea
//!
//! Sequence one genome and you have one gene repertoire. Sequence a second and
//! some of its genes are new. The question a pangenome paper opens with is
//! whether that ever stops: does the total repertoire flatten off, in which
//! case the species has a pangenome you can finish counting, or does it keep
//! climbing, in which case every new isolate still brings genes nobody has seen.
//! The core, the genes every genome has, falls at the same time and for the
//! same reason.
//!
//! Both curves depend on the order the genomes went in, which is an accident of
//! how they were listed. So the order is shuffled many times over and the
//! spread across those shuffles is drawn along with the middle of them. A curve
//! without that spread is one arbitrary ordering presented as a result.
//!
//! # What the x axis is
//!
//! The number of genomes, not a position. Put the figure over
//! `Region::new("genomes", 0, n)` and an
//! [`AxisTrack`](crate::AxisTrack) with
//! [`center_on_bases`](crate::AxisTrack::center_on_bases) will number it
//! correctly, since a count of `k` genomes sits in the `k`th column.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Track};

/// One curve, and every observation behind each of its points.
#[derive(Debug, Clone, PartialEq)]
pub struct AccumulationCurve {
    /// What the curve counts, drawn in the legend.
    pub name: String,
    /// `samples[k - 1]` holds the counts seen at `k` genomes, one per shuffle.
    pub samples: Vec<Vec<f64>>,
}

impl AccumulationCurve {
    /// A curve named `name`.
    pub fn new(name: impl Into<String>, samples: Vec<Vec<f64>>) -> Self {
        AccumulationCurve {
            name: name.into(),
            samples,
        }
    }

    /// How many genome counts the curve covers.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Whether the curve has no points.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// The quantile of the observations at `k` genomes, `k` counting from one.
    ///
    /// `quantile(k, 0.5)` is the median, which is what the line is drawn
    /// through. `None` when nothing was observed at that count.
    pub fn quantile(&self, genomes: usize, q: f64) -> Option<f64> {
        let observations = self.samples.get(genomes.checked_sub(1)?)?;
        let mut sorted: Vec<f64> = observations
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .collect();
        if sorted.is_empty() {
            return None;
        }
        sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN survives the filter"));
        let position = q.clamp(0.0, 1.0) * (sorted.len() - 1) as f64;
        let below = position.floor() as usize;
        let above = position.ceil() as usize;
        let weight = position - below as f64;
        Some(sorted[below] * (1.0 - weight) + sorted[above] * weight)
    }
}

/// Pangenome and core accumulation curves.
///
/// ```
/// use karyon::{AccumulationTrack, AxisTrack, Figure, Region};
///
/// // Four genomes: three share two genes, the fourth brings one of its own.
/// let genomes = vec![
///     vec![true, true, false],
///     vec![true, true, false],
///     vec![true, true, false],
///     vec![true, true, true],
/// ];
///
/// let track = AccumulationTrack::from_presence(&genomes, 20, 7);
/// let svg = Figure::new(Region::new("genomes", 0, 4).unwrap())
///     .push(track.label("gene families"))
///     .push(AxisTrack::new().center_on_bases(true))
///     .to_svg();
/// assert!(svg.contains("pangenome"));
/// ```
#[derive(Debug, Clone)]
pub struct AccumulationTrack {
    curves: Vec<AccumulationCurve>,
    label: Option<String>,
    height: f64,
    colors: Vec<String>,
    band: (f64, f64),
    show_band: bool,
    show_points: bool,
    show_legend: bool,
    show_scale: bool,
    max: Option<f64>,
    floor: Option<f64>,
}

impl AccumulationTrack {
    /// A track over curves you have already computed.
    pub fn new(curves: impl Into<Vec<AccumulationCurve>>) -> Self {
        AccumulationTrack {
            curves: curves.into(),
            label: None,
            height: 150.0,
            colors: Vec::new(),
            band: (0.05, 0.95),
            show_band: true,
            show_points: true,
            show_legend: true,
            show_scale: true,
            max: None,
            floor: None,
        }
    }

    /// The pangenome and core curves of a presence and absence matrix.
    ///
    /// One row per genome, one column per gene family, `true` where the genome
    /// carries it. The genome order is shuffled `permutations` times from
    /// `seed`, and every shuffle contributes one observation at each count, so
    /// the spread drawn around the curves is the spread the ordering causes.
    ///
    /// A hundred permutations is the usual number and is cheap here: the whole
    /// thing is one pass per genome per shuffle.
    ///
    /// A presence matrix held as [`MatrixRow`](crate::MatrixRow)s, which is
    /// what the pangenome matrix plot takes, converts in a line:
    ///
    /// ```
    /// use karyon::{AccumulationTrack, MatrixRow};
    ///
    /// let rows = vec![MatrixRow::new("ERR001", vec![1.0, 0.0, 1.0])];
    /// let genomes: Vec<Vec<bool>> = rows
    ///     .iter()
    ///     .map(|row| row.values.iter().map(|v| *v > 0.5).collect())
    ///     .collect();
    /// assert_eq!(AccumulationTrack::from_presence(&genomes, 10, 1).curves().len(), 2);
    /// ```
    pub fn from_presence(genomes: &[Vec<bool>], permutations: usize, seed: u64) -> Self {
        let count = genomes.len();
        let families = genomes.iter().map(|row| row.len()).max().unwrap_or(0);
        if count == 0 || families == 0 {
            return AccumulationTrack::new(Vec::new());
        }

        let mut pan: Vec<Vec<f64>> = vec![Vec::new(); count];
        let mut core: Vec<Vec<f64>> = vec![Vec::new(); count];
        let mut order: Vec<usize> = (0..count).collect();
        let mut rng = Lcg::new(seed);

        for permutation in 0..permutations.max(1) {
            // The first ordering is the one the caller handed in, so a single
            // permutation reproduces their own list rather than a shuffle of it.
            if permutation > 0 {
                for index in (1..count).rev() {
                    order.swap(index, (rng.next() % (index as u64 + 1)) as usize);
                }
            }
            let mut seen = vec![0usize; families];
            let mut present = 0usize;
            for (added, genome) in order.iter().enumerate() {
                let row = &genomes[*genome];
                for (family, count) in seen.iter_mut().enumerate() {
                    // A row shorter than the widest one is missing those
                    // families, not carrying them.
                    if row.get(family).copied().unwrap_or(false) {
                        if *count == 0 {
                            present += 1;
                        }
                        *count += 1;
                    }
                }
                let shared = seen.iter().filter(|count| **count == added + 1).count();
                pan[added].push(present as f64);
                core[added].push(shared as f64);
            }
        }

        AccumulationTrack::new(vec![
            AccumulationCurve::new("pangenome", pan),
            AccumulationCurve::new("core", core),
        ])
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(30.0);
        self
    }

    /// Overrides the colours, one per curve, cycled if there are too few.
    pub fn colors(mut self, colors: impl Into<Vec<String>>) -> Self {
        self.colors = colors.into();
        self
    }

    /// Sets the quantiles the shaded band runs between.
    pub fn band(mut self, low: f64, high: f64) -> Self {
        self.band = (low.clamp(0.0, 1.0), high.clamp(0.0, 1.0));
        self
    }

    /// Draws or hides the shaded spread across permutations.
    ///
    /// Leave it on. Without it the figure is one arbitrary genome ordering
    /// drawn as though it were the answer.
    pub fn show_band(mut self, show: bool) -> Self {
        self.show_band = show;
        self
    }

    /// Draws or hides a marker at each genome count.
    pub fn show_points(mut self, show: bool) -> Self {
        self.show_points = show;
        self
    }

    /// Draws or hides the legend naming the curves.
    pub fn show_legend(mut self, show: bool) -> Self {
        self.show_legend = show;
        self
    }

    /// Draws or hides the count axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Pins the top of the count axis.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self
    }

    /// Pins the bottom of the count axis.
    ///
    /// `floor(0.0)` anchors it at zero, which is the honest default for a bar
    /// chart and the wrong one here: see [`AccumulationTrack::range`].
    pub fn floor(mut self, floor: f64) -> Self {
        self.floor = Some(floor);
        self
    }

    /// The curves.
    pub fn curves(&self) -> &[AccumulationCurve] {
        &self.curves
    }

    /// How many genomes the widest curve covers.
    pub fn genomes(&self) -> usize {
        self.curves
            .iter()
            .map(|curve| curve.len())
            .max()
            .unwrap_or(0)
    }

    /// The counts the axis spans, low first.
    ///
    /// Fitted to the curves rather than anchored at zero. A pangenome and its
    /// core differ from each other by far less than either differs from
    /// nothing, so a zero-anchored axis squashes both into the top of the band
    /// and hides the one thing the figure is about, which is the shape of the
    /// climb and the fall. Both ends are labelled, so the axis is not being
    /// quiet about what it did, and [`AccumulationTrack::floor`] puts it back
    /// on zero for anyone who wants it there.
    pub fn range(&self) -> (f64, f64) {
        let (low, high) = self.band;
        let extremes = |q: f64| {
            self.curves
                .iter()
                .flat_map(move |curve| (1..=curve.len()).filter_map(move |k| curve.quantile(k, q)))
        };
        let bottom = extremes(low.min(0.5)).fold(f64::INFINITY, f64::min);
        let top = extremes(high.max(0.5)).fold(f64::NEG_INFINITY, f64::max);
        if !bottom.is_finite() || !top.is_finite() {
            return (self.floor.unwrap_or(0.0), self.max.unwrap_or(1.0));
        }
        let pad = ((top - bottom) * 0.08).max(1.0);
        (
            self.floor.unwrap_or((bottom - pad).max(0.0)),
            self.max.unwrap_or(top + pad),
        )
    }

    /// Top of the count axis.
    pub fn ceiling(&self) -> f64 {
        self.range().1
    }

    fn color(&self, index: usize, theme: &Theme) -> String {
        if self.colors.is_empty() {
            return theme.color(index).to_string();
        }
        self.colors[index % self.colors.len()].clone()
    }
}

impl Track for AccumulationTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.curves.is_empty() {
            return 0.0;
        }
        let (bottom, top) = self.range();
        let widest = [bottom, top]
            .iter()
            .map(|value| text_width(&group_thousands(*value), theme.font_size - 1.0))
            .fold(0.0f64, f64::max);
        widest + 8.0
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

        if self.curves.is_empty() || self.genomes() == 0 {
            return;
        }
        let (bottom, top) = self.range();
        // Leave room for the legend, or the pangenome curve runs straight
        // through it: the axis is fitted to the curves, so the top one comes
        // within a few per cent of the band's ceiling by construction.
        let legend_room = if self.show_legend {
            ctx.theme.font_size + 4.0
        } else {
            0.0
        };
        let plot_height = (band.h - legend_room).max(4.0);
        let y_of = |value: f64| {
            baseline - ((value - bottom) / (top - bottom)).clamp(0.0, 1.0) * plot_height
        };
        // A count of k genomes belongs in the kth column, which is index k - 1.
        let x_of = |genomes: usize| ctx.scale.x_center(genomes as u64 - 1);
        let (low, high) = self.band;

        for (index, curve) in self.curves.iter().enumerate() {
            let color = self.color(index, ctx.theme);

            if self.show_band && high > low {
                // Up one edge and back down the other, so the two quantiles
                // close into a shape rather than being two lines to compare.
                let mut top: Vec<(f64, f64)> = Vec::new();
                let mut bottom: Vec<(f64, f64)> = Vec::new();
                for genomes in 1..=curve.len() {
                    if let (Some(hi), Some(lo)) =
                        (curve.quantile(genomes, high), curve.quantile(genomes, low))
                    {
                        top.push((x_of(genomes), y_of(hi)));
                        bottom.push((x_of(genomes), y_of(lo)));
                    }
                }
                if top.len() >= 2 {
                    let mut d = String::with_capacity(top.len() * 24);
                    for (step, (x, y)) in top.iter().chain(bottom.iter().rev()).enumerate() {
                        d.push_str(if step == 0 { "M" } else { " L" });
                        d.push_str(&num(*x));
                        d.push(' ');
                        d.push_str(&num(*y));
                    }
                    d.push_str(" Z");
                    ctx.svg.path(&d, &color, 0.18);
                }
            }

            let middle: Vec<(f64, f64)> = (1..=curve.len())
                .filter_map(|genomes| {
                    curve
                        .quantile(genomes, 0.5)
                        .map(|value| (x_of(genomes), y_of(value)))
                })
                .collect();
            ctx.svg.polyline(&middle, &color, 2.0);

            if self.show_points {
                // A marker only where there is room for one, so a hundred
                // genomes stay a curve rather than becoming a string of beads.
                let spacing = band.w / curve.len().max(1) as f64;
                if spacing >= 9.0 {
                    for (x, y) in &middle {
                        ctx.svg
                            .circle_ringed(*x, *y, 2.4, &color, ctx.theme.surface(), 1.2);
                    }
                }
            }
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            ctx.svg.text(
                right,
                y_of(top) + size * 0.35,
                &group_thousands(top),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
            // Both ends are named, because an axis that does not start at zero
            // and does not say where it does start is a misleading axis.
            ctx.svg.text(
                right,
                baseline - size * 0.22,
                &group_thousands(bottom),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }

        if self.show_legend {
            let size = ctx.theme.font_size - 1.0;
            let mut x = band.right();
            for (index, curve) in self.curves.iter().enumerate().rev() {
                let color = self.color(index, ctx.theme);
                x -= text_width(&curve.name, size);
                ctx.svg.text(
                    x,
                    band.y + size,
                    &curve.name,
                    &ctx.theme.muted,
                    size,
                    Anchor::Start,
                );
                x -= 6.0;
                ctx.svg.circle(x, band.y + size * 0.65, 3.0, &color);
                x -= 14.0;
            }
        }
    }
}

/// `12345` as `12,345`, rounded to a whole count.
fn group_thousands(value: f64) -> String {
    let digits = format!("{}", value.max(0.0).round() as u64);
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, c) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// A linear congruential generator, so a figure is reproducible from its seed
/// and the crate stays free of dependencies.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    /// Four genomes. Two families in all of them, one family each in the last
    /// two, so the core is two and the pangenome is four.
    fn genomes() -> Vec<Vec<bool>> {
        vec![
            vec![true, true, false, false],
            vec![true, true, false, false],
            vec![true, true, true, false],
            vec![true, true, false, true],
        ]
    }

    fn region() -> Region {
        Region::new("genomes", 0, 4).unwrap()
    }

    #[test]
    fn the_pangenome_climbs_and_the_core_falls() {
        let track = AccumulationTrack::from_presence(&genomes(), 50, 11);
        let pan = &track.curves()[0];
        let core = &track.curves()[1];
        assert_eq!(pan.name, "pangenome");
        assert_eq!(core.name, "core");

        // Every ordering ends at the same place: all four families seen, two of
        // them in every genome.
        assert_eq!(pan.quantile(4, 0.5), Some(4.0));
        assert_eq!(core.quantile(4, 0.5), Some(2.0));

        // Neither curve ever goes the wrong way.
        for genomes in 2..=4 {
            assert!(pan.quantile(genomes, 0.5) >= pan.quantile(genomes - 1, 0.5));
            assert!(core.quantile(genomes, 0.5) <= core.quantile(genomes - 1, 0.5));
        }
    }

    #[test]
    fn one_genome_is_its_own_core_and_its_own_pangenome() {
        let track = AccumulationTrack::from_presence(&genomes(), 30, 3);
        let pan = track.curves()[0].quantile(1, 0.5).unwrap();
        let core = track.curves()[1].quantile(1, 0.5).unwrap();
        assert_eq!(pan, core, "with one genome there is nothing else to share");
    }

    #[test]
    fn the_ordering_is_what_the_spread_measures() {
        let track = AccumulationTrack::from_presence(&genomes(), 200, 99);
        let pan = &track.curves()[0];
        // At two genomes the answer depends on which two: two of the first pair
        // gives three families between them, the last two give four.
        let lo = pan.quantile(2, 0.0).unwrap();
        let hi = pan.quantile(2, 1.0).unwrap();
        assert!(hi > lo, "every ordering gave the same answer: {lo} to {hi}");
    }

    #[test]
    fn the_same_seed_gives_the_same_figure() {
        let one = AccumulationTrack::from_presence(&genomes(), 40, 2024);
        let two = AccumulationTrack::from_presence(&genomes(), 40, 2024);
        assert_eq!(one.curves(), two.curves());
    }

    #[test]
    fn the_first_ordering_is_the_one_that_was_handed_in() {
        // A single permutation reproduces the caller's own genome list rather
        // than a shuffle of it, so the curve can be checked by hand.
        let track = AccumulationTrack::from_presence(&genomes(), 1, 0);
        let pan: Vec<f64> = (1..=4)
            .map(|k| track.curves()[0].quantile(k, 0.5).unwrap())
            .collect();
        assert_eq!(pan, vec![2.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn a_quantile_interpolates_between_observations() {
        let curve = AccumulationCurve::new("c", vec![vec![10.0, 20.0, 30.0]]);
        assert_eq!(curve.quantile(1, 0.0), Some(10.0));
        assert_eq!(curve.quantile(1, 0.5), Some(20.0));
        assert_eq!(curve.quantile(1, 1.0), Some(30.0));
        assert_eq!(curve.quantile(1, 0.25), Some(15.0));
        assert_eq!(curve.quantile(2, 0.5), None, "past the end of the curve");
        assert_eq!(curve.quantile(0, 0.5), None, "zero genomes is no curve");
    }

    #[test]
    fn a_curve_of_nothing_but_missing_values_has_no_quantile() {
        let curve = AccumulationCurve::new("c", vec![vec![f64::NAN]]);
        assert_eq!(curve.quantile(1, 0.5), None);
    }

    #[test]
    fn an_empty_matrix_gives_an_empty_track_rather_than_a_panic() {
        assert!(AccumulationTrack::from_presence(&[], 10, 1)
            .curves()
            .is_empty());
        assert!(AccumulationTrack::from_presence(&[vec![], vec![]], 10, 1)
            .curves()
            .is_empty());
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(AccumulationTrack::new(Vec::new()).label("none"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn ragged_rows_are_read_as_absences_rather_than_dropped() {
        let ragged = vec![vec![true, true, true], vec![true]];
        let track = AccumulationTrack::from_presence(&ragged, 5, 1);
        assert_eq!(track.curves()[0].quantile(2, 0.5), Some(3.0));
        assert_eq!(track.curves()[1].quantile(2, 0.5), Some(1.0));
    }

    #[test]
    fn both_curves_and_their_spread_are_drawn() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&genomes(), 60, 5).label("families"))
            .to_svg();
        assert_eq!(svg.matches("<polyline").count(), 2, "one line per curve");
        assert!(svg.contains("<path"), "the spread is a shape");
        assert!(svg.contains(">pangenome</text>"));
        assert!(svg.contains(">core</text>"));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_spread_can_be_turned_off_but_the_lines_stay() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&genomes(), 60, 5).show_band(false))
            .to_svg();
        assert!(!svg.contains("<path"));
        assert_eq!(svg.matches("<polyline").count(), 2);
    }

    #[test]
    fn markers_go_away_once_the_genomes_are_packed_too_tight() {
        let many: Vec<Vec<bool>> = (0..400).map(|i| vec![true, i % 2 == 0]).collect();
        let region = Region::new("genomes", 0, 400).unwrap();
        let crowded = Figure::new(region)
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&many, 5, 1).show_legend(false))
            .to_svg();
        assert!(
            !crowded.contains("<circle"),
            "four hundred beads on a string"
        );

        // Few enough genomes and every count is worth marking.
        let few: Vec<Vec<bool>> = (0..6).map(|i| vec![true, i % 2 == 0]).collect();
        let roomy = Figure::new(Region::new("genomes", 0, 6).unwrap())
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&few, 5, 1).show_legend(false))
            .to_svg();
        assert!(roomy.contains("<circle"));
    }

    #[test]
    fn the_axis_is_only_reserved_when_it_is_wanted() {
        let theme = Theme::light();
        let track = AccumulationTrack::from_presence(&genomes(), 10, 1);
        assert!(track.y_axis_width(&theme) > 0.0);
        assert_eq!(track.clone().show_scale(false).y_axis_width(&theme), 0.0);
        assert_eq!(AccumulationTrack::new(Vec::new()).y_axis_width(&theme), 0.0);
    }

    #[test]
    fn the_axis_fits_the_curves_rather_than_starting_at_zero() {
        // A core of three thousand and a pangenome of four is a figure about
        // the gap between them. Anchored at zero both curves land in the top
        // fifth of the band and the shape of neither can be read.
        let big: Vec<Vec<bool>> = (0..6)
            .map(|genome| {
                (0..4_000)
                    .map(|family| family < 3_000 || family % 6 == genome)
                    .collect()
            })
            .collect();
        let track = AccumulationTrack::from_presence(&big, 20, 1);
        let (bottom, top) = track.range();
        assert!(bottom > 2_000.0, "the empty bottom is not drawn: {bottom}");
        assert!(top > 3_000.0 && top < 4_500.0);

        // And it can be put back on zero for anyone who wants it there.
        assert_eq!(track.clone().floor(0.0).range().0, 0.0);
    }

    #[test]
    fn both_ends_of_the_count_axis_are_named() {
        // An axis that does not start at zero and does not say where it starts
        // is a misleading axis.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&genomes(), 20, 1).show_legend(false))
            .to_svg();
        let (bottom, top) = AccumulationTrack::from_presence(&genomes(), 20, 1).range();
        assert!(svg.contains(&format!(">{}</text>", group_thousands(bottom))));
        assert!(svg.contains(&format!(">{}</text>", group_thousands(top))));
    }

    #[test]
    fn the_legend_gets_room_rather_than_sitting_on_the_curves() {
        // The axis is fitted to the curves, so the top one reaches within a few
        // per cent of the band's ceiling. A legend written there is a legend
        // with a line through it.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(AccumulationTrack::from_presence(&genomes(), 20, 1).height(120.0))
            .to_svg();

        let legend_baseline = {
            let at = svg.find(">pangenome</text>").unwrap();
            let head = &svg[..at];
            let y = head.rfind("y=\"").unwrap() + 3;
            head[y..].split('"').next().unwrap().parse::<f64>().unwrap()
        };
        // The curves themselves, which the legend swatches are not.
        let highest_point = svg
            .match_indices("<polyline points=\"")
            .flat_map(|(at, key)| {
                let rest = &svg[at + key.len()..];
                let list = &rest[..rest.find('"').unwrap()];
                list.split(' ')
                    .filter_map(|pair| pair.split_once(','))
                    .filter_map(|(_, y)| y.parse::<f64>().ok())
                    .collect::<Vec<f64>>()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(
            highest_point > legend_baseline,
            "a curve reaches {highest_point}, above the legend baseline {legend_baseline}"
        );
    }

    #[test]
    fn thousands_are_grouped_in_the_axis_label() {
        assert_eq!(group_thousands(0.0), "0");
        assert_eq!(group_thousands(999.4), "999");
        assert_eq!(group_thousands(12_345.0), "12,345");
    }
}
