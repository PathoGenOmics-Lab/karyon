//! Association statistics along the sequence: the Manhattan plot.
//!
//! One point per test, height by significance, and a line where significance
//! starts. The plot is named for the skyline that appears when a real signal
//! stacks a run of neighbouring markers into a tower, and it is read as a
//! texture rather than as a set of markers.
//!
//! # A tower is only a tower against a crowd
//!
//! Which settles the drawing. Points are small by default, and a point above
//! the threshold is given a ring rather than a larger disc: the ring keeps it
//! findable where the crowd is densest without swallowing the neighbours it has
//! to be seen against.
//!
//! # What a pixel holds
//!
//! A scan is usually far denser than the page: a hundred thousand tests over a
//! megabase put more than a hundred on every pixel column. A point drawn over
//! one of its own shape and colour on the same pixel is the same ink twice, so
//! each pixel keeps one point of each look, and the document grows with the
//! figure rather than with the scan. The point kept is the last one drawn
//! there, the one that was on top, and every point left out is less than a
//! pixel away, in each direction, from one of its own look that is drawn.
//!
//! # Where significance starts is not decided here
//!
//! [`ManhattanTrack::threshold`] takes a number from the caller,
//! [`ManhattanTrack::significant`] returns nothing at all until it has one, and
//! no point is coloured as a hit without it. What is being plotted is the
//! caller's too: [`Association::from_p_value`] converts a p-value into the
//! conventional `-log10(p)`, but a test statistic or a Bayes factor plots the
//! same way and only [`ManhattanTrack::unit`] tells the axis which of them it
//! is showing.

use std::collections::BTreeSet;

use crate::scale::Scale;
use crate::style::{legible_ticks, Emphasis, LinePattern, QuantitativeAxis, Symbol};
use crate::svg::Anchor;
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// One tested position and what the test said about it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Association {
    /// Position, 0-based.
    pub pos: u64,
    /// The plotted statistic, conventionally `-log10(p)`.
    pub value: f64,
}

impl Association {
    /// A point at a 0-based position.
    pub fn new(pos: u64, value: f64) -> Self {
        Association { pos, value }
    }

    /// A point from a p-value, converted to `-log10(p)`.
    ///
    /// A p-value of zero has no logarithm, so it is floored at the smallest
    /// positive double rather than becoming infinity. Anything at that floor is
    /// already far above any threshold you would draw.
    pub fn from_p_value(pos: u64, p: f64) -> Self {
        let p = if p.is_finite() {
            p.max(f64::MIN_POSITIVE)
        } else {
            1.0
        };
        Association {
            pos,
            value: -p.log10(),
        }
    }
}

/// Points scored by significance, with a threshold line.
///
/// The x axis is genomic, so this draws one chromosome or one region of one.
/// A genome-wide plot that lays every chromosome end to end is a different
/// coordinate system, and this crate does not pretend otherwise: give it a
/// region spanning the sequence you are testing.
///
/// Each pixel keeps one point of each look however many tests land on it, so
/// the SVG grows with the figure rather than with the scan (see
/// [the module](crate::track::manhattan)).
///
/// ```
/// use karyon::{Association, Figure, ManhattanTrack, Region};
///
/// let hits = vec![
///     Association::from_p_value(1_000, 0.4),
///     Association::from_p_value(2_000, 3e-9),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 0, 3_000).unwrap())
///     .push(ManhattanTrack::new(hits).genome_wide_threshold().label("GWAS"))
///     .to_svg();
/// assert!(svg.contains("<circle"));
/// ```
#[derive(Debug, Clone)]
pub struct ManhattanTrack {
    points: Vec<Association>,
    label: Option<String>,
    height: f64,
    radius: f64,
    max: Option<f64>,
    threshold: Option<f64>,
    threshold_label: Option<String>,
    color: Option<String>,
    significant_color: Option<String>,
    unit: String,
    title: Option<String>,
    show_scale: bool,
    bands: Vec<u64>,
    axis: QuantitativeAxis,
}

impl ManhattanTrack {
    /// A track holding `points`.
    pub fn new(points: impl Into<Vec<Association>>) -> Self {
        ManhattanTrack {
            points: points.into(),
            label: None,
            height: 90.0,
            radius: 2.2,
            max: None,
            threshold: None,
            threshold_label: None,
            color: None,
            significant_color: None,
            unit: String::new(),
            title: None,
            show_scale: true,
            bands: Vec::new(),
            axis: QuantitativeAxis::new(),
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(8.0);
        self
    }

    /// Alternates the point colour at each of these positions.
    ///
    /// For a figure whose axis is several sequences laid end to end: hand it
    /// [`Genome::boundaries`](crate::Genome::boundaries) and the chromosomes
    /// come out in alternating shades, which is the only thing telling a
    /// reader of a genome-wide plot where one ends and the next begins. A
    /// point past the last boundary keeps the last shade.
    pub fn bands(mut self, boundaries: impl Into<Vec<u64>>) -> Self {
        self.bands = boundaries.into();
        self.bands.sort_unstable();
        self
    }

    /// Which band a position falls in, counting from zero.
    pub fn band_of(&self, position: u64) -> usize {
        self.bands
            .iter()
            .rposition(|start| *start <= position)
            .unwrap_or(0)
    }

    /// Sets the radius of a point.
    ///
    /// Points are small on purpose: a Manhattan plot is read as a texture with
    /// towers in it, and a marker big enough to admire individually turns the
    /// texture into a smear.
    pub fn radius(mut self, radius: f64) -> Self {
        self.radius = radius.max(0.3);
        self
    }

    /// Pins the top of the y axis.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self.axis.max = Some(max);
        self
    }

    /// Uses a shared quantitative-axis contract for range, ticks, units and
    /// reference lines.
    pub fn axis(mut self, axis: QuantitativeAxis) -> Self {
        self.max = axis.max;
        self.unit = axis.unit.clone();
        self.axis = axis;
        self
    }

    /// Draws a horizontal line at `value` and colours what rises above it.
    ///
    /// The line is labelled with its value, as the axis writes its numbers;
    /// [`ManhattanTrack::threshold_label`] says something else.
    pub fn threshold(mut self, value: f64) -> Self {
        self.threshold = Some(value);
        self.threshold_label = None;
        self
    }

    /// Draws the line at the p-value `p`, on a track whose points are
    /// `-log10 p`, and labels it with the p-value, as `p = 1e-5`.
    ///
    /// A `p` outside nought to one is no p-value, and draws no line.
    pub fn p_value_threshold(mut self, p: f64) -> Self {
        if p > 0.0 && p <= 1.0 {
            self = self.threshold(-p.log10());
            self.threshold_label = Some(format!("p = {}", p_text(p)));
        }
        self
    }

    /// The usual genome-wide threshold, p = 5e-8, drawn at `-log10(5e-8)`,
    /// about 7.3.
    ///
    /// It is a Bonferroni correction for a million independent tests, which is
    /// the convention in human GWAS and frequently the wrong number everywhere
    /// else: what it should be follows from how many independent tests were
    /// really run, and a shorter genome or stronger linkage leaves far fewer
    /// than a million. Set your own if you know it.
    pub fn genome_wide_threshold(self) -> Self {
        self.p_value_threshold(5e-8)
    }

    /// Sets the words written on the threshold line, in place of its value;
    /// an empty string writes none.
    pub fn threshold_label(mut self, label: impl Into<String>) -> Self {
        self.threshold_label = Some(label.into());
        self
    }

    /// What the value axis measures, written under the track's name, such as
    /// `-log10 p`.
    pub fn axis_title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Sets the colour of points below the threshold.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets the colour of points above it.
    pub fn significant_color(mut self, color: impl Into<String>) -> Self {
        self.significant_color = Some(color.into());
        self
    }

    /// Sets a unit suffix for the top tick's label.
    ///
    /// What the axis measures, such as `-log10 p`, is a title rather than a
    /// unit, and reads as one under the track's name:
    /// [`ManhattanTrack::axis_title`].
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self.axis.unit = self.unit.clone();
        self
    }

    /// Draws or hides the value axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// The points in the track.
    pub fn points(&self) -> &[Association] {
        &self.points
    }

    /// Points at or above the threshold, in input order.
    ///
    /// Empty when no threshold was set, since without one no point is
    /// significant and calling some of them so would be an invention.
    pub fn significant(&self) -> Vec<Association> {
        let Some(threshold) = self.threshold else {
            return Vec::new();
        };
        self.points
            .iter()
            .filter(|p| p.value.is_finite() && p.value >= threshold)
            .copied()
            .collect()
    }

    /// Top of the axis: the pinned maximum, or the tallest point with a little
    /// headroom, or the threshold when every point is under it.
    fn ceiling(&self) -> f64 {
        if let Some(max) = self.axis.max.or(self.max) {
            return max.max(1e-9);
        }
        let tallest = self
            .points
            .iter()
            .map(|p| p.value)
            .filter(|v| v.is_finite())
            .fold(0.0f64, f64::max);
        let with_threshold = self.threshold.map_or(tallest, |t| tallest.max(t));
        (with_threshold * 1.08).max(1e-9)
    }

    /// The floor and ceiling points are plotted against, with a free ceiling
    /// rounded up to a value worth labelling.
    ///
    /// The rounding is the headroom, so it starts from the tallest point and
    /// the threshold themselves rather than from the ceiling [`Self::ceiling`]
    /// has already padded: padding and then rounding stacked the two, and an
    /// association peaking at ten was drawn against fifteen.
    fn range(&self) -> (f64, f64) {
        if self.axis.max.or(self.max).is_some() {
            return self.axis.resolve(0.0, self.ceiling());
        }
        let (floor, ceiling) = self.axis.resolve(0.0, self.tallest().max(1e-9));
        self.axis.nice(floor, ceiling)
    }

    /// The tallest point, or the threshold if that is higher.
    fn tallest(&self) -> f64 {
        let tallest = self
            .points
            .iter()
            .map(|p| p.value)
            .filter(|v| v.is_finite())
            .fold(0.0f64, f64::max);
        self.threshold.map_or(tallest, |t| tallest.max(t))
    }
}

impl Track for ManhattanTrack {
    fn noun(&self) -> &str {
        "an association scan"
    }

    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn axis_title(&self) -> Option<&str> {
        (self.show_scale && !self.points.is_empty())
            .then_some(self.title.as_deref())
            .flatten()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.points.is_empty() {
            return 0.0;
        }
        let (floor, ceiling) = self.range();
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

        if self.points.is_empty() {
            return;
        }
        let (floor, ceiling) = self.range();
        let y_of =
            |value: f64| baseline - ((value - floor) / (ceiling - floor)).clamp(0.0, 1.0) * band.h;
        let size = ctx.theme.font_size - 1.0;
        let ticks = self.axis.values(floor, ceiling);
        let shown = legible_ticks(&ticks, y_of, size);

        for &value in shown.iter().filter(|&&value| value > floor) {
            ctx.svg.line(
                band.x,
                y_of(value),
                band.right(),
                y_of(value),
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
        }

        let plain = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.muted.clone());
        let significant = self
            .significant_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());

        if let Some(threshold) = self.threshold {
            if threshold >= floor && threshold <= ceiling {
                let y = y_of(threshold);
                ctx.svg.line_pattern(
                    band.x,
                    y,
                    band.right(),
                    y,
                    &significant,
                    ctx.theme.tokens.stroke,
                    LinePattern::Dashed,
                );
            }
        }

        for reference in &self.axis.references {
            if !reference.value.is_finite() || reference.value < floor || reference.value > ceiling
            {
                continue;
            }
            let style = ctx.theme.mark_style(reference.emphasis);
            let color = if reference.emphasis == Emphasis::Alert {
                &significant
            } else {
                &ctx.theme.muted
            };
            ctx.svg.line_pattern(
                band.x,
                y_of(reference.value),
                band.right(),
                y_of(reference.value),
                color,
                style.stroke_width,
                reference.pattern,
            );
        }

        // A point drawn over one of its own look on the same pixel is the same
        // ink twice, which adds an element to the document and nothing to the
        // picture, and a scan is usually far denser than its pixels. So each
        // pixel keeps one point of each look, and the one kept is the last,
        // since the last was on top: keeping the first would bring a hit in
        // front of a miss that was drawn over it. Measured on a hundred
        // thousand tests over a megabase, 101,385 elements and 5.9 MB became
        // 10,996 elements and 0.6 MB.
        let mut taken: BTreeSet<(i64, i64, bool, u8, bool)> = BTreeSet::new();
        let mut kept: Vec<(f64, f64, bool, Symbol, bool)> = Vec::new();
        for point in self.points.iter().rev() {
            if !ctx.region.contains(point.pos) || !point.value.is_finite() {
                continue;
            }
            let (x, y) = (ctx.scale.x_center(point.pos), y_of(point.value));
            let above = self.threshold.is_some_and(|t| point.value >= t);
            // A hit looks the same in every band. A miss takes the shape of its
            // band, and every other sequence a shade lighter, which is what
            // separates one chromosome from the next when the axis is all of
            // them.
            let (symbol, lighter) = if above {
                (Symbol::Diamond, false)
            } else {
                let nth = self.band_of(point.pos);
                (ctx.theme.symbol(nth), self.bands.len() > 1 && nth % 2 == 1)
            };
            let pixel = (x.round() as i64, y.round() as i64);
            if taken.insert((pixel.0, pixel.1, above, symbol as u8, lighter)) {
                kept.push((x, y, above, symbol, lighter));
            }
        }

        let shaded = mix(&plain, ctx.theme.surface(), 0.42);
        let radius = self.radius * ctx.visual_scale;
        for &(x, y, above, symbol, lighter) in kept.iter().rev() {
            if above {
                // A hit is worth a ring, so it stays a point where the texture
                // around it is densest.
                ctx.svg.symbol_ringed(
                    x,
                    y,
                    radius + ctx.theme.tokens.hairline,
                    symbol,
                    &significant,
                    ctx.theme.surface(),
                    ctx.theme.tokens.hairline,
                );
            } else {
                let color = if lighter { &shaded } else { &plain };
                ctx.svg.symbol(x, y, radius, symbol, color);
            }
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let right = ctx.axis.right() - 4.0;
            for (value, label) in ticks.iter().zip(self.axis.labels(&ticks)) {
                if !shown.contains(value) {
                    continue;
                }
                ctx.svg.text(
                    right,
                    (y_of(*value) + size * 0.35)
                        .max(band.y + size)
                        .min(baseline - size * 0.22),
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

        // Last, so no point is drawn over the words. A dashed line with
        // nothing on it was a line the reader had to be told the meaning of.
        if let (Some(threshold), Some(text)) = (self.threshold, self.threshold_text()) {
            if threshold >= floor && threshold <= ceiling {
                // Over the line where a line of text fits there, and under it
                // where the line runs along the top of the band, rather than
                // pushed down onto the line and struck through by it.
                let y = y_of(threshold);
                let gap = ctx.theme.tokens.row_gap;
                let baseline = if y - gap - size * 0.8 >= band.y {
                    y - gap
                } else {
                    y + gap + size * 0.8
                };
                ctx.svg.text(
                    band.x + ctx.theme.tokens.label_gap,
                    baseline,
                    &text,
                    &significant,
                    size,
                    Anchor::Start,
                );
            }
        }
    }
}

impl ManhattanTrack {
    /// The words on the threshold line: what was set, or the value as the
    /// axis writes its numbers, and none for an empty label.
    fn threshold_text(&self) -> Option<String> {
        let value = self.threshold?;
        match &self.threshold_label {
            Some(label) if label.is_empty() => None,
            Some(label) => Some(label.clone()),
            None => Some(self.axis.format.format(value)),
        }
    }
}

/// A p-value as it is usually written: `0.05` and `0.001` as they are, and
/// smaller ones in scientific notation, `5e-8`.
fn p_text(p: f64) -> String {
    if p >= 1e-3 {
        format!("{p}")
    } else {
        format!("{p:e}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 0, 10_000).unwrap()
    }

    #[test]
    fn a_p_value_becomes_minus_log10_of_itself() {
        assert!((Association::from_p_value(0, 0.1).value - 1.0).abs() < 1e-12);
        assert!((Association::from_p_value(0, 1e-8).value - 8.0).abs() < 1e-12);
        assert_eq!(Association::from_p_value(0, 1.0).value, 0.0);
    }

    #[test]
    fn a_p_value_of_zero_stays_finite() {
        let point = Association::from_p_value(0, 0.0);
        assert!(point.value.is_finite());
        assert!(
            point.value > 300.0,
            "and is still enormous: {}",
            point.value
        );
    }

    #[test]
    fn the_genome_wide_threshold_is_the_conventional_one() {
        let track = ManhattanTrack::new(Vec::new()).genome_wide_threshold();
        let threshold = track.threshold.unwrap();
        assert!((threshold - 7.30103).abs() < 1e-5, "got {threshold}");
    }

    #[test]
    fn banding_alternates_at_each_boundary() {
        let track = ManhattanTrack::new(Vec::new()).bands(vec![0u64, 1_000, 1_600]);
        assert_eq!(track.band_of(0), 0);
        assert_eq!(track.band_of(999), 0);
        assert_eq!(track.band_of(1_000), 1);
        assert_eq!(track.band_of(1_599), 1);
        assert_eq!(track.band_of(9_000), 2, "past the last boundary");
        // Boundaries arrive in whatever order they are given.
        assert_eq!(
            ManhattanTrack::new(Vec::new())
                .bands(vec![1_600u64, 0, 1_000])
                .band_of(1_200),
            1
        );
    }

    #[test]
    fn banded_sequences_come_out_in_two_shades() {
        let points = vec![Association::new(10, 3.0), Association::new(1_200, 3.0)];
        let plain = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points.clone()))
            .to_svg();
        let banded = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).bands(vec![0u64, 1_000]))
            .to_svg();
        assert_ne!(plain, banded, "nothing alternated");
        let base = Theme::light().muted;
        assert!(banded.contains(&mix(&base, Theme::light().surface(), 0.42)));
    }

    #[test]
    fn one_boundary_is_not_a_banding() {
        // A single sequence has nothing to alternate with.
        let points = vec![Association::new(10, 3.0)];
        let one = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points.clone()).bands(vec![0u64]))
            .to_svg();
        let none = Figure::new(Region::new("genome", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(ManhattanTrack::new(points))
            .to_svg();
        assert_eq!(one, none);
    }

    /// The text of every `<text>` in a document, in order.
    fn words(svg: &str) -> Vec<String> {
        svg.split("<text")
            .skip(1)
            .filter_map(|piece| {
                let body = &piece[piece.find('>')? + 1..];
                Some(body[..body.find("</text>")?].to_string())
            })
            .collect()
    }

    /// A dashed line with nothing on it had to be explained to every reader,
    /// and `-log10 p` after the top tick read as ten minus something.
    #[test]
    fn the_threshold_says_its_value_and_the_axis_says_what_it_measures() {
        let points = vec![Association::new(100, 12.0), Association::new(900, 1.0)];
        let drawn = |track: ManhattanTrack| {
            words(&Figure::new(region()).push(track.label("scan")).to_svg())
        };
        let genome_wide = drawn(
            ManhattanTrack::new(points.clone())
                .genome_wide_threshold()
                .axis_title("-log10 p"),
        );
        assert!(
            genome_wide.contains(&"p = 5e-8".to_string()),
            "{genome_wide:?}"
        );
        assert!(
            genome_wide.contains(&"-log10 p".to_string()),
            "{genome_wide:?}"
        );
        assert!(
            genome_wide.iter().all(|word| !word.ends_with(" -log10 p")),
            "a tick carries the title: {genome_wide:?}"
        );
        let own = drawn(ManhattanTrack::new(points.clone()).p_value_threshold(1e-5));
        assert!(own.contains(&"p = 1e-5".to_string()), "{own:?}");
        let loose = drawn(ManhattanTrack::new(points.clone()).p_value_threshold(0.05));
        assert!(loose.contains(&"p = 0.05".to_string()), "{loose:?}");
        // In the axis's own units, a line says the number it is at.
        let plain = drawn(ManhattanTrack::new(points.clone()).threshold(6.0));
        assert!(plain.contains(&"6".to_string()), "{plain:?}");
        // Unless it is told to say something else, or nothing.
        let named = drawn(
            ManhattanTrack::new(points.clone())
                .threshold(6.0)
                .threshold_label("FDR 5%"),
        );
        assert!(named.contains(&"FDR 5%".to_string()), "{named:?}");
        let quiet = drawn(
            ManhattanTrack::new(points.clone())
                .genome_wide_threshold()
                .threshold_label(""),
        );
        assert!(
            !quiet.iter().any(|word| word.starts_with("p =")),
            "{quiet:?}"
        );
        // A p outside nought to one draws no line.
        assert_eq!(
            ManhattanTrack::new(points.clone())
                .p_value_threshold(0.0)
                .threshold,
            None
        );
        assert_eq!(
            ManhattanTrack::new(points.clone())
                .p_value_threshold(2.0)
                .threshold,
            None
        );
        // No scale, no title for it.
        let hidden = ManhattanTrack::new(points)
            .axis_title("-log10 p")
            .show_scale(false);
        assert_eq!(Track::axis_title(&hidden), None);
    }

    /// A line along the top of the band has no room over it, and a label
    /// pushed down onto it was struck through.
    #[test]
    fn a_threshold_at_the_top_is_labelled_under_its_line() {
        let track =
            ManhattanTrack::new(vec![Association::new(100, 7.0), Association::new(900, 1.0)])
                .genome_wide_threshold();
        let svg = Figure::new(region()).push(track).to_svg();
        let line_y: f64 = {
            let at = svg.find("stroke-dasharray").expect("a dashed line");
            let tag = &svg[svg[..at].rfind("<line").expect("the line")..at];
            let y1 = &tag[tag.find("y1=\"").unwrap() + 4..];
            y1[..y1.find('"').unwrap()].parse().unwrap()
        };
        let label_y: f64 = {
            let at = svg.find(">p = 5e-8<").expect("the label");
            let tag = &svg[svg[..at].rfind("<text").unwrap()..at];
            let y = &tag[tag.find(" y=\"").unwrap() + 4..];
            y[..y.find('"').unwrap()].parse().unwrap()
        };
        assert!(
            label_y > line_y,
            "the label at {label_y} is over the line at {line_y}"
        );
    }

    #[test]
    fn nothing_is_significant_without_a_threshold() {
        let points = vec![Association::new(1, 50.0), Association::new(2, 0.1)];
        let track = ManhattanTrack::new(points.clone());
        assert!(track.significant().is_empty());
        assert_eq!(
            ManhattanTrack::new(points)
                .threshold(1.0)
                .significant()
                .len(),
            1
        );
    }

    #[test]
    fn the_ceiling_leaves_headroom_but_a_pin_does_not() {
        let points = vec![Association::new(1, 10.0)];
        let fitted = ManhattanTrack::new(points.clone());
        let pinned = ManhattanTrack::new(points).max(10.0);
        assert!(fitted.ceiling() > 10.0);
        assert_eq!(pinned.ceiling(), 10.0);
    }

    #[test]
    fn a_threshold_above_every_point_still_fits_on_the_axis() {
        let track = ManhattanTrack::new(vec![Association::new(1, 2.0)]).threshold(7.3);
        assert!(track.ceiling() >= 7.3, "the line would be off the top");
    }

    #[test]
    fn an_empty_track_asks_for_no_axis_and_draws_a_baseline() {
        let track = ManhattanTrack::new(Vec::new());
        assert_eq!(track.y_axis_width(&Theme::light()), 0.0);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("<line"));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn hits_and_misses_are_drawn_in_different_colours() {
        let points = vec![Association::new(1_000, 1.0), Association::new(2_000, 9.0)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                ManhattanTrack::new(points)
                    .threshold(5.0)
                    .color("#111111")
                    .significant_color("#ff0000"),
            )
            .to_svg();
        assert!(svg.contains("#111111"), "the quiet point");
        assert!(svg.contains("#ff0000"), "the hit and the threshold line");
    }

    #[test]
    fn non_finite_values_are_skipped_rather_than_drawn_at_zero() {
        let points = vec![
            Association::new(1_000, f64::NAN),
            Association::new(2_000, 4.0),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).show_scale(false))
            .to_svg();
        assert_eq!(svg.matches("<circle").count(), 1);
    }

    #[test]
    fn points_outside_the_region_are_not_drawn() {
        let points = vec![Association::new(50_000, 4.0)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).show_scale(false))
            .to_svg();
        assert!(!svg.contains("<circle"));
    }

    #[test]
    fn marks_are_bounded_by_the_pixels_and_not_by_the_points() {
        // Every test was drawn before this, so the band here held a hundred
        // thousand circles, more than six for every pixel it has, and a scan
        // of that size through the command line was 5.9 MB. A point drawn over
        // one of its own look on the same pixel adds an element to the
        // document and nothing to the picture.
        let points: Vec<Association> = (0..100_000u64)
            .map(|i| Association::new(i * 10, (i * 7_919 % 1_000) as f64 / 100.0))
            .collect();
        let (width, height) = (400.0, 40.0);
        let svg = Figure::new(Region::new("chr1", 0, 1_000_000).unwrap())
            .width(width)
            .show_region_label(false)
            .push(ManhattanTrack::new(points).height(height).show_scale(false))
            .to_svg();

        // The band is no wider than the figure and exactly as tall as asked,
        // so this is every pixel a mark can be centred on.
        let pixels = ((width + 1.0) * (height + 1.0)) as usize;
        let marks = svg.matches("<circle").count();
        assert!(
            marks <= pixels,
            "{marks} marks for a band of at most {pixels} pixels"
        );
    }

    #[test]
    fn a_pixel_keeps_one_mark_of_each_look_and_the_one_on_top() {
        // A miss, a hit and a second miss on one pixel. The hit and the misses
        // look different, so both looks are drawn; the second miss was drawn
        // over the first and over the hit, so it is the miss kept, and it
        // still goes on after the hit. Keeping the first of each look instead
        // would bring the hit in front of what used to cover it.
        let at = 5_000;
        let points = vec![
            Association::new(at, 5.0 - 1e-9),
            Association::new(at, 5.0),
            Association::new(at, 5.0 - 2e-9),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(ManhattanTrack::new(points).threshold(5.0).show_scale(false))
            .to_svg();

        assert_eq!(svg.matches("<circle").count(), 1, "one miss on one pixel");
        assert_eq!(
            svg.matches("<polygon").count(),
            2,
            "the hit, a diamond and its ring"
        );
        let hit = svg.rfind("<polygon").unwrap();
        let miss = svg.rfind("<circle").unwrap();
        assert!(
            miss > hit,
            "the hit was brought in front of the miss over it"
        );
    }

    #[test]
    fn a_pixel_holds_one_mark_per_look_and_not_one_per_sequence() {
        // Four sequences a base long after the first, so a point on the first,
        // the second and the fifth all land on one pixel. The first and the
        // fifth are drawn alike, a circle in the plain shade, and the second
        // is a square a shade lighter, which is a mark of its own.
        let points = vec![
            Association::new(999, 3.0),
            Association::new(1_000, 3.0),
            Association::new(1_003, 3.0),
        ];
        let svg = Figure::new(Region::new("genome", 0, 100_000).unwrap())
            .show_region_label(false)
            .push(
                ManhattanTrack::new(points)
                    .bands(vec![0u64, 1_000, 1_001, 1_002, 1_003])
                    .show_scale(false),
            )
            .to_svg();

        let plain = Theme::light().muted;
        let lighter = mix(&plain, Theme::light().surface(), 0.42);
        assert_eq!(
            svg.matches(&format!("fill=\"{plain}\"")).count(),
            1,
            "two circles in one shade on one pixel"
        );
        assert_eq!(
            svg.matches(&format!("fill=\"{lighter}\"")).count(),
            1,
            "the lighter square went with them"
        );
    }
}
