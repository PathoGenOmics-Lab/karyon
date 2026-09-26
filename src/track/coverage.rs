//! Per-base quantitative signal: read depth, GC content, mappability.
//!
//! Dense values from a start position are what [`CoverageTrack::new`] takes;
//! scattered points come in through [`CoverageTrack::from_pairs`], and spans
//! as a bedGraph states them through [`CoverageTrack::from_spans`]. The profile
//! is drawn upwards from the floor of the band, so the quantity has to be one
//! whose zero is the bottom in fact and not by convention: a signed or centred
//! statistic belongs in [`WindowTrack`](crate::WindowTrack).
//!
//! # What a pixel column throws away
//!
//! A region wider than the figure has more bases than columns, so each column
//! is a summary of the bases beneath it. Which summary is [`Aggregate`], and
//! the default is [`Aggregate::Max`], so a single base spike survives being
//! drawn at a megabase.
//!
//! The cost is the mirror image, and it is quiet: a maximum cannot show a hole.
//! Across a 5 Mb view one pixel spans thousands of bases, and a 300 bp deletion
//! inside an otherwise well covered gene contributes nothing to the maximum of
//! its column, so the profile is drawn flat over it. [`Aggregate::Min`] takes
//! the floor of each column instead, and the dropout appears. A column with
//! nothing under it at all is skipped rather than drawn at zero, so a gap in
//! the input stays a gap in the profile.
//!
//! # The ceiling moves unless it is pinned
//!
//! Left alone, the top of the band is the largest value on screen, so the
//! profile fills the band whatever the depth is. [`CoverageTrack::max`] fixes
//! the top and is taken literally; an automatic ceiling is lifted by six per
//! cent instead, so the tallest point reads as a peak rather than as something
//! that ran out of room. [`CoverageTrack::log_scale`] is the other way to make
//! a wide range fit into one band.

use crate::region::Region;
use crate::scale::Scale;
use crate::style::{legible_ticks, Emphasis, QuantitativeAxis};
use crate::svg::{num, text_width, Anchor};
use crate::theme::Theme;
use crate::track::{unbroken, DrawContext, Track};

/// How a coverage track is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageStyle {
    /// Filled area under the profile. The default, and the one that survives
    /// being printed small.
    Area,
    /// Outline only, for overlaying several samples.
    Line,
    /// One bar per pixel column, closest to what a genome browser shows.
    Bars,
}

/// How values are combined when one pixel covers many bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aggregate {
    /// Highest value in the column. Keeps narrow spikes visible, which is
    /// usually what you want when hunting for duplications, and for the same
    /// reason hides a dropout narrower than a column: the covered bases beside
    /// it still set the column's height. [`Aggregate::Min`] is the one that
    /// shows it.
    Max,
    /// Mean of the column. Truer to the overall level, hides single-base
    /// spikes.
    Mean,
    /// Lowest value in the column. The one to use when looking for dropouts
    /// that a mean would smooth away.
    Min,
}

/// A quantitative signal sampled once per base.
///
/// Values are kept from a start position as runs of bases holding one value,
/// so a stretch of one depth costs one entry however long it is. When a pixel
/// covers more than one base the column is reduced with [`Aggregate`]; the SVG
/// therefore has at most one point per pixel however wide the region is.
///
/// ```
/// use karyon::{CoverageTrack, Figure, Region};
///
/// let depth: Vec<f64> = (0..1000).map(|i| 30.0 + (i % 17) as f64).collect();
/// let svg = Figure::new(Region::parse("chr1:1-1000").unwrap())
///     .push(CoverageTrack::new(0, depth).label("depth"))
///     .to_svg();
/// assert!(svg.contains("depth"));
/// ```
#[derive(Debug, Clone)]
pub struct CoverageTrack {
    start: u64,
    /// The values as runs of bases holding one value, each ending where the
    /// next begins, the first beginning at `start`. A value per base was a
    /// double for every base of the region: two gigabytes to draw a whole
    /// chromosome from a bedGraph of a few thousand windows.
    runs: Vec<Run>,
    label: Option<String>,
    height: f64,
    color: Option<String>,
    style: CoverageStyle,
    aggregate: Aggregate,
    max: Option<f64>,
    log_scale: bool,
    fill_opacity: Option<f64>,
    show_max: bool,
    axis: QuantitativeAxis,
    title: Option<String>,
}

impl CoverageTrack {
    /// A track whose `values[i]` describes base `start + i`, 0-based.
    ///
    /// Values need not cover the whole region: anything outside is not
    /// drawn, and non-finite values are treated as missing.
    pub fn new(start: u64, values: impl Into<Vec<f64>>) -> Self {
        let mut runs = Vec::new();
        let mut at = start;
        for value in values.into() {
            at = at.saturating_add(1);
            push_run(&mut runs, at, value);
        }
        CoverageTrack {
            start,
            runs,
            label: None,
            height: 60.0,
            color: None,
            style: CoverageStyle::Area,
            aggregate: Aggregate::Max,
            max: None,
            log_scale: false,
            fill_opacity: None,
            show_max: true,
            axis: QuantitativeAxis::new(),
            title: None,
        }
    }

    /// A track built from sparse `(position, value)` pairs.
    ///
    /// Positions are 0-based. The profile spans `region`, and memory follows
    /// the changes of value in it rather than its length. Pairs outside the
    /// region are ignored, and positions not listed stay at zero.
    pub fn from_pairs(region: &Region, pairs: impl IntoIterator<Item = (u64, f64)>) -> Self {
        Self::from_spans(
            region,
            pairs.into_iter().map(|(pos, value)| (pos, pos + 1, value)),
        )
    }

    /// A profile from half-open `(start, end, value)` spans over `region`.
    ///
    /// What a bedGraph states, taken as it states it. A row covering a hundred
    /// thousand bases is one span here and not a hundred thousand pairs, which
    /// is the difference between a kilobyte of input costing a kilobyte and it
    /// costing six gigabytes.
    ///
    /// A base no span covers stays at nought, the same as
    /// [`CoverageTrack::from_pairs`], because that is what a depth of nought
    /// means and what a bedGraph leaves out.
    pub fn from_spans(region: &Region, spans: impl IntoIterator<Item = (u64, u64, f64)>) -> Self {
        let mut track = CoverageTrack::new(region.start(), Vec::new());
        if region.end() > region.start() {
            track.runs.push(Run {
                end: region.end(),
                value: 0.0,
            });
        }
        for (start, end, value) in spans {
            track.paint(start, end, value);
        }
        track
    }

    /// Writes `value` over the bases from `start` up to `end`, clamped to what
    /// this track covers.
    ///
    /// This is what [`CoverageTrack::from_spans`] does to each span in turn,
    /// separated out so a reader can hand over one span at a time and never
    /// hold the whole list. `samtools depth` writes a line per base, and over
    /// ten million bases collecting them first cost 231 MB for a track of 76.
    pub(crate) fn paint(&mut self, start: u64, end: u64, value: f64) {
        let (lo, hi) = (start.max(self.start), end.min(self.end()));
        if lo >= hi {
            return;
        }
        // A file in order paints inside the last run, the stretch it has not
        // reached yet, so that run is split where it lies, with no search and
        // nothing allocated. `samtools depth` paints ten million bases this
        // way, and a search for each one made that figure almost twice as
        // slow as the value per base these runs replaced.
        let tail = self.runs.len() - 1;
        let tail_start = self.run_start(tail);
        if lo >= tail_start {
            let Run { end, value: under } = self.runs[tail];
            self.runs.pop();
            if tail_start < lo {
                push_run(&mut self.runs, lo, under);
            }
            push_run(&mut self.runs, hi, value);
            if hi < end {
                push_run(&mut self.runs, end, under);
            }
            return;
        }
        // Anywhere else, the run holding `lo` and the run holding the base
        // before `hi` are replaced by what the span leaves of them.
        let first = self.runs.partition_point(|run| run.end <= lo);
        let last = self.runs.partition_point(|run| run.end < hi);
        let mut replacement = Vec::with_capacity(3);
        let before = self.run_start(first);
        if before < lo {
            replacement.push(Run {
                end: lo,
                value: self.runs[first].value,
            });
        }
        replacement.push(Run { end: hi, value });
        if hi < self.runs[last].end {
            replacement.push(Run {
                end: self.runs[last].end,
                value: self.runs[last].value,
            });
        }
        self.runs.splice(first..=last, replacement);
        // Runs of one value either side of the new one become one.
        let from = first.saturating_sub(1);
        let mut index = from;
        while index + 1 < self.runs.len() && index <= from + 3 {
            if same(self.runs[index].value, self.runs[index + 1].value) {
                self.runs[index].end = self.runs[index + 1].end;
                self.runs.remove(index + 1);
            } else {
                index += 1;
            }
        }
    }

    /// The value at a 0-based position, or `None` outside what this track holds.
    ///
    /// A track laid over a region holds a value for every base of it, so a
    /// `Some(0.0)` here is a nought the file stated or a base it left out, and
    /// `None` is a position this track does not reach at all.
    pub fn at(&self, pos: u64) -> Option<f64> {
        if pos < self.start {
            return None;
        }
        let index = self.runs.partition_point(|run| run.end <= pos);
        self.runs.get(index).map(|run| run.value)
    }

    /// One past the last base this track holds a value for.
    fn end(&self) -> u64 {
        self.runs.last().map_or(self.start, |run| run.end)
    }

    /// Where run `index` begins, which is where the one before it ends.
    fn run_start(&self, index: usize) -> u64 {
        match index {
            0 => self.start,
            _ => self.runs[index - 1].end,
        }
    }

    /// The runs over `[lo, hi)`, each cut to the part inside, as
    /// `(from, to, value)`.
    fn overlapping(&self, lo: u64, hi: u64) -> impl Iterator<Item = (u64, u64, f64)> + '_ {
        let first = self.runs.partition_point(|run| run.end <= lo);
        (first..self.runs.len()).map_while(move |index| {
            let from = self.run_start(index);
            (from < hi).then(|| {
                let run = self.runs[index];
                (from.max(lo), run.end.min(hi), run.value)
            })
        })
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(1.0);
        self
    }

    /// Overrides the colour, which defaults to the theme accent.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// What the value axis measures, written under the track's name, such
    /// as `cM/Mb` for a recombination rate.
    pub fn axis_title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Chooses area, line or bars.
    pub fn style(mut self, style: CoverageStyle) -> Self {
        self.style = style;
        self
    }

    /// Chooses how a pixel column summarises the bases under it.
    pub fn aggregate(mut self, aggregate: Aggregate) -> Self {
        self.aggregate = aggregate;
        self
    }

    /// Pins the top of the y axis, instead of taking the visible maximum.
    ///
    /// Pin it when comparing samples side by side, or the eye will read two
    /// different scales as the same one.
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

    /// Plots `log10(1 + value)` so that a 1000x depth range fits in one band.
    ///
    /// The axis label still reports the value in the original units.
    pub fn log_scale(mut self, enabled: bool) -> Self {
        self.log_scale = enabled;
        self
    }

    /// Sets the fill opacity, between 0 and 1.
    ///
    /// Left alone, each style picks what suits it: an area is a wash under a
    /// drawn line, so it fills at a fifth and lets whatever is behind it show
    /// through, while bars are the mark itself and stay solid.
    pub fn fill_opacity(mut self, opacity: f64) -> Self {
        self.fill_opacity = Some(opacity.clamp(0.0, 1.0));
        self
    }

    /// Shows or hides the small maximum-value annotation.
    pub fn show_max(mut self, show: bool) -> Self {
        self.show_max = show;
        self
    }

    /// Largest finite value inside `region`, or `None` when nothing overlaps.
    pub fn visible_max(&self, region: &Region) -> Option<f64> {
        let (lo, hi) = self.visible_slice(region)?;
        self.overlapping(lo, hi)
            .map(|(_, _, value)| value)
            .filter(|v| v.is_finite())
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            })
    }

    /// The stretch of positions this track holds that `region` overlaps.
    fn visible_slice(&self, region: &Region) -> Option<(u64, u64)> {
        let end = self.end();
        if region.end() <= self.start || region.start() >= end {
            return None;
        }
        let lo = region.start().max(self.start);
        let hi = region.end().min(end);
        (hi > lo).then_some((lo, hi))
    }

    /// Reduces the bases in `[lo, hi)` to the single value a pixel shows.
    fn sample(&self, lo: f64, hi: f64) -> Option<f64> {
        let origin = self.start as f64;
        let first = (lo - origin).floor().max(0.0);
        let last = (hi - origin).ceil().min((self.end() - self.start) as f64);
        if last <= first {
            return None;
        }
        let (from, to) = (self.start + first as u64, self.start + last as u64);
        let mut count = 0u64;
        let mut acc = match self.aggregate {
            Aggregate::Max => f64::NEG_INFINITY,
            Aggregate::Min => f64::INFINITY,
            Aggregate::Mean => 0.0,
        };
        for (run_from, run_to, value) in self.overlapping(from, to) {
            if !value.is_finite() {
                continue;
            }
            // Every base of the run counts once, as it did one at a time.
            let bases = run_to - run_from;
            count += bases;
            acc = match self.aggregate {
                Aggregate::Max => acc.max(value),
                Aggregate::Min => acc.min(value),
                Aggregate::Mean => acc + value * bases as f64,
            };
        }
        if count == 0 {
            return None;
        }
        Some(match self.aggregate {
            Aggregate::Mean => acc / count as f64,
            _ => acc,
        })
    }

    fn transform(&self, value: f64) -> f64 {
        if self.log_scale {
            (1.0 + value.max(0.0)).log10()
        } else {
            value
        }
    }
}

/// A stretch of bases holding one value, from where the run before it ends to
/// `end`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Run {
    end: u64,
    value: f64,
}

/// Whether two values are one: the same bits, so two missing values are one
/// value and a nought is not a missing one.
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// Adds a run ending at `end`, joined to the one before it where it holds the
/// same value.
fn push_run(runs: &mut Vec<Run>, end: u64, value: f64) {
    match runs.last_mut() {
        Some(last) if same(last.value, value) => last.end = end,
        _ => runs.push(Run { end, value }),
    }
}

/// How many samples a band of `width` pixels is reduced to, and how far apart.
///
/// One per pixel column, which is what the profile is drawn from, up to a
/// ceiling. A width is a number a caller supplies, and turning it straight
/// into a column count turned a large finite width into an allocation that
/// aborted the render: `1e30` pixels is not an image, but it is a `f64` a
/// builder will accept. Past the ceiling the samples are spread over the whole
/// band instead of stopping part way across it, so the profile still spans the
/// band it was given. Below it the step is exactly one pixel and nothing about
/// the drawing changes.
fn column_grid(width: f64) -> (usize, f64) {
    const MAX_COLUMNS: usize = 100_000;
    let wanted = width.max(1.0).ceil();
    // A NaN width casts to zero, so the floor is applied after the cast too.
    let columns = (wanted as usize).clamp(1, MAX_COLUMNS);
    let step = (wanted / columns as f64).max(1.0);
    (columns, step)
}

impl Track for CoverageTrack {
    fn noun(&self) -> &str {
        "a coverage profile"
    }

    fn axis_title(&self) -> Option<&str> {
        self.show_max.then_some(self.title.as_deref()).flatten()
    }

    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_max {
            return 0.0;
        }
        // Room for the widest label this track could print. The data it holds
        // bounds every region it can be shown over, so the ends it would
        // round to over all of it bound the labels over any part of it.
        let size = theme.font_size - 1.0;
        let whole = self
            .runs
            .iter()
            .map(|run| run.value)
            .filter(|v| v.is_finite())
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            });
        let widest = match self
            .axis
            .max
            .or(self.max)
            .or(whole)
            .filter(|m| m.is_finite() && *m > 0.0)
        {
            Some(ceiling) => {
                let (floor, ceiling, _) = self.ends(ceiling);
                self.tick_values(floor, ceiling)
                    .iter()
                    .zip(self.tick_labels(floor, ceiling))
                    .map(|(_, label)| text_width(&label, size))
                    .fold(0.0f64, f64::max)
            }
            None => text_width("0", size),
        };
        widest + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());

        let baseline = band.bottom();
        // Half a stroke inside the band. A one pixel rule centred on the clip
        // edge keeps only half its ink and renders at half the weight of the
        // ceiling rule drawn with identical parameters, which AxisTrack
        // already worked around with the same half pixel.
        ctx.svg.line(
            band.x,
            baseline - 0.5,
            band.right(),
            baseline - 0.5,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );

        let data_ceiling = self
            .axis
            .max
            .or(self.max)
            .or_else(|| self.visible_max(ctx.region))
            .filter(|m| m.is_finite() && *m > 0.0);
        let Some(data_ceiling) = data_ceiling else {
            return;
        };
        let (floor, ceiling, visual_ceiling) = self.ends(data_ceiling);
        let transformed_floor = self.transform(floor);
        let span = self.transform(visual_ceiling) - transformed_floor;
        if span <= 0.0 {
            return;
        }
        let y_of = |value: f64| {
            baseline - ((self.transform(value) - transformed_floor) / span).clamp(0.0, 1.0) * band.h
        };

        let (columns, step) = column_grid(band.w);
        let mut points: Vec<Option<(f64, f64)>> = Vec::with_capacity(columns);
        for column in 0..columns {
            let x = band.x + column as f64 * step;
            let lo = ctx.scale.pos_at_x(x);
            let hi = ctx.scale.pos_at_x(x + step);
            let Some(value) = self.sample(lo, hi) else {
                points.push(None);
                continue;
            };
            let y = y_of(value);
            match self.style {
                CoverageStyle::Bars => {
                    ctx.svg.rect_opacity(
                        x,
                        y,
                        step,
                        baseline - y,
                        &color,
                        self.fill_opacity.unwrap_or(1.0),
                    );
                }
                _ => points.push(Some((x, y))),
            }
        }

        // The profile is drawn once per unbroken run of columns. A column with
        // nothing under it is a gap, and drawn as one line the profile ran
        // straight across it, from the last value before the gap to the first
        // one after, at a depth nothing measured.
        if self.style != CoverageStyle::Bars {
            for mut run in unbroken(points) {
                // A column standing alone between two gaps still holds a
                // value, and a line needs two ends, so it is drawn across its
                // own width rather than dropped.
                if let [(x, y)] = run[..] {
                    run.push((x + step, y));
                }
                if self.style == CoverageStyle::Area {
                    let mut d = String::with_capacity(run.len() * 14);
                    d.push('M');
                    d.push_str(&num(run[0].0));
                    d.push(' ');
                    d.push_str(&num(baseline));
                    for (x, y) in &run {
                        d.push_str(" L");
                        d.push_str(&num(*x));
                        d.push(' ');
                        d.push_str(&num(*y));
                    }
                    let last = run[run.len() - 1].0;
                    d.push_str(" L");
                    d.push_str(&num(last));
                    d.push(' ');
                    d.push_str(&num(baseline));
                    d.push('Z');
                    // A fade under a drawn line, rather than a saturated block.
                    // The line is what carries the shape; the fill only says
                    // which side of it is under the curve, and it is strongest
                    // up against the line and gone by the baseline, which is
                    // where there is least to say. An opacity asked for is a
                    // flat wash of exactly that, since whoever asked for it
                    // wanted that number.
                    match self.fill_opacity {
                        Some(opacity) => ctx.svg.path(&d, &color, opacity),
                        None => {
                            let top = (ctx.theme.tokens.area_opacity * 3.0).min(1.0);
                            let paint = ctx.svg.fade_down(&color, top, 0.03);
                            ctx.svg.path(&d, &paint, 1.0);
                        }
                    }
                }
                ctx.svg
                    .polyline(&run, &color, ctx.theme.tokens.stroke * 1.25);
            }
        }

        if self.show_max {
            self.draw_axis(ctx, floor, ceiling, visual_ceiling);
        }
    }
}

impl CoverageTrack {
    /// The floor, the labelled ceiling and the ceiling the band is scaled to.
    ///
    /// A free ceiling is rounded up to a value worth labelling, which is also
    /// the headroom that keeps the tallest point a peak rather than something
    /// that ran out of band. A pinned maximum is taken literally, because that
    /// is the whole reason for pinning one.
    fn ends(&self, data_ceiling: f64) -> (f64, f64, f64) {
        let (floor, ceiling) = self.axis.resolve(0.0, data_ceiling);
        let pinned = self.axis.max.is_some() || self.max.is_some();
        if pinned {
            return (floor, ceiling, ceiling);
        }
        if self.log_scale {
            return (floor, ceiling, floor + (ceiling - floor) * 1.06);
        }
        let (floor, rounded) = self.axis.nice(floor, ceiling);
        (floor, rounded, rounded)
    }

    /// Where the value axis puts its ticks.
    ///
    /// A log scale spends its height on the small values, so round linear
    /// steps would pile up at the top; it is ticked at powers of ten instead.
    fn tick_values(&self, floor: f64, ceiling: f64) -> Vec<f64> {
        if !self.log_scale {
            return self.axis.values(floor, ceiling);
        }
        let mut values = vec![floor.max(0.0)];
        let mut power = 1.0;
        while power < ceiling * 0.8 {
            if power > floor {
                values.push(power);
            }
            power *= 10.0;
        }
        values.push(ceiling);
        values
    }

    fn tick_labels(&self, floor: f64, ceiling: f64) -> Vec<String> {
        let values = self.tick_values(floor, ceiling);
        let mut labels = self.axis.labels(&values);
        if self.log_scale {
            if let Some(last) = labels.last_mut() {
                last.push_str(" log");
            }
        }
        labels
    }

    /// Draws the value axis in the strip the figure reserved for it.
    ///
    /// Round values from nought to the rounded ceiling, as many as the band
    /// has room to label, each with a hairline across the plot. A coverage
    /// track is read for its shape and its order of magnitude, so a short band
    /// is left with its two ends rather than a ladder of gridlines that would
    /// be more ink than the profile it is measuring.
    fn draw_axis(&self, ctx: &mut DrawContext<'_>, floor: f64, ceiling: f64, visual_ceiling: f64) {
        let band = ctx.band;
        let size = ctx.theme.font_size - 1.0;
        let baseline = band.bottom();
        let transformed_floor = self.transform(floor);
        let span = self.transform(visual_ceiling) - transformed_floor;
        let y_of =
            |value: f64| baseline - ((self.transform(value) - transformed_floor) / span) * band.h;
        let all = self.tick_values(floor, ceiling);
        let labels = self.tick_labels(floor, ceiling);
        let shown = legible_ticks(&all, y_of, size);
        for &value in &shown {
            if value <= floor {
                continue;
            }
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
            if !reference.value.is_finite() || reference.value < floor || reference.value > ceiling
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

        if ctx.axis.w <= 0.0 {
            return;
        }
        let right = ctx.axis.right() - 4.0;
        for (value, label) in all.iter().zip(&labels) {
            if !shown.contains(value) {
                continue;
            }
            let y = y_of(*value);
            ctx.svg.text(
                right,
                (y + size * 0.35).max(band.y + size * 0.78).min(baseline),
                label,
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

/// Compact label for the top of the y axis.
///
/// A depth of 72.04 is labelled `72`, not `72.0`: the extra digit is noise from
/// however the values were computed, and it reads as precision the data does
/// not have.
#[cfg(test)]
fn format_value(value: f64) -> String {
    if value >= 1_000_000.0 {
        format!("{}M", crate::svg::text_rounded(value / 1e6, 1))
    } else if value >= 1_000.0 {
        format!("{}k", crate::svg::text_rounded(value / 1e3, 1))
    } else {
        crate::svg::text_rounded(value, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> Region {
        Region::new("chr1", 0, 100).unwrap()
    }

    #[test]
    fn sample_takes_the_maximum_of_a_binned_column() {
        let track = CoverageTrack::new(0, vec![1.0, 9.0, 2.0, 3.0]);
        assert_eq!(track.sample(0.0, 4.0), Some(9.0));
    }

    #[test]
    fn sample_can_average_or_take_the_minimum() {
        let values = vec![1.0, 9.0, 2.0, 4.0];
        let mean = CoverageTrack::new(0, values.clone()).aggregate(Aggregate::Mean);
        assert_eq!(mean.sample(0.0, 4.0), Some(4.0));
        let min = CoverageTrack::new(0, values).aggregate(Aggregate::Min);
        assert_eq!(min.sample(0.0, 4.0), Some(1.0));
    }

    #[test]
    fn sample_returns_one_base_when_zoomed_past_single_base_resolution() {
        let track = CoverageTrack::new(10, vec![5.0, 6.0, 7.0]);
        assert_eq!(track.sample(11.2, 11.6), Some(6.0));
    }

    #[test]
    fn sample_outside_the_data_is_missing_rather_than_zero() {
        let track = CoverageTrack::new(10, vec![5.0, 6.0]);
        assert_eq!(track.sample(0.0, 5.0), None);
        assert_eq!(track.sample(50.0, 60.0), None);
    }

    #[test]
    fn non_finite_values_are_ignored_not_propagated() {
        let track = CoverageTrack::new(0, vec![f64::NAN, 4.0, f64::INFINITY]);
        assert_eq!(track.sample(0.0, 3.0), Some(4.0));
        assert_eq!(track.visible_max(&region()), Some(4.0));
    }

    #[test]
    fn a_fully_missing_column_has_no_maximum() {
        let track = CoverageTrack::new(0, vec![f64::NAN, f64::NAN]);
        assert_eq!(track.sample(0.0, 2.0), None);
        assert_eq!(track.visible_max(&region()), None);
    }

    #[test]
    fn visible_max_only_looks_inside_the_region() {
        let track = CoverageTrack::new(0, vec![1.0, 2.0, 100.0, 3.0]);
        let window = Region::new("chr1", 0, 2).unwrap();
        assert_eq!(track.visible_max(&window), Some(2.0));
    }

    #[test]
    fn data_disjoint_from_the_region_has_no_visible_slice() {
        let track = CoverageTrack::new(1000, vec![1.0, 2.0]);
        assert_eq!(track.visible_slice(&region()), None);
        assert_eq!(track.visible_max(&region()), None);
    }

    #[test]
    fn from_pairs_places_values_at_their_positions() {
        let window = Region::new("chr1", 100, 105).unwrap();
        let track = CoverageTrack::from_pairs(&window, [(101, 7.0), (104, 3.0), (900, 99.0)]);
        let values: Vec<f64> = (100..105).map(|pos| track.at(pos).unwrap()).collect();
        assert_eq!(values, vec![0.0, 7.0, 0.0, 0.0, 3.0]);
        assert_eq!(track.visible_max(&window), Some(7.0));
    }

    /// Spans painted in order, as a file writes them, or in any order at all,
    /// leave every base holding what a value per base would, and no two
    /// neighbouring runs holding one value.
    #[test]
    fn painting_spans_leaves_what_a_value_per_base_would() {
        let window = Region::new("chr1", 100, 400).unwrap();
        // A generator of its own, since the crate takes no dependencies.
        let mut seed = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = |below: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % below
        };
        for round in 0..200 {
            let mut track = CoverageTrack::from_spans(&window, std::iter::empty());
            let mut expected = [0.0f64; 300];
            let mut reached = 90;
            for _ in 0..40 {
                let (start, end) = if round % 2 == 0 {
                    let start = reached + next(4);
                    reached = start + 1 + next(12);
                    (start, reached)
                } else {
                    let start = 90 + next(330);
                    (start, start + 1 + next(60))
                };
                let value = [0.0, 1.0, 2.0, f64::NAN][next(4) as usize];
                track.paint(start, end, value);
                for pos in start.max(100)..end.min(400) {
                    expected[(pos - 100) as usize] = value;
                }
            }
            for pos in 100..400 {
                let (held, wanted) = (track.at(pos).unwrap(), expected[(pos - 100) as usize]);
                assert!(
                    same(held, wanted),
                    "round {round}, base {pos}: {held} for {wanted}"
                );
            }
            assert_eq!(track.end(), 400, "round {round}");
            assert!(
                track
                    .runs
                    .windows(2)
                    .all(|pair| pair[0].end < pair[1].end && !same(pair[0].value, pair[1].value)),
                "round {round}: {:?}",
                track.runs
            );
        }
    }

    #[test]
    fn log_scale_compresses_the_range() {
        let track = CoverageTrack::new(0, vec![0.0]).log_scale(true);
        assert!((track.transform(9.0) - 1.0).abs() < 1e-12);
        assert_eq!(track.transform(0.0), 0.0);
    }

    #[test]
    fn the_axis_labels_a_round_ceiling_and_the_floor() {
        use crate::figure::Figure;
        let depth: Vec<f64> = (0..500).map(|i| (i % 87) as f64).collect();
        let svg = Figure::new(Region::parse("chr1:1-500").unwrap())
            .show_region_label(false)
            .push(CoverageTrack::new(0, depth).label("depth"))
            .to_svg();
        assert!(svg.contains(">100</text>"), "the ceiling rounds up to 100");
        assert!(svg.contains(">50</text>"), "with a round step under it");
        assert!(svg.contains(">0</text>"), "and the floor is labelled");
        assert!(!svg.contains(">86</text>"), "not the tallest sample");
    }

    #[test]
    fn the_shared_axis_draws_formatted_ticks_and_labelled_reference_lines() {
        use crate::style::{AxisFormat, ReferenceLine};
        let axis = QuantitativeAxis::new()
            .range(10.0, 30.0)
            .ticks(3)
            .unit("x")
            .format(AxisFormat::Fixed(1))
            .reference(ReferenceLine::new(20.0).label("target"));
        let svg = crate::Figure::new(region())
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![20.0; 100]).axis(axis))
            .to_svg();
        assert!(svg.contains(">30.0x</text>"), "the unit on the top tick");
        assert!(svg.contains(">10.0</text>"), "and nowhere else");
        assert!(!svg.contains(">10.0x</text>"), "not repeated on every tick");
        assert!(svg.contains(">target</text>"));
        assert!(svg.contains("stroke-dasharray"));
    }

    #[test]
    fn a_pinned_maximum_is_what_the_axis_says() {
        use crate::figure::Figure;
        let svg = Figure::new(Region::parse("chr1:1-500").unwrap())
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 500]).max(250.0).label("d"))
            .to_svg();
        assert!(svg.contains(">250</text>"));
        assert!(!svg.contains(">10</text>"));
    }

    #[test]
    fn value_labels_shorten_large_numbers() {
        assert_eq!(format_value(42.0), "42");
        assert_eq!(format_value(1500.0), "1.5k");
        assert_eq!(format_value(2_500_000.0), "2.5M");
        assert_eq!(format_value(0.5), "0.5");
    }

    #[test]
    fn value_labels_drop_a_decimal_that_rounds_to_zero() {
        assert_eq!(format_value(72.04), "72");
        assert_eq!(format_value(72.44), "72.4");
        assert_eq!(format_value(2_000_400.0), "2M");
    }

    #[test]
    fn data_starting_at_the_top_of_the_coordinate_range_is_off_screen_not_a_panic() {
        use crate::figure::Figure;
        // `start + values.len()` used to be computed before the disjointness
        // test below it, so a track at u64::MAX aborted the render.
        let track = CoverageTrack::new(u64::MAX, vec![7.0]);
        assert_eq!(track.visible_slice(&region()), None);
        assert_eq!(track.visible_max(&region()), None);
        let svg = Figure::new(region()).push(track).to_svg();
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn a_band_wider_than_any_image_is_still_sampled_across_its_whole_width() {
        // One sample per pixel up to the ceiling; past it the samples are
        // spread over the band rather than turned into an allocation. 1e30
        // pixels used to abort with a capacity overflow.
        assert_eq!(column_grid(800.4), (801, 1.0));
        assert_eq!(column_grid(0.0), (1, 1.0));
        let (columns, step) = column_grid(1e30);
        assert_eq!(columns, 100_000);
        assert!((columns as f64 * step - 1e30).abs() <= 1e15);
    }

    #[test]
    fn a_figure_a_thousand_times_wider_than_any_screen_still_renders() {
        use crate::figure::Figure;
        let svg = Figure::new(Region::new("chr1", 0, 1000).unwrap())
            .width(1e30)
            .push(CoverageTrack::new(0, vec![1.0, 2.0, 3.0]))
            .to_svg();
        assert!(svg.contains("</svg>"));
    }

    /// A depth over three hundred bases whose middle hundred were never
    /// measured.
    fn holed() -> Vec<f64> {
        (0..300)
            .map(|i| {
                if (100..200).contains(&i) {
                    f64::NAN
                } else {
                    20.0
                }
            })
            .collect()
    }

    /// `track` alone in a figure over `region`, as a document.
    fn drawn(region: Region, track: CoverageTrack) -> String {
        crate::Figure::new(region)
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    /// The vertices of every filled area in a document, one list per area.
    fn areas(svg: &str) -> Vec<Vec<(f64, f64)>> {
        svg.split("<path d=\"")
            .skip(1)
            .filter_map(|rest| rest.split('"').next())
            .map(|d| {
                let numbers: Vec<f64> = d
                    .replace(['M', 'L', 'Z'], " ")
                    .split_whitespace()
                    .filter_map(|piece| piece.parse().ok())
                    .collect();
                numbers
                    .chunks_exact(2)
                    .map(|pair| (pair[0], pair[1]))
                    .collect()
            })
            .collect()
    }

    /// Fails on the first pair of neighbouring vertices more than a pixel
    /// column apart. Each vertex of the profile is a column, so a profile that
    /// stays with the data never steps further than that.
    fn assert_no_column_skipped(shapes: &[Vec<(f64, f64)>]) {
        for shape in shapes {
            for pair in shape.windows(2) {
                assert!(
                    pair[1].0 - pair[0].0 <= 1.0 + 1e-9,
                    "{:?} is joined to {:?}, {} pixels on, over columns that hold nothing",
                    pair[0],
                    pair[1],
                    pair[1].0 - pair[0].0
                );
            }
        }
    }

    #[test]
    fn a_stretch_of_missing_values_is_a_gap_in_the_line() {
        // The hundred bases in the middle have no value, and the line used to
        // run straight from the last column before them to the first one
        // after, drawing a depth across the whole stretch that nothing
        // measured.
        let svg = drawn(
            Region::new("chr1", 0, 300).unwrap(),
            CoverageTrack::new(0, holed()).style(CoverageStyle::Line),
        );
        let lines = crate::track::polylines(&svg);
        assert_no_column_skipped(&lines);
        assert_eq!(lines.len(), 2, "one line either side of the gap");
    }

    #[test]
    fn a_stretch_of_missing_values_is_a_gap_in_the_area() {
        // The fill has the same edge as the line on top of it, so it bridged
        // the same stretch, and filled everything under the bridge as well.
        let svg = drawn(
            Region::new("chr1", 0, 300).unwrap(),
            CoverageTrack::new(0, holed()).style(CoverageStyle::Area),
        );
        let fills = areas(&svg);
        assert_no_column_skipped(&fills);
        assert_eq!(fills.len(), 2, "one area either side of the gap");
        let lines = crate::track::polylines(&svg);
        assert_no_column_skipped(&lines);
        assert_eq!(lines.len(), 2, "and one edge on each");
    }

    #[test]
    fn a_column_standing_alone_is_drawn_rather_than_dropped() {
        // One base in a megabase, so one pixel column holds a value and every
        // other column holds nothing. A line needs two ends, and the one value
        // on screen used to be drawn as no line at all, in both styles, which
        // is a band that says there is nothing there.
        let region = Region::new("chr1", 0, 1_000_000).unwrap();
        for style in [CoverageStyle::Line, CoverageStyle::Area] {
            let svg = drawn(
                region.clone(),
                CoverageTrack::new(500_000, vec![20.0]).style(style),
            );
            let lines = crate::track::polylines(&svg);
            assert_eq!(lines.len(), 1, "{style:?} draws the one value on screen");
            assert!(
                (lines[0][1].0 - lines[0][0].0 - 1.0).abs() < 1e-9,
                "across its own column: {:?}",
                lines[0]
            );
            if style == CoverageStyle::Area {
                assert_eq!(areas(&svg).len(), 1, "and fills under it");
            }
        }
    }
}
