//! A statistic computed in windows, drawn against a line it may fall below.
//!
//! A [`Window`] is a half-open span and one number; the track is the run of
//! them. pN/pS along a gene, Tajima's D along a chromosome, GC skew around a
//! plasmid: all of them arrive in that shape, and what none of them arrive with
//! is where the line they are read against belongs.
//!
//! # Why this is not a coverage track
//!
//! [`CoverageTrack`](crate::CoverageTrack) draws upwards from the floor of its
//! band, which is right for a read depth: depth cannot go negative, and zero is
//! the floor in fact and not by convention. A selection statistic has neither
//! property. pN/pS is centred on one, and a half is as far below it as a two is
//! above; Tajima's D and GC skew are signed and cross zero wherever the thing
//! they measure changes direction. Drawn up from the bottom of a band, all of
//! those lose the one thing they were computed to say, which is which side of
//! the line a window fell on.
//!
//! So this track puts the line where the statistic says it belongs and lets a
//! window hang below it, in a colour of its own.
//!
//! # The constructor is where the statistic is decided
//!
//! Where the line goes and whether the values are transformed on the way in are
//! the two things easiest to get wrong here and hardest to catch afterwards, so
//! each named constructor settles them rather than leaving them to the caller.
//!
//! A ratio needs one more step than a signed statistic. On a linear axis a
//! pN/pS of 0.5 sits half a unit under the line and a 2.0 sits a whole unit
//! over it, so the same twofold departure looks twice as big in one direction
//! as the other. [`WindowTrack::ratios`] plots log2 of the ratio instead, which
//! puts the two at equal distances from the line, where they belong.
//!
//! A composition needs the line moved instead. [`WindowTrack::gc_content`] puts
//! it at the GC of the sequence it was handed rather than at a half;
//! [`WindowTrack::gc_skew`] leaves it at zero, where the sign of the statistic
//! already means what it says. [`WindowTrack::new`] takes the values as they
//! came, and then [`WindowTrack::baseline`] is yours to place.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Track};

/// One window and the value computed in it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Window {
    /// First base of the window, 0-based.
    pub start: u64,
    /// One past the last base of the window.
    pub end: u64,
    /// The statistic. A non-finite value is a window with no answer, usually
    /// because it held nothing to compute one from, and is left blank.
    pub value: f64,
}

impl Window {
    /// A window over `[start, end)` holding `value`.
    pub fn new(start: u64, end: u64, value: f64) -> Self {
        Window { start, end, value }
    }
}

/// How a windowed statistic is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowStyle {
    /// A block per window, from the baseline out to the value. The default,
    /// and the honest one: a window is an interval, and a block says so.
    Steps,
    /// A line through the windows, for a statistic read as a curve. GC skew is
    /// the usual case, where the shape between the windows is the point.
    Line,
}

/// A statistic in windows along the sequence, drawn either side of a baseline.
///
/// ```
/// use karyon::{Figure, Region, Window, WindowTrack};
///
/// // Selection along a gene, in 300 bp windows.
/// let windows: Vec<Window> = (0..10)
///     .map(|i| Window::new(i * 300, (i + 1) * 300, if i == 4 { 2.4 } else { 0.3 }))
///     .collect();
///
/// let svg = Figure::new(Region::new("Rv0667", 0, 3_000).unwrap())
///     .push(WindowTrack::ratios(windows).label("pN/pS"))
///     .to_svg();
/// assert!(svg.contains("pN/pS"));
/// ```
#[derive(Debug, Clone)]
pub struct WindowTrack {
    windows: Vec<Window>,
    label: Option<String>,
    height: f64,
    baseline: f64,
    style: WindowStyle,
    above_color: Option<String>,
    below_color: Option<String>,
    extent: Option<f64>,
    symmetric: bool,
    show_scale: bool,
    unit: String,
}

impl WindowTrack {
    /// A track over `windows`, with the baseline at zero.
    pub fn new(windows: impl Into<Vec<Window>>) -> Self {
        WindowTrack {
            windows: windows.into(),
            label: None,
            height: 56.0,
            baseline: 0.0,
            style: WindowStyle::Steps,
            above_color: None,
            below_color: None,
            extent: None,
            symmetric: true,
            show_scale: true,
            unit: String::new(),
        }
    }

    /// A track over ratios, plotted as log2 so that both directions read alike.
    ///
    /// Each value becomes `log2(value)`, so a pN/pS of 1 lands on the baseline,
    /// a 2 sits one unit above it and a 0.5 one unit below. A ratio of zero or
    /// less has no logarithm and is dropped, which is the right answer: a gene
    /// with no synonymous change has no ratio, not a ratio of zero.
    pub fn ratios(windows: impl Into<Vec<Window>>) -> Self {
        let logged: Vec<Window> = windows
            .into()
            .into_iter()
            .map(|window| Window {
                value: if window.value > 0.0 {
                    window.value.log2()
                } else {
                    f64::NAN
                },
                ..window
            })
            .collect();
        WindowTrack::new(logged).unit(" log2")
    }

    /// GC skew, `(G - C) / (G + C)`, in windows of `window` bases.
    ///
    /// The classic use is finding the origin and terminus of replication, where
    /// the cumulative skew turns over. That works on anything replicated in two
    /// directions from one origin: a chromosome, a plasmid, an organelle genome. This is the
    /// per-window skew rather than the cumulative one; sum it yourself if the
    /// turning point is what you are after.
    pub fn gc_skew(start: u64, seq: &[u8], window: u64) -> Self {
        WindowTrack::new(scan(start, seq, window, |counts| {
            let gc = counts.g + counts.c;
            if gc == 0 {
                f64::NAN
            } else {
                (counts.g as f64 - counts.c as f64) / gc as f64
            }
        }))
    }

    /// GC content as a fraction, in windows of `window` bases.
    ///
    /// The baseline is put at the GC content of the whole sequence handed in,
    /// so a window reads as richer or poorer than its own genome rather than
    /// against an arbitrary half. In a genome that sits at 65% GC, as
    /// *M. tuberculosis* does, a 50% line would put every window on one side of
    /// it and say nothing.
    pub fn gc_content(start: u64, seq: &[u8], window: u64) -> Self {
        let whole = counts_of(seq);
        let total = whole.g + whole.c + whole.a + whole.t;
        let mean = if total == 0 {
            0.5
        } else {
            (whole.g + whole.c) as f64 / total as f64
        };
        WindowTrack::new(scan(start, seq, window, |counts| {
            let total = counts.a + counts.c + counts.g + counts.t;
            if total == 0 {
                f64::NAN
            } else {
                (counts.g + counts.c) as f64 / total as f64
            }
        }))
        .baseline(mean)
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

    /// Moves the line the statistic is read against.
    ///
    /// Zero for a signed statistic, one for a raw ratio, the genome mean for a
    /// composition. It is the whole point of the track, so it is worth being
    /// deliberate about: a baseline in the wrong place turns "unusual" into
    /// "usual" without any warning that it has.
    pub fn baseline(mut self, baseline: f64) -> Self {
        self.baseline = baseline;
        self
    }

    /// Chooses blocks or a line.
    pub fn style(mut self, style: WindowStyle) -> Self {
        self.style = style;
        self
    }

    /// Sets the colours for windows above and below the baseline.
    pub fn colors(mut self, above: impl Into<String>, below: impl Into<String>) -> Self {
        self.above_color = Some(above.into());
        self.below_color = Some(below.into());
        self
    }

    /// Pins how far the axis reaches either side of the baseline.
    ///
    /// Pin it when putting two of these side by side, or the eye will read two
    /// different scales as one.
    pub fn extent(mut self, extent: f64) -> Self {
        self.extent = Some(extent.abs());
        self
    }

    /// Whether the axis reaches equally far above and below the baseline.
    ///
    /// On by default, because a departure of the same size should look the same
    /// size whichever way it went. Turn it off when the data are lopsided
    /// enough that half the band would be empty. Either way the baseline itself
    /// stays on the axis, since a band that does not show the line cannot say
    /// which side of it a window fell on.
    pub fn symmetric(mut self, symmetric: bool) -> Self {
        self.symmetric = symmetric;
        self
    }

    /// Shows or hides the value axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Sets a suffix printed after each axis number, such as `" log2"`.
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self
    }

    /// The windows.
    pub fn windows(&self) -> &[Window] {
        &self.windows
    }

    /// Where the baseline sits.
    ///
    /// Worth asking for after [`WindowTrack::gc_content`], which works it out
    /// from the sequence rather than taking it from the caller.
    pub fn baseline_value(&self) -> f64 {
        self.baseline
    }

    /// The values the axis spans, low first.
    ///
    /// Taken over every window handed in rather than only the visible ones, so
    /// that panning across a figure does not quietly rescale the axis under
    /// you. A statistic is read by comparing one stretch of sequence with
    /// another, and that comparison needs the axis to hold still.
    pub fn range(&self) -> (f64, f64) {
        if let Some(extent) = self.extent {
            return (self.baseline - extent, self.baseline + extent);
        }
        let mut lo = self.baseline;
        let mut hi = self.baseline;
        for window in &self.windows {
            if window.value.is_finite() {
                lo = lo.min(window.value);
                hi = hi.max(window.value);
            }
        }
        if self.symmetric {
            let extent = (hi - self.baseline).max(self.baseline - lo) * 1.06;
            let extent = if extent > 0.0 { extent } else { 1.0 };
            return (self.baseline - extent, self.baseline + extent);
        }
        let pad = (hi - lo) * 0.06;
        if hi - lo <= 0.0 {
            (self.baseline - 1.0, self.baseline + 1.0)
        } else {
            (lo - pad, hi + pad)
        }
    }

    /// Per pixel column, the lowest and highest value in it.
    ///
    /// Once the region is wide enough that several windows fall in one column,
    /// a statistic that swings both ways inside that column has to show both
    /// swings: reducing it to a single number would hide the very thing the
    /// sign is there to report.
    fn columns(&self, scale: &Scale, region: &Region, x0: f64, width: f64) -> Vec<Column> {
        let count = width.max(1.0).ceil() as usize;
        let mut columns = vec![Column::default(); count];
        for window in &self.windows {
            if !window.value.is_finite()
                || window.end <= region.start()
                || window.start >= region.end()
            {
                continue;
            }
            let left = scale.x(window.start).max(x0);
            let right = scale.x(window.end).min(x0 + width);
            if right < left {
                continue;
            }
            // Rounded rather than floored and ceiled, so that a column belongs
            // to the window covering most of it and two neighbours never share
            // one. Sharing would put a one pixel tick of the opposite colour at
            // every join, which reads as a pattern the data does not have.
            let first = ((left - x0).round().max(0.0) as usize).min(count - 1);
            let last = ((right - x0).round().max(0.0) as usize).clamp(first + 1, count);
            for column in &mut columns[first..last] {
                column.add(window.value);
            }
        }
        columns
    }
}

/// What one pixel column of the track holds.
#[derive(Debug, Clone, Copy)]
struct Column {
    lo: f64,
    hi: f64,
    sum: f64,
    count: usize,
}

impl Default for Column {
    fn default() -> Self {
        Column {
            lo: f64::INFINITY,
            hi: f64::NEG_INFINITY,
            sum: 0.0,
            count: 0,
        }
    }
}

impl Column {
    fn add(&mut self, value: f64) {
        self.lo = self.lo.min(value);
        self.hi = self.hi.max(value);
        self.sum += value;
        self.count += 1;
    }

    fn mean(&self) -> Option<f64> {
        (self.count > 0).then(|| self.sum / self.count as f64)
    }
}

impl Track for WindowTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.windows.is_empty() {
            return 0.0;
        }
        let (lo, hi) = self.range();
        let widest = [lo, hi, self.baseline]
            .iter()
            .map(|value| {
                text_width(
                    &format!("{}{}", trim(*value), self.unit),
                    theme.font_size - 1.0,
                )
            })
            .fold(0.0f64, f64::max);
        widest + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let (lo, hi) = self.range();
        let y_of = |value: f64| band.bottom() - ((value - lo) / (hi - lo)).clamp(0.0, 1.0) * band.h;
        let base_y = y_of(self.baseline);

        // The baseline is drawn first and always, even with no data, because it
        // is what the band means. An empty track that shows its line says "no
        // windows here"; an empty rectangle says nothing at all.
        ctx.svg
            .line(band.x, base_y, band.right(), base_y, &ctx.theme.rule, 1.0);

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            let mut label = |y: f64, value: f64| {
                ctx.svg.text(
                    right,
                    y,
                    &format!("{}{}", trim(value), self.unit),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            };
            label(band.y + size, hi);
            label(base_y + size * 0.35, self.baseline);
            // Lifted by a descender's worth: the band is what the track is
            // clipped to, and a baseline sitting exactly on it would have the
            // tail of a "g" cut off by the clip.
            label(band.bottom() - size * 0.22, lo);
        }

        if self.windows.is_empty() {
            return;
        }

        let above = self
            .above_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(0).to_string());
        let below = self
            .below_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());
        let columns = self.columns(ctx.scale, ctx.region, band.x, band.w);

        match self.style {
            WindowStyle::Steps => {
                let mut index = 0;
                while index < columns.len() {
                    let column = columns[index];
                    if column.count == 0 {
                        index += 1;
                        continue;
                    }
                    // Every column inside one window carries the same pair, so
                    // the run is emitted as a single block. Drawn one pixel at
                    // a time instead, each sliver lands on a fractional
                    // boundary and the renderer leaves a seam at every join,
                    // which reads as a bar chart of the wrong bars.
                    let mut run = index + 1;
                    while run < columns.len()
                        && columns[run].count > 0
                        && columns[run].lo == column.lo
                        && columns[run].hi == column.hi
                    {
                        run += 1;
                    }
                    let x = band.x + index as f64;
                    let width = (run - index) as f64;
                    if column.hi > self.baseline {
                        let top = y_of(column.hi);
                        ctx.svg.rect(x, top, width, base_y - top, &above);
                    }
                    if column.lo < self.baseline {
                        let bottom = y_of(column.lo);
                        ctx.svg.rect(x, base_y, width, bottom - base_y, &below);
                    }
                    index = run;
                }
            }
            WindowStyle::Line => {
                // A line cannot show both ends of a column that swings either
                // way, so it shows the mean of the column and the baseline
                // underneath it says which side that fell on.
                let points: Vec<(f64, f64, bool)> = columns
                    .iter()
                    .enumerate()
                    .filter_map(|(index, column)| {
                        column.mean().map(|mean| {
                            (
                                band.x + index as f64 + 0.5,
                                y_of(mean),
                                mean >= self.baseline,
                            )
                        })
                    })
                    .collect();
                // Broken into runs either side of the baseline, so the line
                // wears the same two colours the blocks do. One colour for the
                // whole trace would say every window was above the line.
                let mut run: Vec<(f64, f64)> = Vec::new();
                let mut side = points.first().map(|point| point.2).unwrap_or(true);
                for (x, y, up) in &points {
                    if *up != side && run.len() > 1 {
                        ctx.svg
                            .polyline(&run, if side { &above } else { &below }, 1.6);
                        run = vec![*run.last().expect("the run is not empty")];
                    }
                    side = *up;
                    run.push((*x, *y));
                }
                if run.len() > 1 {
                    ctx.svg
                        .polyline(&run, if side { &above } else { &below }, 1.6);
                }
            }
        }
    }
}

/// A, C, G and T in a stretch of sequence, everything else ignored.
#[derive(Debug, Clone, Copy, Default)]
struct Counts {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
}

fn counts_of(seq: &[u8]) -> Counts {
    let mut counts = Counts::default();
    for base in seq {
        match base.to_ascii_uppercase() {
            b'A' => counts.a += 1,
            b'C' => counts.c += 1,
            b'G' => counts.g += 1,
            b'T' | b'U' => counts.t += 1,
            _ => {}
        }
    }
    counts
}

/// Walks `seq` in windows and turns each one into a value.
///
/// The last window is whatever is left over rather than being padded or
/// dropped, so the track covers the sequence it was given exactly.
fn scan(start: u64, seq: &[u8], window: u64, value: impl Fn(Counts) -> f64) -> Vec<Window> {
    let step = window.max(1) as usize;
    seq.chunks(step)
        .enumerate()
        .map(|(index, chunk)| {
            let from = start + (index * step) as u64;
            Window::new(from, from + chunk.len() as u64, value(counts_of(chunk)))
        })
        .collect()
}

/// Two decimals at most, and none when they would be zeros.
fn trim(value: f64) -> String {
    if !value.is_finite() {
        return "0".to_string();
    }
    let rounded = (value * 100.0).round() / 100.0;
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

    fn region() -> Region {
        Region::new("chr1", 0, 3_000).unwrap()
    }

    fn windows() -> Vec<Window> {
        vec![
            Window::new(0, 1_000, 1.0),
            Window::new(1_000, 2_000, -0.5),
            Window::new(2_000, 3_000, 0.25),
        ]
    }

    #[test]
    fn a_signed_statistic_reaches_both_ways_from_the_baseline() {
        let track = WindowTrack::new(windows());
        let (lo, hi) = track.range();
        assert!(lo < 0.0 && hi > 0.0);
        // Symmetric by default: a departure of a given size looks the same size
        // whichever way it went.
        assert!((lo + hi).abs() < 1e-9, "{lo} and {hi} are not a pair");
    }

    #[test]
    fn turning_symmetry_off_gives_the_empty_half_of_the_band_back() {
        let lopsided = vec![Window::new(0, 100, 10.0), Window::new(100, 200, 9.0)];
        let (lo, hi) = WindowTrack::new(lopsided.clone()).symmetric(false).range();
        assert!(hi < 11.0 && lo > -1.0, "the axis stays near the data");
        // The baseline is still on it, whichever way symmetry is set: a band
        // that does not show the line cannot say which side of it anything is.
        assert!(lo <= 0.0 && hi >= 0.0);

        let (lo, hi) = WindowTrack::new(lopsided).range();
        assert!(hi > 10.0 && (lo + hi).abs() < 1e-9, "symmetric pairs up");
    }

    #[test]
    fn a_pinned_extent_is_taken_literally() {
        let track = WindowTrack::new(windows()).baseline(1.0).extent(3.0);
        assert_eq!(track.range(), (-2.0, 4.0));
    }

    #[test]
    fn a_flat_track_still_gets_an_axis_to_draw_on() {
        let flat = vec![Window::new(0, 100, 0.0)];
        let (lo, hi) = WindowTrack::new(flat).range();
        assert!(hi > lo, "a zero range would divide by zero");
    }

    #[test]
    fn ratios_become_logs_so_that_halving_matches_doubling() {
        let track = WindowTrack::ratios(vec![
            Window::new(0, 100, 2.0),
            Window::new(100, 200, 0.5),
            Window::new(200, 300, 1.0),
        ]);
        let values: Vec<f64> = track.windows().iter().map(|w| w.value).collect();
        assert_eq!(values, vec![1.0, -1.0, 0.0]);
    }

    #[test]
    fn a_ratio_with_no_denominator_is_dropped_rather_than_flattened() {
        // A gene with no synonymous change has no pN/pS, and calling that zero
        // would draw it as the most conserved thing on the chromosome.
        let track = WindowTrack::ratios(vec![Window::new(0, 100, 0.0)]);
        assert!(track.windows()[0].value.is_nan());
    }

    #[test]
    fn gc_skew_is_the_difference_over_the_sum() {
        let track = WindowTrack::gc_skew(0, b"GGGGCCCCGGGGGGGG", 8);
        let values: Vec<f64> = track.windows().iter().map(|w| w.value).collect();
        // Four G against four C is no skew; eight G against none is all of it.
        assert_eq!(values, vec![0.0, 1.0]);
    }

    #[test]
    fn a_window_with_nothing_in_it_has_no_answer() {
        let track = WindowTrack::gc_skew(0, b"AAAATTTT", 4);
        assert!(track.windows().iter().all(|w| w.value.is_nan()));
    }

    #[test]
    fn gc_content_is_read_against_the_genome_it_came_from() {
        // Two thirds GC overall, so that is where the line goes: against a half
        // every window of this sequence would sit on the same side.
        let track = WindowTrack::gc_content(0, b"GGCCAAGGCCAAGGCCAA", 6);
        assert!((track.baseline - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn the_last_window_is_the_remainder_rather_than_a_padded_one() {
        let track = WindowTrack::gc_skew(100, b"GGGGGGGGGG", 4);
        let spans: Vec<(u64, u64)> = track.windows().iter().map(|w| (w.start, w.end)).collect();
        assert_eq!(spans, vec![(100, 104), (104, 108), (108, 110)]);
    }

    #[test]
    fn windows_are_painted_on_the_side_they_fell() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(WindowTrack::new(windows()).colors("#111111", "#222222"))
            .to_svg();
        assert!(svg.contains("#111111"), "nothing above the line");
        assert!(svg.contains("#222222"), "nothing below the line");
    }

    #[test]
    fn a_column_holding_a_swing_both_ways_shows_both() {
        // Two thousand windows over nine hundred pixels puts several in every
        // column, and a column that swings both ways has to say so.
        let windows: Vec<Window> = (0..2_000)
            .map(|i| Window::new(i, i + 1, if i % 2 == 0 { 1.0 } else { -1.0 }))
            .collect();
        let track = WindowTrack::new(windows).colors("#111111", "#222222");
        let region = Region::new("chr1", 0, 2_000).unwrap();
        let scale = Scale::new(&region, 0.0, 900.0);
        let columns = track.columns(&scale, &region, 0.0, 900.0);
        assert!(columns.iter().any(|c| c.lo < 0.0 && c.hi > 0.0));

        let svg = Figure::new(region)
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("#111111") && svg.contains("#222222"));
    }

    #[test]
    fn one_window_is_one_block_rather_than_a_row_of_slivers() {
        // Three windows over a thousand bases is hundreds of pixels each. Drawn
        // a pixel at a time every join lands on a fractional boundary and the
        // renderer leaves a seam, so a window comes out as one rectangle.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(WindowTrack::new(windows()).colors("#111111", "#222222"))
            .to_svg();
        let blocks = svg.matches("#111111").count() + svg.matches("#222222").count();
        assert_eq!(blocks, 3, "one rectangle per window, not one per pixel");
    }

    #[test]
    fn an_empty_track_still_draws_its_baseline() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(WindowTrack::new(Vec::new()).label("pN/pS"))
            .to_svg();
        assert!(svg.contains("<line"), "the line is what the band means");
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_axis_is_only_reserved_when_it_is_wanted() {
        let theme = Theme::light();
        assert!(WindowTrack::new(windows()).y_axis_width(&theme) > 0.0);
        assert_eq!(
            WindowTrack::new(windows())
                .show_scale(false)
                .y_axis_width(&theme),
            0.0
        );
        assert_eq!(WindowTrack::new(Vec::new()).y_axis_width(&theme), 0.0);
    }

    #[test]
    fn windows_outside_the_region_are_not_drawn() {
        let track = WindowTrack::new(vec![Window::new(90_000, 91_000, 5.0)]);
        let region = region();
        let scale = Scale::new(&region, 0.0, 900.0);
        let columns = track.columns(&scale, &region, 0.0, 900.0);
        assert!(columns.iter().all(|c| c.count == 0));
    }

    #[test]
    fn a_line_style_draws_a_line() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(WindowTrack::new(windows()).style(WindowStyle::Line))
            .to_svg();
        assert!(svg.contains("<polyline"));
    }
}
