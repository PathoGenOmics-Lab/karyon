//! Raw nanopore current, before it was ever a base.
//!
//! A read starts life as a few hundred thousand current measurements, in
//! picoamperes, sampled thousands of times a second while the strand ratchets
//! through the pore. Basecalling turns that into letters and throws the rest
//! away. When a basecall is in doubt, or a modification is the thing being
//! measured, the current is the evidence and the letters are the summary: a
//! methylated cytosine does not change the base, it changes the current.
//!
//! The y axis is that current as it was recorded, which is comparable within a
//! read and not between reads. [`SquiggleTrack::normalized`] is what puts two
//! reads on one axis, at the price of an axis that is no longer in picoamperes.
//!
//! # The x axis is time, not position
//!
//! It counts samples, so put the figure over `Region::new("read", 0, samples)`.
//! Nothing in the band is measured in bases until [`SquiggleTrack::moves`]
//! attaches the basecaller's move table; with it the trace picks up a tint per
//! called base, letters where there is room for them, and the dwell times
//! [`SquiggleTrack::dwells`] reports.
//!
//! # Reading it at a distance
//!
//! A read is far longer than a figure is wide, so above one sample per pixel
//! each column is drawn as the range of the samples underneath it, which is how
//! an oscilloscope and an audio editor draw the same problem. The trace becomes
//! an envelope: the extremes are honest and the shape between them is not
//! there. Zoom in far enough and the samples are drawn as themselves.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Track};

/// Where a called base starts in the signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Move {
    /// Index of the first sample belonging to this base.
    pub sample: usize,
    /// The base that was called.
    pub base: u8,
}

impl Move {
    /// A base beginning at `sample`.
    pub fn new(sample: usize, base: u8) -> Self {
        Move { sample, base }
    }
}

/// Raw signal from one read.
///
/// ```
/// use karyon::{Figure, Region, SquiggleTrack};
///
/// let signal: Vec<f64> = (0..2_000)
///     .map(|i| 95.0 + 12.0 * ((i as f64) / 40.0).sin())
///     .collect();
///
/// let svg = Figure::new(Region::new("read", 0, 2_000).unwrap())
///     .push(SquiggleTrack::new(0, signal).label("pA"))
///     .to_svg();
/// assert!(svg.contains("<path") || svg.contains("<rect"));
/// ```
#[derive(Debug, Clone)]
pub struct SquiggleTrack {
    start: usize,
    signal: Vec<f64>,
    moves: Vec<Move>,
    label: Option<String>,
    height: f64,
    color: Option<String>,
    range: Option<(f64, f64)>,
    show_scale: bool,
    show_bases: bool,
    unit: String,
    point_threshold: f64,
}

impl SquiggleTrack {
    /// A track whose `signal[i]` is sample `start + i`.
    pub fn new(start: usize, signal: impl Into<Vec<f64>>) -> Self {
        SquiggleTrack {
            start,
            signal: signal.into(),
            moves: Vec::new(),
            label: None,
            height: 90.0,
            color: None,
            range: None,
            show_scale: true,
            show_bases: true,
            unit: " pA".to_string(),
            point_threshold: 4.0,
        }
    }

    /// A track over signal rescaled the way a nanopore tool rescales it.
    ///
    /// Subtracts the median and divides by the median absolute deviation, which
    /// is what puts two reads from two pores on one axis: the open pore current
    /// drifts between pores and over the life of a flow cell, so raw
    /// picoamperes from two reads are not comparable and normalised ones are.
    /// The unit becomes a count of deviations rather than a current.
    ///
    /// A signal with no spread at all is returned untouched, since dividing by
    /// its deviation would divide by zero.
    pub fn normalized(start: usize, signal: impl Into<Vec<f64>>) -> Self {
        let signal = signal.into();
        let median = median_of(&signal);
        let deviations: Vec<f64> = signal
            .iter()
            .filter(|value| value.is_finite())
            .map(|value| (value - median).abs())
            .collect();
        let mad = median_of(&deviations);
        let scaled = if mad > 0.0 {
            signal
                .iter()
                .map(|value| (value - median) / mad)
                .collect::<Vec<f64>>()
        } else {
            signal
        };
        SquiggleTrack::new(start, scaled).unit(" MAD")
    }

    /// Attaches the base boundaries, so the trace can be read against letters.
    ///
    /// One entry per called base, holding the first sample that belongs to it.
    /// This is what a basecaller's move table gives you, and it is the only
    /// thing in the plot that connects time to sequence.
    pub fn moves(mut self, moves: impl Into<Vec<Move>>) -> Self {
        self.moves = moves.into();
        self.moves.sort_by_key(|step| step.sample);
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(12.0);
        self
    }

    /// Overrides the trace colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Pins the current axis, instead of taking it from the signal.
    pub fn range(mut self, low: f64, high: f64) -> Self {
        self.range = Some((low.min(high), low.max(high)));
        self
    }

    /// Sets the unit printed after the axis numbers.
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = unit.into();
        self
    }

    /// Draws or hides the current axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Draws or hides the called bases over the trace.
    pub fn show_bases(mut self, show: bool) -> Self {
        self.show_bases = show;
        self
    }

    /// Sets how many pixels a sample needs before samples are drawn as points.
    pub fn point_threshold(mut self, pixels: f64) -> Self {
        self.point_threshold = pixels.max(0.0);
        self
    }

    /// The signal.
    pub fn signal(&self) -> &[f64] {
        &self.signal
    }

    /// The base boundaries.
    pub fn base_moves(&self) -> &[Move] {
        &self.moves
    }

    /// The sample at an index, if the track carries it.
    pub fn sample(&self, index: usize) -> Option<f64> {
        index
            .checked_sub(self.start)
            .and_then(|offset| self.signal.get(offset))
            .copied()
            .filter(|value| value.is_finite())
    }

    /// How many samples the read holds.
    pub fn len(&self) -> usize {
        self.signal.len()
    }

    /// Whether the read is empty.
    pub fn is_empty(&self) -> bool {
        self.signal.is_empty()
    }

    /// The current values the axis spans, low first.
    pub fn extent(&self) -> (f64, f64) {
        if let Some(range) = self.range {
            return range;
        }
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for value in self.signal.iter().filter(|value| value.is_finite()) {
            low = low.min(*value);
            high = high.max(*value);
        }
        if !low.is_finite() || high <= low {
            return (0.0, 1.0);
        }
        let pad = (high - low) * 0.06;
        (low - pad, high + pad)
    }

    /// How long each called base lasts, in samples.
    ///
    /// The dwell time is the measurement a homopolymer breaks: ten adenines in
    /// a row look like one long adenine, and the only thing that says there
    /// were ten is how long the current stayed put.
    pub fn dwells(&self) -> Vec<(u8, usize)> {
        let end = self.start + self.signal.len();
        self.moves
            .iter()
            .enumerate()
            .map(|(index, step)| {
                let next = self
                    .moves
                    .get(index + 1)
                    .map_or(end, |following| following.sample);
                (step.base, next.saturating_sub(step.sample))
            })
            .collect()
    }
}

impl Track for SquiggleTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.signal.is_empty() {
            return 0.0;
        }
        let (low, high) = self.extent();
        let widest = [low, high]
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
        if self.signal.is_empty() {
            return;
        }
        let (low, high) = self.extent();
        let y_of =
            |value: f64| band.bottom() - ((value - low) / (high - low)).clamp(0.0, 1.0) * band.h;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.accent.clone());

        if self.show_bases && !self.moves.is_empty() {
            self.draw_base_tints(ctx, band.h);
        }

        let per_px = ctx.scale.px_per_bp();
        if per_px >= self.point_threshold.max(0.01) {
            // Zoomed in far enough that a sample is a thing you can point at,
            // so it is drawn as one rather than as part of an envelope.
            let points: Vec<(f64, f64)> = self
                .visible(ctx.region)
                .filter_map(|index| {
                    self.sample(index)
                        .map(|value| (ctx.scale.x_center(index as u64), y_of(value)))
                })
                .collect();
            ctx.svg.polyline(&points, &color, 1.4);
            if per_px >= self.point_threshold * 2.5 {
                for (x, y) in &points {
                    ctx.svg
                        .circle_ringed(*x, *y, 1.8, &color, ctx.theme.surface(), 0.8);
                }
            }
        } else {
            // One bar per column from the lowest sample in it to the highest,
            // which is the envelope an oscilloscope draws and the only honest
            // summary of a column holding a hundred measurements.
            let columns = band.w.max(1.0).ceil() as usize;
            for column in 0..columns {
                let x = band.x + column as f64;
                let from = ctx.scale.pos_at_x(x).floor().max(0.0) as usize;
                let to = ctx.scale.pos_at_x(x + 1.0).ceil().max(0.0) as usize;
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for index in from..to.max(from + 1) {
                    if let Some(value) = self.sample(index) {
                        lo = lo.min(value);
                        hi = hi.max(value);
                    }
                }
                if !lo.is_finite() {
                    continue;
                }
                let top = y_of(hi);
                ctx.svg.rect(x, top, 1.0, (y_of(lo) - top).max(0.8), &color);
            }
        }

        // The letters go over the trace rather than under it: the current
        // climbs into the top of the band wherever the base is a T, and a
        // letter behind the line it labels is a letter nobody can read.
        if self.show_bases && !self.moves.is_empty() {
            self.draw_base_letters(ctx);
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            ctx.svg.text(
                right,
                band.y + size,
                &format!("{}{}", trim(high), self.unit),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
            ctx.svg.text(
                right,
                band.bottom() - size * 0.22,
                &format!("{}{}", trim(low), self.unit),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }
    }
}

impl SquiggleTrack {
    /// Sample indices inside the region, and inside the read.
    fn visible(&self, region: &Region) -> impl Iterator<Item = usize> + '_ {
        let first = (region.start() as usize).max(self.start);
        let last = (region.end() as usize).min(self.start + self.signal.len());
        (first..last.max(first)).collect::<Vec<_>>().into_iter()
    }

    /// Alternating tints behind the samples of each base.
    fn draw_base_tints(&self, ctx: &mut DrawContext<'_>, height: f64) {
        let band = ctx.band;
        let tint = crate::theme::mix(ctx.theme.surface(), &ctx.theme.rule, 0.22);
        for (index, (x0, x1, _)) in self.base_spans(ctx).into_iter().enumerate() {
            // Every other base, so the boundaries read as boundaries rather
            // than as a wash of colour over the whole trace.
            if index % 2 == 0 {
                ctx.svg.rect(x0, band.y, x1 - x0, height, &tint);
            }
        }
    }

    /// The called base over each of its samples, where there is room.
    fn draw_base_letters(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let size = ctx.theme.font_size - 1.0;
        for (x0, x1, base) in self.base_spans(ctx) {
            if x1 - x0 >= size * 1.2 {
                ctx.svg.text(
                    (x0 + x1) / 2.0,
                    band.y + size,
                    &(base as char).to_ascii_uppercase().to_string(),
                    &ctx.theme.foreground,
                    size,
                    Anchor::Middle,
                );
            }
        }
    }

    /// Where each called base begins and ends on the page.
    fn base_spans(&self, ctx: &DrawContext<'_>) -> Vec<(f64, f64, u8)> {
        let band = ctx.band;
        let end = self.start + self.signal.len();
        self.moves
            .iter()
            .enumerate()
            .filter_map(|(index, step)| {
                let next = self
                    .moves
                    .get(index + 1)
                    .map_or(end, |following| following.sample);
                if next <= step.sample {
                    return None;
                }
                let x0 = ctx.scale.x(step.sample as u64);
                let x1 = ctx.scale.x(next as u64);
                (x1 >= band.x && x0 <= band.right()).then_some((x0, x1, step.base))
            })
            .collect()
    }
}

/// The middle value, or zero when there is nothing to take a middle of.
fn median_of(values: &[f64]) -> f64 {
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if sorted.is_empty() {
        return 0.0;
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN survives the filter"));
    let middle = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    }
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

    fn signal() -> Vec<f64> {
        (0..1_000)
            .map(|i| 95.0 + 12.0 * ((i as f64) / 30.0).sin())
            .collect()
    }

    fn region(samples: u64) -> Region {
        Region::new("read", 0, samples).unwrap()
    }

    #[test]
    fn a_sample_is_looked_up_by_its_index_in_the_read() {
        let track = SquiggleTrack::new(500, vec![1.0, 2.0, 3.0]);
        assert_eq!(track.sample(500), Some(1.0));
        assert_eq!(track.sample(502), Some(3.0));
        assert_eq!(track.sample(499), None, "before the first sample held");
        assert_eq!(track.sample(503), None, "past the last");
        assert_eq!(track.len(), 3);
    }

    #[test]
    fn a_missing_sample_is_missing_rather_than_zero() {
        let track = SquiggleTrack::new(0, vec![1.0, f64::NAN, 3.0]);
        assert_eq!(track.sample(1), None);
        let (low, high) = track.extent();
        assert!(low < 1.0 && high > 3.0, "and does not poison the axis");
    }

    #[test]
    fn normalising_puts_two_pores_on_one_axis() {
        // The same shape at two different open pore currents comes out the
        // same, which is the whole point: raw picoamperes are not comparable
        // between pores or across the life of a flow cell.
        let one: Vec<f64> = vec![90.0, 100.0, 110.0, 100.0];
        let two: Vec<f64> = one.iter().map(|value| value + 40.0).collect();
        let a = SquiggleTrack::normalized(0, one);
        let b = SquiggleTrack::normalized(0, two);
        assert_eq!(a.signal(), b.signal());
        // Centred on zero, in deviations rather than picoamperes.
        assert_eq!(a.signal()[1], 0.0);
        assert!(a.unit.contains("MAD"));
    }

    #[test]
    fn a_flat_signal_survives_normalising_rather_than_dividing_by_zero() {
        let track = SquiggleTrack::normalized(0, vec![100.0; 20]);
        assert!(track.signal().iter().all(|value| value.is_finite()));
        assert_eq!(track.signal()[0], 100.0, "returned untouched");
    }

    #[test]
    fn the_median_is_the_middle_either_way_round() {
        assert_eq!(median_of(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median_of(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median_of(&[]), 0.0);
        assert_eq!(median_of(&[f64::NAN]), 0.0);
    }

    #[test]
    fn dwell_times_are_the_gap_to_the_next_base() {
        let track = SquiggleTrack::new(0, vec![0.0; 100]).moves(vec![
            Move::new(0, b'A'),
            Move::new(10, b'C'),
            Move::new(40, b'G'),
        ]);
        // Ten, thirty, and then the rest of the read.
        assert_eq!(track.dwells(), vec![(b'A', 10), (b'C', 30), (b'G', 60)]);
    }

    #[test]
    fn moves_are_sorted_however_they_arrive() {
        let track = SquiggleTrack::new(0, vec![0.0; 100])
            .moves(vec![Move::new(40, b'G'), Move::new(0, b'A')]);
        assert_eq!(track.base_moves()[0].sample, 0);
    }

    #[test]
    fn a_long_read_is_drawn_as_an_envelope() {
        // A hundred thousand samples across nine hundred pixels is a hundred
        // measurements a column, and the honest summary of a column is the
        // range of what is in it.
        let long: Vec<f64> = (0..100_000)
            .map(|i| 95.0 + 12.0 * ((i as f64) / 7.0).sin())
            .collect();
        let svg = Figure::new(region(100_000))
            .show_region_label(false)
            .push(SquiggleTrack::new(0, long))
            .to_svg();
        assert!(!svg.contains("<polyline"), "not one point per sample");
        assert!(svg.matches("<rect").count() > 100, "one bar per column");
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn a_short_read_is_drawn_as_itself() {
        let svg = Figure::new(region(80))
            .show_region_label(false)
            .push(SquiggleTrack::new(0, signal()[..80].to_vec()))
            .to_svg();
        assert!(svg.contains("<polyline"), "the trace is a trace");
        assert!(svg.contains("<circle"), "and the samples are visible");
    }

    #[test]
    fn base_boundaries_tint_alternate_bases_and_letter_them() {
        let moves: Vec<Move> = (0..10).map(|i| Move::new(i * 8, b"ACGT"[i % 4])).collect();
        let svg = Figure::new(region(80))
            .show_region_label(false)
            .push(SquiggleTrack::new(0, signal()[..80].to_vec()).moves(moves))
            .to_svg();
        assert!(svg.contains(">A</text>"));
        assert!(svg.contains(">C</text>"));
        // Every other one, so the boundaries read as boundaries.
        assert_eq!(svg.matches(">A</text>").count(), 3);
    }

    #[test]
    fn bases_can_be_turned_off() {
        let moves = vec![Move::new(0, b'A'), Move::new(20, b'C')];
        let svg = Figure::new(region(80))
            .show_region_label(false)
            .push(
                SquiggleTrack::new(0, signal()[..80].to_vec())
                    .moves(moves)
                    .show_bases(false),
            )
            .to_svg();
        assert!(!svg.contains(">A</text>"));
    }

    #[test]
    fn the_axis_is_only_reserved_when_it_is_wanted() {
        let theme = Theme::light();
        let track = SquiggleTrack::new(0, signal());
        assert!(track.y_axis_width(&theme) > 0.0);
        assert_eq!(track.clone().show_scale(false).y_axis_width(&theme), 0.0);
        assert_eq!(SquiggleTrack::new(0, Vec::new()).y_axis_width(&theme), 0.0);
    }

    #[test]
    fn a_pinned_range_is_taken_literally() {
        let track = SquiggleTrack::new(0, signal()).range(0.0, 200.0);
        assert_eq!(track.extent(), (0.0, 200.0));
    }

    #[test]
    fn an_empty_read_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(region(10))
            .show_region_label(false)
            .push(SquiggleTrack::new(0, Vec::new()).label("pA"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
