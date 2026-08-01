//! The six reading frames, with their stops and their open stretches.
//!
//! # The shape
//!
//! Six lanes: three above a line for the frames read left to right, three below
//! for the frames read right to left off the other strand. A stop codon is a
//! tick across its lane and an open reading frame is the bar between two of
//! them. Nothing outside molecular biology looks like this, and it is what you
//! look at to decide whether an unannotated stretch of sequence is coding, and
//! which way round.
//!
//! # What counts as open
//!
//! A run of codons with no stop in it, at least [`OrfTrack::min_codons`] long.
//! Whether it starts at a methionine is a separate question and a separate
//! switch: [`OrfTrack::require_start`] asks for one. `ATG`, `GTG` and `TTG` all
//! count, since alternative starts are ordinary in bacteria, archaea and
//! organelle genomes, and the default is off because the first thing you want
//! from a six frame map is where the stops are not.
//!
//! # The frames are not the strands
//!
//! Frame `-1` is read from the far end of the sequence given, not from the far
//! end of the chromosome. Hand this the sequence of the region on display and
//! the reverse frames are the reverse frames of that region, which is what a
//! reader assumes and what every ORF finder does.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, wash, Theme};
use crate::track::{DrawContext, Track};

/// One open reading frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Orf {
    /// First base, 0-based, in the coordinates the track was given.
    pub start: u64,
    /// One past the last base.
    pub end: u64,
    /// Which frame it is in: 1, 2 or 3 forward, -1, -2 or -3 reverse.
    pub frame: i8,
}

impl Orf {
    /// How many codons long it is.
    pub fn codons(&self) -> u64 {
        (self.end - self.start) / 3
    }

    /// Whether it is read off the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.frame < 0
    }
}

/// The six reading frames of a stretch of sequence.
///
/// ```
/// use karyon::{Figure, OrfTrack, Region};
///
/// // A stretch with a long open frame in it.
/// let mut seq = b"ATG".to_vec();
/// seq.extend(std::iter::repeat(b'A').take(300));
/// seq.extend(b"TAA");
///
/// let track = OrfTrack::new(0, seq).min_codons(20);
/// assert!(!track.orfs().is_empty());
///
/// let svg = Figure::new(Region::new("plasmid", 0, 306).unwrap())
///     .push(track.label("frames"))
///     .to_svg();
/// assert!(svg.contains("<rect"));
/// ```
#[derive(Debug, Clone)]
pub struct OrfTrack {
    start: u64,
    seq: Vec<u8>,
    label: Option<String>,
    lane_height: f64,
    lane_gap: f64,
    min_codons: u64,
    require_start: bool,
    forward_color: Option<String>,
    reverse_color: Option<String>,
    show_frames: bool,
}

impl OrfTrack {
    /// A track over `seq`, whose first base is at `start`.
    pub fn new(start: u64, seq: impl Into<Vec<u8>>) -> Self {
        OrfTrack {
            start,
            seq: seq.into(),
            label: None,
            lane_height: 9.0,
            lane_gap: 2.0,
            min_codons: 30,
            require_start: false,
            forward_color: None,
            reverse_color: None,
            show_frames: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one frame lane.
    pub fn lane_height(mut self, height: f64) -> Self {
        self.lane_height = height.max(2.0);
        self
    }

    /// Sets the gap between lanes.
    pub fn lane_gap(mut self, gap: f64) -> Self {
        self.lane_gap = gap.max(0.0);
        self
    }

    /// Sets how many codons a stretch needs before it counts as open.
    ///
    /// Thirty by default. In random sequence a stop turns up every twenty-one
    /// codons or so, so a shorter floor fills all six lanes with stretches that
    /// are open by chance and say nothing.
    pub fn min_codons(mut self, codons: u64) -> Self {
        self.min_codons = codons.max(1);
        self
    }

    /// Requires an open frame to begin at a start codon.
    ///
    /// `ATG`, `GTG` and `TTG`. `ATG` is the common one everywhere, and the other
    /// two begin real genes often enough in bacteria, archaea and organelle
    /// genomes to be worth counting. Off by default:
    /// the first thing a six frame map is read for is where the stops are not,
    /// and insisting on a start hides a frame that is open but whose beginning
    /// is off the left of the view.
    pub fn require_start(mut self, require: bool) -> Self {
        self.require_start = require;
        self
    }

    /// Sets the colours of the forward and reverse lanes.
    pub fn colors(mut self, forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        self.forward_color = Some(forward.into());
        self.reverse_color = Some(reverse.into());
        self
    }

    /// Draws or hides the frame numbers in the left of each lane.
    pub fn show_frames(mut self, show: bool) -> Self {
        self.show_frames = show;
        self
    }

    /// The sequence.
    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    /// Every stop codon, as its position and its frame.
    pub fn stops(&self) -> Vec<(u64, i8)> {
        let mut stops = Vec::new();
        for offset in 0..3usize {
            for (position, codon) in self.codons(offset, false) {
                if is_stop(&codon) {
                    stops.push((position, (offset + 1) as i8));
                }
            }
            for (position, codon) in self.codons(offset, true) {
                if is_stop(&codon) {
                    stops.push((position, -((offset + 1) as i8)));
                }
            }
        }
        stops.sort_unstable();
        stops
    }

    /// Every open reading frame, in every frame.
    pub fn orfs(&self) -> Vec<Orf> {
        let mut orfs = Vec::new();
        for offset in 0..3usize {
            for reverse in [false, true] {
                let frame = if reverse {
                    -((offset + 1) as i8)
                } else {
                    (offset + 1) as i8
                };
                orfs.extend(self.frame_orfs(offset, reverse, frame));
            }
        }
        orfs.sort_by_key(|orf| (orf.start, orf.frame));
        orfs
    }

    /// The open stretches of one frame.
    fn frame_orfs(&self, offset: usize, reverse: bool, frame: i8) -> Vec<Orf> {
        let codons = self.codons(offset, reverse);
        let mut orfs = Vec::new();
        let mut open: Option<u64> = None;
        let mut started = !self.require_start;

        for (position, codon) in &codons {
            if is_stop(codon) {
                if let Some(from) = open.take() {
                    if started {
                        let (start, end) = span(from, *position, reverse);
                        push_orf(&mut orfs, start, end, frame, self.min_codons);
                    }
                }
                started = !self.require_start;
                continue;
            }
            if open.is_none() && (!self.require_start || is_start(codon)) {
                open = Some(*position);
                started = true;
            }
        }
        // A frame still open at the end of the sequence is still an answer: the
        // gene runs off the edge of what was handed in.
        if let (Some(from), true) = (open, started) {
            let last = codons.last().map_or(from, |(position, _)| *position);
            let (start, end) = span(from, last + 3, reverse);
            push_orf(&mut orfs, start, end, frame, self.min_codons);
        }
        orfs
    }

    /// The codons of one frame, each as the position of its first base in the
    /// coordinates the track was given.
    ///
    /// A reverse frame is walked from the far end, so its codons come back in
    /// descending order and each position is where the codon starts on the
    /// forward strand.
    fn codons(&self, offset: usize, reverse: bool) -> Vec<(u64, [u8; 3])> {
        let len = self.seq.len();
        let mut out = Vec::new();
        if len < offset + 3 {
            return out;
        }
        if reverse {
            let mut end = len - offset;
            while end >= 3 {
                let at = end - 3;
                let codon = [
                    complement(self.seq[end - 1]),
                    complement(self.seq[end - 2]),
                    complement(self.seq[at]),
                ];
                out.push((self.start + at as u64, codon));
                end -= 3;
            }
        } else {
            let mut at = offset;
            while at + 3 <= len {
                let codon = [self.seq[at], self.seq[at + 1], self.seq[at + 2]];
                out.push((self.start + at as u64, codon));
                at += 3;
            }
        }
        out
    }

    /// Where a frame's lane sits inside the band, top first.
    fn lane(&self, frame: i8, band_y: f64) -> f64 {
        // Frames run +3, +2, +1 down to the midline and then -1, -2, -3 below
        // it, so the two strands mirror each other about the line.
        let index = match frame {
            3 => 0,
            2 => 1,
            1 => 2,
            -1 => 3,
            -2 => 4,
            _ => 5,
        } as f64;
        band_y + index * (self.lane_height + self.lane_gap)
    }
}

/// A codon that ends translation.
fn is_stop(codon: &[u8; 3]) -> bool {
    matches!(
        [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase()
        ],
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
    )
}

/// A codon a gene may begin at, in the genomes that use more than `ATG`.
fn is_start(codon: &[u8; 3]) -> bool {
    matches!(
        [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase()
        ],
        [b'A', b'T', b'G'] | [b'G', b'T', b'G'] | [b'T', b'T', b'G']
    )
}

/// The complement of a base, leaving anything unexpected alone.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        other => other,
    }
}

/// The forward-strand span of a stretch walked in either direction.
fn span(from: u64, to: u64, reverse: bool) -> (u64, u64) {
    if reverse {
        (to.min(from), from.max(to) + 3)
    } else {
        (from, to)
    }
}

fn push_orf(orfs: &mut Vec<Orf>, start: u64, end: u64, frame: i8, min_codons: u64) {
    if end > start && (end - start) / 3 >= min_codons {
        orfs.push(Orf { start, end, frame });
    }
}

impl Track for OrfTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        6.0 * self.lane_height + 5.0 * self.lane_gap
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_frames {
            return 0.0;
        }
        let size = (theme.font_size - 3.0).min(self.lane_height);
        text_width("+3", size) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let forward = self
            .forward_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(0).to_string());
        let reverse = self
            .reverse_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());
        let rule = mix(ctx.theme.surface(), &ctx.theme.rule, 0.75);

        // A lane of its own for every frame, empty or not: six lanes is the
        // shape, and hiding the empty ones would say there were fewer frames.
        for frame in [3i8, 2, 1, -1, -2, -3] {
            let top = self.lane(frame, band.y);
            ctx.svg.rect(
                band.x,
                top + self.lane_height / 2.0 - 0.5,
                band.w,
                1.0,
                &rule,
            );
        }

        // Stops first, so an open frame is drawn between them rather than over
        // them: the ticks are what bound it.
        for (position, frame) in self.stops() {
            let top = self.lane(frame, band.y);
            let x = ctx.scale.x_center(position + 1);
            let color = if frame < 0 { &reverse } else { &forward };
            ctx.svg.rect(
                x - 0.5,
                top,
                1.0,
                self.lane_height,
                &mix(ctx.theme.surface(), color, 0.55),
            );
        }

        for orf in self.orfs() {
            let top = self.lane(orf.frame, band.y);
            let x0 = ctx.scale.x(orf.start);
            let x1 = ctx.scale.x(orf.end).max(x0 + 1.0);
            let color = if orf.is_reverse() { &reverse } else { &forward };
            ctx.svg.rect_rounded(
                x0,
                top,
                x1 - x0,
                self.lane_height,
                ctx.theme.corner_radius,
                &wash(color, ctx.theme),
            );
            ctx.svg.rect_outline(
                x0 + 0.4,
                top + 0.4,
                (x1 - x0 - 0.8).max(0.2),
                (self.lane_height - 0.8).max(0.2),
                color,
                0.9,
            );
        }

        if self.show_frames && ctx.axis.w > 0.0 {
            let size = (ctx.theme.font_size - 3.0).min(self.lane_height);
            for frame in [3i8, 2, 1, -1, -2, -3] {
                let top = self.lane(frame, band.y);
                let text = if frame > 0 {
                    format!("+{frame}")
                } else {
                    format!("{frame}")
                };
                // In the strip rather than in the plotting area. Written at the
                // left of the band they sat on top of whatever open frame
                // reached the edge of the view, which is most of them.
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    top + self.lane_height / 2.0 + size * 0.35,
                    &text,
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    /// A start, ninety-nine codons of alanine, then a stop, in frame +1.
    fn coding() -> Vec<u8> {
        let mut seq = b"ATG".to_vec();
        seq.extend(std::iter::repeat(b'G').take(297));
        seq.extend(b"TAA");
        seq
    }

    fn region(len: u64) -> Region {
        Region::new("plasmid", 0, len).unwrap()
    }

    #[test]
    fn a_stop_is_found_in_the_frame_it_is_in() {
        // TAA at offset 300, which is frame +1 since 300 is a multiple of three.
        let track = OrfTrack::new(0, coding());
        let stops = track.stops();
        assert!(stops.contains(&(300, 1)), "{stops:?}");
        assert!(!stops.iter().any(|(_, frame)| *frame == 2));
    }

    #[test]
    fn every_stop_codon_counts() {
        for stop in [b"TAA", b"TAG", b"TGA"] {
            let mut seq = b"AAA".to_vec();
            seq.extend(stop);
            let track = OrfTrack::new(0, seq);
            assert!(
                track.stops().iter().any(|(position, _)| *position == 3),
                "{} was not a stop",
                String::from_utf8_lossy(stop)
            );
        }
        // And a codon that is not one is not one.
        assert!(OrfTrack::new(0, b"AAATTT".to_vec())
            .stops()
            .iter()
            .all(|(_, frame)| *frame != 1));
    }

    #[test]
    fn an_open_frame_is_the_stretch_between_two_stops() {
        let track = OrfTrack::new(0, coding()).min_codons(50);
        let orfs = track.orfs();
        let plus_one: Vec<&Orf> = orfs.iter().filter(|orf| orf.frame == 1).collect();
        assert_eq!(plus_one.len(), 1, "{orfs:?}");
        assert_eq!(plus_one[0].start, 0);
        assert_eq!(plus_one[0].end, 300);
        assert_eq!(plus_one[0].codons(), 100);
        assert!(!plus_one[0].is_reverse());
    }

    #[test]
    fn a_short_stretch_is_not_an_open_frame() {
        // In random sequence a stop turns up every twenty-one codons or so, so
        // a low floor fills all six lanes with stretches open by chance.
        let track = OrfTrack::new(0, coding()).min_codons(200);
        assert!(track.orfs().iter().all(|orf| orf.frame != 1));
        assert!(!OrfTrack::new(0, coding()).min_codons(50).orfs().is_empty());
    }

    #[test]
    fn requiring_a_start_changes_where_a_frame_opens() {
        // GGG for twenty codons, then ATG, then more, then a stop.
        let mut seq = std::iter::repeat(b'G').take(60).collect::<Vec<u8>>();
        seq.extend(b"ATG");
        seq.extend(std::iter::repeat(b'C').take(150));
        seq.extend(b"TGA");

        let loose = OrfTrack::new(0, seq.clone()).min_codons(10);
        let strict = OrfTrack::new(0, seq).min_codons(10).require_start(true);
        let first = |track: &OrfTrack| {
            track
                .orfs()
                .into_iter()
                .find(|orf| orf.frame == 1)
                .map(|orf| orf.start)
        };
        assert_eq!(first(&loose), Some(0), "opens where the sequence does");
        assert_eq!(first(&strict), Some(60), "opens at the methionine");
    }

    #[test]
    fn every_bacterial_start_codon_counts() {
        for start in [b"ATG", b"GTG", b"TTG"] {
            let mut seq = start.to_vec();
            seq.extend(std::iter::repeat(b'C').take(150));
            seq.extend(b"TAA");
            let track = OrfTrack::new(0, seq).min_codons(10).require_start(true);
            assert!(
                track
                    .orfs()
                    .iter()
                    .any(|orf| orf.frame == 1 && orf.start == 0),
                "{} did not open a frame",
                String::from_utf8_lossy(start)
            );
        }
    }

    #[test]
    fn the_reverse_frames_are_read_off_the_other_strand() {
        // TTA on the forward strand is TAA read backwards on the reverse one.
        let mut seq = b"TTA".to_vec();
        seq.extend(std::iter::repeat(b'C').take(300));
        let track = OrfTrack::new(0, seq);
        assert!(
            track.stops().iter().any(|(_, frame)| *frame < 0),
            "no stop on the reverse strand"
        );
        // And nothing on the forward one, since TTA and CCC are not stops.
        assert!(track.stops().iter().all(|(_, frame)| *frame < 0));
    }

    #[test]
    fn a_frame_still_open_at_the_end_is_still_an_answer() {
        // No stop at all: the gene runs off the edge of what was handed in.
        let seq: Vec<u8> = std::iter::repeat(b'A').take(600).collect();
        let track = OrfTrack::new(0, seq).min_codons(50);
        assert!(track.orfs().iter().any(|orf| orf.frame == 1));
    }

    #[test]
    fn the_track_is_six_lanes_whether_or_not_they_hold_anything() {
        let scale = Scale::new(&region(300), 0.0, 800.0);
        let track = OrfTrack::new(0, coding()).lane_height(10.0).lane_gap(2.0);
        assert_eq!(track.height(&scale), 6.0 * 10.0 + 5.0 * 2.0);

        let svg = Figure::new(region(306))
            .show_region_label(false)
            .push(OrfTrack::new(0, coding()))
            .to_svg();
        for frame in ["+1", "+2", "+3", "-1", "-2", "-3"] {
            assert!(svg.contains(&format!(">{frame}</text>")), "missing {frame}");
        }
    }

    #[test]
    fn the_frame_numbers_take_a_strip_rather_than_the_plot() {
        // Written at the left of the band they sat on top of whatever open
        // frame reached the edge of the view, which is most of them.
        let theme = Theme::light();
        let track = OrfTrack::new(0, coding());
        assert!(track.y_axis_width(&theme) > 0.0);
        assert_eq!(track.clone().show_frames(false).y_axis_width(&theme), 0.0);
    }

    #[test]
    fn the_two_strands_mirror_each_other_about_the_line() {
        let track = OrfTrack::new(0, coding());
        let top = |frame: i8| track.lane(frame, 0.0);
        // Forward frames count down to the line, reverse frames away from it.
        assert!(top(3) < top(2) && top(2) < top(1));
        assert!(top(1) < top(-1));
        assert!(top(-1) < top(-2) && top(-2) < top(-3));
    }

    #[test]
    fn an_empty_sequence_draws_the_lanes_and_nothing_else() {
        let svg = Figure::new(region(100))
            .show_region_label(false)
            .push(OrfTrack::new(0, Vec::new()).label("frames"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
        assert!(OrfTrack::new(0, Vec::new()).orfs().is_empty());
        assert!(OrfTrack::new(0, b"AT".to_vec()).stops().is_empty());
    }
}
