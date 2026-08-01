//! Reads that align to more than one place, drawn as one molecule.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Strand, Track};

/// One alignment of one read.
///
/// The reference span is where it landed; the read span is where those bases sat
/// on the molecule, which is what puts the segments in order.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitSegment {
    start: u64,
    end: u64,
    strand: Strand,
    read: Option<(u64, u64)>,
    mapq: Option<u8>,
}

impl SplitSegment {
    /// An alignment covering the reference span `start..end`, 0-based half-open.
    pub fn new(start: u64, end: u64, strand: Strand) -> Self {
        let (start, end) = if end >= start {
            (start, end)
        } else {
            (end, start)
        };
        SplitSegment {
            start,
            end,
            strand,
            read: None,
            mapq: None,
        }
    }

    /// Where on the read these bases sat, 0-based half-open.
    ///
    /// From a SAM `SA` tag this is the soft-clip either side of the CIGAR. It is
    /// optional: without it the segments are taken to run back to back in the
    /// order they were given, which is what a caller who has already sorted them
    /// means anyway.
    pub fn read_span(mut self, from: u64, to: u64) -> Self {
        self.read = Some(if to >= from { (from, to) } else { (to, from) });
        self
    }

    /// Mapping quality, drawn as a fading of the segment.
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = Some(mapq);
        self
    }

    /// Reference start, 0-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Reference end, exclusive.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Which strand it aligned to.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Reference position where the molecule enters this segment.
    ///
    /// The low coordinate on the forward strand and the high one on the reverse,
    /// because a reverse alignment is read backwards along the reference.
    pub fn entry(&self) -> u64 {
        if self.strand == Strand::Reverse {
            self.end
        } else {
            self.start
        }
    }

    /// Reference position where the molecule leaves this segment.
    pub fn exit(&self) -> u64 {
        if self.strand == Strand::Reverse {
            self.start
        } else {
            self.end
        }
    }
}

/// One molecule and every place it aligned, in read order.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitRead {
    name: Option<String>,
    segments: Vec<SplitSegment>,
}

impl SplitRead {
    /// A read from its alignments, which must already be in read order.
    pub fn new(segments: Vec<SplitSegment>) -> Self {
        SplitRead {
            name: None,
            segments,
        }
    }

    /// Sets the read name shown at the left of its row.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The alignments, in read order.
    pub fn segments(&self) -> &[SplitSegment] {
        &self.segments
    }

    /// Lowest and highest reference coordinate the molecule touched.
    pub fn extent(&self) -> Option<(u64, u64)> {
        let first = self.segments.first()?;
        Some(
            self.segments
                .iter()
                .fold((first.start, first.end), |(low, high), segment| {
                    (low.min(segment.start), high.max(segment.end))
                }),
        )
    }

    /// Whether consecutive segments ever move backwards along the reference.
    ///
    /// True for an inversion and for a transposon that landed upstream of where
    /// it came from, false for a plain deletion junction.
    pub fn goes_backwards(&self) -> bool {
        self.segments
            .windows(2)
            .any(|pair| pair[1].entry() < pair[0].exit())
    }
}

/// Reads that align in pieces, one row per molecule.
///
/// A [`PileupTrack`](crate::PileupTrack) read is one start, one CIGAR and one
/// strand, so a molecule that visits three places in the reference cannot be
/// written down in it at all. A [`StructuralTrack`](crate::StructuralTrack) arc
/// starts from a finished two-breakpoint call, so it cannot either: by the time
/// there is an arc, the evidence has been summarised away.
///
/// This is the evidence. One row is one molecule, one bar is one alignment, and
/// the connectors between the bars say in what order and in which orientation
/// that single piece of DNA visited those coordinates. A connector that runs
/// backwards is drawn under the row rather than over it, so a read that crosses
/// an inversion looks different from a read that crosses a deletion instead of
/// merely being annotated differently.
///
/// Colour ramps along the read from its 5' end, which is the other half of the
/// claim: it is what distinguishes a molecule that went A then B then C from one
/// that went C then B then A across the same three places.
///
/// A transposition is three segments and not an arc, which is the kind of figure
/// this exists for.
///
/// ```
/// use karyon::{Figure, Region, SplitRead, SplitReadTrack, SplitSegment, Strand};
///
/// let read = SplitRead::new(vec![
///     SplitSegment::new(1_000, 1_600, Strand::Forward),
///     SplitSegment::new(9_000, 9_400, Strand::Reverse),
///     SplitSegment::new(1_600, 2_100, Strand::Forward),
/// ])
/// .name("m64011_1");
/// assert!(read.goes_backwards());
///
/// let svg = Figure::new(Region::new("NC_000962.3", 0, 10_000).unwrap())
///     .push(SplitReadTrack::new(vec![read]).label("split reads"))
///     .to_svg();
/// assert!(svg.contains("<polygon"));
/// ```
#[derive(Debug, Clone)]
pub struct SplitReadTrack {
    reads: Vec<SplitRead>,
    row_height: f64,
    row_gap: f64,
    label: Option<String>,
    min_segment: f64,
    show_names: bool,
    ramp: bool,
    color: Option<String>,
}

impl SplitReadTrack {
    /// A track over these molecules, one row each in the order given.
    pub fn new(reads: Vec<SplitRead>) -> Self {
        SplitReadTrack {
            reads,
            row_height: 13.0,
            row_gap: 7.0,
            label: None,
            min_segment: 3.0,
            show_names: true,
            ramp: true,
            color: None,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one molecule's bar.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(3.0);
        self
    }

    /// Sets the gap between rows, which is where a backward connector runs.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(2.0);
        self
    }

    /// Sets the narrowest a segment may be drawn, in pixels.
    ///
    /// A 300 bp segment over a 4 Mb view is a twentieth of a pixel, and a
    /// segment that vanishes takes its junction with it. Three pixels by
    /// default: wide enough to see, narrow enough that it cannot be mistaken
    /// for a real span.
    pub fn min_segment(mut self, pixels: f64) -> Self {
        self.min_segment = pixels.max(0.5);
        self
    }

    /// Whether to print read names beside their rows.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Whether colour ramps along the read from its 5' end.
    ///
    /// On by default. Turning it off makes every segment of a read one colour,
    /// which loses the order and leaves the connectors to carry it alone.
    pub fn read_ramp(mut self, ramp: bool) -> Self {
        self.ramp = ramp;
        self
    }

    /// Overrides the base colour the ramp is built from.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// The molecules, in the order they are drawn.
    pub fn reads(&self) -> &[SplitRead] {
        &self.reads
    }

    /// How many segments fall entirely outside the view.
    pub fn not_in_view(&self, scale: &Scale) -> usize {
        let (from, to) = scale.bounds();
        self.reads
            .iter()
            .flat_map(|read| read.segments.iter())
            .filter(|segment| segment.end <= from || segment.start >= to)
            .count()
    }

    /// Where a segment sits along its read, from 0 at the 5' end to 1 at the 3'.
    ///
    /// Read spans are used when they were given, and the segment order stands in
    /// for them when they were not.
    fn along(&self, read: &SplitRead, index: usize) -> f64 {
        let count = read.segments.len();
        if count <= 1 {
            return 0.0;
        }
        let measured: Option<Vec<(u64, u64)>> =
            read.segments.iter().map(|segment| segment.read).collect();
        match measured {
            Some(spans) => {
                let total = spans.last().map(|span| span.1).unwrap_or(0) as f64;
                if total <= 0.0 {
                    return index as f64 / (count - 1) as f64;
                }
                let (from, to) = spans[index];
                ((from + to) as f64 / 2.0) / total
            }
            None => index as f64 / (count - 1) as f64,
        }
    }
}

impl Track for SplitReadTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        if self.reads.is_empty() {
            return self.row_height;
        }
        let rows = self.reads.len() as f64;
        rows * self.row_height + (rows - 1.0) * self.row_gap + 12.0
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_names {
            return 0.0;
        }
        let widest = self
            .reads
            .iter()
            .filter_map(|read| read.name.as_deref())
            .map(|name| text_width(name, theme.font_size))
            .fold(0.0f64, f64::max);
        if widest <= 0.0 {
            0.0
        } else {
            widest + 8.0
        }
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        if self.reads.is_empty() {
            return;
        }
        let band = ctx.band;
        let theme = ctx.theme;
        let base = self
            .color
            .clone()
            .unwrap_or_else(|| theme.color(0).to_string());
        let pitch = self.row_height + self.row_gap;

        for (row, read) in self.reads.iter().enumerate() {
            let top = band.y + row as f64 * pitch;
            if top > band.bottom() {
                break;
            }
            self.draw_read(ctx, read, top, &base);
        }

        // Off-screen segments are the one thing a reader cannot see for
        // themselves, so they are counted rather than quietly dropped.
        let hidden = self.not_in_view(ctx.scale);
        if hidden > 0 {
            let text = if hidden == 1 {
                "1 segment not in view".to_string()
            } else {
                format!("{hidden} segments not in view")
            };
            ctx.svg.text(
                band.right() - 2.0,
                band.bottom() - 2.0,
                &text,
                &theme.muted,
                theme.font_size,
                Anchor::End,
            );
        }
    }
}

impl SplitReadTrack {
    fn draw_read(&self, ctx: &mut DrawContext<'_>, read: &SplitRead, top: f64, base: &str) {
        let band = ctx.band;
        let theme = ctx.theme;
        let h = self.row_height;
        let mid = top + h / 2.0;

        // The connectors go down first, so a segment always paints over the
        // line that arrives at it rather than the other way round.
        for (index, pair) in read.segments.windows(2).enumerate() {
            let (from, to) = (pair[0].exit(), pair[1].entry());
            let x0 = ctx.scale.x(from);
            let x1 = ctx.scale.x(to);
            if x0.max(x1) < band.x - 4.0 || x0.min(x1) > band.right() + 4.0 {
                continue;
            }
            let backwards = to < from;
            // Forward junctions arch over the row and backward ones dip under
            // it, into the gap the row pitch reserves. The direction of travel
            // is then readable without following the colour ramp.
            let reach = (self.row_gap * 0.7).min(h * 0.9);
            let apex = if backwards {
                mid + h / 2.0 + reach
            } else {
                top - reach
            };
            let ink = if backwards {
                mix(base, "#000000", 0.35)
            } else {
                theme.muted.clone()
            };
            let d = format!(
                "M {} {} Q {} {} {} {}",
                num(x0),
                num(if backwards { mid + h / 2.0 } else { top }),
                num((x0 + x1) / 2.0),
                num(apex + (apex - mid) * 0.15),
                num(x1),
                num(if backwards { mid + h / 2.0 } else { top }),
            );
            ctx.svg
                .path_stroked(&d, &ink, if backwards { 1.4 } else { 1.0 });
            let _ = index;
        }

        for (index, segment) in read.segments.iter().enumerate() {
            let along = if self.ramp {
                self.along(read, index)
            } else {
                0.0
            };
            // Pale at the 5' end, saturated at the 3'. Which end a molecule
            // started at is the claim, so it gets the strongest signal the
            // track has after position itself.
            let fill = mix(&mix(base, theme.surface(), 0.62), base, along);
            let fill = match segment.mapq {
                Some(q) if q < 20 => mix(&fill, theme.surface(), 0.45),
                _ => fill,
            };
            self.draw_segment(ctx, segment, top, h, &fill);
        }

        // Names go in the axis strip rather than beside the leftmost segment.
        // Hung off the reads themselves they step in and out with wherever each
        // molecule happens to start, which reads as raggedness and not as data.
        if self.show_names {
            if let Some(name) = &read.name {
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    mid + theme.font_size * 0.35,
                    name,
                    &theme.muted,
                    theme.font_size,
                    Anchor::End,
                );
            }
        }
    }

    /// One alignment: a bar with a head on the end the molecule leaves by.
    fn draw_segment(
        &self,
        ctx: &mut DrawContext<'_>,
        segment: &SplitSegment,
        top: f64,
        h: f64,
        fill: &str,
    ) {
        let band = ctx.band;
        let x0 = ctx.scale.x(segment.start);
        let x1 = ctx.scale.x(segment.end);
        let (x0, x1) = if x1 - x0 < self.min_segment {
            let mid = (x0 + x1) / 2.0;
            (mid - self.min_segment / 2.0, mid + self.min_segment / 2.0)
        } else {
            (x0, x1)
        };
        if x1 < band.x - 2.0 || x0 > band.right() + 2.0 {
            return;
        }

        let head = ((x1 - x0) * 0.3).min(h * 0.8).min(x1 - x0).max(0.0);
        let bottom = top + h;
        let reverse = segment.strand == Strand::Reverse;
        let points: Vec<(f64, f64)> = if reverse {
            vec![
                (x1, top),
                (x0 + head, top),
                (x0, top + h / 2.0),
                (x0 + head, bottom),
                (x1, bottom),
            ]
        } else {
            vec![
                (x0, top),
                (x1 - head, top),
                (x1, top + h / 2.0),
                (x1 - head, bottom),
                (x0, bottom),
            ]
        };
        ctx.svg.polygon(&points, fill);
    }
}

/// The rectangle a row of this track occupies, for tests.
#[cfg(test)]
fn row_rect(track: &SplitReadTrack, band: crate::track::Rect, row: usize) -> crate::track::Rect {
    let pitch = track.row_height + track.row_gap;
    crate::track::Rect {
        x: band.x,
        y: band.y + row as f64 * pitch,
        w: band.w,
        h: track.row_height,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn scale(start: u64, end: u64) -> Scale {
        Scale::new(&Region::new("chr", start, end).unwrap(), 0.0, 900.0)
    }

    fn transposition() -> SplitRead {
        // A read that leaves the donor, crosses an IS copy and comes back:
        // three alignments, and the middle one on the other strand.
        SplitRead::new(vec![
            SplitSegment::new(1_000, 1_600, Strand::Forward).read_span(0, 600),
            SplitSegment::new(9_000, 9_400, Strand::Reverse).read_span(600, 1_000),
            SplitSegment::new(1_600, 2_100, Strand::Forward).read_span(1_000, 1_500),
        ])
        .name("m64011_1")
    }

    #[test]
    fn a_segment_is_entered_and_left_by_strand() {
        let forward = SplitSegment::new(100, 200, Strand::Forward);
        assert_eq!((forward.entry(), forward.exit()), (100, 200));
        let reverse = SplitSegment::new(100, 200, Strand::Reverse);
        assert_eq!(
            (reverse.entry(), reverse.exit()),
            (200, 100),
            "a reverse alignment is read backwards along the reference"
        );
    }

    #[test]
    fn a_backwards_hop_is_what_separates_an_inversion_from_a_deletion() {
        assert!(transposition().goes_backwards());

        // A plain deletion junction only ever moves forwards.
        let deletion = SplitRead::new(vec![
            SplitSegment::new(1_000, 1_500, Strand::Forward),
            SplitSegment::new(4_000, 4_500, Strand::Forward),
        ]);
        assert!(!deletion.goes_backwards());
    }

    #[test]
    fn the_extent_covers_every_place_the_molecule_touched() {
        assert_eq!(transposition().extent(), Some((1_000, 9_400)));
        assert_eq!(SplitRead::new(vec![]).extent(), None);
    }

    #[test]
    fn the_ramp_follows_read_position_when_it_is_given() {
        let track = SplitReadTrack::new(vec![transposition()]);
        let read = &track.reads()[0];
        let first = track.along(read, 0);
        let middle = track.along(read, 1);
        let last = track.along(read, 2);
        assert!(first < middle && middle < last, "{first} {middle} {last}");
        assert!(first >= 0.0 && last <= 1.0);
    }

    #[test]
    fn without_read_spans_the_order_given_is_the_order_used() {
        let read = SplitRead::new(vec![
            SplitSegment::new(500, 600, Strand::Forward),
            SplitSegment::new(100, 200, Strand::Forward),
            SplitSegment::new(900, 950, Strand::Forward),
        ]);
        let track = SplitReadTrack::new(vec![read]);
        let read = &track.reads()[0];
        assert_eq!(track.along(read, 0), 0.0);
        assert_eq!(track.along(read, 2), 1.0);
        assert!(track.along(read, 1) > 0.0);
    }

    #[test]
    fn one_segment_of_a_read_is_all_at_the_five_prime_end() {
        let read = SplitRead::new(vec![SplitSegment::new(0, 100, Strand::Forward)]);
        let track = SplitReadTrack::new(vec![read]);
        assert_eq!(track.along(&track.reads()[0], 0), 0.0);
    }

    #[test]
    fn segments_off_screen_are_counted_rather_than_dropped() {
        let track = SplitReadTrack::new(vec![transposition()]);
        assert_eq!(track.not_in_view(&scale(0, 10_000)), 0);
        assert_eq!(
            track.not_in_view(&scale(0, 3_000)),
            1,
            "the hop to 9 kb is out of view"
        );
        assert_eq!(track.not_in_view(&scale(500_000, 600_000)), 3);
    }

    #[test]
    fn rows_do_not_overlap_each_other() {
        let track = SplitReadTrack::new(vec![transposition(), transposition()]);
        let band = crate::track::Rect {
            x: 0.0,
            y: 0.0,
            w: 900.0,
            h: track.height(&scale(0, 10_000)),
        };
        let first = row_rect(&track, band, 0);
        let second = row_rect(&track, band, 1);
        assert!(second.y >= first.bottom(), "{second:?} under {first:?}");
        assert!(second.bottom() <= band.bottom());
    }

    #[test]
    fn a_short_segment_is_still_wide_enough_to_see() {
        // 200 bp over 4.4 Mb is a twentieth of a pixel, and a junction whose
        // segment vanished is a junction that did not happen.
        let read = SplitRead::new(vec![
            SplitSegment::new(1_000, 1_200, Strand::Forward),
            SplitSegment::new(3_000_000, 3_400_000, Strand::Forward),
        ]);
        let svg = Figure::new(Region::new("chr", 0, 4_400_000).unwrap())
            .width(900.0)
            .push(SplitReadTrack::new(vec![read]).min_segment(4.0))
            .to_svg();
        let widths: Vec<f64> = svg
            .match_indices("<polygon points=\"")
            .map(|(at, key)| {
                let rest = &svg[at + key.len()..];
                let head = &rest[..rest.find('"').unwrap()];
                let xs: Vec<f64> = head
                    .split(' ')
                    .filter_map(|pair| pair.split(',').next()?.parse().ok())
                    .collect();
                xs.iter().cloned().fold(f64::MIN, f64::max)
                    - xs.iter().cloned().fold(f64::MAX, f64::min)
            })
            .collect();
        assert!(!widths.is_empty());
        assert!(
            widths.iter().all(|w| *w >= 3.9),
            "every segment is drawn wide enough: {widths:?}"
        );
    }

    #[test]
    fn a_backward_connector_dips_below_its_row_and_a_forward_one_arches_over() {
        let track = SplitReadTrack::new(vec![transposition()]).row_height(14.0);
        let svg = Figure::new(Region::new("chr", 0, 10_000).unwrap())
            .push(track)
            .to_svg();
        let paths: Vec<&str> = svg
            .match_indices("<path d=\"")
            .map(|(at, key)| {
                let rest = &svg[at + key.len()..];
                &rest[..rest.find('"').unwrap()]
            })
            .collect();
        assert_eq!(paths.len(), 2, "one connector per junction: {paths:?}");

        let control_y = |d: &str| -> f64 {
            let after = d.split('Q').nth(1).unwrap().trim();
            after.split(' ').nth(1).unwrap().parse().unwrap()
        };
        let start_y = |d: &str| -> f64 {
            d.split('M')
                .nth(1)
                .unwrap()
                .trim()
                .split(' ')
                .nth(1)
                .unwrap()
                .parse()
                .unwrap()
        };
        // The hop out to 9 kb goes forwards, the hop back to 1,600 goes
        // backwards, so one control point is above its row and one below.
        assert!(
            control_y(paths[0]) < start_y(paths[0]),
            "forward arches over"
        );
        assert!(
            control_y(paths[1]) > start_y(paths[1]),
            "backward dips under"
        );
    }

    #[test]
    fn an_empty_track_draws_nothing_and_still_has_a_height() {
        let track = SplitReadTrack::new(vec![]);
        assert!(track.height(&scale(0, 100)) > 0.0);
        let svg = Figure::new(Region::new("chr", 0, 100).unwrap())
            .push(track)
            .to_svg();
        assert!(!svg.contains("<polygon"));
    }

    #[test]
    fn a_reverse_segment_points_the_other_way() {
        let forward = SplitReadTrack::new(vec![SplitRead::new(vec![SplitSegment::new(
            100,
            900,
            Strand::Forward,
        )])]);
        let reverse = SplitReadTrack::new(vec![SplitRead::new(vec![SplitSegment::new(
            100,
            900,
            Strand::Reverse,
        )])]);
        let render = |track: SplitReadTrack| {
            Figure::new(Region::new("chr", 0, 1_000).unwrap())
                .push(track)
                .to_svg()
        };
        assert_ne!(
            render(forward),
            render(reverse),
            "strand has to change the glyph, not just the colour"
        );
    }
}
