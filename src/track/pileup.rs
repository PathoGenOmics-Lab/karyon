//! Aligned reads, stacked the way a genome browser stacks them.
//!
//! This is the track you open when a variant call looks wrong. It draws the
//! reads themselves: where they start, which strand they came off, where they
//! disagree with the reference, and where the aligner had to open a gap. Most
//! of the work is in the CIGAR, because a read is not an interval. It is a
//! sequence of operations against the reference, and only some of them advance
//! along it.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::Theme;
use crate::track::feature::Strand;
use crate::track::{DrawContext, Rect, Track};

/// One CIGAR operation.
///
/// The names follow the SAM specification, and so does what each one consumes.
/// A conversion from a BAM record is one arm per letter:
///
/// | SAM | Here | Reference | Query |
/// |:----|:-----|:----------|:------|
/// | `M`, `=`, `X` | [`Match`](CigarOp::Match) | yes | yes |
/// | `I` | [`Insertion`](CigarOp::Insertion) | no | yes |
/// | `D` | [`Deletion`](CigarOp::Deletion) | yes | no |
/// | `N` | [`Skip`](CigarOp::Skip) | yes | no |
/// | `S` | [`SoftClip`](CigarOp::SoftClip) | no | yes |
/// | `H` | [`HardClip`](CigarOp::HardClip) | no | no |
///
/// `P` has no entry because it draws as nothing; drop it on the way in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// Aligned bases, whether or not they agree with the reference. `M`, `=`
    /// and `X` all land here; the track compares the sequences itself rather
    /// than trusting the operation.
    Match(u32),
    /// Bases in the read that are not in the reference.
    Insertion(u32),
    /// Bases in the reference that are not in the read.
    Deletion(u32),
    /// A reference gap the read jumps, as a spliced alignment does.
    Skip(u32),
    /// Bases kept in the record but not aligned.
    SoftClip(u32),
    /// Bases removed from the record entirely.
    HardClip(u32),
}

impl CigarOp {
    /// How far along the reference this operation moves.
    pub fn reference_len(self) -> u64 {
        match self {
            CigarOp::Match(n) | CigarOp::Deletion(n) | CigarOp::Skip(n) => n as u64,
            _ => 0,
        }
    }

    /// How far along the read this operation moves.
    pub fn query_len(self) -> u64 {
        match self {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::SoftClip(n) => n as u64,
            _ => 0,
        }
    }
}

/// A piece of a read, already placed on the reference.
///
/// Produced by [`Read::segments`], which is where the CIGAR stops being a list
/// of operations and starts being geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Segment {
    /// Bases aligned to the reference, starting at `query` in the read.
    Aligned {
        /// First reference position, 0-based.
        start: u64,
        /// How many bases.
        len: u64,
        /// Offset of the first base in the read sequence.
        query: usize,
    },
    /// Reference bases the read does not have.
    Deletion {
        /// First reference position, 0-based.
        start: u64,
        /// How many bases.
        len: u64,
    },
    /// A reference gap the read jumps over.
    Skip {
        /// First reference position, 0-based.
        start: u64,
        /// How many bases.
        len: u64,
    },
    /// Read bases with nowhere to go, sitting between two reference positions.
    Insertion {
        /// Reference position the insertion sits before.
        at: u64,
        /// How many bases were inserted.
        len: u32,
    },
}

/// One aligned read.
///
/// Positions are 0-based, so a SAM `POS` becomes `POS - 1`. The sequence is
/// optional: without it the read still draws, it just cannot show mismatches.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    /// Leftmost aligned reference position, 0-based.
    pub start: u64,
    /// The alignment, operation by operation.
    pub cigar: Vec<CigarOp>,
    /// Query bases, including any soft clipped ones, as they appear in the
    /// record. Empty when the sequence is not to hand.
    pub sequence: Vec<u8>,
    /// Which strand the read came off.
    pub strand: Strand,
    /// Mapping quality, used to fade unreliable reads when asked.
    pub mapping_quality: Option<u8>,
}

impl Read {
    /// A read aligned at `start` with the given CIGAR.
    pub fn new(start: u64, cigar: impl Into<Vec<CigarOp>>) -> Self {
        Read {
            start,
            cigar: cigar.into(),
            sequence: Vec::new(),
            strand: Strand::Unknown,
            mapping_quality: None,
        }
    }

    /// A read that is one ungapped block, the common case.
    pub fn aligned(start: u64, len: u32) -> Self {
        Read::new(start, vec![CigarOp::Match(len)])
    }

    /// Attaches the query sequence, which is what mismatches are read from.
    pub fn sequence(mut self, sequence: impl Into<Vec<u8>>) -> Self {
        self.sequence = sequence.into();
        self
    }

    /// Sets the strand.
    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Sets the mapping quality.
    pub fn mapping_quality(mut self, quality: u8) -> Self {
        self.mapping_quality = Some(quality);
        self
    }

    /// How much reference the read covers, gaps and skips included.
    pub fn reference_span(&self) -> u64 {
        self.cigar.iter().map(|op| op.reference_len()).sum()
    }

    /// One past the last reference position the read touches.
    pub fn end(&self) -> u64 {
        self.start + self.reference_span()
    }

    /// Whether the read touches a region.
    pub fn overlaps(&self, region: &Region) -> bool {
        self.start < region.end() && self.end() > region.start()
    }

    /// Walks the CIGAR and places every piece on the reference.
    ///
    /// Soft and hard clips produce nothing: they advance the read but the track
    /// has nowhere to draw them, so they are consumed silently rather than
    /// shifting everything downstream.
    pub fn segments(&self) -> Vec<Segment> {
        let mut segments = Vec::with_capacity(self.cigar.len());
        let mut reference = self.start;
        let mut query = 0usize;

        for op in &self.cigar {
            match *op {
                CigarOp::Match(len) => {
                    if len > 0 {
                        segments.push(Segment::Aligned {
                            start: reference,
                            len: len as u64,
                            query,
                        });
                    }
                    reference += len as u64;
                    query += len as usize;
                }
                CigarOp::Deletion(len) => {
                    if len > 0 {
                        segments.push(Segment::Deletion {
                            start: reference,
                            len: len as u64,
                        });
                    }
                    reference += len as u64;
                }
                CigarOp::Skip(len) => {
                    if len > 0 {
                        segments.push(Segment::Skip {
                            start: reference,
                            len: len as u64,
                        });
                    }
                    reference += len as u64;
                }
                CigarOp::Insertion(len) => {
                    if len > 0 {
                        segments.push(Segment::Insertion { at: reference, len });
                    }
                    query += len as usize;
                }
                CigarOp::SoftClip(len) => query += len as usize,
                CigarOp::HardClip(_) => {}
            }
        }

        segments
    }

    /// The read base at a reference position, if the read has one there.
    ///
    /// `None` covers three different things: the position is outside the read,
    /// it falls in a deletion or a skip, or no sequence was attached.
    pub fn base_at(&self, position: u64) -> Option<u8> {
        if self.sequence.is_empty() {
            return None;
        }
        for segment in self.segments() {
            if let Segment::Aligned { start, len, query } = segment {
                if position >= start && position < start + len {
                    let offset = query + (position - start) as usize;
                    return self.sequence.get(offset).copied();
                }
            }
        }
        None
    }
}

/// What decides the colour of a read body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadColoring {
    /// One colour for every read, the quietest background for mismatches to
    /// stand out against. This is what IGV does by default and why.
    Uniform,
    /// Forward and reverse reads take different colours, which is how strand
    /// bias becomes visible at a glance.
    Strand,
}

/// Aligned reads, packed into rows.
///
/// ```
/// use karyon::{CigarOp, Figure, PileupTrack, Read, Region, Strand};
///
/// let reads = vec![
///     Read::new(100, vec![CigarOp::Match(50)])
///         .sequence(vec![b'A'; 50])
///         .strand(Strand::Forward),
///     Read::new(120, vec![CigarOp::Match(20), CigarOp::Deletion(5), CigarOp::Match(25)])
///         .strand(Strand::Reverse),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 90, 220).unwrap())
///     .push(PileupTrack::new(reads).reference(90, vec![b'C'; 130]).label("reads"))
///     .to_svg();
/// assert!(svg.contains("<rect"));
/// ```
#[derive(Debug, Clone)]
pub struct PileupTrack {
    reads: Vec<Read>,
    reference_start: u64,
    reference: Vec<u8>,
    label: Option<String>,
    read_height: f64,
    row_gap: f64,
    max_rows: Option<usize>,
    coloring: ReadColoring,
    color: Option<String>,
    forward_color: Option<String>,
    reverse_color: Option<String>,
    fade_by_quality: bool,
    show_mismatches: bool,
    mismatch_threshold: f64,
    letter_threshold: f64,
}

impl PileupTrack {
    /// A pileup of `reads`.
    pub fn new(reads: impl Into<Vec<Read>>) -> Self {
        PileupTrack {
            reads: reads.into(),
            reference_start: 0,
            reference: Vec::new(),
            label: None,
            read_height: 9.0,
            row_gap: 2.0,
            max_rows: Some(40),
            coloring: ReadColoring::Uniform,
            color: None,
            forward_color: None,
            reverse_color: None,
            fade_by_quality: false,
            show_mismatches: true,
            mismatch_threshold: 0.2,
            letter_threshold: 7.0,
        }
    }

    /// Attaches the reference, without which no mismatch can be found.
    ///
    /// `reference[i]` is the base at 0-based position `start + i`, the same
    /// shape [`SequenceTrack`](crate::SequenceTrack) takes.
    pub fn reference(mut self, start: u64, sequence: impl Into<Vec<u8>>) -> Self {
        self.reference_start = start;
        self.reference = sequence.into();
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of a single read in pixels.
    pub fn read_height(mut self, height: f64) -> Self {
        self.read_height = height.max(1.0);
        self
    }

    /// Sets the vertical gap between rows of reads.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Caps how many rows are drawn, or removes the cap with `None`.
    ///
    /// A pileup at a thousandfold depth is a thousand rows tall and useful to
    /// nobody, so the default stops at forty. Reads that do not fit are counted
    /// and reported on the figure rather than disappearing quietly.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|r| r.max(1));
        self
    }

    /// Chooses what decides a read's colour.
    pub fn coloring(mut self, coloring: ReadColoring) -> Self {
        self.coloring = coloring;
        self
    }

    /// Sets the colour of a read under [`ReadColoring::Uniform`].
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets the two colours used under [`ReadColoring::Strand`].
    pub fn strand_colors(mut self, forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        self.forward_color = Some(forward.into());
        self.reverse_color = Some(reverse.into());
        self
    }

    /// Fades reads in proportion to their mapping quality.
    ///
    /// A read with quality zero could have come from anywhere, and drawing it
    /// as solidly as a uniquely placed one is how a repeat becomes a variant
    /// call. Reads with no quality attached are drawn solid.
    pub fn fade_by_quality(mut self, fade: bool) -> Self {
        self.fade_by_quality = fade;
        self
    }

    /// Draws or hides mismatched bases.
    pub fn show_mismatches(mut self, show: bool) -> Self {
        self.show_mismatches = show;
        self
    }

    /// Sets how many pixels a base needs before mismatches are looked for.
    ///
    /// Below this a mismatch would be a sliver too thin to see, and finding it
    /// would still mean walking every base of every read, so the search is
    /// skipped entirely.
    pub fn mismatch_threshold(mut self, pixels: f64) -> Self {
        self.mismatch_threshold = pixels;
        self
    }

    /// The reads in the track.
    pub fn reads(&self) -> &[Read] {
        &self.reads
    }

    /// The reference base at a position, if the track was given one.
    pub fn reference_base(&self, position: u64) -> Option<u8> {
        if position < self.reference_start {
            return None;
        }
        self.reference
            .get((position - self.reference_start) as usize)
            .copied()
    }

    /// Packs the reads overlapping `region` into rows.
    ///
    /// Returns the row of each packed read, paired with its index in
    /// [`PileupTrack::reads`], plus how many reads did not fit under the row
    /// cap. Only reads on screen take part, so the packing is as tight as the
    /// view allows.
    pub fn layout(&self, scale: &Scale, region: &Region) -> PileupLayout {
        let mut visible: Vec<usize> = (0..self.reads.len())
            .filter(|i| self.reads[*i].overlaps(region))
            .collect();
        visible.sort_by_key(|i| (self.reads[*i].start, self.reads[*i].end()));

        // Two reads sharing a row need a visible gap between them, measured in
        // pixels so it survives any zoom.
        let padding = 3.0;
        let mut row_ends: Vec<f64> = Vec::new();
        let mut placed: Vec<(usize, usize)> = Vec::with_capacity(visible.len());
        let mut hidden = 0usize;

        for index in visible {
            let read = &self.reads[index];
            let left = scale.x(read.start);
            let right = scale.x(read.end()).max(left + 1.0);

            let row = match row_ends.iter().position(|end| *end + padding <= left) {
                Some(row) => {
                    row_ends[row] = right;
                    row
                }
                None => {
                    row_ends.push(right);
                    row_ends.len() - 1
                }
            };

            match self.max_rows {
                Some(cap) if row >= cap => hidden += 1,
                _ => placed.push((index, row)),
            }
        }

        let rows = match self.max_rows {
            Some(cap) => row_ends.len().min(cap),
            None => row_ends.len(),
        };

        PileupLayout {
            placed,
            rows: rows.max(1),
            hidden,
        }
    }

    fn read_color(&self, read: &Read, theme: &Theme) -> String {
        match self.coloring {
            ReadColoring::Uniform => self.color.clone().unwrap_or_else(|| "#b3bcc6".to_string()),
            ReadColoring::Strand => match read.strand {
                Strand::Reverse => self
                    .reverse_color
                    .clone()
                    .unwrap_or_else(|| theme.color(0).to_string()),
                _ => self
                    .forward_color
                    .clone()
                    .unwrap_or_else(|| theme.color(1).to_string()),
            },
        }
    }

    fn opacity(&self, read: &Read) -> f64 {
        if !self.fade_by_quality {
            return 1.0;
        }
        match read.mapping_quality {
            // Zero means "this could be anywhere"; 30 and up is as good as it
            // gets from most aligners, so that is where the ramp tops out.
            Some(quality) => 0.25 + 0.75 * (quality as f64 / 30.0).min(1.0),
            None => 1.0,
        }
    }
}

/// The result of packing a pileup into rows.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupLayout {
    /// Each drawn read as its index in the track paired with its row.
    pub placed: Vec<(usize, usize)>,
    /// How many rows the pileup occupies, at least one.
    pub rows: usize,
    /// Reads that overlapped the region but fell past the row cap.
    pub hidden: usize,
}

impl Track for PileupTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        // Height cannot depend on the region, which the trait does not offer
        // here, so the cap is what fixes it. Without a cap the track sizes
        // itself to the deepest possible stack instead.
        let rows = match self.max_rows {
            Some(cap) => cap,
            None => self.reads.len().max(1),
        };
        let rows = rows.min(self.reads.len().max(1));
        rows as f64 * self.read_height + rows.saturating_sub(1) as f64 * self.row_gap
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let layout = self.layout(ctx.scale, ctx.region);
        let px_per_bp = ctx.scale.px_per_bp();
        let hunt_mismatches = self.show_mismatches
            && !self.reference.is_empty()
            && px_per_bp >= self.mismatch_threshold;
        let draw_letters = px_per_bp >= self.letter_threshold;

        for (index, row) in &layout.placed {
            let read = &self.reads[*index];
            let top = band.y + *row as f64 * (self.read_height + self.row_gap);
            if top + self.read_height > band.bottom() + 0.5 {
                continue;
            }
            let color = self.read_color(read, ctx.theme);
            let opacity = self.opacity(read);
            let middle = top + self.read_height / 2.0;

            for segment in read.segments() {
                match segment {
                    Segment::Aligned { start, len, .. } => {
                        let left = ctx.scale.x(start);
                        let block = Rect {
                            x: left,
                            y: top,
                            w: ctx.scale.x(start + len) - left,
                            h: self.read_height,
                        };
                        self.draw_body(ctx, read, block, &color, opacity);
                    }
                    Segment::Deletion { start, len } => {
                        // A gap in the read: a line across, so the reads either
                        // side still read as one alignment.
                        ctx.svg.line(
                            ctx.scale.x(start),
                            middle,
                            ctx.scale.x(start + len),
                            middle,
                            &color,
                            1.5,
                        );
                    }
                    Segment::Skip { start, len } => {
                        ctx.svg.line(
                            ctx.scale.x(start),
                            middle,
                            ctx.scale.x(start + len),
                            middle,
                            &ctx.theme.rule,
                            1.0,
                        );
                    }
                    Segment::Insertion { at, .. } => {
                        let x = ctx.scale.x(at);
                        ctx.svg
                            .rect(x - 0.75, top, 1.5, self.read_height, "#8e44ad");
                    }
                }
            }

            if hunt_mismatches {
                self.draw_mismatches(ctx, read, top, draw_letters);
            }
        }

        if layout.hidden > 0 {
            let text = format!("+{} reads not shown", layout.hidden);
            let size = ctx.theme.font_size - 1.0;
            let width = text_width(&text, size) + 6.0;
            ctx.svg.rect_opacity(
                band.right() - width,
                band.bottom() - size - 3.0,
                width,
                size + 3.0,
                &ctx.theme.background,
                0.8,
            );
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 2.0,
                &text,
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }
    }
}

impl PileupTrack {
    /// Draws one aligned block, with a point on the end the read runs towards.
    fn draw_body(
        &self,
        ctx: &mut DrawContext<'_>,
        read: &Read,
        block: Rect,
        color: &str,
        opacity: f64,
    ) {
        let (left, right, top) = (block.x, block.right(), block.y);
        let width = right - left;
        if width <= 0.0 {
            return;
        }
        let bottom = top + self.read_height;
        let middle = top + self.read_height / 2.0;
        // Only the outermost block gets the arrowhead, and only when there is
        // room for one that still leaves a body behind it.
        let head = (width * 0.4).min(4.0);
        let pointed = width > 5.0
            && match read.strand {
                Strand::Forward => (right - ctx.scale.x(read.end())).abs() < 0.01,
                Strand::Reverse => (left - ctx.scale.x(read.start)).abs() < 0.01,
                Strand::Unknown => false,
            };

        if !pointed {
            ctx.svg
                .rect_opacity(left, top, width.max(0.6), self.read_height, color, opacity);
            return;
        }

        let points = match read.strand {
            Strand::Forward => vec![
                (left, top),
                (right - head, top),
                (right, middle),
                (right - head, bottom),
                (left, bottom),
            ],
            _ => vec![
                (right, top),
                (left + head, top),
                (left, middle),
                (left + head, bottom),
                (right, bottom),
            ],
        };
        if opacity < 1.0 {
            // The polygon writer takes no opacity, so a faded read is drawn as
            // a plain block rather than silently coming out solid.
            ctx.svg
                .rect_opacity(left, top, width, self.read_height, color, opacity);
        } else {
            ctx.svg.polygon(&points, color);
        }
    }

    /// Paints the bases where a read disagrees with the reference.
    fn draw_mismatches(
        &self,
        ctx: &mut DrawContext<'_>,
        read: &Read,
        top: f64,
        draw_letters: bool,
    ) {
        if read.sequence.is_empty() {
            return;
        }
        let first = ctx.region.start();
        let last = ctx.region.end();

        for segment in read.segments() {
            let Segment::Aligned { start, len, query } = segment else {
                continue;
            };
            let from = start.max(first);
            let to = (start + len).min(last);
            for position in from..to {
                let offset = query + (position - start) as usize;
                let Some(base) = read.sequence.get(offset).copied() else {
                    continue;
                };
                let Some(reference) = self.reference_base(position) else {
                    continue;
                };
                if base.eq_ignore_ascii_case(&reference) {
                    continue;
                }

                let x = ctx.scale.x(position);
                let width = ctx.scale.x(position + 1) - x;
                let color = ctx.theme.bases.of(base).to_string();
                // A mismatch is the reason the track exists, so it never gets
                // to be thinner than a pixel. The floor only bites at the wide
                // end of the zoom, where a base is sub-pixel anyway.
                ctx.svg
                    .rect(x, top, (width + 0.2).max(1.0), self.read_height, &color);

                if draw_letters && self.read_height >= 7.0 {
                    let size = (width * 0.75).min(self.read_height * 0.9);
                    ctx.svg.glyph(
                        x,
                        top + self.read_height - (self.read_height - size) / 2.0,
                        width,
                        size / ctx.theme.cap_height_ratio.max(0.1),
                        &(base as char).to_ascii_uppercase().to_string(),
                        "#ffffff",
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;

    fn region() -> Region {
        Region::new("chr1", 0, 200).unwrap()
    }

    fn scale(region: &Region) -> Scale {
        Scale::new(region, 0.0, 800.0)
    }

    #[test]
    fn a_cigar_knows_what_it_consumes() {
        assert_eq!(CigarOp::Match(10).reference_len(), 10);
        assert_eq!(CigarOp::Match(10).query_len(), 10);
        assert_eq!(CigarOp::Insertion(3).reference_len(), 0);
        assert_eq!(CigarOp::Insertion(3).query_len(), 3);
        assert_eq!(CigarOp::Deletion(4).reference_len(), 4);
        assert_eq!(CigarOp::Deletion(4).query_len(), 0);
        assert_eq!(CigarOp::Skip(100).reference_len(), 100);
        assert_eq!(CigarOp::SoftClip(5).reference_len(), 0);
        assert_eq!(CigarOp::SoftClip(5).query_len(), 5);
        assert_eq!(CigarOp::HardClip(5).reference_len(), 0);
        assert_eq!(CigarOp::HardClip(5).query_len(), 0);
    }

    #[test]
    fn the_reference_span_counts_only_what_advances_along_it() {
        let read = Read::new(
            100,
            vec![
                CigarOp::SoftClip(5),
                CigarOp::Match(20),
                CigarOp::Insertion(3),
                CigarOp::Match(10),
                CigarOp::Deletion(4),
                CigarOp::Match(10),
            ],
        );
        assert_eq!(read.reference_span(), 44);
        assert_eq!(read.end(), 144);
    }

    #[test]
    fn segments_place_every_piece_on_the_reference() {
        let read = Read::new(
            100,
            vec![
                CigarOp::SoftClip(5),
                CigarOp::Match(10),
                CigarOp::Deletion(2),
                CigarOp::Match(5),
            ],
        );
        assert_eq!(
            read.segments(),
            vec![
                // The soft clip advances the query but places nothing, so the
                // first aligned base is read offset 5 at reference 100.
                Segment::Aligned {
                    start: 100,
                    len: 10,
                    query: 5
                },
                Segment::Deletion { start: 110, len: 2 },
                Segment::Aligned {
                    start: 112,
                    len: 5,
                    query: 15
                },
            ]
        );
    }

    #[test]
    fn an_insertion_advances_the_query_but_not_the_reference() {
        let read = Read::new(
            50,
            vec![CigarOp::Match(4), CigarOp::Insertion(3), CigarOp::Match(4)],
        );
        assert_eq!(
            read.segments(),
            vec![
                Segment::Aligned {
                    start: 50,
                    len: 4,
                    query: 0
                },
                Segment::Insertion { at: 54, len: 3 },
                Segment::Aligned {
                    start: 54,
                    len: 4,
                    query: 7
                },
            ]
        );
        assert_eq!(read.reference_span(), 8);
    }

    #[test]
    fn a_spliced_read_jumps_the_intron() {
        let read = Read::new(
            1000,
            vec![CigarOp::Match(50), CigarOp::Skip(2000), CigarOp::Match(50)],
        );
        assert_eq!(read.reference_span(), 2100);
        assert!(matches!(
            read.segments()[1],
            Segment::Skip { len: 2000, .. }
        ));
    }

    #[test]
    fn zero_length_operations_produce_no_segments() {
        let read = Read::new(
            10,
            vec![CigarOp::Match(0), CigarOp::Deletion(0), CigarOp::Match(5)],
        );
        assert_eq!(read.segments().len(), 1);
    }

    #[test]
    fn base_at_walks_the_cigar_to_find_the_right_offset() {
        // Query ACGT TTT GGGG, aligned as 4M 3I 4M from reference 100.
        let read = Read::new(
            100,
            vec![CigarOp::Match(4), CigarOp::Insertion(3), CigarOp::Match(4)],
        )
        .sequence(b"ACGTTTTGGGG".to_vec());

        assert_eq!(read.base_at(100), Some(b'A'));
        assert_eq!(read.base_at(103), Some(b'T'));
        // Reference 104 is the first base after the insertion, read offset 7.
        assert_eq!(read.base_at(104), Some(b'G'));
        assert_eq!(read.base_at(107), Some(b'G'));
        assert_eq!(read.base_at(108), None);
        assert_eq!(read.base_at(99), None);
    }

    #[test]
    fn a_deleted_position_has_no_base() {
        let read = Read::new(
            100,
            vec![CigarOp::Match(4), CigarOp::Deletion(2), CigarOp::Match(4)],
        )
        .sequence(b"ACGTACGT".to_vec());
        assert_eq!(read.base_at(103), Some(b'T'));
        assert_eq!(read.base_at(104), None, "inside the deletion");
        assert_eq!(read.base_at(105), None, "inside the deletion");
        assert_eq!(read.base_at(106), Some(b'A'));
    }

    #[test]
    fn a_read_without_a_sequence_has_no_bases() {
        let read = Read::aligned(100, 10);
        assert_eq!(read.base_at(105), None);
    }

    #[test]
    fn overlapping_reads_take_separate_rows() {
        let reads = vec![
            Read::aligned(0, 50),
            Read::aligned(10, 50),
            Read::aligned(20, 50),
        ];
        let track = PileupTrack::new(reads);
        let layout = track.layout(&scale(&region()), &region());
        assert_eq!(layout.rows, 3);
        assert_eq!(layout.hidden, 0);
        let rows: Vec<usize> = layout.placed.iter().map(|(_, row)| *row).collect();
        assert_eq!(rows, vec![0, 1, 2]);
    }

    #[test]
    fn disjoint_reads_share_a_row() {
        let reads = vec![
            Read::aligned(0, 20),
            Read::aligned(60, 20),
            Read::aligned(120, 20),
        ];
        let track = PileupTrack::new(reads);
        let layout = track.layout(&scale(&region()), &region());
        assert_eq!(layout.rows, 1);
    }

    #[test]
    fn reads_outside_the_region_are_left_out_of_the_packing() {
        let reads = vec![Read::aligned(0, 50), Read::aligned(10_000, 50)];
        let track = PileupTrack::new(reads);
        let layout = track.layout(&scale(&region()), &region());
        assert_eq!(layout.placed.len(), 1);
        assert_eq!(layout.placed[0].0, 0);
    }

    #[test]
    fn a_deep_pileup_is_capped_and_says_so() {
        let reads: Vec<Read> = (0..100).map(|i| Read::aligned(i, 50)).collect();
        let track = PileupTrack::new(reads).max_rows(Some(10));
        let layout = track.layout(&scale(&region()), &region());
        assert_eq!(layout.rows, 10);
        assert!(layout.hidden > 0);
        assert_eq!(layout.placed.len() + layout.hidden, 100);
    }

    #[test]
    fn the_cap_can_be_lifted() {
        let reads: Vec<Read> = (0..60).map(|i| Read::aligned(i, 50)).collect();
        let track = PileupTrack::new(reads).max_rows(None);
        let layout = track.layout(&scale(&region()), &region());
        assert_eq!(layout.hidden, 0);
        assert!(layout.rows > 40);
    }

    #[test]
    fn the_hidden_count_is_reported_on_the_figure() {
        let reads: Vec<Read> = (0..40).map(|i| Read::aligned(i, 60)).collect();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(PileupTrack::new(reads).max_rows(Some(5)))
            .to_svg();
        assert!(svg.contains("reads not shown"), "no note in the output");
    }

    #[test]
    fn mismatches_are_painted_and_matches_are_not() {
        let reference = vec![b'A'; 200];
        // One read of all A, one of all G, over a reference of A.
        let reads = vec![
            Read::aligned(10, 20).sequence(vec![b'A'; 20]),
            Read::aligned(40, 20).sequence(vec![b'G'; 20]),
        ];
        let svg = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .show_region_label(false)
            .push(PileupTrack::new(reads).reference(0, reference))
            .to_svg();

        let guanine = crate::theme::BaseColors::default().g;
        let adenine = crate::theme::BaseColors::default().a;
        assert_eq!(
            svg.matches(&guanine).count(),
            20,
            "every mismatched base should be painted"
        );
        assert!(
            !svg.contains(&adenine),
            "matching bases should not be painted"
        );
    }

    #[test]
    fn mismatch_hunting_is_case_insensitive() {
        let reads = vec![Read::aligned(10, 5).sequence(b"acgta".to_vec())];
        let svg = Figure::new(Region::new("chr1", 0, 40).unwrap())
            .show_region_label(false)
            .push(PileupTrack::new(reads).reference(10, b"ACGTA".to_vec()))
            .to_svg();
        for color in [
            crate::theme::BaseColors::default().c,
            crate::theme::BaseColors::default().g,
            crate::theme::BaseColors::default().t,
        ] {
            assert!(
                !svg.contains(&color),
                "soft masked bases are not mismatches"
            );
        }
    }

    #[test]
    fn a_mismatch_is_never_thinner_than_a_pixel() {
        // Twenty thousand bases across the figure puts a base at a fraction of
        // a pixel, but a mismatch still has to be visible.
        let reads = vec![Read::aligned(100, 30).sequence(vec![b'G'; 30])];
        let svg = Figure::new(Region::new("chr1", 0, 20_000).unwrap())
            .show_region_label(false)
            .push(
                PileupTrack::new(reads)
                    .reference(100, vec![b'A'; 30])
                    .mismatch_threshold(0.01),
            )
            .to_svg();

        let guanine = crate::theme::BaseColors::default().g;
        let widths: Vec<f64> = svg
            .match_indices(r#"<rect x="#)
            .filter_map(|(index, _)| {
                let element = &svg[index..svg[index..].find("/>").map(|e| index + e)?];
                if !element.contains(&guanine) {
                    return None;
                }
                let start = element.find(r#"width=""#)? + 7;
                let rest = &element[start..];
                rest[..rest.find('"')?].parse().ok()
            })
            .collect();

        assert!(!widths.is_empty(), "no mismatch was drawn");
        assert!(
            widths.iter().all(|w| *w >= 1.0),
            "a mismatch came out at {widths:?}"
        );
    }

    #[test]
    fn without_a_reference_no_mismatch_can_be_found() {
        let reads = vec![Read::aligned(10, 20).sequence(vec![b'G'; 20])];
        let svg = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .show_region_label(false)
            .push(PileupTrack::new(reads))
            .to_svg();
        assert!(!svg.contains(&crate::theme::BaseColors::default().g));
    }

    #[test]
    fn a_wide_view_does_not_go_looking_for_mismatches() {
        // Ten megabases across the figure: a mismatch would be invisible, and
        // finding one would mean walking every base of every read.
        let reads = vec![Read::aligned(1_000, 150).sequence(vec![b'G'; 150])];
        let svg = Figure::new(Region::new("chr1", 0, 10_000_000).unwrap())
            .show_region_label(false)
            .push(PileupTrack::new(reads).reference(1_000, vec![b'A'; 150]))
            .to_svg();
        assert!(!svg.contains(&crate::theme::BaseColors::default().g));
    }

    #[test]
    fn strand_colouring_uses_two_colours() {
        let reads = vec![
            Read::aligned(10, 20).strand(Strand::Forward),
            Read::aligned(50, 20).strand(Strand::Reverse),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(PileupTrack::new(reads).coloring(ReadColoring::Strand))
            .to_svg();
        let theme = Theme::light();
        assert!(svg.contains(theme.color(0)));
        assert!(svg.contains(theme.color(1)));
    }

    #[test]
    fn low_quality_reads_are_faded_when_asked() {
        let track = PileupTrack::new(Vec::new()).fade_by_quality(true);
        let solid = Read::aligned(0, 10).mapping_quality(60);
        let poor = Read::aligned(0, 10).mapping_quality(0);
        let unknown = Read::aligned(0, 10);
        assert_eq!(track.opacity(&solid), 1.0);
        assert!(track.opacity(&poor) < 0.5);
        assert_eq!(track.opacity(&unknown), 1.0);

        let off = PileupTrack::new(Vec::new());
        assert_eq!(off.opacity(&poor), 1.0);
    }

    #[test]
    fn an_empty_pileup_renders_nothing_and_does_not_panic() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(PileupTrack::new(Vec::new()).label("reads"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn a_deletion_draws_a_line_rather_than_a_block() {
        let reads = vec![Read::new(
            10,
            vec![
                CigarOp::Match(20),
                CigarOp::Deletion(10),
                CigarOp::Match(20),
            ],
        )];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(PileupTrack::new(reads))
            .to_svg();
        assert!(svg.contains("<line"), "the gap should be spanned");
        assert!(svg.contains("<rect"), "the aligned blocks should be drawn");
    }

    #[test]
    fn an_insertion_draws_its_own_mark() {
        let reads = vec![Read::new(
            10,
            vec![
                CigarOp::Match(20),
                CigarOp::Insertion(4),
                CigarOp::Match(20),
            ],
        )];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(PileupTrack::new(reads))
            .to_svg();
        assert!(svg.contains("#8e44ad"), "insertions get their own colour");
    }

    #[test]
    fn track_height_follows_the_row_cap() {
        let reads: Vec<Read> = (0..100).map(|i| Read::aligned(i, 50)).collect();
        let capped = PileupTrack::new(reads.clone()).max_rows(Some(10));
        let shallow = PileupTrack::new(vec![Read::aligned(0, 10)]).max_rows(Some(10));
        let s = scale(&region());
        assert_eq!(capped.height(&s), 10.0 * 9.0 + 9.0 * 2.0);
        // A track with one read does not reserve ten rows of empty space.
        assert_eq!(shallow.height(&s), 9.0);
    }
}
