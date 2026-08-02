//! A ruler in codons, so a coding sequence can be read in protein coordinates.
//!
//! Every other ruler in the crate counts bases. This one counts residues,
//! because that is the unit a coding change is named and argued about in.
//! [`CodonTrack`] is given the span of a coding sequence and the strand it is
//! read off, and everything here, the numbering and the translation, follows
//! from that pair.
//!
//! # The numbering is a partition, not a division by three
//!
//! A codon number belongs to the coding sequence and not to the base underneath
//! it. [`CodonTrack::codon_of`] takes a base to the codon holding it and
//! [`CodonTrack::span_of`] takes a codon back to the three bases it owns, and
//! both of them count from the far end of the span when the strand is reverse.
//! A length that is not a multiple of three keeps its whole codons and drops
//! the remainder, at the low end on the reverse strand and the high end on the
//! forward one.
//!
//! # Translating is a table, and the table is a choice
//!
//! Bases are optional. Hand them to [`CodonTrack::sequence`] and every codon
//! with room for a letter carries one; leave them out and the ruler is numbers
//! alone. Which letter a codon gets is [`CodonTrack::genetic_code`]'s to
//! decide, and it is asked for rather than assumed, because nothing about a
//! figure translated with the wrong table looks wrong.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::axis::{group_thousands, nice_step};
use crate::track::feature::span_label;
use crate::track::{strand_color, DrawContext, Strand, Track};

/// NCBI translation table 1, the standard code, indexed by
/// `base * 16 + base * 4 + base` with T, C, A and G reading 0 to 3.
///
/// Table 11, which bacteria, archaea and plastids use, differs from it only in
/// which codons may start a protein, so the residues are identical. The
/// mitochondrial and ciliate tables do reassign residues, and
/// [`CodonTrack::genetic_code`] takes one of those in the same form.
const RESIDUES: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// A ruler counting codons across a coding sequence.
///
/// A variant in a coding sequence is named by residue: BRAF V600E, TP53 R175H,
/// rpoB S450L. A figure drawn in bases cannot be pointed at with any of those
/// names. This is the sibling of [`AxisTrack`](crate::AxisTrack) that can: it
/// partitions the sequence into codons, numbers them, and translates them when
/// there is room for a letter.
///
/// The partition is itself the claim. Two changes at different bases of one
/// codon are competing alleles at one residue rather than a double mutant, and
/// two changes in neighbouring codons are two substitutions however few bases
/// apart they are. Neither statement can be made on a ruler of bases.
///
/// On the reverse strand codon 1 sits at the **highest** coordinate and the
/// numbering runs right to left, which is the whole reason this is a track and
/// not a division by three. Roughly half the coding sequences in any annotation
/// run backwards, and getting their numbering wrong is silent: the figure still
/// draws, it just names the wrong residue.
///
/// ```
/// use karyon::{CodonTrack, Figure, Region, Strand};
///
/// // rpoB, forward strand. The resistance determining region is codons 426 to
/// // 452, so S450 is the codon holding base 761,154.
/// let track = CodonTrack::new(759_806, 763_325, Strand::Forward);
/// assert_eq!(track.codon_of(761_154), Some(450));
/// assert_eq!(track.span_of(450), Some((761_153, 761_156)));
///
/// let svg = Figure::new(Region::parse("NC_000962.3:761140-761180").unwrap())
///     .push(track.label("rpoB"))
///     .to_svg();
/// assert!(svg.contains("450"));
/// ```
#[derive(Debug, Clone)]
pub struct CodonTrack {
    start: u64,
    end: u64,
    strand: Strand,
    seq: Option<(u64, Vec<u8>)>,
    code: [u8; 64],
    cell_height: f64,
    label: Option<String>,
    color: Option<String>,
    numbering: Option<u64>,
    show_start: bool,
}

impl CodonTrack {
    /// A ruler over the coding sequence `start..end`, 0-based half-open.
    ///
    /// A GFF or GenBank CDS is 1-based inclusive, so its `start` is `start - 1`
    /// here and its `end` is unchanged.
    pub fn new(start: u64, end: u64, strand: Strand) -> Self {
        let (start, end) = if end >= start {
            (start, end)
        } else {
            (end, start)
        };
        CodonTrack {
            start,
            end,
            strand,
            seq: None,
            code: *RESIDUES,
            cell_height: 13.0,
            label: None,
            color: None,
            numbering: None,
            show_start: true,
        }
    }

    /// Attaches the reference bases, so codons can be translated.
    ///
    /// `start` is the 0-based position of the first base of `bases`, the same
    /// convention as [`SequenceTrack`](crate::SequenceTrack), and the slice does
    /// not have to cover the whole coding sequence: a codon whose bases are
    /// missing is drawn without a letter rather than guessed at.
    ///
    /// On the reverse strand the bases are complemented and read backwards, so
    /// the caller hands in the reference as it is and never the reverse
    /// complement.
    pub fn sequence(mut self, start: u64, bases: impl Into<Vec<u8>>) -> Self {
        self.seq = Some((start, bases.into()));
        self
    }

    /// Translates with a genetic code other than the standard one.
    ///
    /// The table is the sixty-four residues in NCBI order, `TTT` `TTC` `TTA`
    /// `TTG` `TCT` and so on to `GGG`, which is the `AAs` line of an NCBI
    /// translation table copied across unchanged.
    ///
    /// Table 1 is the default and table 11, which bacteria, archaea and plastids
    /// use, gives the same residues, so this is only needed for the codes that
    /// reassign one. Vertebrate mitochondria read `TGA` as tryptophan and `AGA`
    /// and `AGG` as stops; several ciliates read `TAA` and `TAG` as glutamine.
    /// Translating those with the standard code produces a plausible protein
    /// that is wrong, which is worse than refusing, so the table is a knob
    /// rather than an assumption.
    ///
    /// ```
    /// use karyon::{CodonTrack, Strand};
    ///
    /// // Vertebrate mitochondrial, NCBI table 2.
    /// let mito = CodonTrack::new(0, 3, Strand::Forward)
    ///     .sequence(0, b"TGA".to_vec())
    ///     .genetic_code(b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG");
    /// assert_eq!(mito.residue_of(1), Some(b'W'));
    /// ```
    pub fn genetic_code(mut self, residues: &[u8; 64]) -> Self {
        self.code = *residues;
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of the codon cells in pixels.
    pub fn cell_height(mut self, height: f64) -> Self {
        self.cell_height = height.max(4.0);
        self
    }

    /// Overrides the colour, which is otherwise the crate's strand colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Numbers every `step` codons instead of at a spacing chosen for the zoom.
    pub fn numbering(mut self, step: u64) -> Self {
        self.numbering = Some(step.max(1));
        self
    }

    /// Whether to mark the first codon with a chevron pointing the way the
    /// sequence is read. On by default, and the only thing on a reverse strand
    /// ruler that says so at a glance.
    pub fn show_start(mut self, show: bool) -> Self {
        self.show_start = show;
        self
    }

    /// How many whole codons the sequence holds.
    ///
    /// A coding sequence whose length is not a multiple of three is not
    /// rejected: the trailing partial codon is left out of the count and out of
    /// the ruler, since a third of a residue is not a residue.
    pub fn codons(&self) -> u64 {
        (self.end - self.start) / 3
    }

    /// Which strand the sequence is read off.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// The 1-based codon holding the base at 0-based `pos`, if any.
    pub fn codon_of(&self, pos: u64) -> Option<u64> {
        if pos < self.start || pos >= self.end {
            return None;
        }
        let offset = if self.strand == Strand::Reverse {
            self.end - 1 - pos
        } else {
            pos - self.start
        };
        if offset >= self.codons() * 3 {
            return None;
        }
        Some(offset / 3 + 1)
    }

    /// The reference span of the 1-based codon `codon`, 0-based half-open.
    ///
    /// The span is in reference order whichever strand the sequence is on, so
    /// it can be handed straight to a [`Scale`] without a second thought.
    pub fn span_of(&self, codon: u64) -> Option<(u64, u64)> {
        if codon == 0 || codon > self.codons() {
            return None;
        }
        Some(if self.strand == Strand::Reverse {
            (self.end - 3 * codon, self.end - 3 * codon + 3)
        } else {
            (self.start + 3 * (codon - 1), self.start + 3 * codon)
        })
    }

    /// The one-letter residue of the 1-based codon `codon`.
    ///
    /// `Some(b'*')` is a stop. `None` means either no sequence was attached, the
    /// bases for this codon are not in the slice, or one of them is ambiguous.
    pub fn residue_of(&self, codon: u64) -> Option<u8> {
        let (start, end) = self.span_of(codon)?;
        let (at, bases) = self.seq.as_ref()?;
        let mut index = 0usize;
        for offset in 0..3u64 {
            // Reading order: the last reference base first on the reverse
            // strand, and complemented.
            let pos = if self.strand == Strand::Reverse {
                end - 1 - offset
            } else {
                start + offset
            };
            let base = bases.get(pos.checked_sub(*at)? as usize)?;
            let base = if self.strand == Strand::Reverse {
                complement(*base)
            } else {
                *base
            };
            index = index * 4 + code_of(base)? as usize;
        }
        Some(self.code[index])
    }

    /// Codons overlapping the region on screen, as an inclusive 1-based range.
    fn visible(&self, scale: &Scale) -> Option<(u64, u64)> {
        let total = self.codons();
        if total == 0 {
            return None;
        }
        let (view_start, view_end) = scale.bounds();
        if view_end <= self.start || view_start >= self.end {
            return None;
        }
        let at = |pos: u64| self.codon_of(pos.clamp(self.start, self.end - 1));
        let left = at(view_start)?;
        let right = at(view_end.saturating_sub(1).max(view_start))?;
        Some((left.min(right), left.max(right)))
    }

    /// How many codons to leave between printed numbers at this zoom.
    fn number_step(&self, scale: &Scale, theme: &Theme) -> u64 {
        if let Some(step) = self.numbering {
            return step;
        }
        // Room for the widest number the ruler will print, plus a gap.
        let widest = text_width(&format!("{}", self.codons()), theme.font_size) + 14.0;
        let per_codon = scale.px_per_bp() * 3.0;
        if per_codon <= 0.0 {
            return 1;
        }
        nice_step(widest / per_codon, 1.0)
    }
}

impl Track for CodonTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.cell_height + 4.0 + 12.0
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let theme = ctx.theme;
        let font = theme.font_size;
        let ink = self
            .color
            .clone()
            .unwrap_or_else(|| strand_color(self.strand, theme).to_string());

        let cells_y = band.y + 1.0;
        let cells_h = self.cell_height.min(band.h - font - 4.0).max(2.0);
        let baseline = cells_y + cells_h + font;

        let Some((first, last)) = self.visible(ctx.scale) else {
            return;
        };
        let per_codon = ctx.scale.px_per_bp() * 3.0;

        // Three levels of detail, the same ladder SequenceTrack climbs: letters
        // while they fit, then alternating cells that still show where one
        // codon ends, then a plain rule once a codon is thinner than a hair.
        let letters = self.seq.is_some() && per_codon >= font * 0.95;
        let cells = per_codon >= 2.5;

        let pale = mix(&ink, theme.surface(), 0.86);
        let paler = mix(&ink, theme.surface(), 0.94);

        if cells {
            for codon in first..=last {
                let Some((from, to)) = self.span_of(codon) else {
                    continue;
                };
                let x0 = ctx.scale.x(from);
                let x1 = ctx.scale.x(to);
                let residue = self.residue_of(codon);
                // A cell is drawn only once a codon is two and a half pixels
                // wide, which is the same test that decides it is a thing on
                // the page rather than a stripe of the rule, so every cell
                // drawn here is worth naming.
                ctx.svg
                    .begin_titled(&self.cell_title(codon, from, to, residue));
                let fill = if codon % 2 == 0 { &pale } else { &paler };
                ctx.svg.rect(x0, cells_y, (x1 - x0).max(0.4), cells_h, fill);
                if letters {
                    let glyph = residue.unwrap_or(b'-');
                    let text = (glyph as char).to_string();
                    // A stop is the one residue worth reading as a warning.
                    let color = if glyph == b'*' {
                        theme.bases.t.clone()
                    } else {
                        theme.foreground.clone()
                    };
                    ctx.svg.text(
                        (x0 + x1) / 2.0,
                        cells_y + cells_h - (cells_h - font * 0.72) / 2.0,
                        &text,
                        &color,
                        font,
                        Anchor::Middle,
                    );
                }
                ctx.svg.end_group();
            }
        } else {
            // Below cell resolution the extent of the sequence is still worth
            // drawing, so the numbers underneath have something to belong to.
            let x0 = ctx.scale.x(self.start).max(band.x);
            let x1 = ctx.scale.x(self.end).min(band.right());
            ctx.svg
                .rect(x0, cells_y, (x1 - x0).max(1.0), cells_h, &pale);
        }

        // The rule sits under the cells, a hair inside the band so a clip on the
        // boundary cannot take half its width.
        let rule_y = cells_y + cells_h + 0.5;
        ctx.svg
            .line(band.x, rule_y, band.right(), rule_y, &theme.rule, 1.0);

        let step = self.number_step(ctx.scale, theme);
        let mut codon = first.div_ceil(step) * step;
        if codon == 0 {
            codon = step;
        }
        while codon <= last {
            if let Some((from, to)) = self.span_of(codon) {
                let x = (ctx.scale.x(from) + ctx.scale.x(to)) / 2.0;
                ctx.svg.line(x, rule_y, x, rule_y + 3.0, &theme.rule, 1.0);
                let text = format!("{codon}");
                let half = text_width(&text, font) / 2.0;
                // A number half off the edge of the figure is a number that
                // cannot be read, so it is dropped rather than clipped.
                if x - half >= band.x && x + half <= band.right() {
                    ctx.svg
                        .text(x, baseline, &text, &theme.muted, font, Anchor::Middle);
                }
            }
            codon += step;
        }

        if self.show_start {
            self.draw_start(ctx, cells_y, cells_h, &ink);
        }
    }
}

impl CodonTrack {
    /// What one codon cell is, for its tooltip.
    ///
    /// The span is written 1-based inclusive, the form the ruler under it
    /// prints and the form a variant is quoted in, and the residue comes last
    /// because it is the thing the cell already shows when there is room.
    fn cell_title(&self, codon: u64, from: u64, to: u64, residue: Option<u8>) -> String {
        let span = span_label(from, to);
        // The codon number is a count like any other number in a tooltip, so
        // a gene long enough to reach four figures reads `codon 1,234`.
        let at = group_thousands(codon);
        match residue {
            Some(glyph) => format!("codon {at}, {span}, {}", glyph as char),
            None => format!("codon {at}, {span}"),
        }
    }

    /// A chevron on the first codon, pointing the way the sequence is read.
    fn draw_start(&self, ctx: &mut DrawContext<'_>, y: f64, h: f64, ink: &str) {
        let Some((from, to)) = self.span_of(1) else {
            return;
        };
        let band = ctx.band;
        let (x0, x1) = (ctx.scale.x(from), ctx.scale.x(to));
        if x1 < band.x || x0 > band.right() {
            return;
        }
        let reverse = self.strand == Strand::Reverse;
        // The chevron points out of the sequence, away from codon 2, so on the
        // reverse strand it points right from the highest coordinate.
        let tip = if reverse {
            x1.min(band.right()) + 6.0
        } else {
            x0.max(band.x) - 6.0
        };
        let root = if reverse {
            x1.min(band.right())
        } else {
            x0.max(band.x)
        };
        let mid = y + h / 2.0;
        ctx.svg.polygon(
            &[
                (root, y + 1.0),
                (tip, mid),
                (root, y + h - 1.0),
                (root + if reverse { -3.0 } else { 3.0 }, mid),
            ],
            ink,
        );
    }
}

/// Complement of one base, leaving anything unrecognised alone so it fails the
/// lookup rather than silently becoming an A.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        other => other,
    }
}

/// T, C, A and G as 0 to 3, the ordering [`RESIDUES`] is written in.
fn code_of(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::Region;

    fn scale(start: u64, end: u64) -> Scale {
        let region = Region::new("chr", start, end).unwrap();
        Scale::new(&region, 0.0, 900.0)
    }

    #[test]
    fn a_forward_codon_is_three_bases_from_the_start() {
        let track = CodonTrack::new(1_000, 1_030, Strand::Forward);
        assert_eq!(track.codons(), 10);
        assert_eq!(track.span_of(1), Some((1_000, 1_003)));
        assert_eq!(track.span_of(10), Some((1_027, 1_030)));
        assert_eq!(track.codon_of(1_000), Some(1));
        assert_eq!(track.codon_of(1_002), Some(1));
        assert_eq!(track.codon_of(1_003), Some(2));
    }

    #[test]
    fn a_reverse_codon_is_counted_from_the_far_end() {
        // This is the case the track exists for: katG and pncA both run
        // backwards, so codon 1 is at the highest coordinate.
        let track = CodonTrack::new(1_000, 1_030, Strand::Reverse);
        assert_eq!(track.span_of(1), Some((1_027, 1_030)));
        assert_eq!(track.span_of(10), Some((1_000, 1_003)));
        assert_eq!(track.codon_of(1_029), Some(1));
        assert_eq!(track.codon_of(1_027), Some(1));
        assert_eq!(track.codon_of(1_026), Some(2));
        assert_eq!(track.codon_of(1_000), Some(10));
    }

    #[test]
    fn every_codon_maps_back_to_itself_on_both_strands() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let track = CodonTrack::new(500, 599, strand);
            for codon in 1..=track.codons() {
                let (from, to) = track.span_of(codon).unwrap();
                assert_eq!(to - from, 3);
                for pos in from..to {
                    assert_eq!(track.codon_of(pos), Some(codon), "{strand:?} {codon}");
                }
            }
        }
    }

    #[test]
    fn a_trailing_partial_codon_is_not_a_codon() {
        let forward = CodonTrack::new(0, 11, Strand::Forward);
        assert_eq!(forward.codons(), 3);
        // Bases 9 and 10 are two thirds of a residue, so they have no number.
        assert_eq!(forward.codon_of(8), Some(3));
        assert_eq!(forward.codon_of(9), None);

        // On the reverse strand the leftover is at the low end instead.
        let reverse = CodonTrack::new(0, 11, Strand::Reverse);
        assert_eq!(reverse.codon_of(10), Some(1));
        assert_eq!(reverse.codon_of(2), Some(3));
        assert_eq!(reverse.codon_of(1), None);
    }

    #[test]
    fn outside_the_sequence_there_is_no_codon() {
        let track = CodonTrack::new(100, 130, Strand::Forward);
        assert_eq!(track.codon_of(99), None);
        assert_eq!(track.codon_of(130), None);
        assert_eq!(track.span_of(0), None);
        assert_eq!(track.span_of(11), None);
    }

    #[test]
    fn the_genetic_code_is_the_standard_one() {
        let track = CodonTrack::new(0, 12, Strand::Forward).sequence(0, b"ATGTGGTAATTT".to_vec());
        assert_eq!(track.residue_of(1), Some(b'M'));
        assert_eq!(track.residue_of(2), Some(b'W'));
        assert_eq!(track.residue_of(3), Some(b'*'));
        assert_eq!(track.residue_of(4), Some(b'F'));
        assert_eq!(RESIDUES.len(), 64);
    }

    #[test]
    fn a_different_genetic_code_reassigns_residues() {
        // Vertebrate mitochondrial, NCBI table 2: TGA is tryptophan rather than
        // a stop, ATA is methionine, and AGA and AGG are stops.
        const TABLE_2: &[u8; 64] =
            b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
        let seq = b"TGAATAAGA".to_vec();

        let standard = CodonTrack::new(0, 9, Strand::Forward).sequence(0, seq.clone());
        assert_eq!(standard.residue_of(1), Some(b'*'));
        assert_eq!(standard.residue_of(2), Some(b'I'));
        assert_eq!(standard.residue_of(3), Some(b'R'));

        let mito = CodonTrack::new(0, 9, Strand::Forward)
            .sequence(0, seq)
            .genetic_code(TABLE_2);
        assert_eq!(mito.residue_of(1), Some(b'W'));
        assert_eq!(mito.residue_of(2), Some(b'M'));
        assert_eq!(mito.residue_of(3), Some(b'*'));
    }

    #[test]
    fn table_eleven_gives_the_same_residues_as_the_standard_code() {
        // Bacteria and archaea use table 11, which differs only in which codons
        // may start a protein, so it is not a separate table here.
        assert_eq!(
            RESIDUES,
            b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        );
    }

    #[test]
    fn a_reverse_sequence_is_complemented_and_read_backwards() {
        // The reference reads TTACCACAT. Backwards and complemented that is
        // ATGTGGTAA, so the reverse strand codons are M, W and a stop.
        let track = CodonTrack::new(0, 9, Strand::Reverse).sequence(0, b"TTACCACAT".to_vec());
        assert_eq!(track.residue_of(1), Some(b'M'));
        assert_eq!(track.residue_of(2), Some(b'W'));
        assert_eq!(track.residue_of(3), Some(b'*'));
    }

    #[test]
    fn an_ambiguous_or_missing_base_translates_to_nothing() {
        let track = CodonTrack::new(0, 6, Strand::Forward).sequence(0, b"ATGNNN".to_vec());
        assert_eq!(track.residue_of(1), Some(b'M'));
        assert_eq!(track.residue_of(2), None, "N is not a base");

        // A slice that does not reach the codon gives no letter rather than a
        // wrong one.
        let short = CodonTrack::new(0, 6, Strand::Forward).sequence(0, b"ATG".to_vec());
        assert_eq!(short.residue_of(2), None);

        // And neither does one that starts after it.
        let late = CodonTrack::new(0, 6, Strand::Forward).sequence(3, b"TGG".to_vec());
        assert_eq!(late.residue_of(1), None);
        assert_eq!(late.residue_of(2), Some(b'W'));
    }

    #[test]
    fn only_the_codons_on_screen_are_drawn() {
        let track = CodonTrack::new(1_000, 1_300, Strand::Forward);
        assert_eq!(track.visible(&scale(1_000, 1_300)), Some((1, 100)));
        assert_eq!(track.visible(&scale(1_030, 1_060)), Some((11, 20)));
        assert_eq!(track.visible(&scale(0, 900)), None, "not in view");
        assert_eq!(track.visible(&scale(2_000, 2_900)), None);
    }

    #[test]
    fn the_reverse_ruler_still_reports_a_low_to_high_codon_range() {
        let track = CodonTrack::new(1_000, 1_300, Strand::Reverse);
        let (first, last) = track.visible(&scale(1_030, 1_060)).unwrap();
        assert!(first <= last, "{first} to {last}");
        assert_eq!((first, last), (81, 90));
    }

    #[test]
    fn numbers_thin_out_as_the_view_widens() {
        let track = CodonTrack::new(0, 3_000, Strand::Forward);
        let theme = Theme::light();
        let close = track.number_step(&scale(0, 90), &theme);
        let far = track.number_step(&scale(0, 3_000), &theme);
        assert!(close < far, "{close} then {far}");
        assert!(close >= 1);
    }

    #[test]
    fn the_ruler_renders_and_names_its_codons() {
        use crate::figure::Figure;

        let bases: Vec<u8> = b"ATG".iter().cycle().take(300).copied().collect();
        let svg = Figure::new(Region::new("chr", 1_000, 1_060).unwrap())
            .push(
                CodonTrack::new(1_000, 1_300, Strand::Forward)
                    .sequence(1_000, bases)
                    .label("rpoB"),
            )
            .to_svg();
        assert!(svg.contains(">M<"), "codons are translated");
        assert!(svg.contains("rpoB"));
    }

    #[test]
    fn a_codon_cell_names_its_number_its_span_and_its_residue() {
        use crate::figure::Figure;

        let bases: Vec<u8> = b"ATG".iter().cycle().take(300).copied().collect();
        let svg = Figure::new(Region::new("chr", 1_000, 1_060).unwrap())
            .show_region_label(false)
            .push(
                CodonTrack::new(1_000, 1_300, Strand::Forward)
                    .sequence(1_000, bases)
                    .label("rpoB"),
            )
            .to_svg();
        // 1-based inclusive, the form the ruler prints and a variant is quoted
        // in, and the residue the cell itself is showing.
        assert!(
            svg.contains("<title>codon 1, 1,001 to 1,003, M</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>codon 20, 1,058 to 1,060, M</title>"));
        // Twenty codons on screen, one group each and all of them closed.
        assert_eq!(svg.matches("<title>").count(), 20);
        assert_eq!(svg.matches("<g>").count(), 20);
    }

    #[test]
    fn a_zero_length_codon_cell_still_reads_forwards() {
        // The span formatter is shared with every other track precisely so
        // that the degenerate case is guarded in one place: written by hand
        // this came out `1,001 to 1,000`.
        let track = CodonTrack::new(1_000, 1_060, Strand::Forward);
        assert_eq!(
            track.cell_title(7, 1_000, 1_000, Some(b'M')),
            "codon 7, 1,001 to 1,001, M"
        );
        // And the codon number is grouped like every other count.
        assert_eq!(
            track.cell_title(1_234, 1_000, 1_003, None),
            "codon 1,234, 1,001 to 1,003"
        );
    }

    #[test]
    fn a_codon_with_no_bases_behind_it_is_still_named() {
        use crate::figure::Figure;

        let svg = Figure::new(Region::new("chr", 1_000, 1_060).unwrap())
            .show_region_label(false)
            .push(CodonTrack::new(1_000, 1_300, Strand::Reverse))
            .to_svg();
        // Reverse strand: codon 1 sits at the far end, so this window holds
        // the high eighties.
        assert!(
            svg.contains("<title>codon 81, 1,058 to 1,060</title>"),
            "{svg}"
        );
    }

    #[test]
    fn a_codon_thinner_than_a_cell_is_not_named() {
        use crate::figure::Figure;

        // Ninety kilobases across the figure: a codon is a hundredth of a
        // pixel, no cell is drawn, and naming one would name nothing.
        let svg = Figure::new(Region::new("chr", 0, 90_000).unwrap())
            .show_region_label(false)
            .push(CodonTrack::new(0, 90_000, Strand::Forward))
            .to_svg();
        assert!(!svg.contains("<title>"), "{svg}");
    }

    #[test]
    fn a_sequence_shorter_than_one_codon_draws_nothing() {
        let track = CodonTrack::new(10, 12, Strand::Forward);
        assert_eq!(track.codons(), 0);
        assert_eq!(track.visible(&scale(0, 100)), None);
    }
}
