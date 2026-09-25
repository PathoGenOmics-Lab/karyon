//! The reference bases themselves.
//!
//! What a base can look like depends entirely on how much room it has, so this
//! track has three regimes rather than one. Wider than
//! [`SequenceTrack::letter_threshold`] pixels per base there is space for a
//! letter and the letter is drawn. Narrower than that each base is a coloured
//! block, which still carries the composition of the stretch even though no
//! individual base is legible. Narrower than
//! [`SequenceTrack::block_threshold`] the track draws nothing at all and prints
//! a line of text asking for more zoom.
//!
//! # The floor is not politeness
//!
//! A megabase of blocks is a million rectangles in the file and tens of
//! megabytes of SVG, spent on marks thinner than the pixels they are drawn on.
//! Above the floor the cost follows the width of the figure rather than the
//! length of the sequence handed in, since only the bases inside the region are
//! visited: a whole chromosome can be given to a track showing two hundred
//! bases of it.
//!
//! # Nothing is left uncoloured, and case is lost
//!
//! Colour comes from the theme's [`BaseColors`](crate::BaseColors), where `U`
//! is the colour of `T` and everything the palette does not name, an `N` or a
//! gap, takes one neutral colour rather than being dropped. Lower case is drawn
//! as upper case, so soft masking survives as sequence but not as a mark: a
//! repeat-masked region has to be a [`FeatureTrack`](crate::FeatureTrack) over
//! the top if it is to be visible.

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{finite_within, Anchor};
use crate::theme::Theme;
use crate::track::legend::Legend;
use crate::track::{DrawContext, Track};

/// How many pixels a base needs before its letter is drawn, in the reference
/// and in a pileup's mismatches alike.
pub(crate) const LETTER_PX: f64 = 7.0;

/// Nucleotides drawn as coloured blocks, with letters once there is room.
///
/// The track adapts to the zoom level the way a genome browser does: letters
/// when a base is at least `letter_threshold` pixels wide, plain blocks when it
/// is narrower, and a hint to zoom in when the bases would be thinner than a
/// pixel. That last case matters, because drawing five million one-pixel
/// rectangles produces a file no viewer will open.
///
/// ```
/// use karyon::{Figure, Region, SequenceTrack};
///
/// let svg = Figure::new(Region::parse("chr1:1-20").unwrap())
///     .push(SequenceTrack::new(0, b"ACGTACGTACGTACGTACGT").label("ref"))
///     .to_svg();
/// assert!(svg.contains(">A</text>"));
/// ```
#[derive(Debug, Clone)]
pub struct SequenceTrack {
    start: u64,
    seq: Vec<u8>,
    label: Option<String>,
    height: f64,
    letter_threshold: f64,
    block_threshold: f64,
}

impl SequenceTrack {
    /// A track whose `seq[i]` is the base at 0-based position `start + i`.
    pub fn new(start: u64, seq: impl Into<Vec<u8>>) -> Self {
        SequenceTrack {
            start,
            seq: seq.into(),
            label: None,
            height: 18.0,
            letter_threshold: LETTER_PX,
            block_threshold: 0.6,
        }
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

    /// Sets how many pixels a base needs before its letter is drawn.
    pub fn letter_threshold(mut self, pixels: f64) -> Self {
        self.letter_threshold = finite_within(pixels, 0.0, 200.0, 7.0);
        self
    }

    /// Sets how many pixels a base needs before any block is drawn at all.
    ///
    /// Below this the track prints a hint instead. Raise it to keep dense
    /// figures light, lower it to force blocks at a wider zoom.
    pub fn block_threshold(mut self, pixels: f64) -> Self {
        self.block_threshold = finite_within(pixels, 0.0, 200.0, 1.5);
        self
    }

    /// The base at a 0-based position, if the track carries it.
    pub fn base_at(&self, pos: u64) -> Option<u8> {
        if pos < self.start {
            return None;
        }
        self.seq.get((pos - self.start) as usize).copied()
    }
}

impl Track for SequenceTrack {
    fn noun(&self) -> &str {
        "the reference sequence"
    }

    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    /// The bases' colours, while they are drawn as blocks: a letter names
    /// itself, and under the floor nothing is drawn to name.
    fn key(&self, region: &Region, px_per_bp: f64, theme: &Theme) -> Option<Legend> {
        if px_per_bp < self.block_threshold || px_per_bp >= self.letter_threshold {
            return None;
        }
        let first = region.start().max(self.start);
        let last = region
            .end()
            .min(self.start.saturating_add(self.seq.len() as u64));
        if last <= first {
            return None;
        }
        // An N or an ambiguity code is one neutral colour, named only where
        // the window holds one.
        let other = (first..last)
            .filter_map(|pos| self.base_at(pos))
            .any(|base| !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U'));
        let key = theme.bases.legend();
        Some(if other {
            key.key("other", theme.bases.other.clone())
        } else {
            key
        })
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let px_per_bp = ctx.scale.px_per_bp();

        if px_per_bp < self.block_threshold {
            ctx.svg.text(
                band.x + 3.0,
                band.mid_y() + ctx.theme.font_size * 0.35,
                "zoom in to see bases",
                &ctx.theme.muted,
                ctx.theme.font_size - 1.0,
                Anchor::Start,
            );
            return;
        }

        // Only the bases actually on screen, so the cost follows the figure
        // width rather than the length of the sequence handed in.
        let first = ctx.region.start().max(self.start);
        let last = ctx
            .region
            .end()
            .min(self.start.saturating_add(self.seq.len() as u64));
        if last <= first {
            return;
        }

        let draw_letters = px_per_bp >= self.letter_threshold;
        let font_size = (px_per_bp * 0.72).min(band.h * 0.8);
        let baseline = band.mid_y() + font_size * 0.35;

        for pos in first..last {
            let Some(base) = self.base_at(pos) else {
                continue;
            };
            let x = ctx.scale.x(pos);
            let width = ctx.scale.x(pos.saturating_add(1)) - x;
            let color = ctx.theme.bases.of(base).to_string();

            if draw_letters {
                let letter = (base as char).to_ascii_uppercase().to_string();
                ctx.svg.text_bold(
                    x + width / 2.0,
                    baseline,
                    &letter,
                    &color,
                    font_size,
                    Anchor::Middle,
                );
            } else {
                // Overdraw by a hair so neighbouring blocks do not show a
                // hairline gap from antialiasing.
                ctx.svg.rect(x, band.y, width + 0.25, band.h, &color);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;
    use crate::theme::BaseColors;

    #[test]
    fn base_at_maps_positions_through_the_offset() {
        let track = SequenceTrack::new(100, b"ACGT");
        assert_eq!(track.base_at(100), Some(b'A'));
        assert_eq!(track.base_at(103), Some(b'T'));
        assert_eq!(track.base_at(104), None);
        assert_eq!(track.base_at(99), None);
    }

    /// A block of colour names its base only to a reader who knows the
    /// palette, and a letter needs no key.
    #[test]
    fn the_bases_are_keyed_only_while_they_are_blocks() {
        let theme = Theme::light();
        let track = SequenceTrack::new(0, b"ACGTACGTAC".repeat(100));
        let region = Region::new("chr1", 0, 1_000).unwrap();
        let key = track.key(&region, 2.0, &theme).expect("blocks are keyed");
        let labels: Vec<String> = key
            .items()
            .iter()
            .map(|item| match item {
                crate::track::legend::LegendItem::Key { label, .. } => label.clone(),
                other => panic!("not a key: {other:?}"),
            })
            .collect();
        assert_eq!(labels, ["A", "C", "G", "T"]);
        assert_eq!(
            track.key(&region, 8.0, &theme),
            None,
            "letters name themselves"
        );
        assert_eq!(track.key(&region, 0.1, &theme), None, "nothing is drawn");
        // An N is named, as the neutral colour, where the window holds one.
        let masked = SequenceTrack::new(0, b"ACGTNNACGT".repeat(100));
        let key = masked.key(&region, 2.0, &theme).unwrap();
        assert_eq!(key.len(), 5);
    }

    #[test]
    fn wide_regions_print_a_hint_instead_of_millions_of_rectangles() {
        let region = Region::parse("chr1:1-1000000").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(SequenceTrack::new(0, vec![b'A'; 1_000_000]))
            .to_svg();
        assert!(svg.contains("zoom in to see bases"));
        assert!(
            !svg.contains(&BaseColors::default().a),
            "no base should have been painted"
        );
        assert!(svg.len() < 4_000, "output grew to {} bytes", svg.len());
    }

    #[test]
    fn medium_zoom_draws_blocks_without_letters() {
        let region = Region::parse("chr1:1-200").unwrap();
        let svg = Figure::new(region)
            .width(600.0)
            .show_region_label(false)
            .push(SequenceTrack::new(0, vec![b'C'; 200]))
            .to_svg();
        assert_eq!(svg.matches(&BaseColors::default().c).count(), 200);
        assert!(!svg.contains("</text>"));
    }

    #[test]
    fn sequence_shorter_than_the_region_draws_only_what_it_has() {
        let region = Region::parse("chr1:1-40").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(SequenceTrack::new(0, b"ACGT"))
            .to_svg();
        assert_eq!(svg.matches("</text>").count(), 4);
    }

    #[test]
    fn sequence_outside_the_region_draws_nothing() {
        let region = Region::parse("chr1:1-40").unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(SequenceTrack::new(10_000, b"ACGT"))
            .to_svg();
        assert!(!svg.contains("</text>"));
    }
}
