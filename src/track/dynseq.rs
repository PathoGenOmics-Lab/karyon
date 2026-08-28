//! Per-base model attribution: one reference base carrying one signed number.
//!
//! A model that predicts something from a sequence can be asked which bases it
//! used, and the answer is one number per base with a sign: this base pushed
//! the prediction up, that one pulled it down. The convention for drawing it,
//! which BPNet, DeepSHAP and the dynseq browser track all use, is to draw the
//! base itself at a height proportional to its score, hanging below the line
//! where the score is negative. A motif then appears as a word.
//!
//! # Why this is not a logo
//!
//! A [`LogoTrack`](crate::LogoTrack) normalises within a column. Handed one
//! symbol carrying one weight it has nothing to divide by but that weight, so
//! the symbol takes the whole column whatever the number was. Measured, with
//! four bases and one symbol per column: `0.1` and `0.9` both come out at
//! height `1.0` under `Probability` and both at `8.65` under `LogOdds`. And
//! `LogoColumn::add` clamps a negative weight to zero before any score is
//! chosen, so a base that pulled the prediction down draws as nothing at all.
//!
//! The magnitude and the sign are the whole measurement, and a logo destroys
//! both. Its below-the-line half is not a sign channel: it is depletion against
//! a background, worked out downstream from probabilities that are all
//! positive.
//!
//! # Why this is not a window track either, which is the closer call
//!
//! [`WindowTrack`](crate::WindowTrack) draws a signed statistic against a line
//! and reduces each pixel column to its lowest and highest value, which is what
//! this does below one pixel per base. Two things separate them, and neither is
//! a setting.
//!
//! The first is the datum. A window is an interval carrying a statistic, and a
//! base is not an interval: a window track over a megabase of per-base scores
//! is a million `Window` values, each one base long, to say what a sequence and
//! a vector of numbers say. The second is that a base has an identity. Above
//! about five pixels a base the letter is drawn, and the letter is the reason
//! the figure is read at all: a reader looking for `CAGGTG` is looking for a
//! word, not for six bars.
//!
//! # Three regimes, and the zoom picks
//!
//! Letters where a letter fits, bars where the box is at least a pixel wide,
//! and below that an envelope of the extremes in one neutral ink. The envelope
//! is never in a base colour, because a column spanning forty bases has no base
//! and colouring it green would claim one.
//!
//! There is no aggregate to choose here and that is a decision rather than an
//! omission. Every answer [`Aggregate`](crate::Aggregate) offers is wrong for a
//! signed score: a maximum hides a strong negative, a minimum hides a strong
//! positive, and a mean cancels a `+2` against a `-2` into a nought that says
//! the model ignored the place. Both extremes, or nothing.
//!
//! # A base with no score is not a base scoring nought
//!
//! The rule is not one line across the band. It is one line under each run of
//! bases that carry a score, so a stretch the model was never run over has no
//! rule under it. Without that, a score of exactly nought and a base nobody
//! scored would both draw no glyph and look identical, which is this crate's
//! named mistake exactly.

use crate::scale::Scale;
use crate::style::QuantitativeAxis;
use crate::svg::{finite_within, text_rounded, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Rect, Track};
use crate::Region;

/// Per-base attribution scores drawn as the bases themselves.
///
/// ```
/// use karyon::{plot, DynseqTrack};
///
/// // A weak base, a strong one, one the model pushed the other way, and one
/// // it was never run over.
/// let scores = vec![0.1, 0.9, -0.6, f64::NAN];
/// let track = DynseqTrack::new(100, b"ACGT".to_vec(), scores).label("contribution");
///
/// assert_eq!(track.score_at(100), Some(0.1));
/// assert_eq!(track.score_at(103), None, "an unscored base became a number");
///
/// let svg = plot("chr1:101-104").expect("a region").add_track(track).to_svg();
/// assert!(svg.contains("contribution"));
/// ```
#[derive(Debug, Clone)]
pub struct DynseqTrack {
    start: u64,
    seq: Vec<u8>,
    /// A score per base, with a non-finite value meaning the base was not
    /// scored. Private, and reached only through [`DynseqTrack::score_at`],
    /// which answers with `None`, the way `MatrixRow::value` already does for
    /// a sample nobody typed.
    scores: Vec<f64>,
    label: Option<String>,
    height: f64,
    symmetric: bool,
    max_extent: Option<f64>,
    axis: QuantitativeAxis,
    letter_threshold: f64,
    show_scale: bool,
}

impl DynseqTrack {
    /// `seq[i]` and `scores[i]` describe the base at `start + i`, 0-based.
    ///
    /// A score that is not a number is a base that was not scored, and it is
    /// drawn as one: no glyph, and no rule under it either.
    pub fn new(start: u64, seq: impl Into<Vec<u8>>, scores: impl Into<Vec<f64>>) -> Self {
        DynseqTrack {
            start,
            seq: seq.into(),
            scores: scores.into(),
            label: None,
            height: 70.0,
            symmetric: true,
            max_extent: None,
            axis: QuantitativeAxis::new(),
            letter_threshold: 5.0,
            show_scale: true,
        }
    }

    /// Lays scattered `(position, score)` pairs over a sequence.
    ///
    /// A base no pair mentions stays unscored. It does not become a nought,
    /// which is the difference between a model that looked and found nothing
    /// and a model that never looked.
    pub fn from_pairs(
        start: u64,
        seq: impl Into<Vec<u8>>,
        pairs: impl IntoIterator<Item = (u64, f64)>,
    ) -> Self {
        Self::from_spans(
            start,
            seq,
            pairs.into_iter().map(|(pos, score)| (pos, pos + 1, score)),
        )
    }

    /// Lays half-open `(start, end, score)` spans over a sequence.
    ///
    /// What a bedGraph states, taken as it states it. A row covering a hundred
    /// bases is one span here and not a hundred pairs, which is the difference
    /// between a file of a few rows staying small and a file of a few rows
    /// across a whole window becoming gigabytes.
    ///
    /// A base no span covers stays unscored. It does not become a nought,
    /// which is the difference between a model that looked and found nothing
    /// and a model that never looked.
    pub fn from_spans(
        start: u64,
        seq: impl Into<Vec<u8>>,
        spans: impl IntoIterator<Item = (u64, u64, f64)>,
    ) -> Self {
        let seq = seq.into();
        let mut scores = vec![f64::NAN; seq.len()];
        for (from, to, score) in spans {
            let first = from.saturating_sub(start);
            let last = to.saturating_sub(start);
            let Ok(first) = usize::try_from(first) else {
                continue;
            };
            let last = usize::try_from(last)
                .unwrap_or(scores.len())
                .min(scores.len());
            for slot in scores.iter_mut().take(last).skip(first) {
                *slot = score;
            }
        }
        Self::new(start, seq, scores)
    }

    /// Sets the name in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, pixels: f64) -> Self {
        self.height = finite_within(pixels, 4.0, 400.0, 70.0);
        self
    }

    /// How wide a base has to be, in pixels, before its letter is drawn.
    pub fn letter_threshold(mut self, pixels: f64) -> Self {
        self.letter_threshold = finite_within(pixels, 0.0, 200.0, 5.0);
        self
    }

    /// Draws or hides the score axis in the left strip.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Puts the rule in the middle and reaches the same distance either way.
    ///
    /// The default, because the sign is the claim: a scan whose positive half
    /// is drawn twice the scale of its negative half says the pushes were
    /// bigger than the pulls, and it says so even when they were not.
    pub fn symmetric(mut self, yes: bool) -> Self {
        self.symmetric = yes;
        self
    }

    /// Pins how far the axis reaches either side, in the model's own units.
    ///
    /// Two panels of the same height are two different rulers, and nothing on
    /// the page says so. Pin both to compare them.
    pub fn max_extent(mut self, extent: f64) -> Self {
        if extent.is_finite() && extent > 0.0 {
            self.max_extent = Some(extent);
        }
        self
    }

    /// Replaces the value axis.
    pub fn axis(mut self, axis: QuantitativeAxis) -> Self {
        self.axis = axis;
        self
    }

    /// The score at a 0-based position, or `None` where that base is unscored.
    pub fn score_at(&self, pos: u64) -> Option<f64> {
        let index = usize::try_from(pos.checked_sub(self.start)?).ok()?;
        self.scores.get(index).copied().filter(|s| s.is_finite())
    }

    /// The reference base at a 0-based position.
    pub fn base_at(&self, pos: u64) -> Option<u8> {
        let index = usize::try_from(pos.checked_sub(self.start)?).ok()?;
        self.seq.get(index).copied()
    }

    /// The first base this track carries.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// How many bases it carries.
    pub fn len(&self) -> usize {
        self.seq.len().max(self.scores.len())
    }

    /// Whether it carries none.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The largest score either way inside `region`, or `None` where nothing in
    /// it was scored.
    pub fn visible_extent(&self, region: &Region) -> Option<f64> {
        let mut widest = 0.0f64;
        let mut any = false;
        for pos in self.window(region) {
            if let Some(score) = self.score_at(pos) {
                widest = widest.max(score.abs());
                any = true;
            }
        }
        any.then_some(widest.max(f64::MIN_POSITIVE))
    }

    /// How many bases inside `region` carry a letter and no score.
    pub fn unscored(&self, region: &Region) -> usize {
        self.window(region)
            .filter(|pos| self.score_at(*pos).is_none())
            .count()
    }

    /// The positions this track holds that fall inside `region`.
    ///
    /// A range and not a vector. A region may span 2^28 bases, and collecting
    /// one `u64` per base to walk it is two gigabytes to read a few thousand
    /// scores. The property net found that by being killed for it.
    fn window(&self, region: &Region) -> std::ops::Range<u64> {
        let first = self.start.max(region.start());
        let last = self
            .start
            .saturating_add(self.len() as u64)
            .min(region.end());
        first..last.max(first)
    }

    /// How far the axis reaches each way, resolved against the data.
    fn extent(&self, region: &Region) -> (f64, f64) {
        if let Some(pinned) = self.max_extent {
            return (-pinned, pinned);
        }
        let mut lowest = 0.0f64;
        let mut highest = 0.0f64;
        for pos in self.window(region) {
            if let Some(score) = self.score_at(pos) {
                lowest = lowest.min(score);
                highest = highest.max(score);
            }
        }
        if self.symmetric {
            let reach = lowest.abs().max(highest).max(f64::MIN_POSITIVE);
            return self.axis.resolve(-reach, reach);
        }
        self.axis.resolve(lowest.min(0.0), highest.max(0.0))
    }
}

impl Track for DynseqTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale {
            return 0.0;
        }
        // Measured over every label that will be drawn, not just the nought.
        // The clip covers the strip, so a wider label loses its left end, and
        // a score of minus twelve thousand printed with its minus sign cut off
        // is a positive number an order of magnitude too small.
        let size = theme.font_size - 1.0;
        // The widest score the track carries anywhere, not the widest inside
        // the window, because the figure asks for this width before it says
        // which window. A visible extent is never larger than that one.
        let reach = self.max_extent.unwrap_or_else(|| {
            self.scores
                .iter()
                .filter(|score| score.is_finite())
                .fold(0.0f64, |widest, score| widest.max(score.abs()))
                .max(1.0)
        });
        [reach, -reach, 0.0]
            .into_iter()
            .map(|value| crate::svg::text_width(&self.axis.label(value), size))
            .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let (lo, hi) = self.extent(ctx.region);
        let span = (hi - lo).max(f64::MIN_POSITIVE);
        let y_of = |score: f64| band.bottom() - ((score - lo) / span).clamp(0.0, 1.0) * band.h;
        let rule_y = y_of(0.0);
        let quiet = mix(&ctx.theme.foreground, ctx.theme.surface(), 0.45);
        let per_base = ctx.scale.px_per_bp();

        // No numbers beside a band nobody scored. The extent would be the
        // symmetric default, and printing it puts a quantitative scale on the
        // page for a measurement that was never made.
        let measured = self.visible_extent(ctx.region).is_some();
        if self.show_scale && measured && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            // Ordered before clamping. On a band shorter than the text is
            // tall the two bounds cross, and `clamp` panics rather than
            // choosing: a track four pixels high is allowed by `height`, so
            // this was reachable.
            let top = band.y + size * 0.8;
            let bottom = band.bottom() - size * 0.15;
            for value in [hi, 0.0, lo] {
                let y = (y_of(value) + size * 0.35).clamp(top.min(bottom), top.max(bottom));
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    y,
                    &self.axis.label(value),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }

        // One rule per run of scored bases, never one across the band. A base
        // scoring exactly nought draws no glyph, and so does a base nobody
        // scored, so the rule under it is the only thing that tells them apart.
        let mut run: Option<(u64, u64)> = None;
        let mut runs: Vec<(u64, u64)> = Vec::new();
        for pos in self.window(ctx.region) {
            match (self.score_at(pos).is_some(), &mut run) {
                (true, Some(open)) if open.1 == pos => open.1 = pos + 1,
                (true, slot) => {
                    if let Some(done) = slot.take() {
                        runs.push(done);
                    }
                    *slot = Some((pos, pos + 1));
                }
                (false, slot) => {
                    if let Some(done) = slot.take() {
                        runs.push(done);
                    }
                }
            }
        }
        if let Some(done) = run.take() {
            runs.push(done);
        }
        // Runs closer together than a pixel are joined before anything is
        // written. Every gap is a separate `<line>`, and a file scoring every
        // other base over a megabase produced half a million of them in a nine
        // hundred pixel band: forty-two megabytes of SVG to draw a rule.
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (from, to) in runs {
            match merged.last_mut() {
                Some(open) if ctx.scale.x(from) - ctx.scale.x(open.1) < 1.0 => open.1 = to,
                _ => merged.push((from, to)),
            }
        }
        for (from, to) in &merged {
            ctx.svg.line(
                ctx.scale.x(*from),
                rule_y,
                ctx.scale.x(*to),
                rule_y,
                &ctx.theme.rule,
                ctx.theme.tokens.stroke,
            );
        }

        if per_base >= 1.0 {
            let letters = per_base >= self.letter_threshold;
            for pos in self.window(ctx.region) {
                let Some(score) = self.score_at(pos) else {
                    continue;
                };
                let x = ctx.scale.x(pos);
                let width = (ctx.scale.x(pos + 1) - x).max(0.6);
                let y = y_of(score);
                let top = y.min(rule_y);
                let depth = (y - rule_y).abs();
                if depth <= 0.2 {
                    continue;
                }
                let base = self.base_at(pos).unwrap_or(b'N');
                let color = ctx.theme.bases.of(base).to_string();

                ctx.svg.begin_titled(&format!(
                    "{} at {}, {}",
                    base as char,
                    crate::track::axis::group_thousands(pos.saturating_add(1)),
                    text_rounded(score, 4)
                ));
                if letters {
                    // The glyph baseline is the far edge of the box either way,
                    // which is what `LogoTrack::draw_symbol` uses: a letter
                    // hanging below the rule then sits on the bottom of its own
                    // box and grows up towards the rule.
                    let cell = Rect {
                        x,
                        y: top,
                        w: width,
                        h: depth,
                    };
                    let font = cell.h / ctx.theme.cap_height_ratio.max(0.1);
                    ctx.svg.glyph(
                        cell.x,
                        cell.bottom(),
                        cell.w,
                        font,
                        &(base as char).to_string(),
                        &color,
                    );
                } else {
                    ctx.svg.rect(x, top, width, depth, &color);
                }
                ctx.svg.end_group();
            }
            return;
        }

        // Below a pixel a base, an envelope of the extremes in one neutral ink.
        // A column spanning forty bases has no base, and painting it green
        // would claim one.
        let columns = (band.w.ceil() as usize).max(1);
        let mut lowest = vec![f64::MAX; columns];
        let mut highest = vec![f64::MIN; columns];
        for pos in self.window(ctx.region) {
            let Some(score) = self.score_at(pos) else {
                continue;
            };
            let at = ctx.scale.x_center(pos) - band.x;
            if at < 0.0 {
                continue;
            }
            let index = (at as usize).min(columns - 1);
            lowest[index] = lowest[index].min(score);
            highest[index] = highest[index].max(score);
        }
        for (index, (low, high)) in lowest.iter().zip(&highest).enumerate() {
            if low > high {
                continue;
            }
            let top = y_of(high.max(0.0));
            let bottom = y_of(low.min(0.0));
            ctx.svg.rect(
                band.x + index as f64,
                top,
                1.0,
                (bottom - top).max(0.6),
                &quiet,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::theme::Theme;

    fn region(width: u64) -> Region {
        Region::new("chr1", 0, width).unwrap()
    }

    fn drawn(track: DynseqTrack, width: u64, pixels: f64) -> String {
        Figure::new(region(width))
            .width(pixels)
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    fn rules(svg: &str) -> usize {
        svg.matches("stroke=\"#d7dce2\"").count()
    }

    #[test]
    fn a_base_scoring_nought_and_a_base_nobody_scored_are_told_apart() {
        // Both draw no glyph, because a glyph of no height is no glyph. The
        // rule under one of them is the whole difference, and without it this
        // is the crate's named mistake.
        let zero = DynseqTrack::new(0, b"ACGT".to_vec(), vec![0.5, 0.0, 0.0, 0.5]);
        let absent = DynseqTrack::new(0, b"ACGT".to_vec(), vec![0.5, f64::NAN, f64::NAN, 0.5]);
        assert_eq!(zero.score_at(1), Some(0.0));
        assert_eq!(absent.score_at(1), None);
        assert_eq!(
            rules(&drawn(zero, 4, 400.0)),
            1,
            "a run of scored bases was broken by a nought"
        );
        assert_eq!(
            rules(&drawn(absent, 4, 400.0)),
            2,
            "the rule ran straight over bases nobody scored"
        );
    }

    #[test]
    fn the_zoom_picks_the_regime_and_nothing_else_does() {
        let scores: Vec<f64> = (0..400).map(|i| ((i % 9) as f64 - 4.0) / 4.0).collect();
        let seq: Vec<u8> = (0..400).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let track = || DynseqTrack::new(0, seq.clone(), scores.clone());

        // Wide: letters.
        let letters = drawn(track(), 20, 900.0);
        assert!(letters.contains("textLength"), "no letters where they fit");

        // Middle: boxes, the identity still true and merely illegible.
        let bars = drawn(track(), 300, 900.0);
        assert!(
            !bars.contains("textLength"),
            "letters drawn at a sliver wide"
        );
        assert!(bars.contains("<rect "), "no bars in the middle regime");

        // Narrow: one neutral envelope, never a base colour.
        let envelope = drawn(track(), 400, 300.0);
        assert!(!envelope.contains("textLength"));
        for base in ["#3d9970", "#e6194b"] {
            assert!(
                !envelope.contains(base),
                "a column spanning many bases was painted as one of them"
            );
        }
    }

    #[test]
    fn the_envelope_shows_both_directions_rather_than_an_average() {
        // A +2 and a -2 in one column average to a nought, which says the model
        // ignored the place. Both extremes or nothing.
        let seq: Vec<u8> = (0..400).map(|_| b'A').collect();
        let scores: Vec<f64> = (0..400)
            .map(|i| if i % 2 == 0 { 2.0 } else { -2.0 })
            .collect();
        let svg = drawn(
            DynseqTrack::new(0, seq, scores).label("both ways"),
            400,
            300.0,
        );
        // Every column holds both, so every envelope rectangle has to cross the
        // rule rather than sit on it.
        let tall = svg
            .split("<rect ")
            .skip(1)
            .filter_map(|piece| piece.split("/>").next())
            .filter(|head| head.contains("width=\"1\""))
            .count();
        assert!(tall > 10, "the envelope collapsed: {tall} columns drawn");
    }

    #[test]
    fn the_two_sides_reach_the_same_distance_by_default() {
        let track = DynseqTrack::new(0, b"AC".to_vec(), vec![3.0, -0.5]);
        let (lo, hi) = track.extent(&region(2));
        assert!((lo + hi).abs() < 1e-9, "{lo} and {hi} are not a mirror");

        let free = track.clone().symmetric(false);
        let (lo, hi) = free.extent(&region(2));
        assert!(hi > 0.0 && lo < 0.0 && (lo + hi).abs() > 1e-9);

        let pinned = DynseqTrack::new(0, b"AC".to_vec(), vec![3.0, -0.5]).max_extent(10.0);
        assert_eq!(pinned.extent(&region(2)), (-10.0, 10.0));
    }

    #[test]
    fn a_position_no_pair_mentions_stays_unscored() {
        // Not nought. The difference is a model that looked and found nothing
        // against one that never looked.
        let track = DynseqTrack::from_pairs(100, b"ACGT".to_vec(), [(100, 0.4), (102, -0.2)]);
        assert_eq!(track.score_at(100), Some(0.4));
        assert_eq!(track.score_at(101), None);
        assert_eq!(track.score_at(102), Some(-0.2));
        assert_eq!(track.unscored(&Region::new("chr1", 100, 104).unwrap()), 2);
    }

    #[test]
    fn a_score_that_is_not_a_number_never_reaches_the_figure() {
        let track = DynseqTrack::new(
            0,
            b"ACGT".to_vec(),
            vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.5],
        );
        let svg = drawn(track, 4, 400.0);
        for spelling in ["NaN", "nan", "inf", "Infinity"] {
            assert!(!svg.contains(spelling), "{spelling} reached the output");
        }
    }

    #[test]
    fn a_track_with_nothing_scored_in_it_draws_no_rule_at_all() {
        // Rather than a rule across the band, which would say the model was run
        // here and found nothing.
        let track = DynseqTrack::new(0, b"ACGT".to_vec(), vec![f64::NAN; 4]);
        assert_eq!(track.visible_extent(&region(4)), None);
        assert_eq!(rules(&drawn(track, 4, 400.0)), 0);
    }

    #[test]
    fn a_band_shorter_than_its_own_labels_does_not_panic() {
        // The two bounds of the label clamp cross once the band is shorter than
        // the text is tall, and `clamp` panics rather than choosing. `height`
        // allows four pixels, so this was reachable without asking for it.
        for pixels in [4.0, 8.0, 10.0, 10.4, 10.5, 40.0] {
            let track = DynseqTrack::new(0, b"ACGT".to_vec(), vec![0.5; 4]).height(pixels);
            let svg = Figure::new(region(4))
                .show_region_label(false)
                .push(track)
                .to_svg();
            assert!(svg.starts_with("<svg "), "{pixels} px");
        }
    }

    #[test]
    fn the_strip_is_wide_enough_for_the_labels_that_go_in_it() {
        // Sized from the label of nought alone, the extremes were cut off on
        // the left by the clip, and a score of minus twelve thousand printed
        // without its minus sign is a positive number an order of magnitude
        // too small.
        let big = DynseqTrack::new(0, b"AC".to_vec(), vec![12_345.678, -12_345.678]);
        let small = DynseqTrack::new(0, b"AC".to_vec(), vec![0.5, -0.5]);
        let theme = Theme::light();
        assert!(
            big.y_axis_width(&theme) > small.y_axis_width(&theme),
            "the strip is the same width whatever the numbers are"
        );
        let size = theme.font_size - 1.0;
        let widest = crate::svg::text_width(&big.axis.label(-12_345.678), size);
        assert!(
            big.y_axis_width(&theme) >= widest,
            "the widest label will not fit"
        );
    }

    #[test]
    fn the_rule_is_bounded_by_the_pixels_and_not_by_the_bases() {
        // Every other base scored is half a million runs over a megabase, and
        // one line each was forty-two megabytes of SVG to draw a rule.
        let bases: Vec<u8> = (0..40_000).map(|_| b'A').collect();
        let scores: Vec<f64> = (0..40_000)
            .map(|i| if i % 2 == 0 { 0.5 } else { f64::NAN })
            .collect();
        let svg = drawn(DynseqTrack::new(0, bases, scores), 40_000, 900.0);
        let lines = svg.matches("<line ").count();
        assert!(lines < 200, "{lines} rule segments in a 900 pixel band");
    }

    #[test]
    fn a_letter_is_coloured_by_its_base_and_named_in_a_tooltip() {
        let svg = drawn(
            DynseqTrack::new(10, b"ACGT".to_vec(), vec![0.5, 0.4, 0.3, 0.2]),
            20,
            600.0,
        );
        assert!(svg.contains(">A</text>"), "no letter drawn");
        // 1-based in the tooltip, the way the ruler counts.
        assert!(svg.contains("A at 11,"), "{svg}");
    }
}
