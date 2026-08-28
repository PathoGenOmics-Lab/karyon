//! Copy number along the sequence, as the segments a caller fitted.
//!
//! A segmentation is a set of intervals, each carrying how many copies were
//! called over it, and the good ones carry two numbers rather than one: how
//! many copies in total, and how many of them came from the quieter of the two
//! alleles. [`CopyNumberTrack`] draws both, because the relation between them
//! is the finding.
//!
//! Copies are continuous here rather than counted. A subclonal call and a
//! purity-adjusted one are both fractional, and a ploidy of 2.6 is a number a
//! caller fits. The rungs of the ladder are at whole copies because whole
//! copies are where the interpretable states are, not because a call has to
//! land on one.
//!
//! # Why this is not a window track
//!
//! [`WindowTrack`](crate::WindowTrack) can hold two numbers over one span,
//! since nothing says one row per place, so the obvious objection to a second
//! type does not hold on the count of numbers. It holds on what the two tracks
//! draw where nothing happened.
//!
//! A window track fills from its baseline out to the value, so it draws a
//! window only when `hi` is above the line or `lo` is below it. A segment
//! called at exactly the ploidy is neither, and it draws nothing. That is not
//! an edge case here: a balanced segment is most of a genome, and a copy number
//! figure in which the quiet arms are blank is a figure in which the quiet arms
//! and the arms nobody called look the same. Measured, on four windows at two,
//! no call, seven and nought copies with the baseline at two: two marks come
//! out, and the two that vanish are the balanced one and the missing one.
//!
//! So a level here is a bar drawn at the level, not a fill from a line to it.
//! Every called segment is a mark, and a segment nobody called is the only
//! blank.
//!
//! The second half is loss of heterozygosity. A minor allele of nought with
//! copies still present is a finding, and it must not look like the absence of
//! one, so the absence is not representable: [`CopyNumber::Total`] has no field
//! to put a minor allele in, and [`CopyNumber::minor`] answers with `None`.
//! Writing `minor().unwrap_or(0.0)` anywhere is this bug put back.
//!
//! The third is the edges. A window's edges are a grid the analyst chose, so
//! rounding them to whole pixel columns loses nothing. A segment's edges are
//! breakpoints somebody inferred, so they are drawn where they were reported.
//!
//! # Where balanced sits is the caller's to say
//!
//! [`CopyNumberTrack::at_ploidy`] takes the number and has no default, because
//! this crate does not know what it is looking at. Two is right for a human
//! autosome, one for most bacteria, and neither for a polyploid or for a tumour
//! whose ploidy the caller fitted. A rule in the wrong place does not merely
//! mis-scale the figure: it swaps gain for loss, and says so confidently.
//!
//! # What a pixel column throws away, and what is done about it
//!
//! No average is offered, and that is deliberate. The mean of a call of one and
//! a call of three is two, which is a level nobody called, and on a diploid
//! ladder it lands exactly on the rule that means unchanged. So every segment
//! is drawn at its own level, at a pixel wide where it is narrower than one,
//! and a hairline riser joins the lowest and the highest level any segment
//! reached in each column. A focal amplification inside a quiet arm is then a
//! one pixel bar with a stalk down to the arm, rather than a bump in an
//! average. Where two segments land in the same column both are drawn, in the
//! order they were given, so the picture is a function of the file.
//!
//! Nothing is ever drawn between two segments. A caller reports intervals, not
//! a polyline, and a line across a gap asserts a breakpoint at a coordinate
//! nobody reported.

use crate::scale::Scale;
use crate::style::LinePattern;
use crate::svg::{finite_within, text_rounded, Anchor};
use crate::theme::{mix, Theme};
use crate::track::feature::span_label;
use crate::track::legend::Legend;
use crate::track::{DrawContext, Track};

/// How many copies were called, and whether the alleles were told apart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CopyNumber {
    /// Total copies, with nothing said about the two alleles.
    Total(f64),
    /// An allele-specific call, `major` being the larger of the pair.
    Allelic {
        /// Copies of the commoner allele.
        major: f64,
        /// Copies of the rarer one.
        minor: f64,
    },
}

impl CopyNumber {
    /// An allele-specific call, ordered, so a caller holding the pair the other
    /// way round is right anyway.
    ///
    /// A pair with a value that is not a number is not a call at all, and in
    /// particular it is not a total. A caller reporting three major copies and
    /// `NA` for the minor has not said there are three copies, it has said the
    /// total is three plus something nobody measured. Keeping the surviving
    /// number as a total would put that missing something at nought, which is
    /// this crate's named mistake with an extra step in front of it.
    pub fn allelic(a: f64, b: f64) -> CopyNumber {
        if !a.is_finite() || !b.is_finite() {
            return CopyNumber::Total(f64::NAN);
        }
        CopyNumber::Allelic {
            major: a.max(b),
            minor: a.min(b),
        }
    }

    /// Total copies: the number itself, or the two added.
    pub fn total(self) -> f64 {
        match self {
            CopyNumber::Total(total) => total,
            CopyNumber::Allelic { major, minor } => major + minor,
        }
    }

    /// Copies of the rarer allele, or `None` where the alleles were not told
    /// apart.
    ///
    /// Never nought standing in for a call nobody made. That substitution turns
    /// "we did not look" into loss of heterozygosity, which is the one mistake
    /// this whole module is arranged around.
    pub fn minor(self) -> Option<f64> {
        match self {
            CopyNumber::Total(_) => None,
            CopyNumber::Allelic { minor, .. } => Some(minor),
        }
    }

    /// Copies of the commoner allele, or `None` where they were not told apart.
    pub fn major(self) -> Option<f64> {
        match self {
            CopyNumber::Total(_) => None,
            CopyNumber::Allelic { major, .. } => Some(major),
        }
    }

    /// Whether every number in the call is one a figure can be drawn from.
    pub fn is_called(self) -> bool {
        match self {
            CopyNumber::Total(total) => total.is_finite(),
            CopyNumber::Allelic { major, minor } => major.is_finite() && minor.is_finite(),
        }
    }
}

/// One interval a caller fitted, and what it called over it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CopyNumberSegment {
    /// First base, 0-based.
    pub start: u64,
    /// One past the last base.
    pub end: u64,
    /// What was called over it.
    pub copy: CopyNumber,
}

impl CopyNumberSegment {
    /// A segment with a total and no allele split.
    pub fn total(start: u64, end: u64, total: f64) -> Self {
        CopyNumberSegment {
            start: start.min(end),
            end: start.max(end),
            copy: CopyNumber::Total(total),
        }
    }

    /// A segment with an allele-specific call.
    pub fn allelic(start: u64, end: u64, a: f64, b: f64) -> Self {
        CopyNumberSegment {
            start: start.min(end),
            end: start.max(end),
            copy: CopyNumber::allelic(a, b),
        }
    }

    /// How many bases it covers.
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Whether it covers no bases.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether heterozygosity was lost here, or `None` where the alleles were
    /// not told apart.
    ///
    /// `None` rather than `false`, so that counting the segments that lost it
    /// cannot quietly count the ones nobody looked at as having kept it.
    pub fn loh(&self) -> Option<bool> {
        let minor = self.copy.minor()?;
        Some(minor == 0.0 && self.copy.total() > 0.0)
    }
}

/// The most copies a ladder is drawn to.
///
/// A caller does not report a hundred thousand copies, and a file that says so
/// says it by mistake. The number is a bound on the drawing rather than a claim
/// about biology: past it the rungs are closer than a pixel and the ladder
/// stops being one.
const CEILING: f64 = 100_000.0;

/// The lane along the foot of the band, saying what the alleles did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Allele {
    /// The alleles were not told apart here.
    Unresolved,
    /// Both alleles are present.
    Both,
    /// One allele is gone and copies remain.
    Lost,
    /// No copies at all.
    Absent,
}

/// Segmented copy number, drawn on a ladder of whole copies.
///
/// ```
/// use karyon::{plot, CopyNumberSegment, CopyNumberTrack};
///
/// let segments = vec![
///     CopyNumberSegment::allelic(0, 40_000, 1.0, 1.0),
///     // Two copies, both from one allele: the total bar sits exactly on the
///     // balanced rule, and the lane along the bottom is what says the other
///     // allele went.
///     CopyNumberSegment::allelic(40_000, 90_000, 2.0, 0.0),
///     CopyNumberSegment::total(90_000, 120_000, 7.0),
/// ];
///
/// let svg = plot("chr8:1-120,000")
///     .expect("a region")
///     .add_copy_number(segments, 2.0)
///     .label("copy number")
///     .to_svg();
///
/// assert!(svg.contains("heterozygosity lost"));
/// assert!(svg.contains("allele split not called"));
/// ```
#[derive(Debug, Clone)]
pub struct CopyNumberTrack {
    segments: Vec<CopyNumberSegment>,
    ploidy: f64,
    label: Option<String>,
    height: f64,
    cap: Option<f64>,
    gain: Option<String>,
    loss: Option<String>,
    neutral: Option<String>,
    loh: Option<String>,
    show_scale: bool,
    show_alleles: bool,
}

impl CopyNumberTrack {
    /// A track over `segments`, with `ploidy` copies as the balanced state.
    ///
    /// There is no default for the ploidy and there is not going to be one.
    /// This crate is agnostic about what it is drawing, so it cannot know
    /// whether two copies is normal, and a rule in the wrong place turns every
    /// gain into a loss without saying anything.
    pub fn at_ploidy(segments: impl Into<Vec<CopyNumberSegment>>, ploidy: f64) -> Self {
        CopyNumberTrack {
            segments: segments.into(),
            ploidy: if ploidy.is_finite() && ploidy >= 0.0 {
                ploidy
            } else {
                2.0
            },
            label: None,
            height: 74.0,
            cap: None,
            gain: None,
            loss: None,
            neutral: None,
            loh: None,
            show_scale: true,
            show_alleles: true,
        }
    }

    /// Two copies as the balanced state.
    pub fn diploid(segments: impl Into<Vec<CopyNumberSegment>>) -> Self {
        Self::at_ploidy(segments, 2.0)
    }

    /// One copy as the balanced state.
    pub fn haploid(segments: impl Into<Vec<CopyNumberSegment>>) -> Self {
        Self::at_ploidy(segments, 1.0)
    }

    /// Sets the name in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        if height.is_finite() {
            self.height = height.max(24.0);
        }
        self
    }

    /// Pins the top rung of the ladder, taken literally.
    ///
    /// For comparing two figures whose amplifications reach different heights,
    /// where letting each pick its own ceiling makes the quieter one look as
    /// dramatic as the louder.
    pub fn cap(mut self, copies: f64) -> Self {
        if copies.is_finite() && copies > 0.0 {
            self.cap = Some(copies);
        }
        self
    }

    /// Sets the inks for gain, loss and the balanced state.
    pub fn colors(
        mut self,
        gain: impl Into<String>,
        loss: impl Into<String>,
        neutral: impl Into<String>,
    ) -> Self {
        self.gain = Some(gain.into());
        self.loss = Some(loss.into());
        self.neutral = Some(neutral.into());
        self
    }

    /// Sets the ink of the lost-heterozygosity mark.
    pub fn loh_color(mut self, color: impl Into<String>) -> Self {
        self.loh = Some(color.into());
        self
    }

    /// Draws or hides the copy scale in the left strip.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Draws or hides the lane along the foot of the band.
    ///
    /// Worth hiding only where no segment carries an allele split. Copy-neutral
    /// loss of heterozygosity puts the total bar exactly on the balanced rule,
    /// so what says the allele went is the minor bar down at nought copies and
    /// this lane, and the lane is the one a reader finds without looking for
    /// it.
    pub fn show_alleles(mut self, show: bool) -> Self {
        self.show_alleles = show;
        self
    }

    /// The segments, in the order they were given.
    pub fn segments(&self) -> &[CopyNumberSegment] {
        &self.segments
    }

    /// Where balanced sits.
    pub fn ploidy(&self) -> f64 {
        self.ploidy
    }

    /// The top rung of the ladder.
    ///
    /// At least two copies above the ploidy even for a flat genome, so that a
    /// figure with nothing gained in it still has the room to show that.
    pub fn ceiling(&self) -> f64 {
        if let Some(cap) = self.cap {
            return cap;
        }
        let highest = self
            .segments
            .iter()
            .filter(|segment| segment.copy.is_called())
            .map(|segment| segment.copy.total())
            .filter(|total| total.is_finite())
            .fold(0.0f64, f64::max);
        // Bounded, and not only for tidiness. A ceiling of infinity makes the
        // spacing between rungs nought, and the loop that walks the ladder to
        // the top then never reaches it. The property net found that by being
        // killed for the memory.
        let asked = (highest * 1.06).max(self.ploidy + 2.0).ceil();
        finite_within(asked, 1.0, CEILING, self.ploidy + 2.0)
    }

    /// The spans where one allele is gone and copies remain.
    pub fn loh_spans(&self) -> Vec<(u64, u64)> {
        self.segments
            .iter()
            .filter(|segment| segment.loh() == Some(true))
            .map(|segment| (segment.start, segment.end))
            .collect()
    }

    /// The spans with no copies at all.
    pub fn homozygous_deletions(&self) -> Vec<(u64, u64)> {
        self.segments
            .iter()
            .filter(|segment| segment.copy.is_called() && segment.copy.total() == 0.0)
            .map(|segment| (segment.start, segment.end))
            .collect()
    }

    /// How many segments carry a total and no allele split.
    ///
    /// The number that says whether the lane along the bottom is mostly a
    /// finding or mostly a gap in the file.
    pub fn without_allele_call(&self) -> usize {
        self.segments
            .iter()
            .filter(|segment| segment.copy.minor().is_none())
            .count()
    }

    /// A key holding only the marks this data actually used.
    pub fn legend(&self, theme: &Theme) -> Legend {
        // The theme is taken rather than the colours written down, because the
        // inks come from it unless the caller overrode them, and a key holding
        // literals names the wrong colour the first time a figure is drawn
        // dark.
        let gain = self
            .gain
            .clone()
            .unwrap_or_else(|| theme.color(1).to_string());
        let loss = self
            .loss
            .clone()
            .unwrap_or_else(|| theme.color(0).to_string());
        let lost = self
            .loh
            .clone()
            .unwrap_or_else(|| theme.color(2).to_string());

        let called = |f: fn(&CopyNumberSegment, f64) -> bool| {
            self.segments
                .iter()
                .any(|segment| segment.copy.is_called() && f(segment, self.ploidy))
        };
        let has = |state: Allele| self.segments.iter().any(|s| self.state(s) == state);

        let mut legend = Legend::new();
        if called(|s, ploidy| s.copy.total() > ploidy) {
            legend = legend.key("gain", gain);
        }
        if called(|s, ploidy| s.copy.total() < ploidy) {
            legend = legend.key("loss", loss.clone());
        }
        if has(Allele::Lost) {
            legend = legend.key("heterozygosity lost", lost);
        }
        if has(Allele::Absent) {
            legend = legend.key("no copies", loss);
        }
        if has(Allele::Unresolved) {
            legend = legend.key("allele split not called", theme.rule.clone());
        }
        legend
    }

    /// What the lane says over one segment.
    fn state(&self, segment: &CopyNumberSegment) -> Allele {
        if !segment.copy.is_called() {
            return Allele::Unresolved;
        }
        match segment.copy.minor() {
            None => Allele::Unresolved,
            Some(_) if segment.copy.total() == 0.0 => Allele::Absent,
            Some(0.0) => Allele::Lost,
            Some(_) => Allele::Both,
        }
    }

    /// What a pointer hovering one segment is told.
    fn tooltip(&self, segment: &CopyNumberSegment) -> String {
        let where_it_is = span_label(segment.start, segment.end);
        if !segment.copy.is_called() {
            return format!("{where_it_is}, no call");
        }
        let total = text_rounded(segment.copy.total(), 2);
        match (segment.copy.major(), segment.copy.minor()) {
            (Some(major), Some(minor)) => {
                let split = format!("{} + {}", text_rounded(major, 2), text_rounded(minor, 2));
                if minor == 0.0 && segment.copy.total() > 0.0 {
                    format!("{where_it_is}, {total} copies ({split}), heterozygosity lost")
                } else {
                    format!("{where_it_is}, {total} copies ({split})")
                }
            }
            _ => format!("{where_it_is}, {total} copies, allele split not called"),
        }
    }

    /// The rungs of the ladder, thinned until they have room to be read.
    fn rungs(&self, band_height: f64) -> Vec<f64> {
        let ceiling = self.ceiling();
        let spacing = band_height / ceiling.max(1.0);
        // One, two, five, ten and then the same run again a decade up, rather
        // than stopping at ten: a ceiling of a thousand copies with a step of
        // ten is a rung every two thirds of a pixel, which paints the whole
        // band in the rule colour and is no longer a ladder.
        let mut step = 1.0f64;
        while spacing * step < 9.0 && step < ceiling {
            let decade = 10f64.powf(step.log10().floor());
            step = match (step / decade).round() as u32 {
                1 => decade * 2.0,
                2 => decade * 5.0,
                _ => decade * 10.0,
            };
        }
        let mut rungs = Vec::new();
        let mut copies = 0.0;
        // Counted rather than compared, so no ceiling and no step can make this
        // a loop that does not end.
        let steps = ((ceiling / step).floor() as usize).min(4_096);
        for _ in 0..=steps {
            rungs.push(copies);
            copies += step;
        }
        rungs
    }
}

impl Track for CopyNumberTrack {
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
        let widest = self
            .rungs(self.height)
            .into_iter()
            .chain(std::iter::once(self.ploidy))
            .map(|copies| crate::svg::text_width(&text_rounded(copies, 2), theme.font_size - 1.0))
            .fold(0.0f64, f64::max);
        widest + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let lane = if self.show_alleles { ctx.px(5.0) } else { 0.0 };
        let floor = band.bottom() - lane - if lane > 0.0 { ctx.px(2.0) } else { 0.0 };
        let ladder = (floor - band.y).max(1.0);
        let ceiling = self.ceiling();
        let y_of = |copies: f64| floor - (copies / ceiling).clamp(0.0, 1.0) * ladder;

        let gain = self
            .gain
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());
        let loss = self
            .loss
            .clone()
            .unwrap_or_else(|| ctx.theme.color(0).to_string());
        let neutral = self
            .neutral
            .clone()
            .unwrap_or_else(|| mix(&ctx.theme.foreground, ctx.theme.surface(), 0.45));
        let lost = self
            .loh
            .clone()
            .unwrap_or_else(|| ctx.theme.color(2).to_string());

        // Whole copies only. An evenly divided axis would print three and a
        // half copies, which is a number of copies nobody has.
        for copies in self.rungs(ladder) {
            let y = y_of(copies);
            if (copies - self.ploidy).abs() < f64::EPSILON {
                continue;
            }
            ctx.svg.line(
                band.x,
                y,
                band.right(),
                y,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
        }

        // Drawn before any data test and even with nothing to draw, because it
        // is what the band means: an empty track showing its rule says nothing
        // was called here, and an empty rectangle says nothing at all.
        let balanced = y_of(self.ploidy);
        ctx.svg.line(
            band.x,
            balanced,
            band.right(),
            balanced,
            &ctx.theme.rule,
            ctx.theme.tokens.stroke,
        );

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            // The rungs may be closer together than a label is tall, and a
            // column of numbers touching each other is a column nobody reads.
            // The lines stay where they are; only the labels thin out.
            let rungs = self.rungs(ladder);
            let apart = rungs
                .windows(2)
                .map(|pair| (y_of(pair[0]) - y_of(pair[1])).abs())
                .fold(f64::MAX, f64::min);
            let every = if apart >= size + 3.0 {
                1
            } else {
                (((size + 3.0) / apart.max(0.1)).ceil() as usize).max(1)
            };
            for copies in rungs.into_iter().step_by(every) {
                // Nudged down where the top rung sits on the edge of the band.
                // The clip covers the strip along with the band, so a label
                // whose ascenders reach above `band.y` loses them.
                // Held inside the band at both ends. With the allele lane
                // hidden the floor is the bottom of the band, and the nought
                // label then sat below it and lost its lower third to the clip.
                let top = band.y + size * 0.8;
                let bottom = band.bottom() - size * 0.15;
                let y = (y_of(copies) + size * 0.35).clamp(top.min(bottom), top.max(bottom));
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    y,
                    &text_rounded(copies, 2),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }

        let visible: Vec<&CopyNumberSegment> = self
            .segments
            .iter()
            .filter(|segment| {
                segment.end > ctx.region.start()
                    && segment.start < ctx.region.end()
                    && !segment.is_empty()
            })
            .collect();

        // One entry per pixel column, holding the lowest and the highest level
        // any segment reached over it. The riser between the two is what keeps
        // a focal event from being swallowed by the arm it sits in, without
        // averaging the two into a level nobody called.
        let columns = (band.w.ceil() as usize).max(1);
        let mut lo = vec![f64::MAX; columns];
        let mut hi = vec![f64::MIN; columns];
        let mut touched = false;

        for segment in &visible {
            if !segment.copy.is_called() {
                continue;
            }
            let total = segment.copy.total();
            let x0 = ctx.scale.x_at(segment.start as f64).max(band.x);
            let x1 = ctx.scale.x_at(segment.end as f64).min(band.right());
            if x1 <= x0 {
                continue;
            }
            let first = ((x0 - band.x).floor().max(0.0)) as usize;
            let last = ((x1 - band.x).ceil().min(columns as f64)) as usize;
            for slot in lo.iter_mut().zip(hi.iter_mut()).take(last).skip(first) {
                *slot.0 = slot.0.min(total);
                *slot.1 = slot.1.max(total);
                touched = true;
            }
        }

        if touched {
            for (index, (low, high)) in lo.iter().zip(&hi).enumerate() {
                if *low >= *high {
                    continue;
                }
                let x = band.x + index as f64 + 0.5;
                ctx.svg.line(
                    x,
                    y_of(*low),
                    x,
                    y_of(*high),
                    &neutral,
                    ctx.theme.tokens.hairline,
                );
            }
        }

        let bar = ctx.px(3.0);
        let thin = ctx.px(2.0);

        for segment in &visible {
            let x0 = ctx.scale.x_at(segment.start as f64).max(band.x);
            let x1 = ctx.scale.x_at(segment.end as f64).min(band.right());
            let width = (x1 - x0).max(1.0);

            ctx.svg.begin_titled(&self.tooltip(segment));
            if segment.copy.is_called() {
                let total = segment.copy.total();
                let y = y_of(total);
                let ink = if total > self.ploidy {
                    &gain
                } else if total < self.ploidy {
                    &loss
                } else {
                    &neutral
                };

                // The wash is what makes a departure readable at a glance. The
                // bar alone is a thin line, and a filled column up from the
                // floor would bury the minor allele bar inside it.
                let top = y.min(balanced);
                let depth = (y - balanced).abs();
                if depth > 0.5 {
                    ctx.svg
                        .rect_opacity(x0, top, width, depth, ink, ctx.theme.tokens.area_opacity);
                }
                ctx.svg.rect(x0, y - bar / 2.0, width, bar, ink);

                if let Some(minor) = segment.copy.minor() {
                    let my = y_of(minor);
                    ctx.svg.rect(x0, my - thin / 2.0, width, thin, &neutral);
                }
            }
            ctx.svg.end_group();
        }

        if lane > 0.0 {
            let top = band.bottom() - lane;
            for segment in &visible {
                let x0 = ctx.scale.x_at(segment.start as f64).max(band.x);
                let x1 = ctx.scale.x_at(segment.end as f64).min(band.right());
                let width = (x1 - x0).max(1.0);
                match self.state(segment) {
                    // Dotted, and a rule rather than a block, so that a stretch
                    // nobody resolved cannot be read as a state that was.
                    Allele::Unresolved => ctx.svg.line_pattern(
                        x0,
                        top + lane / 2.0,
                        x0 + width,
                        top + lane / 2.0,
                        &ctx.theme.rule,
                        ctx.theme.tokens.hairline,
                        LinePattern::Dotted,
                    ),
                    Allele::Both => ctx.svg.line(
                        x0,
                        top + lane / 2.0,
                        x0 + width,
                        top + lane / 2.0,
                        &neutral,
                        ctx.theme.tokens.hairline,
                    ),
                    Allele::Lost => ctx.svg.rect(x0, top, width, lane, &lost),
                    Allele::Absent => ctx.svg.rect(x0, top, width, lane, &loss),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;
    use crate::theme::Theme;

    fn region() -> Region {
        Region::new("chr8", 0, 10_000).unwrap()
    }

    fn drawn(track: CopyNumberTrack) -> String {
        Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    #[test]
    fn a_segment_at_the_ploidy_is_a_mark_and_a_segment_nobody_called_is_not() {
        // The whole reason this is not a window track. There, a value equal to
        // the baseline draws nothing, which makes the quiet arms of a genome
        // pixel-identical to the arms nobody called.
        let called = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, 2.0)]);
        let missing = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, f64::NAN)]);
        let with = drawn(called).matches("<rect ").count();
        let without = drawn(missing).matches("<rect ").count();
        assert!(
            with > without,
            "a balanced segment drew as little as an uncalled one: {with} against {without}"
        );
    }

    #[test]
    fn a_minor_allele_of_nought_is_not_the_same_as_no_allele_call() {
        let lost = CopyNumberSegment::allelic(0, 5_000, 2.0, 0.0);
        let unresolved = CopyNumberSegment::total(0, 5_000, 2.0);
        assert_eq!(lost.loh(), Some(true));
        assert_eq!(unresolved.loh(), None, "an absent call answered a question");
        assert_eq!(lost.copy.minor(), Some(0.0));
        assert_eq!(unresolved.copy.minor(), None);

        let one = drawn(CopyNumberTrack::diploid(vec![lost]));
        let other = drawn(CopyNumberTrack::diploid(vec![unresolved]));
        assert!(one.contains("heterozygosity lost"), "{one}");
        assert!(other.contains("allele split not called"), "{other}");
        assert_ne!(
            one.matches("<rect ").count(),
            other.matches("<rect ").count(),
            "the two drew the same marks"
        );
    }

    #[test]
    fn copy_neutral_loss_of_heterozygosity_is_only_in_the_lane() {
        // Two copies, both from one allele. The ladder puts the bar exactly on
        // the balanced rule, where a normal diploid arm sits, so the lane along
        // the foot is the one thing in the figure that says the allele went.
        let track = CopyNumberTrack::diploid(vec![CopyNumberSegment::allelic(0, 5_000, 2.0, 0.0)]);
        assert_eq!(track.loh_spans(), [(0, 5_000)]);
        let svg = drawn(track.clone());
        assert!(svg.contains("heterozygosity lost"));

        let hidden = drawn(track.show_alleles(false));
        assert!(
            svg.matches("<rect ").count() > hidden.matches("<rect ").count(),
            "hiding the lane removed nothing, so it was drawing nothing"
        );
    }

    #[test]
    fn the_pair_is_ordered_and_a_half_read_pair_is_not_a_pair() {
        assert_eq!(
            CopyNumber::allelic(1.0, 4.0),
            CopyNumber::Allelic {
                major: 4.0,
                minor: 1.0
            }
        );
        // One number missing is not an allele split with a zero in it, and it
        // is not a total of the number that survived either: three major copies
        // and an unmeasured minor is a total of three plus something nobody
        // counted, and calling that three puts the something at nought.
        let half = CopyNumber::allelic(3.0, f64::NAN);
        assert!(!half.is_called(), "a half-read pair became a call");
        assert_eq!(half.minor(), None);
        let segment = CopyNumberSegment::allelic(0, 5_000, 3.0, f64::NAN);
        assert_eq!(segment.loh(), None);
        let svg = drawn(CopyNumberTrack::diploid(vec![segment]));
        assert!(svg.contains("no call"), "{svg}");
        assert!(!svg.contains("3 copies"), "an unmeasured total was drawn");
    }

    #[test]
    fn the_ladder_leaves_room_for_a_gain_even_when_nothing_gained() {
        let flat = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, 2.0)]);
        assert!(
            flat.ceiling() >= 4.0,
            "a flat genome has no room to show a gain: {}",
            flat.ceiling()
        );
        let tall = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, 12.0)]);
        assert!(tall.ceiling() >= 12.0);
        assert_eq!(tall.clone().cap(20.0).ceiling(), 20.0);
    }

    #[test]
    fn the_rungs_are_whole_copies_and_thin_out_rather_than_crowd() {
        let track = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, 40.0)]);
        for copies in track.rungs(60.0) {
            assert_eq!(
                copies.fract(),
                0.0,
                "{copies} copies is not a number of copies"
            );
        }
        let roomy = track.rungs(400.0).len();
        let cramped = track.rungs(40.0).len();
        assert!(
            cramped < roomy,
            "the rungs did not thin out: {cramped} against {roomy}"
        );
    }

    #[test]
    fn the_rungs_keep_thinning_past_ten_rather_than_painting_the_band() {
        // A step that stops at ten is a rung every two thirds of a pixel once
        // the ceiling passes about six hundred and seventy, and the whole
        // plotting area is then filled with the rule colour.
        let track = CopyNumberTrack::diploid(vec![CopyNumberSegment::total(0, 5_000, 5_000.0)]);
        let rungs = track.rungs(70.0);
        assert!(
            rungs.len() < 20,
            "{} rungs in a seventy pixel band",
            rungs.len()
        );
        let svg = drawn(track);
        assert!(
            svg.matches("<line ").count() < 20,
            "the ladder painted the band"
        );
    }

    #[test]
    fn the_key_names_the_inks_that_were_drawn_and_no_others() {
        let theme = Theme::dark();
        let track = CopyNumberTrack::diploid(vec![
            CopyNumberSegment::allelic(0, 5_000, 6.0, 1.0),
            CopyNumberSegment::allelic(5_000, 9_000, 2.0, 0.0),
        ]);
        let keys = format!("{:?}", track.legend(&theme).items());
        // The theme's inks, not the light theme's literals.
        assert!(keys.contains(theme.color(1)), "{keys}");
        // Nothing here is below the ploidy, so there is no loss to name.
        assert!(!keys.contains("loss"), "{keys}");
    }

    #[test]
    fn the_balanced_rule_is_drawn_over_an_empty_track() {
        // An empty band showing its rule says nothing was called here. An empty
        // rectangle says nothing at all.
        let svg = drawn(CopyNumberTrack::diploid(Vec::new()).label("nothing"));
        assert!(svg.contains("<line "), "no rule on an empty track");
        assert!(svg.contains("nothing"));
    }

    #[test]
    fn a_ploidy_that_is_not_a_number_does_not_become_the_figure() {
        let track =
            CopyNumberTrack::at_ploidy(vec![CopyNumberSegment::total(0, 100, 2.0)], f64::NAN);
        assert!(track.ploidy().is_finite());
        let svg = drawn(track);
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn a_segment_narrower_than_a_pixel_still_reaches_the_figure() {
        // A focal amplification inside a quiet arm. Averaged into its column it
        // would become a level nobody called, so it is drawn at its own level
        // with a riser down to the arm.
        let segments = vec![
            CopyNumberSegment::total(0, 4_990, 2.0),
            CopyNumberSegment::total(4_990, 4_991, 9.0),
            CopyNumberSegment::total(4_991, 10_000, 2.0),
        ];
        let svg = drawn(CopyNumberTrack::diploid(segments));
        assert!(svg.contains("9 copies"), "the focal event vanished");
    }
}
