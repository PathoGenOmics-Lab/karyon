//! How many gene families sit in how many genomes.
//!
//! # The idea
//!
//! Count, for every gene family in a pangenome, how many of the genomes carry
//! it, then count how many families got each answer. The result is the same
//! shape in every bacterial species: a tall bar at "all of them", a tall bar at
//! "one of them", and a long shallow trough between. That U is the pangenome.
//! The left tower is what the species is still picking up and losing, the right
//! tower is what it cannot do without, and the trough is everything whose
//! distribution is a question rather than an answer.
//!
//! It is the companion to
//! [`AccumulationTrack`](crate::AccumulationTrack), and the two say different
//! things. The accumulation curve says whether the pangenome closes; this says
//! what is in it. A collection can have a flat accumulation curve and still be
//! mostly accessory, and the curve alone would not tell you.
//!
//! # Where the lines are drawn
//!
//! Core, shell and cloud are conventions with numbers attached, and the numbers
//! are arguments rather than facts. The defaults here are the ones Roary uses,
//! core at 99% of genomes and cloud under 15%, and
//! [`FrequencyTrack::thresholds`] moves them. Whichever you use belongs in the
//! caption, because a family at 96% is core under one convention and shell
//! under another.
//!
//! # What the x axis is
//!
//! The number of genomes a family is in, not a position. Put the figure over
//! `Region::new("genomes", 0, n)`.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, wash, Theme};
use crate::track::{DrawContext, Track};

/// Which part of the pangenome a family belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Frequency {
    /// In nearly every genome. What the species cannot do without.
    Core,
    /// In a middling number of them. The part whose distribution is a question.
    Shell,
    /// In a handful. What the species is still picking up and losing.
    Cloud,
}

impl Frequency {
    /// The name this class goes by in a caption.
    pub fn name(self) -> &'static str {
        match self {
            Frequency::Core => "core",
            Frequency::Shell => "shell",
            Frequency::Cloud => "cloud",
        }
    }
}

/// The gene frequency spectrum of a pangenome.
///
/// ```
/// use karyon::{Figure, FrequencyTrack, Region};
///
/// // Four genomes: two families in all of them, one in a single genome.
/// let genomes = vec![
///     vec![true, true, false],
///     vec![true, true, false],
///     vec![true, true, false],
///     vec![true, true, true],
/// ];
///
/// let track = FrequencyTrack::from_presence(&genomes);
/// assert_eq!(track.families_in(4), 2);
/// assert_eq!(track.families_in(1), 1);
///
/// let svg = Figure::new(Region::new("genomes", 0, 4).unwrap())
///     .push(track.label("gene families"))
///     .to_svg();
/// assert!(svg.contains("<rect"));
/// ```
#[derive(Debug, Clone)]
pub struct FrequencyTrack {
    counts: Vec<usize>,
    label: Option<String>,
    height: f64,
    core: f64,
    cloud: f64,
    colors: Option<(String, String, String)>,
    show_scale: bool,
    max: Option<usize>,
    log_scale: bool,
}

impl FrequencyTrack {
    /// A spectrum from counts you have already worked out.
    ///
    /// `counts[k - 1]` is how many families sit in exactly `k` genomes.
    pub fn new(counts: impl Into<Vec<usize>>) -> Self {
        FrequencyTrack {
            counts: counts.into(),
            label: None,
            height: 150.0,
            core: 0.99,
            cloud: 0.15,
            colors: None,
            show_scale: true,
            max: None,
            log_scale: true,
        }
    }

    /// The spectrum of a presence and absence matrix.
    ///
    /// One row per genome, one column per gene family, `true` where the genome
    /// carries it. A family in no genome at all is not counted, since it is not
    /// in the pangenome.
    pub fn from_presence(genomes: &[Vec<bool>]) -> Self {
        let families = genomes.iter().map(|row| row.len()).max().unwrap_or(0);
        let mut counts = vec![0usize; genomes.len()];
        for family in 0..families {
            let carried = genomes
                .iter()
                .filter(|row| row.get(family).copied().unwrap_or(false))
                .count();
            if carried > 0 {
                counts[carried - 1] += 1;
            }
        }
        FrequencyTrack::new(counts)
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(30.0);
        self
    }

    /// Moves the lines between core, shell and cloud.
    ///
    /// Both are fractions of the collection. The defaults are Roary's, core at
    /// 0.99 and cloud under 0.15. They are conventions rather than facts, so
    /// whichever you use belongs in the caption.
    pub fn thresholds(mut self, core: f64, cloud: f64) -> Self {
        self.core = core.clamp(0.0, 1.0);
        self.cloud = cloud.clamp(0.0, 1.0);
        self
    }

    /// Sets the colours of the three classes.
    pub fn colors(
        mut self,
        core: impl Into<String>,
        shell: impl Into<String>,
        cloud: impl Into<String>,
    ) -> Self {
        self.colors = Some((core.into(), shell.into(), cloud.into()));
        self
    }

    /// Draws or hides the count axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Plots the counts on a logarithmic axis.
    ///
    /// On by default, and it has to be. The core is thousands of families and
    /// a shell bin is tens, so on a linear axis the two towers of the U are the
    /// only things with any height and the trough between them, which is most
    /// of the accessory genome, is a flat line along the bottom. The axis is
    /// ticked at each power of ten so that nobody has to guess it is a log.
    pub fn log_scale(mut self, log: bool) -> Self {
        self.log_scale = log;
        self
    }

    /// Pins the top of the count axis.
    pub fn max(mut self, families: usize) -> Self {
        self.max = Some(families);
        self
    }

    /// How many genomes the spectrum covers.
    pub fn genomes(&self) -> usize {
        self.counts.len()
    }

    /// Whether the spectrum is empty.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// The counts, `counts()[k - 1]` being the families in exactly `k` genomes.
    pub fn counts(&self) -> &[usize] {
        &self.counts
    }

    /// How many families sit in exactly `genomes` of them.
    pub fn families_in(&self, genomes: usize) -> usize {
        genomes
            .checked_sub(1)
            .and_then(|index| self.counts.get(index))
            .copied()
            .unwrap_or(0)
    }

    /// Which class a family in `genomes` of them belongs to.
    pub fn class(&self, genomes: usize) -> Frequency {
        if self.counts.is_empty() {
            return Frequency::Shell;
        }
        let fraction = genomes as f64 / self.counts.len() as f64;
        if fraction >= self.core {
            Frequency::Core
        } else if fraction < self.cloud {
            Frequency::Cloud
        } else {
            Frequency::Shell
        }
    }

    /// How many families fall in one class.
    ///
    /// These are the numbers a pangenome paper quotes, so they come off the
    /// same thresholds the figure was drawn with rather than being counted
    /// again somewhere else with a different convention.
    pub fn total(&self, class: Frequency) -> usize {
        (1..=self.counts.len())
            .filter(|genomes| self.class(*genomes) == class)
            .map(|genomes| self.families_in(genomes))
            .sum()
    }

    /// Every family, of any class.
    pub fn pangenome(&self) -> usize {
        self.counts.iter().sum()
    }

    /// The colour of one class.
    pub fn color(&self, class: Frequency, theme: &Theme) -> String {
        match &self.colors {
            Some((core, shell, cloud)) => match class {
                Frequency::Core => core.clone(),
                Frequency::Shell => shell.clone(),
                Frequency::Cloud => cloud.clone(),
            },
            None => match class {
                Frequency::Core => theme.color(0).to_string(),
                Frequency::Shell => theme.color(2).to_string(),
                Frequency::Cloud => theme.color(1).to_string(),
            },
        }
    }

    /// How far up the band a count sits, from 0 at the floor to 1 at the top.
    pub fn fraction(&self, families: usize) -> f64 {
        let ceiling = self.ceiling() as f64;
        if families == 0 {
            return 0.0;
        }
        if self.log_scale {
            // log(1 + x), not log(x). A bin holding exactly one family is the
            // commonest bin there is, and log(1) is zero, so the plain
            // logarithm draws the whole left tower as nothing at all.
            (1.0 + families as f64).ln() / (1.0 + ceiling).ln()
        } else {
            families as f64 / ceiling
        }
        .clamp(0.0, 1.0)
    }

    /// The powers of ten the axis is ticked at.
    pub fn ticks(&self) -> Vec<usize> {
        if !self.log_scale {
            return vec![self.ceiling()];
        }
        let mut ticks = Vec::new();
        let mut at = 1usize;
        while at <= self.ceiling() {
            ticks.push(at);
            at = match at.checked_mul(10) {
                Some(next) => next,
                None => break,
            };
        }
        ticks
    }

    /// Top of the count axis.
    pub fn ceiling(&self) -> usize {
        self.max
            .unwrap_or_else(|| self.counts.iter().copied().max().unwrap_or(1))
            .max(1)
    }
}

impl Track for FrequencyTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.counts.is_empty() {
            return 0.0;
        }
        self.ticks()
            .into_iter()
            .map(|tick| text_width(&group_thousands(tick), theme.font_size - 1.0))
            .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let baseline = band.bottom();
        ctx.svg.line(
            band.x,
            baseline - 0.5,
            band.right(),
            baseline - 0.5,
            &ctx.theme.rule,
            1.0,
        );
        if self.counts.is_empty() {
            return;
        }

        let width = (band.w / self.counts.len() as f64).max(1.0);
        for genomes in 1..=self.counts.len() {
            let families = self.families_in(genomes);
            if families == 0 {
                continue;
            }
            let class = self.class(genomes);
            let color = self.color(class, ctx.theme);
            let full = self.fraction(families) * band.h;
            let x = band.x + (genomes - 1) as f64 * width;
            // A bar wide enough to have an edge gets one, in its own hue: a row
            // of abutting bars of one colour is otherwise one block, and how
            // many genomes each step covers is the x axis of the figure.
            let inset = if width > 4.0 { 0.6 } else { 0.0 };
            ctx.svg.rect(
                x + inset,
                baseline - full,
                (width - inset * 2.0).max(0.6),
                full,
                &wash(&color, ctx.theme),
            );
            if inset > 0.0 {
                ctx.svg.rect_outline(
                    x + inset,
                    baseline - full,
                    (width - inset * 2.0).max(0.6),
                    full,
                    &color,
                    1.0,
                );
            }
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            // One label per power of ten, and a faint rule at each, so that
            // nobody has to guess the axis is a log.
            for tick in self.ticks() {
                let y = baseline - self.fraction(tick) * band.h;
                ctx.svg.line(
                    band.x,
                    y,
                    band.right(),
                    y,
                    &mix(ctx.theme.surface(), &ctx.theme.rule, 0.55),
                    0.6,
                );
                ctx.svg.text(
                    right,
                    y + size * 0.35,
                    &group_thousands(tick),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }
    }
}

/// `12345` as `12,345`.
fn group_thousands(value: usize) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, c) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    /// Ten genomes: a core of three, a shell family in half of them, and one
    /// family per genome.
    fn genomes() -> Vec<Vec<bool>> {
        (0..10)
            .map(|genome| {
                (0..14)
                    .map(|family| match family {
                        0..=2 => true,
                        3 => genome < 5,
                        _ => family - 4 == genome,
                    })
                    .collect()
            })
            .collect()
    }

    fn region() -> Region {
        Region::new("genomes", 0, 10).unwrap()
    }

    #[test]
    fn the_spectrum_counts_families_by_how_many_carry_them() {
        let track = FrequencyTrack::from_presence(&genomes());
        assert_eq!(track.genomes(), 10);
        assert_eq!(track.families_in(10), 3, "the core");
        assert_eq!(track.families_in(5), 1, "the shell family");
        assert_eq!(track.families_in(1), 10, "one of its own per genome");
        assert_eq!(track.families_in(7), 0);
        assert_eq!(track.pangenome(), 14);
    }

    #[test]
    fn a_family_in_no_genome_is_not_in_the_pangenome() {
        let empty_column = vec![vec![true, false], vec![true, false]];
        let track = FrequencyTrack::from_presence(&empty_column);
        assert_eq!(track.pangenome(), 1);
    }

    #[test]
    fn the_three_classes_split_on_their_thresholds() {
        let track = FrequencyTrack::from_presence(&genomes());
        // Roary's numbers: core at 99% of ten genomes is all ten.
        assert_eq!(track.class(10), Frequency::Core);
        assert_eq!(track.class(9), Frequency::Shell);
        assert_eq!(track.class(5), Frequency::Shell);
        assert_eq!(track.class(1), Frequency::Cloud);

        assert_eq!(track.total(Frequency::Core), 3);
        assert_eq!(track.total(Frequency::Shell), 1);
        assert_eq!(track.total(Frequency::Cloud), 10);
        // Every family lands in exactly one class.
        assert_eq!(
            track.total(Frequency::Core)
                + track.total(Frequency::Shell)
                + track.total(Frequency::Cloud),
            track.pangenome()
        );
    }

    #[test]
    fn the_thresholds_are_arguments_rather_than_facts() {
        let relaxed = FrequencyTrack::from_presence(&genomes()).thresholds(0.9, 0.15);
        // At ninety per cent, nine of ten genomes is core rather than shell.
        assert_eq!(relaxed.class(9), Frequency::Core);
        let strict = FrequencyTrack::from_presence(&genomes()).thresholds(0.99, 0.55);
        assert_eq!(strict.class(5), Frequency::Cloud);
    }

    #[test]
    fn an_empty_pangenome_has_no_classes_and_does_not_panic() {
        let track = FrequencyTrack::from_presence(&[]);
        assert!(track.is_empty());
        assert_eq!(track.pangenome(), 0);
        assert_eq!(track.ceiling(), 1, "and no zero to divide by");
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track.label("families"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_u_is_drawn_in_three_colours() {
        let theme = Theme::light();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(FrequencyTrack::from_presence(&genomes()))
            .to_svg();
        let track = FrequencyTrack::from_presence(&genomes());
        for class in [Frequency::Core, Frequency::Shell, Frequency::Cloud] {
            let color = track.color(class, &theme);
            assert!(svg.contains(&color), "no {} bars", class.name());
        }
    }

    #[test]
    fn a_single_family_still_has_a_bar_on_the_log_axis() {
        // The commonest bin in any pangenome is "in exactly one genome", and
        // log(1) is zero: the plain logarithm draws the left tower as nothing.
        let track = FrequencyTrack::new(vec![1usize, 0, 0, 900]);
        assert!(track.fraction(1) > 0.0);
        assert!(track.fraction(900) > track.fraction(1));
        assert_eq!(track.fraction(0), 0.0, "and nothing is still nothing");
    }

    #[test]
    fn the_axis_is_ticked_at_every_power_of_ten() {
        let track = FrequencyTrack::new(vec![0usize, 0, 4_512]);
        assert_eq!(track.ticks(), vec![1, 10, 100, 1_000]);
        // Linear has one number on it and no decades to walk.
        assert_eq!(track.clone().log_scale(false).ticks(), vec![4_512]);
    }

    #[test]
    fn a_count_of_zero_draws_no_bar_at_all() {
        // The trough of the U is real: nothing sits in seven of the ten.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(FrequencyTrack::from_presence(&genomes()).show_scale(false))
            .to_svg();
        // Three bars with counts, each washed and edged, plus the page and the
        // clip path's own rectangle. Nothing sits in seven of the ten, and
        // nothing is drawn there.
        assert_eq!(svg.matches("<rect").count(), 3 * 2 + 2);
    }

    #[test]
    fn the_axis_is_only_reserved_when_it_is_wanted() {
        let theme = Theme::light();
        let track = FrequencyTrack::from_presence(&genomes());
        assert!(track.y_axis_width(&theme) > 0.0);
        assert_eq!(track.clone().show_scale(false).y_axis_width(&theme), 0.0);
        assert_eq!(FrequencyTrack::new(Vec::new()).y_axis_width(&theme), 0.0);
    }

    #[test]
    fn thousands_are_grouped() {
        assert_eq!(group_thousands(0), "0");
        assert_eq!(group_thousands(999), "999");
        assert_eq!(group_thousands(4_512), "4,512");
    }
}
