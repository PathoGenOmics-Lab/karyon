//! Site-wise molecular-selection evidence along a coding sequence.
//!
//! A selection scan has two answers which should not be collapsed into one
//! colour: how strong the evidence is, and in which direction the rates point.
//! [`SelectionTrack`] therefore draws an evidence skyline above a signed
//! `log2(omega)` effect strip.  The two panels share genomic x coordinates, so
//! domains, variants and codons can be stacked without re-labelling sites.
//!
//! Missing estimates are omitted rather than drawn as zero.  Exact rates and
//! evidence remain in SVG tooltips; visual values are capped only to keep one
//! very large estimate from flattening the rest of the scan, and a ratio held
//! at a cap is drawn open beside an end written as a bound, `ω ≥ 8`.  The
//! sites past the threshold are named over their marks while there are few
//! enough to read.

use crate::scale::Scale;
use crate::style::{LinePattern, Symbol};
use crate::svg::{text_exact, text_rounded, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::{DrawContext, Track};

/// How far over its mark an episodic site's capsule is centred, and how far
/// over the mark its top stands: the lift, half its height and its ring.
const CAPSULE_LIFT: f64 = 7.0;
const CAPSULE_RISE: f64 = CAPSULE_LIFT + 2.1 + 0.8;

/// The most sites past the threshold the evidence panel names over their
/// marks. Past it the names are a wall of numbers, and the tooltips name
/// each site still.
const NAMED_SITES: usize = 30;

/// Statistical quantity shown in the evidence half of a [`SelectionTrack`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SelectionEvidence {
    /// Frequentist p-values, drawn as `-log10(p)`.
    #[default]
    PValue,
    /// Posterior probability of selection, drawn from zero to one.
    Posterior,
}

/// One tested coding position.
///
/// Positions are 0-based genomic coordinates.  `dS` and `dN` may be replaced
/// by synonymous and nonsynonymous rate parameters from codon models; the
/// ratio is the quantity encoded visually.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectionSite {
    /// Tested position, 0-based.
    pub pos: u64,
    synonymous: Option<f64>,
    nonsynonymous: Option<f64>,
    p_value: Option<f64>,
    posterior: Option<f64>,
    episodic: Option<EpisodicRates>,
    label: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct EpisodicRates {
    beta_minus: f64,
    beta_plus: f64,
    positive_weight: f64,
}

impl SelectionSite {
    /// Starts a site at 0-based genomic position `pos`.
    pub fn new(pos: u64) -> Self {
        SelectionSite {
            pos,
            synonymous: None,
            nonsynonymous: None,
            p_value: None,
            posterior: None,
            episodic: None,
            label: None,
        }
    }

    /// Adds the synonymous and nonsynonymous rate estimates.
    pub fn rates(mut self, synonymous: f64, nonsynonymous: f64) -> Self {
        self.synonymous = finite_nonnegative(synonymous);
        self.nonsynonymous = finite_nonnegative(nonsynonymous);
        self
    }

    /// Adds a frequentist p-value.
    pub fn p_value(mut self, value: f64) -> Self {
        self.p_value = if value.is_finite() && value >= 0.0 {
            Some(value.min(1.0))
        } else {
            None
        };
        self
    }

    /// Adds a posterior probability.
    pub fn posterior(mut self, value: f64) -> Self {
        self.posterior = if value.is_finite() {
            Some(value.clamp(0.0, 1.0))
        } else {
            None
        };
        self
    }

    /// Adds the two nonsynonymous-rate classes used by an episodic model.
    ///
    /// `positive_weight` is the fitted proportion assigned to `beta_plus`; the
    /// remainder is assigned to `beta_minus`.  These values are reported in
    /// tooltips and as a small two-part capsule on the evidence marker.
    pub fn episodic_rates(mut self, beta_minus: f64, beta_plus: f64, positive_weight: f64) -> Self {
        if beta_minus.is_finite()
            && beta_minus >= 0.0
            && beta_plus.is_finite()
            && beta_plus >= 0.0
            && positive_weight.is_finite()
        {
            self.episodic = Some(EpisodicRates {
                beta_minus,
                beta_plus,
                positive_weight: positive_weight.clamp(0.0, 1.0),
            });
        }
        self
    }

    /// Adds a human-readable site, residue or hypothesis label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Returns `dN/dS`, preserving an infinite estimate when `dS = 0`.
    pub fn omega(&self) -> Option<f64> {
        match (self.synonymous, self.nonsynonymous) {
            (Some(ds), Some(dn)) => rate_ratio(ds, dn),
            _ => None,
        }
    }

    /// The omega of each episodic rate class, answered the way
    /// [`omega`](SelectionSite::omega) answers it.
    ///
    /// `None` is a class that cannot be placed on the omega scale at all: a
    /// site that kept no synonymous rate has no denominator to divide by, and
    /// a synonymous rate of nought under a nonsynonymous rate of nought is a
    /// ratio nobody can form.  The capsule used to stand a missing `dS` in as
    /// 1 and paint the classes as though the site had been measured, which is
    /// a value given for the absence of a value.
    fn class_omegas(&self, rates: EpisodicRates) -> (Option<f64>, Option<f64>) {
        let ratio = |beta| self.synonymous.and_then(|alpha| rate_ratio(alpha, beta));
        (ratio(rates.beta_minus), ratio(rates.beta_plus))
    }

    /// Returns the stored p-value, when supplied.
    pub fn p(&self) -> Option<f64> {
        self.p_value
    }

    /// Returns the stored posterior probability, when supplied.
    pub fn probability(&self) -> Option<f64> {
        self.posterior
    }
}

/// A two-tier selection atlas: statistical evidence above, rate effect below.
///
/// ```
/// use karyon::{Figure, Region, SelectionEvidence, SelectionSite, SelectionTrack};
///
/// let sites = vec![
///     SelectionSite::new(41).rates(0.3, 0.05).p_value(0.4),
///     SelectionSite::new(86).rates(0.2, 1.4).p_value(0.002).label("surface loop"),
/// ];
/// let svg = Figure::new(Region::new("gene", 0, 120).unwrap())
///     .push(SelectionTrack::new(sites)
///         .evidence(SelectionEvidence::PValue)
///         .label("FEL"))
///     .to_svg();
/// assert!(svg.contains("surface loop"));
/// ```
#[derive(Debug, Clone)]
pub struct SelectionTrack {
    sites: Vec<SelectionSite>,
    label: Option<String>,
    evidence: SelectionEvidence,
    height: f64,
    p_threshold: f64,
    posterior_threshold: f64,
    neutral_lower: f64,
    neutral_upper: f64,
    saturation: f64,
    show_evidence: bool,
    show_effect: bool,
}

impl SelectionTrack {
    /// Creates a selection atlas for `sites`.
    pub fn new(sites: impl Into<Vec<SelectionSite>>) -> Self {
        SelectionTrack {
            sites: sites.into(),
            label: None,
            evidence: SelectionEvidence::PValue,
            height: 150.0,
            p_threshold: 0.05,
            posterior_threshold: 0.90,
            neutral_lower: 0.90,
            neutral_upper: 1.10,
            saturation: 8.0,
            show_evidence: true,
            show_effect: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Chooses p-values or posterior probability for the evidence skyline.
    pub fn evidence(mut self, evidence: SelectionEvidence) -> Self {
        self.evidence = evidence;
        self
    }

    /// Sets the frequentist significance threshold.
    pub fn p_threshold(mut self, threshold: f64) -> Self {
        if threshold.is_finite() && threshold > 0.0 && threshold <= 1.0 {
            self.p_threshold = threshold;
        }
        self
    }

    /// Sets the posterior-probability threshold.
    pub fn posterior_threshold(mut self, threshold: f64) -> Self {
        if threshold.is_finite() {
            self.posterior_threshold = threshold.clamp(0.0, 1.0);
        }
        self
    }

    /// Sets the omega interval treated as visually neutral.
    pub fn neutral_band(mut self, lower: f64, upper: f64) -> Self {
        if lower.is_finite() && upper.is_finite() && lower > 0.0 && upper >= lower {
            self.neutral_lower = lower;
            self.neutral_upper = upper;
        }
        self
    }

    /// Sets the positive omega value at which effect height and colour cap.
    pub fn saturation(mut self, omega: f64) -> Self {
        if omega.is_finite() && omega > 1.0 {
            self.saturation = omega;
        }
        self
    }

    /// Sets the total track height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(54.0);
        self
    }

    /// Shows or hides the evidence skyline.
    pub fn show_evidence(mut self, show: bool) -> Self {
        self.show_evidence = show;
        self
    }

    /// Shows or hides the signed omega effect strip.
    pub fn show_effect(mut self, show: bool) -> Self {
        self.show_effect = show;
        self
    }

    /// Returns the sites in input order.
    pub fn sites(&self) -> &[SelectionSite] {
        &self.sites
    }

    /// Number of episodic rate classes that have no omega, and so are counted
    /// rather than drawn.
    ///
    /// A class is undrawable when its site kept no synonymous rate to divide
    /// by, or when both of its rates are nought.  The rates themselves stay in
    /// the tooltip; what is missing is the ratio the capsule paints, and a
    /// capsule painted from an assumed denominator is a measurement nobody
    /// made.  Sites outside the drawn region are not in here: those are a
    /// choice about what to show, this is data that cannot be shown.
    pub fn undrawable_rate_class_count(&self) -> usize {
        self.sites
            .iter()
            .filter_map(|site| site.episodic.map(|rates| site.class_omegas(rates)))
            .map(|(minus, plus)| usize::from(minus.is_none()) + usize::from(plus.is_none()))
            .sum()
    }

    /// Returns sites that cross the selected evidence threshold.
    pub fn selected_sites(&self) -> Vec<&SelectionSite> {
        self.sites
            .iter()
            .filter(|site| self.is_selected(site))
            .collect()
    }

    fn is_selected(&self, site: &SelectionSite) -> bool {
        match self.evidence {
            SelectionEvidence::PValue => site.p_value.is_some_and(|p| p <= self.p_threshold),
            SelectionEvidence::Posterior => site
                .posterior
                .is_some_and(|p| p >= self.posterior_threshold),
        }
    }

    fn evidence_value(&self, site: &SelectionSite) -> Option<f64> {
        match self.evidence {
            SelectionEvidence::PValue => {
                site.p_value.map(evidence_height).filter(|v| v.is_finite())
            }
            SelectionEvidence::Posterior => site.posterior,
        }
    }

    fn evidence_threshold(&self) -> f64 {
        match self.evidence {
            SelectionEvidence::PValue => evidence_height(self.p_threshold),
            SelectionEvidence::Posterior => self.posterior_threshold,
        }
    }
}

impl Track for SelectionTrack {
    fn noun(&self) -> &str {
        "selection tests by site"
    }

    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, _theme: &Theme) -> f64 {
        48.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        if !self.show_evidence && !self.show_effect {
            return;
        }

        let header = ctx.px(21.0).min(ctx.band.h * 0.28);
        let top = ctx.band.y + header;
        let available = (ctx.band.bottom() - top).max(8.0);
        let evidence_h = match (self.show_evidence, self.show_effect) {
            (true, true) => available * 0.58,
            (true, false) => available,
            _ => 0.0,
        };
        let gap = if self.show_evidence && self.show_effect {
            ctx.px(7.0).min(available * 0.12)
        } else {
            0.0
        };
        let effect_top = top + evidence_h + gap;
        let effect_h = (ctx.band.bottom() - effect_top).max(0.0);

        self.draw_header(ctx);
        if self.show_evidence {
            self.draw_evidence(ctx, top, evidence_h.max(4.0));
        }
        if self.show_effect && effect_h > 3.0 {
            self.draw_effect(ctx, effect_top, effect_h);
        }
    }
}

impl SelectionTrack {
    /// A key to the colours and the diamond. The track's own name is in the
    /// gutter, where every track's is; a title here read `selection atlas`
    /// over every figure of one, whatever it held.
    fn draw_header(&self, ctx: &mut DrawContext<'_>) {
        let y = ctx.band.y + ctx.px(12.5);
        let size = ctx.theme.font_size * 0.78;
        let mut x = ctx.band.x + ctx.px(4.0);
        for (label, omega) in [("ω < 1", 0.3), ("ω ≈ 1", 1.0), ("ω > 1", 4.0)] {
            let color = selection_color(
                ctx.theme,
                omega,
                self.neutral_lower,
                self.neutral_upper,
                self.saturation,
            );
            ctx.svg.circle(x, y - ctx.px(3.0), ctx.px(2.7), &color);
            ctx.svg.text(
                x + ctx.px(5.0),
                y,
                label,
                &ctx.theme.muted,
                size,
                Anchor::Start,
            );
            x += ctx.px(5.0) + text_width(label, size) + ctx.px(12.0);
        }

        // The legend states the rule the track selects sites by, so it prints
        // the threshold and not a rounding of it: three decimals turn a
        // genome-wide 5e-8 into `p ≤ 0`, a claim no site can satisfy, and two
        // turn a posterior threshold of 0.995 into `PP ≥ 1`.
        let evidence = match self.evidence {
            SelectionEvidence::PValue => format!("p ≤ {}", text_exact(self.p_threshold)),
            SelectionEvidence::Posterior => {
                format!("PP ≥ {}", text_exact(self.posterior_threshold))
            }
        };
        ctx.svg.symbol(
            x,
            y - ctx.px(3.0),
            ctx.px(3.0),
            Symbol::Diamond,
            ctx.theme.color(1),
        );
        ctx.svg.text(
            x + ctx.px(6.0),
            y,
            &evidence,
            &ctx.theme.muted,
            size,
            Anchor::Start,
        );
        let right = x + ctx.px(6.0) + text_width(&evidence, size);
        ctx.svg
            .keep_clear(ctx.band.x, y - size, right - ctx.band.x, size * 1.3);
    }

    fn draw_evidence(&self, ctx: &mut DrawContext<'_>, top: f64, height: f64) {
        let bottom = top + height - ctx.px(4.0);
        let threshold = self.evidence_threshold();
        let max = match self.evidence {
            SelectionEvidence::PValue => self
                .sites
                .iter()
                .filter(|site| ctx.region.contains(site.pos))
                .filter_map(|site| self.evidence_value(site))
                .fold(threshold.max(1.0), f64::max)
                .min(16.0)
                .max(threshold * 1.12),
            SelectionEvidence::Posterior => 1.0,
        };
        // The sites past the threshold are named over their marks while
        // there are few enough to read, and the panel keeps a line at its top
        // for the names of the tallest. Nine diamonds over a protein with no
        // number on any had a reader working the sites out with awk.
        let selected = self
            .sites
            .iter()
            .filter(|site| ctx.region.contains(site.pos) && self.is_selected(site))
            .filter(|site| self.evidence_value(site).is_some())
            .count();
        let named = (1..=NAMED_SITES).contains(&selected);
        let name_size = ctx.theme.font_size * 0.78;
        // A name goes over the capsule of an episodic site, which stands
        // over its mark, so the room is taller by a capsule where one could
        // be named.
        let capsules = named
            && self
                .sites
                .iter()
                .any(|site| site.episodic.is_some() && self.is_selected(site));
        let room = match (named, capsules) {
            (false, _) => 0.0,
            (true, false) => name_size + ctx.px(3.0),
            (true, true) => name_size + ctx.px(3.0 + CAPSULE_RISE),
        };
        let usable = (height - ctx.px(8.0) - room).max(2.0);
        let map_y = |value: f64| bottom - value.clamp(0.0, max) / max * usable;
        let mut marks: Vec<(f64, f64, u64)> = Vec::new();

        ctx.svg.line(
            ctx.band.x,
            bottom,
            ctx.band.right(),
            bottom,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );
        let threshold_y = map_y(threshold);
        ctx.svg.line_pattern(
            ctx.band.x,
            threshold_y,
            ctx.band.right(),
            threshold_y,
            &mix(&ctx.theme.rule, ctx.theme.color(1), 0.45),
            ctx.theme.tokens.hairline,
            LinePattern::Dashed,
        );
        self.axis_text(
            ctx,
            top + ctx.px(7.0),
            &evidence_axis_top(self.evidence, max),
        );
        self.axis_text(
            ctx,
            threshold_y + ctx.px(3.0),
            &evidence_axis_threshold(self),
        );

        for site in &self.sites {
            if !ctx.region.contains(site.pos) {
                continue;
            }
            let Some(value) = self.evidence_value(site) else {
                continue;
            };
            let x = ctx.scale.x_center(site.pos);
            let y = map_y(value);
            let selected = self.is_selected(site);
            let color = site
                .omega()
                .map(|omega| {
                    selection_color(
                        ctx.theme,
                        omega,
                        self.neutral_lower,
                        self.neutral_upper,
                        self.saturation,
                    )
                })
                .unwrap_or_else(|| ctx.theme.accent.clone());
            ctx.svg.begin_titled(&site_title(site));
            ctx.svg.line(
                x,
                bottom,
                x,
                y,
                &mix(ctx.theme.surface(), &color, 0.78),
                ctx.theme.tokens.hairline,
            );
            let symbol = if selected {
                Symbol::Diamond
            } else {
                Symbol::Circle
            };
            let radius = ctx.px(if selected { 3.7 } else { 2.3 });
            ctx.svg.symbol_ringed(
                x,
                y,
                radius,
                symbol,
                &color,
                ctx.theme.surface(),
                ctx.px(1.15),
            );
            if let Some(episodic) = site.episodic {
                draw_episodic_capsule(ctx, x, y - ctx.px(CAPSULE_LIFT), site, episodic, self);
            }
            ctx.svg.end_group();
            if named && selected {
                let top = if site.episodic.is_some() {
                    y - ctx.px(CAPSULE_RISE)
                } else {
                    y - radius
                };
                marks.push((x, top - ctx.px(2.5), site.pos));
            }
        }
        if named {
            for (x, y, text) in site_names(marks, name_size, ctx.band.x, ctx.band.right()) {
                let width = text_width(&text, name_size);
                ctx.svg
                    .keep_clear(x - width / 2.0, y - name_size, width, name_size * 1.3);
                ctx.svg.text(
                    x,
                    y,
                    &text,
                    &ctx.theme.foreground,
                    name_size,
                    Anchor::Middle,
                );
            }
        }
    }

    fn draw_effect(&self, ctx: &mut DrawContext<'_>, top: f64, height: f64) {
        let mid = top + height / 2.0;
        let half = (height / 2.0 - ctx.px(4.0)).max(1.0);
        let log_cap = self.saturation.log2().max(1.0);
        let cap = log_cap.exp2();
        // A ratio past either end of the strip stands at the end, where it
        // reads as that end: sites of ω 8.80 and 8.91 stood at `ω 8`, and 44
        // under an eighth at `1/8`, as if each were exactly that. Such a mark
        // is drawn open, and the end it is held at is written as a bound.
        let past_high = |omega: f64| omega.is_infinite() || omega > cap;
        let past_low = |omega: f64| omega <= 0.0 || omega < 1.0 / cap;
        let shown = || {
            self.sites
                .iter()
                .filter(|site| ctx.region.contains(site.pos))
                .filter_map(SelectionSite::omega)
        };
        let (any_high, any_low) = (shown().any(past_high), shown().any(past_low));
        let neutral_top = omega_y(self.neutral_upper, mid, half, log_cap);
        let neutral_bottom = omega_y(self.neutral_lower, mid, half, log_cap);
        ctx.svg.rect_opacity(
            ctx.band.x,
            neutral_top.min(neutral_bottom),
            ctx.band.w,
            (neutral_bottom - neutral_top).abs().max(ctx.px(2.0)),
            &ctx.theme.rule,
            0.28,
        );
        ctx.svg.line(
            ctx.band.x,
            mid,
            ctx.band.right(),
            mid,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );
        let bound = |past: bool, at: &'static str| if past { at } else { "" };
        self.axis_text(
            ctx,
            top + ctx.px(6.0),
            &format!(
                "ω {}{}",
                bound(any_high, "≥ "),
                text_rounded(self.saturation, 1)
            ),
        );
        self.axis_text(ctx, mid + ctx.px(3.0), "ω 1");
        self.axis_text(
            ctx,
            top + height,
            &format!(
                "{}1/{}",
                bound(any_low, "≤ "),
                text_rounded(self.saturation, 1)
            ),
        );

        for site in &self.sites {
            if !ctx.region.contains(site.pos) {
                continue;
            }
            let Some(omega) = site.omega() else {
                continue;
            };
            let x = ctx.scale.x_center(site.pos);
            let y = omega_y(omega, mid, half, log_cap);
            let selected = self.is_selected(site);
            let color = selection_color(
                ctx.theme,
                omega,
                self.neutral_lower,
                self.neutral_upper,
                self.saturation,
            );
            ctx.svg.begin_titled(&site_title(site));
            ctx.svg.line(
                x,
                mid,
                x,
                y,
                &color,
                ctx.px(if selected { 1.8 } else { 1.05 }),
            );
            let symbol = if selected {
                Symbol::Diamond
            } else {
                Symbol::Circle
            };
            let radius = ctx.px(if selected { 3.6 } else { 2.15 });
            if past_high(omega) || past_low(omega) {
                // Open: an edge in its colour round the page, as wide as the
                // filled mark it stands for.
                let edge = ctx.px(1.2);
                ctx.svg.symbol_ringed(
                    x,
                    y,
                    radius + ctx.px(1.0) - edge,
                    symbol,
                    ctx.theme.surface(),
                    &color,
                    edge,
                );
            } else {
                ctx.svg.symbol_ringed(
                    x,
                    y,
                    radius,
                    symbol,
                    &color,
                    ctx.theme.surface(),
                    ctx.px(1.0),
                );
            }
            ctx.svg.end_group();
        }
    }

    fn axis_text(&self, ctx: &mut DrawContext<'_>, y: f64, text: &str) {
        if ctx.axis.w <= 0.0 {
            return;
        }
        ctx.svg.text(
            ctx.axis.right() - ctx.px(4.0),
            y,
            text,
            &ctx.theme.muted,
            ctx.theme.font_size * 0.68,
            Anchor::End,
        );
    }
}

/// The names written over the marks of the sites past the threshold, each
/// `(x, baseline, text)`, from marks `(x, top, position)`.
///
/// Sites whose names would touch share one, written as their numbers with
/// runs of three or more as ranges, `58, 59, 63-65`, centred over them and
/// kept inside the band. Positions are written from one, as the ruler writes
/// them.
fn site_names(
    mut marks: Vec<(f64, f64, u64)>,
    size: f64,
    left: f64,
    right: f64,
) -> Vec<(f64, f64, String)> {
    marks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let gap = size * 0.5;
    // Each group: its leftmost and rightmost mark, its highest top, and the
    // positions it names.
    let mut groups: Vec<(f64, f64, f64, Vec<u64>)> = Vec::new();
    let reach = |group: &(f64, f64, f64, Vec<u64>)| -> (f64, f64) {
        let width = text_width(&listed(&group.3), size);
        let centre = ((group.0 + group.1) / 2.0)
            .max(left + width / 2.0)
            .min(right - width / 2.0);
        (centre - width / 2.0, centre + width / 2.0)
    };
    for (x, top, pos) in marks {
        groups.push((x, x, top, vec![pos.saturating_add(1)]));
        // A group that grew may reach the one before it, so merging walks
        // back until two neighbours stand apart.
        while groups.len() > 1 {
            let last = groups.len() - 1;
            if reach(&groups[last - 1]).1 + gap <= reach(&groups[last]).0 {
                break;
            }
            let (from, to, top, positions) = groups.pop().expect("two groups");
            let before = groups.last_mut().expect("two groups");
            before.0 = before.0.min(from);
            before.1 = before.1.max(to);
            before.2 = before.2.min(top);
            before.3.extend(positions);
        }
    }
    groups
        .iter()
        .map(|group| {
            let (from, to) = reach(group);
            ((from + to) / 2.0, group.2, listed(&group.3))
        })
        .collect()
}

/// Positions as a reader lists them: in order, a run of three or more as its
/// two ends.
fn listed(positions: &[u64]) -> String {
    let mut sorted = positions.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut parts: Vec<String> = Vec::new();
    let mut at = 0;
    while at < sorted.len() {
        let mut end = at;
        while end + 1 < sorted.len() && sorted[end + 1] == sorted[end] + 1 {
            end += 1;
        }
        if end - at >= 2 {
            parts.push(format!(
                "{}-{}",
                group_thousands(sorted[at]),
                group_thousands(sorted[end])
            ));
        } else {
            parts.extend(sorted[at..=end].iter().map(|pos| group_thousands(*pos)));
        }
        at = end + 1;
    }
    parts.join(", ")
}

fn finite_nonnegative(value: f64) -> Option<f64> {
    (value.is_finite() && value >= 0.0).then_some(value)
}

/// The three answers a rate ratio has, in the one place that knows them.
///
/// A denominator of nought is not a small denominator: `beta / 0` is either the
/// largest ratio there is or no ratio at all, and which of the two it is depends
/// on the numerator.  The capsule used to answer this in two adjacent lines,
/// once with nought and once with infinity, so a site whose two rate classes
/// held the same number was painted deepest purifying beside strongest
/// positive.  Everything that divides a nonsynonymous rate by a synonymous one
/// now gets the same answer.
fn rate_ratio(alpha: f64, beta: f64) -> Option<f64> {
    if alpha > 0.0 {
        Some(beta / alpha)
    } else if beta > 0.0 {
        Some(f64::INFINITY)
    } else {
        None
    }
}

/// The smallest p-value the evidence axis can tell apart from nought.
const EVIDENCE_FLOOR: f64 = 1e-16;

/// The height a p-value reaches on the evidence axis.
///
/// `-log10 0` is not a height, so every p-value is read against the floor
/// before it becomes one.  The significance rule is read against it too: it
/// used to be taken raw, so a threshold below the floor drew its rule above
/// every mark the axis can reach and put every significant site on the wrong
/// side of its own significance line.
fn evidence_height(p: f64) -> f64 {
    -p.max(EVIDENCE_FLOOR).log10()
}

fn omega_y(omega: f64, mid: f64, half: f64, log_cap: f64) -> f64 {
    let signed = if omega.is_infinite() {
        log_cap
    } else if omega > 0.0 {
        omega.log2().clamp(-log_cap, log_cap)
    } else {
        -log_cap
    };
    mid - signed / log_cap * half
}

fn selection_color(
    theme: &Theme,
    omega: f64,
    neutral_lower: f64,
    neutral_upper: f64,
    saturation: f64,
) -> String {
    if omega >= neutral_lower && omega <= neutral_upper {
        return mix(theme.surface(), &theme.muted, 0.55);
    }
    if omega > neutral_upper {
        let strength = if omega.is_infinite() {
            1.0
        } else {
            (omega / neutral_upper).ln() / (saturation / neutral_upper).ln()
        }
        .clamp(0.18, 1.0);
        return mix(&theme.rule, theme.color(1), strength);
    }
    let floor = 1.0 / saturation;
    let strength = if omega <= 0.0 {
        1.0
    } else {
        (neutral_lower / omega).ln() / (neutral_lower / floor).ln()
    }
    .clamp(0.18, 1.0);
    mix(&theme.rule, theme.color(0), strength)
}

fn draw_episodic_capsule(
    ctx: &mut DrawContext<'_>,
    x: f64,
    y: f64,
    site: &SelectionSite,
    rates: EpisodicRates,
    track: &SelectionTrack,
) {
    let colour = |omega: Option<f64>| {
        omega.map(|omega| {
            selection_color(
                ctx.theme,
                omega,
                track.neutral_lower,
                track.neutral_upper,
                track.saturation,
            )
        })
    };
    let (minus_omega, plus_omega) = site.class_omegas(rates);
    let (minus, plus) = (colour(minus_omega), colour(plus_omega));
    if minus.is_none() && plus.is_none() {
        // Neither class has a denominator, so there is no colour here that
        // would not be invented.  The rates stay in the tooltip and the site
        // keeps its marker; the capsule is what cannot be drawn, and
        // `undrawable_rate_class_count` is where the reader is told how many.
        return;
    }
    let width = ctx.px(11.0);
    let height = ctx.px(4.2);
    let left_width = width * (1.0 - rates.positive_weight);
    let right_width = width - left_width;
    ctx.svg.rect_rounded(
        x - width / 2.0 - ctx.px(0.8),
        y - height / 2.0 - ctx.px(0.8),
        width + ctx.px(1.6),
        height + ctx.px(1.6),
        height / 2.0,
        ctx.theme.surface(),
    );
    // A class with no omega and a class with no width are both left out; the
    // rectangle writer drops a width of nought on its own.
    if let Some(fill) = &minus {
        ctx.svg
            .rect(x - width / 2.0, y - height / 2.0, left_width, height, fill);
    }
    if let Some(fill) = &plus {
        ctx.svg.rect(
            x - width / 2.0 + left_width,
            y - height / 2.0,
            right_width,
            height,
            fill,
        );
    }
}

fn site_title(site: &SelectionSite) -> String {
    let mut parts = vec![format!("position {}", site.pos.saturating_add(1))];
    if let Some(label) = &site.label {
        parts.push(label.clone());
    }
    if let (Some(ds), Some(dn)) = (site.synonymous, site.nonsynonymous) {
        // Exact, because this is the only place the rates survive and a rounded
        // rate contradicts the ratio printed beside it: a dN of 1e-7 written as
        // `0` over a dS of `0` asks a reader to believe that nought divided by
        // nought is infinity.  The ratio itself stays rounded, since both of
        // its inputs are here for anyone who wants it to the last bit.
        parts.push(format!("dS {}", text_exact(ds)));
        parts.push(format!("dN {}", text_exact(dn)));
        let omega = match site.omega() {
            Some(value) if value.is_infinite() => "infinity".into(),
            Some(value) => text_rounded(value, 4),
            None => "not estimable".into(),
        };
        parts.push(format!("omega {omega}"));
    }
    if let Some(p) = site.p_value {
        // Exact for the same reason, and one the skyline makes visible: the
        // mark is placed on a log scale, so p 1e-8 and p 1e-12 stand at
        // different heights while six decimals print both of them as `0`, and
        // a reader who hovers to find out why is told they are one number.
        parts.push(format!("p {}", text_exact(p)));
    }
    if let Some(posterior) = site.posterior {
        parts.push(format!("posterior {}", text_exact(posterior)));
    }
    if let Some(rates) = site.episodic {
        parts.push(format!(
            "beta- {}; beta+ {}; positive class weight {}",
            text_exact(rates.beta_minus),
            text_exact(rates.beta_plus),
            text_exact(rates.positive_weight)
        ));
    }
    parts.join(" | ")
}

fn evidence_axis_top(evidence: SelectionEvidence, max: f64) -> String {
    match evidence {
        SelectionEvidence::PValue => format!("-log10 p {}", text_rounded(max, 1)),
        SelectionEvidence::Posterior => "PP 1".into(),
    }
}

fn evidence_axis_threshold(track: &SelectionTrack) -> String {
    match track.evidence {
        // This one labels a height rather than the rule the legend states, so
        // it reports the p-value the dashed line is drawn at: a threshold below
        // the axis floor is drawn at the floor, and printing the asked for
        // number there would label a height with a p-value no mark on the axis
        // can reach.  Exact, for the reason the legend is.
        SelectionEvidence::PValue => {
            format!("p {}", text_exact(track.p_threshold.max(EVIDENCE_FLOOR)))
        }
        SelectionEvidence::Posterior => {
            format!("PP {}", text_exact(track.posterior_threshold))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Figure, Region};

    #[test]
    fn evidence_and_effect_have_exact_tooltips() {
        let track = SelectionTrack::new(vec![SelectionSite::new(10)
            .rates(0.2, 1.4)
            .p_value(0.004)
            .episodic_rates(0.05, 2.8, 0.22)
            .label("surface loop")]);
        let svg = Figure::new(Region::new("gene", 0, 30).unwrap())
            .push(track)
            .to_svg();
        assert!(svg.contains("surface loop"));
        assert!(svg.contains("omega 7"));
        assert!(svg.contains("positive class weight 0.22"));
        assert!(svg.contains("<polygon"), "selected sites need a shape cue");
        assert!(!svg.contains("NaN"));
        assert!(!svg.contains("inf\""));
    }

    /// The sites past the threshold are named over their marks while there
    /// are thirty or fewer, sites whose names would touch sharing one: nine
    /// diamonds with no number on any sent a reader to awk for 58, 59, 63,
    /// 64, 65, 181, 183, 184 and 185.
    #[test]
    fn the_sites_past_the_threshold_are_named_while_few() {
        let site = |at: u64, p: f64| SelectionSite::new(at - 1).rates(0.2, 1.4).p_value(p);
        let mut sites: Vec<SelectionSite> = (1..=300).map(|at| site(at, 0.6)).collect();
        for at in [58, 59, 63, 64, 65, 181, 183, 184, 185] {
            sites[at as usize - 1] = site(at, 0.01);
        }
        let region = Region::new("site", 0, 300).unwrap();
        let svg = Figure::new(region.clone())
            .push(SelectionTrack::new(sites.clone()))
            .to_svg();
        assert!(svg.contains(">58, 59, 63-65</text>"), "{svg}");
        assert!(svg.contains(">181, 183-185</text>"), "{svg}");
        // Thirty-one are a wall of numbers, and none is written.
        let many: Vec<SelectionSite> = (1..=300)
            .map(|at| site(at, if at % 9 == 0 { 0.01 } else { 0.6 }))
            .collect();
        let svg = Figure::new(region).push(SelectionTrack::new(many)).to_svg();
        assert!(!svg.contains(">9, "), "{svg}");
        assert!(!svg.contains(">9</text>"), "{svg}");

        assert_eq!(listed(&[65, 58, 63, 59, 64]), "58, 59, 63-65");
        assert_eq!(listed(&[1_200, 1_201]), "1,200, 1,201");
        assert_eq!(listed(&[761_155, 761_156, 761_157]), "761,155-761,157");
        // A name is kept inside the band at its ends.
        let names = site_names(vec![(1.0, 50.0, 0)], 10.0, 0.0, 100.0);
        assert!(names[0].0 >= text_width("1", 10.0) / 2.0, "{names:?}");
    }

    /// The name of an episodic site goes over its capsule, which stands over
    /// its mark: written over the mark, `46` sat on the capsule.
    #[test]
    fn a_name_goes_over_an_episodic_capsule() {
        let svg = Figure::new(Region::new("site", 0, 60).unwrap())
            .push(SelectionTrack::new(vec![
                SelectionSite::new(10).rates(1.0, 0.5).p_value(0.5),
                SelectionSite::new(45)
                    .rates(0.2, 1.4)
                    .p_value(0.004)
                    .episodic_rates(0.05, 2.8, 0.22),
            ]))
            .to_svg();
        let at = svg.find(">46</text>").expect("the name");
        let open = svg[..at].rfind("<text").unwrap();
        let baseline: f64 = svg[open..at]
            .split(" y=\"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        // The capsule's backing, the one rounded rect of the page colour in
        // the site's group, and its top.
        let group = &svg[svg.find("<title>position 46").unwrap()..];
        let backing = group.split("<rect ").nth(1).unwrap();
        let top: f64 = backing
            .split(" y=\"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert!(
            baseline < top,
            "the name at {baseline} reaches the capsule at {top}"
        );
    }

    /// A ratio past an end of the effect strip stands at the end, and read
    /// as exactly that: ω 8.80 and 8.91 at `ω 8`. It is drawn open and the
    /// end is written as a bound, only where a mark in view is held there.
    #[test]
    fn a_ratio_held_at_an_end_of_the_strip_is_open_and_the_end_a_bound() {
        let theme = Theme::light();
        let drawn = |sites: Vec<SelectionSite>| {
            Figure::new(Region::new("site", 0, 40).unwrap())
                .push(SelectionTrack::new(sites))
                .to_svg()
        };
        let past = drawn(vec![
            SelectionSite::new(5).rates(0.1, 0.88).p_value(0.5),
            SelectionSite::new(20).rates(1.0, 0.05).p_value(0.5),
            SelectionSite::new(30).rates(1.0, 1.2).p_value(0.5),
        ]);
        assert!(past.contains(">ω ≥ 8</text>"), "{past}");
        assert!(past.contains(">≤ 1/8</text>"), "{past}");
        // Two open marks, each a disc of the page inside an edge.
        let open = format!("fill=\"{}\"", theme.surface());
        let effect = past.split(">ω ≥ 8</text>").nth(1).unwrap();
        let discs = effect
            .split("<g><title>")
            .filter(|mark| mark.matches("<circle").count() == 2)
            .filter(|mark| mark.split("<circle").nth(2).unwrap().contains(&open))
            .count();
        assert_eq!(discs, 2, "{effect}");

        let within = drawn(vec![
            SelectionSite::new(5).rates(0.1, 0.8).p_value(0.5),
            SelectionSite::new(20).rates(1.0, 0.125).p_value(0.5),
        ]);
        assert!(within.contains(">ω 8</text>"), "{within}");
        assert!(within.contains(">1/8</text>"), "{within}");
    }

    #[test]
    fn posterior_mode_uses_its_own_threshold() {
        let track = SelectionTrack::new(vec![
            SelectionSite::new(2).posterior(0.89),
            SelectionSite::new(3).posterior(0.97),
        ])
        .evidence(SelectionEvidence::Posterior)
        .posterior_threshold(0.95);
        assert_eq!(track.selected_sites().len(), 1);
        let svg = Figure::new(Region::new("gene", 0, 10).unwrap())
            .push(track)
            .to_svg();
        assert!(svg.contains("PP ≥ 0.95"));
    }

    #[test]
    fn zero_synonymous_rate_never_emits_nonfinite_geometry() {
        let site = SelectionSite::new(4).rates(0.0, 1.0).p_value(0.01);
        assert!(site.omega().unwrap().is_infinite());
        let svg = Figure::new(Region::new("gene", 0, 10).unwrap())
            .push(SelectionTrack::new(vec![site]))
            .to_svg();
        assert!(svg.contains("omega infinity"));
        assert!(!svg.contains("NaN"));
        assert!(!svg.contains("Infinity"));
    }

    /// Fills of the rectangles inside one site's tooltipped group: the capsule
    /// surround, then one rectangle per rate class the drawing could place.
    fn capsule_fills<'a>(svg: &'a str, site: &str) -> Vec<&'a str> {
        svg.split("<g><title>")
            .find(|group| group.starts_with(site))
            .and_then(|group| group.split("</g>").next())
            .unwrap_or_default()
            .split("<rect ")
            .skip(1)
            .filter_map(|piece| piece.split("fill=\"").nth(1))
            .filter_map(|piece| piece.split('"').next())
            .collect()
    }

    #[test]
    fn two_equal_rate_classes_over_a_zero_synonymous_rate_share_one_colour() {
        // Both halves of the capsule put the same question to the same number,
        // so they cannot come back with opposite answers. A dS of nought used
        // to make the negative class nought and the positive class infinity in
        // adjacent lines, and one site was painted deepest purifying beside
        // strongest positive from a single pair of identical rates.
        let svg = Figure::new(Region::new("gene", 0, 10).unwrap())
            .push(SelectionTrack::new(vec![SelectionSite::new(4)
                .rates(0.0, 1.0)
                .p_value(0.01)
                .episodic_rates(2.0, 2.0, 0.5)]))
            .to_svg();
        let fills = capsule_fills(&svg, "position 5");
        assert_eq!(fills.len(), 3, "a surround and two classes: {fills:?}");
        assert_eq!(
            fills[1], fills[2],
            "two rate classes of 2 are painted {} and {}",
            fills[1], fills[2]
        );
    }

    #[test]
    fn rate_classes_with_no_denominator_are_counted_rather_than_drawn() {
        // A site with no synonymous rate has no ratio, and the capsule used to
        // divide by 1 and paint it as though the site had been measured. The
        // rates stay in the tooltip, the colours are not invented, and the
        // count says how many classes were left out.
        let track = SelectionTrack::new(vec![
            SelectionSite::new(4)
                .p_value(0.01)
                .episodic_rates(0.05, 2.8, 0.3),
            SelectionSite::new(6)
                .rates(0.0, 0.0)
                .p_value(0.01)
                .episodic_rates(0.0, 3.0, 0.3),
        ]);
        assert_eq!(track.undrawable_rate_class_count(), 3);
        let svg = Figure::new(Region::new("gene", 0, 10).unwrap())
            .push(track)
            .to_svg();
        assert!(
            capsule_fills(&svg, "position 5").is_empty(),
            "a capsule drawn from an assumed dS"
        );
        assert!(svg.contains("beta- 0.05"), "the rates are still reported");
        assert_eq!(
            capsule_fills(&svg, "position 7").len(),
            2,
            "the class with a ratio is still drawn"
        );
    }

    #[test]
    fn a_threshold_below_the_axis_floor_keeps_its_sites_on_their_own_side() {
        // The marks and the rule are two readings of one axis, so they are
        // floored alike. Taken raw, a threshold of 1e-300 drew its rule at 300
        // where no mark can reach, and every significant site landed below its
        // own significance line.
        let track = SelectionTrack::new(vec![
            SelectionSite::new(3).p_value(0.0),
            SelectionSite::new(5).p_value(1e-20),
        ])
        .p_threshold(1e-300);
        let rule = track.evidence_threshold();
        assert!(rule.is_finite(), "the rule is drawn at {rule}");
        // A mark on the rule makes no claim either way, which is where a site
        // whose p-value the axis cannot separate from the threshold belongs.
        for site in track.sites() {
            let value = track.evidence_value(site).expect("a p-value has a height");
            if track.is_selected(site) {
                assert!(
                    value >= rule,
                    "a selected site sits at {value}, under a rule at {rule}"
                );
            } else {
                assert!(
                    value <= rule,
                    "an unselected site sits at {value}, over a rule at {rule}"
                );
            }
        }
    }

    #[test]
    fn a_tooltip_keeps_the_numbers_the_drawing_was_given() {
        // Rounding loses the one thing the tooltip is for. Two p-values four
        // orders of magnitude apart both printed `p 0` while the skyline drew
        // them at different heights, and a dN of 1e-7 printed as `0` sat beside
        // the `omega infinity` it had produced.
        let svg = Figure::new(Region::new("gene", 0, 30).unwrap())
            .push(SelectionTrack::new(vec![
                SelectionSite::new(4).rates(0.0, 1e-7).p_value(1e-8),
                SelectionSite::new(9)
                    .rates(2.0, 1.0)
                    .p_value(1e-12)
                    .posterior(1e-9),
            ]))
            .to_svg();
        assert!(svg.contains("dN 1e-7"), "a rate rounded into a nought");
        assert!(svg.contains("omega infinity"));
        assert!(svg.contains("p 1e-8") && svg.contains("p 1e-12"));
        assert!(svg.contains("posterior 1e-9"));
        assert!(!svg.contains("p 0<") && !svg.contains("p 0 "));
    }

    #[test]
    fn a_genome_wide_threshold_is_printed_rather_than_rounded_away() {
        // Three decimals turn 5e-8 into a rule reading `p ≤ 0`, which nothing
        // can satisfy, on a legend that is the only statement of which sites
        // the track picked out.
        let svg = Figure::new(Region::new("gene", 0, 30).unwrap())
            .push(
                SelectionTrack::new(vec![SelectionSite::new(4).rates(0.2, 1.4).p_value(1e-9)])
                    .p_threshold(5e-8),
            )
            .to_svg();
        assert!(svg.contains("p ≤ 5e-8"));
        assert!(svg.contains("p 5e-8"), "the axis labels the rule's height");
        assert!(!svg.contains("p ≤ 0<") && !svg.contains("p ≤ 0.0<"));
    }

    #[test]
    fn missing_and_outside_sites_do_not_become_zeroes() {
        let svg = Figure::new(Region::new("gene", 0, 10).unwrap())
            .push(SelectionTrack::new(vec![
                SelectionSite::new(3),
                SelectionSite::new(30).rates(1.0, 5.0).p_value(1e-8),
            ]))
            .to_svg();
        assert!(!svg.contains("position 31"));
        assert!(!svg.contains("omega 0"));
    }
}
