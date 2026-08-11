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
//! very large estimate from flattening the rest of the scan.

use crate::scale::Scale;
use crate::style::{LinePattern, Symbol};
use crate::svg::{text_rounded, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

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
            (Some(ds), Some(dn)) if ds > 0.0 => Some(dn / ds),
            (Some(0.0), Some(dn)) if dn > 0.0 => Some(f64::INFINITY),
            _ => None,
        }
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
            SelectionEvidence::PValue => site
                .p_value
                .map(|p| -p.max(1e-16).log10())
                .filter(|v| v.is_finite()),
            SelectionEvidence::Posterior => site.posterior,
        }
    }

    fn evidence_threshold(&self) -> f64 {
        match self.evidence {
            SelectionEvidence::PValue => -self.p_threshold.log10(),
            SelectionEvidence::Posterior => self.posterior_threshold,
        }
    }
}

impl Track for SelectionTrack {
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
    fn draw_header(&self, ctx: &mut DrawContext<'_>) {
        let y = ctx.band.y + ctx.px(12.5);
        let size = ctx.theme.font_size * 0.78;
        let mut x = ctx.band.x + ctx.px(4.0);
        ctx.svg.text_bold(
            x,
            y,
            "selection atlas",
            &ctx.theme.foreground,
            size,
            Anchor::Start,
        );
        x += ctx.px(78.0);
        for (label, omega) in [("omega<1", 0.3), ("~1", 1.0), (">1", 4.0)] {
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
            x += ctx.px(if label == "omega<1" { 49.0 } else { 30.0 });
        }

        let evidence = match self.evidence {
            SelectionEvidence::PValue => format!("p<={}", text_rounded(self.p_threshold, 3)),
            SelectionEvidence::Posterior => {
                format!("PP>={}", text_rounded(self.posterior_threshold, 2))
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
        let usable = (height - ctx.px(8.0)).max(2.0);
        let map_y = |value: f64| bottom - value.clamp(0.0, max) / max * usable;

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
                draw_episodic_capsule(ctx, x, y - ctx.px(7.0), episodic, site.synonymous, self);
            }
            ctx.svg.end_group();
        }
    }

    fn draw_effect(&self, ctx: &mut DrawContext<'_>, top: f64, height: f64) {
        let mid = top + height / 2.0;
        let half = (height / 2.0 - ctx.px(4.0)).max(1.0);
        let log_cap = self.saturation.log2().max(1.0);
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
        self.axis_text(
            ctx,
            top + ctx.px(6.0),
            &format!("omega {}", text_rounded(self.saturation, 1)),
        );
        self.axis_text(ctx, mid + ctx.px(3.0), "omega 1");
        self.axis_text(
            ctx,
            top + height,
            &format!("1/{}", text_rounded(self.saturation, 1)),
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
            ctx.svg.symbol_ringed(
                x,
                y,
                ctx.px(if selected { 3.6 } else { 2.15 }),
                if selected {
                    Symbol::Diamond
                } else {
                    Symbol::Circle
                },
                &color,
                ctx.theme.surface(),
                ctx.px(1.0),
            );
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

fn finite_nonnegative(value: f64) -> Option<f64> {
    (value.is_finite() && value >= 0.0).then_some(value)
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
    rates: EpisodicRates,
    synonymous: Option<f64>,
    track: &SelectionTrack,
) {
    let width = ctx.px(11.0);
    let height = ctx.px(4.2);
    let left_width = width * (1.0 - rates.positive_weight);
    let right_width = width - left_width;
    let alpha = synonymous.unwrap_or(1.0);
    let minus_omega = if alpha > 0.0 {
        rates.beta_minus / alpha
    } else {
        0.0
    };
    let plus_omega = if alpha > 0.0 && rates.beta_plus.is_finite() {
        rates.beta_plus / alpha
    } else {
        f64::INFINITY
    };
    let minus = selection_color(
        ctx.theme,
        minus_omega,
        track.neutral_lower,
        track.neutral_upper,
        track.saturation,
    );
    let plus = selection_color(
        ctx.theme,
        plus_omega,
        track.neutral_lower,
        track.neutral_upper,
        track.saturation,
    );
    ctx.svg.rect_rounded(
        x - width / 2.0 - ctx.px(0.8),
        y - height / 2.0 - ctx.px(0.8),
        width + ctx.px(1.6),
        height + ctx.px(1.6),
        height / 2.0,
        ctx.theme.surface(),
    );
    if left_width > 0.0 {
        ctx.svg.rect(
            x - width / 2.0,
            y - height / 2.0,
            left_width,
            height,
            &minus,
        );
    }
    if right_width > 0.0 {
        ctx.svg.rect(
            x - width / 2.0 + left_width,
            y - height / 2.0,
            right_width,
            height,
            &plus,
        );
    }
}

fn site_title(site: &SelectionSite) -> String {
    let mut parts = vec![format!("position {}", site.pos.saturating_add(1))];
    if let Some(label) = &site.label {
        parts.push(label.clone());
    }
    if let (Some(ds), Some(dn)) = (site.synonymous, site.nonsynonymous) {
        parts.push(format!("dS {}", text_rounded(ds, 5)));
        parts.push(format!("dN {}", text_rounded(dn, 5)));
        let omega = match site.omega() {
            Some(value) if value.is_infinite() => "infinity".into(),
            Some(value) => text_rounded(value, 4),
            None => "not estimable".into(),
        };
        parts.push(format!("omega {omega}"));
    }
    if let Some(p) = site.p_value {
        parts.push(format!("p {}", text_rounded(p, 6)));
    }
    if let Some(posterior) = site.posterior {
        parts.push(format!("posterior {}", text_rounded(posterior, 4)));
    }
    if let Some(rates) = site.episodic {
        parts.push(format!(
            "beta- {}; beta+ {}; positive class weight {}",
            text_rounded(rates.beta_minus, 4),
            text_rounded(rates.beta_plus, 4),
            text_rounded(rates.positive_weight, 4)
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
        SelectionEvidence::PValue => format!("p {}", text_rounded(track.p_threshold, 3)),
        SelectionEvidence::Posterior => {
            format!("PP {}", text_rounded(track.posterior_threshold, 2))
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
        assert!(svg.contains("PP&gt;=0.95") || svg.contains("PP>=0.95"));
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
