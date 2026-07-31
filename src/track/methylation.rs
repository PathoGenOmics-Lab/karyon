//! Per-site methylation, one lane per strand.
//!
//! # The idea
//!
//! A methylation call is not a variant. The base is the same base; what changed
//! is a chemical group on it, and the measurement is a fraction of reads rather
//! than a genotype. Two things follow, and both of them are why this is not a
//! [`VariantTrack`](crate::VariantTrack).
//!
//! The first is strand. Methylation is a property of one strand of a duplex,
//! and the two strands of the same site are two measurements, not one. In
//! bacteria the asymmetry is the finding: a fully methylated GATC has both
//! strands done, a hemimethylated one has just been replicated, and a plot that
//! averages the pair cannot tell you which. So each strand gets a lane.
//!
//! The second is coverage. A site called from four reads and a site called from
//! four hundred are drawn the same size by anything that plots the fraction
//! alone, and the first is noise. Sites under
//! [`MethylationTrack::min_coverage`] are dropped, and the rest fade with how
//! well covered they are, so a thin call looks thin.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Strand, Track};

/// One methylation call at one position on one strand.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MethylSite {
    /// Position of the modified base, 0-based.
    pub pos: u64,
    /// Which strand carries the modification.
    pub strand: Strand,
    /// Fraction of reads calling it modified, between 0 and 1.
    pub fraction: f64,
    /// How many reads the call was made from.
    pub coverage: u32,
}

impl MethylSite {
    /// A call at `pos`.
    pub fn new(pos: u64, strand: Strand, fraction: f64, coverage: u32) -> Self {
        MethylSite {
            pos,
            strand,
            fraction: fraction.clamp(0.0, 1.0),
            coverage,
        }
    }

    /// Whether the call is on the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.strand == Strand::Reverse
    }
}

/// Methylation calls, forward strand above the line and reverse below.
///
/// ```
/// use karyon::{Figure, MethylSite, MethylationTrack, Region, Strand};
///
/// let sites = vec![
///     MethylSite::new(1_010, Strand::Forward, 0.95, 40),
///     MethylSite::new(1_010, Strand::Reverse, 0.12, 38),
/// ];
///
/// let svg = Figure::new(Region::new("NC_000962.3", 1_000, 1_100).unwrap())
///     .push(MethylationTrack::new(sites).label("6mA"))
///     .to_svg();
/// assert!(svg.contains("6mA"));
/// ```
#[derive(Debug, Clone)]
pub struct MethylationTrack {
    sites: Vec<MethylSite>,
    label: Option<String>,
    height: f64,
    forward_color: Option<String>,
    reverse_color: Option<String>,
    min_coverage: u32,
    saturating_coverage: u32,
    radius: f64,
    show_scale: bool,
    show_stems: bool,
}

impl MethylationTrack {
    /// A track over `sites`.
    pub fn new(sites: impl Into<Vec<MethylSite>>) -> Self {
        MethylationTrack {
            sites: sites.into(),
            label: None,
            height: 76.0,
            forward_color: None,
            reverse_color: None,
            min_coverage: 5,
            saturating_coverage: 30,
            radius: 2.6,
            show_scale: true,
            show_stems: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(20.0);
        self
    }

    /// Sets the colours of the two strands.
    pub fn colors(mut self, forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        self.forward_color = Some(forward.into());
        self.reverse_color = Some(reverse.into());
        self
    }

    /// Drops sites called from fewer than this many reads.
    ///
    /// The default of five is a floor, not a recommendation. A methylation
    /// fraction from four reads takes five values and none of them mean much;
    /// what the right number is depends on the depth of the run and on how
    /// small a difference the study is trying to see.
    pub fn min_coverage(mut self, reads: u32) -> Self {
        self.min_coverage = reads;
        self
    }

    /// Sets the coverage at which a site is drawn at full strength.
    pub fn saturating_coverage(mut self, reads: u32) -> Self {
        self.saturating_coverage = reads.max(1);
        self
    }

    /// Sets the radius of one site marker.
    pub fn radius(mut self, radius: f64) -> Self {
        self.radius = radius.max(0.5);
        self
    }

    /// Draws or hides the fraction axis.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Draws or hides the stem joining each marker to the midline.
    pub fn show_stems(mut self, show: bool) -> Self {
        self.show_stems = show;
        self
    }

    /// The calls.
    pub fn sites(&self) -> &[MethylSite] {
        &self.sites
    }

    /// Calls that clear the coverage floor.
    pub fn confident(&self) -> impl Iterator<Item = &MethylSite> {
        self.sites
            .iter()
            .filter(|site| site.coverage >= self.min_coverage)
    }

    /// How many calls the coverage floor dropped.
    pub fn discarded(&self) -> usize {
        self.sites.len() - self.confident().count()
    }

    /// Positions where the two strands disagree by more than `by`.
    ///
    /// Hemimethylation: one strand modified and the other not. In an organism
    /// that methylates symmetrically it means the site was caught between
    /// replication and maintenance, and finding those is usually the reason for
    /// drawing the two strands apart in the first place.
    pub fn hemimethylated(&self, by: f64) -> Vec<u64> {
        let mut found = Vec::new();
        for site in self.confident().filter(|site| !site.is_reverse()) {
            let partner = self
                .confident()
                .find(|other| other.pos == site.pos && other.is_reverse());
            if let Some(partner) = partner {
                if (site.fraction - partner.fraction).abs() > by {
                    found.push(site.pos);
                }
            }
        }
        found.sort_unstable();
        found.dedup();
        found
    }

    /// How solid a marker is drawn, from how many reads are behind it.
    fn weight(&self, coverage: u32) -> f64 {
        let ratio = coverage as f64 / self.saturating_coverage as f64;
        0.35 + 0.65 * ratio.clamp(0.0, 1.0)
    }
}

impl Track for MethylationTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.sites.is_empty() {
            return 0.0;
        }
        text_width("100%", theme.font_size - 1.0) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let middle = band.mid_y();
        let half = band.h / 2.0;

        // The midline is zero for both lanes, and it is drawn whether or not
        // there is data: it is what says the band has two halves.
        ctx.svg
            .line(band.x, middle, band.right(), middle, &ctx.theme.rule, 1.0);

        let forward = self
            .forward_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(0).to_string());
        let reverse = self
            .reverse_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());
        let stem = mix(ctx.theme.surface(), &ctx.theme.rule, 0.7);

        for site in self.confident() {
            if !ctx.region.contains(site.pos) {
                continue;
            }
            let x = ctx.scale.x_center(site.pos);
            let reach = site.fraction * (half - self.radius - 1.0).max(1.0);
            let y = if site.is_reverse() {
                middle + reach
            } else {
                middle - reach
            };
            let hue = if site.is_reverse() {
                &reverse
            } else {
                &forward
            };
            // Faded towards the page by how few reads are behind it, so a call
            // from six reads does not look like a call from six hundred.
            let color = mix(ctx.theme.surface(), hue, self.weight(site.coverage));

            if self.show_stems {
                ctx.svg.line(x, middle, x, y, &stem, 0.8);
            }
            // The ring is what keeps two sites a base apart reading as two.
            ctx.svg.circle_ringed(
                x,
                y,
                self.radius,
                &color,
                ctx.theme.surface(),
                self.radius * 0.45,
            );
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            for (y, text) in [
                (band.y + size, "100%"),
                (middle + size * 0.35, "0"),
                (band.bottom() - size * 0.22, "100%"),
            ] {
                ctx.svg
                    .text(right, y, text, &ctx.theme.muted, size, Anchor::End);
            }
        }

        let hidden = self.discarded();
        if hidden > 0 {
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 2.0,
                &format!("{hidden} under {}x", self.min_coverage),
                &ctx.theme.muted,
                ctx.theme.font_size - 2.0,
                Anchor::End,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn sites() -> Vec<MethylSite> {
        vec![
            MethylSite::new(1_010, Strand::Forward, 0.95, 40),
            MethylSite::new(1_010, Strand::Reverse, 0.92, 38),
            MethylSite::new(1_040, Strand::Forward, 0.90, 44),
            MethylSite::new(1_040, Strand::Reverse, 0.05, 41),
            MethylSite::new(1_070, Strand::Forward, 0.50, 2),
        ]
    }

    fn region() -> Region {
        Region::new("NC_000962.3", 1_000, 1_100).unwrap()
    }

    #[test]
    fn a_fraction_outside_its_range_is_brought_back_in() {
        assert_eq!(MethylSite::new(1, Strand::Forward, 4.0, 10).fraction, 1.0);
        assert_eq!(MethylSite::new(1, Strand::Forward, -1.0, 10).fraction, 0.0);
    }

    #[test]
    fn a_thin_call_is_dropped_and_counted() {
        let track = MethylationTrack::new(sites());
        assert_eq!(track.confident().count(), 4);
        assert_eq!(track.discarded(), 1, "the two read call");

        // And the floor can be moved.
        let strict = MethylationTrack::new(sites()).min_coverage(42);
        assert_eq!(strict.confident().count(), 1);
        assert_eq!(strict.discarded(), 4);
    }

    #[test]
    fn hemimethylation_is_the_two_strands_disagreeing() {
        let track = MethylationTrack::new(sites());
        // 1,040 is 0.90 on one strand and 0.05 on the other; 1,010 is neither.
        assert_eq!(track.hemimethylated(0.5), vec![1_040]);
        assert!(track.hemimethylated(0.99).is_empty());
    }

    #[test]
    fn a_site_with_no_partner_is_not_hemimethylated() {
        // One strand alone is one measurement, not a disagreement.
        let lonely = vec![MethylSite::new(1_010, Strand::Forward, 1.0, 40)];
        assert!(MethylationTrack::new(lonely).hemimethylated(0.1).is_empty());
    }

    #[test]
    fn a_thin_call_cannot_make_a_pair_look_hemimethylated() {
        let noisy = vec![
            MethylSite::new(1_010, Strand::Forward, 1.0, 40),
            MethylSite::new(1_010, Strand::Reverse, 0.0, 1),
        ];
        assert!(MethylationTrack::new(noisy).hemimethylated(0.5).is_empty());
    }

    #[test]
    fn coverage_decides_how_solid_a_marker_is() {
        let track = MethylationTrack::new(sites()).saturating_coverage(30);
        assert!(track.weight(60) > track.weight(10));
        assert_eq!(track.weight(30), 1.0, "saturated");
        assert_eq!(track.weight(60), 1.0, "and no more than saturated");
        assert!(track.weight(0) > 0.0, "still visible");
    }

    #[test]
    fn the_two_strands_go_either_side_of_the_line() {
        let deep: Vec<MethylSite> = sites()
            .iter()
            .map(|site| MethylSite {
                coverage: 60,
                ..*site
            })
            .collect();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(deep).colors("#111111", "#222222"))
            .to_svg();
        assert!(svg.contains("#111111"), "nothing on the forward strand");
        assert!(svg.contains("#222222"), "nothing on the reverse strand");
    }

    #[test]
    fn a_thinly_covered_call_is_drawn_paler_than_a_deep_one() {
        let paint = |coverage: u32| {
            Figure::new(region())
                .show_region_label(false)
                .push(
                    MethylationTrack::new(vec![MethylSite::new(
                        1_010,
                        Strand::Forward,
                        1.0,
                        coverage,
                    )])
                    .colors("#000000", "#000000")
                    .show_scale(false),
                )
                .to_svg()
        };
        // Full strength at the saturating coverage, and washed out well below.
        assert!(paint(60).contains("#000000"));
        assert!(!paint(6).contains("#000000"));
    }

    #[test]
    fn the_dropped_calls_are_reported_rather_than_hidden() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(sites()))
            .to_svg();
        assert!(svg.contains("under 5x"), "a silent filter is a lie");
    }

    #[test]
    fn sites_outside_the_region_are_not_drawn() {
        let far = vec![MethylSite::new(90_000, Strand::Forward, 1.0, 50)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(far).show_scale(false))
            .to_svg();
        assert!(!svg.contains("<circle"));
    }

    #[test]
    fn the_axis_names_both_halves() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(sites()))
            .to_svg();
        // A hundred per cent either way from a zero in the middle.
        assert_eq!(svg.matches(">100%</text>").count(), 2);
        assert!(svg.contains(">0</text>"));
    }

    #[test]
    fn an_empty_track_still_draws_its_midline() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(Vec::new()).label("6mA"))
            .to_svg();
        assert!(svg.contains("<line"), "the line is what makes two lanes");
        assert!(!svg.contains("NaN"));
    }
}
