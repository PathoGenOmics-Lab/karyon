//! Per-site methylation, one lane per strand.
//!
//! A methylation call is not a variant. The base is the same base; what changed
//! is a chemical group on it, and the measurement is a fraction of reads rather
//! than a genotype. Two things follow from that, and between them they are why
//! this is not a [`VariantTrack`](crate::VariantTrack).
//!
//! # The two strands are two measurements
//!
//! Methylation is a property of one strand of a duplex, so the two strands of
//! one site are two numbers and not one. The asymmetry is often the finding: a
//! fully methylated palindromic site has both strands done, a hemimethylated one
//! has just been replicated and the second strand has not caught up, and a plot
//! that averages the pair cannot tell you which. Each strand gets a lane, above
//! and below a midline that is zero for both.
//!
//! [`MethylationTrack::hemimethylated`] is that same comparison in numbers,
//! pairing each forward call with the nearest reverse one within
//! [`pair_within`](MethylationTrack::pair_within), which is one base by default
//! and worth checking against how your caller reports the two strands.
//!
//! # What a fraction leaves out
//!
//! How many reads it came from. A site called from four reads and a site called
//! from four hundred are drawn the same size by anything that plots the fraction
//! alone, and the first is noise. Calls under
//! [`min_coverage`](MethylationTrack::min_coverage) are dropped and how many
//! were dropped is printed on the figure, since a filter nobody can see is
//! worse than no filter; the rest fade towards the page by how thin they are, so
//! a call from six reads does not look like a call from six hundred.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::feature::strand_color;
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
    pair_within: u64,
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
            pair_within: 1,
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

    /// How far apart the two strands' calls for one site may sit.
    ///
    /// The two modified bases of a palindrome are not at the same coordinate. A
    /// `GATC` carries its 6mA at the second base reading one strand and at the
    /// third reading the other, so the two calls are one base apart; the two
    /// cytosines of a `CpG` are one apart for the same reason. One base by
    /// default, which covers both. Set it to zero if you register both strands
    /// at the site's anchor rather than at the base that was modified, and
    /// wider for a longer recognition sequence.
    pub fn pair_within(mut self, bases: u64) -> Self {
        self.pair_within = bases;
        self
    }

    /// Positions where the two strands disagree by more than `by`.
    ///
    /// Hemimethylation: one strand modified and the other not. In an organism
    /// that methylates symmetrically it means the site was caught between
    /// replication and maintenance, and finding those is usually the reason for
    /// drawing the two strands apart in the first place.
    ///
    /// The position returned is the forward strand's. Its partner is the
    /// nearest reverse call within [`MethylationTrack::pair_within`], so a
    /// palindrome whose two modified bases are a base apart still pairs.
    pub fn hemimethylated(&self, by: f64) -> Vec<u64> {
        let mut found = Vec::new();
        for site in self.confident().filter(|site| !site.is_reverse()) {
            let partner = self
                .confident()
                .filter(|other| other.is_reverse())
                .filter(|other| other.pos.abs_diff(site.pos) <= self.pair_within)
                .min_by_key(|other| other.pos.abs_diff(site.pos));
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
        text_width("fwd 100%", theme.font_size - 1.0) + 8.0
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
            .unwrap_or_else(|| strand_color(Strand::Forward, ctx.theme).to_string());
        let reverse = self
            .reverse_color
            .clone()
            .unwrap_or_else(|| strand_color(Strand::Reverse, ctx.theme).to_string());
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

            // A marker is a fraction and the fraction is the height, so the one
            // thing the geometry cannot show is how many reads it came from.
            ctx.svg.begin_titled(&site_tooltip(site));
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
            ctx.svg.end_group();
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            // The numbers alone would leave the reader guessing which lane is
            // which strand, and the figure has nothing else that says so.
            // Written where a full fraction actually lands, not at the band
            // edges: the markers stop a radius short of those, so labelling
            // the edges overstates the lane by a tenth.
            let full = (half - self.radius - 1.0).max(1.0);
            for (y, text) in [
                (middle - full + size * 0.35, "fwd 100%"),
                (middle + size * 0.35, "0"),
                (middle + full + size * 0.35, "rev 100%"),
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

/// What a reader hovering one call is told.
///
/// The noun leads, as it does in every other tooltip in the crate: a mark that
/// opened on a bare `3,924,562` was the one glyph in a stack of bands whose
/// answer had to be decoded from its position on the page.
///
/// The percentage is written as a percentage rather than as a bare fraction,
/// and it is followed by the depth it was computed from, because those are the
/// two numbers the marker's height cannot separate: a nine tenths from six
/// reads and a nine tenths from six hundred sit at exactly the same place, and
/// only the wash of the colour tells them apart.
fn site_tooltip(site: &MethylSite) -> String {
    format!(
        "methylation site, {}, {}, {}% modified in {} {}",
        group_thousands(site.pos.saturating_add(1)),
        // The track has two lanes and everything that is not reverse is drawn
        // in the forward one, so this names the lane the marker is in.
        if site.is_reverse() {
            "reverse"
        } else {
            "forward"
        },
        (site.fraction * 100.0).round() as i64,
        // A deep site is a five figure read count, and a count in a tooltip is
        // grouped wherever it appears.
        group_thousands(site.coverage as u64),
        if site.coverage == 1 { "read" } else { "reads" }
    )
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
    fn a_palindrome_pairs_across_the_one_base_the_two_strands_differ_by() {
        // The two 6mA of a GATC are one base apart, and so are the two 5mC of a
        // CpG. Requiring the same coordinate made the caller lie about one of
        // them, which is what the E. coli oriC example ran into.
        let gatc = vec![
            MethylSite::new(1_001, Strand::Forward, 0.95, 40),
            MethylSite::new(1_002, Strand::Reverse, 0.04, 40),
        ];
        assert_eq!(
            MethylationTrack::new(gatc.clone()).hemimethylated(0.5),
            vec![1_001],
            "the forward position is the one reported"
        );
        assert!(
            MethylationTrack::new(gatc)
                .pair_within(0)
                .hemimethylated(0.5)
                .is_empty(),
            "a caller who anchors both strands together can still say so"
        );
    }

    #[test]
    fn a_partner_further_off_than_the_window_is_not_a_partner() {
        let apart = vec![
            MethylSite::new(1_001, Strand::Forward, 0.95, 40),
            MethylSite::new(1_010, Strand::Reverse, 0.04, 40),
        ];
        assert!(MethylationTrack::new(apart).hemimethylated(0.5).is_empty());
    }

    #[test]
    fn the_nearest_reverse_call_is_the_partner() {
        let crowded = vec![
            MethylSite::new(1_001, Strand::Forward, 0.95, 40),
            MethylSite::new(1_002, Strand::Reverse, 0.90, 40),
            MethylSite::new(1_000, Strand::Reverse, 0.02, 40),
        ];
        // Both reverse calls are one base away, so the first one found wins and
        // the pair agrees; what matters is that it does not silently take the
        // one that would make a better story.
        let track = MethylationTrack::new(crowded).pair_within(1);
        assert!(track.hemimethylated(0.5).is_empty());
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
        // A hundred per cent either way from a zero in the middle, and each
        // end says which strand it belongs to.
        assert!(svg.contains(">fwd 100%</text>"));
        assert!(svg.contains(">rev 100%</text>"));
        assert!(svg.contains(">0</text>"));
    }

    #[test]
    fn a_site_is_named_with_its_position_its_strand_and_its_fraction() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(sites()))
            .to_svg();
        // The fraction says what it is a fraction of, because the height of
        // the marker cannot tell six reads from six hundred.
        assert!(
            svg.contains(
                "<title>methylation site, 1,011, forward, 95% modified in 40 reads</title>"
            ),
            "{svg}"
        );
        assert!(
            svg.contains(
                "<title>methylation site, 1,041, reverse, 5% modified in 41 reads</title>"
            ),
            "{svg}"
        );
        // The call from two reads was dropped, so it is not named either.
        assert!(!svg.contains("1,071, forward"), "{svg}");
    }

    #[test]
    fn one_read_is_a_read_and_not_one_reads() {
        let lone = vec![MethylSite::new(1_010, Strand::Forward, 1.0, 1)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(lone).min_coverage(1))
            .to_svg();
        assert!(
            svg.contains(
                "<title>methylation site, 1,011, forward, 100% modified in 1 read</title>"
            ),
            "{svg}"
        );
    }

    #[test]
    fn every_group_a_methylation_track_opens_is_closed_again() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MethylationTrack::new(sites()))
            .to_svg();
        let open = svg.matches("<g>").count() + svg.matches("<g ").count();
        assert_eq!(open, svg.matches("</g>").count(), "{svg}");
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
