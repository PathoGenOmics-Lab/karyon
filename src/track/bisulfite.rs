//! Methylation one molecule at a time: filled and open circles.
//!
//! A [`Molecule`] is one read and what it said at each of the track's sites,
//! and the track is the grid of them. Everything here follows from keeping the
//! reads apart rather than pooling them into a number per site.
//!
//! # What this says that a fraction cannot
//!
//! [`MethylationTrack`](crate::MethylationTrack) gives a fraction per site: of
//! the reads covering this cytosine, this many were modified. Half the reads
//! methylated at every site in a region is one number, and it has two very
//! different explanations. Either every molecule is methylated at about half
//! its sites, scattered, or half the molecules are methylated at all of them
//! and half at none. The first is a region being modified loosely; the second
//! is two populations of cells, or an allele-specific pattern, and the site
//! fractions are identical in both cases.
//!
//! One row per molecule and one column per site tells them apart at a glance:
//! the first is confetti, the second is stripes.
//!
//! # Filled, open, and nothing
//!
//! Filled for modified, open for not, nothing at all where the molecule did not
//! cover the site. Open and absent must not look the same, which is why an
//! unmethylated call gets a ring and a missing one gets no mark: "measured and
//! not methylated" and "not measured" are different statements.
//!
//! # The columns sit at real distances
//!
//! The x axis is the figure's own, so the columns fall at the real distances
//! between the cytosines. Two sites four bases apart look four bases apart,
//! which matters when the question is whether a CpG island is uniformly
//! modified.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// One sequenced molecule and what it said at each site.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    /// Read or clone name, drawn beside the row.
    pub name: String,
    /// One call per site, in the same order as the sites. `None` is a site the
    /// molecule did not cover, which is not the same as an unmethylated call.
    pub calls: Vec<Option<bool>>,
}

impl Molecule {
    /// A molecule named `name`.
    pub fn new(name: impl Into<String>, calls: impl Into<Vec<Option<bool>>>) -> Self {
        Molecule {
            name: name.into(),
            calls: calls.into(),
        }
    }

    /// What this molecule said at one site.
    pub fn call(&self, site: usize) -> Option<bool> {
        self.calls.get(site).copied().flatten()
    }

    /// How many sites it covered, and how many of those were modified.
    pub fn covered(&self) -> (usize, usize) {
        let covered = self.calls.iter().filter(|call| call.is_some()).count();
        let modified = self
            .calls
            .iter()
            .filter(|call| **call == Some(true))
            .count();
        (covered, modified)
    }

    /// The share of its covered sites that were modified.
    ///
    /// `None` when it covered none, since nothing over nothing is not zero.
    pub fn fraction(&self) -> Option<f64> {
        let (covered, modified) = self.covered();
        (covered > 0).then(|| modified as f64 / covered as f64)
    }
}

/// Methylation per molecule, as the classic grid of filled and open circles.
///
/// ```
/// use karyon::{BisulfiteTrack, Figure, Molecule, Region};
///
/// let sites = vec![1_010u64, 1_024, 1_031, 1_058];
/// let molecules = vec![
///     Molecule::new("read_1", vec![Some(true), Some(true), Some(true), Some(true)]),
///     Molecule::new("read_2", vec![Some(false), Some(false), None, Some(false)]),
/// ];
///
/// let track = BisulfiteTrack::new(sites, molecules);
/// assert_eq!(track.molecules()[0].fraction(), Some(1.0));
///
/// let svg = Figure::new(Region::new("chr1", 1_000, 1_100).unwrap())
///     .push(track.label("CpG"))
///     .to_svg();
/// assert!(svg.contains("<circle"));
/// ```
#[derive(Debug, Clone)]
pub struct BisulfiteTrack {
    sites: Vec<u64>,
    molecules: Vec<Molecule>,
    label: Option<String>,
    row_height: f64,
    radius: f64,
    color: Option<String>,
    show_names: bool,
    show_rule: bool,
    max_rows: Option<usize>,
}

impl BisulfiteTrack {
    /// A grid over `sites`, one row per molecule.
    pub fn new(sites: impl Into<Vec<u64>>, molecules: impl Into<Vec<Molecule>>) -> Self {
        BisulfiteTrack {
            sites: sites.into(),
            molecules: molecules.into(),
            label: None,
            row_height: 12.0,
            radius: 3.6,
            color: None,
            show_names: true,
            show_rule: true,
            max_rows: Some(40),
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one molecule's row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(4.0);
        self
    }

    /// Sets the radius of one circle.
    pub fn radius(mut self, radius: f64) -> Self {
        self.radius = radius.max(1.0);
        self
    }

    /// Sets the colour the circles are drawn in.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Draws or hides the molecule names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Draws or hides the faint line along each molecule.
    ///
    /// It is what makes a row a molecule rather than a scatter of circles that
    /// happen to be level with each other.
    pub fn show_rule(mut self, show: bool) -> Self {
        self.show_rule = show;
        self
    }

    /// Caps how many molecules are drawn.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|count| count.max(1));
        self
    }

    /// The site positions.
    pub fn sites(&self) -> &[u64] {
        &self.sites
    }

    /// The molecules.
    pub fn molecules(&self) -> &[Molecule] {
        &self.molecules
    }

    /// How many rows are drawn and how many the cap hid.
    pub fn visible_rows(&self) -> (usize, usize) {
        match self.max_rows {
            Some(cap) if self.molecules.len() > cap => (cap, self.molecules.len() - cap),
            _ => (self.molecules.len(), 0),
        }
    }

    /// The share of covered calls at one site that were modified.
    ///
    /// `None` when no molecule covered it.
    pub fn site_fraction(&self, site: usize) -> Option<f64> {
        let calls: Vec<bool> = self
            .molecules
            .iter()
            .filter_map(|molecule| molecule.call(site))
            .collect();
        (!calls.is_empty())
            .then(|| calls.iter().filter(|call| **call).count() as f64 / calls.len() as f64)
    }

    /// How far the molecules are from all being the same.
    ///
    /// Zero when every molecule is methylated at the same share of its sites,
    /// and larger the more they differ. It is the number that separates the two
    /// readings of one site fraction: loose modification of every molecule
    /// gives a small spread, two populations give a large one.
    ///
    /// `None` when fewer than two molecules covered anything.
    pub fn discordance(&self) -> Option<f64> {
        let fractions: Vec<f64> = self
            .molecules
            .iter()
            .filter_map(|molecule| molecule.fraction())
            .collect();
        if fractions.len() < 2 {
            return None;
        }
        let mean = fractions.iter().sum::<f64>() / fractions.len() as f64;
        let variance = fractions
            .iter()
            .map(|value| (value - mean).powi(2))
            .sum::<f64>()
            / fractions.len() as f64;
        Some(variance.sqrt())
    }
}

impl Track for BisulfiteTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let (rows, _) = self.visible_rows();
        rows.max(1) as f64 * self.row_height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_names || self.molecules.is_empty() {
            return 0.0;
        }
        let size = (theme.font_size - 2.0).min(self.row_height);
        let (rows, _) = self.visible_rows();
        self.molecules
            .iter()
            .take(rows)
            .map(|molecule| text_width(&molecule.name, size))
            .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let (rows, hidden) = self.visible_rows();
        if self.sites.is_empty() || rows == 0 {
            return;
        }
        let ink = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.foreground.clone());
        let rule = mix(ctx.theme.surface(), &ctx.theme.rule, 0.8);
        let name_size = (ctx.theme.font_size - 2.0).min(self.row_height);

        let first = self.sites.iter().copied().min().unwrap_or(0);
        let last = self.sites.iter().copied().max().unwrap_or(0);

        for (row, molecule) in self.molecules.iter().take(rows).enumerate() {
            let y = band.y + (row as f64 + 0.5) * self.row_height;

            if self.show_rule {
                // The line is what makes a row a molecule rather than a scatter
                // of circles that happen to be level with each other.
                ctx.svg
                    .line(ctx.scale.x(first), y, ctx.scale.x(last), y, &rule, 0.8);
            }

            for (site, position) in self.sites.iter().enumerate() {
                let Some(modified) = molecule.call(site) else {
                    // Not covered. Nothing is drawn, because an open circle
                    // already means "measured and not methylated" and the two
                    // must not look the same.
                    continue;
                };
                let x = ctx.scale.x_center(*position);
                if modified {
                    ctx.svg.circle(x, y, self.radius, &ink);
                } else {
                    // Filled with the page rather than left hollow, so the rule
                    // under it does not run through the middle of the ring.
                    ctx.svg.circle(x, y, self.radius, ctx.theme.surface());
                    ctx.svg
                        .circle_ringed(x, y, self.radius - 0.9, ctx.theme.surface(), &ink, 0.9);
                }
            }

            if self.show_names && ctx.axis.w > 0.0 {
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    y + name_size * 0.35,
                    &molecule.name,
                    &ctx.theme.muted,
                    name_size,
                    Anchor::End,
                );
            }
        }

        if hidden > 0 {
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 2.0,
                &format!("+{hidden} more"),
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

    fn sites() -> Vec<u64> {
        vec![1_010u64, 1_024, 1_031, 1_058]
    }

    /// Two populations: two molecules methylated throughout, two at nothing.
    fn stripes() -> Vec<Molecule> {
        vec![
            Molecule::new("read_1", vec![Some(true); 4]),
            Molecule::new("read_2", vec![Some(true); 4]),
            Molecule::new("read_3", vec![Some(false); 4]),
            Molecule::new("read_4", vec![Some(false); 4]),
        ]
    }

    /// Loose modification: every molecule half methylated, scattered.
    fn confetti() -> Vec<Molecule> {
        vec![
            Molecule::new(
                "read_1",
                vec![Some(true), Some(false), Some(true), Some(false)],
            ),
            Molecule::new(
                "read_2",
                vec![Some(false), Some(true), Some(false), Some(true)],
            ),
            Molecule::new(
                "read_3",
                vec![Some(true), Some(false), Some(false), Some(true)],
            ),
            Molecule::new(
                "read_4",
                vec![Some(false), Some(true), Some(true), Some(false)],
            ),
        ]
    }

    fn region() -> Region {
        Region::new("chr1", 1_000, 1_100).unwrap()
    }

    #[test]
    fn a_site_the_molecule_did_not_cover_is_not_an_unmethylated_call() {
        let molecule = Molecule::new("r", vec![Some(true), None, Some(false)]);
        assert_eq!(molecule.call(0), Some(true));
        assert_eq!(molecule.call(1), None, "not covered");
        assert_eq!(molecule.call(2), Some(false), "covered and not methylated");
        assert_eq!(molecule.covered(), (2, 1));
        assert_eq!(molecule.fraction(), Some(0.5));
    }

    #[test]
    fn a_molecule_covering_nothing_has_no_fraction() {
        // Nothing over nothing is not zero.
        assert_eq!(Molecule::new("r", vec![None, None]).fraction(), None);
        assert_eq!(Molecule::new("r", Vec::new()).covered(), (0, 0));
    }

    #[test]
    fn two_patterns_with_one_site_fraction_are_told_apart() {
        // Both are half methylated at every site. That is the whole point: the
        // site fractions cannot separate them and the molecules can.
        let striped = BisulfiteTrack::new(sites(), stripes());
        let scattered = BisulfiteTrack::new(sites(), confetti());
        for site in 0..4 {
            assert_eq!(striped.site_fraction(site), Some(0.5));
            assert_eq!(scattered.site_fraction(site), Some(0.5));
        }
        // And the discordance separates them: all or nothing per molecule
        // against half of each.
        assert_eq!(striped.discordance(), Some(0.5));
        assert_eq!(scattered.discordance(), Some(0.0));
    }

    #[test]
    fn discordance_needs_two_molecules_to_be_anything() {
        let lonely = BisulfiteTrack::new(sites(), vec![Molecule::new("r", vec![Some(true); 4])]);
        assert_eq!(lonely.discordance(), None);
        let uncovered = BisulfiteTrack::new(sites(), vec![Molecule::new("r", vec![None; 4]); 3]);
        assert_eq!(uncovered.discordance(), None);
    }

    #[test]
    fn a_site_nothing_covered_has_no_fraction() {
        let sparse = vec![Molecule::new("r", vec![Some(true), None, None, None])];
        let track = BisulfiteTrack::new(sites(), sparse);
        assert_eq!(track.site_fraction(0), Some(1.0));
        assert_eq!(track.site_fraction(1), None);
        assert_eq!(track.site_fraction(9), None, "past the last site");
    }

    #[test]
    fn a_modified_call_is_filled_and_an_unmodified_one_is_a_ring() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), stripes()).color("#111111"))
            .to_svg();
        // Eight filled circles, and eight rings drawn as two circles each.
        assert_eq!(svg.matches("<circle").count(), 8 + 8 * 3);
        assert!(svg.contains("#111111"));
    }

    #[test]
    fn a_site_nobody_covered_draws_nothing_at_all() {
        // An open circle already means "measured and not methylated", so a
        // missing call cannot be one too.
        let with = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), stripes()).show_names(false))
            .to_svg();
        let gapped: Vec<Molecule> = stripes()
            .into_iter()
            .map(|mut molecule| {
                molecule.calls[1] = None;
                molecule
            })
            .collect();
        let without = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), gapped).show_names(false))
            .to_svg();
        assert!(without.matches("<circle").count() < with.matches("<circle").count());
    }

    #[test]
    fn each_row_carries_a_line_so_it_reads_as_a_molecule() {
        let with = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), stripes()))
            .to_svg();
        let without = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), stripes()).show_rule(false))
            .to_svg();
        assert_eq!(with.matches("<line").count(), 4);
        assert_eq!(without.matches("<line").count(), 0);
    }

    #[test]
    fn the_columns_sit_at_the_real_distances_between_the_sites() {
        // Two sites seven bases apart have to look seven bases apart, which is
        // the question when a CpG island is being read.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(sites(), stripes()).show_names(false))
            .to_svg();
        let xs: Vec<f64> = svg
            .match_indices("<circle cx=\"")
            .filter_map(|(at, key)| {
                let rest = &svg[at + key.len()..];
                rest[..rest.find('"')?].parse().ok()
            })
            .collect();
        let mut unique: Vec<f64> = xs.clone();
        unique.sort_by(|a, b| a.partial_cmp(b).unwrap());
        unique.dedup();
        assert_eq!(unique.len(), 4, "one column per site");
        // 1,024 to 1,031 is seven bases; 1,031 to 1,058 is twenty-seven.
        let gaps: Vec<f64> = unique.windows(2).map(|pair| pair[1] - pair[0]).collect();
        assert!(gaps[2] > gaps[1] * 3.0, "{gaps:?}");
    }

    #[test]
    fn a_deep_pile_of_molecules_is_capped_and_says_so() {
        let many: Vec<Molecule> = (0..80)
            .map(|index| Molecule::new(format!("r{index}"), vec![Some(true); 4]))
            .collect();
        let track = BisulfiteTrack::new(sites(), many).max_rows(Some(10));
        assert_eq!(track.visible_rows(), (10, 70));
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("+70 more"));
    }

    #[test]
    fn an_empty_track_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(BisulfiteTrack::new(Vec::new(), Vec::new()).label("CpG"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
