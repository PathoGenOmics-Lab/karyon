//! Several loci from several genomes, with their homologies drawn between.
//!
//! # The idea
//!
//! A genomic island, an operon, a prophage or a biosynthetic cluster is a
//! stretch of a few dozen genes, and the question asked of it is almost never
//! "what is in it" but "what is in it that the other one has not". Answering
//! that needs the loci one under another, the genes drawn as genes, and a line
//! from each gene to whatever it matches in the row below. What the reader is
//! then looking for is the gaps: the ribbon that goes nowhere, the arrow that
//! points the other way, the block that arrived as a unit.
//!
//! # What the x axis is
//!
//! The figure's own, shared with every other track, so a kilobase is a kilobase
//! in every row and the loci can be compared for length as well as for content.
//! Give each [`Locus`] genes in whatever coordinates they came in, and shift a
//! row with [`Locus::offset`] to line it up with its neighbour.
//!
//! Homologies are between neighbouring rows only. That is not a limitation of
//! the drawing so much as of the reading: a ribbon that skips a row crosses one
//! it has nothing to do with, and a figure of those is a figure of crossings.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::feature::strand_color;
use crate::track::{DrawContext, Feature, Strand, Track};

/// One locus, from one genome.
#[derive(Debug, Clone, PartialEq)]
pub struct Locus {
    /// Which genome this row came from.
    pub name: String,
    /// The genes, in the coordinates they arrived in.
    pub genes: Vec<Feature>,
}

impl Locus {
    /// A locus named `name`.
    pub fn new(name: impl Into<String>, genes: impl Into<Vec<Feature>>) -> Self {
        Locus {
            name: name.into(),
            genes: genes.into(),
        }
    }

    /// Shifts every gene in the locus by `offset` bases.
    ///
    /// For lining one row up with the row above it, which is usually what turns
    /// a tangle of crossing ribbons into a set of parallel ones.
    pub fn offset(mut self, offset: i64) -> Self {
        for gene in &mut self.genes {
            gene.start = gene.start.saturating_add_signed(offset);
            gene.end = gene.end.saturating_add_signed(offset);
        }
        self
    }

    /// First and last base the locus covers.
    pub fn span(&self) -> Option<(u64, u64)> {
        let start = self.genes.iter().map(|gene| gene.start).min()?;
        let end = self.genes.iter().map(|gene| gene.end).max()?;
        Some((start, end))
    }
}

/// A gene in one row matching a gene in the row below it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Homology {
    /// Index of the upper row. The lower one is always the next.
    pub row: usize,
    /// Index of the gene in the upper row.
    pub from: usize,
    /// Index of the gene in the lower row.
    pub to: usize,
    /// How alike they are, between 0 and 1, which sets how dark the ribbon is.
    pub identity: f64,
}

impl Homology {
    /// A match between `from` in row `row` and `to` in the row under it.
    pub fn new(row: usize, from: usize, to: usize, identity: f64) -> Self {
        Homology {
            row,
            from,
            to,
            identity: identity.clamp(0.0, 1.0),
        }
    }
}

/// Loci stacked with their homologies between them.
///
/// ```
/// use karyon::{Feature, Figure, Homology, Locus, LocusTrack, Region, Strand};
///
/// let loci = vec![
///     Locus::new("H37Rv", vec![
///         Feature::new(0, 1_200).name("espA").strand(Strand::Forward),
///         Feature::new(1_300, 2_100).name("espC").strand(Strand::Forward),
///     ]),
///     Locus::new("CDC1551", vec![
///         Feature::new(0, 1_200).name("espA").strand(Strand::Forward),
///     ]),
/// ];
///
/// let svg = Figure::new(Region::new("ESX-1", 0, 2_400).unwrap())
///     .push(LocusTrack::new(loci).links(vec![Homology::new(0, 0, 0, 0.98)]))
///     .to_svg();
/// assert!(svg.contains("H37Rv"));
/// ```
#[derive(Debug, Clone)]
pub struct LocusTrack {
    loci: Vec<Locus>,
    links: Vec<Homology>,
    label: Option<String>,
    gene_height: f64,
    link_height: f64,
    show_names: bool,
    show_gene_names: bool,
    color: Option<String>,
    reverse_color: Option<String>,
    min_gene_width: f64,
    identity_range: (f64, f64),
}

impl LocusTrack {
    /// A track over `loci`, drawn top to bottom in the order given.
    pub fn new(loci: impl Into<Vec<Locus>>) -> Self {
        LocusTrack {
            loci: loci.into(),
            links: Vec::new(),
            label: None,
            gene_height: 16.0,
            link_height: 26.0,
            show_names: true,
            show_gene_names: true,
            color: None,
            reverse_color: None,
            min_gene_width: 2.0,
            identity_range: (0.7, 1.0),
        }
    }

    /// Adds the homologies between neighbouring rows.
    ///
    /// A link naming a row or a gene that is not there is dropped rather than
    /// drawn somewhere arbitrary.
    pub fn links(mut self, links: impl Into<Vec<Homology>>) -> Self {
        self.links = links.into();
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets how tall one gene is drawn.
    pub fn gene_height(mut self, height: f64) -> Self {
        self.gene_height = height.max(4.0);
        self
    }

    /// Sets how much room the ribbons get between two rows.
    pub fn link_height(mut self, height: f64) -> Self {
        self.link_height = height.max(0.0);
        self
    }

    /// Draws or hides the genome names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Draws or hides the gene names, where there is room for them.
    pub fn show_gene_names(mut self, show: bool) -> Self {
        self.show_gene_names = show;
        self
    }

    /// Sets the colours of the forward and reverse genes.
    ///
    /// Only used for a gene that carries no colour of its own. Giving a gene
    /// family one colour across every row is what makes a rearrangement
    /// visible without following a single ribbon.
    pub fn colors(mut self, forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        self.color = Some(forward.into());
        self.reverse_color = Some(reverse.into());
        self
    }

    /// Sets how narrow a gene may be drawn, in pixels.
    pub fn min_gene_width(mut self, pixels: f64) -> Self {
        self.min_gene_width = pixels.max(0.5);
        self
    }

    /// Sets the identities that map to the palest and darkest ribbon.
    ///
    /// Real homologies do not run from nought to one. A set of orthologues sits
    /// between about seventy and a hundred per cent identical, so a ramp spread
    /// over the whole range would draw all of them the same shade and say
    /// nothing. Narrow it to the range the data actually occupies.
    pub fn identity_range(mut self, low: f64, high: f64) -> Self {
        let (low, high) = (low.clamp(0.0, 1.0), high.clamp(0.0, 1.0));
        self.identity_range = if high > low { (low, high) } else { (0.0, 1.0) };
        self
    }

    /// How dark the ribbon for a given identity is drawn.
    fn shade(&self, identity: f64) -> f64 {
        let (low, high) = self.identity_range;
        let fraction = ((identity - low) / (high - low)).clamp(0.0, 1.0);
        // Kept pale on purpose. A ribbon is context for the genes, and at half
        // the ink of the page it stops being context and becomes the subject.
        0.08 + 0.24 * fraction
    }

    /// The loci.
    pub fn loci(&self) -> &[Locus] {
        &self.loci
    }

    /// The homologies.
    pub fn homologies(&self) -> &[Homology] {
        &self.links
    }

    /// A gene in a row, if both exist.
    pub fn gene(&self, row: usize, index: usize) -> Option<&Feature> {
        self.loci.get(row)?.genes.get(index)
    }

    /// Genes in `row` that no homology reaches.
    ///
    /// These are the answer to the question the figure was drawn to ask: what
    /// this locus has that its neighbours do not. Rows at the two ends are
    /// judged against their one neighbour, so a gene is unmatched only when
    /// nothing it could have matched did.
    pub fn unmatched(&self, row: usize) -> Vec<usize> {
        let Some(locus) = self.loci.get(row) else {
            return Vec::new();
        };
        (0..locus.genes.len())
            .filter(|index| {
                !self.links.iter().any(|link| {
                    (link.row == row && link.from == *index)
                        || (link.row + 1 == row && link.to == *index)
                })
            })
            .collect()
    }

    /// Top of a row inside the band.
    fn row_top(&self, row: usize) -> f64 {
        row as f64 * (self.gene_height + self.link_height)
    }
}

impl Track for LocusTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.loci.len().max(1) as f64;
        rows * self.gene_height + (rows - 1.0).max(0.0) * self.link_height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_names || self.loci.is_empty() {
            return 0.0;
        }
        let size = (theme.font_size - 1.0).min(self.gene_height);
        self.loci
            .iter()
            .map(|locus| text_width(&locus.name, size))
            .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        if self.loci.is_empty() {
            return;
        }

        // Ribbons first, so a gene is never half hidden under one.
        for link in &self.links {
            let (Some(upper), Some(lower)) = (
                self.gene(link.row, link.from),
                self.gene(link.row + 1, link.to),
            ) else {
                continue;
            };
            let top = band.y + self.row_top(link.row) + self.gene_height;
            let bottom = band.y + self.row_top(link.row + 1);
            if bottom <= top {
                continue;
            }
            // A grey ramp rather than a hue, because the colour of a gene is
            // already carrying which family it belongs to and two colour
            // scales in one plot are one too many.
            let shade = mix(
                ctx.theme.surface(),
                &ctx.theme.foreground,
                self.shade(link.identity),
            );
            let points = [
                (ctx.scale.x(upper.start), top),
                (ctx.scale.x(upper.end), top),
                (ctx.scale.x(lower.end), bottom),
                (ctx.scale.x(lower.start), bottom),
            ];
            ctx.svg.polygon(&points, &shade);
        }

        for (row, locus) in self.loci.iter().enumerate() {
            let top = band.y + self.row_top(row);
            for gene in &locus.genes {
                self.draw_gene(ctx, gene, top);
            }
            if self.show_names && ctx.axis.w > 0.0 {
                let size = (ctx.theme.font_size - 1.0).min(self.gene_height);
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    top + self.gene_height / 2.0 + size * 0.35,
                    &locus.name,
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }
    }
}

impl LocusTrack {
    /// One gene, as an arrow pointing the way it is transcribed.
    fn draw_gene(&self, ctx: &mut DrawContext<'_>, gene: &Feature, top: f64) {
        let x0 = ctx.scale.x(gene.start);
        let x1 = ctx.scale.x(gene.end).max(x0 + self.min_gene_width);
        let height = self.gene_height;
        let color = gene.color.clone().unwrap_or_else(|| {
            if gene.strand == Strand::Reverse {
                self.reverse_color
                    .clone()
                    .unwrap_or_else(|| strand_color(Strand::Reverse, ctx.theme).to_string())
            } else {
                self.color
                    .clone()
                    .unwrap_or_else(|| strand_color(Strand::Forward, ctx.theme).to_string())
            }
        });

        // The head takes a third of the gene, or less when the gene is short:
        // a head longer than its body reads as an arrowhead on its own and
        // says nothing about which way the gene points.
        let head = (x1 - x0).min(height * 0.7).max(0.0);
        let body = (x1 - x0 - head).max(0.0);
        let mid = top + height / 2.0;
        let points: Vec<(f64, f64)> = match gene.strand {
            Strand::Reverse => vec![
                (x1, top),
                (x0 + head, top),
                (x0, mid),
                (x0 + head, top + height),
                (x1, top + height),
            ],
            _ => vec![
                (x0, top),
                (x0 + body, top),
                (x1, mid),
                (x0 + body, top + height),
                (x0, top + height),
            ],
        };
        ctx.svg.polygon(&points, &color);

        if self.show_gene_names {
            if let Some(name) = &gene.name {
                let size = (height * 0.62).min(ctx.theme.font_size);
                if text_width(name, size) < body.max(1.0) - 2.0 {
                    ctx.svg.text(
                        (x0 + x1) / 2.0,
                        mid + size * 0.35,
                        name,
                        contrast_ink(&color),
                        size,
                        Anchor::Middle,
                    );
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

    fn loci() -> Vec<Locus> {
        vec![
            Locus::new(
                "H37Rv",
                vec![
                    Feature::new(0, 1_200).name("espA").strand(Strand::Forward),
                    Feature::new(1_300, 2_100)
                        .name("espC")
                        .strand(Strand::Forward),
                    Feature::new(2_200, 3_000)
                        .name("espD")
                        .strand(Strand::Reverse),
                ],
            ),
            Locus::new(
                "CDC1551",
                vec![
                    Feature::new(0, 1_200).name("espA").strand(Strand::Forward),
                    Feature::new(1_400, 2_200)
                        .name("espD")
                        .strand(Strand::Reverse),
                ],
            ),
        ]
    }

    fn links() -> Vec<Homology> {
        vec![Homology::new(0, 0, 0, 0.99), Homology::new(0, 2, 1, 0.91)]
    }

    fn region() -> Region {
        Region::new("ESX-1", 0, 3_200).unwrap()
    }

    fn scale() -> Scale {
        Scale::new(&region(), 0.0, 800.0)
    }

    #[test]
    fn a_locus_knows_what_it_covers() {
        assert_eq!(loci()[0].span(), Some((0, 3_000)));
        assert_eq!(Locus::new("empty", Vec::new()).span(), None);
    }

    #[test]
    fn an_offset_moves_every_gene_together() {
        let shifted = loci()[0].clone().offset(5_000);
        assert_eq!(shifted.span(), Some((5_000, 8_000)));
        // And back, without running off the bottom of the coordinate system.
        let clamped = Locus::new("a", vec![Feature::new(10, 20)]).offset(-1_000);
        assert_eq!(clamped.span(), Some((0, 0)));
    }

    #[test]
    fn an_identity_outside_its_range_is_brought_back_in() {
        assert_eq!(Homology::new(0, 0, 0, 4.0).identity, 1.0);
        assert_eq!(Homology::new(0, 0, 0, -1.0).identity, 0.0);
    }

    #[test]
    fn height_follows_the_row_count_and_the_room_between_them() {
        let track = LocusTrack::new(loci());
        assert_eq!(track.height(&scale()), 2.0 * 16.0 + 26.0);
        assert_eq!(
            LocusTrack::new(loci()).link_height(0.0).height(&scale()),
            32.0
        );
        // An empty track still has a row of height rather than none.
        assert_eq!(LocusTrack::new(Vec::new()).height(&scale()), 16.0);
    }

    #[test]
    fn an_unmatched_gene_is_what_the_figure_is_for() {
        let track = LocusTrack::new(loci()).links(links());
        // espC is in H37Rv and matched by nothing below it.
        assert_eq!(track.unmatched(0), vec![1]);
        // Both genes of the second row are matched from above.
        assert!(track.unmatched(1).is_empty());
    }

    #[test]
    fn everything_is_unmatched_without_any_homologies() {
        let track = LocusTrack::new(loci());
        assert_eq!(track.unmatched(0), vec![0, 1, 2]);
        assert_eq!(track.unmatched(1), vec![0, 1]);
        assert!(track.unmatched(9).is_empty(), "a row that is not there");
    }

    #[test]
    fn a_link_to_a_gene_that_is_not_there_is_dropped() {
        let track = LocusTrack::new(loci()).links(vec![
            Homology::new(0, 99, 0, 1.0),
            Homology::new(5, 0, 0, 1.0),
        ]);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        // Five genes and no ribbons: the arrows are the only polygons.
        assert_eq!(svg.matches("<polygon").count(), 5);
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn ribbons_go_under_the_genes_rather_than_over_them() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(links()))
            .to_svg();
        // Two ribbons then five genes, in that order.
        assert_eq!(svg.matches("<polygon").count(), 7);
        let first_gene = svg.find(Theme::light().color(0)).unwrap();
        let first_ribbon = svg.find("<polygon").unwrap();
        assert!(first_ribbon < first_gene, "a gene half under a ribbon");
    }

    #[test]
    fn the_identity_ramp_covers_the_range_the_data_uses() {
        // Orthologues sit between about seventy and a hundred per cent, so a
        // ramp over the whole of nought to one draws them all the same shade.
        let track = LocusTrack::new(loci());
        assert!(
            track.shade(1.0) - track.shade(0.7) > 0.2,
            "the default range"
        );
        // Everything below the floor is the palest shade, not a negative one.
        assert_eq!(track.shade(0.0), track.shade(0.7));
        // A ribbon stays context: never more than a third of the page's ink.
        assert!(track.shade(1.0) < 0.34);

        // Spread over the whole range instead, and a seventy per cent match is
        // already most of the way to solid, leaving nothing for the rest.
        let wide = LocusTrack::new(loci()).identity_range(0.0, 1.0);
        assert!(wide.shade(0.7) > track.shade(0.7));
        assert!(wide.shade(1.0) - wide.shade(0.7) < 0.1);
        // A range the wrong way round falls back rather than dividing by zero.
        let broken = LocusTrack::new(loci()).identity_range(0.9, 0.1);
        assert!(broken.shade(0.5).is_finite());
    }

    #[test]
    fn a_closer_match_is_a_darker_ribbon() {
        let dark = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(vec![Homology::new(0, 0, 0, 1.0)]))
            .to_svg();
        let pale = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(vec![Homology::new(0, 0, 0, 0.0)]))
            .to_svg();
        assert_ne!(dark, pale);
    }

    #[test]
    fn a_gene_points_the_way_it_is_transcribed() {
        let forward = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 2_000).strand(Strand::Forward)],
            )]))
            .to_svg();
        let reverse = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 2_000).strand(Strand::Reverse)],
            )]))
            .to_svg();
        assert_ne!(forward, reverse, "the arrow is the strand");
    }

    #[test]
    fn a_gene_keeps_its_own_colour_when_it_has_one() {
        // One colour per gene family across every row is what makes a
        // rearrangement visible without following a single ribbon.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 900).color("#123456")],
            )]))
            .to_svg();
        assert!(svg.contains("#123456"));
    }

    #[test]
    fn names_are_drawn_and_size_the_strip() {
        let theme = Theme::light();
        let track = LocusTrack::new(loci());
        assert!(track.y_axis_width(&theme) > 0.0);
        assert_eq!(track.clone().show_names(false).y_axis_width(&theme), 0.0);

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()))
            .to_svg();
        assert!(svg.contains(">H37Rv</text>"));
        assert!(svg.contains(">CDC1551</text>"));
        assert!(svg.contains(">espA</text>"), "and the gene names too");
    }

    #[test]
    fn a_gene_name_wider_than_its_gene_is_left_out() {
        let cramped = Figure::new(Region::new("x", 0, 400_000).unwrap())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 900).name("a_very_long_gene_name")],
            )]))
            .to_svg();
        assert!(!cramped.contains("a_very_long_gene_name</text>"));
    }

    #[test]
    fn an_empty_track_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(Vec::new()).label("none"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
