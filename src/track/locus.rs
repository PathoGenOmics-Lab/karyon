//! Several loci from several genomes, with their homologies drawn between.
//!
//! # The idea
//!
//! A locus of a few dozen features, whether that is a gene cluster, an operon,
//! a viral genome or a syntenic block, is rarely looked at on its own: the
//! question asked of it is almost never "what is in it" but "what is in it that
//! the other one has not". Answering
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
use crate::svg::{num, text_width, Anchor};
use crate::theme::{contrast_ink, mix, wash, Theme};
use crate::track::feature::strand_color;

/// How much of a gene's height the shaft of an arrow takes, the rest being the
/// overhang of its head.
const SHAFT: f64 = 0.62;
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

/// The outline a gene is drawn with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GeneShape {
    /// A shaft with a head that overhangs it, the way an arrow is drawn.
    ///
    /// The default. The head is a shape in its own right rather than a taper
    /// on the end of a bar, so which way a gene points is legible at a glance
    /// and stays legible when the gene is long.
    #[default]
    Arrow,
    /// A block with one end brought to a point.
    ///
    /// Squarer, and worth having when the genes are short: a head that
    /// overhangs needs a shaft to overhang, and a gene a few pixels wide has
    /// none.
    Pointed,
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
    mark_unmatched: bool,
    soft_fills: bool,
    link_inset: f64,
    shape: GeneShape,
}

impl LocusTrack {
    /// A track over `loci`, drawn top to bottom in the order given.
    pub fn new(loci: impl Into<Vec<Locus>>) -> Self {
        LocusTrack {
            loci: loci.into(),
            links: Vec::new(),
            label: None,
            gene_height: 22.0,
            link_height: 34.0,
            show_names: true,
            show_gene_names: true,
            color: None,
            reverse_color: None,
            min_gene_width: 2.0,
            identity_range: (0.7, 1.0),
            mark_unmatched: true,
            soft_fills: true,
            link_inset: 4.0,
            shape: GeneShape::Arrow,
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

    /// Chooses the outline a gene is drawn with.
    pub fn shape(mut self, shape: GeneShape) -> Self {
        self.shape = shape;
        self
    }

    /// Whether a gene is a wash of its colour edged in the colour itself.
    ///
    /// On by default. A gene arrow is a large filled shape, and a large shape
    /// at full saturation shouts: eight of them and the figure is a colour
    /// chart. Putting the hue in the edge and a wash of it in the body keeps
    /// the identity and gives the page back. It also lets the name inside be
    /// dark ink, which reads better than white on a saturated block.
    ///
    /// Turn it off to have the colour used exactly as given.
    pub fn soft_fills(mut self, soft: bool) -> Self {
        self.soft_fills = soft;
        self
    }

    /// Sets the gap between a gene and the ribbon leaving it, in pixels.
    ///
    /// Ribbons that touch the arrows read as one shape with them. A few pixels
    /// of page is what makes a ribbon a connection between two things.
    pub fn link_inset(mut self, inset: f64) -> Self {
        self.link_inset = inset.max(0.0);
        self
    }

    /// Outlines the genes no homology reaches.
    ///
    /// On by default, because it is the question the figure exists to answer.
    /// The missing ribbon says it too, but only to a reader who thought to
    /// look for an absence, and an absence is the hardest thing to notice.
    pub fn mark_unmatched(mut self, mark: bool) -> Self {
        self.mark_unmatched = mark;
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

    /// The two ends of the ribbon shading, palest first.
    ///
    /// For handing to a [`Legend`](crate::Legend), so the key and the ribbons
    /// cannot drift apart: a legend that names its own colours is a legend
    /// that goes stale the first time the ramp is touched.
    pub fn ramp_ends(&self, theme: &Theme) -> (String, String) {
        let (low, high) = self.identity_range;
        (
            mix(theme.surface(), &theme.foreground, self.shade(low)),
            mix(theme.surface(), &theme.foreground, self.shade(high)),
        )
    }

    /// How dark the ribbon for a given identity is drawn.
    fn shade(&self, identity: f64) -> f64 {
        let (low, high) = self.identity_range;
        let fraction = ((identity - low) / (high - low)).clamp(0.0, 1.0);
        // Spread widely enough that two identities a few per cent apart are
        // two different greys, and capped low enough that the darkest is still
        // background. A perfect match is already the widest ribbon on the
        // page; making it the heaviest as well buries the gaps, and the gaps
        // are what the figure is about.
        0.05 + 0.17 * fraction
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
            let top = band.y + self.row_top(link.row) + self.gene_height + self.link_inset;
            let bottom = band.y + self.row_top(link.row + 1) - self.link_inset;
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
            // Curved sides rather than straight ones. A quadrilateral between
            // two rows reads as a block of colour; a ribbon that leaves each
            // gene vertically and arrives vertically reads as a connection,
            // and tells the eye which end goes with which without being
            // traced corner to corner.
            let (ax0, ax1) = (ctx.scale.x(upper.start), ctx.scale.x(upper.end));
            let (bx0, bx1) = (ctx.scale.x(lower.start), ctx.scale.x(lower.end));
            let waist = (top + bottom) / 2.0;
            let d = format!(
                "M{} {}L{} {}C{} {} {} {} {} {}L{} {}C{} {} {} {} {} {}Z",
                num(ax0),
                num(top),
                num(ax1),
                num(top),
                num(ax1),
                num(waist),
                num(bx1),
                num(waist),
                num(bx1),
                num(bottom),
                num(bx0),
                num(bottom),
                num(bx0),
                num(waist),
                num(ax0),
                num(waist),
                num(ax0),
                num(top),
            );
            ctx.svg.path(&d, &shade, 1.0);
            // A hairline edge, and no more. At these tints a dark outline
            // fights the genes, and the genes are the subject.
            ctx.svg.path_stroked(&d, &mix(&shade, "#000000", 0.10), 0.5);
        }

        for (row, locus) in self.loci.iter().enumerate() {
            let top = band.y + self.row_top(row);
            // The sequence the genes sit on, drawn from the first to the last
            // of them. Without it a row is a line of floating arrows and the
            // space between two genes is nothing at all; with it that space is
            // intergenic sequence, which is what it actually is.
            if let Some((from, to)) = locus.span() {
                let middle = top + self.gene_height / 2.0;
                ctx.svg.line(
                    ctx.scale.x(from),
                    middle,
                    ctx.scale.x(to),
                    middle,
                    &ctx.theme.rule,
                    1.0,
                );
            }
            let orphans = if self.mark_unmatched {
                self.unmatched(row)
            } else {
                Vec::new()
            };
            // A name written under a gene has to clear the one before it.
            // Two adjacent short genes would otherwise print their names on
            // top of each other, which is worse than printing neither.
            let mut label_right = f64::NEG_INFINITY;
            for (index, gene) in locus.genes.iter().enumerate() {
                let ink = ctx.theme.foreground.clone();
                let unmatched = orphans.contains(&index);
                label_right = self.draw_gene(ctx, gene, top, unmatched, ink, label_right);
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
    fn draw_gene(
        &self,
        ctx: &mut DrawContext<'_>,
        gene: &Feature,
        top: f64,
        unmatched: bool,
        ink: String,
        label_right: f64,
    ) -> f64 {
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

        // A head is a third of the gene, capped at a little over its height.
        // Fixed to the height alone it is a nick on the end of a long gene,
        // which reads as a bar rather than as an arrow; unbounded it swallows
        // a short gene and says nothing about which way that one points.
        let head = ((x1 - x0) * 0.32).min(height * 1.15).min(x1 - x0).max(0.0);
        let body = (x1 - x0 - head).max(0.0);
        let mid = top + height / 2.0;
        let bottom = top + height;
        // Half the shaft, which the head overhangs on both sides. A gene too
        // short to have a shaft is drawn as the head alone.
        let shaft = height * SHAFT / 2.0;
        // A head that overhangs needs a shaft long enough to be seen as one.
        // Below that the arrow is all head, and a two pixel gene with a one
        // pixel shaft is a smudge with a notch in it.
        let barbed = self.shape == GeneShape::Arrow && body >= height * 0.4;
        let points: Vec<(f64, f64)> = match (barbed, gene.strand) {
            (true, Strand::Reverse) => vec![
                (x1, mid - shaft),
                (x0 + head, mid - shaft),
                (x0 + head, top),
                (x0, mid),
                (x0 + head, bottom),
                (x0 + head, mid + shaft),
                (x1, mid + shaft),
            ],
            (true, _) => vec![
                (x0, mid - shaft),
                (x0 + body, mid - shaft),
                (x0 + body, top),
                (x1, mid),
                (x0 + body, bottom),
                (x0 + body, mid + shaft),
                (x0, mid + shaft),
            ],
            (false, Strand::Reverse) => vec![
                (x1, top),
                (x0 + head, top),
                (x0, mid),
                (x0 + head, bottom),
                (x1, bottom),
            ],
            (false, _) => vec![
                (x0, top),
                (x0 + body, top),
                (x1, mid),
                (x0 + body, bottom),
                (x0, bottom),
            ],
        };
        // The hue goes in the edge and a wash of it in the body. It costs a
        // little ink and buys the shape: two abutting genes of one family stop
        // being one long blob, an arrowhead against a ribbon of similar weight
        // stays an arrowhead, and a large area stops shouting.
        let (fill, edge_color) = if self.soft_fills {
            (wash(&color, ctx.theme), color.clone())
        } else {
            (color.clone(), mix(&color, "#000000", 0.3))
        };
        ctx.svg.polygon(&points, &fill);
        let mut edge: Vec<(f64, f64)> = points.clone();
        edge.push(points[0]);
        ctx.svg.polyline(&edge, &edge_color, 1.1);
        if unmatched {
            // Nothing in the neighbouring rows matched this one, which is the
            // finding. Outlined rather than recoloured, so it keeps whatever
            // family colour it came in with.
            let mut outline: Vec<(f64, f64)> = points.clone();
            outline.push(points[0]);
            ctx.svg.polyline(&outline, &ink, 1.6);
        }

        if self.show_gene_names {
            if let Some(name) = &gene.name {
                let room = if barbed { height * SHAFT } else { height };
                let size = (room * 0.72).min(ctx.theme.font_size);
                if text_width(name, size) < body.max(1.0) - 4.0 {
                    ctx.svg.text(
                        (x0 + x1) / 2.0,
                        mid + size * 0.35,
                        name,
                        contrast_ink(&fill),
                        size,
                        Anchor::Middle,
                    );
                } else if self.link_height >= size + 4.0 {
                    // Under the gene rather than nowhere. A short gene is
                    // exactly the one whose name a reader needs, and dropping
                    // it leaves the figure quietly incomplete. It still has to
                    // clear the last one written down there.
                    let width = text_width(name, size);
                    let left = (x0 + x1) / 2.0 - width / 2.0;
                    if left > label_right + 3.0 {
                        ctx.svg.text(
                            (x0 + x1) / 2.0,
                            top + height + size,
                            name,
                            &ctx.theme.muted,
                            size,
                            Anchor::Middle,
                        );
                        return left + width;
                    }
                }
            }
        }
        label_right
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
        assert_eq!(track.height(&scale()), 2.0 * 22.0 + 34.0);
        assert_eq!(
            LocusTrack::new(loci()).link_height(0.0).height(&scale()),
            44.0
        );
        // An empty track still has a row of height rather than none.
        assert_eq!(LocusTrack::new(Vec::new()).height(&scale()), 22.0);
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
        // Five genes as polygons, and the two ribbons as curved paths, each
        // filled and then edged.
        assert_eq!(svg.matches("<polygon").count(), 5);
        assert_eq!(svg.matches("<path").count(), 4);
        let first_gene = svg.find("<polygon").unwrap();
        let first_ribbon = svg.find("<path").unwrap();
        assert!(first_ribbon < first_gene, "a gene half under a ribbon");
    }

    #[test]
    fn the_identity_ramp_covers_the_range_the_data_uses() {
        // Orthologues sit between about seventy and a hundred per cent, so a
        // ramp over the whole of nought to one draws them all the same shade.
        let track = LocusTrack::new(loci());
        assert!(
            track.shade(1.0) - track.shade(0.7) > 0.15,
            "the default range"
        );
        // Everything below the floor is the palest shade, not a negative one.
        assert_eq!(track.shade(0.0), track.shade(0.7));
        // A ribbon stays context: never more than a third of the page's ink.
        assert!(track.shade(1.0) < 0.25, "still context, not the subject");

        // Spread over the whole range instead, and a seventy per cent match is
        // already most of the way to solid, leaving nothing for the rest.
        let wide = LocusTrack::new(loci()).identity_range(0.0, 1.0);
        assert!(wide.shade(0.7) > track.shade(0.7));
        assert!(
            wide.shade(1.0) - wide.shade(0.7) < (track.shade(1.0) - track.shade(0.7)) / 2.0,
            "the useful range gets less than half the ramp"
        );
        // A range the wrong way round falls back rather than dividing by zero.
        let broken = LocusTrack::new(loci()).identity_range(0.9, 0.1);
        assert!(broken.shade(0.5).is_finite());
    }

    #[test]
    fn the_legend_can_take_the_ramp_from_the_track() {
        // A legend that names its own colours goes stale the first time the
        // ramp is touched.
        let theme = Theme::light();
        let track = LocusTrack::new(loci());
        let (pale, dark) = track.ramp_ends(&theme);
        assert_ne!(pale, dark);

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(vec![Homology::new(0, 0, 0, 1.0)]))
            .to_svg();
        assert!(
            svg.contains(&dark),
            "the ramp's dark end is not on the page"
        );
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
    fn a_gene_nothing_matched_is_outlined() {
        // The missing ribbon says it too, but only to a reader who thought to
        // look for an absence.
        let marked = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(links()))
            .to_svg();
        let plain = Figure::new(region())
            .show_region_label(false)
            .push(LocusTrack::new(loci()).links(links()).mark_unmatched(false))
            .to_svg();
        // Five genes, each with an edge of its own, and espC in the top row
        // carries a second outline because nothing below it matched.
        assert_eq!(marked.matches("<polyline").count(), 6);
        assert_eq!(plain.matches("<polyline").count(), 5);
        let ink = Theme::light().foreground;
        assert!(marked.contains(&format!("stroke=\"{ink}\"")));
        assert!(!plain.contains(&format!("stroke=\"{ink}\"")));
    }

    #[test]
    fn an_arrow_has_a_head_that_overhangs_its_shaft() {
        let paint = |shape: GeneShape| {
            Figure::new(region())
                .show_region_label(false)
                .push(
                    LocusTrack::new(vec![Locus::new(
                        "a",
                        vec![Feature::new(0, 2_400).strand(Strand::Forward)],
                    )])
                    .shape(shape),
                )
                .to_svg()
        };
        let corners = |svg: &str| {
            let at = svg.find("<polygon points=\"").unwrap() + 17;
            svg[at..].split('"').next().unwrap().split(' ').count()
        };
        // Seven corners for a shaft with a head on it, five for a block
        // brought to a point.
        assert_eq!(corners(&paint(GeneShape::Arrow)), 7);
        assert_eq!(corners(&paint(GeneShape::Pointed)), 5);
    }

    #[test]
    fn a_gene_too_short_for_a_shaft_is_drawn_as_the_head_alone() {
        // A head that overhangs needs a shaft to overhang, and a gene a few
        // pixels wide has none.
        let svg = Figure::new(Region::new("x", 0, 400_000).unwrap())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 200).strand(Strand::Forward)],
            )]))
            .to_svg();
        let at = svg.find("<polygon points=\"").unwrap() + 17;
        let corners = svg[at..].split('"').next().unwrap().split(' ').count();
        assert_eq!(corners, 5, "no room for a shaft");
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
    fn a_gene_name_wider_than_its_gene_goes_underneath_it() {
        // A short gene is exactly the one whose name a reader needs, so it is
        // moved rather than dropped.
        let cramped = Figure::new(Region::new("x", 0, 400_000).unwrap())
            .show_region_label(false)
            .push(LocusTrack::new(vec![Locus::new(
                "a",
                vec![Feature::new(0, 900).name("a_very_long_gene_name")],
            )]))
            .to_svg();
        assert!(cramped.contains(">a_very_long_gene_name</text>"));
        // In the muted ink under the gene, not the contrast ink inside it.
        assert!(cramped.contains(&format!("fill=\"{}\"", Theme::light().muted)));

        // With no room between the rows there is nowhere to put it.
        let flat = Figure::new(Region::new("x", 0, 400_000).unwrap())
            .show_region_label(false)
            .push(
                LocusTrack::new(vec![Locus::new(
                    "a",
                    vec![Feature::new(0, 900).name("a_very_long_gene_name")],
                )])
                .link_height(0.0),
            )
            .to_svg();
        assert!(!flat.contains(">a_very_long_gene_name</text>"));
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
