//! Variable sites only, evenly spaced.
//!
//! # The idea
//!
//! An alignment of closely related genomes is almost entirely agreement. Thirty
//! kilobases of it might carry twenty differences, and a plot that draws all
//! thirty thousand columns spends every pixel on the twenty-nine thousand nine
//! hundred and eighty that say nothing. So this track throws the invariant
//! columns away and spaces what is left evenly, which turns a smear into
//! twenty legible columns.
//!
//! # What that costs
//!
//! The x axis stops being linear in the genome. Two adjacent columns here may
//! be nine bases apart or nine kilobases apart, and nothing about the spacing
//! says which. That is the trade, and it is why every column carries its own
//! position underneath rather than relying on a ruler: a ruler under this
//! panel would be a lie. Put an [`AxisTrack`](crate::AxisTrack) under something
//! else.
//!
//! The figure's region is therefore the site index space: a panel of twenty
//! sites is `Region::new("sites", 0, 20)`.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::msa::{is_gap, MsaSequence};
use crate::track::tree::{draw_tree, leaf_order, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// One position where the sequences disagree.
#[derive(Debug, Clone, PartialEq)]
pub struct SnpSite {
    /// Where the site is in whatever coordinates the caller is using.
    pub position: u64,
    /// The reference residue.
    pub reference: u8,
    /// What each sample has here, in row order.
    pub alleles: Vec<u8>,
}

impl SnpSite {
    /// A site at `position`.
    pub fn new(position: u64, reference: u8, alleles: impl Into<Vec<u8>>) -> Self {
        SnpSite {
            position,
            reference,
            alleles: alleles.into(),
        }
    }

    /// What sample `row` has here.
    pub fn allele(&self, row: usize) -> Option<u8> {
        self.alleles.get(row).copied()
    }

    /// Whether sample `row` differs from the reference here.
    ///
    /// A gap is a difference, since a deletion is a real observation, but an
    /// absent row is not: past the end of the allele list there is nothing to
    /// compare.
    pub fn differs(&self, row: usize) -> bool {
        match self.allele(row) {
            Some(allele) => !allele.eq_ignore_ascii_case(&self.reference),
            None => false,
        }
    }

    /// How many samples differ from the reference here.
    pub fn count(&self) -> usize {
        (0..self.alleles.len())
            .filter(|row| self.differs(*row))
            .count()
    }
}

/// A panel of variable sites, one row per sample.
///
/// ```
/// use karyon::{Figure, MsaSequence, Region, SnpTrack};
///
/// let rows = vec![
///     MsaSequence::new("reference", b"ACGTACGTAC".to_vec()),
///     MsaSequence::new("sample_1", b"ACGTTCGTAC".to_vec()),
///     MsaSequence::new("sample_2", b"ACGTTCGTGC".to_vec()),
/// ];
///
/// // Two columns of ten differ, so the panel is two wide.
/// let panel = SnpTrack::from_alignment(0, &rows);
/// assert_eq!(panel.sites().len(), 2);
///
/// let svg = Figure::new(Region::new("sites", 0, 2).unwrap())
///     .push(panel.label("isolates"))
///     .to_svg();
/// assert!(svg.contains("sample_1"));
/// ```
#[derive(Debug, Clone)]
pub struct SnpTrack {
    names: Vec<String>,
    sites: Vec<SnpSite>,
    reference_name: String,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    show_names: bool,
    show_reference: bool,
    show_positions: bool,
    show_counts: bool,
    zebra: bool,
    letter_threshold: f64,
    max_rows: Option<usize>,
    match_color: Option<String>,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
}

impl SnpTrack {
    /// A panel from sites you have already worked out.
    pub fn new(names: impl Into<Vec<String>>, sites: impl Into<Vec<SnpSite>>) -> Self {
        SnpTrack {
            names: names.into(),
            sites: sites.into(),
            reference_name: "reference".to_string(),
            label: None,
            row_height: 15.0,
            row_gap: 2.0,
            show_names: true,
            show_reference: true,
            show_positions: true,
            show_counts: true,
            zebra: true,
            letter_threshold: 7.0,
            max_rows: Some(40),
            match_color: None,
            tree: None,
            tree_width: 90.0,
            tree_shape: TreeShape::Phylogram,
        }
    }

    /// A panel from an alignment, keeping only the columns that vary.
    ///
    /// Row `reference` is the one everything is compared against and is left
    /// out of the sample rows. A column counts as variable when any other row
    /// disagrees with it, gaps included, since a deletion is an observation
    /// too. Positions are alignment column indices; use
    /// [`SnpTrack::offset`] when the alignment starts somewhere other than
    /// zero, or build the sites yourself when gaps make the mapping
    /// interesting.
    pub fn from_alignment(reference: usize, sequences: &[MsaSequence]) -> Self {
        let Some(reference_row) = sequences.get(reference) else {
            return SnpTrack::new(Vec::new(), Vec::new());
        };
        let samples: Vec<&MsaSequence> = sequences
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != reference)
            .map(|(_, row)| row)
            .collect();

        let columns = sequences
            .iter()
            .map(|row| row.residues.len())
            .max()
            .unwrap_or(0);

        let mut sites = Vec::new();
        for column in 0..columns {
            let Some(base) = reference_row.residue(column) else {
                continue;
            };
            let alleles: Vec<u8> = samples
                .iter()
                .map(|row| row.residue(column).unwrap_or(b'-'))
                .collect();
            let varies = alleles
                .iter()
                .any(|allele| !allele.eq_ignore_ascii_case(&base));
            if varies {
                sites.push(SnpSite::new(column as u64, base, alleles));
            }
        }

        let mut panel = SnpTrack::new(
            samples
                .iter()
                .map(|row| row.name.clone())
                .collect::<Vec<_>>(),
            sites,
        );
        panel.reference_name = reference_row.name.clone();
        panel
    }

    /// Shifts every site position by `offset`, for an alignment that starts
    /// somewhere other than zero.
    pub fn offset(mut self, offset: u64) -> Self {
        for site in &mut self.sites {
            site.position += offset;
        }
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(2.0);
        self
    }

    /// Sets the gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Names the reference row.
    pub fn reference_name(mut self, name: impl Into<String>) -> Self {
        self.reference_name = name.into();
        self
    }

    /// Draws or hides the sample names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Draws or hides the reference row along the top.
    pub fn show_reference(mut self, show: bool) -> Self {
        self.show_reference = show;
        self
    }

    /// Draws or hides the position under each column.
    ///
    /// They are the only thing saying where these sites actually are, since the
    /// spacing deliberately does not.
    pub fn show_positions(mut self, show: bool) -> Self {
        self.show_positions = show;
        self
    }

    /// Draws or hides the per-sample difference count on the right.
    pub fn show_counts(mut self, show: bool) -> Self {
        self.show_counts = show;
        self
    }

    /// Draws or hides the alternating column tint.
    ///
    /// It is there to carry the eye across a wide panel of sparse cells; turn
    /// it off for a narrow one where it is just decoration.
    pub fn zebra(mut self, show: bool) -> Self {
        self.zebra = show;
        self
    }

    /// Caps how many sample rows are drawn.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|r| r.max(1));
        self
    }

    /// Sets the colour of a cell that matches the reference.
    pub fn match_color(mut self, color: impl Into<String>) -> Self {
        self.match_color = Some(color.into());
        self
    }

    /// Draws a phylogeny beside the panel and sorts the rows to match it.
    ///
    /// This is the pairing the tree exists for. Sorted by name, a clade's
    /// shared substitutions are scattered down the panel and look like noise;
    /// sorted by descent they line up into a block, and the block is the
    /// finding. The tree goes in the strip the figure already reserves for row
    /// names, so it costs no new machinery.
    ///
    /// Rows are matched to leaves by name. A sample the tree does not mention
    /// keeps its place at the bottom rather than disappearing, because a row
    /// silently dropped from a figure is worse than a row out of order.
    pub fn tree(mut self, tree: Tree) -> Self {
        let sorted = leaf_order(&tree, &self.names);
        self.names = sorted
            .iter()
            .map(|index| self.names[*index].clone())
            .collect();
        for site in &mut self.sites {
            site.alleles = sorted
                .iter()
                .map(|index| site.alleles.get(*index).copied().unwrap_or(b'-'))
                .collect();
        }
        self.tree = Some(tree);
        self
    }

    /// Sets how much of the strip the tree gets, in pixels.
    pub fn tree_width(mut self, width: f64) -> Self {
        self.tree_width = width.max(0.0);
        self
    }

    /// Chooses a phylogram or a cladogram for the tree beside the panel.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// The tree beside the panel, if there is one.
    pub fn attached_tree(&self) -> Option<&Tree> {
        self.tree.as_ref()
    }

    /// The variable sites, in position order.
    pub fn sites(&self) -> &[SnpSite] {
        &self.sites
    }

    /// The sample names.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// How many sites sample `row` differs at.
    pub fn differences(&self, row: usize) -> usize {
        self.sites.iter().filter(|site| site.differs(row)).count()
    }

    /// How many rows are drawn and how many the cap hid.
    pub fn visible_rows(&self) -> (usize, usize) {
        match self.max_rows {
            Some(cap) if self.names.len() > cap => (cap, self.names.len() - cap),
            _ => (self.names.len(), 0),
        }
    }

    /// Height of the strip under the panel holding the position labels.
    fn position_strip(&self, theme: &Theme) -> f64 {
        if !self.show_positions || self.sites.is_empty() {
            return 0.0;
        }
        let size = theme.font_size - 2.0;
        let widest = self
            .sites
            .iter()
            .map(|site| text_width(&site.position.to_string(), size))
            .fold(0.0f64, f64::max);
        widest + 6.0
    }

    /// How many rows of cells, reference included.
    fn drawn_rows(&self) -> usize {
        let (rows, _) = self.visible_rows();
        rows + usize::from(self.show_reference)
    }
}

impl Track for SnpTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.drawn_rows().max(1) as f64;
        rows * self.row_height
            + (rows - 1.0).max(0.0) * self.row_gap
            + self.position_strip(&Theme::default())
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let tree = if self.tree.is_some() {
            self.tree_width
        } else {
            0.0
        };
        if !self.show_names {
            return tree;
        }
        let size = (theme.font_size - 1.0).min(self.row_height);
        let (rows, _) = self.visible_rows();
        let widest = self
            .names
            .iter()
            .take(rows)
            .chain(std::iter::once(&self.reference_name))
            .map(|name| text_width(name, size))
            .fold(0.0f64, f64::max);
        if widest <= 0.0 {
            tree
        } else {
            widest + 8.0 + tree
        }
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        if self.sites.is_empty() {
            return;
        }
        let (rows, hidden) = self.visible_rows();
        let name_size = (ctx.theme.font_size - 1.0).min(self.row_height);
        let cell_width = ctx.scale.x(1) - ctx.scale.x(0);
        let letters = cell_width >= self.letter_threshold;
        let matched = self
            .match_color
            .clone()
            .unwrap_or_else(|| mix(&ctx.theme.background, &ctx.theme.rule, 0.22));

        // A tint on every other column, so the eye can cross a wide panel of
        // sparse cells without losing its place.
        if self.zebra {
            let tint = mix(&ctx.theme.background, &ctx.theme.rule, 0.13);
            let strip = self.drawn_rows() as f64 * (self.row_height + self.row_gap);
            for index in (0..self.sites.len()).step_by(2) {
                ctx.svg
                    .rect(ctx.scale.x(index as u64), band.y, cell_width, strip, &tint);
            }
        }

        // The tree takes the left of the strip and the names the right of it,
        // so a leaf, its name and its row of cells are all on one line.
        if let Some(tree) = &self.tree {
            let reference_offset = if self.show_reference {
                self.row_height + self.row_gap
            } else {
                0.0
            };
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: band.y + reference_offset,
                w: (self.tree_width - 6.0).max(1.0),
                h: band.h - reference_offset,
            };
            draw_tree(
                ctx.svg,
                tree,
                area,
                self.row_height + self.row_gap,
                area.y + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &ctx.theme.foreground,
                    width: 1.1,
                },
            );
        }

        let mut top = band.y;
        if self.show_reference {
            for (index, site) in self.sites.iter().enumerate() {
                let x = ctx.scale.x(index as u64);
                self.paint_cell(
                    ctx,
                    Rect {
                        x,
                        y: top,
                        w: cell_width,
                        h: self.row_height,
                    },
                    site.reference,
                    letters,
                    true,
                );
            }
            self.paint_name(ctx, top, name_size, &self.reference_name.clone(), true);
            top += self.row_height + self.row_gap;
        }

        for row in 0..rows {
            for (index, site) in self.sites.iter().enumerate() {
                let x = ctx.scale.x(index as u64);
                if !site.differs(row) {
                    // An agreement is a quiet bar. The whole point of the panel
                    // is that the disagreements are the only thing worth ink.
                    ctx.svg.rect_rounded(
                        x + 1.0,
                        top + self.row_height * 0.32,
                        (cell_width - 2.0).max(1.0),
                        self.row_height * 0.36,
                        ctx.theme.corner_radius,
                        &matched,
                    );
                    continue;
                }
                let Some(allele) = site.allele(row) else {
                    continue;
                };
                self.paint_cell(
                    ctx,
                    Rect {
                        x,
                        y: top,
                        w: cell_width,
                        h: self.row_height,
                    },
                    allele,
                    letters,
                    false,
                );
            }

            if let Some(name) = self.names.get(row) {
                self.paint_name(ctx, top, name_size, &name.clone(), false);
            }
            if self.show_counts {
                let count = self.differences(row);
                ctx.svg.text(
                    band.right() - 3.0,
                    top + self.row_height / 2.0 + name_size * 0.35,
                    &count.to_string(),
                    &ctx.theme.muted,
                    name_size - 1.0,
                    Anchor::End,
                );
            }
            top += self.row_height + self.row_gap;
        }

        if self.show_positions {
            let size = ctx.theme.font_size - 2.0;
            for (index, site) in self.sites.iter().enumerate() {
                // Standing the label on end is what lets a five digit position
                // sit under a column narrower than it is.
                ctx.svg.text_rotated(
                    (
                        ctx.scale.x(index as u64) + cell_width / 2.0 + size * 0.35,
                        top + 2.0,
                    ),
                    -90.0,
                    &site.position.to_string(),
                    &ctx.theme.muted,
                    size,
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
                ctx.theme.font_size - 1.0,
                Anchor::End,
            );
        }
    }
}

impl SnpTrack {
    /// Paints one residue cell, with its letter when there is room.
    fn paint_cell(
        &self,
        ctx: &mut DrawContext<'_>,
        cell: Rect,
        residue: u8,
        letters: bool,
        reference: bool,
    ) {
        let (x, top, width) = (cell.x, cell.y, cell.w);
        let color = if is_gap(residue) {
            mix(&ctx.theme.background, &ctx.theme.foreground, 0.35)
        } else if reference {
            mix(&ctx.theme.background, &ctx.theme.rule, 0.55)
        } else {
            ctx.theme.bases.of(residue).to_string()
        };
        ctx.svg.rect_rounded(
            x + 1.0,
            top,
            (width - 2.0).max(1.0),
            self.row_height,
            ctx.theme.corner_radius,
            &color,
        );
        if letters && !is_gap(residue) {
            let size = (self.row_height * 0.62).min(width * 0.7);
            ctx.svg.glyph(
                x + (width - size * 0.75) / 2.0,
                top + self.row_height - (self.row_height - size) / 2.0,
                size * 0.75,
                size / ctx.theme.cap_height_ratio.max(0.1),
                &(residue as char).to_ascii_uppercase().to_string(),
                contrast_ink(&color),
            );
        }
    }

    /// Draws one row name in the axis strip.
    fn paint_name(
        &self,
        ctx: &mut DrawContext<'_>,
        top: f64,
        size: f64,
        name: &str,
        reference: bool,
    ) {
        if !self.show_names || ctx.axis.w <= 0.0 {
            return;
        }
        let color = if reference {
            ctx.theme.foreground.clone()
        } else {
            ctx.theme.muted.clone()
        };
        ctx.svg.text(
            ctx.axis.right() - 4.0,
            top + self.row_height / 2.0 + size * 0.35,
            name,
            &color,
            size,
            Anchor::End,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn alignment() -> Vec<MsaSequence> {
        vec![
            MsaSequence::new("reference", b"ACGTACGTAC".to_vec()),
            MsaSequence::new("sample_1", b"ACGTTCGTAC".to_vec()),
            MsaSequence::new("sample_2", b"ACGTTCGTGC".to_vec()),
            MsaSequence::new("sample_3", b"ACGTACG-AC".to_vec()),
        ]
    }

    #[test]
    fn only_the_columns_that_vary_survive() {
        let panel = SnpTrack::from_alignment(0, &alignment());
        let positions: Vec<u64> = panel.sites().iter().map(|s| s.position).collect();
        // Columns 4, 7 and 8 differ; the other seven agree and are dropped.
        assert_eq!(positions, vec![4, 7, 8]);
    }

    #[test]
    fn the_reference_row_is_not_one_of_the_samples() {
        let panel = SnpTrack::from_alignment(0, &alignment());
        assert_eq!(panel.names(), ["sample_1", "sample_2", "sample_3"]);
        assert_eq!(panel.reference_name, "reference");
        assert_eq!(panel.sites()[0].alleles.len(), 3);
    }

    #[test]
    fn any_row_can_be_the_reference() {
        let panel = SnpTrack::from_alignment(1, &alignment());
        assert_eq!(panel.reference_name, "sample_1");
        assert!(panel.names().contains(&"reference".to_string()));
    }

    #[test]
    fn a_gap_is_a_difference_because_a_deletion_is_an_observation() {
        let panel = SnpTrack::from_alignment(0, &alignment());
        let site = panel.sites().iter().find(|s| s.position == 7).unwrap();
        // Only sample_3 has the gap.
        assert_eq!(site.count(), 1);
        assert!(site.differs(2));
        assert!(!site.differs(0));
    }

    #[test]
    fn agreement_is_case_insensitive() {
        let rows = vec![
            MsaSequence::new("reference", b"ACGT".to_vec()),
            MsaSequence::new("soft_masked", b"acgt".to_vec()),
        ];
        assert!(
            SnpTrack::from_alignment(0, &rows).sites().is_empty(),
            "a soft masked copy is not a pile of SNPs"
        );
    }

    #[test]
    fn a_missing_reference_row_gives_an_empty_panel() {
        let panel = SnpTrack::from_alignment(99, &alignment());
        assert!(panel.sites().is_empty());
        assert!(panel.names().is_empty());
    }

    #[test]
    fn an_alignment_that_agrees_everywhere_has_no_panel_to_draw() {
        let rows = vec![
            MsaSequence::new("a", b"ACGT".to_vec()),
            MsaSequence::new("b", b"ACGT".to_vec()),
        ];
        assert!(SnpTrack::from_alignment(0, &rows).sites().is_empty());
    }

    #[test]
    fn a_short_row_is_padded_with_gaps_rather_than_dropped() {
        let rows = vec![
            MsaSequence::new("reference", b"ACGT".to_vec()),
            MsaSequence::new("truncated", b"AC".to_vec()),
        ];
        let panel = SnpTrack::from_alignment(0, &rows);
        // Columns 2 and 3 are missing from the short row and count as gaps.
        assert_eq!(panel.sites().len(), 2);
        assert_eq!(panel.differences(0), 2);
    }

    #[test]
    fn counts_are_per_sample_and_per_site() {
        let panel = SnpTrack::from_alignment(0, &alignment());
        assert_eq!(panel.differences(0), 1, "sample_1 has one");
        assert_eq!(panel.differences(1), 2, "sample_2 has two");
        assert_eq!(panel.differences(2), 1, "sample_3 has the gap");
        assert_eq!(panel.sites()[0].count(), 2, "column 4 differs in two rows");
    }

    #[test]
    fn an_offset_moves_every_position_together() {
        let panel = SnpTrack::from_alignment(0, &alignment()).offset(10_000);
        let positions: Vec<u64> = panel.sites().iter().map(|s| s.position).collect();
        assert_eq!(positions, vec![10_004, 10_007, 10_008]);
    }

    #[test]
    fn the_panel_is_indexed_by_site_not_by_position() {
        // Two sites nine thousand bases apart sit next to each other, which is
        // the whole trade this track makes.
        let sites = vec![
            SnpSite::new(100, b'A', b"C".to_vec()),
            SnpSite::new(9_100, b'G', b"T".to_vec()),
        ];
        let panel = SnpTrack::new(vec!["s".to_string()], sites);
        let svg = Figure::new(Region::new("sites", 0, 2).unwrap())
            .show_region_label(false)
            .push(panel)
            .to_svg();
        // Both positions are written out, since the spacing cannot say them.
        assert!(svg.contains(">100</text>"));
        assert!(svg.contains(">9100</text>"));
        assert!(svg.contains("rotate(-90)"));
    }

    #[test]
    fn a_deep_panel_is_capped_and_says_so() {
        let names: Vec<String> = (0..60).map(|i| format!("s{i}")).collect();
        let sites = vec![SnpSite::new(1, b'A', vec![b'C'; 60])];
        let panel = SnpTrack::new(names, sites).max_rows(Some(8));
        assert_eq!(panel.visible_rows(), (8, 52));

        let svg = Figure::new(Region::new("sites", 0, 1).unwrap())
            .show_region_label(false)
            .push(panel)
            .to_svg();
        assert!(svg.contains("more"));
    }

    #[test]
    fn names_include_the_reference_when_sizing_the_strip() {
        let theme = Theme::light();
        let short = SnpTrack::from_alignment(0, &alignment()).reference_name("r");
        let long =
            SnpTrack::from_alignment(0, &alignment()).reference_name("a_very_long_reference_name");
        assert!(long.y_axis_width(&theme) > short.y_axis_width(&theme));
    }

    #[test]
    fn the_reference_row_can_be_turned_off_and_the_panel_shrinks() {
        let scale = Scale::new(&Region::new("sites", 0, 3).unwrap(), 0.0, 100.0);
        let with = SnpTrack::from_alignment(0, &alignment());
        let without = SnpTrack::from_alignment(0, &alignment()).show_reference(false);
        assert!(with.height(&scale) > without.height(&scale));
    }

    #[test]
    fn the_position_strip_only_takes_room_when_it_is_used() {
        let theme = Theme::light();
        let with = SnpTrack::from_alignment(0, &alignment());
        let without = SnpTrack::from_alignment(0, &alignment()).show_positions(false);
        assert!(with.position_strip(&theme) > 0.0);
        assert_eq!(without.position_strip(&theme), 0.0);
    }

    #[test]
    fn an_empty_panel_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(Region::new("sites", 0, 1).unwrap())
            .show_region_label(false)
            .push(SnpTrack::new(Vec::new(), Vec::new()).label("none"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
