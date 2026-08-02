//! Variable sites only, evenly spaced.
//!
//! This is the panel a phylogeny is usually read beside: samples down the side,
//! the sites that tell them apart across the top, and nothing in between. What
//! it buys is legibility and what it costs is the x axis, and both halves of
//! that trade need to be understood before the panel goes into a figure.
//!
//! # Throwing away the columns that agree
//!
//! An alignment of closely related genomes is almost entirely agreement. Thirty
//! kilobases of it might carry twenty differences, and a plot that draws all
//! thirty thousand columns spends every pixel on the twenty-nine thousand nine
//! hundred and eighty that say nothing. So this track throws the invariant
//! columns away and spaces what is left evenly, which turns a smear into
//! twenty legible columns.
//!
//! # The axis is an index, not a coordinate
//!
//! Spacing stops being linear in the genome. Two adjacent columns here may be
//! nine bases apart or nine kilobases apart, and nothing about the spacing says
//! which. That is the trade, and it is why every column carries its own
//! position underneath rather than relying on a ruler: a ruler under this
//! panel would be a lie. Put an [`AxisTrack`](crate::AxisTrack) under something
//! else.
//!
//! The figure's region is therefore the site index space: a panel of twenty
//! sites is `Region::new("sites", 0, 20)`.
//!
//! It follows that the panel places its own columns rather than asking the
//! figure's scale where they go, as most tracks here do. Two things make that
//! necessary. An index axis means nothing to the track above or below, so there
//! is no alignment with a neighbour left to preserve, and the per-sample
//! difference counts need a strip of the band held back on the right, which is
//! width the shared scale knows nothing about and would give to the cells.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::axis::group_thousands;
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

    /// Width of the strip on the right holding the per-sample counts.
    ///
    /// The counts need room of their own. Drawn at the right edge of the
    /// plotting area instead, they land on top of the last site's cells, and a
    /// muted number over a saturated base colour is a number nobody reads.
    fn counts_strip(&self, theme: &Theme) -> f64 {
        if !self.show_counts || self.sites.is_empty() || self.names.is_empty() {
            return 0.0;
        }
        let widest = (0..self.visible_rows().0)
            .map(|row| self.differences(row))
            .max()
            .unwrap_or(0);
        text_width(&widest.to_string(), theme.font_size - 2.0) + 8.0
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
        // The panel lays itself out rather than going through the shared scale.
        // Its x axis is the site index, which means nothing to any neighbour,
        // and the counts on the right need a strip the scale knows nothing of.
        let strip = self.counts_strip(ctx.theme);
        let cell_width = (band.w - strip).max(1.0) / self.sites.len() as f64;
        let x_of = |index: usize| band.x + index as f64 * cell_width;
        let letters = cell_width >= self.letter_threshold;
        let matched = self
            .match_color
            .clone()
            .unwrap_or_else(|| mix(ctx.theme.surface(), &ctx.theme.rule, 0.22));

        // A tint on every other column, so the eye can cross a wide panel of
        // sparse cells without losing its place.
        if self.zebra {
            let tint = mix(ctx.theme.surface(), &ctx.theme.rule, 0.13);
            let height = self.drawn_rows() as f64 * (self.row_height + self.row_gap);
            for index in (0..self.sites.len()).step_by(2) {
                ctx.svg.rect(x_of(index), band.y, cell_width, height, &tint);
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
                    mirror: false,
                },
            );
        }

        let mut top = band.y;
        if self.show_reference {
            for (index, site) in self.sites.iter().enumerate() {
                let x = x_of(index);
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
                let x = x_of(index);
                // The cell is what a pointer lands on, so the cell is what
                // carries the site. Naming only the column label left the
                // panel itself unanswerable, and left `show_positions(false)`
                // with no tooltips at all.
                //
                // Only a cell that differs is named. An agreement is the quiet
                // state this panel exists to see past, it is roughly half the
                // grid, and "matches" is the one thing a reader can already
                // tell from the colour without asking. Naming them cost a
                // third of the file for the least it could say.
                let named = site.differs(row);
                if named {
                    ctx.svg.begin_titled(&self.cell_tooltip(site, row));
                }
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
                    if named {
                        ctx.svg.end_group();
                    }
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
                if named {
                    ctx.svg.end_group();
                }
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
                // The column is the datum here, not the cell: one column is one
                // site, and the cells in it are one site seen once per sample.
                // So the site is named once, on the label that identifies the
                // column, rather than once per cell. The cells cannot carry it
                // between them because they are written row by row, and a
                // column is therefore not a run of elements anything can be
                // wrapped around without reordering the document.
                ctx.svg.begin_titled(&self.site_tooltip(site));
                // Standing the label on end is what lets a five digit position
                // sit under a column narrower than it is.
                ctx.svg.text_rotated(
                    (x_of(index) + cell_width / 2.0 + size * 0.35, top + 2.0),
                    -90.0,
                    &site.position.to_string(),
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
                ctx.svg.end_group();
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
    /// What a reader hovering one column is told.
    ///
    /// What it is, where the site is, what the reference has there, and how
    /// much of the panel departs from it. The last of those is the only one a
    /// reader cannot get by looking: counting coloured cells down a column is
    /// exactly the work a tooltip is for.
    ///
    /// The position is written **exactly as the column's own label writes it**,
    /// ungrouped and without the shift to 1-based coordinates the rest of the
    /// crate applies. This is the crate's one exception to
    /// [`group_thousands`](crate::track::axis::group_thousands) like every
    /// other number the crate writes into a tooltip. The drawn column label
    /// stays ungrouped, which is not a disagreement but the same split
    /// [`VariantTrack`](crate::VariantTrack) makes between `format_value` and
    /// its exact form: a label stood on end in a column narrower than itself
    /// pays for every character, and a tooltip pays for none, so the label
    /// abbreviates and the tooltip is exact.
    ///
    /// The count clause takes the participle every other `N of M` clause in the
    /// crate takes. `differs` and `differ` are a finite verb agreeing with a
    /// number, which put two grammars in one slot for no gain: the tooltip is a
    /// label, and a label does not conjugate.
    fn site_tooltip(&self, site: &SnpSite) -> String {
        let differing = site.count();
        let total = site.alleles.len();
        let mut text = format!("site, {}", group_thousands(site.position));
        text.push_str(", reference ");
        text.push((site.reference as char).to_ascii_uppercase());
        if total > 0 {
            text.push_str(&format!(
                ", {} of {} sample{} differing",
                group_thousands(differing as u64),
                group_thousands(total as u64),
                if total == 1 { "" } else { "s" }
            ));
        }
        text
    }

    /// What a reader hovering one cell is told.
    ///
    /// The cell is the only thing in this panel a pointer can land on, so it
    /// carries the site and then says what this one sample has there, which is
    /// the question a cell exists to answer. A cell that agrees says so rather
    /// than repeating the reference base, since agreement is the quiet state
    /// and naming it twice would bury the difference.
    fn cell_tooltip(&self, site: &SnpSite, row: usize) -> String {
        let mut text = format!("site, {}", group_thousands(site.position));
        text.push_str(", reference ");
        text.push((site.reference as char).to_ascii_uppercase());
        let sample = self.names.get(row).map(String::as_str).unwrap_or("sample");
        match site.allele(row) {
            Some(allele) if site.differs(row) => {
                text.push_str(&format!(
                    ", {sample} has {}",
                    (allele as char).to_ascii_uppercase()
                ));
            }
            _ => text.push_str(&format!(", {sample} matches")),
        }
        text
    }

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
            mix(ctx.theme.surface(), &ctx.theme.foreground, 0.35)
        } else if reference {
            mix(ctx.theme.surface(), &ctx.theme.rule, 0.55)
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
    fn the_counts_get_a_strip_rather_than_landing_on_the_last_column() {
        // Drawn at the right edge of the plotting area, every count sat on top
        // of the last site's cells, and the one over a saturated base colour
        // was unreadable.
        let panel = SnpTrack::from_alignment(0, &alignment());
        let theme = Theme::light();
        assert!(panel.counts_strip(&theme) > 0.0);
        assert_eq!(
            panel.clone().show_counts(false).counts_strip(&theme),
            0.0,
            "no counts, no strip"
        );

        let svg = Figure::new(Region::new("sites", 0, 3).unwrap())
            .show_region_label(false)
            .push(SnpTrack::from_alignment(0, &alignment()).show_positions(false))
            .to_svg();

        // The rightmost cell now stops short of the rightmost count.
        let rights: Vec<f64> = svg
            .match_indices("<rect x=\"")
            .filter_map(|(at, key)| {
                let rest = &svg[at + key.len()..];
                let x: f64 = rest[..rest.find('"')?].parse().ok()?;
                let after = &rest[rest.find("width=\"")? + 7..];
                let w: f64 = after[..after.find('"')?].parse().ok()?;
                (w < 400.0).then_some(x + w)
            })
            .collect();
        let cell_edge = rights.iter().copied().fold(0.0f64, f64::max);
        let count_at = 900.0 - 16.0 - 3.0;
        assert!(
            cell_edge < count_at,
            "cells reach {cell_edge} and the counts start at {count_at}"
        );
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

    /// Every group a track opens has to be closed by exactly one `end_group`.
    fn groups_balance(svg: &str) -> bool {
        svg.matches("<g ").count() + svg.matches("<g>").count() == svg.matches("</g>").count()
    }

    fn panel_svg(panel: SnpTrack) -> String {
        let sites = panel.sites().len().max(1) as u64;
        Figure::new(Region::new("sites", 0, sites).unwrap())
            .show_region_label(false)
            .push(panel)
            .to_svg()
    }

    #[test]
    fn a_cell_says_where_the_site_is_and_what_this_sample_has() {
        let svg = panel_svg(SnpTrack::from_alignment(0, &alignment()).offset(1_472_000));
        // Column 4 of the alignment, where sample_1 and sample_2 both differ.
        // The cell is what a pointer lands on, so the cell is what answers.
        assert!(
            svg.contains("<title>site, 1,472,004, reference A, sample_1 has T</title>"),
            "{svg}"
        );
        assert!(groups_balance(&svg));
    }

    #[test]
    fn only_a_cell_that_differs_is_named() {
        // Agreement is the quiet state this panel exists to see past, it is
        // about half the grid, and the colour already says it. Naming it cost
        // a third of the file for the least it could say.
        let svg = panel_svg(SnpTrack::from_alignment(0, &alignment()).offset(1_472_000));
        assert!(!svg.contains("matches</title>"), "{svg}");
        assert_eq!(svg.matches(" has ").count(), 4, "four cells differ: {svg}");
    }

    #[test]
    fn the_count_clause_takes_a_participle_and_never_conjugates() {
        let svg = panel_svg(SnpTrack::from_alignment(0, &alignment()).offset(1_472_000));
        // Column 7, where only sample_3 carries the gap. One and two read the
        // same way, which is the whole point of a participle: `N of M noun
        // participle` is the idiom every count clause in the crate uses.
        assert!(
            svg.contains("<title>site, 1,472,007, reference T, 1 of 3 samples differing</title>"),
            "{svg}"
        );
        assert!(!svg.contains("differs"), "{svg}");
    }

    #[test]
    fn every_cell_can_be_pointed_at_even_with_the_labels_turned_off() {
        // The tooltip used to hang on the rotated column label, so the panel
        // itself answered nothing and turning the labels off took every
        // tooltip in the track with them.
        let with_labels = panel_svg(SnpTrack::from_alignment(0, &alignment()));
        let without = panel_svg(SnpTrack::from_alignment(0, &alignment()).show_positions(false));

        // Four of the nine cells differ, and each answers for itself, plus one
        // label per site when the labels are drawn.
        assert_eq!(without.matches("<title>").count(), 4, "{without}");
        assert_eq!(with_labels.matches("<title>").count(), 7, "{with_labels}");
        assert!(!without.contains("rotate(-90)"), "no labels are drawn");
        assert!(without.contains(" has "), "the cells still answer");
        assert!(groups_balance(&without));
    }

    #[test]
    fn the_tooltip_is_exact_where_the_drawn_label_abbreviates() {
        // A label stood on end in a narrow column pays for every character and
        // a tooltip pays for none, which is the same split VariantTrack makes.
        let sites = vec![SnpSite::new(9_100, b'G', b"T".to_vec())];
        let svg = panel_svg(SnpTrack::new(vec!["s".to_string()], sites));
        assert!(svg.contains(">9100</text>"), "the label the reader sees");
        assert!(
            svg.contains("<title>site, 9,100, reference G, s has T</title>"),
            "{svg}"
        );
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
