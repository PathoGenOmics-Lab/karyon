//! Stretches of DNA that belong to a clade rather than to a sample.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, wash, Theme};
use crate::track::tree::{draw_tree, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// A contiguous stretch of reference carried by a named set of taxa.
#[derive(Debug, Clone, PartialEq)]
pub struct CladeBlock {
    start: u64,
    end: u64,
    taxa: Vec<String>,
    name: Option<String>,
    color: Option<String>,
}

impl CladeBlock {
    /// A block covering `start..end`, 0-based half-open, carried by `taxa`.
    pub fn new(start: u64, end: u64, taxa: impl IntoIterator<Item = impl Into<String>>) -> Self {
        let (start, end) = if end >= start {
            (start, end)
        } else {
            (end, start)
        };
        CladeBlock {
            start,
            end,
            taxa: taxa.into_iter().map(Into::into).collect(),
            name: None,
            color: None,
        }
    }

    /// Sets the name drawn on the block when it is wide enough to hold it.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Overrides the fill colour, which is otherwise taken from the palette.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Reference start, 0-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Reference end, exclusive.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// How much sequence the block covers.
    pub fn span(&self) -> u64 {
        self.end - self.start
    }

    /// The taxa carrying it.
    pub fn taxa(&self) -> &[String] {
        &self.taxa
    }
}

/// Where a block sits once its taxa have been found in the tree.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Placed {
    first: usize,
    last: usize,
    carriers: usize,
    missing: usize,
}

impl Placed {
    /// Rows between the first and last carrier that do not carry the block.
    fn rows(&self) -> usize {
        self.last - self.first + 1
    }

    /// Whether the carriers are every leaf between the first and the last, which
    /// is what makes the block one event on one branch.
    fn is_clade(&self) -> bool {
        self.rows() == self.carriers
    }
}

/// Genomic intervals painted onto a phylogeny, one block per inferred event.
///
/// A [`MatrixTrack`](crate::MatrixTrack) cell is one base wide and cells never
/// merge, so a matrix can only say that these six samples each carry something
/// here. That is six observations. This says one: a rectangle whose vertical
/// extent covers a whole clade asserts a single acquisition or loss on the branch
/// below which every carrier sits, and the same six samples drawn as six separate
/// rows assert six independent events. The difference between those two claims is
/// most of what a comparative genomics figure is arguing about.
///
/// Horizontal extent is a real coordinate span, so the questions a reader asks of
/// it are coordinate questions: does the block cover *katG*, do two blocks on
/// different branches share an endpoint, does the deletion stop where an IS copy
/// sits.
///
/// The track can only lie in one way, and it refuses to. When the carriers are
/// not every leaf under their common ancestor, the block is still drawn across
/// the rows it spans, but every row inside it that does not carry the block is
/// cut out, so a paraphyletic set can never pass for a clade.
///
/// ```
/// use karyon::{CladeBlock, CladeTrack, Figure, Region};
/// use karyon::tree::Tree;
///
/// let tree = Tree::parse_newick("(((A:1,B:1):1,C:2):1,D:3);").unwrap();
/// let track = CladeTrack::new(tree, vec![
///     CladeBlock::new(1_000, 4_000, ["A", "B"]).name("RD1"),
/// ]);
/// assert!(track.is_clade(0), "A and B are sisters");
///
/// let svg = Figure::new(Region::new("chr", 0, 10_000).unwrap())
///     .push(track.label("lineages"))
///     .to_svg();
/// assert!(svg.contains("RD1"));
/// ```
#[derive(Debug, Clone)]
pub struct CladeTrack {
    tree: Tree,
    blocks: Vec<CladeBlock>,
    order: Vec<usize>,
    placed: Vec<Option<Placed>>,
    row_height: f64,
    row_gap: f64,
    label: Option<String>,
    tree_width: f64,
    tree_shape: TreeShape,
    show_names: bool,
    show_block_names: bool,
    min_block: f64,
}

impl CladeTrack {
    /// A track over `tree`, painting `blocks` onto it.
    pub fn new(tree: Tree, blocks: Vec<CladeBlock>) -> Self {
        let rows = row_of_leaf(&tree);
        let placed: Vec<Option<Placed>> = blocks
            .iter()
            .map(|block| place(block, &rows, tree.leaf_count()))
            .collect();

        // Widest clade first, so a block on a tip paints over the one on the
        // branch that carries it rather than under it. StructuralTrack sorts its
        // arcs by the same rule and for the same reason.
        let mut order: Vec<usize> = (0..blocks.len()).collect();
        order.sort_by(|a, b| {
            let rows_of = |index: usize| placed[index].map(|p| p.rows()).unwrap_or(0);
            rows_of(*b)
                .cmp(&rows_of(*a))
                .then(blocks[*b].span().cmp(&blocks[*a].span()))
                .then(a.cmp(b))
        });

        CladeTrack {
            tree,
            blocks,
            order,
            placed,
            row_height: 13.0,
            row_gap: 6.0,
            label: None,
            tree_width: 100.0,
            tree_shape: TreeShape::Phylogram,
            show_names: true,
            show_block_names: true,
            min_block: 2.0,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one taxon's row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(3.0);
        self
    }

    /// Sets the gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Sets how much room the tree gets, in pixels.
    pub fn tree_width(mut self, width: f64) -> Self {
        self.tree_width = width.max(0.0);
        self
    }

    /// Draws the tree as a cladogram or a phylogram.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// Whether to print taxon names between the tree and the blocks.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Whether to print block names on the blocks that are wide enough.
    pub fn show_block_names(mut self, show: bool) -> Self {
        self.show_block_names = show;
        self
    }

    /// Sets the narrowest a block may be drawn, in pixels.
    pub fn min_block(mut self, pixels: f64) -> Self {
        self.min_block = pixels.max(0.5);
        self
    }

    /// The blocks, in the order they were given.
    pub fn blocks(&self) -> &[CladeBlock] {
        &self.blocks
    }

    /// The tree the blocks are painted onto.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Whether block `index` is carried by every leaf under its ancestor.
    ///
    /// False means the carriers are paraphyletic, so the block is either two
    /// events, or one event followed by a loss, or a mistake. The figure says so
    /// by cutting the rows that do not carry it out of the rectangle.
    pub fn is_clade(&self, index: usize) -> bool {
        self.placed
            .get(index)
            .and_then(|p| p.as_ref())
            .map(|p| p.is_clade())
            .unwrap_or(false)
    }

    /// Rows inside block `index` that do not carry it.
    pub fn cut_rows(&self, index: usize) -> usize {
        self.placed
            .get(index)
            .and_then(|p| p.as_ref())
            .map(|p| p.rows() - p.carriers)
            .unwrap_or(0)
    }

    /// Taxa named by a block that the tree does not have.
    ///
    /// A name that matched nothing would otherwise shrink a block silently,
    /// which is the one failure a reader could never spot.
    pub fn unmatched(&self) -> usize {
        self.placed
            .iter()
            .map(|placed| placed.map(|p| p.missing).unwrap_or(0))
            .sum()
    }

    /// Blocks whose taxa were all unknown, so nothing could be drawn for them.
    pub fn unplaced(&self) -> usize {
        self.placed.iter().filter(|placed| placed.is_none()).count()
    }

    fn row_centre(&self, band: Rect, row: usize) -> f64 {
        band.y + row as f64 * (self.row_height + self.row_gap) + self.row_height / 2.0
    }
}

impl Track for CladeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.tree.leaf_count().max(1) as f64;
        rows * self.row_height + (rows - 1.0) * self.row_gap + 12.0
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let names = if self.show_names {
            let size = (theme.font_size - 1.0).min(self.row_height);
            self.tree
                .leaf_names()
                .iter()
                .map(|name| text_width(name, size))
                .fold(0.0f64, f64::max)
                + 8.0
        } else {
            0.0
        };
        self.tree_width + names
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let theme = ctx.theme;
        let rows = self.tree.leaf_count();
        if rows == 0 {
            return;
        }
        let pitch = self.row_height + self.row_gap;
        let name_size = (theme.font_size - 1.0).min(self.row_height);

        if self.tree_width > 0.0 {
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: band.y,
                w: (self.tree_width - 6.0).max(1.0),
                h: band.h,
            };
            draw_tree(
                ctx.svg,
                &self.tree,
                area,
                pitch,
                band.y + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &theme.foreground,
                    width: 1.1,
                    mirror: false,
                },
            );
        }

        if self.show_names {
            for (row, name) in self.tree.leaf_names().iter().enumerate() {
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    self.row_centre(band, row) + name_size * 0.35,
                    name,
                    &theme.muted,
                    name_size,
                    Anchor::End,
                );
            }
        }

        // A hairline per row, so a taxon carrying nothing is still a taxon
        // rather than a gap in the figure.
        let wire = mix(theme.surface(), &theme.rule, 0.5);
        for row in 0..rows {
            let y = self.row_centre(band, row);
            ctx.svg.line(band.x, y, band.right(), y, &wire, 0.6);
        }

        for (slot, index) in self.order.iter().enumerate() {
            let Some(placed) = self.placed[*index] else {
                continue;
            };
            let block = &self.blocks[*index];
            let color = block
                .color
                .clone()
                .unwrap_or_else(|| theme.color(slot).to_string());
            self.draw_block(ctx, block, placed, &color);
        }

        let mut notes: Vec<String> = Vec::new();
        if self.unmatched() > 0 {
            notes.push(format!("{} taxa not in the tree", self.unmatched()));
        }
        if self.unplaced() > 0 {
            notes.push(format!("{} blocks unplaced", self.unplaced()));
        }
        if !notes.is_empty() {
            ctx.svg.text(
                band.right() - 2.0,
                band.bottom() - 2.0,
                &notes.join(", "),
                &theme.muted,
                theme.font_size,
                Anchor::End,
            );
        }
    }
}

impl CladeTrack {
    fn draw_block(
        &self,
        ctx: &mut DrawContext<'_>,
        block: &CladeBlock,
        placed: Placed,
        color: &str,
    ) {
        let band = ctx.band;
        let theme = ctx.theme;
        let x0 = ctx.scale.x(block.start);
        let x1 = ctx.scale.x(block.end);
        let (x0, x1) = if x1 - x0 < self.min_block {
            let mid = (x0 + x1) / 2.0;
            (mid - self.min_block / 2.0, mid + self.min_block / 2.0)
        } else {
            (x0, x1)
        };
        if x1 < band.x - 2.0 || x0 > band.right() + 2.0 {
            return;
        }

        let top = self.row_centre(band, placed.first) - self.row_height / 2.0;
        let bottom = self.row_centre(band, placed.last) + self.row_height / 2.0;
        let fill = wash(color, theme);

        ctx.svg.rect(x0, top, x1 - x0, bottom - top, &fill);

        // The rows that do not carry it are cut out rather than covered, which
        // is the only thing keeping the rectangle honest. Consecutive cut rows
        // become one hole: cutting them one at a time leaves the row gaps
        // between them filled, and four missing taxa then read as stripes rather
        // than as an absence.
        let mut row = placed.first;
        while row <= placed.last {
            if self.carries(block, row) {
                row += 1;
                continue;
            }
            let first = row;
            while row <= placed.last && !self.carries(block, row) {
                row += 1;
            }
            let last = row - 1;
            let reach = self.row_height / 2.0 + self.row_gap / 2.0;
            let y = self.row_centre(band, first) - reach;
            let h = self.row_centre(band, last) + reach - y;
            ctx.svg.rect(x0, y, x1 - x0, h, theme.surface());
            ctx.svg.line(x0, y, x1, y, color, 0.7);
            ctx.svg.line(x0, y + h, x1, y + h, color, 0.7);
        }

        ctx.svg
            .rect_outline(x0, top, x1 - x0, bottom - top, color, 1.0);

        if self.show_block_names {
            if let Some(name) = &block.name {
                let size = theme.font_size;
                let width = text_width(name, size);
                let (left, right) = (x0.max(band.x), x1.min(band.right()));
                let mid = ((left + right) / 2.0).clamp(
                    band.x + width / 2.0,
                    (band.right() - width / 2.0).max(band.x + width / 2.0),
                );
                if right - left >= width + 8.0 {
                    ctx.svg.text(
                        mid,
                        (top + bottom) / 2.0 + size * 0.35,
                        name,
                        &theme.foreground,
                        size,
                        Anchor::Middle,
                    );
                } else {
                    // A ten kilobase deletion on a four megabase axis is two
                    // pixels wide, so at genome scale no block ever holds its
                    // own name. Beside it the name still points at it, and the
                    // alternative is a figure with every label missing.
                    let above = top - 2.0;
                    let below = bottom + size;
                    let y = if above - size >= band.y {
                        Some(above)
                    } else if below <= band.bottom() {
                        Some(below)
                    } else {
                        None
                    };
                    if let Some(y) = y {
                        ctx.svg
                            .text(mid, y, name, &theme.muted, size, Anchor::Middle);
                    }
                }
            }
        }
    }

    /// Whether the leaf on `row` is one of the block's taxa.
    fn carries(&self, block: &CladeBlock, row: usize) -> bool {
        self.tree
            .leaf_names()
            .get(row)
            .map(|name| block.taxa.iter().any(|taxon| taxon == name))
            .unwrap_or(false)
    }
}

/// Row of every leaf, by name, in the order the tree lays them out.
fn row_of_leaf(tree: &Tree) -> Vec<String> {
    tree.leaf_names()
}

/// The row span a block covers, or `None` when the tree names none of its taxa.
fn place(block: &CladeBlock, rows: &[String], leaves: usize) -> Option<Placed> {
    let mut first = usize::MAX;
    let mut last = 0usize;
    let mut carriers = 0usize;
    let mut missing = 0usize;

    for taxon in &block.taxa {
        match rows.iter().position(|name| name == taxon) {
            Some(row) if row < leaves => {
                first = first.min(row);
                last = last.max(row);
                carriers += 1;
            }
            _ => missing += 1,
        }
    }
    if carriers == 0 {
        return None;
    }
    Some(Placed {
        first,
        last,
        carriers,
        missing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn tree() -> Tree {
        // ((A,B),(C,D)),(E,F) with A and B sisters, C and D sisters.
        Tree::parse_newick("((((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):2):1,G:4);")
            .expect("the tree in this test is well formed")
    }

    fn scale() -> Scale {
        Scale::new(&Region::new("chr", 0, 10_000).unwrap(), 0.0, 900.0)
    }

    #[test]
    fn sisters_make_a_clade_and_a_scattered_pair_does_not() {
        let track = CladeTrack::new(
            tree(),
            vec![
                CladeBlock::new(0, 100, ["A", "B"]),
                CladeBlock::new(0, 100, ["A", "D"]),
                CladeBlock::new(0, 100, ["A", "B", "C", "D"]),
            ],
        );
        assert!(track.is_clade(0));
        assert!(!track.is_clade(1), "A and D have B and C between them");
        assert!(track.is_clade(2));
    }

    #[test]
    fn the_rows_a_block_does_not_carry_are_counted() {
        let track = CladeTrack::new(tree(), vec![CladeBlock::new(0, 100, ["A", "D"])]);
        assert_eq!(track.cut_rows(0), 2, "B and C sit inside the span");
        let clean = CladeTrack::new(tree(), vec![CladeBlock::new(0, 100, ["A", "B"])]);
        assert_eq!(clean.cut_rows(0), 0);
    }

    #[test]
    fn a_taxon_the_tree_does_not_have_is_reported() {
        let track = CladeTrack::new(
            tree(),
            vec![CladeBlock::new(0, 100, ["A", "B", "ERR_not_here"])],
        );
        assert_eq!(track.unmatched(), 1);
        assert_eq!(track.unplaced(), 0);
        assert!(track.is_clade(0), "the two that matched are still sisters");
    }

    #[test]
    fn a_block_nobody_carries_is_not_drawn_and_is_counted() {
        let track = CladeTrack::new(tree(), vec![CladeBlock::new(0, 100, ["nobody"])]);
        assert_eq!(track.unplaced(), 1);
        assert!(!track.is_clade(0));
        let svg = Figure::new(Region::new("chr", 0, 10_000).unwrap())
            .push(track)
            .to_svg();
        assert!(svg.contains("1 blocks unplaced"));
    }

    #[test]
    fn one_carrier_is_a_block_one_row_tall() {
        let track = CladeTrack::new(tree(), vec![CladeBlock::new(0, 100, ["E"])]);
        assert!(track.is_clade(0));
        assert_eq!(track.cut_rows(0), 0);
    }

    #[test]
    fn the_widest_clade_is_drawn_first_so_a_tip_block_sits_on_top() {
        let track = CladeTrack::new(
            tree(),
            vec![
                CladeBlock::new(0, 100, ["A"]).name("tip"),
                CladeBlock::new(0, 100, ["A", "B", "C", "D"]).name("deep"),
            ],
        );
        // Order is by rows covered, so the four-taxon block paints first.
        assert_eq!(track.order, vec![1, 0]);
    }

    #[test]
    fn the_cut_out_rows_actually_reach_the_svg() {
        let scattered = CladeTrack::new(tree(), vec![CladeBlock::new(2_000, 8_000, ["A", "D"])]);
        let solid = CladeTrack::new(
            tree(),
            vec![CladeBlock::new(2_000, 8_000, ["A", "B", "C", "D"])],
        );
        let render = |track: CladeTrack| {
            Figure::new(Region::new("chr", 0, 10_000).unwrap())
                .push(track)
                .to_svg()
        };
        let (a, b) = (render(scattered), render(solid));
        assert_ne!(a, b, "a paraphyletic block cannot look like a clade");
        assert!(
            a.matches("<rect").count() > b.matches("<rect").count(),
            "the cut rows are extra rectangles"
        );
    }

    #[test]
    fn a_block_is_placed_on_true_coordinates() {
        let track = CladeTrack::new(tree(), vec![CladeBlock::new(2_000, 3_000, ["A", "B"])]);
        let svg = Figure::new(Region::new("chr", 0, 10_000).unwrap())
            .width(1_000.0)
            .push(track.tree_width(0.0).show_names(false))
            .to_svg();
        // A fifth of the way across a 10 kb view, a tenth of it wide. If the
        // block were drawn on an ordinal axis this would be nowhere near.
        assert!(svg.contains("<rect"));
        let track = CladeTrack::new(tree(), vec![CladeBlock::new(2_000, 3_000, ["A", "B"])]);
        assert_eq!(track.blocks()[0].span(), 1_000);
    }

    #[test]
    fn a_narrow_block_is_still_visible() {
        let track = CladeTrack::new(
            tree(),
            vec![CladeBlock::new(1_000_000, 1_000_050, ["A", "B"])],
        )
        .min_block(4.0);
        let svg = Figure::new(Region::new("chr", 0, 4_400_000).unwrap())
            .width(900.0)
            .push(track)
            .to_svg();
        assert!(svg.contains("<rect"));
    }

    #[test]
    fn a_block_off_screen_draws_nothing_of_itself() {
        let render = |blocks: Vec<CladeBlock>| {
            Figure::new(Region::new("chr", 0, 10_000).unwrap())
                .push(
                    CladeTrack::new(tree(), blocks)
                        .tree_width(0.0)
                        .show_names(false),
                )
                .to_svg()
        };
        // The figure's own page rectangle is always there, so the test is
        // whether the block added any of its own.
        let empty = render(vec![]).matches("<rect").count();
        let away = render(vec![CladeBlock::new(900_000, 950_000, ["A", "B"])]);
        let here = render(vec![CladeBlock::new(2_000, 3_000, ["A", "B"])]);
        assert_eq!(away.matches("<rect").count(), empty, "nothing is in view");
        assert!(here.matches("<rect").count() > empty);
    }

    #[test]
    fn height_follows_the_number_of_leaves() {
        let seven = CladeTrack::new(tree(), vec![]);
        let two = CladeTrack::new(Tree::parse_newick("(A:1,B:1);").unwrap(), vec![]);
        assert!(seven.height(&scale()) > two.height(&scale()));
    }

    #[test]
    fn the_names_strip_grows_with_the_longest_name() {
        let theme = Theme::light();
        let short = CladeTrack::new(tree(), vec![]);
        let long = CladeTrack::new(
            Tree::parse_newick("(a_very_long_isolate_name_indeed:1,B:1);").unwrap(),
            vec![],
        );
        assert!(long.y_axis_width(&theme) > short.y_axis_width(&theme));
    }
}
