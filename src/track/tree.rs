//! A phylogeny drawn beside its rows.
//!
//! Like the ideogram, this track does not use the shared horizontal scale: its
//! x axis is evolutionary distance, not genomic position, and the two have
//! nothing to do with each other. Unlike the ideogram, its y axis means
//! something to its neighbours, because a leaf is a row.
//!
//! That is the whole point. [`SnpTrack::tree`](crate::SnpTrack::tree) puts one
//! of these in the strip beside a panel of variable sites and sorts the rows to
//! match, so a clade's shared substitutions line up into a block instead of
//! being scattered down the panel in whatever order the samples were listed.

use crate::scale::Scale;
use crate::svg::text_width;
use crate::theme::Theme;
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// How to draw the horizontal extent of a tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeShape {
    /// Branch lengths carry distance, so tip positions mean something.
    Phylogram,
    /// Branches are counted rather than measured, and every tip lines up on
    /// the right. Use it when the lengths are missing or not to be trusted.
    Cladogram,
}

/// The order `names` should be in for its rows to line up with `tree`.
///
/// Returns a permutation: index `i` of the result is the row that belongs on
/// line `i`. Matching is by name, and a row the tree never names keeps its
/// place at the end rather than being dropped, because a row silently missing
/// from a figure is worse than a row out of order. A duplicate name is matched
/// once, so two rows called the same thing both survive.
pub fn leaf_order(tree: &Tree, names: &[String]) -> Vec<usize> {
    let mut taken = vec![false; names.len()];
    let mut order: Vec<usize> = Vec::with_capacity(names.len());

    for leaf in tree.leaf_names() {
        if let Some(index) = names
            .iter()
            .enumerate()
            .position(|(index, name)| *name == leaf && !taken[index])
        {
            taken[index] = true;
            order.push(index);
        }
    }
    for (index, used) in taken.iter().enumerate() {
        if !used {
            order.push(index);
        }
    }
    order
}

/// How the branches of a tree are drawn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TreeStyle<'a> {
    /// Phylogram or cladogram.
    pub shape: TreeShape,
    /// Branch colour.
    pub color: &'a str,
    /// Branch width in pixels.
    pub width: f64,
    /// Whether the root is on the right and the tips on the left.
    ///
    /// For the second tree of a tanglegram, which faces the first.
    pub mirror: bool,
}

/// Draws a tree into a rectangle, rows evenly spaced down it.
///
/// Shared by the standalone [`TreeTrack`] and by the row tracks that draw a
/// tree in their own strip, so both put a branch in exactly the same place.
pub fn draw_tree(
    svg: &mut crate::svg::SvgWriter,
    tree: &Tree,
    area: Rect,
    row_pitch: f64,
    first_row_centre: f64,
    style: TreeStyle<'_>,
) {
    let (shape, color, width) = (style.shape, style.color, style.width);
    let style = &style;
    let cladogram = shape == TreeShape::Cladogram;
    let layout = tree.layout(cladogram);
    let span = tree.max_depth(cladogram);
    let x_of = |depth: f64| {
        let fraction = if span <= 0.0 { 0.0 } else { depth / span };
        if style.mirror {
            area.right() - fraction * area.w
        } else {
            area.x + fraction * area.w
        }
    };
    let y_of = |row: f64| first_row_centre + row * row_pitch;

    for placement in &layout {
        let node = &tree.nodes()[placement.node];
        let Some(parent) = node.parent else {
            continue;
        };
        let parent_placement = layout[parent];
        let x0 = x_of(parent_placement.depth);
        let x1 = x_of(placement.depth);
        let y = y_of(placement.row);

        // Rectangular branches: along to the child's depth, then the parent's
        // riser joins its children. Diagonals would imply the tree says
        // something about the space between two rows, and it does not.
        svg.line(x0, y, x1, y, color, width);
    }

    for placement in &layout {
        let node = &tree.nodes()[placement.node];
        if node.is_leaf() || node.children.is_empty() {
            continue;
        }
        let rows: Vec<f64> = node
            .children
            .iter()
            .map(|child| layout[*child].row)
            .collect();
        let (top, bottom) = rows.iter().fold((f64::MAX, f64::MIN), |(lo, hi), row| {
            (lo.min(*row), hi.max(*row))
        });
        let x = x_of(placement.depth);
        svg.line(x, y_of(top), x, y_of(bottom), color, width);
    }
}

/// A phylogeny as a track of its own.
///
/// ```
/// use karyon::{Figure, Region, TreeTrack};
/// use karyon::tree::Tree;
///
/// let tree = Tree::parse_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
/// let svg = Figure::new(Region::new("tree", 0, 1).unwrap())
///     .push(TreeTrack::new(tree).label("phylogeny"))
///     .to_svg();
/// assert!(svg.contains("<line"));
/// ```
#[derive(Debug, Clone)]
pub struct TreeTrack {
    tree: Tree,
    label: Option<String>,
    row_height: f64,
    shape: TreeShape,
    color: Option<String>,
    line_width: f64,
    show_tips: bool,
}

impl TreeTrack {
    /// A track drawing `tree`.
    pub fn new(tree: Tree) -> Self {
        TreeTrack {
            tree,
            label: None,
            row_height: 15.0,
            shape: TreeShape::Phylogram,
            color: None,
            line_width: 1.2,
            show_tips: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the vertical pitch of one leaf.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(2.0);
        self
    }

    /// Chooses a phylogram or a cladogram.
    pub fn shape(mut self, shape: TreeShape) -> Self {
        self.shape = shape;
        self
    }

    /// Sets the branch colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets the branch width.
    pub fn line_width(mut self, width: f64) -> Self {
        self.line_width = width.max(0.2);
        self
    }

    /// Draws or hides the tip names.
    pub fn show_tips(mut self, show: bool) -> Self {
        self.show_tips = show;
        self
    }

    /// The tree.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Width the tip names need.
    fn tip_width(&self, theme: &Theme) -> f64 {
        if !self.show_tips {
            return 0.0;
        }
        self.tree
            .leaf_names()
            .iter()
            .map(|name| text_width(name, theme.font_size - 1.0))
            .fold(0.0f64, f64::max)
            + 6.0
    }
}

impl Track for TreeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        (self.tree.leaf_count().max(1) as f64) * self.row_height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let color = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.foreground.clone());
        let tips = self.tip_width(ctx.theme);
        let area = Rect {
            x: band.x,
            y: band.y,
            w: (band.w - tips).max(1.0),
            h: band.h,
        };

        draw_tree(
            ctx.svg,
            &self.tree,
            area,
            self.row_height,
            band.y + self.row_height / 2.0,
            TreeStyle {
                shape: self.shape,
                color: &color,
                width: self.line_width,
                mirror: false,
            },
        );

        if self.show_tips {
            let size = ctx.theme.font_size - 1.0;
            for (row, name) in self.tree.leaf_names().iter().enumerate() {
                ctx.svg.text(
                    area.right() + 4.0,
                    band.y + self.row_height / 2.0 + row as f64 * self.row_height + size * 0.35,
                    name,
                    &ctx.theme.muted,
                    size,
                    crate::svg::Anchor::Start,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn tree() -> Tree {
        Tree::parse_newick("((A:0.1,B:0.2)0.9:0.3,(C:0.15,D:0.05):0.2);").unwrap()
    }

    fn region() -> Region {
        Region::new("tree", 0, 1).unwrap()
    }

    #[test]
    fn height_follows_the_leaf_count() {
        let scale = Scale::new(&region(), 0.0, 100.0);
        assert_eq!(TreeTrack::new(tree()).height(&scale), 4.0 * 15.0);
        assert_eq!(
            TreeTrack::new(tree()).row_height(20.0).height(&scale),
            4.0 * 20.0
        );
    }

    #[test]
    fn every_branch_and_every_riser_is_drawn() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()).show_tips(false))
            .to_svg();
        // Six branches, since the root has none, and three risers: one for
        // each pair of tips and one joining those two pairs at the root.
        assert_eq!(svg.matches("<line").count(), 9);
    }

    #[test]
    fn tip_names_are_drawn_when_asked_for() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()))
            .to_svg();
        for tip in ["A", "B", "C", "D"] {
            assert!(svg.contains(&format!(">{tip}</text>")), "missing {tip}");
        }
    }

    #[test]
    fn a_cladogram_lines_every_tip_up_and_a_phylogram_does_not() {
        let x_of_tips = |shape: TreeShape| {
            let svg = Figure::new(region())
                .show_region_label(false)
                .push(TreeTrack::new(tree()).shape(shape).show_tips(false))
                .to_svg();
            let mut ends: Vec<String> = svg
                .match_indices(r#"x2=""#)
                .map(|(index, prefix)| {
                    let rest = &svg[index + prefix.len()..];
                    rest[..rest.find('"').unwrap()].to_string()
                })
                .collect();
            ends.sort();
            ends.dedup();
            ends
        };
        // A cladogram has fewer distinct branch ends, because the tips share one.
        assert!(x_of_tips(TreeShape::Cladogram).len() < x_of_tips(TreeShape::Phylogram).len());
    }

    #[test]
    fn a_tree_with_no_lengths_still_draws_as_a_phylogram() {
        let flat = Tree::parse_newick("((A,B),(C,D));").unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(flat).show_tips(false))
            .to_svg();
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn a_single_leaf_draws_without_dividing_by_zero() {
        let one = Tree::parse_newick("A;").unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(one))
            .to_svg();
        assert!(!svg.contains("NaN"));
        assert!(svg.contains(">A</text>"));
    }
}
