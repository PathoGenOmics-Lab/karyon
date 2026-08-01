//! Two trees face to face, with the same tips joined.
//!
//! # What it is for
//!
//! Two trees over the same taxa that disagree, and the question of where. A
//! gene tree against a species tree, a core tree against an accessory tree,
//! two methods over one alignment. Drawn side by side the disagreement is
//! something you have to hold in your head; drawn facing each other with each
//! tip joined to its twin, the disagreement is the crossings, and a crossing is
//! a thing you can point at.
//!
//! # What a crossing means
//!
//! That two taxa are in one order on the left and the other order on the right.
//! [`TanglegramTrack::crossings`] counts them, which is a summary of how much
//! the two trees disagree and is worth putting in the caption.
//!
//! It is not, however, a statistic. The count depends on how each tree happened
//! to rotate its clades, and a clade can be rotated freely without changing
//! what the tree says. Two trees that agree completely can be drawn with a
//! great many crossings by an unlucky rotation. Untangling is a separate
//! problem this does not solve; the count is what the drawing shows, and the
//! drawing is one of many.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::tree::{draw_tree, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// Two trees drawn facing each other, tips joined.
///
/// ```
/// use karyon::tree::Tree;
/// use karyon::{Figure, Region, TanglegramTrack};
///
/// let genes = Tree::parse_newick("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();
/// let species = Tree::parse_newick("((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);").unwrap();
///
/// let track = TanglegramTrack::new(genes, species);
/// assert!(track.crossings() > 0, "the two disagree");
///
/// let svg = Figure::new(Region::new("taxa", 0, 4).unwrap())
///     .push(track.names("gene tree", "species tree"))
///     .to_svg();
/// assert!(svg.contains("gene tree"));
/// ```
#[derive(Debug, Clone)]
pub struct TanglegramTrack {
    left: Tree,
    right: Tree,
    label: Option<String>,
    left_name: Option<String>,
    right_name: Option<String>,
    row_height: f64,
    tree_width: f64,
    shape: TreeShape,
    color: Option<String>,
    tie_color: Option<String>,
    crossing_color: Option<String>,
    show_tips: bool,
}

impl TanglegramTrack {
    /// A pair of trees, the first on the left.
    pub fn new(left: Tree, right: Tree) -> Self {
        TanglegramTrack {
            left,
            right,
            label: None,
            left_name: None,
            right_name: None,
            row_height: 16.0,
            tree_width: 0.3,
            shape: TreeShape::Phylogram,
            color: None,
            tie_color: None,
            crossing_color: None,
            show_tips: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Names the two trees, drawn over them.
    pub fn names(mut self, left: impl Into<String>, right: impl Into<String>) -> Self {
        self.left_name = Some(left.into());
        self.right_name = Some(right.into());
        self
    }

    /// Sets the vertical pitch of one tip.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(4.0);
        self
    }

    /// Sets how much of the width each tree gets, as a fraction.
    ///
    /// The rest is the middle, where the ties run. A third each leaves a third
    /// for them, which is enough for a crossing to be a crossing rather than a
    /// kink.
    pub fn tree_width(mut self, fraction: f64) -> Self {
        self.tree_width = fraction.clamp(0.05, 0.45);
        self
    }

    /// Chooses a phylogram or a cladogram for both trees.
    pub fn shape(mut self, shape: TreeShape) -> Self {
        self.shape = shape;
        self
    }

    /// Sets the branch colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets the colours of a tie that runs straight and one that crosses.
    ///
    /// Two colours because the crossings are the finding, and a tanglegram in
    /// one colour makes the reader trace every line to find them.
    pub fn tie_colors(mut self, straight: impl Into<String>, crossing: impl Into<String>) -> Self {
        self.tie_color = Some(straight.into());
        self.crossing_color = Some(crossing.into());
        self
    }

    /// Draws or hides the tip names down the middle.
    pub fn show_tips(mut self, show: bool) -> Self {
        self.show_tips = show;
        self
    }

    /// The two trees.
    pub fn trees(&self) -> (&Tree, &Tree) {
        (&self.left, &self.right)
    }

    /// The tips both trees have, in the order the left one draws them.
    ///
    /// A tip in one tree and not the other has nothing to be joined to and is
    /// left out of the ties, though it is still drawn on the tree that has it:
    /// a taxon missing from one analysis is a fact about the analysis.
    pub fn shared(&self) -> Vec<String> {
        let right = self.right.leaf_names();
        self.left
            .leaf_names()
            .into_iter()
            .filter(|name| !name.is_empty() && right.contains(name))
            .collect()
    }

    /// Tips the two trees do not have in common.
    pub fn unshared(&self) -> Vec<String> {
        let shared = self.shared();
        let mut odd: Vec<String> = self
            .left
            .leaf_names()
            .into_iter()
            .chain(self.right.leaf_names())
            .filter(|name| !name.is_empty() && !shared.contains(name))
            .collect();
        odd.sort();
        odd.dedup();
        odd
    }

    /// Each shared tip as its row in the left tree and its row in the right.
    pub fn ties(&self) -> Vec<(String, usize, usize)> {
        let left = self.left.leaf_names();
        let right = self.right.leaf_names();
        self.shared()
            .into_iter()
            .filter_map(|name| {
                let from = left.iter().position(|leaf| *leaf == name)?;
                let to = right.iter().position(|leaf| *leaf == name)?;
                Some((name, from, to))
            })
            .collect()
    }

    /// How many pairs of ties cross.
    ///
    /// The count depends on how each tree happened to rotate its clades, and a
    /// clade rotates freely without changing what the tree says. It is what
    /// this drawing shows, not a property of the two trees.
    pub fn crossings(&self) -> usize {
        let ties = self.ties();
        let mut count = 0usize;
        for (index, (_, from, to)) in ties.iter().enumerate() {
            for (other_from, other_to) in ties[index + 1..].iter().map(|(_, a, b)| (a, b)) {
                if (from < other_from) != (to < other_to) {
                    count += 1;
                }
            }
        }
        count
    }

    /// Whether one tie crosses any other.
    fn crosses(&self, ties: &[(String, usize, usize)], index: usize) -> bool {
        let (_, from, to) = &ties[index];
        ties.iter()
            .enumerate()
            .any(|(other, (_, a, b))| other != index && ((from < a) != (to < b)))
    }

    /// How many rows the taller of the two trees has.
    fn rows(&self) -> usize {
        self.left.leaf_count().max(self.right.leaf_count()).max(1)
    }

    /// Room above the trees for their names.
    ///
    /// The theme reaches the height only through the font size, and the default
    /// is the one the figure will use unless the caller changed it, as in the
    /// feature track.
    fn header(&self, theme: &Theme) -> f64 {
        if self.left_name.is_some() {
            theme.font_size - 2.0 + 4.0
        } else {
            0.0
        }
    }
}

impl Track for TanglegramTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        // The names sit above the trees, so they are part of the height. Left
        // out, the last tip of an eight tip pair fell past the bottom of the
        // band and the clip took it.
        self.rows() as f64 * self.row_height + self.header(&Theme::default())
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
        let straight = self
            .tie_color
            .clone()
            .unwrap_or_else(|| mix(ctx.theme.surface(), &ctx.theme.rule, 0.9));
        let crossing = self
            .crossing_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());

        let side = band.w * self.tree_width;
        let name_size = ctx.theme.font_size - 2.0;
        let header = if self.left_name.is_some() {
            name_size + 4.0
        } else {
            0.0
        };
        let top = band.y + header;
        let tip_first = top + self.row_height / 2.0;

        // Ties first, under both trees.
        let ties = self.ties();
        for (index, (name, from, to)) in ties.iter().enumerate() {
            let y0 = tip_first + *from as f64 * self.row_height;
            let y1 = tip_first + *to as f64 * self.row_height;
            let crosses = self.crosses(&ties, index);
            let ink = if crosses { &crossing } else { &straight };
            let x0 = band.x + side;
            let x1 = band.right() - side;
            // Leaving each tip horizontally so the line reads as belonging to
            // that tip rather than pointing at whatever is above it.
            let d = format!(
                "M{} {}C{} {} {} {} {} {}",
                num(x0),
                num(y0),
                num(x0 + (x1 - x0) * 0.35),
                num(y0),
                num(x1 - (x1 - x0) * 0.35),
                num(y1),
                num(x1),
                num(y1)
            );
            ctx.svg
                .path_stroked(&d, ink, if crosses { 1.3 } else { 0.9 });

            if self.show_tips {
                // At the left tip's row, not the middle of the tie. On the
                // midpoint two crossing ties put their names in the same place
                // and neither can be read.
                ctx.svg.text(
                    x0 + (x1 - x0) * 0.22,
                    y0 + name_size * 0.35,
                    name,
                    &ctx.theme.muted,
                    name_size,
                    Anchor::Start,
                );
            }
        }

        for (tree, mirror) in [(&self.left, false), (&self.right, true)] {
            let area = Rect {
                x: if mirror { band.right() - side } else { band.x },
                y: top,
                w: side,
                h: band.h - header,
            };
            draw_tree(
                ctx.svg,
                tree,
                area,
                self.row_height,
                tip_first,
                TreeStyle {
                    shape: self.shape,
                    color: &color,
                    width: 1.2,
                    mirror,
                },
            );
        }

        if let (Some(left), Some(right)) = (&self.left_name, &self.right_name) {
            ctx.svg.text(
                band.x,
                band.y + name_size,
                left,
                &ctx.theme.muted,
                name_size,
                Anchor::Start,
            );
            ctx.svg.text(
                band.right(),
                band.y + name_size,
                right,
                &ctx.theme.muted,
                name_size,
                Anchor::End,
            );
        }

        let odd = self.unshared().len();
        if odd > 0 {
            // A tip on one tree and not the other has no tie, and a missing
            // line is exactly the thing a reader will not notice.
            let text = format!("{odd} not in both");
            ctx.svg.text(
                band.right() - text_width(&text, name_size - 1.0) / 2.0 - 3.0,
                band.bottom() - 2.0,
                &text,
                &ctx.theme.muted,
                name_size - 1.0,
                Anchor::Middle,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn agreeing() -> TanglegramTrack {
        let one = Tree::parse_newick("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();
        let two = Tree::parse_newick("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();
        TanglegramTrack::new(one, two)
    }

    fn disagreeing() -> TanglegramTrack {
        let genes = Tree::parse_newick("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);").unwrap();
        let species = Tree::parse_newick("((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);").unwrap();
        TanglegramTrack::new(genes, species)
    }

    fn region() -> Region {
        Region::new("taxa", 0, 4).unwrap()
    }

    #[test]
    fn two_trees_in_the_same_order_have_no_crossings() {
        assert_eq!(agreeing().crossings(), 0);
        assert_eq!(agreeing().shared().len(), 4);
    }

    #[test]
    fn a_disagreement_is_a_crossing() {
        // A B C D against A C B D: B and C swap, which is one crossing.
        let track = disagreeing();
        assert_eq!(track.crossings(), 1);
        assert_eq!(
            track.ties(),
            vec![
                ("A".to_string(), 0, 0),
                ("B".to_string(), 1, 2),
                ("C".to_string(), 2, 1),
                ("D".to_string(), 3, 3),
            ]
        );
    }

    #[test]
    fn a_reversed_tree_crosses_everything() {
        let forward = Tree::parse_newick("(((A,B),C),D);").unwrap();
        let backward = Tree::parse_newick("(D,(C,(B,A)));").unwrap();
        let track = TanglegramTrack::new(forward, backward);
        // Four tips reversed is every pair of them, which is six.
        assert_eq!(track.crossings(), 6);
    }

    #[test]
    fn a_tip_in_only_one_tree_has_nothing_to_join() {
        let one = Tree::parse_newick("((A,B),(C,D));").unwrap();
        let two = Tree::parse_newick("((A,B),(C,E));").unwrap();
        let track = TanglegramTrack::new(one, two);
        assert_eq!(track.shared(), vec!["A", "B", "C"]);
        assert_eq!(track.unshared(), vec!["D", "E"]);
        assert_eq!(track.ties().len(), 3);
    }

    #[test]
    fn the_odd_tips_are_reported_rather_than_left_silent() {
        let one = Tree::parse_newick("((A,B),(C,D));").unwrap();
        let two = Tree::parse_newick("((A,B),(C,E));").unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TanglegramTrack::new(one, two))
            .to_svg();
        assert!(svg.contains("2 not in both"), "a missing line goes unseen");
    }

    #[test]
    fn a_crossing_tie_is_drawn_differently_from_a_straight_one() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(disagreeing().tie_colors("#111111", "#222222"))
            .to_svg();
        assert!(svg.contains("#111111"), "no straight ties");
        assert!(svg.contains("#222222"), "no crossing ties");

        let calm = Figure::new(region())
            .show_region_label(false)
            .push(agreeing().tie_colors("#111111", "#222222"))
            .to_svg();
        assert!(!calm.contains("#222222"), "nothing crosses here");
    }

    #[test]
    fn both_trees_are_drawn_and_the_second_faces_the_first() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(disagreeing().show_tips(false))
            .to_svg();
        // Two trees of four tips: six branches and three risers each.
        assert_eq!(svg.matches("<line").count(), 18);
        // The mirrored one reaches the right edge; the other starts at the left.
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn height_follows_the_taller_of_the_two() {
        let scale = Scale::new(&region(), 0.0, 400.0);
        let lopsided = TanglegramTrack::new(
            Tree::parse_newick("((A,B),(C,D));").unwrap(),
            Tree::parse_newick("(A,B);").unwrap(),
        );
        assert_eq!(lopsided.height(&scale), 4.0 * 16.0);
        // Naming the trees puts a line above them, and that line is part of
        // the height: left out, the last tip fell past the band.
        let named = TanglegramTrack::new(
            Tree::parse_newick("((A,B),(C,D));").unwrap(),
            Tree::parse_newick("((A,B),(C,D));").unwrap(),
        )
        .names("one", "two");
        assert!(named.height(&scale) > 4.0 * 16.0);
    }

    #[test]
    fn the_trees_can_be_named() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(disagreeing().names("gene tree", "species tree"))
            .to_svg();
        assert!(svg.contains(">gene tree</text>"));
        assert!(svg.contains(">species tree</text>"));
    }

    #[test]
    fn two_trees_with_nothing_in_common_draw_without_panicking() {
        let one = Tree::parse_newick("((A,B),(C,D));").unwrap();
        let two = Tree::parse_newick("((W,X),(Y,Z));").unwrap();
        let track = TanglegramTrack::new(one, two);
        assert!(track.shared().is_empty());
        assert_eq!(track.crossings(), 0);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(!svg.contains("NaN"));
        assert!(svg.contains("8 not in both"));
    }
}
