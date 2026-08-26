//! Two trees face to face, with the same tips joined.
//!
//! [`TanglegramTrack`] takes two trees over the same taxa, draws the first
//! facing the second, and runs a tie from each tip to the tip of the same name
//! opposite. The ties are the figure; the trees are there to put the tips in
//! the order each of them argues for.
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
//! [`TanglegramTrack::crossings`] counts them, which is worth putting in a
//! caption and is not a statistic: the count belongs to this drawing rather
//! than to the two trees, and the drawing is one of many.
//! [`TanglegramTrack::untangle`] greedily rotates free clades on both sides and
//! keeps only strict improvements. It never changes topology or makes the
//! drawing worse, but it is a deterministic local heuristic rather than a
//! claim to the global minimum.
//!
//! # Neither axis is the figure's
//!
//! Like [`TreeTrack`](crate::TreeTrack), this track never reads the shared
//! scale. Each tree draws its own axis of evolutionary distance into the share
//! of the band [`TanglegramTrack::tree_width`] gives it, mirrored on the right.
//! The [`Region`](crate::Region) a figure is built on is still required and
//! still ignored: nothing here is measured in bases.
//!
//! The vertical axis is not shared either, and that is the subject rather than
//! an oversight. Row three on the left and row three on the right are two
//! different taxa wherever the trees disagree, so there is no one row order a
//! neighbouring track could be sorted into. A tanglegram is read on its own.

use crate::scale::Scale;
use std::collections::BTreeMap;

use crate::style::LinePattern;
use crate::svg::{fit_text, num, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::tree::{draw_tree_titled, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// Geometry used to join matching taxa between the two trees.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TangleTieStyle {
    /// A cubic curve leaves and reaches each tip horizontally.
    #[default]
    Curved,
    /// A direct segment, useful for compact or low-crossing comparisons.
    Straight,
    /// A translucent band that remains visible in large exported figures.
    Ribbon,
}

/// Which terminal names are written in the comparison corridor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TangleLabels {
    /// Keep names in branch tooltips only.
    None,
    /// Write one name beside the left tree.
    Left,
    /// Write one name beside the right tree.
    Right,
    /// Write names at both ends of every tie, the clearest tracing view.
    #[default]
    Both,
}

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
    tie_style: TangleTieStyle,
    labels: TangleLabels,
    label_width: f64,
    tie_width: f64,
    crossing_width: f64,
    color_by: Option<String>,
    show_summary: bool,
    initial_crossings: usize,
}

impl TanglegramTrack {
    /// A pair of trees, the first on the left.
    pub fn new(left: Tree, right: Tree) -> Self {
        let initial_crossings = count_crossings(&left, &right);
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
            tie_style: TangleTieStyle::Curved,
            labels: TangleLabels::Both,
            label_width: 72.0,
            tie_width: 0.9,
            crossing_width: 1.5,
            color_by: None,
            show_summary: true,
            initial_crossings,
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

    /// Chooses curved lines, direct segments or translucent ribbons.
    pub fn tie_style(mut self, style: TangleTieStyle) -> Self {
        self.tie_style = style;
        self
    }

    /// Sets the ordinary and crossing tie widths in pixels.
    pub fn tie_widths(mut self, ordinary: f64, crossing: f64) -> Self {
        if ordinary.is_finite() {
            self.tie_width = ordinary.clamp(0.2, 12.0);
        }
        if crossing.is_finite() {
            self.crossing_width = crossing.clamp(self.tie_width, 16.0);
        }
        self
    }

    /// Chooses where matching terminal names are written.
    pub fn labels(mut self, labels: TangleLabels) -> Self {
        self.labels = labels;
        self
    }

    /// Sets the maximum room reserved for a terminal name at each side.
    pub fn label_width(mut self, width: f64) -> Self {
        if width.is_finite() {
            self.label_width = width.clamp(24.0, 180.0);
        }
        self
    }

    /// Colours ties by a terminal annotation shared by the two trees.
    ///
    /// Equal values receive one categorical colour. A disagreement is drawn
    /// with the crossing colour and a dashed centre line, while endpoint marks
    /// retain the value from each tree. Tooltips report both exact values.
    pub fn color_by(mut self, key: impl Into<String>) -> Self {
        self.color_by = Some(key.into());
        self
    }

    /// Draws or hides the compact crossing and matching summary.
    pub fn show_summary(mut self, show: bool) -> Self {
        self.show_summary = show;
        self
    }

    /// Greedily rotates free clades on both sides to reduce crossings.
    ///
    /// Rotation never changes a clade or branch length. The heuristic keeps a
    /// rotation only when it strictly lowers the crossing count, alternates
    /// between the two trees and therefore cannot make the drawing worse.
    pub fn untangle(self) -> Self {
        self.untangle_passes(12)
    }

    /// As [`TanglegramTrack::untangle`], with a bounded number of passes.
    pub fn untangle_passes(mut self, passes: usize) -> Self {
        for _ in 0..passes.clamp(1, 64) {
            let right = improve_rotations(&mut self.right, &self.left, false);
            let left = improve_rotations(&mut self.left, &self.right, true);
            if !right && !left {
                break;
            }
        }
        self
    }

    /// Draws or hides the tip names down the middle.
    pub fn show_tips(mut self, show: bool) -> Self {
        self.labels = if show {
            TangleLabels::Both
        } else {
            TangleLabels::None
        };
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
        count_crossings(&self.left, &self.right)
    }

    /// Crossing count before any call to [`TanglegramTrack::untangle`].
    pub fn initial_crossings(&self) -> usize {
        self.initial_crossings
    }

    /// Number of crossings removed by automatic clade rotation.
    pub fn crossing_reduction(&self) -> usize {
        self.initial_crossings.saturating_sub(self.crossings())
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
        if self.left_name.is_some() || self.right_name.is_some() || self.show_summary {
            theme.font_size - 2.0 + 8.0
        } else {
            0.0
        }
    }

    fn annotation(&self, tree: &Tree, name: &str) -> Option<String> {
        let key = self.color_by.as_deref()?;
        let node = tree.node_named(name)?;
        tree.annotation(node, key)
            .or_else(|| {
                tree.ancestors(node)
                    .into_iter()
                    .find_map(|ancestor| tree.annotation(ancestor, key))
            })
            .map(ToString::to_string)
    }

    fn annotation_categories(&self) -> BTreeMap<String, usize> {
        let mut categories = BTreeMap::new();
        for name in self.shared() {
            for value in [
                self.annotation(&self.left, &name),
                self.annotation(&self.right, &name),
            ]
            .into_iter()
            .flatten()
            {
                categories.insert(value, 0);
            }
        }
        for (index, value) in categories.values_mut().enumerate() {
            *value = index;
        }
        categories
    }
}

fn count_crossings(left: &Tree, right: &Tree) -> usize {
    let left = left.leaf_names();
    let right = right.leaf_names();
    let ties: Vec<(usize, usize)> = left
        .iter()
        .enumerate()
        .filter(|(_, name)| !name.is_empty())
        .filter_map(|(from, name)| {
            right
                .iter()
                .position(|leaf| leaf == name)
                .map(|to| (from, to))
        })
        .collect();
    let mut count = 0usize;
    for (index, (from, to)) in ties.iter().enumerate() {
        for (other_from, other_to) in ties[index + 1..].iter() {
            if (from < other_from) != (to < other_to) {
                count += 1;
            }
        }
    }
    count
}

fn improve_rotations(candidate: &mut Tree, fixed: &Tree, candidate_is_left: bool) -> bool {
    let mut internal: Vec<usize> = candidate
        .nodes()
        .iter()
        .enumerate()
        .filter(|(_, node)| node.children.len() > 1)
        .map(|(node, _)| node)
        .collect();
    internal.sort_by_key(|node| std::cmp::Reverse(candidate.ancestors(*node).len()));
    let mut improved = false;
    for node in internal {
        let before = if candidate_is_left {
            count_crossings(candidate, fixed)
        } else {
            count_crossings(fixed, candidate)
        };
        candidate.rotate(node);
        let after = if candidate_is_left {
            count_crossings(candidate, fixed)
        } else {
            count_crossings(fixed, candidate)
        };
        if after < before {
            improved = true;
        } else {
            candidate.rotate(node);
        }
    }
    improved
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
        let categories = self.annotation_categories();

        let side = band.w * self.tree_width;
        let name_size = ctx.theme.font_size - 2.0;
        let header = self.header(ctx.theme);
        let top = band.y + header;
        let tip_first = top + self.row_height / 2.0;
        let x0 = band.x + side;
        let x1 = band.right() - side;
        let corridor = (x1 - x0).max(1.0);
        let label_room = self.label_width.min(corridor * 0.34);
        let left_labels = matches!(self.labels, TangleLabels::Left | TangleLabels::Both);
        let right_labels = matches!(self.labels, TangleLabels::Right | TangleLabels::Both);
        let tie_start = x0 + if left_labels { label_room } else { 0.0 };
        let tie_end = x1 - if right_labels { label_room } else { 0.0 };

        // Quiet tree panels make the comparison corridor read as a third,
        // deliberate region rather than unused space between two plots.
        let panel = mix(ctx.theme.surface(), &ctx.theme.rule, 0.16);
        ctx.svg.rect_rounded(
            band.x,
            top,
            side,
            (band.h - header).max(1.0),
            ctx.theme.corner_radius,
            &panel,
        );
        ctx.svg.rect_rounded(
            band.right() - side,
            top,
            side,
            (band.h - header).max(1.0),
            ctx.theme.corner_radius,
            &panel,
        );

        // Ties first, under both trees.
        let ties = self.ties();
        for (index, (name, from, to)) in ties.iter().enumerate() {
            let y0 = tip_first + *from as f64 * self.row_height;
            let y1 = tip_first + *to as f64 * self.row_height;
            let crosses = self.crosses(&ties, index);
            let left_value = self.annotation(&self.left, name);
            let right_value = self.annotation(&self.right, name);
            let mismatch = self.color_by.is_some() && left_value != right_value;
            let annotation_color = left_value
                .as_ref()
                .and_then(|value| categories.get(value))
                .map(|index| ctx.theme.color(*index).to_string());
            let right_annotation_color = right_value
                .as_ref()
                .and_then(|value| categories.get(value))
                .map(|index| ctx.theme.color(*index).to_string());
            let ink = if mismatch {
                crossing.clone()
            } else if let Some(color) = &annotation_color {
                color.clone()
            } else if crosses {
                crossing.clone()
            } else {
                straight.clone()
            };
            // A tie is one taxon drawn as one curve, and whether it crosses is
            // the finding. Saying so in words is what a reader who cannot tell
            // the two colours apart has instead.
            let titled = !name.is_empty();
            if titled {
                let state = if crosses { "crossing" } else { "straight" };
                let mut title = format!("{name}, {state}");
                if let Some(key) = &self.color_by {
                    title.push_str(&format!(
                        "; left {key} {}; right {key} {}",
                        left_value.as_deref().unwrap_or("missing"),
                        right_value.as_deref().unwrap_or("missing")
                    ));
                    if mismatch {
                        title.push_str("; annotation mismatch");
                    }
                }
                ctx.svg.begin_titled(&title);
            }
            let span = (tie_end - tie_start).max(1.0);
            let centre = match self.tie_style {
                TangleTieStyle::Curved | TangleTieStyle::Ribbon => format!(
                    "M{} {}C{} {} {} {} {} {}",
                    num(tie_start),
                    num(y0),
                    num(tie_start + span * 0.35),
                    num(y0),
                    num(tie_end - span * 0.35),
                    num(y1),
                    num(tie_end),
                    num(y1)
                ),
                TangleTieStyle::Straight => format!(
                    "M{} {}L{} {}",
                    num(tie_start),
                    num(y0),
                    num(tie_end),
                    num(y1)
                ),
            };
            let width = if crosses {
                self.crossing_width
            } else {
                self.tie_width
            };
            let pattern = if crosses || mismatch {
                LinePattern::Dashed
            } else {
                LinePattern::Solid
            };
            if self.tie_style == TangleTieStyle::Ribbon {
                let half = width.max(1.2) * 0.75;
                let ribbon = format!(
                    "M{} {}C{} {} {} {} {} {}L{} {}C{} {} {} {} {} {}Z",
                    num(tie_start),
                    num(y0 - half),
                    num(tie_start + span * 0.35),
                    num(y0 - half),
                    num(tie_end - span * 0.35),
                    num(y1 - half),
                    num(tie_end),
                    num(y1 - half),
                    num(tie_end),
                    num(y1 + half),
                    num(tie_end - span * 0.35),
                    num(y1 + half),
                    num(tie_start + span * 0.35),
                    num(y0 + half),
                    num(tie_start),
                    num(y0 + half)
                );
                ctx.svg
                    .path(&ribbon, &ink, if crosses { 0.38 } else { 0.24 });
                if crosses || mismatch {
                    ctx.svg
                        .path_stroked_pattern(&centre, &ink, ctx.theme.tokens.hairline, pattern);
                }
            } else {
                ctx.svg.path_stroked_pattern(&centre, &ink, width, pattern);
            }

            if left_labels {
                let visible = fit_text(name, (label_room - 9.0).max(1.0), name_size);
                ctx.svg.text(
                    x0 + 4.0,
                    y0 + name_size * 0.35,
                    &visible,
                    &ctx.theme.muted,
                    name_size,
                    Anchor::Start,
                );
            }
            if right_labels {
                let visible = fit_text(name, (label_room - 9.0).max(1.0), name_size);
                ctx.svg.text(
                    x1 - 4.0,
                    y1 + name_size * 0.35,
                    &visible,
                    &ctx.theme.muted,
                    name_size,
                    Anchor::End,
                );
            }
            let left_mark = annotation_color.as_deref().unwrap_or(&ink);
            let right_mark = right_annotation_color.as_deref().unwrap_or(&ink);
            ctx.svg.circle_ringed(
                tie_start,
                y0,
                1.8,
                left_mark,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
            ctx.svg.circle_ringed(
                tie_end,
                y1,
                1.8,
                right_mark,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
            if titled {
                ctx.svg.end_group();
            }
        }

        for (tree, mirror) in [(&self.left, false), (&self.right, true)] {
            let area = Rect {
                x: if mirror { band.right() - side } else { band.x },
                y: top,
                w: side,
                h: band.h - header,
            };
            draw_tree_titled(
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
                matches!(
                    (mirror, self.labels),
                    (_, TangleLabels::None)
                        | (false, TangleLabels::Right)
                        | (true, TangleLabels::Left)
                ),
            );
        }

        if let Some(left) = &self.left_name {
            ctx.svg.text_bold(
                band.x,
                band.y + name_size,
                left,
                &ctx.theme.foreground,
                name_size,
                Anchor::Start,
            );
        }
        if let Some(right) = &self.right_name {
            ctx.svg.text_bold(
                band.right(),
                band.y + name_size,
                right,
                &ctx.theme.foreground,
                name_size,
                Anchor::End,
            );
        }

        if self.show_summary {
            let current = self.crossings();
            // One crossing is a crossing. The summary is the shortest line on
            // the figure and the one a reader takes the result from, so it is
            // the last place to be sloppy about it.
            let plural = if current == 1 {
                "crossing"
            } else {
                "crossings"
            };
            let crossing_text = if current == self.initial_crossings {
                format!("{current} {plural}")
            } else {
                format!("{} → {current} {plural}", self.initial_crossings)
            };
            let odd = self.unshared().len();
            let mut summary = format!("{crossing_text} · {} linked", ties.len());
            if odd > 0 {
                summary.push_str(&format!(" · {odd} unmatched"));
            }
            let width = text_width(&summary, name_size - 1.0) + 14.0;
            let x = band.x + band.w / 2.0;
            ctx.svg.rect_rounded(
                x - width / 2.0,
                band.y,
                width,
                name_size + 4.0,
                (name_size + 4.0) / 2.0,
                &panel,
            );
            ctx.svg.text_bold(
                x,
                band.y + name_size,
                &summary,
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
        assert!(svg.contains("2 unmatched"), "a missing line goes unseen");
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
    fn a_tie_names_its_taxon_and_says_whether_it_crosses() {
        // A B C D against A C B D: B and C swap and the other two do not.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(disagreeing())
            .to_svg();
        assert!(svg.contains("<title>A, straight</title>"), "{svg}");
        assert!(svg.contains("<title>B, crossing</title>"), "{svg}");
        assert!(svg.contains("<title>C, crossing</title>"));
        assert!(svg.contains("<title>D, straight</title>"));
        // Which is what a reader who cannot tell the two colours apart has
        // instead of the colour.
        assert_eq!(svg.matches("crossing</title>").count(), 2);
    }

    #[test]
    fn every_tip_of_both_trees_is_named_too() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(disagreeing().show_tips(false))
            .to_svg();
        // Four taxa, on both trees, plus the four ties.
        assert_eq!(svg.matches("<title>A</title>").count(), 2, "{svg}");
        assert_eq!(svg.matches("<title>").count(), 12);
        assert_eq!(svg.matches("<g>").count(), 12);
    }

    #[test]
    fn height_follows_the_taller_of_the_two() {
        let scale = Scale::new(&region(), 0.0, 400.0);
        let lopsided = TanglegramTrack::new(
            Tree::parse_newick("((A,B),(C,D));").unwrap(),
            Tree::parse_newick("(A,B);").unwrap(),
        );
        assert_eq!(
            lopsided.clone().show_summary(false).height(&scale),
            4.0 * 16.0
        );
        assert!(lopsided.height(&scale) > 4.0 * 16.0);
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
        assert!(svg.contains("8 unmatched"));
    }

    #[test]
    fn free_clade_rotations_reduce_crossings_without_reordering_taxa_by_hand() {
        let left = Tree::parse_newick("((A,B),(C,D));").unwrap();
        let right = Tree::parse_newick("((B,A),(D,C));").unwrap();
        let track = TanglegramTrack::new(left, right);
        assert_eq!(track.initial_crossings(), 2);

        let track = track.untangle();
        assert_eq!(track.crossings(), 0);
        assert_eq!(track.crossing_reduction(), 2);
        assert_eq!(track.shared(), vec!["A", "B", "C", "D"]);

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("2 → 0 crossings"), "{svg}");
    }

    #[test]
    fn annotations_colour_ties_and_report_disagreement_exactly() {
        let left = Tree::parse_annotated_newick(
            "((A[&country=Peru],B[&country=Spain]),(C[&country=Peru],D[&country=Spain]));",
        )
        .unwrap();
        let right = Tree::parse_annotated_newick(
            "((A[&country=Peru],C[&country=Peru]),(B[&country=Kenya],D[&country=Spain]));",
        )
        .unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TanglegramTrack::new(left, right).color_by("country"))
            .to_svg();
        assert!(svg.contains("left country Spain; right country Kenya"));
        assert!(svg.contains("annotation mismatch"));
        assert!(svg.contains("stroke-dasharray"));
    }

    #[test]
    fn ribbons_and_one_sided_labels_are_available_for_dense_comparisons() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                disagreeing()
                    .tie_style(TangleTieStyle::Ribbon)
                    .labels(TangleLabels::Left),
            )
            .to_svg();
        assert!(svg.contains("fill-opacity"), "the tie is not a ribbon");
        assert_eq!(svg.matches(">A</text>").count(), 1, "{svg}");
        assert!(!svg.contains("NaN"));
    }
}
