//! A phylogeny drawn beside its rows.
//!
//! Like the ideogram, this track does not use the shared horizontal scale: its
//! x axis is evolutionary distance, not genomic position, and the two have
//! nothing to do with each other. Unlike the ideogram, its y axis means
//! something to its neighbours, because a leaf is a row.
//!
//! # The rows are sorted to match the tree
//!
//! That is the whole point. [`SnpTrack::tree`](crate::SnpTrack::tree) puts one
//! of these in the strip beside a panel of variable sites and sorts the rows to
//! match, so a clade's shared substitutions line up into a block instead of
//! being scattered down the panel in whatever order the samples were listed.
//! [`leaf_order`] is that sort, and it never drops a row.
//!
//! # One function draws every tree in the crate
//!
//! [`draw_tree`] is a free function rather than a method on [`TreeTrack`]: the
//! standalone track, the tracks that carry a tree in a strip of their own and
//! both halves of a tanglegram all go through it. What it draws is rectangular
//! rather than diagonal, because a diagonal would imply the tree says something
//! about the space between two rows, and it says nothing about it.
//!
//! The tracks whose subject is the tree itself take the same drawing with its
//! branches named, so a clade can be pointed at for its support. A tree
//! standing beside a panel of rows does not, because the rows are named down
//! the side already and a title on every branch would be that same string a
//! second time. A tip is named on its branch only when its label is not drawn,
//! for exactly the same reason.

use std::collections::{BTreeMap, BTreeSet};

use crate::scale::Scale;
use crate::svg::{fit_text, num, text_width};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::{AnnotationValue, Placement, TimeDirection, Tree};

/// How to draw the horizontal extent of a tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeShape {
    /// Branch lengths carry distance, so tip positions mean something.
    Phylogram,
    /// Branches are counted rather than measured, and every tip lines up on
    /// the right. Use it when the lengths are missing or not to be trusted.
    Cladogram,
}

/// How a phylogenetic trait column maps values to colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitScale {
    /// Each distinct value receives a categorical palette colour.
    Categorical,
    /// Numeric values form one continuous muted-to-accent ramp.
    Continuous,
}

/// One metadata column drawn beside the terminal taxa of a tree.
#[derive(Debug, Clone, PartialEq)]
pub struct TraitColumn {
    key: String,
    label: String,
    scale: TraitScale,
    width: f64,
    show_values: bool,
}

impl TraitColumn {
    /// Builds a categorical column from annotation `key`.
    pub fn categorical(key: impl Into<String>) -> Self {
        let key = key.into();
        TraitColumn {
            label: key.clone(),
            key,
            scale: TraitScale::Categorical,
            width: 56.0,
            show_values: true,
        }
    }

    /// Builds a continuous column from numeric annotation `key`.
    pub fn continuous(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.scale = TraitScale::Continuous;
        column
    }

    /// Replaces the visible column heading without changing its metadata key.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Sets the cell width in pixels.
    pub fn width(mut self, width: f64) -> Self {
        self.width = if width.is_finite() {
            width.max(12.0)
        } else {
            56.0
        };
        self
    }

    /// Draws or hides the value text inside each cell.
    pub fn show_values(mut self, show: bool) -> Self {
        self.show_values = show;
        self
    }

    /// The annotation key read from each terminal taxon.
    pub fn key(&self) -> &str {
        &self.key
    }

    /// The colour mapping used by this column.
    pub fn scale(&self) -> TraitScale {
        self.scale
    }
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
    draw(
        svg,
        tree,
        area,
        row_pitch,
        first_row_centre,
        style,
        Titles {
            nodes: false,
            leaves: false,
        },
    );
}

/// The same drawing, with the branches named.
///
/// An internal node's riser always carries its support, since there is nothing
/// else on the page that does. A leaf's branch carries the leaf's name only
/// when `name_leaves` is set, which is the caller saying it has not drawn that
/// name anywhere else.
///
/// That switch is the whole reason this takes an argument. A tip label is
/// drawn four pixels from the branch it belongs to, at a width the track
/// reserved for it so it is never clipped, and a tooltip repeating it is a
/// pointer answering with what the reader is already looking at. Suppressed,
/// the hover falls through to nothing, which is the honest result: there is
/// no second thing to say about a tip whose name is on the page. With tips
/// hidden the title is the only way to read the tree at all, so it comes back.
pub(crate) fn draw_tree_titled(
    svg: &mut crate::svg::SvgWriter,
    tree: &Tree,
    area: Rect,
    row_pitch: f64,
    first_row_centre: f64,
    style: TreeStyle<'_>,
    name_leaves: bool,
) {
    draw(
        svg,
        tree,
        area,
        row_pitch,
        first_row_centre,
        style,
        Titles {
            nodes: true,
            leaves: name_leaves,
        },
    );
}

/// Which parts of a tree name themselves.
#[derive(Debug, Clone, Copy)]
struct Titles {
    /// Whether an internal node's riser carries its support.
    nodes: bool,
    /// Whether a leaf's branch carries the leaf's name.
    leaves: bool,
}

#[allow(clippy::too_many_arguments)]
fn draw(
    svg: &mut crate::svg::SvgWriter,
    tree: &Tree,
    area: Rect,
    row_pitch: f64,
    first_row_centre: f64,
    style: TreeStyle<'_>,
    titles: Titles,
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

        // A leaf is named on the branch that ends at it, the one piece of the
        // drawing that belongs to it alone.
        let name = if titles.leaves && node.is_leaf() {
            node.name.as_deref().unwrap_or_default()
        } else {
            ""
        };
        if !name.is_empty() {
            svg.begin_titled(name);
        }
        // Rectangular branches: along to the child's depth, then the parent's
        // riser joins its children. Diagonals would imply the tree says
        // something about the space between two rows, and it does not.
        svg.line(x0, y, x1, y, color, width);
        if !name.is_empty() {
            svg.end_group();
        }
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
        // An internal node has nothing to be called, so its riser is named by
        // the one number it does carry, and a node without one gets no group
        // rather than an empty one.
        let support = match (titles.nodes, node.support) {
            (true, Some(support)) if support.is_finite() => {
                Some(format!("clade support {}", num(support)))
            }
            _ => None,
        };
        if let Some(text) = &support {
            svg.begin_titled(text);
        }
        svg.line(x, y_of(top), x, y_of(bottom), color, width);
        if support.is_some() {
            svg.end_group();
        }
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
    time: Option<TimeAxis>,
    color_by: Option<String>,
    collapsed: BTreeSet<usize>,
    show_nodes: bool,
    trait_columns: Vec<TraitColumn>,
}

#[derive(Debug, Clone)]
struct TimeAxis {
    key: String,
    direction: TimeDirection,
    unit: Option<String>,
    show_axis: bool,
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
            time: None,
            color_by: None,
            collapsed: BTreeSet::new(),
            show_nodes: false,
            trait_columns: Vec::new(),
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

    /// Places the tree on a numeric annotation such as a decimal `date`.
    ///
    /// Every tip must carry the annotation. Missing internal values are
    /// inferred from child values and branch lengths.
    pub fn time(mut self, key: impl Into<String>) -> Self {
        self.time = Some(TimeAxis {
            key: key.into(),
            direction: TimeDirection::Increasing,
            unit: None,
            show_axis: true,
        });
        self
    }

    /// Chooses whether time values increase or decrease from root to tips.
    pub fn time_direction(mut self, direction: TimeDirection) -> Self {
        if let Some(time) = &mut self.time {
            time.direction = direction;
        }
        self
    }

    /// Adds a unit after temporal axis values.
    pub fn time_unit(mut self, unit: impl Into<String>) -> Self {
        if let Some(time) = &mut self.time {
            time.unit = Some(unit.into());
        }
        self
    }

    /// Draws or hides the temporal axis created by [`TreeTrack::time`].
    pub fn show_time_axis(mut self, show: bool) -> Self {
        if let Some(time) = &mut self.time {
            time.show_axis = show;
        }
        self
    }

    /// Colours each incoming branch by one node annotation.
    pub fn color_by(mut self, key: impl Into<String>) -> Self {
        self.color_by = Some(key.into());
        self
    }

    /// Collapses one internal node visually while preserving the source tree.
    pub fn collapse(mut self, node: usize) -> Self {
        if self
            .tree
            .nodes()
            .get(node)
            .is_some_and(|clade| !clade.is_leaf())
        {
            self.collapsed.insert(node);
        }
        self
    }

    /// Draws or hides a point at every visible internal node.
    pub fn show_nodes(mut self, show: bool) -> Self {
        self.show_nodes = show;
        self
    }

    /// Adds one metadata strip aligned to the visible terminal taxa.
    pub fn trait_column(mut self, column: TraitColumn) -> Self {
        self.trait_columns.push(column);
        self
    }

    /// Adds a categorical metadata strip.
    pub fn trait_categorical(self, key: impl Into<String>) -> Self {
        self.trait_column(TraitColumn::categorical(key))
    }

    /// Adds a continuous numeric metadata strip.
    pub fn trait_continuous(self, key: impl Into<String>) -> Self {
        self.trait_column(TraitColumn::continuous(key))
    }

    /// The tree.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Width the tip names need.
    fn tip_width(&self, theme: &Theme, scene: &TreeScene) -> f64 {
        if !self.show_tips {
            return 0.0;
        }
        scene
            .terminals
            .iter()
            .map(|node| {
                text_width(
                    &terminal_label(&self.tree, *node, &self.collapsed),
                    theme.font_size - 1.0,
                )
            })
            .fold(0.0f64, f64::max)
            + 6.0
    }

    fn axis_room(&self, theme: &Theme) -> f64 {
        self.time
            .as_ref()
            .filter(|time| time.show_axis)
            .map_or(0.0, |_| theme.font_size + theme.tokens.tick_length + 5.0)
    }

    fn trait_width(&self, theme: &Theme) -> f64 {
        if self.trait_columns.is_empty() {
            0.0
        } else {
            self.trait_columns
                .iter()
                .map(|column| column.width)
                .sum::<f64>()
                + theme.tokens.legend_gap * (self.trait_columns.len().saturating_sub(1) as f64)
                + theme.tokens.label_gap
        }
    }

    fn trait_header_room(&self) -> f64 {
        if self.trait_columns.is_empty() {
            0.0
        } else {
            18.0
        }
    }
}

impl Track for TreeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = visible_terminals(&self.tree, &self.collapsed).len().max(1) as f64;
        rows * self.row_height
            + self
                .time
                .as_ref()
                .filter(|time| time.show_axis)
                .map_or(0.0, |_| 22.0)
            + self.trait_header_room()
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
        let scene = TreeScene::new(&self.tree, self.shape, self.time.as_ref(), &self.collapsed);
        let tips = self.tip_width(ctx.theme, &scene);
        let axis_room = self.axis_room(ctx.theme);
        let traits = self.trait_width(ctx.theme);
        let header_room = self.trait_header_room();
        let area = Rect {
            x: band.x,
            y: band.y + header_room,
            w: (band.w - tips - traits).max(1.0),
            h: (band.h - axis_room - header_room).max(1.0),
        };

        draw_tree_scene(
            ctx,
            &self.tree,
            &scene,
            area,
            self.row_height,
            &color,
            self.line_width,
            self.color_by.as_deref(),
            self.show_nodes,
            !self.show_tips,
        );

        if self.show_tips {
            let size = ctx.theme.font_size - 1.0;
            for (row, node) in scene.terminals.iter().enumerate() {
                let name = terminal_label(&self.tree, *node, &self.collapsed);
                ctx.svg.text(
                    area.right() + 4.0,
                    area.y + self.row_height / 2.0 + row as f64 * self.row_height + size * 0.35,
                    &name,
                    &ctx.theme.muted,
                    size,
                    crate::svg::Anchor::Start,
                );
            }
        }
        draw_trait_columns(
            ctx,
            &self.tree,
            &scene,
            &self.collapsed,
            area,
            tips,
            &self.trait_columns,
            self.row_height,
        );
        if let Some(time) = self.time.as_ref().filter(|time| time.show_axis) {
            draw_time_axis(ctx, &scene, area, time);
        }
    }
}

struct TreeScene {
    placements: Vec<Option<Placement>>,
    source_placements: Vec<Placement>,
    terminals: Vec<usize>,
    minimum: f64,
    maximum: f64,
    temporal: bool,
    direction: TimeDirection,
}

impl TreeScene {
    fn new(
        tree: &Tree,
        shape: TreeShape,
        time: Option<&TimeAxis>,
        collapsed: &BTreeSet<usize>,
    ) -> Self {
        let timed = time.and_then(|axis| tree.time_layout(&axis.key, axis.direction));
        let temporal = timed.is_some();
        let source_placements = timed.unwrap_or_else(|| tree.layout(shape == TreeShape::Cladogram));
        let direction = time.map_or(TimeDirection::Increasing, |axis| axis.direction);
        let terminals = visible_terminals(tree, collapsed);
        let mut rows = vec![None; tree.nodes().len()];
        for (row, node) in terminals.iter().enumerate() {
            rows[*node] = Some(row as f64);
        }
        let visible = visible_nodes(tree, collapsed);
        for node in postorder_nodes(tree) {
            if !visible[node] || rows[node].is_some() {
                continue;
            }
            let children: Vec<f64> = tree.nodes()[node]
                .children
                .iter()
                .filter(|child| visible[**child])
                .filter_map(|child| rows[*child])
                .collect();
            if !children.is_empty() {
                rows[node] = Some(children.iter().sum::<f64>() / children.len() as f64);
            }
        }
        let placements: Vec<Option<Placement>> = source_placements
            .iter()
            .map(|placement| {
                visible[placement.node].then_some(Placement {
                    row: rows[placement.node].unwrap_or(0.0),
                    ..*placement
                })
            })
            .collect();
        let (minimum, maximum) = source_placements
            .iter()
            .map(|placement| placement.depth)
            .filter(|value| value.is_finite())
            .fold((f64::MAX, f64::MIN), |(minimum, maximum), value| {
                (minimum.min(value), maximum.max(value))
            });
        let (minimum, maximum) = if minimum.is_finite() && maximum.is_finite() {
            (minimum, maximum)
        } else {
            (0.0, 1.0)
        };
        TreeScene {
            placements,
            source_placements,
            terminals,
            minimum,
            maximum,
            temporal,
            direction,
        }
    }

    fn x(&self, area: Rect, value: f64) -> f64 {
        let span = self.maximum - self.minimum;
        let fraction = if span <= 0.0 {
            0.0
        } else {
            match self.direction {
                TimeDirection::Increasing => (value - self.minimum) / span,
                TimeDirection::Decreasing if self.temporal => (self.maximum - value) / span,
                TimeDirection::Decreasing => (value - self.minimum) / span,
            }
        };
        area.x + fraction.clamp(0.0, 1.0) * area.w
    }
}

fn visible_nodes(tree: &Tree, collapsed: &BTreeSet<usize>) -> Vec<bool> {
    let mut visible = vec![false; tree.nodes().len()];
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        visible[node] = true;
        if collapsed.contains(&node) {
            continue;
        }
        for child in tree.nodes()[node].children.iter().rev() {
            stack.push(*child);
        }
    }
    visible
}

fn visible_terminals(tree: &Tree, collapsed: &BTreeSet<usize>) -> Vec<usize> {
    let mut terminals = Vec::new();
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        let clade = &tree.nodes()[node];
        if clade.is_leaf() || collapsed.contains(&node) {
            terminals.push(node);
            continue;
        }
        for child in clade.children.iter().rev() {
            stack.push(*child);
        }
    }
    terminals
}

fn postorder_nodes(tree: &Tree) -> Vec<usize> {
    let mut order = Vec::with_capacity(tree.nodes().len());
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        order.push(node);
        stack.extend(tree.nodes()[node].children.iter().copied());
    }
    order.reverse();
    order
}

fn terminal_label(tree: &Tree, node: usize, collapsed: &BTreeSet<usize>) -> String {
    let name = tree.nodes()[node].name.as_deref().unwrap_or("clade");
    if collapsed.contains(&node) {
        format!("{} ({} tips)", name, tree.clade_size(node))
    } else {
        name.to_string()
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_tree_scene(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    scene: &TreeScene,
    area: Rect,
    row_pitch: f64,
    default_color: &str,
    width: f64,
    color_by: Option<&str>,
    show_nodes: bool,
    name_leaves: bool,
) {
    let colors = branch_colors(tree, scene, color_by, ctx.theme, default_color);
    let y_of = |row: f64| area.y + row_pitch / 2.0 + row * row_pitch;

    for placement in scene.placements.iter().flatten() {
        let node = &tree.nodes()[placement.node];
        let Some(parent) = node.parent else {
            continue;
        };
        let Some(parent_placement) = scene.placements[parent] else {
            continue;
        };
        let x0 = scene.x(area, parent_placement.depth);
        let x1 = scene.x(area, placement.depth);
        let y = y_of(placement.row);
        let title = branch_title(tree, placement.node, color_by, name_leaves, false);
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        ctx.svg.line(x0, y, x1, y, &colors[placement.node], width);
        if title.is_some() {
            ctx.svg.end_group();
        }
    }

    for placement in scene.placements.iter().flatten() {
        let node = &tree.nodes()[placement.node];
        if node.is_leaf() || node.children.is_empty() || scene.terminals.contains(&placement.node) {
            continue;
        }
        let rows: Vec<f64> = node
            .children
            .iter()
            .filter_map(|child| scene.placements[*child].map(|value| value.row))
            .collect();
        if rows.is_empty() {
            continue;
        }
        let (top, bottom) = rows.iter().fold((f64::MAX, f64::MIN), |(lo, hi), row| {
            (lo.min(*row), hi.max(*row))
        });
        let x = scene.x(area, placement.depth);
        let title = branch_title(tree, placement.node, color_by, false, true);
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        ctx.svg.line(
            x,
            y_of(top),
            x,
            y_of(bottom),
            &colors[placement.node],
            width,
        );
        if title.is_some() {
            ctx.svg.end_group();
        }
        if show_nodes {
            ctx.svg.circle_ringed(
                x,
                y_of(placement.row),
                ctx.theme.tokens.marker_radius * 0.65,
                &colors[placement.node],
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
        }
    }

    for node in &scene.terminals {
        if !tree.nodes()[*node].is_leaf() {
            let placement = scene.placements[*node].unwrap();
            let start = scene.x(area, placement.depth);
            let far = tree
                .descendants(*node)
                .into_iter()
                .map(|descendant| scene.x(area, scene.source_placements[descendant].depth))
                .fold(start, f64::max);
            let y = y_of(placement.row);
            let half = row_pitch * 0.34;
            let d = format!(
                "M {} {} L {} {} L {} {} Z",
                num(start),
                num(y),
                num(far.max(start + 2.0)),
                num(y - half),
                num(far.max(start + 2.0)),
                num(y + half)
            );
            let title = format!(
                "{} ({} tips)",
                tree.nodes()[*node].name.as_deref().unwrap_or("clade"),
                tree.clade_size(*node)
            );
            ctx.svg.begin_titled(&title);
            ctx.svg.path(&d, &colors[*node], 0.28);
            ctx.svg.end_group();
        }
    }
}

fn branch_colors(
    tree: &Tree,
    scene: &TreeScene,
    key: Option<&str>,
    theme: &Theme,
    default_color: &str,
) -> Vec<String> {
    let mut colors = vec![default_color.to_string(); tree.nodes().len()];
    let Some(key) = key else {
        return colors;
    };
    let values: Vec<Option<&AnnotationValue>> = (0..tree.nodes().len())
        .map(|node| inherited_annotation(tree, node, key))
        .collect();
    let numeric: Vec<f64> = values
        .iter()
        .filter_map(|value| value.and_then(AnnotationValue::as_number))
        .collect();
    let all_numeric = numeric.len() == values.iter().filter(|value| value.is_some()).count()
        && !numeric.is_empty();
    if all_numeric {
        let minimum = numeric.iter().copied().fold(f64::MAX, f64::min);
        let maximum = numeric.iter().copied().fold(f64::MIN, f64::max);
        for node in scene
            .placements
            .iter()
            .flatten()
            .map(|placement| placement.node)
        {
            if let Some(value) = values[node].and_then(AnnotationValue::as_number) {
                let fraction = if maximum <= minimum {
                    1.0
                } else {
                    (value - minimum) / (maximum - minimum)
                };
                colors[node] = mix(&theme.muted, &theme.accent, fraction);
            }
        }
    } else {
        let mut categories = BTreeMap::new();
        for node in scene
            .placements
            .iter()
            .flatten()
            .map(|placement| placement.node)
        {
            let Some(value) = values[node] else {
                continue;
            };
            let value = value.to_string();
            let next = categories.len();
            let index = *categories.entry(value).or_insert(next);
            colors[node] = theme.color(index).to_string();
        }
    }
    colors
}

fn inherited_annotation<'a>(tree: &'a Tree, node: usize, key: &str) -> Option<&'a AnnotationValue> {
    tree.annotation(node, key).or_else(|| {
        tree.ancestors(node)
            .into_iter()
            .find_map(|ancestor| tree.annotation(ancestor, key))
    })
}

fn branch_title(
    tree: &Tree,
    node: usize,
    color_by: Option<&str>,
    name_leaf: bool,
    include_support: bool,
) -> Option<String> {
    let clade = &tree.nodes()[node];
    let mut parts = Vec::new();
    if name_leaf && clade.is_leaf() {
        if let Some(name) = &clade.name {
            if !name.is_empty() {
                parts.push(name.clone());
            }
        }
    }
    if include_support {
        if let Some(support) = clade.support.filter(|value| value.is_finite()) {
            parts.push(format!("clade support {}", num(support)));
        }
    }
    if let Some(key) = color_by {
        if let Some(value) = inherited_annotation(tree, node, key) {
            parts.push(format!("{key} {value}"));
        }
    }
    (!parts.is_empty()).then(|| parts.join("; "))
}

#[allow(clippy::too_many_arguments)]
fn draw_trait_columns(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    scene: &TreeScene,
    collapsed: &BTreeSet<usize>,
    area: Rect,
    tip_width: f64,
    columns: &[TraitColumn],
    row_pitch: f64,
) {
    if columns.is_empty() {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let mut x = area.right() + tip_width + ctx.theme.tokens.label_gap;

    for column in columns {
        let values: Vec<Option<&AnnotationValue>> = scene
            .terminals
            .iter()
            .map(|node| inherited_annotation(tree, *node, &column.key))
            .collect();
        let categories: BTreeMap<String, usize> = values
            .iter()
            .filter_map(|value| value.map(ToString::to_string))
            .fold(BTreeMap::new(), |mut categories, value| {
                let next = categories.len();
                categories.entry(value).or_insert(next);
                categories
            });
        let numeric: Vec<f64> = values
            .iter()
            .filter_map(|value| value.and_then(AnnotationValue::as_number))
            .filter(|value| value.is_finite())
            .collect();
        let minimum = numeric.iter().copied().fold(f64::MAX, f64::min);
        let maximum = numeric.iter().copied().fold(f64::MIN, f64::max);

        let heading = fit_text(&column.label, column.width, size);
        ctx.svg.text(
            x + column.width / 2.0,
            area.y - 5.0,
            &heading,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );

        for (row, node) in scene.terminals.iter().enumerate() {
            let y = area.y + row as f64 * row_pitch + 1.0;
            let height = (row_pitch - 2.0).max(1.0);
            let name = terminal_label(tree, *node, collapsed);
            let value = values[row];
            let fill = match (column.scale, value) {
                (TraitScale::Categorical, Some(value)) => {
                    let value = value.to_string();
                    categories
                        .get(&value)
                        .map(|index| ctx.theme.color(*index).to_string())
                }
                (TraitScale::Continuous, Some(value)) => value.as_number().and_then(|value| {
                    value.is_finite().then(|| {
                        let fraction = if maximum <= minimum {
                            1.0
                        } else {
                            (value - minimum) / (maximum - minimum)
                        };
                        mix(&ctx.theme.muted, &ctx.theme.accent, fraction)
                    })
                }),
                _ => None,
            };
            let displayed = value.map(ToString::to_string);
            let title = match &displayed {
                Some(value) => format!("{name}; {} {value}", column.key),
                None => format!("{name}; {} missing", column.key),
            };
            ctx.svg.begin_titled(&title);
            if let Some(fill) = &fill {
                ctx.svg.rect_rounded(
                    x,
                    y,
                    column.width,
                    height,
                    ctx.theme.corner_radius.min(2.0),
                    fill,
                );
                if column.show_values {
                    if let Some(value) = &displayed {
                        let value = fit_text(value, column.width - 4.0, size);
                        ctx.svg.text(
                            x + column.width / 2.0,
                            y + height / 2.0 + size * 0.35,
                            &value,
                            contrast_ink(fill),
                            size,
                            crate::svg::Anchor::Middle,
                        );
                    }
                }
            } else {
                ctx.svg.rect_outline(
                    x,
                    y,
                    column.width,
                    height,
                    &ctx.theme.rule,
                    ctx.theme.tokens.hairline,
                );
                if column.show_values {
                    ctx.svg.text(
                        x + column.width / 2.0,
                        y + height / 2.0 + size * 0.35,
                        "—",
                        &ctx.theme.muted,
                        size,
                        crate::svg::Anchor::Middle,
                    );
                }
            }
            ctx.svg.end_group();
        }
        x += column.width + ctx.theme.tokens.legend_gap;
    }
}

fn draw_time_axis(ctx: &mut DrawContext<'_>, scene: &TreeScene, area: Rect, time: &TimeAxis) {
    if !scene.temporal {
        return;
    }
    let y = area.bottom() + 2.0;
    ctx.svg.line(
        area.x,
        y,
        area.right(),
        y,
        &ctx.theme.foreground,
        ctx.theme.tokens.hairline,
    );
    let size = ctx.theme.font_size - 1.0;
    for index in 0..=2 {
        let fraction = index as f64 / 2.0;
        let value = scene.minimum + fraction * (scene.maximum - scene.minimum);
        let x = scene.x(area, value);
        ctx.svg.line(
            x,
            y,
            x,
            y + ctx.theme.tokens.tick_length,
            &ctx.theme.foreground,
            ctx.theme.tokens.hairline,
        );
        let label = match &time.unit {
            Some(unit) => format!("{} {unit}", num(value)),
            None => num(value),
        };
        ctx.svg.text(
            x,
            y + ctx.theme.tokens.tick_length + size,
            &label,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );
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
    fn a_leaf_is_named_on_its_own_branch_when_its_label_is_not_drawn() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()).show_tips(false))
            .to_svg();
        for tip in ["A", "B", "C", "D"] {
            assert!(svg.contains(&format!("<title>{tip}</title>")), "{svg}");
        }
        // One of the two clades carries a support value and the other does
        // not, so only one of them opens a group.
        assert!(svg.contains("<title>clade support 0.9</title>"), "{svg}");
        assert_eq!(svg.matches("clade support").count(), 1);
        assert_eq!(svg.matches("<title>").count(), 5);
        assert_eq!(svg.matches("<g>").count(), 5);
    }

    #[test]
    fn a_tip_whose_label_is_drawn_is_not_named_a_second_time() {
        // The label sits four pixels from the branch, at a width the track
        // reserved for it, so it is never clipped. A tooltip carrying that
        // same string is the pointer answering with what is already on screen.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()))
            .to_svg();
        for tip in ["A", "B", "C", "D"] {
            assert!(svg.contains(&format!(">{tip}</text>")), "the label, {svg}");
            assert!(
                !svg.contains(&format!("<title>{tip}</title>")),
                "the same string twice, {svg}"
            );
        }
        // The clade support is the one thing no label carries, so it stays.
        assert!(svg.contains("<title>clade support 0.9</title>"), "{svg}");
        assert_eq!(svg.matches("<title>").count(), 1);
    }

    #[test]
    fn an_unnamed_leaf_opens_no_group() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(Tree::parse_newick("((,),(,));").unwrap()).show_tips(false))
            .to_svg();
        assert!(!svg.contains("<title>"), "{svg}");
    }

    #[test]
    fn a_tree_drawn_beside_rows_is_left_unnamed() {
        // The panel tracks name their rows down the side already, so the plain
        // drawing has to stay plain: a title on every branch there would be
        // the same string twice.
        let mut svg = crate::svg::SvgWriter::new();
        draw_tree(
            &mut svg,
            &tree(),
            Rect {
                x: 0.0,
                y: 0.0,
                w: 100.0,
                h: 60.0,
            },
            15.0,
            7.5,
            TreeStyle {
                shape: TreeShape::Phylogram,
                color: "#111111",
                width: 1.0,
                mirror: false,
            },
        );
        let out = svg.finish(100.0, 60.0, "none", "sans-serif");
        assert!(!out.contains("<title"), "{out}");
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
    fn a_time_tree_draws_calendar_values_on_its_axis() {
        let tree = Tree::parse_annotated_newick(
            "((A[&date=2024]:1,B[&date=2025]:2)AB:1,C[&date=2023]:3);",
        )
        .unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree).time("date").time_unit("year"))
            .to_svg();
        for label in ["2021 year", "2023 year", "2025 year"] {
            assert!(svg.contains(&format!(">{label}</text>")), "{label}: {svg}");
        }
    }

    #[test]
    fn branch_annotations_drive_colour_and_accessible_text() {
        let tree = Tree::parse_annotated_newick(
            "((A[&country=Peru]:1,B[&country=Chile]:1)[&country=Peru]:1,C[&country=Chile]:2);",
        )
        .unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree).color_by("country"))
            .to_svg();
        assert!(svg.contains("country Peru"), "{svg}");
        assert!(svg.contains("country Chile"), "{svg}");
        assert!(svg.contains("#0072b2"), "first categorical colour: {svg}");
        assert!(svg.contains("#d55e00"), "second categorical colour: {svg}");
    }

    #[test]
    fn visual_collapse_keeps_the_source_tree_and_names_the_triangle() {
        let tree = Tree::parse_newick("((A:1,B:1)outbreak:1,C:2);").unwrap();
        let outbreak = tree.node_named("outbreak").unwrap();
        let track = TreeTrack::new(tree).collapse(outbreak);
        assert_eq!(track.tree().leaf_names(), ["A", "B", "C"]);
        assert_eq!(track.height(&Scale::new(&region(), 0.0, 100.0)), 30.0);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("outbreak (2 tips)"), "{svg}");
        assert!(svg.contains("fill-opacity=\"0.28\""), "{svg}");
    }

    #[test]
    fn internal_node_points_are_optional() {
        let plain = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()))
            .to_svg();
        let marked = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()).show_nodes(true))
            .to_svg();
        assert!(marked.matches("<circle").count() > plain.matches("<circle").count());
    }

    #[test]
    fn trait_columns_align_exact_metadata_with_terminal_taxa() {
        let tree = Tree::parse_annotated_newick(
            "((A[&country=Peru,coverage=18]:1,B[&country=Chile]:1):1,C[&country=Peru,coverage=42]:2);",
        )
        .unwrap();
        let track = TreeTrack::new(tree)
            .trait_column(TraitColumn::categorical("country").label("Country"))
            .trait_column(TraitColumn::continuous("coverage").label("Depth"));
        assert_eq!(
            track.height(&Scale::new(&region(), 0.0, 100.0)),
            3.0 * 15.0 + 18.0
        );
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        for text in [
            ">Country</text>",
            ">Depth</text>",
            ">Peru</text>",
            ">Chile</text>",
        ] {
            assert!(svg.contains(text), "{text}: {svg}");
        }
        for title in [
            "A; country Peru",
            "A; coverage 18",
            "B; country Chile",
            "B; coverage missing",
            "C; coverage 42",
        ] {
            assert!(svg.contains(&format!("<title>{title}</title>")), "{svg}");
        }
        assert!(svg.contains("#4b5563"), "continuous minimum: {svg}");
        assert!(svg.contains("#0072b2"), "continuous maximum: {svg}");
        assert!(svg.contains(">—</text>"), "missing value: {svg}");
    }

    #[test]
    fn trait_column_builders_expose_their_mapping() {
        let categorical = TraitColumn::categorical("lineage");
        let continuous = TraitColumn::continuous("clock_rate");
        assert_eq!(categorical.key(), "lineage");
        assert_eq!(categorical.scale(), TraitScale::Categorical);
        assert_eq!(continuous.key(), "clock_rate");
        assert_eq!(continuous.scale(), TraitScale::Continuous);
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
