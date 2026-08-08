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

/// Coordinate projection used by [`TreeTrack`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TreeProjection {
    /// Root on the left and terminal taxa in rows on the right.
    #[default]
    Rectangular,
    /// Root and tips arranged on concentric radii.
    Circular,
}

/// Direction in which branches radiate in a circular tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RadialDirection {
    /// Root nearest the centre and terminal taxa towards the circumference.
    #[default]
    Outward,
    /// Root at the circumference and terminal taxa towards the centre.
    Inward,
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
    ring_width: f64,
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
            ring_width: 10.0,
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

    /// Sets the thickness of this trait when drawn as a circular ring.
    pub fn ring_width(mut self, width: f64) -> Self {
        self.ring_width = if width.is_finite() {
            width.clamp(2.0, 24.0)
        } else {
            10.0
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
    projection: TreeProjection,
    radial: RadialLayout,
    color: Option<String>,
    line_width: f64,
    show_tips: bool,
    time: Option<TimeAxis>,
    color_by: Option<String>,
    collapsed: BTreeSet<usize>,
    show_nodes: bool,
    trait_columns: Vec<TraitColumn>,
}

#[derive(Debug, Clone, Copy)]
struct RadialLayout {
    start_degrees: f64,
    sweep_degrees: f64,
    direction: RadialDirection,
    inner_radius: f64,
    size: f64,
}

impl Default for RadialLayout {
    fn default() -> Self {
        RadialLayout {
            start_degrees: -90.0,
            sweep_degrees: 360.0,
            direction: RadialDirection::Outward,
            inner_radius: 0.08,
            size: 440.0,
        }
    }
}

fn finite_between(value: f64, minimum: f64, maximum: f64, fallback: f64) -> f64 {
    if value.is_finite() {
        value.clamp(minimum, maximum)
    } else {
        fallback
    }
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
            projection: TreeProjection::Rectangular,
            radial: RadialLayout::default(),
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

    /// Chooses rectangular or circular coordinates.
    pub fn projection(mut self, projection: TreeProjection) -> Self {
        self.projection = projection;
        self
    }

    /// Draws a complete circular tree radiating outwards by default.
    pub fn circular(mut self) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.sweep_degrees = 360.0;
        self
    }

    /// Draws a circular fan covering `sweep_degrees` clockwise.
    pub fn fan(mut self, sweep_degrees: f64) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.sweep_degrees = finite_between(sweep_degrees, 10.0, 359.0, 240.0);
        self
    }

    /// Sets the angle where a circular tree begins, in clockwise degrees.
    ///
    /// Zero is three o'clock and -90 is twelve o'clock.
    pub fn radial_start(mut self, degrees: f64) -> Self {
        self.projection = TreeProjection::Circular;
        if degrees.is_finite() {
            self.radial.start_degrees = degrees;
        }
        self
    }

    /// Sets the clockwise angular span of a circular tree in degrees.
    pub fn radial_sweep(mut self, degrees: f64) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.sweep_degrees = finite_between(degrees, 10.0, 360.0, 360.0);
        self
    }

    /// Chooses whether tips point away from or towards the centre.
    pub fn radial_direction(mut self, direction: RadialDirection) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.direction = direction;
        self
    }

    /// Sets the central gap as a fraction of the tree radius.
    pub fn inner_radius(mut self, fraction: f64) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.inner_radius = finite_between(fraction, 0.0, 0.85, 0.08);
        self
    }

    /// Sets the requested height of the circular drawing in pixels.
    pub fn radial_size(mut self, size: f64) -> Self {
        self.radial.size = if size.is_finite() {
            size.max(120.0)
        } else {
            440.0
        };
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

    fn draw_rectangular(&self, ctx: &mut DrawContext<'_>) {
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

impl Track for TreeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        match self.projection {
            TreeProjection::Rectangular => {
                let rows = visible_terminals(&self.tree, &self.collapsed).len().max(1) as f64;
                rows * self.row_height
                    + self
                        .time
                        .as_ref()
                        .filter(|time| time.show_axis)
                        .map_or(0.0, |_| 22.0)
                    + self.trait_header_room()
            }
            TreeProjection::Circular => self.radial.size + self.trait_header_room(),
        }
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        match self.projection {
            TreeProjection::Rectangular => self.draw_rectangular(ctx),
            TreeProjection::Circular => {
                draw_radial_track(self, ctx);
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct RadialGeometry {
    cx: f64,
    cy: f64,
    tree_inner: f64,
    tree_outer: f64,
    ring_outer: f64,
    label_radius: f64,
    start: f64,
    sweep: f64,
    terminals: usize,
    direction: RadialDirection,
}

impl RadialGeometry {
    fn new(track: &TreeTrack, theme: &Theme, scene: &TreeScene, area: Rect) -> Self {
        let label_extent = if track.show_tips {
            scene
                .terminals
                .iter()
                .map(|node| {
                    text_width(
                        &terminal_label(&track.tree, *node, &track.collapsed),
                        theme.font_size - 1.0,
                    )
                })
                .fold(0.0f64, f64::max)
                + 6.0
        } else {
            4.0
        };
        let gap = theme.tokens.legend_gap.clamp(1.0, 4.0);
        let ring_room = if track.trait_columns.is_empty() {
            0.0
        } else {
            track
                .trait_columns
                .iter()
                .map(|column| column.ring_width)
                .sum::<f64>()
                + gap * track.trait_columns.len() as f64
        };
        let half = (area.w.min(area.h) / 2.0 - 4.0).max(2.0);
        let (tree_outer, tree_inner, ring_outer, label_radius) = match track.radial.direction {
            RadialDirection::Outward => {
                let ring_outer = (half - label_extent).max(4.0);
                let tree_outer = (ring_outer - ring_room).max(2.0);
                let tree_inner = tree_outer * track.radial.inner_radius;
                (tree_outer, tree_inner, ring_outer, ring_outer + 4.0)
            }
            RadialDirection::Inward => {
                let ring_outer = half;
                let tree_outer = (ring_outer - ring_room).max(2.0);
                let requested = tree_outer * track.radial.inner_radius;
                let tree_inner = if track.show_tips {
                    requested.max(label_extent + 4.0).min(tree_outer * 0.9)
                } else {
                    requested
                };
                (
                    tree_outer,
                    tree_inner,
                    ring_outer,
                    (tree_inner - 4.0).max(0.0),
                )
            }
        };
        RadialGeometry {
            cx: area.x + area.w / 2.0,
            cy: area.y + area.h / 2.0,
            tree_inner,
            tree_outer,
            ring_outer,
            label_radius,
            start: track.radial.start_degrees.to_radians(),
            sweep: track.radial.sweep_degrees.to_radians(),
            terminals: scene.terminals.len().max(1),
            direction: track.radial.direction,
        }
    }

    fn angle(&self, row: f64) -> f64 {
        if self.terminals == 1 {
            return if self.full_circle() {
                self.start
            } else {
                self.start + self.sweep / 2.0
            };
        }
        let denominator = if self.full_circle() {
            self.terminals as f64
        } else {
            (self.terminals - 1) as f64
        };
        self.start + self.sweep * row / denominator
    }

    fn angular_step(&self) -> f64 {
        if self.terminals <= 1 {
            self.sweep
        } else if self.full_circle() {
            self.sweep / self.terminals as f64
        } else {
            self.sweep / (self.terminals - 1) as f64
        }
    }

    fn full_circle(&self) -> bool {
        self.sweep >= std::f64::consts::TAU - 1e-6
    }

    fn point(&self, radius: f64, angle: f64) -> (f64, f64) {
        (
            self.cx + angle.cos() * radius,
            self.cy + angle.sin() * radius,
        )
    }

    fn radius(&self, scene: &TreeScene, value: f64) -> f64 {
        let fraction = scene.fraction(value);
        match self.direction {
            RadialDirection::Outward => {
                self.tree_inner + fraction * (self.tree_outer - self.tree_inner)
            }
            RadialDirection::Inward => {
                self.tree_outer - fraction * (self.tree_outer - self.tree_inner)
            }
        }
    }

    fn terminal_boundary(&self) -> f64 {
        match self.direction {
            RadialDirection::Outward => self.tree_outer,
            RadialDirection::Inward => self.tree_inner,
        }
    }
}

fn draw_radial_track(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    let color = track
        .color
        .clone()
        .unwrap_or_else(|| ctx.theme.foreground.clone());
    let scene = TreeScene::new(
        &track.tree,
        track.shape,
        track.time.as_ref(),
        &track.collapsed,
    );
    let header_room = track.trait_header_room();
    let area = Rect {
        x: ctx.band.x,
        y: ctx.band.y + header_room,
        w: ctx.band.w,
        h: (ctx.band.h - header_room).max(1.0),
    };
    let geometry = RadialGeometry::new(track, ctx.theme, &scene, area);
    let colors = branch_colors(
        &track.tree,
        &scene,
        track.color_by.as_deref(),
        ctx.theme,
        &color,
    );

    if let Some(time) = track.time.as_ref().filter(|time| time.show_axis) {
        draw_radial_time_axis(ctx, &scene, &geometry, time);
    }
    draw_radial_padding(track, ctx, &scene, &geometry);
    draw_radial_branches(track, ctx, &scene, &geometry, &colors);
    draw_radial_collapsed(track, ctx, &scene, &geometry, &colors);
    draw_trait_rings(track, ctx, &scene, &geometry);
    draw_radial_labels(track, ctx, &scene, &geometry);
    draw_trait_ring_headings(track, ctx);
}

fn draw_radial_padding(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if !track.show_tips && track.trait_columns.is_empty() {
        return;
    }
    let boundary = geometry.terminal_boundary();
    for (row, node) in scene.terminals.iter().enumerate() {
        if !track.tree.nodes()[*node].is_leaf() {
            continue;
        }
        let placement = scene.placements[*node].unwrap();
        let radius = geometry.radius(scene, placement.depth);
        if (radius - boundary).abs() <= 0.5 {
            continue;
        }
        let angle = geometry.angle(row as f64);
        let (x0, y0) = geometry.point(radius, angle);
        let (x1, y1) = geometry.point(boundary, angle);
        ctx.svg
            .line(x0, y0, x1, y1, &ctx.theme.rule, ctx.theme.tokens.hairline);
    }
}

fn draw_radial_branches(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    colors: &[String],
) {
    for placement in scene.placements.iter().flatten() {
        let node = &track.tree.nodes()[placement.node];
        let Some(parent) = node.parent else {
            continue;
        };
        let Some(parent_placement) = scene.placements[parent] else {
            continue;
        };
        let angle = geometry.angle(placement.row);
        let (x0, y0) = geometry.point(geometry.radius(scene, parent_placement.depth), angle);
        let (x1, y1) = geometry.point(geometry.radius(scene, placement.depth), angle);
        let title = branch_title(
            &track.tree,
            placement.node,
            track.color_by.as_deref(),
            !track.show_tips,
            false,
        );
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        ctx.svg
            .line(x0, y0, x1, y1, &colors[placement.node], track.line_width);
        if title.is_some() {
            ctx.svg.end_group();
        }
    }

    for placement in scene.placements.iter().flatten() {
        let node = &track.tree.nodes()[placement.node];
        if node.is_leaf() || scene.terminals.contains(&placement.node) {
            continue;
        }
        let angles: Vec<f64> = node
            .children
            .iter()
            .filter_map(|child| scene.placements[*child])
            .map(|child| geometry.angle(child.row))
            .collect();
        let (Some(start), Some(end)) = (angles.first(), angles.last()) else {
            continue;
        };
        let radius = geometry.radius(scene, placement.depth);
        if radius > 0.5 && (end - start).abs() > 1e-9 {
            let title = branch_title(
                &track.tree,
                placement.node,
                track.color_by.as_deref(),
                false,
                true,
            );
            if let Some(title) = &title {
                ctx.svg.begin_titled(title);
            }
            ctx.svg.path_stroked(
                &radial_arc_path(geometry, radius, *start, *end),
                &colors[placement.node],
                track.line_width,
            );
            if title.is_some() {
                ctx.svg.end_group();
            }
        }
        if track.show_nodes {
            let angle = geometry.angle(placement.row);
            let (x, y) = geometry.point(radius, angle);
            ctx.svg.circle_ringed(
                x,
                y,
                ctx.theme.tokens.marker_radius * 0.65,
                &colors[placement.node],
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
        }
    }
}

fn draw_radial_labels(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if !track.show_tips {
        return;
    }
    let size = ctx.theme.font_size - 1.0;
    for (row, node) in scene.terminals.iter().enumerate() {
        let angle = geometry.angle(row as f64);
        let (x, y) = geometry.point(geometry.label_radius, angle);
        let degrees = angle.to_degrees().rem_euclid(360.0);
        let right = angle.cos() >= 0.0;
        let (rotation, anchor) = match (geometry.direction, right) {
            (RadialDirection::Outward, true) => (degrees, crate::svg::Anchor::Start),
            (RadialDirection::Outward, false) => (degrees + 180.0, crate::svg::Anchor::End),
            (RadialDirection::Inward, true) => (degrees, crate::svg::Anchor::End),
            (RadialDirection::Inward, false) => (degrees + 180.0, crate::svg::Anchor::Start),
        };
        ctx.svg.text_rotated(
            (x, y + size * 0.32),
            rotation,
            &terminal_label(&track.tree, *node, &track.collapsed),
            &ctx.theme.muted,
            size,
            anchor,
        );
    }
}

fn radial_arc_path(geometry: &RadialGeometry, radius: f64, start: f64, end: f64) -> String {
    let delta = (end - start).abs();
    let (x0, y0) = geometry.point(radius, start);
    if delta >= std::f64::consts::TAU - 1e-6 {
        let middle = start + std::f64::consts::PI;
        let (xm, ym) = geometry.point(radius, middle);
        return format!(
            "M {} {} A {} {} 0 0 1 {} {} A {} {} 0 0 1 {} {}",
            num(x0),
            num(y0),
            num(radius),
            num(radius),
            num(xm),
            num(ym),
            num(radius),
            num(radius),
            num(x0),
            num(y0)
        );
    }
    let (x1, y1) = geometry.point(radius, end);
    format!(
        "M {} {} A {} {} 0 {} 1 {} {}",
        num(x0),
        num(y0),
        num(radius),
        num(radius),
        usize::from(delta > std::f64::consts::PI),
        num(x1),
        num(y1)
    )
}

fn draw_radial_time_axis(
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    time: &TimeAxis,
) {
    if !scene.temporal {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    for index in 0..=2 {
        let fraction = index as f64 / 2.0;
        let value = scene.minimum + fraction * (scene.maximum - scene.minimum);
        let radius = geometry.radius(scene, value);
        if radius > 0.5 {
            ctx.svg.path_stroked(
                &radial_arc_path(
                    geometry,
                    radius,
                    geometry.start,
                    geometry.start + geometry.sweep,
                ),
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
        }
        let label = match &time.unit {
            Some(unit) => format!("{} {unit}", num(value)),
            None => num(value),
        };
        let (x, y) = geometry.point((radius - 4.0).max(0.0), geometry.start);
        let rotation = upright_tangent(geometry.start);
        ctx.svg.text_rotated(
            (x, y - 2.0),
            rotation,
            &label,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );
    }
}

fn draw_radial_collapsed(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    colors: &[String],
) {
    for (row, node) in scene.terminals.iter().enumerate() {
        if track.tree.nodes()[*node].is_leaf() {
            continue;
        }
        let placement = scene.placements[*node].unwrap();
        let start_radius = geometry.radius(scene, placement.depth);
        let descendant_radii: Vec<f64> = track
            .tree
            .descendants(*node)
            .into_iter()
            .map(|descendant| geometry.radius(scene, scene.source_placements[descendant].depth))
            .collect();
        let far_radius = match geometry.direction {
            RadialDirection::Outward => descendant_radii
                .into_iter()
                .fold(start_radius + 2.0, f64::max),
            RadialDirection::Inward => descendant_radii
                .into_iter()
                .fold((start_radius - 2.0).max(0.0), f64::min),
        };
        let angle = geometry.angle(row as f64);
        let half = (geometry.angular_step() * 0.34).min(std::f64::consts::PI * 0.24);
        let (tip_x, tip_y) = geometry.point(start_radius, angle);
        let (left_x, left_y) = geometry.point(far_radius, angle - half);
        let (right_x, right_y) = geometry.point(far_radius, angle + half);
        let d = format!(
            "M {} {} L {} {} A {} {} 0 0 1 {} {} Z",
            num(tip_x),
            num(tip_y),
            num(left_x),
            num(left_y),
            num(far_radius),
            num(far_radius),
            num(right_x),
            num(right_y)
        );
        let title = format!(
            "{} ({} tips)",
            track.tree.nodes()[*node].name.as_deref().unwrap_or("clade"),
            track.tree.clade_size(*node)
        );
        ctx.svg.begin_titled(&title);
        ctx.svg.path(&d, &colors[*node], 0.28);
        ctx.svg.end_group();
    }
}

fn draw_trait_rings(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if track.trait_columns.is_empty() {
        return;
    }
    let gap = ctx.theme.tokens.legend_gap.clamp(1.0, 4.0);
    let mut inner = geometry.tree_outer + gap;
    for column in &track.trait_columns {
        let outer = (inner + column.ring_width).min(geometry.ring_outer);
        let values: Vec<Option<&AnnotationValue>> = scene
            .terminals
            .iter()
            .map(|node| inherited_annotation(&track.tree, *node, &column.key))
            .collect();
        let categories: BTreeMap<String, usize> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| inherited_annotation(&track.tree, placement.node, &column.key))
            .map(ToString::to_string)
            .fold(BTreeMap::new(), |mut categories, value| {
                let next = categories.len();
                categories.entry(value).or_insert(next);
                categories
            });
        let numeric: Vec<f64> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| inherited_annotation(&track.tree, placement.node, &column.key))
            .filter_map(AnnotationValue::as_number)
            .filter(|value| value.is_finite())
            .collect();
        let minimum = numeric.iter().copied().fold(f64::MAX, f64::min);
        let maximum = numeric.iter().copied().fold(f64::MIN, f64::max);
        for (row, node) in scene.terminals.iter().enumerate() {
            let angle = geometry.angle(row as f64);
            let gap_angle = if outer > 0.0 { 0.8 / outer } else { 0.0 };
            let half = (geometry.angular_step() / 2.0 - gap_angle)
                .max(geometry.angular_step() * 0.12)
                .min(std::f64::consts::PI * 0.45);
            let start = if geometry.full_circle() {
                angle - half
            } else {
                (angle - half).max(geometry.start)
            };
            let end = if geometry.full_circle() {
                angle + half
            } else {
                (angle + half).min(geometry.start + geometry.sweep)
            };
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
            let name = terminal_label(&track.tree, *node, &track.collapsed);
            let title = match &displayed {
                Some(value) => format!("{name}; {} {value}", column.key),
                None => format!("{name}; {} missing", column.key),
            };
            let path = radial_sector_path(geometry, inner, outer, start, end);
            ctx.svg.begin_titled(&title);
            if let Some(fill) = &fill {
                ctx.svg.path(&path, fill, 1.0);
            } else {
                ctx.svg
                    .path_stroked(&path, &ctx.theme.rule, ctx.theme.tokens.hairline);
            }
            if column.show_values {
                let text = displayed.as_deref().unwrap_or("—");
                let size = (ctx.theme.font_size - 3.0).max(6.0);
                let middle_radius = (inner + outer) / 2.0;
                let arc_room = middle_radius * (end - start).abs();
                if column.ring_width >= size + 1.0 && arc_room >= text_width(text, size) + 4.0 {
                    let (x, y) = geometry.point(middle_radius, angle);
                    let ink = match &fill {
                        Some(fill) => contrast_ink(fill),
                        None => ctx.theme.muted.as_str(),
                    };
                    ctx.svg.text_rotated(
                        (x, y + size * 0.3),
                        upright_tangent(angle),
                        text,
                        ink,
                        size,
                        crate::svg::Anchor::Middle,
                    );
                }
            }
            ctx.svg.end_group();
        }
        inner = outer + gap;
    }
}

fn radial_sector_path(
    geometry: &RadialGeometry,
    inner: f64,
    outer: f64,
    start: f64,
    end: f64,
) -> String {
    let (x0, y0) = geometry.point(outer, start);
    let (x1, y1) = geometry.point(outer, end);
    let (x2, y2) = geometry.point(inner, end);
    let (x3, y3) = geometry.point(inner, start);
    let large = usize::from((end - start).abs() > std::f64::consts::PI);
    format!(
        "M {} {} A {} {} 0 {} 1 {} {} L {} {} A {} {} 0 {} 0 {} {} Z",
        num(x0),
        num(y0),
        num(outer),
        num(outer),
        large,
        num(x1),
        num(y1),
        num(x2),
        num(y2),
        num(inner),
        num(inner),
        large,
        num(x3),
        num(y3)
    )
}

fn draw_trait_ring_headings(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    if track.trait_columns.is_empty() {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let slot = ctx.band.w / track.trait_columns.len() as f64;
    for (index, column) in track.trait_columns.iter().enumerate() {
        let x = ctx.band.x + index as f64 * slot;
        let visible = fit_text(&column.label, (slot - 18.0).max(1.0), size);
        if visible != column.label {
            ctx.svg.begin_titled(&column.label);
        }
        ctx.svg.circle_ringed(
            x + 7.0,
            ctx.band.y + 7.0,
            4.0,
            &ctx.theme.background,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );
        ctx.svg.text(
            x + 16.0,
            ctx.band.y + size + 2.0,
            &visible,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Start,
        );
        if visible != column.label {
            ctx.svg.end_group();
        }
    }
}

fn upright_tangent(angle: f64) -> f64 {
    let mut degrees = (angle.to_degrees() + 90.0).rem_euclid(360.0);
    if degrees > 90.0 && degrees < 270.0 {
        degrees += 180.0;
    }
    degrees
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
        area.x + self.fraction(value) * area.w
    }

    fn fraction(&self, value: f64) -> f64 {
        let span = self.maximum - self.minimum;
        if span <= 0.0 {
            return 0.0;
        }
        let fraction = match self.direction {
            TimeDirection::Increasing => (value - self.minimum) / span,
            TimeDirection::Decreasing if self.temporal => (self.maximum - value) / span,
            TimeDirection::Decreasing => (value - self.minimum) / span,
        };
        fraction.clamp(0.0, 1.0)
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
    let visible: Vec<usize> = scene
        .placements
        .iter()
        .flatten()
        .map(|placement| placement.node)
        .collect();
    let numeric: Vec<f64> = visible
        .iter()
        .filter_map(|node| values[*node].and_then(AnnotationValue::as_number))
        .collect();
    let all_numeric = numeric.len()
        == visible
            .iter()
            .filter(|node| values[**node].is_some())
            .count()
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
        let categories: BTreeMap<String, usize> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| inherited_annotation(tree, placement.node, &column.key))
            .map(ToString::to_string)
            .fold(BTreeMap::new(), |mut categories, value| {
                let next = categories.len();
                categories.entry(value).or_insert(next);
                categories
            });
        let numeric: Vec<f64> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| inherited_annotation(tree, placement.node, &column.key))
            .filter_map(AnnotationValue::as_number)
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
            match index {
                0 => crate::svg::Anchor::Start,
                2 => crate::svg::Anchor::End,
                _ => crate::svg::Anchor::Middle,
            },
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
        assert!(
            svg.contains("text-anchor=\"start\">2021 year</text>"),
            "{svg}"
        );
        assert!(
            svg.contains("text-anchor=\"end\">2025 year</text>"),
            "{svg}"
        );
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
    fn trait_categories_keep_branch_colours_after_ladderizing_and_collapsing() {
        let mut tree = Tree::parse_annotated_newick(
            "((A[&kind=alpha]:1,B[&kind=alpha]:1)alpha_clade[&kind=alpha]:1,C[&kind=beta]:2);",
        )
        .unwrap();
        let alpha = tree.node_named("alpha_clade").unwrap();
        tree.ladderize(false);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                TreeTrack::new(tree)
                    .color_by("kind")
                    .collapse(alpha)
                    .trait_categorical("kind"),
            )
            .to_svg();
        let beta = svg.find("<title>C; kind beta</title>").unwrap();
        assert!(
            svg[beta..(beta + 180).min(svg.len())].contains("fill=\"#d55e00\""),
            "{svg}"
        );
        let alpha = svg
            .find("<title>alpha_clade (2 tips); kind alpha</title>")
            .unwrap();
        assert!(
            svg[alpha..(alpha + 220).min(svg.len())].contains("fill=\"#0072b2\""),
            "{svg}"
        );
    }

    #[test]
    fn a_circular_tree_preserves_every_branch_and_draws_internal_arcs() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()).circular().show_tips(false))
            .to_svg();
        assert_eq!(svg.matches("<line").count(), 6, "one radial line per edge");
        assert_eq!(svg.matches("<path").count(), 3, "one arc per internal node");
        assert!(!svg.contains("NaN"), "{svg}");
    }

    #[test]
    fn circular_tip_labels_stay_upright_on_both_halves() {
        let svg = Figure::new(region())
            .width(520.0)
            .show_region_label(false)
            .push(TreeTrack::new(tree()).circular().radial_size(360.0))
            .to_svg();
        for tip in ["A", "B", "C", "D"] {
            assert!(svg.contains(&format!(">{tip}</text>")), "{svg}");
        }
        assert!(svg.contains("rotate(0)"), "right-facing label: {svg}");
        assert!(svg.contains("rotate(360)"), "left-facing label: {svg}");
    }

    #[test]
    fn a_fan_and_an_inward_tree_are_distinct_finite_projections() {
        let outward = Figure::new(region())
            .show_region_label(false)
            .push(TreeTrack::new(tree()).fan(220.0).show_tips(false))
            .to_svg();
        let inward = Figure::new(region())
            .show_region_label(false)
            .push(
                TreeTrack::new(tree())
                    .fan(220.0)
                    .radial_direction(RadialDirection::Inward)
                    .inner_radius(0.35)
                    .show_tips(false),
            )
            .to_svg();
        assert_ne!(outward, inward);
        assert!(!outward.contains("NaN"));
        assert!(!inward.contains("NaN"));
    }

    #[test]
    fn circular_time_guides_keep_their_exact_values() {
        let tree = Tree::parse_annotated_newick(
            "((A[&date=2024]:1,B[&date=2025]:2)AB:1,C[&date=2023]:3);",
        )
        .unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                TreeTrack::new(tree)
                    .circular()
                    .time("date")
                    .time_unit("year")
                    .show_tips(false),
            )
            .to_svg();
        for label in ["2021 year", "2023 year", "2025 year"] {
            assert!(svg.contains(&format!(">{label}</text>")), "{svg}");
        }
    }

    #[test]
    fn trait_columns_become_annotated_rings_in_circular_trees() {
        let tree = Tree::parse_annotated_newick(
            "(A[&country=Peru,coverage=18]:1,B[&country=Chile]:1,C[&country=Peru,coverage=42]:1);",
        )
        .unwrap();
        let svg = Figure::new(region())
            .width(560.0)
            .show_region_label(false)
            .push(
                TreeTrack::new(tree)
                    .circular()
                    .trait_column(
                        TraitColumn::categorical("country")
                            .label("Country")
                            .ring_width(12.0),
                    )
                    .trait_column(
                        TraitColumn::continuous("coverage")
                            .label("Depth")
                            .ring_width(12.0),
                    ),
            )
            .to_svg();
        for title in [
            "A; country Peru",
            "B; country Chile",
            "B; coverage missing",
            "C; coverage 42",
        ] {
            assert!(svg.contains(&format!("<title>{title}</title>")), "{svg}");
        }
        for heading in [">Country</text>", ">Depth</text>"] {
            assert!(svg.contains(heading), "{svg}");
        }
    }

    #[test]
    fn circular_collapse_is_a_non_destructive_wedge() {
        let tree = Tree::parse_newick("((A:1,B:1)outbreak:1,C:2);").unwrap();
        let outbreak = tree.node_named("outbreak").unwrap();
        let track = TreeTrack::new(tree)
            .circular()
            .collapse(outbreak)
            .show_tips(false);
        assert_eq!(track.tree().leaf_names(), ["A", "B", "C"]);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("<title>outbreak (2 tips)</title>"), "{svg}");
        assert!(svg.contains("fill-opacity=\"0.28\""), "{svg}");
    }

    #[test]
    fn radial_height_is_explicit_and_independent_of_leaf_count() {
        let scale = Scale::new(&region(), 0.0, 100.0);
        assert_eq!(
            TreeTrack::new(tree())
                .projection(TreeProjection::Circular)
                .radial_size(320.0)
                .height(&scale),
            320.0
        );
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
