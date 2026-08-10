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
//! # One function draws every rectangular tree in the crate
//!
//! [`draw_tree`] is a free function rather than a method on [`TreeTrack`]: the
//! standalone track, the tracks that carry a tree in a strip of their own and
//! both halves of a tanglegram all go through it. What it draws is rectangular
//! rather than diagonal, because a diagonal would imply the tree says something
//! about the space between two rows, and it says nothing about it.
//! A standalone [`TreeTrack`] can instead use [`TreeProjection::Circular`] or
//! [`TreeProjection::Unrooted`]. Circular coordinates retain the rooted depth;
//! unrooted coordinates choose a topology-balanced centre and do not privilege
//! the arbitrary root in the source Newick. Neither aligns to neighbouring rows.
//!
//! The tracks whose subject is the tree itself take the same drawing with its
//! branches named, so a clade can be pointed at for its support. A tree
//! standing beside a panel of rows does not, because the rows are named down
//! the side already and a title on every branch would be that same string a
//! second time. A tip is named on its branch only when its label is not drawn,
//! for exactly the same reason.

use std::collections::{BTreeMap, BTreeSet};

use crate::scale::Scale;
use crate::style::LinePattern;
use crate::svg::{finite_within, fit_text, num, text_rounded, text_width};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::{AnnotationValue, Placement, TimeDirection, Tree};

mod decorate;
mod interactions;
mod radial;
mod rectangular;
mod scale;
mod scene;
mod unrooted;

#[cfg(test)]
mod tests;

use self::decorate::*;
use self::interactions::*;
use self::radial::*;
use self::rectangular::*;
use self::scale::*;
use self::scene::*;
use self::unrooted::*;

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
    /// Topology drawn without assigning the source root a privileged position.
    Unrooted,
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

/// Mark used for one metadata dataset beside or around a tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TraitStyle {
    /// One filled cell or annular sector per terminal taxon.
    #[default]
    Strip,
    /// Numeric value encoded by bar length or radial height.
    Bar,
    /// Boolean or zero/non-zero value encoded by presence of a marker.
    Binary,
    /// Category encoded redundantly by both colour and marker shape.
    Symbol,
}

/// Visible encoding used for internal-node support values.
///
/// Support always remains available in exact SVG tooltips. This setting adds
/// marks or text when the values need to be readable without hovering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SupportStyle {
    /// Keep support in tooltips only.
    #[default]
    None,
    /// Scale an internal-node marker by support.
    Symbols,
    /// Print the original support value beside the node.
    Labels,
    /// Draw both the scaled marker and its value.
    SymbolsAndLabels,
}

/// Graphic attached directly to an annotated phylogenetic node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeGlyphStyle {
    /// One circle whose area follows a numeric annotation.
    Bubble,
    /// A composition drawn as sectors of a filled circle.
    Pie,
    /// A composition drawn as an annulus with a quiet centre.
    Donut,
    /// A composition drawn as one compact horizontal stacked bar.
    StackedBar,
}

/// Which annotated nodes receive a [`NodeGlyph`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NodeGlyphTarget {
    /// Draw every node carrying the requested numeric annotation data.
    #[default]
    All,
    /// Draw internal nodes only, including the root.
    Internal,
    /// Draw terminal taxa only.
    Leaves,
}

/// A branch-wise mixture of fitted omega rate classes and their weights.
///
/// Each rate key is paired with the weight key at the same index.  Values are
/// read directly from the child node that owns the incoming branch and are not
/// inherited.  The visible capsule normalises weights only for geometry; exact
/// supplied weights remain in the tooltip.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchRateMixture {
    rate_keys: Vec<String>,
    weight_keys: Vec<String>,
    label: String,
    width: f64,
    thickness: f64,
    neutral_lower: f64,
    neutral_upper: f64,
    saturation: f64,
}

impl BranchRateMixture {
    /// Pairs `rate_keys` and `weight_keys` in iterator order.
    pub fn new<R, W, RS, WS>(rate_keys: R, weight_keys: W) -> Self
    where
        R: IntoIterator<Item = RS>,
        W: IntoIterator<Item = WS>,
        RS: Into<String>,
        WS: Into<String>,
    {
        BranchRateMixture {
            rate_keys: rate_keys.into_iter().map(Into::into).collect(),
            weight_keys: weight_keys.into_iter().map(Into::into).collect(),
            label: "branch omega mixture".into(),
            width: 24.0,
            thickness: 5.2,
            neutral_lower: 0.95,
            neutral_upper: 1.05,
            saturation: 4.0,
        }
    }

    /// Replaces the visible legend label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Sets the preferred capsule length in pixels.
    pub fn width(mut self, width: f64) -> Self {
        if width.is_finite() {
            self.width = width.max(6.0);
        }
        self
    }

    /// Sets the capsule thickness in pixels.
    pub fn thickness(mut self, thickness: f64) -> Self {
        if thickness.is_finite() {
            self.thickness = thickness.clamp(2.0, 14.0);
        }
        self
    }

    /// Sets the inclusive rate interval rendered as approximately neutral.
    pub fn neutral_band(mut self, lower: f64, upper: f64) -> Self {
        if lower.is_finite() && upper.is_finite() && (0.0..=1.0).contains(&lower) && upper >= 1.0 {
            self.neutral_lower = lower;
            self.neutral_upper = upper;
        }
        self
    }

    /// Sets the positive omega value at which colours saturate.
    pub fn saturation(mut self, omega: f64) -> Self {
        if omega.is_finite() && omega > 1.0 {
            self.saturation = omega;
        }
        self
    }
}

/// Connections between branches carrying the same direct event annotation.
///
/// This is deliberately named for a visual hypothesis rather than a proof:
/// repeated ancestral-state reconstructions can represent convergence,
/// reversal or uncertainty.  The tooltip calls them recurrent events and
/// leaves that interpretation with the analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct HomoplasyLayer {
    key: String,
    label: String,
    minimum_occurrences: usize,
    maximum_connections: usize,
    width: f64,
}

impl HomoplasyLayer {
    /// Groups direct branch annotations stored under `key`.
    pub fn new(key: impl Into<String>) -> Self {
        let key = key.into();
        HomoplasyLayer {
            label: key.clone(),
            key,
            minimum_occurrences: 2,
            maximum_connections: 96,
            width: 1.15,
        }
    }

    /// Replaces the visible legend label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Requires at least this many branches before connecting an event.
    pub fn minimum_occurrences(mut self, minimum: usize) -> Self {
        self.minimum_occurrences = minimum.max(2);
        self
    }

    /// Caps the number of curves emitted by this layer.
    pub fn maximum_connections(mut self, maximum: usize) -> Self {
        self.maximum_connections = maximum.max(1);
        self
    }

    /// Sets the connection width in pixels.
    pub fn width(mut self, width: f64) -> Self {
        if width.is_finite() {
            self.width = width.clamp(0.4, 5.0);
        }
        self
    }
}

/// A data glyph placed on every matching annotated node.
///
/// Bubble glyphs read one numeric annotation. Pie, donut and stacked-bar
/// glyphs read one numeric annotation per supplied key, preserve key order and
/// normalise only the visible geometry; exact values remain in SVG tooltips.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeGlyph {
    keys: Vec<String>,
    label: String,
    style: NodeGlyphStyle,
    target: NodeGlyphTarget,
    size: f64,
    minimum_size: f64,
}

impl NodeGlyph {
    /// Scales circle area by numeric annotation `key`.
    pub fn bubble(key: impl Into<String>) -> Self {
        let key = key.into();
        NodeGlyph {
            label: key.clone(),
            keys: vec![key],
            style: NodeGlyphStyle::Bubble,
            target: NodeGlyphTarget::All,
            size: 9.0,
            minimum_size: 2.5,
        }
    }

    /// Draws a compositional pie from numeric annotation `keys`.
    pub fn pie<I, S>(keys: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::composition(keys, NodeGlyphStyle::Pie)
    }

    /// Draws a compositional donut from numeric annotation `keys`.
    pub fn donut<I, S>(keys: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::composition(keys, NodeGlyphStyle::Donut)
    }

    /// Draws a compact stacked bar from numeric annotation `keys`.
    pub fn stacked_bar<I, S>(keys: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self::composition(keys, NodeGlyphStyle::StackedBar)
    }

    fn composition<I, S>(keys: I, style: NodeGlyphStyle) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let keys: Vec<String> = keys.into_iter().map(Into::into).collect();
        NodeGlyph {
            label: keys.join(" / "),
            keys,
            style,
            target: NodeGlyphTarget::All,
            size: 9.0,
            minimum_size: 2.5,
        }
    }

    /// Replaces the visible legend label without changing annotation keys.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Restricts the glyph to all nodes, internal nodes or leaves.
    pub fn target(mut self, target: NodeGlyphTarget) -> Self {
        self.target = target;
        self
    }

    /// Sets the largest bubble radius or the nominal composition radius.
    pub fn size(mut self, pixels: f64) -> Self {
        self.size = finite_within(pixels, 2.0, 30.0, 9.0);
        self.minimum_size = self.minimum_size.min(self.size);
        self
    }

    /// Sets the smallest positive bubble radius.
    pub fn minimum_size(mut self, pixels: f64) -> Self {
        self.minimum_size = finite_within(pixels, 0.8, self.size, 2.5);
        self
    }

    /// Annotation keys read by this glyph, in visual order.
    pub fn keys(&self) -> &[String] {
        &self.keys
    }

    /// Visual form of this glyph.
    pub fn glyph_style(&self) -> NodeGlyphStyle {
        self.style
    }
}

/// A translucent field identifying one named or indexed clade.
#[derive(Debug, Clone, PartialEq)]
pub struct CladeHighlight {
    node: usize,
    label: Option<String>,
    color: Option<String>,
    opacity: f64,
}

impl CladeHighlight {
    /// Highlights the descendants of `node` without changing the tree.
    pub fn new(node: usize) -> Self {
        CladeHighlight {
            node,
            label: None,
            color: None,
            opacity: 0.12,
        }
    }

    /// Adds visible text to the highlighted field.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the field colour, defaulting to the categorical theme palette.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Sets fill opacity between 0.03 and 0.35.
    pub fn opacity(mut self, opacity: f64) -> Self {
        self.opacity = finite_within(opacity, 0.03, 0.35, 0.12);
        self
    }

    /// Index of the clade root.
    pub fn node(&self) -> usize {
        self.node
    }
}

impl SupportStyle {
    fn symbols(self) -> bool {
        matches!(self, Self::Symbols | Self::SymbolsAndLabels)
    }

    fn labels(self) -> bool {
        matches!(self, Self::Labels | Self::SymbolsAndLabels)
    }
}

/// One metadata column drawn beside the terminal taxa of a tree.
#[derive(Debug, Clone, PartialEq)]
pub struct TraitColumn {
    key: String,
    label: String,
    scale: TraitScale,
    style: TraitStyle,
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
            style: TraitStyle::Strip,
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

    /// Builds a numeric bar column or radial bar ring.
    pub fn bar(key: impl Into<String>) -> Self {
        let mut column = Self::continuous(key);
        column.style = TraitStyle::Bar;
        column.show_values = false;
        column
    }

    /// Builds a boolean presence/absence marker dataset.
    ///
    /// Boolean values and finite numbers are accepted; zero is absent and a
    /// non-zero number is present. Text is left missing rather than guessed.
    pub fn binary(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.style = TraitStyle::Binary;
        column.width = 28.0;
        column.show_values = false;
        column
    }

    /// Builds a categorical dataset encoded by colour and marker shape.
    pub fn symbol(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.style = TraitStyle::Symbol;
        column.width = 32.0;
        column.show_values = false;
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

    /// Replaces the visual mark while retaining the column's value mapping.
    pub fn style(mut self, style: TraitStyle) -> Self {
        self.style = style;
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

    /// The mark used in rectangular, radial and unrooted projections.
    pub fn trait_style(&self) -> TraitStyle {
        self.style
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
                Some(format!("clade support {}", text_rounded(support, 3)))
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
    dnds: Option<DnDsLayer>,
    rate_mixtures: Vec<BranchRateMixture>,
    homoplasy_layers: Vec<HomoplasyLayer>,
    collapsed: BTreeSet<usize>,
    show_nodes: bool,
    show_root: bool,
    support_style: SupportStyle,
    support_threshold: f64,
    branch_labels: Option<BranchLabels>,
    scale_bar: Option<ScaleBar>,
    trait_columns: Vec<TraitColumn>,
    node_glyphs: Vec<NodeGlyph>,
    clade_highlights: Vec<CladeHighlight>,
}

#[derive(Debug, Clone)]
struct BranchLabels {
    key: String,
    size: f64,
}

/// A branch-wise dN/dS encoding centred on the biologically meaningful
/// neutral ratio rather than on the observed minimum and maximum.
#[derive(Debug, Clone)]
struct DnDsLayer {
    key: String,
    label: String,
    neutral_lower: f64,
    neutral_upper: f64,
    saturation: f64,
    significance: Option<DnDsSignificance>,
}

#[derive(Debug, Clone)]
struct DnDsSignificance {
    key: String,
    maximum: f64,
}

impl DnDsLayer {
    fn new(key: impl Into<String>) -> Self {
        DnDsLayer {
            key: key.into(),
            label: "dN/dS (ω)".to_string(),
            neutral_lower: 0.95,
            neutral_upper: 1.05,
            saturation: 4.0,
            significance: None,
        }
    }
}

#[derive(Debug, Clone, Default)]
struct ScaleBar {
    length: Option<f64>,
    unit: Option<String>,
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
            dnds: None,
            rate_mixtures: Vec::new(),
            homoplasy_layers: Vec::new(),
            collapsed: BTreeSet::new(),
            show_nodes: false,
            show_root: false,
            support_style: SupportStyle::None,
            support_threshold: 0.0,
            branch_labels: None,
            scale_bar: None,
            trait_columns: Vec::new(),
            node_glyphs: Vec::new(),
            clade_highlights: Vec::new(),
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

    /// Chooses rectangular, circular or unrooted coordinates.
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

    /// Draws an equal-angle tree around a topology-balanced central node.
    ///
    /// The source root is not used as the centre. Branch lengths are retained
    /// for a phylogram and topology alone is used for a cladogram.
    pub fn unrooted(mut self) -> Self {
        self.projection = TreeProjection::Unrooted;
        self
    }

    /// Draws a circular fan covering `sweep_degrees` clockwise.
    pub fn fan(mut self, sweep_degrees: f64) -> Self {
        self.projection = TreeProjection::Circular;
        self.radial.sweep_degrees = finite_within(sweep_degrees, 10.0, 359.0, 240.0);
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
        self.radial.sweep_degrees = finite_within(degrees, 10.0, 360.0, 360.0);
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
        self.radial.inner_radius = finite_within(fraction, 0.0, 0.85, 0.08);
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

    /// Sets the requested height of an unrooted drawing in pixels.
    pub fn unrooted_size(self, size: f64) -> Self {
        self.radial_size(size)
    }

    /// Rotates the first equal-angle sector of an unrooted tree.
    pub fn unrooted_start(mut self, degrees: f64) -> Self {
        self.projection = TreeProjection::Unrooted;
        if degrees.is_finite() {
            self.radial.start_degrees = degrees;
        }
        self
    }

    /// Reorients the owned tree around internal `node` and marks the new root.
    ///
    /// An invalid index or sampled tip leaves the tree unchanged. Use
    /// [`Tree::reroot`](crate::Tree::reroot) directly when failure must be
    /// handled rather than represented as an unchanged builder.
    pub fn reroot(mut self, node: usize) -> Self {
        if self.tree.reroot(node) {
            self.show_root = true;
        }
        self
    }

    /// Reorients the owned tree around an internal node with this exact name.
    pub fn reroot_named(mut self, name: &str) -> Self {
        if let Some(node) = self.tree.node_named(name) {
            if self.tree.reroot(node) {
                self.show_root = true;
            }
        }
        self
    }

    /// Roots halfway along the edge leading to a monophyletic named outgroup.
    ///
    /// Missing, duplicate, internal or non-monophyletic names leave the tree
    /// unchanged. The new root is inserted without converting an outgroup tip
    /// into an internal node.
    pub fn reroot_outgroup<I, S>(mut self, names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut nodes = Vec::new();
        for name in names {
            let Some(node) = self.tree.node_named(name.as_ref()) else {
                return self;
            };
            nodes.push(node);
        }
        if self.tree.reroot_outgroup(&nodes).is_some() {
            self.show_root = true;
        }
        self
    }

    /// Roots the owned phylogram at the midpoint of its weighted tip diameter.
    ///
    /// Missing, negative or non-finite branch lengths leave the tree unchanged.
    pub fn reroot_midpoint(mut self) -> Self {
        if self.tree.reroot_midpoint().is_some() {
            self.show_root = true;
        }
        self
    }

    /// Draws or hides the selected root marker in rooted projections.
    pub fn show_root(mut self, show: bool) -> Self {
        self.show_root = show;
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
        self.dnds = None;
        self
    }

    /// Colours incoming branches by a direct dN/dS (ω) annotation.
    ///
    /// Values below one use the cool side of a colour-vision-safe diverging
    /// scale, values near one are neutral and values above one use the warm
    /// side. Unlike [`TreeTrack::color_by`], the annotation is never inherited:
    /// a missing branch estimate stays visibly missing. The legend and exact
    /// SVG tooltips describe the biological regimes without treating ω > 1 as
    /// proof of selection by itself.
    pub fn dnds(mut self, key: impl Into<String>) -> Self {
        self.dnds = Some(DnDsLayer::new(key));
        self.color_by = None;
        self
    }

    /// Replaces the visible label of the dN/dS legend.
    ///
    /// This has no effect until [`TreeTrack::dnds`] has selected an annotation.
    pub fn dnds_label(mut self, label: impl Into<String>) -> Self {
        if let Some(dnds) = &mut self.dnds {
            dnds.label = label.into();
        }
        self
    }

    /// Sets the inclusive interval treated as approximately neutral.
    ///
    /// Invalid, negative or reversed bounds leave the current interval
    /// unchanged. The default is `0.95..=1.05`.
    pub fn dnds_neutral_band(mut self, lower: f64, upper: f64) -> Self {
        if lower.is_finite() && upper.is_finite() && (0.0..=1.0).contains(&lower) && upper >= 1.0 {
            if let Some(dnds) = &mut self.dnds {
                dnds.neutral_lower = lower;
                dnds.neutral_upper = upper;
            }
        }
        self
    }

    /// Sets where each side of the logarithmic dN/dS colour scale saturates.
    ///
    /// `4.0`, the default, makes ω ≥ 4 and ω ≤ 1/4 use the strongest warm and
    /// cool colours. Values between them retain continuous differences.
    pub fn dnds_saturation(mut self, fold: f64) -> Self {
        if fold.is_finite() && fold > 1.0 {
            if let Some(dnds) = &mut self.dnds {
                dnds.saturation = fold;
            }
        }
        self
    }

    /// Emphasises branches whose direct test annotation is at most `maximum`.
    ///
    /// This is commonly a p-value or an adjusted p-value. It changes branch
    /// weight, not colour, so effect size (dN/dS) and evidence remain separate
    /// visual channels. Missing or non-numeric test values are not emphasised.
    pub fn dnds_significance(mut self, key: impl Into<String>, maximum: f64) -> Self {
        if maximum.is_finite() && maximum >= 0.0 {
            if let Some(dnds) = &mut self.dnds {
                dnds.significance = Some(DnDsSignificance {
                    key: key.into(),
                    maximum,
                });
            }
        }
        self
    }

    /// Adds a compact, weighted omega-class capsule to matching branches.
    ///
    /// This is useful for branch-site models such as aBSREL where one mean
    /// omega would erase the fitted episodic class. Missing or invalid class
    /// pairs are omitted and weights are normalised only in visible geometry.
    pub fn branch_rate_mixture(mut self, mixture: BranchRateMixture) -> Self {
        if !mixture.rate_keys.is_empty() && mixture.rate_keys.len() == mixture.weight_keys.len() {
            self.rate_mixtures.push(mixture);
        }
        self
    }

    /// Connects branches carrying the same direct event annotation.
    pub fn homoplasy_layer(mut self, layer: HomoplasyLayer) -> Self {
        self.homoplasy_layers.push(layer);
        self
    }

    /// Convenience form of [`TreeTrack::homoplasy_layer`].
    pub fn homoplasy(self, key: impl Into<String>) -> Self {
        self.homoplasy_layer(HomoplasyLayer::new(key))
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

    /// Chooses how internal-node support is made visible.
    ///
    /// Values in either the `0..=1` or `0..=100` convention are recognised.
    /// Their original representation is retained in labels and tooltips.
    pub fn support_style(mut self, style: SupportStyle) -> Self {
        self.support_style = style;
        self
    }

    /// Hides visible support below `minimum`.
    ///
    /// `0.8` and `80.0` both mean eighty percent. Non-finite values reset the
    /// threshold to zero.
    pub fn support_threshold(mut self, minimum: f64) -> Self {
        self.support_threshold = support_fraction(minimum).unwrap_or(0.0);
        self
    }

    /// Labels each incoming branch with its own annotation `key`.
    ///
    /// Unlike [`TreeTrack::color_by`], values are not inherited from ancestor
    /// nodes. This makes the method suitable for mutations, gains, losses and
    /// other events that belong to one branch. Long text is fitted to the
    /// available segment while the complete value remains in its tooltip.
    pub fn branch_labels(mut self, key: impl Into<String>) -> Self {
        self.branch_labels = Some(BranchLabels {
            key: key.into(),
            size: 8.0,
        });
        self
    }

    /// Sets the font size of labels created by [`TreeTrack::branch_labels`].
    pub fn branch_label_size(mut self, size: f64) -> Self {
        if let Some(labels) = &mut self.branch_labels {
            labels.size = finite_within(size, 5.0, 18.0, 8.0);
        }
        self
    }

    /// Adds an automatically sized branch-length scale bar to a phylogram.
    ///
    /// Cladograms and explicitly time-scaled trees omit it because their axes
    /// do not represent evolutionary branch length.
    pub fn scale_bar(mut self) -> Self {
        self.scale_bar.get_or_insert_with(ScaleBar::default);
        self
    }

    /// Draws or removes the branch-length scale bar.
    pub fn show_scale_bar(mut self, show: bool) -> Self {
        if show {
            self.scale_bar.get_or_insert_with(ScaleBar::default);
        } else {
            self.scale_bar = None;
        }
        self
    }

    /// Requests an exact scale-bar length in the tree's branch-length units.
    ///
    /// Values longer than the visible tree span are clamped to that span.
    /// Invalid values fall back to automatic sizing.
    pub fn scale_bar_length(mut self, length: f64) -> Self {
        let bar = self.scale_bar.get_or_insert_with(ScaleBar::default);
        bar.length = (length.is_finite() && length > 0.0).then_some(length);
        self
    }

    /// Adds a unit such as `substitutions/site` to the scale-bar label.
    pub fn scale_bar_unit(mut self, unit: impl Into<String>) -> Self {
        let bar = self.scale_bar.get_or_insert_with(ScaleBar::default);
        let unit = unit.into();
        bar.unit = (!unit.is_empty()).then_some(unit);
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

    /// Adds a numeric bar column or radial bar ring.
    pub fn trait_bar(self, key: impl Into<String>) -> Self {
        self.trait_column(TraitColumn::bar(key))
    }

    /// Adds a boolean presence/absence dataset.
    pub fn trait_binary(self, key: impl Into<String>) -> Self {
        self.trait_column(TraitColumn::binary(key))
    }

    /// Adds a categorical colour-and-shape dataset.
    pub fn trait_symbol(self, key: impl Into<String>) -> Self {
        self.trait_column(TraitColumn::symbol(key))
    }

    /// Adds a node-attached bubble, composition or compact stacked bar.
    pub fn node_glyph(mut self, glyph: NodeGlyph) -> Self {
        if !glyph.keys.is_empty() {
            self.node_glyphs.push(glyph);
        }
        self
    }

    /// Adds a translucent clade field behind branches and node graphics.
    pub fn clade_highlight(mut self, highlight: CladeHighlight) -> Self {
        if self.tree.nodes().get(highlight.node).is_some() {
            self.clade_highlights.push(highlight);
        }
        self
    }

    /// Highlights a clade by its exact internal or terminal name.
    pub fn highlight_named(mut self, name: &str) -> Self {
        if let Some(node) = self.tree.node_named(name) {
            self.clade_highlights.push(CladeHighlight::new(node));
        }
        self
    }

    /// The tree.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    fn branch_scale(&self) -> Option<&ScaleBar> {
        self.scale_bar
            .as_ref()
            .filter(|_| self.shape == TreeShape::Phylogram && self.time.is_none())
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
        let time = self
            .time
            .as_ref()
            .filter(|time| time.show_axis)
            .map_or(0.0, |_| theme.font_size + theme.tokens.tick_length + 5.0);
        let scale = self
            .branch_scale()
            .map_or(0.0, |_| theme.font_size + theme.tokens.tick_length + 7.0);
        time + scale
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

    fn annotation_header_room(&self) -> f64 {
        if self.trait_columns.is_empty()
            && self.node_glyphs.is_empty()
            && self.dnds.is_none()
            && self.rate_mixtures.is_empty()
            && self.homoplasy_layers.is_empty()
        {
            0.0
        } else {
            22.0
        }
    }

    fn rectangular_glyph_padding(&self) -> (f64, f64) {
        let (horizontal, vertical) =
            self.node_glyphs
                .iter()
                .fold(
                    (0.0f64, 0.0f64),
                    |(horizontal, vertical), glyph| match glyph.style {
                        NodeGlyphStyle::Bubble | NodeGlyphStyle::Pie | NodeGlyphStyle::Donut => {
                            (horizontal.max(glyph.size), vertical.max(glyph.size))
                        }
                        NodeGlyphStyle::StackedBar => (
                            horizontal.max(glyph.size * 1.5),
                            vertical.max(glyph.size * 0.39),
                        ),
                    },
                );
        (
            if horizontal > 0.0 {
                horizontal + 2.0
            } else {
                0.0
            },
            (vertical - self.row_height / 2.0).max(0.0),
        )
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
        let header_room = self.annotation_header_room();
        let (glyph_x, glyph_y) = self.rectangular_glyph_padding();
        let area = Rect {
            x: band.x + glyph_x,
            y: band.y + header_room + glyph_y,
            w: (band.w - tips - traits - glyph_x * 2.0).max(1.0),
            h: (band.h - axis_room - header_room - glyph_y * 2.0).max(1.0),
        };

        draw_rectangular_clade_highlights(self, ctx, &scene, area);
        draw_tree_scene(
            ctx,
            &self.tree,
            &scene,
            area,
            self.row_height,
            &color,
            self.line_width,
            self.color_by.as_deref(),
            self.dnds.as_ref(),
            self.show_nodes,
            self.support_style,
            self.support_threshold,
            self.branch_labels.as_ref(),
            &self.rate_mixtures,
            &self.homoplasy_layers,
            !self.show_tips,
        );
        draw_rectangular_node_glyphs(self, ctx, &scene, area);
        if self.show_root {
            if let Some(root) = scene.placements[self.tree.root()] {
                draw_root_marker(
                    ctx,
                    scene.x(area, root.depth) + ctx.theme.tokens.marker_radius * 1.4,
                    area.y + self.row_height / 2.0 + root.row * self.row_height,
                );
            }
        }

        if self.show_tips {
            let size = ctx.theme.font_size - 1.0;
            for (row, node) in scene.terminals.iter().enumerate() {
                let name = terminal_label(&self.tree, *node, &self.collapsed);
                ctx.svg.text(
                    area.right() + glyph_x + 4.0,
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
            tips + glyph_x,
            &self.trait_columns,
            self.row_height,
        );
        if let Some(time) = self.time.as_ref().filter(|time| time.show_axis) {
            draw_time_axis(ctx, &scene, area, time);
        }
        if let Some(bar) = self.branch_scale() {
            draw_rectangular_scale_bar(ctx, &scene, area, bar);
        }
        draw_annotation_legend(self, ctx);
    }
}

impl Track for TreeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        match self.projection {
            TreeProjection::Rectangular => {
                let rows = visible_terminals(&self.tree, &self.collapsed).len().max(1) as f64;
                let (_, glyph_y) = self.rectangular_glyph_padding();
                rows * self.row_height
                    + glyph_y * 2.0
                    + if self.time.as_ref().is_some_and(|time| time.show_axis) {
                        22.0
                    } else {
                        0.0
                    }
                    + if self.branch_scale().is_some() {
                        22.0
                    } else {
                        0.0
                    }
                    + self.annotation_header_room()
            }
            TreeProjection::Circular => self.radial.size + self.annotation_header_room(),
            TreeProjection::Unrooted => self.radial.size + self.annotation_header_room(),
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
            TreeProjection::Unrooted => draw_unrooted_track(self, ctx),
        }
    }
}
