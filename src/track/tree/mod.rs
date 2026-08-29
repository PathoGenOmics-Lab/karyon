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
use crate::track::traits::{binary_state, draw_column, TraitDomain, TraitRow};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::{AnnotationValue, Placement, TimeDirection, Tree};

// The metadata columns beside a phylogeny are the same columns a matrix or an
// alignment puts beside its rows, so they live in one module and are named
// from here for everything that already reaches them through this one.
pub use crate::track::traits::{TraitColumn, TraitScale, TraitStyle};

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

/// Shape of branches in the rectangular projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BranchGeometry {
    /// Horizontal child branches joined by vertical parent risers.
    #[default]
    Orthogonal,
    /// One straight segment from each parent node to each child node.
    Diagonal,
    /// Smooth parent-to-child curves with horizontal tangents at both ends.
    Curved,
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

/// One or more genomic or amino-acid events marked on their owning branch.
///
/// Text, numbers and booleans become one event. A brace-delimited annotated
/// Newick list becomes several ordered event symbols. Values are direct branch
/// data and are never inherited.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchEventLayer {
    key: String,
    label: String,
    maximum_events: usize,
    size: f64,
}

impl BranchEventLayer {
    /// Reads direct events stored under annotation `key`.
    pub fn new(key: impl Into<String>) -> Self {
        let key = key.into();
        BranchEventLayer {
            label: key.clone(),
            key,
            maximum_events: 8,
            size: 3.0,
        }
    }

    /// Replaces the visible legend label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Caps the number of event marks on one branch.
    pub fn maximum_events(mut self, maximum: usize) -> Self {
        self.maximum_events = maximum.max(1);
        self
    }

    /// Sets the event-symbol radius in pixels.
    pub fn size(mut self, pixels: f64) -> Self {
        self.size = finite_within(pixels, 1.4, 8.0, 3.0);
        self
    }
}

/// A branch estimate with lower and upper bounds.
///
/// The estimate is a point on a compact branch-aligned axis and the interval
/// is a whisker. This can carry concordance factors, ancestral transition
/// support, rate uncertainty or any upstream statistic with a meaningful
/// fixed domain.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchIntervalLayer {
    estimate_key: String,
    lower_key: String,
    upper_key: String,
    label: String,
    minimum: f64,
    maximum: f64,
    threshold: Option<f64>,
    width: f64,
}

impl BranchIntervalLayer {
    /// Reads a point estimate and its lower and upper bounds.
    pub fn new(
        estimate_key: impl Into<String>,
        lower_key: impl Into<String>,
        upper_key: impl Into<String>,
    ) -> Self {
        let estimate_key = estimate_key.into();
        BranchIntervalLayer {
            label: estimate_key.clone(),
            estimate_key,
            lower_key: lower_key.into(),
            upper_key: upper_key.into(),
            minimum: 0.0,
            maximum: 1.0,
            threshold: None,
            width: 27.0,
        }
    }

    /// Replaces the visible legend label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Sets the fixed numeric domain used by every interval in the layer.
    pub fn range(mut self, minimum: f64, maximum: f64) -> Self {
        if minimum.is_finite() && maximum.is_finite() && maximum > minimum {
            self.minimum = minimum;
            self.maximum = maximum;
        }
        self
    }

    /// Emphasises point estimates at or above `threshold`.
    pub fn threshold(mut self, threshold: f64) -> Self {
        self.threshold = threshold.is_finite().then_some(threshold);
        self
    }

    /// Sets the preferred compact-axis length in pixels.
    pub fn width(mut self, pixels: f64) -> Self {
        self.width = finite_within(pixels, 8.0, 80.0, 27.0);
        self
    }
}

/// Posterior probabilities for alternative ancestral states.
///
/// One probability key is supplied per state. Internal nodes receive donut
/// glyphs and, optionally, a branch marker when the maximum-posterior state
/// changes between a parent and child above the confidence threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct AncestralStateLayer {
    keys: Vec<String>,
    label: String,
    confidence: f64,
    size: f64,
    show_transitions: bool,
}

impl AncestralStateLayer {
    /// Uses numeric annotation `keys` as an ordered state composition.
    pub fn new<I, S>(keys: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        AncestralStateLayer {
            keys: keys.into_iter().map(Into::into).collect(),
            label: "ancestral state posterior".into(),
            confidence: 0.70,
            size: 8.0,
            show_transitions: true,
        }
    }

    /// Replaces the visible legend label.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Requires both endpoints to exceed this posterior for a transition mark.
    pub fn confidence(mut self, confidence: f64) -> Self {
        if confidence.is_finite() {
            self.confidence = confidence.clamp(0.0, 1.0);
        }
        self
    }

    /// Sets the internal-node donut radius in pixels.
    pub fn size(mut self, pixels: f64) -> Self {
        self.size = finite_within(pixels, 2.0, 30.0, 8.0);
        self
    }

    /// Shows or hides parent-to-child maximum-posterior state changes.
    pub fn show_transitions(mut self, show: bool) -> Self {
        self.show_transitions = show;
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
    branch_geometry: BranchGeometry,
    radial: RadialLayout,
    color: Option<String>,
    line_width: f64,
    show_tips: bool,
    time: Option<TimeAxis>,
    color_by: Option<String>,
    dnds: Option<DnDsLayer>,
    rate_mixtures: Vec<BranchRateMixture>,
    homoplasy_layers: Vec<HomoplasyLayer>,
    branch_event_layers: Vec<BranchEventLayer>,
    branch_interval_layers: Vec<BranchIntervalLayer>,
    ancestral_state_layers: Vec<AncestralStateLayer>,
    collapsed: BTreeSet<usize>,
    max_rows: Option<usize>,
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
    /// The diameter asked for, or `None` to work one out from the tips.
    size: Option<f64>,
}

impl Default for RadialLayout {
    fn default() -> Self {
        RadialLayout {
            start_degrees: -90.0,
            sweep_degrees: 360.0,
            direction: RadialDirection::Outward,
            inner_radius: 0.08,
            size: None,
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

/// The disc a radial or unrooted tree is drawn on when nothing asks for a size
/// and its tips fit. Trees larger than this grow past it; smaller ones stay
/// here, because the number this replaced was a fixed 440 and a figure that
/// suddenly shrinks reads as a bug even when the labels still fit.
const RADIAL_DIAMETER: f64 = 440.0;

impl TreeTrack {
    /// A track drawing `tree`.
    pub fn new(tree: Tree) -> Self {
        TreeTrack {
            tree,
            label: None,
            row_height: 15.0,
            shape: TreeShape::Phylogram,
            projection: TreeProjection::Rectangular,
            branch_geometry: BranchGeometry::Orthogonal,
            radial: RadialLayout::default(),
            color: None,
            line_width: 1.2,
            show_tips: true,
            time: None,
            color_by: None,
            dnds: None,
            rate_mixtures: Vec::new(),
            homoplasy_layers: Vec::new(),
            branch_event_layers: Vec::new(),
            branch_interval_layers: Vec::new(),
            ancestral_state_layers: Vec::new(),
            collapsed: BTreeSet::new(),
            max_rows: None,
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

    /// Chooses orthogonal, diagonal or curved rectangular branches.
    ///
    /// Circular and unrooted projections retain their own geometry.
    pub fn branch_geometry(mut self, geometry: BranchGeometry) -> Self {
        self.branch_geometry = geometry;
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
        self.radial.size = size.is_finite().then(|| size.max(120.0));
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
            self.fold_to_fit();
        }
        self
    }

    /// Reorients the owned tree around an internal node with this exact name.
    pub fn reroot_named(mut self, name: &str) -> Self {
        if let Some(node) = self.tree.node_named(name) {
            if self.tree.reroot(node) {
                self.show_root = true;
                self.fold_to_fit();
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
            self.fold_to_fit();
        }
        self
    }

    /// Roots the owned phylogram at the midpoint of its weighted tip diameter.
    ///
    /// Missing, negative or non-finite branch lengths leave the tree unchanged.
    pub fn reroot_midpoint(mut self) -> Self {
        if self.tree.reroot_midpoint().is_some() {
            self.show_root = true;
            self.fold_to_fit();
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

    /// Adds ordered genomic or amino-acid event symbols to matching branches.
    pub fn branch_event_layer(mut self, layer: BranchEventLayer) -> Self {
        self.branch_event_layers.push(layer);
        self
    }

    /// Convenience form of [`TreeTrack::branch_event_layer`].
    pub fn branch_events(self, key: impl Into<String>) -> Self {
        self.branch_event_layer(BranchEventLayer::new(key))
    }

    /// Adds a branch-aligned point estimate and uncertainty interval.
    pub fn branch_interval(mut self, layer: BranchIntervalLayer) -> Self {
        self.branch_interval_layers.push(layer);
        self
    }

    /// Adds internal ancestral-state posteriors and optional transitions.
    pub fn ancestral_states(mut self, layer: AncestralStateLayer) -> Self {
        if !layer.keys.is_empty() {
            self.node_glyphs.push(
                NodeGlyph::donut(layer.keys.clone())
                    .label(layer.label.clone())
                    .target(NodeGlyphTarget::Internal)
                    .size(layer.size),
            );
            self.ancestral_state_layers.push(layer);
        }
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

    /// Caps how many rows the tree draws, by collapsing clades until it fits.
    ///
    /// A phylogeny is the one track here that laid a row per tip and never
    /// stopped. Sixty thousand tips drew a figure nine hundred thousand pixels
    /// tall, and there was no way to ask for less: `row_height` floors at two,
    /// so twenty thousand tips could not be brought under forty thousand
    /// pixels by any setting.
    ///
    /// Nothing is dropped. A pileup that meets its cap stops opening rows and
    /// counts the reads it left out; a tree cannot, because a tip is not
    /// interchangeable with the tip below it and cutting the list would cut a
    /// clade in half. So it collapses instead, and every tip is inside a
    /// triangle that says how many it holds.
    ///
    /// Smallest first, so the shape survives: collapsing a cherry costs one
    /// row and hides two names, and collapsing near the root costs nothing and
    /// hides the tree. `None` lifts the cap, which is the default, because a
    /// tree of three hundred tips is an ordinary figure and capping it by
    /// default would fold figures nobody asked to fold.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|rows| rows.max(1));
        self.fold_to_fit();
        self
    }

    /// Collapses the smallest clades until the visible terminals fit the cap.
    ///
    /// Called again after every rerooting, because rerooting moves the tips
    /// about and a fold worked out against the old shape would collapse the
    /// wrong clades.
    fn fold_to_fit(&mut self) {
        let Some(cap) = self.max_rows else {
            return;
        };
        let nodes = self.tree.nodes();
        // How many rows each node contributes as things stand, which is one
        // for a leaf or an already collapsed clade and the sum of its children
        // otherwise. Kept as we go, so collapsing a clade whose own children
        // were folded does not count their rows twice.
        let order = postorder_nodes(&self.tree);
        let mut rows = vec![1usize; nodes.len()];
        // Tips below each node, which is what the clades are ranked by. It
        // comes off the same walk because asking the tree for it one node at a
        // time does not: `clade_size` collects the whole subtree into a vector
        // and counts the leaves in it, so ranking every clade that way is a
        // full traversal per clade, and the sort asks more than once each. A
        // hundred thousand tip tree spent 785 ms getting to sixty rows, more
        // than the 69 ms it took to draw all hundred thousand uncapped.
        let mut tips = vec![1usize; nodes.len()];
        for node in &order {
            if nodes[*node].is_leaf() {
                continue;
            }
            tips[*node] = nodes[*node].children.iter().map(|child| tips[*child]).sum();
            if self.collapsed.contains(node) {
                continue;
            }
            rows[*node] = nodes[*node].children.iter().map(|child| rows[*child]).sum();
        }
        let root = self.tree.root();
        let mut total = rows[root];
        if total <= cap {
            return;
        }

        // Smallest clade first, and the index breaks a tie, so the same tree
        // folds the same way every time. Going up in size also means a node's
        // ancestors are always still open when it is reached, so no collapse
        // here can sit inside another.
        let mut candidates: Vec<usize> = (0..nodes.len())
            .filter(|node| !nodes[*node].is_leaf() && *node != root)
            .collect();
        candidates.sort_by_key(|node| (tips[*node], *node));

        // Folds are recorded in a flat vector and only the outermost ones are
        // kept. Folding smallest first means a clade is often folded and then
        // swallowed by an ancestor folded later, and those inner folds are
        // inside something already collapsed, so nothing ever draws them. A
        // million tip tree made a million of them and spent 131 of its 168 ms
        // putting them into an ordered set, to keep about sixty that matter.
        let mut folded = vec![false; nodes.len()];
        for node in candidates {
            if total <= cap {
                break;
            }
            if self.collapsed.contains(&node) {
                continue;
            }
            // Smallest first means every child of this clade has already been
            // through the loop, so its row count is final and this one can be
            // added up here. The rows a fold saves used to be walked off every
            // node above it instead, which is the depth of the tree per fold
            // and was 131 of the 168 ms a million tip tree spent folding.
            let current: usize = nodes[node].children.iter().map(|child| rows[*child]).sum();
            rows[node] = current;
            if current <= 1 {
                continue;
            }
            folded[node] = true;
            rows[node] = 1;
            total -= current - 1;
        }

        // Parents before children, so a fold is kept only when nothing above
        // it was folded too.
        let mut inside = vec![false; nodes.len()];
        for node in order.iter().rev() {
            let above = nodes[*node].parent.is_some_and(|parent| {
                inside[parent] || folded[parent] || self.collapsed.contains(&parent)
            });
            inside[*node] = above;
            if folded[*node] && !above {
                self.collapsed.insert(*node);
            }
        }
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
    /// How much of the band the tip names need.
    ///
    /// The size here and the size the names are drawn at have to be the one
    /// number, or the gutter is held open for text that is not that big. They
    /// were both a flat `font_size - 1.0`, so a row two pixels tall carried an
    /// eleven pixel name: five rows of them through each other, and the last
    /// few sliced off by the track's own clip. Seven tracks here already clamp
    /// a row's text to the row, and the help for `--row-height` already says a
    /// row too short for a name shrinks the name with it.
    fn tip_size(&self, theme: &Theme) -> f64 {
        (theme.font_size - 1.0).min(self.row_height)
    }

    fn tip_width(&self, theme: &Theme, scene: &TreeScene) -> f64 {
        if !self.show_tips {
            return 0.0;
        }
        let size = self.tip_size(theme);
        scene
            .terminals
            .iter()
            .map(|node| text_width(&terminal_label(&self.tree, *node, &self.collapsed), size))
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
            && self.branch_event_layers.is_empty()
            && self.branch_interval_layers.is_empty()
            && self.ancestral_state_layers.is_empty()
        {
            0.0
        } else {
            22.0
        }
    }

    /// How wide across the circular and unrooted projections are drawn.
    ///
    /// A rectangular tree is as tall as its rows, so it grows with its data.
    /// The disc did not: it was a flat 440 pixels whatever the tree, which is
    /// generous for thirty tips and a solid band at six hundred, and the one
    /// number that decided whether it could be read was not reachable from a
    /// command line at all.
    ///
    /// It sizes itself now. The tips sit at equal angles on a circle, so the
    /// space between two neighbouring labels is the circumference at the label
    /// radius divided by the number of them, and measured against the real
    /// drawing that comes to `2.75 * diameter / tips` on a full turn. Turning
    /// that round gives the diameter at which a label of the theme's size
    /// still clears its neighbour. A fan draws the same tips on a fraction of
    /// the turn and needs the radius back in proportion.
    ///
    /// The figure's width is the ceiling, because the disc is drawn inside
    /// `min(width, height)` and a taller band would only add white above and
    /// below a circle the width already decided. A tree with more tips than
    /// that width can separate is drawn dense rather than drawn wrong, and the
    /// two ways out are the reader's: a wider figure, or `--max-rows`.
    fn radial_diameter(&self, scale: &Scale, theme: &Theme) -> f64 {
        if let Some(size) = self.radial.size {
            return size;
        }
        let terminals = visible_terminals(&self.tree, &self.collapsed);
        let tips = terminals.len().max(1) as f64;
        let size = (theme.font_size - 1.0).min(self.row_height.max(1.0));

        // The room the labels take out of the radius, which is what a fitted
        // proportion gets wrong at the small end: it is a fixed number of
        // pixels, so on a small disc it is most of the radius and on a large
        // one it is a rim.
        let extent = if self.show_tips {
            terminals
                .iter()
                .map(|node| text_width(&terminal_label(&self.tree, *node, &self.collapsed), size))
                .fold(0.0f64, f64::max)
                + 6.0
        } else {
            4.0
        };

        // Two neighbours on a circle of radius r, `turn / tips` of a turn
        // apart, are `2 r sin(pi turn / tips)` from each other. Ask for that to
        // be the height of a label and read the radius back out.
        let turn = (self.radial.sweep_degrees.abs() / 360.0).clamp(0.05, 1.0);
        let step = (std::f64::consts::PI * turn / tips).max(1e-6);
        let radius = size / (2.0 * step.sin());
        let wanted = 2.0 * (radius + extent);
        // Only ever larger. This rule exists so a tree with more tips than the
        // old fixed disc could hold gets the room it needs, and a tree with
        // three tips needs less room than that but does not want less: shrunk
        // to what its labels strictly require, its branches got short enough
        // to start ellipsising the labels drawn along them.
        // The ceiling is the figure's own inner width, because the drawing
        // fits the disc into the shorter of the band's two sides and height
        // beyond that is whitespace. It is spelt `x0 + width` rather than
        // `width` because a track stacked underneath can widen the gutter and
        // narrow every band: measured on the same fourteen tip tree, `width`
        // read 866 alone and 808.269 with company while the sum read 882 both
        // times, and a figure that gets shorter when you add a track to it is
        // a track reading its neighbour's data.
        let figure = scale.x0() + scale.width();
        wanted.clamp(RADIAL_DIAMETER, figure.max(RADIAL_DIAMETER))
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
            &self.branch_event_layers,
            &self.branch_interval_layers,
            &self.ancestral_state_layers,
            self.branch_geometry,
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
            let size = self.tip_size(ctx.theme);
            let names_at = area.right() + glyph_x + 4.0;
            for (row, node) in scene.terminals.iter().enumerate() {
                let middle = area.y + self.row_height / 2.0 + row as f64 * self.row_height;
                // A leader from the branch to the name it belongs to. Without
                // one the names all sit flush at the right while the branches
                // end wherever their lengths put them, so on a phylogram the
                // reader has to guess which name is which: measured on forty
                // tips, the median gap was 649 pixels and the widest 840, the
                // whole of the band a branch can occupy.
                //
                // The two other projections have drawn one all along. This is
                // the same hairline and the same threshold: a gap under half a
                // pixel is not a line, it is a smudge on the end of a branch.
                if let Some(placement) = scene.placements[*node] {
                    let ends = scene.x(area, placement.depth);
                    if names_at - ends > 0.5 {
                        ctx.svg.line(
                            ends,
                            middle,
                            names_at,
                            middle,
                            &ctx.theme.rule,
                            ctx.theme.tokens.hairline,
                        );
                    }
                }
                let name = terminal_label(&self.tree, *node, &self.collapsed);
                ctx.svg.text(
                    names_at,
                    middle + size * 0.35,
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
    fn height(&self, scale: &Scale) -> f64 {
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
            TreeProjection::Circular | TreeProjection::Unrooted => {
                self.radial_diameter(scale, &Theme::default()) + self.annotation_header_room()
            }
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
