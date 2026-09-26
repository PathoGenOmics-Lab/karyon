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
use std::sync::OnceLock;

use crate::scale::Scale;
use crate::style::LinePattern;
use crate::svg::{finite_within, fit_text, num, text_rounded, text_width};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::traits::{
    binary_state, draw_column, Dealt, Join, Stretch, TraitDomain, TraitRow, Traits,
};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::{AnnotationValue, NodeRef, Placement, TimeDirection, Tree};

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
pub(crate) use unrooted::unrooted_layout;

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
    /// Keep support in tooltips only: the default.
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
/// A list annotation is read one event at a time, as [`BranchEventLayer`]
/// reads it, so a branch carrying `{S45N,E88K}` is joined to every other
/// branch carrying either change.
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
    /// The clade as it was named, found in the tree when the highlight is
    /// handed to a track.
    wanted: NodeRef,
    node: usize,
    label: Option<String>,
    color: Option<String>,
    opacity: f64,
}

impl CladeHighlight {
    /// Highlights the descendants of a clade without changing the tree: an
    /// index, a name, or a [`NodeRef`] picking it by its tips or by a value.
    pub fn new(node: impl Into<NodeRef>) -> Self {
        let wanted = node.into();
        let node = match wanted {
            NodeRef::Index(node) => node,
            _ => usize::MAX,
        };
        CladeHighlight {
            wanted,
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

/// The tree to draw beside rows: cut to the tips whose row is among the
/// first `drawn` of `rows`, and with how many of its tips have no row at all.
///
/// Drawn whole, a tree with a tip the panel lacks put every row after that
/// tip beside the branch of the tip before it, and the last branch below the
/// last row: the rows were put in the tree's order, but the tree was never
/// cut to the rows. A row hidden under a cap goes from the tree too, so no
/// branch leads off the band to a row that is not there. `None` when no tip
/// has a row that is drawn.
pub(crate) fn tree_beside_rows<'a>(
    tree: &'a Tree,
    rows: &[String],
    drawn: usize,
) -> (Option<std::borrow::Cow<'a, Tree>>, usize) {
    use std::collections::HashSet;
    let all: HashSet<&str> = rows.iter().map(String::as_str).collect();
    let shown: HashSet<&str> = rows.iter().take(drawn).map(String::as_str).collect();
    let leaves = tree.leaf_names();
    let without_row = leaves
        .iter()
        .filter(|leaf| !all.contains(leaf.as_str()))
        .count();
    let kept: Vec<&str> = leaves
        .iter()
        .map(String::as_str)
        .filter(|leaf| shown.contains(leaf))
        .collect();
    if kept.len() == leaves.len() {
        return (Some(std::borrow::Cow::Borrowed(tree)), without_row);
    }
    (
        tree.keep_tips(kept).map(std::borrow::Cow::Owned),
        without_row,
    )
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

/// Names as a sentence lists them: the first three, and how many more.
fn listed(names: &[String]) -> String {
    match names.len() {
        0..=3 => names.join(", "),
        n => format!("{} and {} more", names[..3].join(", "), n - 3),
    }
}

/// The levels of `key` over `tree`, each with the palette colour a
/// [`TreeTrack`] colouring its branches by `key` deals it: in the order the
/// tree meets them, from the palette's first colour.
pub(crate) fn tree_levels(tree: &Tree, key: &str) -> BTreeMap<String, usize> {
    rectangular::tree_domain(tree, key, Dealt::default()).categories
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
///
/// The settings can be written in any order. One written before the method
/// that turns its layer on, such as [`TreeTrack::time_unit`] before
/// [`TreeTrack::time`], counts the same as one written after it.
#[derive(Debug, Clone)]
pub struct TreeTrack {
    // Every setting is kept in a field of its own, written by its own method
    // and by nothing else, and the time axis, the dN/dS colouring, the branch
    // labels and the folds are put together from them when the tree is drawn.
    // They used to be put together when the method that turns them on was
    // called, which started them afresh: a unit written before `time`, a
    // label size written before `branch_labels` or a hidden root written
    // before a reroot was dropped without a word, and the same settings drew
    // two different figures depending on the order they were written in.
    tree: Tree,
    label: Option<String>,
    row_height: f64,
    shape: TreeShape,
    projection: TreeProjection,
    /// Whether a projection was chosen by name, which a radial setting then
    /// leaves as it is rather than turning the tree into a circle.
    projection_chosen: bool,
    branch_geometry: BranchGeometry,
    radial: RadialLayout,
    color: Option<String>,
    line_width: f64,
    show_tips: bool,
    time: Option<String>,
    time_direction: TimeDirection,
    time_unit: Option<String>,
    show_time_axis: bool,
    color_by: Option<String>,
    dnds: Option<String>,
    dnds_label: String,
    dnds_neutral_band: (f64, f64),
    dnds_saturation: f64,
    dnds_significance: Option<DnDsSignificance>,
    rate_mixtures: Vec<BranchRateMixture>,
    homoplasy_layers: Vec<HomoplasyLayer>,
    branch_event_layers: Vec<BranchEventLayer>,
    branch_interval_layers: Vec<BranchIntervalLayer>,
    ancestral_state_layers: Vec<AncestralStateLayer>,
    /// The clades folded by hand, with [`TreeTrack::collapse`].
    collapsed: BTreeSet<usize>,
    /// What a clade folded as the clade of a value is called: the value.
    fold_names: BTreeMap<usize, String>,
    max_rows: Option<usize>,
    /// Every clade drawn folded: the ones folded by hand and the ones the row
    /// cap folds on top of them. Worked out the first time the tree is drawn,
    /// against the tree as it is by then, and emptied by every method that
    /// changes what it depends on.
    folds: OnceLock<BTreeSet<usize>>,
    show_nodes: bool,
    /// `None` until [`TreeTrack::show_root`] is written, and a reroot marks
    /// the root it chose only while it is.
    show_root: Option<bool>,
    rerooted: bool,
    support_style: SupportStyle,
    support_threshold: f64,
    branch_labels: Option<String>,
    branch_label_size: f64,
    scale_bar: ScaleBar,
    show_scale_bar: bool,
    trait_columns: Vec<TraitColumn>,
    node_glyphs: Vec<NodeGlyph>,
    clade_highlights: Vec<CladeHighlight>,
    /// Requests a builder could not carry out, each said in a line.
    refused: Vec<String>,
    /// What the sheet given to [`TreeTrack::traits`] matched and left out.
    joined: Option<Join>,
    /// How a strip of each key of the joined sheet would deal the palette,
    /// which the branches take for that key with or without the strip.
    sheet_dealing: BTreeMap<String, (Vec<String>, Stretch)>,
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

/// How strongly a folded clade's wedge is filled with the colour of its
/// branches. At 0.28 a pale L4 wedge read as a sliver of nothing, and two
/// readers asked what it was.
const FOLD_FILL: f64 = 0.5;

impl TreeTrack {
    /// A track drawing `tree`.
    pub fn new(tree: Tree) -> Self {
        TreeTrack {
            tree,
            label: None,
            row_height: 15.0,
            shape: TreeShape::Phylogram,
            projection: TreeProjection::Rectangular,
            projection_chosen: false,
            branch_geometry: BranchGeometry::Orthogonal,
            radial: RadialLayout::default(),
            color: None,
            line_width: 1.2,
            show_tips: true,
            time: None,
            time_direction: TimeDirection::Increasing,
            time_unit: None,
            show_time_axis: true,
            color_by: None,
            dnds: None,
            dnds_label: "dN/dS (ω)".to_string(),
            dnds_neutral_band: (0.95, 1.05),
            dnds_saturation: 4.0,
            dnds_significance: None,
            rate_mixtures: Vec::new(),
            homoplasy_layers: Vec::new(),
            branch_event_layers: Vec::new(),
            branch_interval_layers: Vec::new(),
            ancestral_state_layers: Vec::new(),
            collapsed: BTreeSet::new(),
            fold_names: BTreeMap::new(),
            max_rows: None,
            folds: OnceLock::new(),
            show_nodes: false,
            show_root: None,
            rerooted: false,
            support_style: SupportStyle::None,
            support_threshold: 0.0,
            branch_labels: None,
            branch_label_size: 8.0,
            // On by default: a phylogram's widths are its branch lengths, and
            // with no rule to read them against they measured nothing a
            // reader could name. A tree that is not a phylogram draws none.
            scale_bar: ScaleBar::default(),
            show_scale_bar: true,
            trait_columns: Vec::new(),
            node_glyphs: Vec::new(),
            clade_highlights: Vec::new(),
            refused: Vec::new(),
            joined: None,
            sheet_dealing: BTreeMap::new(),
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
    ///
    /// A projection chosen here, or with [`TreeTrack::circular`],
    /// [`TreeTrack::fan`] or [`TreeTrack::unrooted`], is the one drawn,
    /// whatever radial setting comes before or after it. A radial setting on
    /// its own draws a circle.
    pub fn projection(mut self, projection: TreeProjection) -> Self {
        self.projection = projection;
        self.projection_chosen = true;
        self
    }

    /// The projection a radial setting implies, unless one was chosen.
    fn imply(&mut self, projection: TreeProjection) {
        if !self.projection_chosen {
            self.projection = projection;
        }
    }

    /// Chooses orthogonal, diagonal or curved rectangular branches.
    ///
    /// Circular and unrooted projections retain their own geometry.
    pub fn branch_geometry(mut self, geometry: BranchGeometry) -> Self {
        self.branch_geometry = geometry;
        self
    }

    /// Draws a circular tree radiating outwards by default: a complete circle,
    /// unless [`TreeTrack::fan`] or [`TreeTrack::radial_sweep`] asks for less,
    /// before this or after it.
    pub fn circular(self) -> Self {
        self.projection(TreeProjection::Circular)
    }

    /// Draws an equal-angle tree around a topology-balanced central node.
    ///
    /// The source root is not used as the centre. Branch lengths are retained
    /// for a phylogram and topology alone is used for a cladogram.
    pub fn unrooted(self) -> Self {
        self.projection(TreeProjection::Unrooted)
    }

    /// Draws a circular fan covering `sweep_degrees` clockwise.
    pub fn fan(mut self, sweep_degrees: f64) -> Self {
        self.radial.sweep_degrees = finite_within(sweep_degrees, 10.0, 359.0, 240.0);
        self.projection(TreeProjection::Circular)
    }

    /// Sets the angle where a circular tree begins, in clockwise degrees.
    ///
    /// Zero is three o'clock and -90 is twelve o'clock.
    pub fn radial_start(mut self, degrees: f64) -> Self {
        self.imply(TreeProjection::Circular);
        if degrees.is_finite() {
            self.radial.start_degrees = degrees;
        }
        self
    }

    /// Sets the clockwise angular span of a circular tree in degrees.
    pub fn radial_sweep(mut self, degrees: f64) -> Self {
        self.imply(TreeProjection::Circular);
        self.radial.sweep_degrees = finite_within(degrees, 10.0, 360.0, 360.0);
        self
    }

    /// Chooses whether tips point away from or towards the centre.
    pub fn radial_direction(mut self, direction: RadialDirection) -> Self {
        self.imply(TreeProjection::Circular);
        self.radial.direction = direction;
        self
    }

    /// Sets the central gap as a fraction of the tree radius.
    pub fn inner_radius(mut self, fraction: f64) -> Self {
        self.imply(TreeProjection::Circular);
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
        self.imply(TreeProjection::Unrooted);
        if degrees.is_finite() {
            self.radial.start_degrees = degrees;
        }
        self
    }

    /// Reorients the owned tree around an internal node and marks the new
    /// root: an index, a name, or a [`NodeRef`] picking the clade by its tips
    /// or by a value.
    ///
    /// A node the tree does not have, or a sampled tip, leaves the tree
    /// unchanged, and the band says so under the tree, as
    /// [`TreeTrack::warnings`] does. Use [`Tree::reroot`](crate::Tree::reroot)
    /// directly when failure must be handled rather than reported.
    pub fn reroot(mut self, node: impl Into<NodeRef>) -> Self {
        let wanted = node.into();
        if let Some(node) = self.clade(&wanted, "not rerooted") {
            self.rerooted_with(|tree| tree.reroot(node));
        }
        self
    }

    /// Reorients the owned tree around an internal node with this exact name,
    /// as [`TreeTrack::reroot`] does given the name.
    pub fn reroot_named(self, name: &str) -> Self {
        self.reroot(NodeRef::named(name))
    }

    /// The internal node `wanted` names, or `None` with the reason said under
    /// the tree after `refused`: a node the tree does not have, or a tip. A
    /// clade that holds tips it was not named for is found, and that is said
    /// too.
    fn clade(&mut self, wanted: &NodeRef, refused: &str) -> Option<usize> {
        match wanted.find(&self.tree) {
            Err(why) => {
                self.refuse(format!("{refused}: {why}"));
                None
            }
            Ok(found) if self.tree.nodes()[found.node].is_leaf() => {
                self.refuse(format!("{refused}: {wanted} is a tip"));
                None
            }
            Ok(found) => {
                if let Some(also) = found.also {
                    self.refuse(format!("{wanted} {also}"));
                }
                Some(found.node)
            }
        }
    }

    /// Roots halfway along the edge leading to a monophyletic named outgroup.
    ///
    /// Missing, duplicate, internal or non-monophyletic names leave the tree
    /// unchanged, and the band says which under the tree. The new root is
    /// inserted without converting an outgroup tip into an internal node.
    pub fn reroot_outgroup<I, S>(mut self, names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let names: Vec<String> = names
            .into_iter()
            .map(|name| name.as_ref().to_string())
            .collect();
        let mut nodes = Vec::new();
        for name in &names {
            match self.tree.node_named(name) {
                None => {
                    self.refuse(format!("not rerooted: no tip is named {name}"));
                    return self;
                }
                Some(node) if !self.tree.nodes()[node].is_leaf() => {
                    self.refuse(format!("not rerooted: {name} is not a tip"));
                    return self;
                }
                Some(node) => nodes.push(node),
            }
        }
        if !self.rerooted_with(|tree| tree.reroot_outgroup(&nodes).is_some()) {
            self.refuse(format!(
                "not rerooted: {} {} not one clade of the tree",
                listed(&names),
                if names.len() == 1 { "is" } else { "are" }
            ));
        }
        self
    }

    /// Roots the owned phylogram at the midpoint of its weighted tip diameter.
    ///
    /// Missing, negative or non-finite branch lengths leave the tree
    /// unchanged, and the band says so under it.
    pub fn reroot_midpoint(mut self) -> Self {
        if !self.rerooted_with(|tree| tree.reroot_midpoint().is_some()) {
            let why = if self.tree.leaves().len() < 2 {
                "not rerooted at the midpoint: the tree has one tip"
            } else {
                "not rerooted at the midpoint: a branch has no length, or a negative one"
            };
            self.refuse(why.to_string());
        }
        self
    }

    /// Records a request this track was given and cannot carry out.
    fn refuse(&mut self, why: String) {
        if !self.refused.contains(&why) {
            self.refused.push(why);
        }
    }

    /// Reroots with `reroot`, and carries every fold and highlight asked for
    /// before it to the same clade after it. Whether it rerooted.
    ///
    /// A clade is its tips. Rerooting keeps each node where it was in the
    /// list and turns edges round, so an index names another clade after it:
    /// a fold of one lineage asked for before a reroot folded twenty-eight
    /// tips of three lineages after it. Each is found again by its tips, and
    /// one the new root splits is said under the tree rather than drawn.
    ///
    /// The root chosen is marked, unless [`TreeTrack::show_root`] says
    /// otherwise before the reroot or after it, and the folds a row cap made
    /// are emptied to be worked out against the new shape.
    fn rerooted_with(&mut self, reroot: impl FnOnce(&mut Tree) -> bool) -> bool {
        let tips = |tree: &Tree, node: usize| -> BTreeSet<usize> {
            if tree.nodes()[node].is_leaf() {
                return [node].into_iter().collect();
            }
            tree.descendants(node)
                .into_iter()
                .filter(|below| tree.nodes()[*below].is_leaf())
                .collect()
        };
        let folds: Vec<(usize, BTreeSet<usize>)> = self
            .collapsed
            .iter()
            .map(|node| (*node, tips(&self.tree, *node)))
            .collect();
        let named = std::mem::take(&mut self.fold_names);
        let fields: Vec<BTreeSet<usize>> = self
            .clade_highlights
            .iter()
            .map(|highlight| tips(&self.tree, highlight.node))
            .collect();
        if !reroot(&mut self.tree) {
            return false;
        }
        self.rerooted = true;
        self.folds = OnceLock::new();
        let found = |tree: &Tree, held: &BTreeSet<usize>| {
            let nodes: Vec<usize> = held.iter().copied().collect();
            tree.mrca(&nodes)
                .filter(|clade| tips(tree, *clade) == *held)
        };
        self.collapsed = BTreeSet::new();
        for (node, held) in folds {
            match found(&self.tree, &held) {
                Some(clade) => {
                    self.collapsed.insert(clade);
                    if let Some(name) = named.get(&node) {
                        self.fold_names.insert(clade, name.clone());
                    }
                }
                None => self.refuse(format!(
                    "not collapsed: the new root splits the clade of node {node}"
                )),
            }
        }
        let highlights = std::mem::take(&mut self.clade_highlights);
        for (mut highlight, held) in highlights.into_iter().zip(fields) {
            match found(&self.tree, &held) {
                Some(clade) => {
                    highlight.node = clade;
                    self.clade_highlights.push(highlight);
                }
                None => self.refuse(format!(
                    "not highlighted: the new root splits the clade of node {}",
                    highlight.node
                )),
            }
        }
        true
    }

    /// Draws or hides the selected root marker in rooted projections.
    ///
    /// A reroot draws it unless this hides it, written before the reroot or
    /// after it.
    pub fn show_root(mut self, show: bool) -> Self {
        self.show_root = Some(show);
        self
    }

    /// Whether the root is marked: as [`TreeTrack::show_root`] says, and where
    /// it says nothing, once a reroot has chosen it.
    fn shows_root(&self) -> bool {
        self.show_root.unwrap_or(self.rerooted)
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
        self.time = Some(key.into());
        self
    }

    /// Chooses whether time values increase or decrease from root to tips.
    pub fn time_direction(mut self, direction: TimeDirection) -> Self {
        self.time_direction = direction;
        self
    }

    /// Names what the time axis counts, as `year` or `years before present`,
    /// written as the axis's title, under its numbers or at the inner end of
    /// the rings, and not after the latest of them.
    pub fn time_unit(mut self, unit: impl Into<String>) -> Self {
        self.time_unit = Some(unit.into());
        self
    }

    /// Draws or hides the temporal axis created by [`TreeTrack::time`].
    pub fn show_time_axis(mut self, show: bool) -> Self {
        self.show_time_axis = show;
        self
    }

    /// The time axis as it is drawn, put together from its settings.
    fn time_axis(&self) -> Option<TimeAxis> {
        Some(TimeAxis {
            key: self.time.clone()?,
            direction: self.time_direction,
            unit: self.time_unit.clone(),
            show_axis: self.show_time_axis,
        })
    }

    /// Colours each incoming branch by one node annotation.
    ///
    /// [`TreeTrack::dnds`] colours the branches too, and where both are set
    /// the dN/dS colouring is drawn, whichever came first, and the band says
    /// this one was not.
    pub fn color_by(mut self, key: impl Into<String>) -> Self {
        self.color_by = Some(key.into());
        self
    }

    /// The key the branches are coloured by, where one is drawn: a dN/dS
    /// colouring takes the branches over.
    fn branch_key(&self) -> Option<&str> {
        self.color_by.as_deref().filter(|_| self.dnds.is_none())
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
        self.dnds = Some(key.into());
        self
    }

    /// Replaces the visible label of the dN/dS legend.
    ///
    /// Like every `dnds_` setting, it is drawn once [`TreeTrack::dnds`] has
    /// selected an annotation, and counts the same written before it or after.
    pub fn dnds_label(mut self, label: impl Into<String>) -> Self {
        self.dnds_label = label.into();
        self
    }

    /// Sets the inclusive interval treated as approximately neutral.
    ///
    /// Invalid, negative or reversed bounds leave the current interval
    /// unchanged. The default is `0.95..=1.05`.
    pub fn dnds_neutral_band(mut self, lower: f64, upper: f64) -> Self {
        if lower.is_finite() && upper.is_finite() && (0.0..=1.0).contains(&lower) && upper >= 1.0 {
            self.dnds_neutral_band = (lower, upper);
        }
        self
    }

    /// Sets where each side of the logarithmic dN/dS colour scale saturates.
    ///
    /// `4.0`, the default, makes ω ≥ 4 and ω ≤ 1/4 use the strongest warm and
    /// cool colours. Values between them retain continuous differences.
    pub fn dnds_saturation(mut self, fold: f64) -> Self {
        if fold.is_finite() && fold > 1.0 {
            self.dnds_saturation = fold;
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
            self.dnds_significance = Some(DnDsSignificance {
                key: key.into(),
                maximum,
            });
        } else {
            self.refuse(format!(
                "no dN/dS significance drawn: {} is no threshold a test can pass",
                text_rounded(maximum, 3)
            ));
        }
        self
    }

    /// The dN/dS colouring as it is drawn, put together from its settings.
    fn dnds_layer(&self) -> Option<DnDsLayer> {
        Some(DnDsLayer {
            key: self.dnds.clone()?,
            label: self.dnds_label.clone(),
            neutral_lower: self.dnds_neutral_band.0,
            neutral_upper: self.dnds_neutral_band.1,
            saturation: self.dnds_saturation,
            significance: self.dnds_significance.clone(),
        })
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

    /// Collapses one clade visually while preserving the source tree: an
    /// index, a name, or a [`NodeRef`] picking it by its tips, as
    /// `NodeRef::mrca(["S01", "S07"])`, or by a value, as
    /// `NodeRef::holding("lineage", "L4")`.
    ///
    /// A node the tree does not have, or a tip, folds nothing and is said
    /// under the tree, and so is a clade that holds tips it was not named for.
    pub fn collapse(mut self, node: impl Into<NodeRef>) -> Self {
        let wanted = node.into();
        if let Some(node) = self.clade(&wanted, "not collapsed") {
            self.collapsed.insert(node);
            self.folds = OnceLock::new();
            // Folded as the clade of a value, it is named by the value where
            // it has no name of its own: `L4 (16 tips)`, where it read
            // `S26 +15 more`, a sample named for a lineage.
            if let NodeRef::Holding { value, .. } = &wanted {
                if self.tree.nodes()[node].name.is_none() {
                    self.fold_names.insert(node, value.clone());
                }
            }
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
        self.folds = OnceLock::new();
        self
    }

    /// Every clade drawn folded, worked out the first time it is asked for.
    fn folded(&self) -> &BTreeSet<usize> {
        self.folds.get_or_init(|| self.fold_to_fit())
    }

    /// Collapses the smallest clades until the visible terminals fit the cap.
    ///
    /// Worked out when the tree is drawn rather than when the cap is written,
    /// because what is written after the cap changes what it has to fold. It
    /// used to fold straight away: a reroot written after it kept the folds of
    /// the shape it replaced, which are the wrong clades once the tips have
    /// moved; a clade collapsed by hand after it was folded on top of a tree
    /// already fitted to the cap; and `max_rows(None)` after it lifted the cap
    /// and kept every fold.
    fn fold_to_fit(&self) -> BTreeSet<usize> {
        let mut collapsed = self.collapsed.clone();
        let Some(cap) = self.max_rows else {
            return collapsed;
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
            return collapsed;
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
                collapsed.insert(*node);
            }
        }
        collapsed
    }

    /// Draws or hides a point at every visible internal node.
    pub fn show_nodes(mut self, show: bool) -> Self {
        self.show_nodes = show;
        self
    }

    /// Chooses how internal-node support is made visible.
    ///
    /// Values in either the `0..=1` or `0..=100` convention are recognised,
    /// once for the whole tree: when any value runs above one, every value is
    /// read out of a hundred, so a clade at 1 on a bootstrap tree is one
    /// percent and not full support. Their original representation is
    /// retained in labels and tooltips.
    pub fn support_style(mut self, style: SupportStyle) -> Self {
        self.support_style = style;
        self
    }

    /// Reads each clade's support from its numeric annotation `key`, as
    /// `posterior` in a BEAST tree or `prob` in a MrBayes one, or
    /// `support_1` for the first of IQ-TREE's values, as
    /// [`Tree::support_from`](crate::Tree::support_from) does.
    ///
    /// A key no internal node carries sets nothing and is said under the
    /// tree.
    pub fn support_from(mut self, key: &str) -> Self {
        if self.tree.support_from(key) == 0 {
            self.refuse(format!("no support read: no clade carries {key}"));
        }
        self
    }

    /// Hides visible support below `minimum`.
    ///
    /// `0.8` and `80.0` both mean eighty percent. Non-finite values reset the
    /// threshold to zero.
    pub fn support_threshold(mut self, minimum: f64) -> Self {
        self.support_threshold = threshold_fraction(minimum).unwrap_or(0.0);
        self
    }

    /// Labels each incoming branch with its own annotation `key`.
    ///
    /// Unlike [`TreeTrack::color_by`], values are not inherited from ancestor
    /// nodes. This makes the method suitable for mutations, gains, losses and
    /// other events that belong to one branch. Long text is fitted to the
    /// available segment while the complete value remains in its tooltip.
    pub fn branch_labels(mut self, key: impl Into<String>) -> Self {
        self.branch_labels = Some(key.into());
        self
    }

    /// Sets the font size of labels created by [`TreeTrack::branch_labels`].
    pub fn branch_label_size(mut self, size: f64) -> Self {
        self.branch_label_size = finite_within(size, 5.0, 18.0, 8.0);
        self
    }

    /// The branch labels as they are drawn, put together from their settings.
    fn branch_label_layer(&self) -> Option<BranchLabels> {
        Some(BranchLabels {
            key: self.branch_labels.clone()?,
            size: self.branch_label_size,
        })
    }

    /// Draws the automatically sized branch-length scale bar, which a
    /// phylogram draws by default: this undoes
    /// [`TreeTrack::show_scale_bar`]`(false)`.
    ///
    /// Cladograms, explicitly time-scaled trees and trees with no branch
    /// lengths omit it, because their widths do not measure evolutionary
    /// branch length.
    pub fn scale_bar(self) -> Self {
        self.show_scale_bar(true)
    }

    /// Draws or removes the branch-length scale bar, which a phylogram
    /// draws by default.
    ///
    /// Hidden, it stays hidden whatever length or unit is set for it before
    /// or after; those used to bring it back.
    pub fn show_scale_bar(mut self, show: bool) -> Self {
        self.show_scale_bar = show;
        self
    }

    /// Requests an exact scale-bar length in the tree's branch-length units.
    ///
    /// Values longer than the visible tree span are clamped to that span.
    /// Invalid values fall back to automatic sizing.
    pub fn scale_bar_length(mut self, length: f64) -> Self {
        self.scale_bar.length = (length.is_finite() && length > 0.0).then_some(length);
        self
    }

    /// Adds a unit such as `substitutions/site` to the scale-bar label.
    pub fn scale_bar_unit(mut self, unit: impl Into<String>) -> Self {
        let unit = unit.into();
        self.scale_bar.unit = (!unit.is_empty()).then_some(unit);
        self
    }

    /// Adds one metadata strip aligned to the visible terminal taxa.
    pub fn trait_column(mut self, column: TraitColumn) -> Self {
        self.trait_columns.push(column);
        self
    }

    /// Draws a sample sheet's columns beside the tips, joined to them by
    /// name, as `--traits` does on the command line:
    ///
    /// ```
    /// use karyon::{Sheet, Traits, Tree, TreeTrack};
    ///
    /// let sheet = Sheet::parse("sample\tlineage\nA\tL1\nB\tL2\nC\tL2\n")?;
    /// let tree = Tree::parse_newick("((A:1,B:1):1,C:2);")?;
    /// let track = TreeTrack::new(tree).traits(Traits::from_sheet(&sheet).strips(["lineage"]));
    /// assert_eq!(track.join().unwrap().matched.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// Each tip the sheet names takes its values as annotations, so the
    /// strips, [`TreeTrack::color_by`] and [`NodeRef::holding`] read them,
    /// and each column the sheet was spread into is drawn, widened to fit its
    /// heading, which a tree writes across the top of its strip: at the width
    /// a matrix gives it, `lineage` came out as `li…`.
    ///
    /// A sheet with no column spread draws nothing, and still gives the tips
    /// their values, for [`TreeTrack::color_by`] or a fold by value.
    ///
    /// A sheet that names none of the tips draws no strip, and a tip it does
    /// not name is counted under the tree, since its cells are drawn empty.
    /// [`TreeTrack::join`] says what was matched and what was left out on
    /// both sides.
    pub fn traits(mut self, traits: Traits) -> Self {
        let leaves = self.tree.leaf_names();
        let join = traits.join(leaves.iter().map(String::as_str));
        if join.matched.is_empty() {
            self.refuse("no strips: the sheet names none of the tips".to_string());
            self.joined = Some(join);
            return self;
        }
        for name in &join.matched {
            let (Some(values), Some(node)) = (traits.values(name), self.tree.node_named(name))
            else {
                continue;
            };
            let values: Vec<(String, AnnotationValue)> = values
                .iter()
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect();
            if let Some(into) = self.tree.annotations_mut(node) {
                for (key, value) in values {
                    into.insert(key, value);
                }
            }
        }
        if !join.without_row.is_empty() {
            let count = join.without_row.len();
            let (tips, have) = if count == 1 {
                ("tip", "has")
            } else {
                ("tips", "have")
            };
            let shown: Vec<&str> = join
                .without_row
                .iter()
                .take(3)
                .map(String::as_str)
                .collect();
            let more = if count > shown.len() {
                format!(" and {} more", count - shown.len())
            } else {
                String::new()
            };
            self.refuse(format!(
                "{count} {tips} {have} no row in the sheet: {}{more}",
                shown.join(", ")
            ));
        }
        for key in traits.sheet_keys() {
            if let Some(dealing) = traits.dealing_of(key) {
                self.sheet_dealing.insert(key.clone(), dealing);
            }
        }
        for column in traits.columns() {
            // The heading is drawn two points under the body size, and a
            // long column name is capped so it cannot eat the tree beside it.
            let heading = text_width(column.heading(), 9.0) + 8.0;
            self.trait_columns
                .push(column.clone().width(heading.clamp(14.0, 72.0)));
        }
        self.joined = Some(join);
        self
    }

    /// What the sheet handed to [`TreeTrack::traits`] matched and what it left
    /// out, where one was.
    pub fn join(&self) -> Option<&Join> {
        self.joined.as_ref()
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
    ///
    /// An index the tree does not have highlights nothing and is said under
    /// the tree.
    pub fn clade_highlight(mut self, mut highlight: CladeHighlight) -> Self {
        match highlight.wanted.find(&self.tree) {
            Ok(found) => {
                if let Some(also) = found.also {
                    self.refuse(format!("{} {also}", highlight.wanted));
                }
                highlight.node = found.node;
                self.clade_highlights.push(highlight);
            }
            Err(why) => self.refuse(format!("not highlighted: {why}")),
        }
        self
    }

    /// Highlights a clade by its exact internal or terminal name.
    ///
    /// A name no node has highlights nothing and is said under the tree.
    pub fn highlight_named(self, name: &str) -> Self {
        self.clade_highlight(CladeHighlight::new(NodeRef::named(name)))
    }

    /// The tree.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// The trait columns resolved the way a figure resolves them: which level
    /// every node is at, and which colour that level was given.
    ///
    /// A caller drawing the strips itself, on a canvas rather than into an SVG,
    /// needs this rather than the sheet. Working the levels out again from the
    /// sheet would be a second opinion about which blue is which, and two
    /// pictures of one tree that disagree about that are worse than one
    /// picture.
    ///
    /// A continuous column has no levels of its own, so it is given sixteen
    /// bands across its range. That is a choice this makes and the SVG does
    /// not, and it is why the bands are named by the range they cover.
    pub fn strips(&self, theme: &Theme) -> Vec<crate::TraitStrip> {
        const BANDS: usize = 16;
        let nodes = self.tree.nodes().len();
        let dealing = self.dealing(theme.palette.len());
        self.trait_columns
            .iter()
            .zip(dealing.columns)
            .map(|(column, dealt)| {
                let values = rectangular::branch_values(&self.tree, &column.key);
                let domain = rectangular::tree_domain(&self.tree, &column.key, dealt);
                let mut levels: Vec<crate::TraitLevel> = Vec::new();
                let mut of: Vec<Option<usize>> = vec![None; nodes];
                match column.scale {
                    TraitScale::Continuous => {
                        for band in 0..BANDS {
                            let low = domain.minimum
                                + (domain.maximum - domain.minimum) * band as f64 / BANDS as f64;
                            let high = domain.minimum
                                + (domain.maximum - domain.minimum) * (band + 1) as f64
                                    / BANDS as f64;
                            levels.push(crate::TraitLevel {
                                value: format!("{low:.3} to {high:.3}"),
                                color: crate::theme::mix(
                                    &theme.muted,
                                    &theme.accent,
                                    band as f64 / (BANDS - 1) as f64,
                                ),
                                symbol: None,
                            });
                        }
                        for (node, value) in values.iter().enumerate() {
                            let Some(fraction) = domain.fraction(*value) else {
                                continue;
                            };
                            let band = ((fraction * BANDS as f64) as usize).min(BANDS - 1);
                            of[node] = Some(band);
                        }
                    }
                    TraitScale::Categorical => {
                        let shaped = column.drawn_style(&domain, theme) == TraitStyle::Symbol;
                        for (value, index) in domain.levels() {
                            levels.push(crate::TraitLevel {
                                value: value.to_string(),
                                color: domain.paint(index, theme),
                                symbol: shaped.then(|| theme.symbol(index)),
                            });
                        }
                        for (node, value) in values.iter().enumerate() {
                            of[node] = domain.category(*value);
                        }
                    }
                }
                crate::TraitStrip {
                    key: column.key.clone(),
                    label: column.label.clone(),
                    levels,
                    of,
                }
            })
            .collect()
    }

    /// A key for the colours this track paints: the annotation
    /// [`TreeTrack::color_by`] colours the branches by, then every trait
    /// column, each level in the colour it is drawn in.
    ///
    /// Read off the count the branches and the strips are painted from, which
    /// runs over the whole tree in its own order. A key built from the sample
    /// sheet instead, with [`Traits::legend`](crate::track::traits::Traits::legend),
    /// numbers the levels in the order the sheet sorts its rows, and that is a
    /// different order: a figure of two countries printed each one's colour
    /// beside the other's name.
    ///
    /// Each level is keyed with the mark its column draws: a box for a strip,
    /// its own shape for a symbol column, and the two dots of a binary one.
    /// A column over the branch key keys the branch colours too, under its
    /// own heading, and a colour is not keyed twice.
    /// [`Figure::key`](crate::Figure::key) gathers this key in the figure's
    /// own theme; where it goes, and whether the figure needs it, is the
    /// caller's decision.
    pub fn legend(&self, theme: &Theme) -> crate::track::legend::Legend {
        // What a column's colours are dealt from and how it marks them: two
        // columns alike in both key the same entries, and a binary column
        // marks its values with two dots of its own whatever it shares.
        let marks = |key: &'_ str, scale: TraitScale, style: TraitStyle| {
            (key.to_string(), scale, style == TraitStyle::Binary)
        };
        let dealing = self.dealing(theme.palette.len());
        let mut legend = crate::track::legend::Legend::new();
        let mut keyed = Vec::new();
        if let Some(key) = self.branch_key() {
            let values = rectangular::branch_values(&self.tree, key);
            let scale = if rectangular::is_continuous(&values) {
                TraitScale::Continuous
            } else {
                TraitScale::Categorical
            };
            // A column over the same key keys these colours itself, in its
            // own marks and under its own heading. Keyed here as boxes, a
            // column of symbols beside them went unkeyed as a repeat.
            let covered = self.trait_columns.iter().any(|column| {
                marks(&column.key, column.scale, column.style)
                    == marks(key, scale, TraitStyle::Strip)
            });
            if !covered {
                let domain = rectangular::tree_domain(&self.tree, key, dealing.branches);
                legend = crate::track::traits::key_entries(
                    legend,
                    key,
                    scale,
                    TraitStyle::Strip,
                    &domain,
                    theme,
                );
                keyed.push(marks(key, scale, TraitStyle::Strip));
            }
        }
        for (column, dealt) in self.trait_columns.iter().zip(&dealing.columns) {
            let these = marks(&column.key, column.scale, column.style);
            if keyed.contains(&these) {
                continue;
            }
            let domain = rectangular::tree_domain(&self.tree, &column.key, *dealt);
            legend = crate::track::traits::key_entries(
                legend,
                &column.label,
                column.scale,
                column.drawn_style(&domain, theme),
                &domain,
                theme,
            );
            keyed.push(these);
        }
        legend
    }

    /// How the branch colour key deals the palette: as a trait column over
    /// the same key does, so the branches and the strip beside them agree,
    /// and otherwise from a stretch of its own.
    fn color_levels(&self) -> Dealt<'_> {
        self.dealing(crate::track::traits::STRIP_LEVELS).branches
    }

    /// How each trait column deals the palette, in the order of the columns,
    /// and how the branch key does.
    ///
    /// A column given a start of its own keeps it, as every column from a
    /// sample sheet does. The others, and the branch key where no column
    /// covers it, take one stretch of the palette each, in the order they
    /// were asked for, the way a sheet deals its columns: two strips of words
    /// both started at the palette's first colour, so each lineage was also a
    /// country. A key some column chose a start for starts there in every
    /// column over it.
    fn dealing(&self, colors: usize) -> Dealing<'_> {
        let colors = colors.max(1);
        let worded = |column: &TraitColumn| {
            column.scale == TraitScale::Categorical && column.style != TraitStyle::Binary
        };
        let covering = |key: &str| {
            self.trait_columns
                .iter()
                .position(|column| column.key == key && column.scale == TraitScale::Categorical)
        };
        let branch_key = self.branch_key().filter(|key| {
            covering(key).is_none()
                && !rectangular::is_continuous(&rectangular::branch_values(&self.tree, key))
        });
        let mut order: Vec<&str> = branch_key.into_iter().collect();
        for column in self.trait_columns.iter().filter(|column| worded(column)) {
            if !order.contains(&column.key.as_str()) {
                order.push(&column.key);
            }
        }
        let stride = (colors / order.len().max(1)).max(1);
        // A start a column chose or was dealt from its sheet, then the one a
        // strip of the joined sheet would take, then a stretch of this tree's
        // own among the keys it colours.
        let start = |key: &str| {
            self.trait_columns
                .iter()
                .filter(|column| column.key == key)
                .find_map(|column| {
                    column
                        .first
                        .or_else(|| column.stretch.map(|stretch| stretch.start(colors)))
                })
                .or_else(|| {
                    self.sheet_dealing
                        .get(key)
                        .map(|(_, stretch)| stretch.start(colors))
                })
                .or_else(|| {
                    order
                        .iter()
                        .position(|named| *named == key)
                        .map(|place| place * stride)
                })
                .unwrap_or(0)
        };
        let columns: Vec<Dealt<'_>> = self
            .trait_columns
            .iter()
            .map(|column| Dealt {
                levels: &column.levels,
                first: start(&column.key),
                colors: &column.colors,
            })
            .collect();
        // Branches coloured by a key of the sheet with no strip of it deal the
        // palette as that strip would, in the order the sheet gives the
        // levels, so a lineage is one colour whether its strip is drawn or
        // not. They took the order the tree meets the levels in, and L1 was
        // blue in one figure of a set and ochre in the next.
        let branches = match self.branch_key() {
            Some(key) => match covering(key) {
                Some(index) => columns[index],
                None => Dealt {
                    levels: self
                        .sheet_dealing
                        .get(key)
                        .map_or(&[][..], |(levels, _)| levels.as_slice()),
                    first: start(key),
                    colors: &[],
                },
            },
            None => Dealt::default(),
        };
        // The node glyphs take the colours no strip and no branch was dealt
        // first, and only then the ones they were. They took the palette from
        // its start, so a pie's first key was the colour of the first lineage
        // beside it, and a key saying `a` in the colour of `L1` reads as `L1`.
        let glyphs = if self.node_glyphs.is_empty() {
            Vec::new()
        } else {
            let palette = colors;
            let mut used = BTreeSet::new();
            for (column, dealt) in self.trait_columns.iter().zip(&columns) {
                if worded(column) {
                    for (_, index) in
                        rectangular::tree_domain(&self.tree, &column.key, *dealt).keyed()
                    {
                        used.insert(index % palette);
                    }
                }
            }
            if let Some(key) = self.branch_key() {
                if !rectangular::is_continuous(&rectangular::branch_values(&self.tree, key)) {
                    for (_, index) in rectangular::tree_domain(&self.tree, key, branches).keyed() {
                        used.insert(index % palette);
                    }
                }
            }
            let mut order: Vec<usize> =
                (0..palette).filter(|index| !used.contains(index)).collect();
            order.extend((0..palette).filter(|index| used.contains(index)));
            order
        };
        Dealing {
            columns,
            branches,
            glyphs,
        }
    }

    fn branch_scale(&self) -> Option<&ScaleBar> {
        // A tree with no branch lengths has nothing to measure, and room held
        // for its bar was an empty strip under every topology.
        let measured = self
            .tree
            .nodes()
            .iter()
            .any(|node| node.branch_length.is_some_and(|length| length > 0.0));
        // A time axis that cannot be drawn leaves the tree drawn by branch
        // length, and that drawing is measured by the bar as any other is.
        let timed = self
            .time_axis()
            .is_some_and(|time| self.tree.time_layout(&time.key, time.direction).is_some());
        (self.show_scale_bar && measured && self.shape == TreeShape::Phylogram && !timed)
            .then_some(&self.scale_bar)
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
            .map(|node| {
                text_width(
                    &terminal_label(&self.tree, *node, self.folded(), &self.fold_names),
                    size,
                )
            })
            .fold(0.0f64, f64::max)
            + 6.0
    }

    /// The line under a time axis its unit is written on, where it has one.
    fn time_title_room(&self) -> f64 {
        match self.time_unit.as_deref() {
            Some(unit) if !unit.trim().is_empty() => 14.0,
            _ => 0.0,
        }
    }

    fn axis_room(&self, theme: &Theme) -> f64 {
        let time = self
            .time_axis()
            .filter(|time| time.show_axis)
            .map_or(0.0, |_| {
                theme.font_size + theme.tokens.tick_length + 5.0 + self.time_title_room()
            });
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

    /// What this track was asked for and does not draw, a line each.
    ///
    /// A builder that cannot do what it was asked leaves the tree as it was,
    /// so a chain of settings never fails halfway: a reroot, a fold or a
    /// highlight naming a node the tree does not have, an outgroup that is not
    /// one clade. A setting can also say nothing about this tree: a key no
    /// node carries, a time axis some tip has no date for, a support
    /// threshold with no support drawn. Each of these used to leave a figure
    /// that looked finished and was not what was asked for. The band says
    /// every one of them under the tree, in these words, and this hands them
    /// to a caller that would rather stop.
    pub fn warnings(&self) -> Vec<String> {
        self.warnings_in(&Theme::default())
    }

    fn warnings_in(&self, theme: &Theme) -> Vec<String> {
        let mut said = self.refused.clone();
        // Two tips of one name are one tip to everything that finds a tip by
        // its name: a sheet joined by name, a row of an alignment, a clade
        // picked by its tips. The file is taken as it is, and said to be.
        let mut seen = std::collections::BTreeMap::<&str, usize>::new();
        for node in self.tree.leaves() {
            if let Some(name) = self.tree.nodes()[node]
                .name
                .as_deref()
                .filter(|name| !name.is_empty())
            {
                *seen.entry(name).or_default() += 1;
            }
        }
        let repeated: Vec<(&str, usize)> =
            seen.into_iter().filter(|(_, count)| *count > 1).collect();
        if let Some((name, count)) = repeated.first() {
            let others = match repeated.len() {
                1 => String::new(),
                2 => ", and 1 other name is repeated".to_string(),
                more => format!(", and {} other names are repeated", more - 1),
            };
            said.push(format!(
                "{count} tips are called {name}{others}: a sheet or a clade picked by name reaches the first"
            ));
        }
        let negative = self
            .tree
            .nodes()
            .iter()
            .filter(|clade| clade.branch_length.is_some_and(|length| length < 0.0))
            .count();
        if negative > 0 {
            let (branches, have) = if negative == 1 {
                ("branch has", "it is")
            } else {
                ("branches have", "they are")
            };
            said.push(format!(
                "{negative} {branches} a negative length, and {have} drawn as nought"
            ));
        }
        let carried = |key: &str| {
            (0..self.tree.nodes().len()).any(|node| self.tree.annotation(node, key).is_some())
        };
        if let (Some(key), Some(dnds)) = (&self.color_by, &self.dnds) {
            said.push(format!(
                "branches coloured by dN/dS ({dnds}), not by {key}: one colouring at a time"
            ));
        }
        if let Some(key) = self.branch_key() {
            if !carried(key) {
                said.push(format!("no branch is coloured: no node carries {key}"));
            } else {
                let values = rectangular::branch_values(&self.tree, key);
                if !rectangular::is_continuous(&values) {
                    let domain = rectangular::tree_domain(
                        &self.tree,
                        key,
                        self.dealing(theme.palette.len()).branches,
                    );
                    if domain.colors_repeat(theme) {
                        said.push(format!(
                            "{key} has {} values and the palette {} colours, so some branches of two values share one",
                            domain.keyed().len(),
                            theme.palette.len()
                        ));
                    }
                }
            }
        }
        if let Some(key) = &self.dnds {
            if !carried(key) {
                said.push(format!(
                    "no branch is coloured by dN/dS: no node carries {key}"
                ));
            }
        }
        if self.support_style == SupportStyle::None && self.support_threshold > 0.0 {
            said.push("no support is drawn: a threshold was set and no support style".to_string());
        }
        if self.support_style != SupportStyle::None
            && self
                .tree
                .nodes()
                .iter()
                .filter_map(|node| node.support)
                .any(|value| value > 100.0)
        {
            said.push("support above 100 is drawn as full support".to_string());
        }
        if let Some(time) = self.time_axis() {
            said.extend(self.time_warning(&time));
        }
        let dealing = self.dealing(theme.palette.len());
        let colors = theme.palette.len().max(1);
        let most = colors * 4;
        // Each colour a filled strip paints, with the column and the level it
        // paints, to find two strips that paint two things one colour.
        let mut painted: Vec<(String, &str, String)> = Vec::new();
        let mut shared: Vec<String> = Vec::new();
        for (column, dealt) in self.trait_columns.iter().zip(&dealing.columns) {
            let domain = rectangular::tree_domain(&self.tree, &column.key, *dealt);
            let levels = domain.keyed().len();
            let style = column.drawn_style(&domain, theme);
            if style == TraitStyle::Symbol && levels > most {
                said.push(format!(
                    "{}: {levels} values, and shapes and colours tell {most} apart",
                    column.label
                ));
            } else if style == TraitStyle::Symbol && column.style == TraitStyle::Strip {
                // Asked for a strip and drawn as shapes, which the figure
                // should say rather than the reader find out.
                said.push(format!(
                    "{}: {levels} values for {colors} colours, so each is a shape as well",
                    column.label
                ));
            }
            if style != TraitStyle::Strip || column.scale != TraitScale::Categorical {
                continue;
            }
            for (level, color) in domain.painted(theme) {
                let other = painted
                    .iter()
                    .find(|(_, named, had)| *had == color && *named != column.label.as_str());
                if let Some((was, named, _)) = other {
                    shared.push(format!("{named} {was} and {} {level}", column.label));
                }
                painted.push((level, &column.label, color));
            }
        }
        if let Some(first) = shared.first() {
            let others = match shared.len() {
                1 => String::new(),
                2 => ", and 1 other pair is too".to_string(),
                more => format!(", and {} other pairs are too", more - 1),
            };
            said.push(format!("{first} are one colour{others}"));
        }
        said
    }

    /// What is wrong with the time axis asked for, if anything.
    fn time_warning(&self, time: &TimeAxis) -> Option<String> {
        let tips = self.tree.leaves();
        let undated = tips
            .iter()
            .filter(|tip| {
                self.tree
                    .time_value(**tip, &time.key)
                    .map_or(true, |value| !value.is_finite())
            })
            .count();
        if undated > 0 {
            return Some(format!(
                "drawn by branch length: {undated} of {} tips have no number or date under {}",
                tips.len(),
                time.key
            ));
        }
        let Some(placed) = self.tree.time_layout(&time.key, time.direction) else {
            return Some(format!(
                "drawn by branch length: the tree cannot be placed on {}",
                time.key
            ));
        };
        let branches = self
            .tree
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(node, clade)| Some((node, clade.parent?)))
            .collect::<Vec<_>>();
        let backwards = branches
            .iter()
            .filter(|(node, parent)| {
                let (at, from) = (placed[*node].depth, placed[*parent].depth);
                let slack = 1e-9 * (at.abs() + from.abs()).max(1.0);
                match time.direction {
                    TimeDirection::Increasing => at < from - slack,
                    TimeDirection::Decreasing => at > from + slack,
                }
            })
            .count();
        (backwards > 0).then(|| {
            format!(
                "{backwards} of {} branches run backwards in {}",
                branches.len(),
                time.key
            )
        })
    }

    /// The lines [`TreeTrack::warnings`] takes across a band `width` pixels
    /// wide, and the size they are set at.
    fn warning_lines(&self, width: f64, theme: &Theme) -> (Vec<String>, f64) {
        let size = (theme.font_size - 1.0).max(6.0);
        let mut lines = Vec::new();
        for warning in self.warnings_in(theme) {
            let mut line = String::new();
            for word in warning.split(' ') {
                let longer = if line.is_empty() {
                    word.to_string()
                } else {
                    format!("{line} {word}")
                };
                if !line.is_empty() && text_width(&longer, size) > width - 8.0 {
                    lines.push(std::mem::replace(&mut line, word.to_string()));
                } else {
                    line = longer;
                }
            }
            lines.push(line);
        }
        (lines, size)
    }

    /// The room under the tree for its warnings.
    fn warning_room(&self, width: f64, theme: &Theme) -> f64 {
        let (lines, size) = self.warning_lines(width, theme);
        if lines.is_empty() {
            0.0
        } else {
            (lines.len() as f64 * (size + 4.0) + 6.0).ceil()
        }
    }

    /// Writes the warnings in the room kept for them at the foot of the band.
    fn draw_warnings(&self, ctx: &mut DrawContext<'_>) {
        let (lines, size) = self.warning_lines(ctx.band.w, ctx.theme);
        if lines.is_empty() {
            return;
        }
        let room = (lines.len() as f64 * (size + 4.0) + 6.0).ceil();
        let mut y = ctx.band.bottom() - room + 4.0;
        for line in &lines {
            y += size + 4.0;
            ctx.svg.text(
                ctx.band.x + 4.0,
                y - 4.0,
                line,
                &ctx.theme.muted,
                size,
                crate::svg::Anchor::Start,
            );
        }
    }

    /// The room across the top of a band `width` pixels wide for the
    /// headings of the columns or rings and the chips that key the layers.
    ///
    /// A rectangular tree's headings sit over its columns, at the right, and
    /// the chips keep to the left of them, so the two share the room. Rings
    /// have their headings in a row of their own across the top, with the
    /// chips under it: drawn in one row, the first chip covered the first
    /// ring's heading.
    fn annotation_header_room(&self, width: f64, theme: &Theme) -> f64 {
        let chips = match chip_rows(self, self.chip_room(width, theme), theme) {
            0 => 0.0,
            // Whole pixels, so a figure is as many pixels tall as it says.
            rows => (rows as f64 * chip_pitch(theme) + 1.5).ceil(),
        };
        let headings = if self.trait_columns.is_empty() {
            0.0
        } else {
            22.0
        };
        match self.projection {
            TreeProjection::Rectangular => chips.max(headings),
            TreeProjection::Circular | TreeProjection::Unrooted if self.ring_headings() => {
                chips + headings
            }
            TreeProjection::Circular | TreeProjection::Unrooted => chips,
        }
    }

    /// Whether the rings of a circle or a cloud are named across the top.
    ///
    /// Only where there are two or more to tell apart. One ring was named by
    /// a swatch in the corner, which two readers took for a key to a colour
    /// drawn nowhere, and the key under the figure names its levels already.
    fn ring_headings(&self) -> bool {
        self.trait_columns.len() > 1
    }

    /// Where the chips start below the top of the band.
    fn chip_top(&self) -> f64 {
        match self.projection {
            TreeProjection::Rectangular => 1.0,
            TreeProjection::Circular | TreeProjection::Unrooted => {
                if self.ring_headings() {
                    23.0
                } else {
                    1.0
                }
            }
        }
    }

    /// How much of a band `width` pixels wide the chips may take: a
    /// rectangular tree's columns keep the right of it for their headings.
    fn chip_room(&self, width: f64, theme: &Theme) -> f64 {
        match self.projection {
            TreeProjection::Rectangular => (width - self.trait_width(theme)).max(0.0),
            TreeProjection::Circular | TreeProjection::Unrooted => width,
        }
    }

    /// Draws the chips that key the layers where the header keeps them.
    fn draw_layer_chips(&self, ctx: &mut DrawContext<'_>) {
        let top = ctx.band.y + self.chip_top();
        let room = self.chip_room(ctx.band.w, ctx.theme);
        draw_annotation_legend(self, ctx, top, room);
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
        let terminals = visible_terminals(&self.tree, self.folded());
        let tips = terminals.len().max(1) as f64;
        let size = (theme.font_size - 1.0).min(self.row_height.max(1.0));

        // The room the labels take out of the radius, which is what a fitted
        // proportion gets wrong at the small end: it is a fixed number of
        // pixels, so on a small disc it is most of the radius and on a large
        // one it is a rim.
        let extent = if self.show_tips {
            terminals
                .iter()
                .map(|node| {
                    text_width(
                        &terminal_label(&self.tree, *node, self.folded(), &self.fold_names),
                        size,
                    )
                })
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
        // The ceiling is the width of the band this track is given, because
        // that is what the drawing fits the disc into: it takes the shorter of
        // the band's two sides, so a height past the band's width is whitespace
        // and a disc narrower than the height reserved for it is a disc that
        // stopped clearing its own labels.
        //
        // This was `x0 + width` for a while, which is the figure's inner right
        // edge and never moves. It reserved a height the drawing then did not
        // use: on the same two hundred tips at 900 px, changing nothing but
        // this track's own gutter label, the height stayed at 905.515 px while
        // the disc went from radius 339.40 to 283.40 and the gap between names
        // closed from 10.66 px to 8.90 against an 11 px body, which is the
        // collision the sizing exists to prevent.
        //
        // The cost is real and is the reason it was written the other way
        // first: the width of a band is shared, so a track stacked underneath
        // that asks for a wider gutter or a wider value axis narrows every band
        // and this tree gets shorter with it. The height of a round tree
        // genuinely follows the width it is given, and the alternative is a
        // figure whose names collide, so the width wins.
        wanted.clamp(RADIAL_DIAMETER, scale.width().max(RADIAL_DIAMETER))
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
        let time = self.time_axis();
        let scene = TreeScene::new(&self.tree, self.shape, time.as_ref(), self.folded());
        let tips = self.tip_width(ctx.theme, &scene);
        let axis_room = self.axis_room(ctx.theme);
        let traits = self.trait_width(ctx.theme);
        let header_room = self.annotation_header_room(band.w, ctx.theme);
        let warning_room = self.warning_room(band.w, ctx.theme);
        let (glyph_x, glyph_y) = self.rectangular_glyph_padding();
        let area = Rect {
            x: band.x + glyph_x,
            y: band.y + header_room + glyph_y,
            w: (band.w - tips - traits - glyph_x * 2.0).max(1.0),
            h: (band.h - axis_room - header_room - glyph_y * 2.0 - warning_room).max(1.0),
        };

        draw_rectangular_clade_highlights(self, ctx, &scene, area);
        draw_tree_scene(
            ctx,
            &self.tree,
            &self.fold_names,
            &scene,
            area,
            self.row_height,
            &color,
            self.line_width,
            self.branch_key(),
            self.color_levels(),
            self.dnds_layer().as_ref(),
            self.show_nodes,
            self.support_style,
            SupportReading::of(&self.tree, self.support_threshold),
            self.branch_label_layer().as_ref(),
            &self.rate_mixtures,
            &self.homoplasy_layers,
            &self.branch_event_layers,
            &self.branch_interval_layers,
            &self.ancestral_state_layers,
            self.branch_geometry,
            !self.show_tips,
        );
        draw_rectangular_node_glyphs(self, ctx, &scene, area);
        if self.shows_root() {
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
                let name = terminal_label(&self.tree, *node, self.folded(), &self.fold_names);
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
            self.folded(),
            &self.fold_names,
            area,
            tips + glyph_x,
            &self.trait_columns,
            &self.dealing(ctx.theme.palette.len()).columns,
            self.row_height,
        );
        if let Some(time) = time.as_ref().filter(|time| time.show_axis) {
            draw_time_axis(ctx, &scene, area, time);
        }
        if let Some(bar) = self.branch_scale() {
            draw_rectangular_scale_bar(ctx, &scene, area, bar);
        }
        self.draw_layer_chips(ctx);
        self.draw_warnings(ctx);
    }
}

/// How a tree's colour keys deal the palette, worked out once a drawing.
struct Dealing<'a> {
    /// One for each trait column, in their order.
    columns: Vec<Dealt<'a>>,
    /// The key the branches are coloured by.
    branches: Dealt<'a>,
    /// The palette in the order the node glyphs take it: the colours nothing
    /// else was dealt first. Empty where there are no glyphs.
    glyphs: Vec<usize>,
}

impl Dealing<'_> {
    /// The palette index the `index`th colour of the node glyphs is.
    fn glyph(&self, index: usize) -> usize {
        match self.glyphs.len() {
            0 => index,
            len => self.glyphs[index % len],
        }
    }
}

impl Track for TreeTrack {
    fn noun(&self) -> &str {
        "a phylogeny"
    }

    /// The key [`TreeTrack::legend`] builds, in the figure's own theme, so
    /// [`Figure::key`](crate::Figure::key) names every colour the tree
    /// paints. Built by hand from a theme the caller passed, the key of a
    /// dark figure named each level in the light palette's colour.
    fn key(
        &self,
        _region: &crate::region::Region,
        _px_per_bp: f64,
        theme: &Theme,
    ) -> Option<crate::track::legend::Legend> {
        let legend = self.legend(theme);
        (!legend.is_empty()).then_some(legend)
    }

    fn height(&self, scale: &Scale) -> f64 {
        match self.projection {
            TreeProjection::Rectangular => {
                let rows = visible_terminals(&self.tree, self.folded()).len().max(1) as f64;
                let (_, glyph_y) = self.rectangular_glyph_padding();
                rows * self.row_height
                    + glyph_y * 2.0
                    + if self.time_axis().is_some_and(|time| time.show_axis) {
                        22.0 + self.time_title_room()
                    } else {
                        0.0
                    }
                    + if self.branch_scale().is_some() {
                        22.0
                    } else {
                        0.0
                    }
                    + self.annotation_header_room(scale.width(), &Theme::default())
                    + self.warning_room(scale.width(), &Theme::default())
            }
            TreeProjection::Circular | TreeProjection::Unrooted => {
                self.radial_diameter(scale, &Theme::default())
                    + self.annotation_header_room(scale.width(), &Theme::default())
                    + self.warning_room(scale.width(), &Theme::default())
            }
        }
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn on_coordinates(&self) -> bool {
        false
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
