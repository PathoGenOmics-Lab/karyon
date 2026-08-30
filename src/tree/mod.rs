//! Phylogenies: reading them, and working out where their branches go.
//!
//! Nothing here draws. The module turns a Newick string into a [`Tree`] and
//! places every node at a depth and a row; five tracks read that placement, and
//! two of them, [`SnpTrack`](crate::SnpTrack) and
//! [`MatrixTrack`](crate::MatrixTrack), want something from a tree that has
//! nothing to do with drawing one.
//!
//! # The leaf order is half of what a tree is for
//!
//! [`Tree::leaves`] walks the clades depth first with the children in the order
//! the file listed them, which is the order the tips come out in when the tree
//! is drawn. Sort the rows of a panel beside it by that rather than by sample
//! name, which is what [`leaf_order`](crate::track::tree::leaf_order) does, and
//! a pattern carried by a clade stops being a speckle spread down the panel and
//! becomes a rectangle.
//!
//! # A Newick label does not say what it is
//!
//! The format writes support values and internal names in the same place, so
//! the parser decides: a label that reads as a number becomes
//! [`Clade::support`], anything else becomes [`Clade::name`], and there is no
//! switch to ask for the other reading. [`Tree::parse_newick`] discards square
//! bracket comments for compatibility, while [`Tree::parse_annotated_newick`]
//! preserves BEAST and NHX values as typed annotations. A `[&R]` or `[&U]`
//! prefix is stored as rootedness rather than mistaken for another node.
//!
//! # One layout, two ways of measuring depth
//!
//! [`Tree::layout`] takes a single flag, and it changes one of the two numbers
//! it produces. Depth is either the branch lengths added up or the branches
//! counted, so the same tree comes out as a phylogram or as a cladogram without
//! anything downstream knowing which it is looking at. Rows do not change
//! between them, and a parent always sits between its children, so a panel
//! sorted by [`Tree::leaves`] lines up either way.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::error::Error;

#[derive(Debug, Clone, Copy)]
struct TreeEdge {
    node: usize,
    length: Option<f64>,
    support: Option<f64>,
}

type EdgeAdjacency = Vec<Vec<TreeEdge>>;

mod mutation;
mod parse;

pub use mutation::{Mutation, Mutations};

#[cfg(test)]
mod tests;

use self::parse::*;

/// One typed value attached to a phylogenetic node or to the tree itself.
///
/// Annotated Newick written by BEAST commonly carries numbers, booleans,
/// strings and brace-delimited lists. Keeping those distinctions here avoids
/// making every renderer parse the same text again.
#[derive(Debug, Clone, PartialEq)]
pub enum AnnotationValue {
    /// Free text, including values whose syntax is not one of the types below.
    Text(String),
    /// A finite numeric value.
    Number(f64),
    /// A true or false flag.
    Boolean(bool),
    /// An ordered brace-delimited collection.
    List(Vec<AnnotationValue>),
}

impl AnnotationValue {
    /// The value as a number when it was numeric.
    pub fn as_number(&self) -> Option<f64> {
        match self {
            AnnotationValue::Number(value) => Some(*value),
            _ => None,
        }
    }

    /// The value as text when it was textual.
    pub fn as_text(&self) -> Option<&str> {
        match self {
            AnnotationValue::Text(value) => Some(value),
            _ => None,
        }
    }

    /// The value as a boolean when it was a flag.
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            AnnotationValue::Boolean(value) => Some(*value),
            _ => None,
        }
    }
}

/// What is drawn where an annotation has no value to show.
///
/// The trait strips already spell an absent annotation this way, so a value
/// that is not a number is spelled the same: both say there is nothing here to
/// read, which is the only thing a reader can act on.
pub const ABSENT: &str = "\u{2014}";

impl fmt::Display for AnnotationValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AnnotationValue::Text(value) => f.write_str(value),
            // Through `text_exact` rather than straight to the formatter. An
            // annotation is written wherever one is shown, and a BEAST or NHX
            // file carries whatever its sampler wrote: a rate of 1e-300 came
            // out of `write!` as three hundred and twenty four digits in a
            // branch event tooltip. This is the fifth place in the crate the
            // same mistake has been found, so it is fixed at the one point
            // that turns an annotation into text rather than at each reader.
            // A number field whose value is not a number is not a measurement
            // of anything, and this is the point every annotation becomes text
            // through, so it would otherwise reach a tooltip, a trait strip and
            // a node label alike reading `NaN`. The tree already spells an
            // annotation that is not there, and an annotation that is not a
            // number is the same fact about the same field.
            //
            // The infinities are not this. A file can carry a number too large
            // for an f64 and it parses to one, so reporting it verbatim tells
            // the reader what the file actually said.
            AnnotationValue::Number(value) if value.is_nan() => f.write_str(ABSENT),
            AnnotationValue::Number(value) => f.write_str(&crate::svg::text_exact(*value)),
            AnnotationValue::Boolean(value) => write!(f, "{value}"),
            AnnotationValue::List(values) => {
                f.write_str("{")?;
                for (index, value) in values.iter().enumerate() {
                    if index > 0 {
                        f.write_str(",")?;
                    }
                    write!(f, "{value}")?;
                }
                f.write_str("}")
            }
        }
    }
}

/// Stable, ordered annotations keyed by their source name.
pub type Annotations = BTreeMap<String, AnnotationValue>;

/// One node of a tree, internal or a leaf.
#[derive(Debug, Clone, PartialEq)]
pub struct Clade {
    /// Leaf name, or an internal label when the file carried one.
    pub name: Option<String>,
    /// Length of the branch leading to this node.
    pub branch_length: Option<f64>,
    /// Support value on the branch leading to this node, when the internal
    /// label parsed as a number.
    pub support: Option<f64>,
    /// Indices of the children, empty for a leaf.
    pub children: Vec<usize>,
    /// Index of the parent, `None` for the root.
    pub parent: Option<usize>,
}

impl Clade {
    /// Whether this node is a leaf.
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

/// Where one node sits once the tree has been laid out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Placement {
    /// Index of the node.
    pub node: usize,
    /// Distance from the root, in branch length units.
    pub depth: f64,
    /// Row of this node: a leaf's own index, or the mean of its children's.
    pub row: f64,
}

/// Direction in which calendar or height values run from root to tips.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeDirection {
    /// Dates increase towards the tips, as decimal calendar years do.
    #[default]
    Increasing,
    /// Values decrease towards the tips, as heights before present do.
    Decreasing,
}

/// A rooted phylogeny.
///
/// Nodes live in one flat list and refer to each other by index, so a tree can
/// be walked without recursion and a deep one cannot overflow the stack.
///
/// ```
/// use karyon::tree::Tree;
///
/// let tree = Tree::parse_newick("((ERR01:0.01,ERR02:0.012)0.98:0.04,ERR03:0.06);").unwrap();
/// assert_eq!(tree.leaf_names(), ["ERR01", "ERR02", "ERR03"]);
/// assert_eq!(tree.leaf_count(), 3);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    nodes: Vec<Clade>,
    root: usize,
    annotations: Vec<Annotations>,
    tree_annotations: Annotations,
    rooted: Option<bool>,
}

impl Tree {
    /// Reads one Newick string.
    ///
    /// Handles nested clades, branch lengths, quoted names and internal labels,
    /// reading an internal label as a support value when it parses as a number
    /// and keeping it as a name when it does not. Trailing semicolon optional.
    /// A doubled quote inside a quoted name is one literal quote, so
    /// `'O''Brien'` is a single tip.
    ///
    /// An empty label is a tip too: `(,,(,));` is four unnamed leaves, and a
    /// branch length written on an empty slot belongs to that leaf rather than
    /// to the clade around it.
    ///
    /// Square bracket comments are skipped wherever they appear, which is what
    /// lets a file straight out of RAxML or BEAST be read: those write a `[&R]`
    /// rootedness marker before the tree and `[&support=...]` annotations
    /// inside it. Nothing in a comment is kept, including NHX annotations.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidNewick`] for unbalanced parentheses, a stray
    /// comma, a branch length that is not a number, or an empty string.
    pub fn parse_newick(input: &str) -> Result<Self, Error> {
        parse_newick_impl(input, false)
    }

    /// Reads Newick while preserving BEAST, NHX and ordinary comments.
    ///
    /// BEAST fields such as `[&height=12.4,location="Peru"]` and NHX fields
    /// such as `[&&NHX:S=human:B=95]` are attached to the node immediately
    /// before them. Prefix markers `[&R]` and `[&U]` set [`Tree::rooted`].
    /// [`Tree::parse_newick`] remains the compatibility parser and deliberately
    /// discards comments.
    pub fn parse_annotated_newick(input: &str) -> Result<Self, Error> {
        parse_newick_impl(input, true)
    }

    /// Reads the first tree from a Nexus `trees` block.
    ///
    /// A `translate` table is applied to leaf labels and annotations on the
    /// Newick expression are preserved. This intentionally reads the portable
    /// tree subset rather than trying to interpret arbitrary Nexus data blocks.
    pub fn parse_nexus(input: &str) -> Result<Self, Error> {
        parse_nexus(input)
    }

    /// The nodes, in the order they were parsed.
    pub fn nodes(&self) -> &[Clade] {
        &self.nodes
    }

    /// Index of the root.
    pub fn root(&self) -> usize {
        self.root
    }

    /// Annotations attached to `node`, or `None` when the index is absent.
    pub fn annotations(&self, node: usize) -> Option<&Annotations> {
        self.annotations.get(node)
    }

    /// Mutable annotations for `node`, for metadata assembled in Rust.
    pub fn annotations_mut(&mut self, node: usize) -> Option<&mut Annotations> {
        self.annotations.get_mut(node)
    }

    /// One annotation attached to `node`.
    pub fn annotation(&self, node: usize, key: &str) -> Option<&AnnotationValue> {
        self.annotations(node)?.get(key)
    }

    /// Annotations carried by the tree rather than by one node.
    pub fn tree_annotations(&self) -> &Annotations {
        &self.tree_annotations
    }

    /// Mutable annotations carried by the whole tree.
    pub fn tree_annotations_mut(&mut self) -> &mut Annotations {
        &mut self.tree_annotations
    }

    /// Whether the source explicitly marked the tree as rooted or unrooted.
    ///
    /// Plain Newick generally carries no such declaration and returns `None`.
    pub fn rooted(&self) -> Option<bool> {
        self.rooted
    }

    /// Finds the first node with the exact `name`.
    pub fn node_named(&self, name: &str) -> Option<usize> {
        self.nodes
            .iter()
            .position(|node| node.name.as_deref() == Some(name))
    }

    /// Leaf indices, left to right as the tree is drawn.
    pub fn leaves(&self) -> Vec<usize> {
        let mut order = Vec::new();
        let mut stack = vec![self.root];
        // Depth first, children in file order, so the drawn order is the order
        // the file meant.
        while let Some(index) = stack.pop() {
            let node = &self.nodes[index];
            if node.is_leaf() {
                order.push(index);
            } else {
                for child in node.children.iter().rev() {
                    stack.push(*child);
                }
            }
        }
        order
    }

    /// Leaf names in drawn order, unnamed leaves coming back empty.
    pub fn leaf_names(&self) -> Vec<String> {
        self.leaves()
            .into_iter()
            .map(|index| self.nodes[index].name.clone().unwrap_or_default())
            .collect()
    }

    /// How many leaves the tree has.
    pub fn leaf_count(&self) -> usize {
        self.leaves().len()
    }

    /// Ancestors of `node`, from its parent back to the root.
    pub fn ancestors(&self, node: usize) -> Vec<usize> {
        let mut ancestors = Vec::new();
        let mut current = self.nodes.get(node).and_then(|clade| clade.parent);
        while let Some(index) = current {
            ancestors.push(index);
            current = self.nodes[index].parent;
        }
        ancestors
    }

    /// Descendants of `node` in drawn order, excluding `node` itself.
    pub fn descendants(&self, node: usize) -> Vec<usize> {
        let Some(clade) = self.nodes.get(node) else {
            return Vec::new();
        };
        let mut descendants = Vec::new();
        let mut stack: Vec<usize> = clade.children.iter().rev().copied().collect();
        while let Some(index) = stack.pop() {
            descendants.push(index);
            for child in self.nodes[index].children.iter().rev() {
                stack.push(*child);
            }
        }
        descendants
    }

    /// Number of leaves below `node`, counting the node when it is a leaf.
    pub fn clade_size(&self, node: usize) -> usize {
        let Some(clade) = self.nodes.get(node) else {
            return 0;
        };
        if clade.is_leaf() {
            1
        } else {
            self.descendants(node)
                .into_iter()
                .filter(|index| self.nodes[*index].is_leaf())
                .count()
        }
    }

    /// Most recent common ancestor of every node in `nodes`.
    ///
    /// Returns `None` for an empty list or an index outside this tree.
    pub fn mrca(&self, nodes: &[usize]) -> Option<usize> {
        let first = *nodes.first()?;
        self.nodes.get(first)?;
        let mut candidates = vec![first];
        candidates.extend(self.ancestors(first));
        candidates.into_iter().find(|candidate| {
            nodes.iter().skip(1).all(|node| {
                self.nodes.get(*node).is_some()
                    && (*node == *candidate || self.ancestors(*node).contains(candidate))
            })
        })
    }

    /// Reverses the child order of one internal node.
    ///
    /// A rotation changes how a tree is drawn but not the clades it contains.
    pub fn rotate(&mut self, node: usize) -> bool {
        let Some(clade) = self.nodes.get_mut(node) else {
            return false;
        };
        if clade.children.len() < 2 {
            return false;
        }
        clade.children.reverse();
        true
    }

    /// Orders every split by the number of descendant leaves.
    ///
    /// With `largest_first`, the largest clade is drawn first. Equal clades
    /// keep their source order, so repeated calls are deterministic.
    pub fn ladderize(&mut self, largest_first: bool) {
        let mut sizes = vec![0usize; self.nodes.len()];
        for node in self.postorder() {
            sizes[node] = if self.nodes[node].is_leaf() {
                1
            } else {
                self.nodes[node]
                    .children
                    .iter()
                    .map(|child| sizes[*child])
                    .sum()
            };
        }
        for clade in &mut self.nodes {
            if largest_first {
                clade
                    .children
                    .sort_by(|left, right| sizes[*right].cmp(&sizes[*left]));
            } else {
                clade
                    .children
                    .sort_by(|left, right| sizes[*left].cmp(&sizes[*right]));
            }
        }
    }

    /// Reorients the tree around internal `node`, preserving every undirected edge.
    ///
    /// Branch lengths and support stay on their edge when its direction
    /// changes. The new root has neither, because it has no incoming branch. A
    /// leaf is refused because turning a sampled tip into an internal root
    /// would silently remove it from the tip set.
    pub fn reroot(&mut self, node: usize) -> bool {
        let Some(target) = self.nodes.get(node) else {
            return false;
        };
        if target.is_leaf() && self.nodes.len() > 1 {
            return false;
        }
        let adjacency = self.edge_adjacency();
        self.orient_from(node, &adjacency);
        true
    }

    /// Roots the tree halfway along the branch leading to a monophyletic outgroup.
    ///
    /// Every index must name a leaf, and those leaves must be exactly the tips
    /// below their current MRCA. A new degree-two root is inserted on that
    /// clade's incoming edge, so no sampled tip is converted into an internal
    /// node. Returns its stable index, or `None` without changing the tree when
    /// the selection is empty, invalid, non-monophyletic or already spans the
    /// whole tree.
    pub fn reroot_outgroup(&mut self, outgroup: &[usize]) -> Option<usize> {
        let selected: BTreeSet<usize> = outgroup.iter().copied().collect();
        if selected.is_empty()
            || selected.len() != outgroup.len()
            || selected
                .iter()
                .any(|node| !self.nodes.get(*node).is_some_and(Clade::is_leaf))
        {
            return None;
        }
        let nodes: Vec<usize> = selected.iter().copied().collect();
        let mrca = self.mrca(&nodes)?;
        if mrca == self.root {
            return None;
        }
        let clade_leaves: BTreeSet<usize> = if self.nodes[mrca].is_leaf() {
            [mrca].into_iter().collect()
        } else {
            self.descendants(mrca)
                .into_iter()
                .filter(|node| self.nodes[*node].is_leaf())
                .collect()
        };
        if clade_leaves != selected {
            return None;
        }
        self.reroot_on_edge(mrca, 0.5)
    }

    /// Roots a fully weighted tree at the midpoint of its longest tip-to-tip path.
    ///
    /// Every branch must have a finite, non-negative length. A midpoint inside
    /// an edge creates one degree-two root and splits that edge exactly; a
    /// midpoint already on an internal node reuses it. Returns the root index,
    /// or `None` without mutation when the required distances are unavailable.
    pub fn reroot_midpoint(&mut self) -> Option<usize> {
        let leaves = self.leaves();
        if leaves.len() < 2 {
            return None;
        }
        let mut adjacency = vec![Vec::<(usize, f64)>::new(); self.nodes.len()];
        for (child, clade) in self.nodes.iter().enumerate() {
            let Some(parent) = clade.parent else {
                continue;
            };
            let length = clade.branch_length?;
            if !length.is_finite() || length < 0.0 {
                return None;
            }
            adjacency[parent].push((child, length));
            adjacency[child].push((parent, length));
        }

        fn distances(
            start: usize,
            adjacency: &[Vec<(usize, f64)>],
        ) -> (Vec<f64>, Vec<Option<usize>>) {
            let mut distance = vec![0.0; adjacency.len()];
            let mut previous = vec![None; adjacency.len()];
            let mut stack = vec![(start, None)];
            while let Some((node, parent)) = stack.pop() {
                for (next, length) in &adjacency[node] {
                    if Some(*next) == parent {
                        continue;
                    }
                    distance[*next] = distance[node] + length;
                    previous[*next] = Some(node);
                    stack.push((*next, Some(node)));
                }
            }
            (distance, previous)
        }

        let farthest_leaf = |distance: &[f64]| {
            leaves.iter().copied().max_by(|left, right| {
                distance[*left]
                    .total_cmp(&distance[*right])
                    .then_with(|| right.cmp(left))
            })
        };
        let (first_distances, _) = distances(leaves[0], &adjacency);
        let first = farthest_leaf(&first_distances)?;
        let (diameter_distances, previous) = distances(first, &adjacency);
        let last = farthest_leaf(&diameter_distances)?;
        let diameter = diameter_distances[last];
        if !diameter.is_finite() || diameter <= 0.0 {
            return None;
        }

        let mut path = vec![last];
        while *path.last()? != first {
            path.push(previous[*path.last()?]?);
        }
        path.reverse();
        let midpoint = diameter / 2.0;
        let tolerance = diameter.max(1.0) * 1e-12;
        let mut walked = 0.0;
        for edge in path.windows(2) {
            let from = edge[0];
            let to = edge[1];
            if (midpoint - walked).abs() <= tolerance {
                return self.reroot(from).then_some(from);
            }
            let length = adjacency[from]
                .iter()
                .find_map(|(node, length)| (*node == to).then_some(*length))?;
            let next = walked + length;
            if (midpoint - next).abs() <= tolerance {
                return self.reroot(to).then_some(to);
            }
            if midpoint > walked && midpoint < next {
                let original_child = if self.nodes[from].parent == Some(to) {
                    from
                } else {
                    to
                };
                let distance_from_child = if original_child == from {
                    midpoint - walked
                } else {
                    next - midpoint
                };
                return self.reroot_on_edge(original_child, distance_from_child / length);
            }
            walked = next;
        }
        None
    }

    fn edge_adjacency(&self) -> EdgeAdjacency {
        let mut adjacency = vec![Vec::new(); self.nodes.len()];
        for (child, clade) in self.nodes.iter().enumerate() {
            if let Some(parent) = clade.parent {
                adjacency[parent].push(TreeEdge {
                    node: child,
                    length: clade.branch_length,
                    support: clade.support,
                });
                adjacency[child].push(TreeEdge {
                    node: parent,
                    length: clade.branch_length,
                    support: clade.support,
                });
            }
        }
        adjacency
    }

    fn orient_from(&mut self, root: usize, adjacency: &EdgeAdjacency) {
        for clade in &mut self.nodes {
            clade.parent = None;
            clade.children.clear();
            clade.branch_length = None;
            clade.support = None;
        }
        let mut stack = vec![(root, None, None, None)];
        while let Some((current, parent, length, support)) = stack.pop() {
            self.nodes[current].parent = parent;
            self.nodes[current].branch_length = length;
            self.nodes[current].support = support;
            if let Some(parent) = parent {
                self.nodes[parent].children.push(current);
            }
            for edge in adjacency[current].iter().rev() {
                if Some(edge.node) != parent {
                    stack.push((edge.node, Some(current), edge.length, edge.support));
                }
            }
        }
        self.root = root;
        self.rooted = Some(true);
    }

    fn reroot_on_edge(&mut self, child: usize, fraction_from_child: f64) -> Option<usize> {
        let parent = self.nodes.get(child)?.parent?;
        if !fraction_from_child.is_finite() || !(0.0..=1.0).contains(&fraction_from_child) {
            return None;
        }
        let length = self.nodes[child].branch_length;
        let support = self.nodes[child].support;
        let child_length = length.map(|value| value * fraction_from_child);
        let parent_length = length.map(|value| value * (1.0 - fraction_from_child));
        let mut adjacency = self.edge_adjacency();
        adjacency[child].retain(|edge| edge.node != parent);
        adjacency[parent].retain(|edge| edge.node != child);

        let root = self.nodes.len();
        self.nodes.push(Clade {
            name: None,
            branch_length: None,
            support: None,
            children: Vec::new(),
            parent: None,
        });
        self.annotations.push(Annotations::new());
        adjacency.push(Vec::new());
        adjacency[root].push(TreeEdge {
            node: child,
            length: child_length,
            support,
        });
        adjacency[child].push(TreeEdge {
            node: root,
            length: child_length,
            support,
        });
        adjacency[root].push(TreeEdge {
            node: parent,
            length: parent_length,
            support: None,
        });
        adjacency[parent].push(TreeEdge {
            node: root,
            length: parent_length,
            support: None,
        });
        self.orient_from(root, &adjacency);
        Some(root)
    }

    /// Copies the clade rooted at `node` into a standalone tree.
    pub fn subtree(&self, node: usize) -> Option<Tree> {
        self.extract(node)
    }

    /// Replaces the descendants of `node` with one terminal clade.
    ///
    /// The selected node keeps its name, annotations and incoming branch.
    /// This is a data operation; [`TreeTrack`](crate::TreeTrack) also supports
    /// non-destructive visual collapsing.
    pub fn collapse(&mut self, node: usize) -> bool {
        let Some(clade) = self.nodes.get_mut(node) else {
            return false;
        };
        if clade.children.is_empty() {
            return false;
        }
        clade.children.clear();
        let Some(compact) = self.extract(self.root) else {
            return false;
        };
        *self = compact;
        true
    }

    /// Every node placed at a depth and a row.
    ///
    /// Depth is the distance from the root along the branches, so a tree with
    /// branch lengths comes out as a phylogram. `cladogram` ignores the lengths
    /// and uses the number of branches instead, which is what to do when the
    /// lengths are missing or meaningless.
    pub fn layout(&self, cladogram: bool) -> Vec<Placement> {
        let leaves = self.leaves();
        let mut rows: Vec<Option<f64>> = vec![None; self.nodes.len()];
        for (row, index) in leaves.iter().enumerate() {
            rows[*index] = Some(row as f64);
        }

        // Depths first, walking down from the root.
        let mut depths = vec![0.0f64; self.nodes.len()];
        let mut stack = vec![self.root];
        while let Some(index) = stack.pop() {
            let parent_depth = self.nodes[index]
                .parent
                .map_or(0.0, |parent| depths[parent]);
            let step = if cladogram {
                if self.nodes[index].parent.is_some() {
                    1.0
                } else {
                    0.0
                }
            } else {
                self.nodes[index].branch_length.unwrap_or(0.0).max(0.0)
            };
            depths[index] = parent_depth + step;
            stack.extend(self.nodes[index].children.iter().copied());
        }

        // Rows next, walking back up so a parent sees finished children.
        let order = self.postorder();
        for index in order {
            if rows[index].is_some() {
                continue;
            }
            let children = &self.nodes[index].children;
            let sum: f64 = children.iter().filter_map(|child| rows[*child]).sum();
            let count = children
                .iter()
                .filter(|child| rows[**child].is_some())
                .count();
            rows[index] = Some(if count > 0 { sum / count as f64 } else { 0.0 });
        }

        (0..self.nodes.len())
            .map(|node| Placement {
                node,
                depth: depths[node],
                row: rows[node].unwrap_or(0.0),
            })
            .collect()
    }

    /// Places nodes by a numeric annotation such as `date` or `height`.
    ///
    /// Every tip must carry `key`. An unannotated internal node is inferred
    /// from its children and their branch lengths, subtracting lengths for
    /// [`TimeDirection::Increasing`] and adding them for `Decreasing`.
    /// Returns `None` when a tip is missing a finite value.
    pub fn time_layout(&self, key: &str, direction: TimeDirection) -> Option<Vec<Placement>> {
        let rows = self.layout(false);
        let mut values: Vec<Option<f64>> = (0..self.nodes.len())
            .map(|node| {
                self.annotation(node, key)
                    .and_then(AnnotationValue::as_number)
            })
            .collect();
        if self
            .leaves()
            .iter()
            .any(|node| values[*node].map_or(true, |value| !value.is_finite()))
        {
            return None;
        }
        for node in self.postorder() {
            if values[node].is_some() || self.nodes[node].is_leaf() {
                continue;
            }
            let estimates: Vec<f64> = self.nodes[node]
                .children
                .iter()
                .filter_map(|child| {
                    let value = values[*child]?;
                    let branch = self.nodes[*child].branch_length.unwrap_or(0.0).max(0.0);
                    Some(match direction {
                        TimeDirection::Increasing => value - branch,
                        TimeDirection::Decreasing => value + branch,
                    })
                })
                .collect();
            if !estimates.is_empty() {
                values[node] = Some(estimates.iter().sum::<f64>() / estimates.len() as f64);
            }
        }
        if values.iter().any(Option::is_none) {
            return None;
        }
        Some(
            rows.into_iter()
                .map(|placement| Placement {
                    depth: values[placement.node].unwrap_or(0.0),
                    ..placement
                })
                .collect(),
        )
    }

    /// The deepest leaf, in whatever units the layout used.
    pub fn max_depth(&self, cladogram: bool) -> f64 {
        self.layout(cladogram)
            .iter()
            .filter(|placement| self.nodes[placement.node].is_leaf())
            .map(|placement| placement.depth)
            .fold(0.0f64, f64::max)
    }

    fn extract(&self, root: usize) -> Option<Tree> {
        self.nodes.get(root)?;
        let mut order = vec![root];
        order.extend(self.descendants(root));
        let mut mapping = vec![None; self.nodes.len()];
        for (new, old) in order.iter().enumerate() {
            mapping[*old] = Some(new);
        }

        let mut nodes = Vec::with_capacity(order.len());
        let mut annotations = Vec::with_capacity(order.len());
        for old in order {
            let source = &self.nodes[old];
            let parent = if old == root {
                None
            } else {
                source.parent.and_then(|parent| mapping[parent])
            };
            nodes.push(Clade {
                name: source.name.clone(),
                branch_length: if old == root {
                    None
                } else {
                    source.branch_length
                },
                support: source.support,
                children: source
                    .children
                    .iter()
                    .filter_map(|child| mapping[*child])
                    .collect(),
                parent,
            });
            annotations.push(self.annotations[old].clone());
        }
        Some(Tree {
            nodes,
            root: 0,
            annotations,
            tree_annotations: self.tree_annotations.clone(),
            rooted: self.rooted,
        })
    }

    /// Node indices with every child before its parent.
    fn postorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![self.root];
        while let Some(index) = stack.pop() {
            out.push(index);
            stack.extend(self.nodes[index].children.iter().copied());
        }
        out.reverse();
        out
    }
}
