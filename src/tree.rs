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

impl fmt::Display for AnnotationValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AnnotationValue::Text(value) => f.write_str(value),
            AnnotationValue::Number(value) => write!(f, "{value}"),
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

fn parse_newick_impl(input: &str, preserve_annotations: bool) -> Result<Tree, Error> {
    let text = input.trim().trim_end_matches(';').trim();
    if text.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "empty tree",
        });
    }

    let mut nodes: Vec<Clade> = Vec::new();
    let mut annotations: Vec<Annotations> = Vec::new();
    let mut tree_annotations = Annotations::new();
    let mut rooted = None;
    let mut stack: Vec<usize> = Vec::new();
    let mut current: Option<usize> = None;
    let mut chars = text.chars().peekable();

    while let Some(c) = chars.next() {
        match c {
            '(' => {
                let parent = stack.last().copied();
                if parent.is_none() && !nodes.is_empty() {
                    return Err(Error::InvalidNewick {
                        reason: "more than one root",
                    });
                }
                let index = add_node(&mut nodes, &mut annotations, parent);
                stack.push(index);
                current = None;
            }
            ')' => {
                if let (None, Some(parent)) = (current, stack.last().copied()) {
                    add_node(&mut nodes, &mut annotations, Some(parent));
                }
                let closed = stack.pop().ok_or(Error::InvalidNewick {
                    reason: "unbalanced parentheses",
                })?;
                current = Some(closed);
            }
            ',' => {
                if stack.is_empty() {
                    return Err(Error::InvalidNewick {
                        reason: "comma outside any clade",
                    });
                }
                if current.is_none() {
                    add_node(&mut nodes, &mut annotations, stack.last().copied());
                }
                current = None;
            }
            ':' => {
                let mut number = String::new();
                while let Some(next) = chars.peek() {
                    if next.is_ascii_digit() || matches!(next, '.' | '-' | '+' | 'e' | 'E') {
                        number.push(*next);
                        chars.next();
                    } else {
                        break;
                    }
                }
                let length = number.parse::<f64>().map_err(|_| Error::InvalidNewick {
                    reason: "branch length is not a number",
                })?;
                let target = match current {
                    Some(index) => index,
                    None => {
                        let parent = stack.last().copied().ok_or(Error::InvalidNewick {
                            reason: "branch length with nothing to attach to",
                        })?;
                        let index = add_node(&mut nodes, &mut annotations, Some(parent));
                        current = Some(index);
                        index
                    }
                };
                nodes[target].branch_length = Some(length);
            }
            '[' => {
                let mut comment = String::new();
                for next in chars.by_ref() {
                    if next == ']' {
                        break;
                    }
                    comment.push(next);
                }
                if preserve_annotations {
                    let (root_marker, fields) = parse_comment(&comment);
                    if let Some(value) = root_marker {
                        rooted = Some(value);
                    }
                    if let Some(node) = current {
                        annotations[node].extend(fields);
                    } else {
                        tree_annotations.extend(fields);
                    }
                }
            }
            c if c.is_whitespace() => {}
            _ => {
                let mut name = String::new();
                let quoted = c == '\'' || c == '"';
                if !quoted {
                    name.push(c);
                }
                let quote = c;
                while let Some(next) = chars.peek() {
                    if quoted {
                        if *next == quote {
                            chars.next();
                            if chars.peek() == Some(&quote) {
                                name.push(quote);
                                chars.next();
                                continue;
                            }
                            break;
                        }
                        name.push(*next);
                        chars.next();
                    } else if matches!(*next, '(' | ')' | ',' | ':' | ';' | '[') {
                        break;
                    } else {
                        name.push(*next);
                        chars.next();
                    }
                }
                let name = name.trim().to_string();

                match current {
                    Some(index) if !nodes[index].is_leaf() => match name.parse::<f64>() {
                        Ok(support) => nodes[index].support = Some(support),
                        Err(_) => nodes[index].name = Some(name),
                    },
                    _ => {
                        let parent = stack.last().copied();
                        if parent.is_none() && !nodes.is_empty() {
                            return Err(Error::InvalidNewick {
                                reason: "more than one root",
                            });
                        }
                        let index = add_node(&mut nodes, &mut annotations, parent);
                        nodes[index].name = Some(name);
                        current = Some(index);
                    }
                }
            }
        }
    }

    if !stack.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "unbalanced parentheses",
        });
    }
    if nodes.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "empty tree",
        });
    }
    Ok(Tree {
        nodes,
        root: 0,
        annotations,
        tree_annotations,
        rooted,
    })
}

fn add_node(
    nodes: &mut Vec<Clade>,
    annotations: &mut Vec<Annotations>,
    parent: Option<usize>,
) -> usize {
    nodes.push(Clade {
        name: None,
        branch_length: None,
        support: None,
        children: Vec::new(),
        parent,
    });
    annotations.push(Annotations::new());
    let index = nodes.len() - 1;
    if let Some(parent) = parent {
        nodes[parent].children.push(index);
    }
    index
}

fn parse_comment(comment: &str) -> (Option<bool>, Annotations) {
    let text = comment.trim();
    if text.eq_ignore_ascii_case("&R") {
        return (Some(true), Annotations::new());
    }
    if text.eq_ignore_ascii_case("&U") {
        return (Some(false), Annotations::new());
    }

    let mut annotations = Annotations::new();
    if let Some(body) = text.strip_prefix("&&NHX:") {
        for field in split_delimited(body, ':') {
            insert_annotation(&mut annotations, &field, '=');
        }
    } else if let Some(body) = text.strip_prefix('&') {
        for field in split_delimited(body, ',') {
            insert_annotation(&mut annotations, &field, '=');
        }
    } else if !text.is_empty() {
        annotations.insert(
            "comment".to_string(),
            AnnotationValue::Text(text.to_string()),
        );
    }
    (None, annotations)
}

fn insert_annotation(annotations: &mut Annotations, field: &str, separator: char) {
    let field = field.trim();
    if field.is_empty() {
        return;
    }
    let (key, value) = field
        .split_once(separator)
        .map_or((field, "true"), |(key, value)| (key.trim(), value.trim()));
    if !key.is_empty() {
        annotations.insert(key.to_string(), parse_annotation_value(value));
    }
}

fn parse_annotation_value(value: &str) -> AnnotationValue {
    let value = value.trim();
    if value.len() >= 2 {
        let first = value.as_bytes()[0];
        let last = value.as_bytes()[value.len() - 1];
        if (first == b'\'' && last == b'\'') || (first == b'"' && last == b'"') {
            return AnnotationValue::Text(value[1..value.len() - 1].to_string());
        }
        if first == b'{' && last == b'}' {
            return AnnotationValue::List(
                split_delimited(&value[1..value.len() - 1], ',')
                    .into_iter()
                    .map(|item| parse_annotation_value(&item))
                    .collect(),
            );
        }
    }
    if value.eq_ignore_ascii_case("true") {
        return AnnotationValue::Boolean(true);
    }
    if value.eq_ignore_ascii_case("false") {
        return AnnotationValue::Boolean(false);
    }
    if let Ok(number) = value.parse::<f64>() {
        if number.is_finite() {
            return AnnotationValue::Number(number);
        }
    }
    AnnotationValue::Text(value.to_string())
}

fn split_delimited(input: &str, delimiter: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quote = None;
    let mut braces = 0usize;
    for character in input.chars() {
        if let Some(active) = quote {
            field.push(character);
            if character == active {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                field.push(character);
            }
            '{' => {
                braces += 1;
                field.push(character);
            }
            '}' => {
                braces = braces.saturating_sub(1);
                field.push(character);
            }
            value if value == delimiter && braces == 0 => {
                fields.push(field.trim().to_string());
                field.clear();
            }
            _ => field.push(character),
        }
    }
    if !field.trim().is_empty() {
        fields.push(field.trim().to_string());
    }
    fields
}

fn parse_nexus(input: &str) -> Result<Tree, Error> {
    let statements = nexus_statements(input);
    let mut translation = BTreeMap::new();
    let mut expression = None;

    for statement in &statements {
        let trimmed = statement.trim();
        let lower = trimmed.to_ascii_lowercase();
        if lower.starts_with("translate") {
            let body = trimmed.get("translate".len()..).unwrap_or_default().trim();
            for entry in split_delimited(body, ',') {
                let split = entry.find(char::is_whitespace).ok_or(Error::InvalidNexus {
                    reason: "a translate entry has no taxon name",
                })?;
                let key = entry[..split].trim();
                let name = unquote(entry[split..].trim());
                if key.is_empty() || name.is_empty() {
                    return Err(Error::InvalidNexus {
                        reason: "an empty translate entry",
                    });
                }
                translation.insert(key.to_string(), name);
            }
        } else if lower.starts_with("tree ") || lower.starts_with("utree ") {
            expression = trimmed
                .split_once('=')
                .map(|(_, tree)| tree.trim().to_string());
            if expression.is_none() {
                return Err(Error::InvalidNexus {
                    reason: "a tree statement has no equals sign",
                });
            }
            break;
        }
    }

    let expression = expression.ok_or(Error::InvalidNexus {
        reason: "no tree statement",
    })?;
    let mut tree = Tree::parse_annotated_newick(&expression)?;
    for leaf in tree.leaves() {
        let Some(name) = tree.nodes[leaf].name.as_deref() else {
            continue;
        };
        if let Some(translated) = translation.get(name) {
            tree.nodes[leaf].name = Some(translated.clone());
        }
    }
    Ok(tree)
}

fn nexus_statements(input: &str) -> Vec<String> {
    let mut statements = Vec::new();
    let mut statement = String::new();
    let mut quote = None;
    let mut bracket_depth = 0usize;
    for character in input.chars() {
        if let Some(active) = quote {
            statement.push(character);
            if character == active {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                statement.push(character);
            }
            '[' => {
                bracket_depth += 1;
                statement.push(character);
            }
            ']' => {
                bracket_depth = bracket_depth.saturating_sub(1);
                statement.push(character);
            }
            ';' if bracket_depth == 0 => {
                if !statement.trim().is_empty() {
                    statements.push(statement.trim().to_string());
                }
                statement.clear();
            }
            _ => statement.push(character),
        }
    }
    if !statement.trim().is_empty() {
        statements.push(statement.trim().to_string());
    }
    statements
}

fn unquote(value: &str) -> String {
    let value = value.trim();
    if value.len() >= 2 {
        let first = value.as_bytes()[0];
        let last = value.as_bytes()[value.len() - 1];
        if (first == b'\'' && last == b'\'') || (first == b'"' && last == b'"') {
            return value[1..value.len() - 1].replace("''", "'");
        }
    }
    value.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flat_tree_parses() {
        let tree = Tree::parse_newick("(A,B,C);").unwrap();
        assert_eq!(tree.leaf_names(), ["A", "B", "C"]);
        assert_eq!(tree.leaf_count(), 3);
        assert_eq!(tree.nodes().len(), 4, "three leaves and a root");
    }

    #[test]
    fn nesting_and_branch_lengths_parse() {
        let tree = Tree::parse_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
        assert_eq!(tree.leaf_names(), ["A", "B", "C"]);
        let a = tree
            .nodes()
            .iter()
            .find(|node| node.name.as_deref() == Some("A"))
            .unwrap();
        assert_eq!(a.branch_length, Some(0.1));
    }

    #[test]
    fn an_internal_number_is_support_and_a_word_is_a_name() {
        let tree = Tree::parse_newick("((A:0.1,B:0.2)0.98:0.3,C:0.4);").unwrap();
        let internal = tree
            .nodes()
            .iter()
            .find(|node| !node.is_leaf() && node.support.is_some())
            .unwrap();
        assert_eq!(internal.support, Some(0.98));
        assert_eq!(internal.branch_length, Some(0.3));

        let named = Tree::parse_newick("((A,B)clade_one,C);").unwrap();
        assert!(named
            .nodes()
            .iter()
            .any(|node| node.name.as_deref() == Some("clade_one") && !node.is_leaf()));
    }

    #[test]
    fn quoted_names_keep_their_punctuation() {
        let tree = Tree::parse_newick("('ERR (one)':0.1,'B,C':0.2);").unwrap();
        assert_eq!(tree.leaf_names(), ["ERR (one)", "B,C"]);
    }

    #[test]
    fn a_doubled_quote_is_one_literal_quote_and_not_a_second_taxon() {
        let tree = Tree::parse_newick("('O''Brien':0.1,B:0.2);").unwrap();
        assert_eq!(tree.leaf_names(), ["O'Brien", "B"]);
        assert_eq!(tree.leaf_count(), 2);
        let named = tree.leaves();
        assert_eq!(tree.nodes()[named[0]].branch_length, Some(0.1));
    }

    #[test]
    fn an_empty_label_is_an_unnamed_leaf_rather_than_a_dropped_tip() {
        // The canonical example of the format is four unnamed leaves.
        let tree = Tree::parse_newick("(,,(,));").unwrap();
        assert_eq!(tree.leaf_count(), 4);
        assert_eq!(tree.nodes().len(), 6, "four leaves, two clades");
        assert_eq!(tree.leaf_names(), ["", "", "", ""]);

        assert_eq!(Tree::parse_newick("((,),(,));").unwrap().leaf_count(), 4);
        assert_eq!(Tree::parse_newick("(A,B,);").unwrap().leaf_count(), 3);
        assert_eq!(Tree::parse_newick("(A,B,);").unwrap().nodes().len(), 4);
    }

    #[test]
    fn a_branch_length_on_an_empty_label_lands_on_that_leaf() {
        let tree = Tree::parse_newick("(:0.1,:0.2,(:0.3,:0.4):0.5);").unwrap();
        assert_eq!(tree.nodes().len(), 6);
        assert_eq!(tree.leaf_count(), 4);
        // Root first, then the five branches in the order the file wrote them.
        let lengths: Vec<Option<f64>> =
            tree.nodes().iter().map(|node| node.branch_length).collect();
        assert_eq!(
            lengths,
            vec![None, Some(0.1), Some(0.2), Some(0.5), Some(0.3), Some(0.4)]
        );
        // The clade keeps its own 0.5 and its children keep theirs.
        assert_eq!(tree.nodes()[3].children, vec![4, 5]);
        assert!((tree.max_depth(false) - 0.9).abs() < 1e-12, "0.5 then 0.4");
    }

    #[test]
    fn scientific_notation_in_a_branch_length_parses() {
        let tree = Tree::parse_newick("(A:1.5e-4,B:2E-3);").unwrap();
        let lengths: Vec<Option<f64>> = tree
            .leaves()
            .iter()
            .map(|index| tree.nodes()[*index].branch_length)
            .collect();
        assert_eq!(lengths, vec![Some(1.5e-4), Some(2e-3)]);
    }

    #[test]
    fn comments_are_skipped_wherever_a_real_file_puts_them() {
        // RAxML and BEAST write a rootedness marker before the tree and
        // annotations inside it. Read as names, the first one alone makes the
        // whole file fail with "more than one root".
        let plain = Tree::parse_newick("((A:0.1,B:0.2)0.98:0.3,C:0.4);").unwrap();
        let annotated =
            Tree::parse_newick("[&R] ((A[&rate=1.2]:0.1,B:0.2)0.98:0.3[&height=0.4],C:0.4);")
                .unwrap();
        assert_eq!(annotated.leaf_names(), plain.leaf_names());
        assert_eq!(annotated.max_depth(false), plain.max_depth(false));

        // An unclosed comment runs out rather than eating the tree.
        assert_eq!(Tree::parse_newick("(A,B)[oops;").unwrap().leaf_count(), 2);
    }

    #[test]
    fn annotated_newick_keeps_beast_values_and_rootedness() {
        let tree = Tree::parse_annotated_newick(
            "[&R] (A[&date=2024.5,location='Lima',flags={1,2},selected=true]:0.1,B:0.2);",
        )
        .unwrap();
        let a = tree.node_named("A").unwrap();
        assert_eq!(tree.rooted(), Some(true));
        assert_eq!(
            tree.annotation(a, "date")
                .and_then(AnnotationValue::as_number),
            Some(2024.5)
        );
        assert_eq!(
            tree.annotation(a, "location")
                .and_then(AnnotationValue::as_text),
            Some("Lima")
        );
        assert_eq!(
            tree.annotation(a, "selected")
                .and_then(AnnotationValue::as_bool),
            Some(true)
        );
        assert!(matches!(
            tree.annotation(a, "flags"),
            Some(AnnotationValue::List(values)) if values.len() == 2
        ));
    }

    #[test]
    fn annotations_can_be_added_after_parsing() {
        let mut tree = Tree::parse_newick("(A,B);").unwrap();
        let a = tree.node_named("A").unwrap();
        tree.annotations_mut(a)
            .unwrap()
            .insert("country".into(), AnnotationValue::Text("Peru".into()));
        tree.tree_annotations_mut()
            .insert("clock".into(), AnnotationValue::Boolean(true));
        assert_eq!(
            tree.annotation(a, "country")
                .and_then(AnnotationValue::as_text),
            Some("Peru")
        );
        assert_eq!(
            tree.tree_annotations()
                .get("clock")
                .and_then(AnnotationValue::as_bool),
            Some(true)
        );
        assert!(tree.annotations_mut(99).is_none());
    }

    #[test]
    fn nhx_annotations_are_typed_too() {
        let tree = Tree::parse_annotated_newick("(A[&&NHX:S=human:B=95],B);").unwrap();
        let a = tree.node_named("A").unwrap();
        assert_eq!(
            tree.annotation(a, "S").and_then(AnnotationValue::as_text),
            Some("human")
        );
        assert_eq!(
            tree.annotation(a, "B").and_then(AnnotationValue::as_number),
            Some(95.0)
        );
    }

    #[test]
    fn compatibility_newick_still_discards_annotations() {
        let tree = Tree::parse_newick("[&R] (A[&date=2024.5],B);").unwrap();
        assert_eq!(tree.rooted(), None);
        assert!(tree
            .annotations(tree.node_named("A").unwrap())
            .unwrap()
            .is_empty());
    }

    #[test]
    fn nexus_translation_and_annotations_reach_the_tree() {
        let nexus = "#NEXUS\nBegin trees;\nTranslate 1 'sample A', 2 sample_B;\n\
                     Tree outbreak = [&R] (1[&country=Peru]:0.1,2:0.2);\nEnd;";
        let tree = Tree::parse_nexus(nexus).unwrap();
        assert_eq!(tree.leaf_names(), ["sample A", "sample_B"]);
        assert_eq!(tree.rooted(), Some(true));
        let a = tree.node_named("sample A").unwrap();
        assert_eq!(tree.annotation(a, "country").unwrap().to_string(), "Peru");
    }

    #[test]
    fn nexus_without_a_tree_says_what_is_missing() {
        assert!(matches!(
            Tree::parse_nexus("#NEXUS\nBegin taxa; End;"),
            Err(Error::InvalidNexus {
                reason: "no tree statement"
            })
        ));
    }

    #[test]
    fn a_semicolon_is_optional_and_whitespace_is_ignored() {
        let with = Tree::parse_newick("(A:0.1,B:0.2);").unwrap();
        let without = Tree::parse_newick("  (A:0.1,\n B:0.2)  ").unwrap();
        assert_eq!(with, without);
    }

    #[test]
    fn malformed_newick_is_rejected_rather_than_guessed() {
        for bad in ["", ";", "((A,B)", "(A,B))", ",A", "(A:x,B);"] {
            assert!(Tree::parse_newick(bad).is_err(), "{bad:?} should not parse");
        }
    }

    #[test]
    fn leaves_come_back_in_the_order_the_file_meant() {
        let tree = Tree::parse_newick("((D,C),(B,A));").unwrap();
        assert_eq!(tree.leaf_names(), ["D", "C", "B", "A"]);
    }

    #[test]
    fn a_phylogram_places_leaves_at_their_distance_from_the_root() {
        let tree = Tree::parse_newick("((A:0.1,B:0.2):0.3,C:0.05);").unwrap();
        let layout = tree.layout(false);
        let depth_of = |name: &str| {
            let index = tree
                .nodes()
                .iter()
                .position(|node| node.name.as_deref() == Some(name))
                .unwrap();
            layout.iter().find(|p| p.node == index).unwrap().depth
        };
        assert!((depth_of("A") - 0.4).abs() < 1e-12, "0.3 then 0.1");
        assert!((depth_of("B") - 0.5).abs() < 1e-12);
        assert!((depth_of("C") - 0.05).abs() < 1e-12);
    }

    #[test]
    fn a_time_layout_uses_tip_dates_and_infers_internal_dates() {
        let tree = Tree::parse_annotated_newick(
            "((A[&date=2024.0]:2,B[&date=2025.0]:3)AB:1,C[&date=2023.0]:4);",
        )
        .unwrap();
        let layout = tree.time_layout("date", TimeDirection::Increasing).unwrap();
        let value = |name: &str| {
            let node = tree.node_named(name).unwrap();
            layout[node].depth
        };
        assert_eq!(value("A"), 2024.0);
        assert_eq!(value("B"), 2025.0);
        assert_eq!(value("AB"), 2022.0, "mean of 2024-2 and 2025-3");
    }

    #[test]
    fn a_time_layout_refuses_a_tip_without_the_requested_value() {
        let tree = Tree::parse_annotated_newick("(A[&date=2024]:1,B:1);").unwrap();
        assert!(tree
            .time_layout("date", TimeDirection::Increasing)
            .is_none());
    }

    #[test]
    fn a_cladogram_counts_branches_instead_of_measuring_them() {
        let tree = Tree::parse_newick("((A:0.1,B:0.2):0.3,C:99.0);").unwrap();
        let layout = tree.layout(true);
        let depth_of = |name: &str| {
            let index = tree
                .nodes()
                .iter()
                .position(|node| node.name.as_deref() == Some(name))
                .unwrap();
            layout.iter().find(|p| p.node == index).unwrap().depth
        };
        assert_eq!(depth_of("A"), 2.0);
        assert_eq!(depth_of("C"), 1.0, "one branch, however long it is");
    }

    #[test]
    fn a_parent_sits_between_its_children() {
        let tree = Tree::parse_newick("((A,B),C);").unwrap();
        let layout = tree.layout(true);
        let row_of = |index: usize| layout.iter().find(|p| p.node == index).unwrap().row;
        let leaves = tree.leaves();
        // A and B are rows 0 and 1, so their parent is at 0.5.
        assert_eq!(row_of(leaves[0]), 0.0);
        assert_eq!(row_of(leaves[1]), 1.0);
        let parent = tree.nodes()[leaves[0]].parent.unwrap();
        assert_eq!(row_of(parent), 0.5);
    }

    #[test]
    fn a_missing_branch_length_counts_as_zero_rather_than_breaking() {
        let tree = Tree::parse_newick("((A,B):0.3,C:0.1);").unwrap();
        let layout = tree.layout(false);
        assert!(layout.iter().all(|p| p.depth.is_finite()));
        assert!((tree.max_depth(false) - 0.3).abs() < 1e-12);
    }

    #[test]
    fn a_single_leaf_is_a_tree() {
        let tree = Tree::parse_newick("A;").unwrap();
        assert_eq!(tree.leaf_names(), ["A"]);
        assert_eq!(tree.max_depth(true), 0.0);
    }

    #[test]
    fn ancestors_descendants_and_mrca_agree() {
        let tree = Tree::parse_newick("(((A,B)AB,C)ABC,D);").unwrap();
        let a = tree.node_named("A").unwrap();
        let b = tree.node_named("B").unwrap();
        let c = tree.node_named("C").unwrap();
        let ab = tree.node_named("AB").unwrap();
        let abc = tree.node_named("ABC").unwrap();
        assert_eq!(tree.mrca(&[a, b]), Some(ab));
        assert_eq!(tree.mrca(&[a, c]), Some(abc));
        assert!(tree.ancestors(a).contains(&ab));
        assert!(tree.descendants(abc).contains(&c));
        assert_eq!(tree.clade_size(abc), 3);
        assert_eq!(tree.mrca(&[]), None);
    }

    #[test]
    fn rotating_and_ladderizing_change_only_leaf_order() {
        let mut tree = Tree::parse_newick("((A,B,C)large,D);").unwrap();
        let large = tree.node_named("large").unwrap();
        assert!(tree.rotate(large));
        assert_eq!(tree.leaf_names(), ["C", "B", "A", "D"]);
        tree.ladderize(false);
        assert_eq!(tree.leaf_names()[0], "D");
        let mut sorted = tree.leaf_names();
        sorted.sort();
        assert_eq!(sorted, ["A", "B", "C", "D"]);
    }

    #[test]
    fn rerooting_preserves_tips_and_pairwise_distance() {
        let mut tree = Tree::parse_newick("((A:1,B:2)AB:3,(C:4,D:5)CD:6);").unwrap();
        let a = tree.node_named("A").unwrap();
        let c = tree.node_named("C").unwrap();
        let before = pair_distance(&tree, a, c);
        let ab = tree.node_named("AB").unwrap();
        assert!(tree.reroot(ab));
        assert_eq!(tree.root(), ab);
        assert_eq!(tree.rooted(), Some(true));
        assert_eq!(pair_distance(&tree, a, c), before);
        let mut leaves = tree.leaf_names();
        leaves.sort();
        assert_eq!(leaves, ["A", "B", "C", "D"]);
        assert!(!tree.reroot(a), "a sampled tip stays a sampled tip");
    }

    #[test]
    fn rerooting_keeps_support_on_its_undirected_edge() {
        let mut tree = Tree::parse_newick("((A:1,B:1)0.99:2,C:3);").unwrap();
        let supported = tree
            .nodes()
            .iter()
            .position(|node| node.support == Some(0.99))
            .unwrap();
        assert!(tree.reroot(supported));
        assert_eq!(tree.nodes()[tree.root()].support, None);
        assert_eq!(
            tree.nodes()
                .iter()
                .filter(|node| node.support == Some(0.99))
                .count(),
            1
        );
    }

    #[test]
    fn a_monophyletic_outgroup_gets_a_new_root_without_changing_distances() {
        let mut tree = Tree::parse_newick("(((A:1,B:1)AB:2,C:3)ING:4,(O1:2,O2:2)OUT:5);").unwrap();
        let a = tree.node_named("A").unwrap();
        let o1 = tree.node_named("O1").unwrap();
        let o2 = tree.node_named("O2").unwrap();
        let out = tree.node_named("OUT").unwrap();
        let before = pair_distance(&tree, a, o1);
        let root = tree.reroot_outgroup(&[o1, o2]).unwrap();
        assert_eq!(tree.root(), root);
        assert_eq!(tree.rooted(), Some(true));
        assert!(tree.nodes()[root].children.contains(&out));
        assert_eq!(pair_distance(&tree, a, o1), before);
        assert_eq!(tree.leaf_count(), 5);
    }

    #[test]
    fn a_non_monophyletic_outgroup_is_refused_without_mutation() {
        let mut tree = Tree::parse_newick("(((A:1,B:1)AB:2,C:3)ING:4,(O1:2,O2:2)OUT:5);").unwrap();
        let before = tree.clone();
        let a = tree.node_named("A").unwrap();
        let o1 = tree.node_named("O1").unwrap();
        assert_eq!(tree.reroot_outgroup(&[a, o1]), None);
        assert_eq!(tree, before);
    }

    #[test]
    fn midpoint_rooting_splits_the_diameter_edge_exactly() {
        let mut tree = Tree::parse_newick("((A:1,B:1)AB:1,C:4);").unwrap();
        let a = tree.node_named("A").unwrap();
        let b = tree.node_named("B").unwrap();
        let c = tree.node_named("C").unwrap();
        let ac = pair_distance(&tree, a, c);
        let bc = pair_distance(&tree, b, c);
        let old_nodes = tree.nodes().len();
        let root = tree.reroot_midpoint().unwrap();
        assert_eq!(root, old_nodes, "the midpoint lies inside the C edge");
        let layout = tree.layout(false);
        let depth = |node: usize| layout.iter().find(|p| p.node == node).unwrap().depth;
        assert!((depth(a) - 3.0).abs() < 1e-12);
        assert!((depth(c) - 3.0).abs() < 1e-12);
        assert_eq!(pair_distance(&tree, a, c), ac);
        assert_eq!(pair_distance(&tree, b, c), bc);
    }

    #[test]
    fn midpoint_rooting_requires_complete_non_negative_lengths() {
        for input in ["((A:1,B)AB:1,C:4);", "((A:1,B:-1)AB:1,C:4);"] {
            let mut tree = Tree::parse_newick(input).unwrap();
            let before = tree.clone();
            assert_eq!(tree.reroot_midpoint(), None);
            assert_eq!(tree, before);
        }
    }

    #[test]
    fn a_subtree_remaps_nodes_and_keeps_annotations() {
        let tree =
            Tree::parse_annotated_newick("((A[&country=Peru]:1,B[&country=Chile]:1)AB:2,C:3);")
                .unwrap();
        let ab = tree.node_named("AB").unwrap();
        let subtree = tree.subtree(ab).unwrap();
        assert_eq!(subtree.root(), 0);
        assert_eq!(subtree.leaf_names(), ["A", "B"]);
        assert_eq!(subtree.nodes()[0].branch_length, None);
        let a = subtree.node_named("A").unwrap();
        assert_eq!(
            subtree.annotation(a, "country").unwrap().to_string(),
            "Peru"
        );
    }

    #[test]
    fn collapsing_a_clade_removes_only_its_descendants() {
        let mut tree = Tree::parse_newick("((A,B)AB,C);").unwrap();
        let ab = tree.node_named("AB").unwrap();
        assert!(tree.collapse(ab));
        assert_eq!(tree.leaf_names(), ["AB", "C"]);
        assert_eq!(tree.nodes().len(), 3, "root and two terminal clades");
        assert!(!tree.collapse(tree.node_named("AB").unwrap()));
    }

    fn pair_distance(tree: &Tree, left: usize, right: usize) -> f64 {
        let layout = tree.layout(false);
        let depth = |node: usize| layout.iter().find(|p| p.node == node).unwrap().depth;
        let ancestor = tree.mrca(&[left, right]).unwrap();
        depth(left) + depth(right) - 2.0 * depth(ancestor)
    }

    #[test]
    fn a_deep_ladder_does_not_overflow_the_stack() {
        // Ten thousand nested clades, which a recursive walk would not survive.
        let mut newick = String::new();
        for _ in 0..10_000 {
            newick.push('(');
        }
        newick.push('A');
        for index in 0..10_000 {
            newick.push_str(&format!(",L{index})"));
        }
        newick.push(';');
        let tree = Tree::parse_newick(&newick).unwrap();
        assert_eq!(tree.leaf_count(), 10_001);
        assert_eq!(tree.layout(true).len(), tree.nodes().len());
    }
}
