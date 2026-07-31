//! Phylogenies: reading them, and working out where their branches go.
//!
//! A tree on its own is a tree, and there are viewers for those. What this one
//! is for is ordering: its leaf order is the row order of a
//! [`SnpTrack`](crate::SnpTrack), so the panel of variable sites beside it is
//! sorted by descent rather than by sample name, and a clade's shared
//! substitutions line up into a block you can see.

use crate::error::Error;

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
}

impl Tree {
    /// Reads one Newick string.
    ///
    /// Handles nested clades, branch lengths, quoted names and internal labels,
    /// reading an internal label as a support value when it parses as a number
    /// and keeping it as a name when it does not. Trailing semicolon optional.
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
        let text = input.trim().trim_end_matches(';').trim();
        if text.is_empty() {
            return Err(Error::InvalidNewick {
                reason: "empty tree",
            });
        }

        let mut nodes: Vec<Clade> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();
        let mut current: Option<usize> = None;
        let mut chars = text.chars().peekable();

        let new_node = |nodes: &mut Vec<Clade>, parent: Option<usize>| -> usize {
            nodes.push(Clade {
                name: None,
                branch_length: None,
                support: None,
                children: Vec::new(),
                parent,
            });
            let index = nodes.len() - 1;
            if let Some(parent) = parent {
                nodes[parent].children.push(index);
            }
            index
        };

        while let Some(c) = chars.next() {
            match c {
                '(' => {
                    let parent = stack.last().copied();
                    if parent.is_none() && !nodes.is_empty() {
                        return Err(Error::InvalidNewick {
                            reason: "more than one root",
                        });
                    }
                    let index = new_node(&mut nodes, parent);
                    stack.push(index);
                    current = None;
                }
                ')' => {
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
                    let target =
                        current
                            .or_else(|| stack.last().copied())
                            .ok_or(Error::InvalidNewick {
                                reason: "branch length with nothing to attach to",
                            })?;
                    nodes[target].branch_length = Some(length);
                }
                '[' => {
                    // A comment. Newick does not nest them, so the first
                    // closing bracket ends it; an unclosed one runs to the end
                    // rather than swallowing the tree as a name.
                    for next in chars.by_ref() {
                        if next == ']' {
                            break;
                        }
                    }
                }
                c if c.is_whitespace() => {}
                _ => {
                    // A name: either a leaf being introduced, or the label of
                    // the clade that just closed.
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
                        // A label following a closed clade is its support or
                        // its name, not a new leaf.
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
                            let index = new_node(&mut nodes, parent);
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
        Ok(Tree { nodes, root: 0 })
    }

    /// The nodes, in the order they were parsed.
    pub fn nodes(&self) -> &[Clade] {
        &self.nodes
    }

    /// Index of the root.
    pub fn root(&self) -> usize {
        self.root
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
        self.nodes.iter().filter(|node| node.is_leaf()).count()
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

    /// The deepest leaf, in whatever units the layout used.
    pub fn max_depth(&self, cladogram: bool) -> f64 {
        self.layout(cladogram)
            .iter()
            .filter(|placement| self.nodes[placement.node].is_leaf())
            .map(|placement| placement.depth)
            .fold(0.0f64, f64::max)
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
