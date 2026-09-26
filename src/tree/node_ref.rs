//! A clade named the ways a reader names one.

use std::fmt;

use super::{AnnotationValue, Tree};

/// A node of a tree, named the ways a reader names one: by its place in
/// [`Tree::nodes`], by a node's name, as the smallest clade holding some tips,
/// or as the clade of the tips that carry one value.
///
/// Every builder of a [`TreeTrack`](crate::TreeTrack) that picks a clade takes
/// one, and an index or a name still reads as it did, since a `usize` and a
/// `&str` are each a `NodeRef`. A clade picked by its tips or by a value is
/// found in the tree as it stands when the builder is called, so a clade
/// folded before a reroot is still the clade that holds those tips after it.
///
/// ```
/// use karyon::{NodeRef, Tree, TreeTrack};
///
/// let tree = Tree::parse_annotated_newick(
///     "((A[&lineage=L4]:1,B[&lineage=L4]:1):1,C[&lineage=L2]:2);",
/// )?;
/// let track = TreeTrack::new(tree)
///     .collapse(NodeRef::holding("lineage", "L4"))
///     .clade_highlight(karyon::CladeHighlight::new(NodeRef::mrca(["A", "B"])));
/// assert!(track.warnings().is_empty(), "{:?}", track.warnings());
/// # Ok::<(), karyon::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NodeRef {
    /// A node by its index in [`Tree::nodes`].
    Index(usize),
    /// The node with this name: a labelled clade, or a tip.
    Named(String),
    /// The smallest clade holding every one of these tips or named clades.
    Mrca(Vec<String>),
    /// The smallest clade holding every tip whose annotation `key` reads
    /// `value`, from the tip itself or, where it has none, from its nearest
    /// ancestor that has, as a strip reads it.
    Holding {
        /// The annotation.
        key: String,
        /// Its value, as it is written.
        value: String,
    },
}

/// Where a [`NodeRef`] lands in a tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Found {
    /// The node.
    pub node: usize,
    /// The tips the clade holds that were not asked for, said in a clause,
    /// where there are any: the smallest clade holding two tips can hold a
    /// hundred more, and a fold or a highlight over it would say otherwise.
    pub also: Option<String>,
}

impl NodeRef {
    /// The node with this name: a labelled clade, or a tip.
    pub fn named(name: impl Into<String>) -> Self {
        NodeRef::Named(name.into())
    }

    /// The smallest clade holding every one of these tips or named clades, as
    /// ggtree's `MRCA` or iTOL's `A|B` name one.
    pub fn mrca<I, S>(tips: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        NodeRef::Mrca(tips.into_iter().map(Into::into).collect())
    }

    /// The smallest clade holding every tip whose annotation `key` reads
    /// `value`, as `holding("lineage", "L4")` names the clade of lineage L4.
    pub fn holding(key: impl Into<String>, value: impl Into<String>) -> Self {
        NodeRef::Holding {
            key: key.into(),
            value: value.into(),
        }
    }

    /// Finds the node in `tree`.
    ///
    /// # Errors
    ///
    /// Says why there is none, in a clause: an index the tree does not have, a
    /// name no node carries, or a value no tip carries.
    pub fn find(&self, tree: &Tree) -> Result<Found, String> {
        match self {
            NodeRef::Index(node) => match tree.nodes().get(*node) {
                Some(_) => Ok(Found {
                    node: *node,
                    also: None,
                }),
                None => Err(format!("the tree has no node {node}")),
            },
            NodeRef::Named(name) => match tree.node_named(name) {
                Some(node) => Ok(Found { node, also: None }),
                None => Err(format!("no node is named {name}")),
            },
            NodeRef::Mrca(names) => {
                if names.is_empty() {
                    return Err("no tips were named".to_string());
                }
                let mut nodes = Vec::with_capacity(names.len());
                for name in names {
                    match tree.node_named(name) {
                        Some(node) => nodes.push(node),
                        None => return Err(format!("no node is named {name}")),
                    }
                }
                let node = tree
                    .mrca(&nodes)
                    .ok_or_else(|| "the named nodes share no ancestor".to_string())?;
                // Every tip under the named nodes was asked for; any other tip
                // the clade holds was not.
                let asked: std::collections::BTreeSet<usize> = nodes
                    .iter()
                    .flat_map(|named| tips_under(tree, *named))
                    .collect();
                let others: Vec<usize> = tips_under(tree, node)
                    .into_iter()
                    .filter(|tip| !asked.contains(tip))
                    .collect();
                Ok(Found {
                    node,
                    also: also_holds(tree, &others, ""),
                })
            }
            NodeRef::Holding { key, value } => {
                let carries = |tip: usize| {
                    inherited(tree, tip, key).is_some_and(|held| held.to_string() == *value)
                };
                let carriers: Vec<usize> = tree
                    .leaves()
                    .into_iter()
                    .filter(|tip| carries(*tip))
                    .collect();
                if carriers.is_empty() {
                    return Err(format!("no tip has {key} {value}"));
                }
                let node = tree
                    .mrca(&carriers)
                    .ok_or_else(|| format!("the tips with {key} {value} share no ancestor"))?;
                let others: Vec<usize> = tips_under(tree, node)
                    .into_iter()
                    .filter(|tip| !carries(*tip))
                    .collect();
                Ok(Found {
                    node,
                    also: also_holds(tree, &others, &format!(" without {key} {value}")),
                })
            }
        }
    }
}

/// The tips at or under `node`.
fn tips_under(tree: &Tree, node: usize) -> Vec<usize> {
    if tree.nodes()[node].is_leaf() {
        return vec![node];
    }
    tree.descendants(node)
        .into_iter()
        .filter(|index| tree.nodes()[*index].is_leaf())
        .collect()
}

/// A tip's value of `key`: its own, or its nearest ancestor's.
fn inherited<'a>(tree: &'a Tree, node: usize, key: &str) -> Option<&'a AnnotationValue> {
    tree.annotation(node, key).or_else(|| {
        tree.ancestors(node)
            .into_iter()
            .find_map(|ancestor| tree.annotation(ancestor, key))
    })
}

/// "also holds A, B and 3 more", for the tips a clade holds that were not
/// asked for, or `None` where there are none.
fn also_holds(tree: &Tree, others: &[usize], what: &str) -> Option<String> {
    if others.is_empty() {
        return None;
    }
    const SHOWN: usize = 3;
    let names: Vec<String> = others
        .iter()
        .take(SHOWN)
        .map(|tip| {
            tree.nodes()[*tip]
                .name
                .clone()
                .unwrap_or_else(|| format!("node {tip}"))
        })
        .collect();
    let listed = match (names.as_slice(), others.len()) {
        ([one], 1) => one.clone(),
        (_, count) if count <= SHOWN => {
            let (last, rest) = names.split_last().expect("at least two names");
            format!("{} and {last}", rest.join(", "))
        }
        (_, count) => format!("{} and {} more", names.join(", "), count - SHOWN),
    };
    let tips = if others.len() == 1 { "tip" } else { "tips" };
    Some(format!(
        "also holds {} {tips}{what}: {listed}",
        others.len()
    ))
}

impl fmt::Display for NodeRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NodeRef::Index(node) => write!(f, "node {node}"),
            NodeRef::Named(name) => f.write_str(name),
            NodeRef::Mrca(names) => write!(f, "the clade of {}", names.join(" and ")),
            NodeRef::Holding { key, value } => write!(f, "the clade of {key} {value}"),
        }
    }
}

impl From<usize> for NodeRef {
    fn from(node: usize) -> Self {
        NodeRef::Index(node)
    }
}

impl From<&str> for NodeRef {
    fn from(name: &str) -> Self {
        NodeRef::Named(name.to_string())
    }
}

impl From<String> for NodeRef {
    fn from(name: String) -> Self {
        NodeRef::Named(name)
    }
}

impl From<&String> for NodeRef {
    fn from(name: &String) -> Self {
        NodeRef::Named(name.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree() -> Tree {
        Tree::parse_annotated_newick(
            "(((A[&lineage=L4]:1,B[&lineage=L4]:1)AB:1,C[&lineage=L2]:1):1,\
             (D:1,E:1)[&lineage=L1]:1);",
        )
        .unwrap()
    }

    #[test]
    fn a_clade_is_found_by_its_index_its_name_its_tips_or_a_value() {
        let tree = tree();
        let ab = tree.node_named("AB").unwrap();
        assert_eq!(NodeRef::from(ab).find(&tree).unwrap().node, ab);
        assert_eq!(NodeRef::named("AB").find(&tree).unwrap().node, ab);
        assert_eq!(NodeRef::mrca(["A", "B"]).find(&tree).unwrap().node, ab);
        assert_eq!(
            NodeRef::holding("lineage", "L4").find(&tree).unwrap().node,
            ab
        );
        // Inherited from the clade above the tips, as a strip reads it.
        let de = NodeRef::holding("lineage", "L1").find(&tree).unwrap();
        assert_eq!(tree.nodes()[de.node].children.len(), 2);
        assert_eq!(de.also, None);
    }

    #[test]
    fn a_clade_holding_tips_that_were_not_asked_for_says_which() {
        let tree = tree();
        let found = NodeRef::mrca(["A", "C"]).find(&tree).unwrap();
        assert_eq!(found.also.as_deref(), Some("also holds 1 tip: B"));
        let found = NodeRef::mrca(["A", "E"]).find(&tree).unwrap();
        assert_eq!(found.node, tree.root());
        assert_eq!(found.also.as_deref(), Some("also holds 3 tips: B, C and D"));
        let mixed = Tree::parse_annotated_newick(
            "((A[&l=x]:1,B[&l=y]:1):1,(C[&l=x]:1,D[&l=y]:1,E[&l=y]:1,F[&l=y]:1):1);",
        )
        .unwrap();
        let found = NodeRef::holding("l", "x").find(&mixed).unwrap();
        assert_eq!(
            found.also.as_deref(),
            Some("also holds 4 tips without l x: B, D, E and 1 more")
        );
    }

    #[test]
    fn a_clade_that_is_not_there_says_why() {
        let tree = tree();
        assert_eq!(
            NodeRef::from(99).find(&tree).unwrap_err(),
            "the tree has no node 99"
        );
        assert_eq!(
            NodeRef::named("Z").find(&tree).unwrap_err(),
            "no node is named Z"
        );
        assert_eq!(
            NodeRef::mrca(["A", "Z"]).find(&tree).unwrap_err(),
            "no node is named Z"
        );
        assert_eq!(
            NodeRef::holding("lineage", "L9").find(&tree).unwrap_err(),
            "no tip has lineage L9"
        );
        assert_eq!(
            NodeRef::mrca(Vec::<String>::new()).find(&tree).unwrap_err(),
            "no tips were named"
        );
    }
}
