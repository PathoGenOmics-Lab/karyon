//! Changes carried on a branch, and an index that can be asked who carries one.
//!
//! A phylogeny of an outbreak is read as much by its mutations as by its shape:
//! the question is rarely "what does this clade look like" and often "who
//! carries S:D614G". Answering that needs the changes to be structure rather
//! than a string, and needs an index, because walking a million branches to
//! find one spelling is a walk per question.
//!
//! # Where they come from
//!
//! An annotated Newick carries them, which this crate already reads:
//!
//! ```text
//! (a[&mutations="A123T,S:D614G"]:0.1,b:0.2)[&mutations="C241T"]:0.3;
//! ```
//!
//! The key is named rather than fixed, because the tools that write these files
//! do not agree on one: `mutations`, `muts` and `aa_muts` are all in the wild,
//! and a file may carry two of them for nucleotide and amino acid changes.
//!
//! # What is a branch's and what is a node's
//!
//! A mutation belongs to the branch leading to a node, and everything below
//! that node inherits it. That is why [`Mutations::carriers`] answers with a
//! subtree and not with one branch: the change happened once and every
//! descendant has it.

use std::collections::BTreeMap;

use super::Tree;

/// One change on a branch: where, from what, to what.
///
/// A bare `A123T` is a nucleotide change and has no gene. `S:D614G` is an amino
/// acid change in a named gene. Both spellings are in use, sometimes in the
/// same file, and both round-trip through [`Mutation::to_string`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Mutation {
    /// The gene or segment, where the spelling names one.
    pub gene: Option<String>,
    /// The position, in whatever coordinates the file is in.
    pub position: u64,
    /// What was there.
    pub from: String,
    /// What replaced it.
    pub to: String,
}

impl Mutation {
    /// Reads one mutation, or nothing if this is not one.
    ///
    /// Accepts `A123T`, `S:D614G`, and either with a `nt:` or `aa:` prefix,
    /// which some writers add and which says nothing the rest of the spelling
    /// does not. A gene name may itself contain a colon, so the gene is taken
    /// as everything before the *last* one.
    ///
    /// What it will not do is guess. A spelling it cannot take apart is
    /// returned as nothing rather than as a mutation at position zero, because
    /// a figure drawn from a misread file looks exactly like a figure drawn
    /// from a good one.
    pub fn parse(text: &str) -> Option<Mutation> {
        let text = text.trim();
        if text.is_empty() {
            return None;
        }
        // `nt:` and `aa:` say what kind of change it is, which the rest of the
        // spelling already says, so they come off before anything else is read.
        // They come off first and not as part of the gene rule, or `aa:S:D614G`
        // reads as a gene called `aa:S`.
        let text = ["nt:", "aa:", "NT:", "AA:"]
            .iter()
            .find_map(|prefix| text.strip_prefix(*prefix))
            .unwrap_or(text);
        let (gene, change) = match text.rfind(':') {
            Some(at) => (Some(text[..at].to_string()), &text[at + 1..]),
            None => (None, text),
        };
        if gene.as_deref() == Some("") {
            return None;
        }

        let bytes = change.as_bytes();
        let mut at = 0;
        while at < bytes.len() && bytes[at].is_ascii_alphabetic() {
            at += 1;
        }
        let from = &change[..at];
        let start = at;
        while at < bytes.len() && bytes[at].is_ascii_digit() {
            at += 1;
        }
        let digits = &change[start..at];
        let to = &change[at..];

        if from.is_empty() || digits.is_empty() || to.is_empty() {
            return None;
        }
        if !to
            .bytes()
            .all(|b| b.is_ascii_alphabetic() || b == b'*' || b == b'-')
        {
            return None;
        }
        Some(Mutation {
            gene,
            position: digits.parse().ok()?,
            from: from.to_ascii_uppercase(),
            to: to.to_ascii_uppercase(),
        })
    }

    /// Reads a list of them, however the file separated it.
    ///
    /// Commas, semicolons, pipes and spaces all appear as separators, so all
    /// four are taken. A piece that is not a mutation is skipped rather than
    /// failing the list, since a writer that adds a note beside the changes
    /// should not cost a reader the changes.
    pub fn parse_list(text: &str) -> Vec<Mutation> {
        text.split(|c: char| c == ',' || c == ';' || c == '|' || c.is_whitespace())
            .filter_map(Mutation::parse)
            .collect()
    }
}

impl std::fmt::Display for Mutation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.gene {
            Some(gene) => write!(f, "{gene}:{}{}{}", self.from, self.position, self.to),
            None => write!(f, "{}{}{}", self.from, self.position, self.to),
        }
    }
}

/// Every mutation a tree carries, and who carries it.
#[derive(Debug, Clone, Default)]
pub struct Mutations {
    /// What happened on the branch leading to each node.
    on_branch: Vec<Vec<Mutation>>,
    /// Where each spelling occurs, as the branches it happened on. A spelling
    /// with more than one branch is a change that happened more than once.
    branches: BTreeMap<String, Vec<usize>>,
}

impl Mutations {
    /// Reads the mutations a tree carries under `key`.
    ///
    /// Direct branch data and never inherited: a node's own annotation is the
    /// change that happened on the branch above it, and the inheritance is what
    /// [`Mutations::carriers`] works out from the shape of the tree.
    pub fn read(tree: &Tree, key: &str) -> Mutations {
        let mut on_branch = vec![Vec::new(); tree.nodes().len()];
        let mut branches: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (node, held) in on_branch.iter_mut().enumerate() {
            let Some(value) = tree.annotation(node, key) else {
                continue;
            };
            let found = Mutation::parse_list(&value.to_string());
            for mutation in &found {
                branches.entry(mutation.to_string()).or_default().push(node);
            }
            *held = found;
        }
        Mutations {
            on_branch,
            branches,
        }
    }

    /// What happened on the branch leading to `node`.
    pub fn on(&self, node: usize) -> &[Mutation] {
        self.on_branch.get(node).map_or(&[], Vec::as_slice)
    }

    /// Whether anything was read at all.
    pub fn is_empty(&self) -> bool {
        self.branches.is_empty()
    }

    /// How many distinct changes the tree carries.
    pub fn distinct(&self) -> usize {
        self.branches.len()
    }

    /// Every spelling the tree carries, in order.
    pub fn spellings(&self) -> impl Iterator<Item = &str> {
        self.branches.keys().map(String::as_str)
    }

    /// The branches a change happened on, which is more than one where it
    /// happened more than once.
    pub fn happened_on(&self, spelling: &str) -> &[usize] {
        self.branches.get(spelling).map_or(&[], Vec::as_slice)
    }

    /// Every node at or below a branch where this change happened.
    ///
    /// The change happened once and everything under it has it, so this is a
    /// subtree and not a branch. A change that happened on two branches gives
    /// both subtrees.
    pub fn carriers(&self, tree: &Tree, spelling: &str) -> Vec<usize> {
        let mut carrying = vec![false; tree.nodes().len()];
        let mut stack: Vec<usize> = self.happened_on(spelling).to_vec();
        while let Some(node) = stack.pop() {
            if carrying[node] {
                continue;
            }
            carrying[node] = true;
            stack.extend(tree.nodes()[node].children.iter().copied());
        }
        (0..carrying.len()).filter(|node| carrying[*node]).collect()
    }

    /// Everything that happened on the way from the root to `node`, in order.
    ///
    /// Which is this node's genotype, as far as the file records it.
    pub fn path(&self, tree: &Tree, node: usize) -> Vec<&Mutation> {
        let mut chain: Vec<usize> = tree.ancestors(node);
        chain.reverse();
        chain.push(node);
        chain
            .into_iter()
            .flat_map(|at| self.on(at).iter())
            .collect()
    }
}
