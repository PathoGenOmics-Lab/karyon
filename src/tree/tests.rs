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
    let lengths: Vec<Option<f64>> = tree.nodes().iter().map(|node| node.branch_length).collect();
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
        Tree::parse_newick("[&R] ((A[&rate=1.2]:0.1,B:0.2)0.98:0.3[&height=0.4],C:0.4);").unwrap();
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
    let tree = Tree::parse_annotated_newick("((A[&country=Peru]:1,B[&country=Chile]:1)AB:2,C:3);")
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

#[test]
fn a_mutation_is_read_however_the_file_spells_it() {
    // The spellings in the wild, and what each one means. A gene name may hold
    // a colon of its own, so the gene is everything before the last one.
    for (text, gene, position, from, to) in [
        ("A123T", None, 123u64, "A", "T"),
        ("S:D614G", Some("S"), 614, "D", "G"),
        ("nt:C241T", None, 241, "C", "T"),
        ("aa:S:D614G", Some("S"), 614, "D", "G"),
        ("ORF1ab:nsp3:P822L", Some("ORF1ab:nsp3"), 822, "P", "L"),
        ("s:d614g", Some("s"), 614, "D", "G"),
        ("N:R203*", Some("N"), 203, "R", "*"),
    ] {
        let read = Mutation::parse(text).unwrap_or_else(|| panic!("{text} should read"));
        assert_eq!(read.gene.as_deref(), gene, "{text}");
        assert_eq!(read.position, position, "{text}");
        assert_eq!(read.from, from, "{text}");
        assert_eq!(read.to, to, "{text}");
    }

    // And what is not a mutation is not guessed at, because a figure drawn from
    // a misread file looks exactly like one drawn from a good file.
    for text in [
        "", "hello", "123", "A123", "123T", "A12.5T", "S:", ":D614G", "A123-4T",
    ] {
        assert!(
            Mutation::parse(text).is_none(),
            "{text:?} should not read as a mutation"
        );
    }

    // A list, however the file separated it, and a note among them costs
    // nothing.
    let list = Mutation::parse_list("A123T, S:D614G;C241T|  N:R203K   and some words");
    assert_eq!(
        list.iter().map(ToString::to_string).collect::<Vec<_>>(),
        ["A123T", "S:D614G", "C241T", "N:R203K"]
    );
}

#[test]
fn a_change_is_carried_by_everything_below_where_it_happened() {
    // The whole point of the index. A mutation is on a branch, and every tip
    // under that branch has it: asking who carries one is a question about the
    // shape of the tree, not a search through the tips.
    let tree = Tree::parse_annotated_newick(concat!(
        "((a[&muts=\"A1T\"]:0.1,b:0.1)[&muts=\"S:D614G\"]:0.1,",
        "(c[&muts=\"S:D614G\"]:0.1,d:0.1):0.1);"
    ))
    .unwrap();
    let found = Mutations::read(&tree, "muts");

    assert_eq!(
        found.distinct(),
        2,
        "two spellings: {:?}",
        found.spellings().collect::<Vec<_>>()
    );

    // S:D614G happened twice, on two separate branches, which is what makes it
    // worth asking about at all.
    assert_eq!(found.happened_on("S:D614G").len(), 2);

    let named = |nodes: Vec<usize>| {
        let mut names: Vec<String> = nodes
            .iter()
            .filter_map(|node| tree.nodes()[*node].name.clone())
            .collect();
        names.sort();
        names
    };
    // a and b are under the first branch it happened on; c has it directly.
    assert_eq!(named(found.carriers(&tree, "S:D614G")), ["a", "b", "c"]);
    assert_eq!(named(found.carriers(&tree, "A1T")), ["a"]);
    assert!(
        found.carriers(&tree, "Z9Z").is_empty(),
        "a change nothing carries"
    );

    // And a tip's own genotype is what happened on the way down to it.
    let a = tree.node_named("a").unwrap();
    let path: Vec<String> = found
        .path(&tree, a)
        .iter()
        .map(ToString::to_string)
        .collect();
    assert_eq!(path, ["S:D614G", "A1T"], "root to tip, in order");
}

/// Laying a tree out without a root is a different walk from laying it out
/// with one, and the difference is the whole point: what a branch hangs from
/// changes, because the tree is re-rooted at its middle.
#[test]
fn an_unrooted_layout_hangs_everything_off_the_middle() {
    // A ladder, so the root sits at one end and the middle is somewhere else.
    // A tree whose root already is its middle would let this pass without
    // checking anything.
    let tree = Tree::parse_newick("((((((a:1,b:1):1,c:1):1,d:1):1,e:1):1,f:1):1,g:1);").unwrap();
    let laid = tree.unrooted(false);
    assert_ne!(
        laid.centre,
        tree.root(),
        "the middle is the root here, so this proves nothing"
    );

    assert_eq!(
        laid.spots.len(),
        tree.nodes().len(),
        "a tree of one piece leaves nothing unplaced"
    );

    let mut toward = vec![None; tree.nodes().len()];
    for spot in &laid.spots {
        toward[spot.node] = Some(spot.toward);
    }
    let loose: Vec<usize> = (0..tree.nodes().len())
        .filter(|node| toward[*node] == Some(None))
        .collect();
    assert_eq!(
        loose,
        vec![laid.centre],
        "one loose end, and it is the middle"
    );

    // Follow what every node hangs from and it ends at the middle, which is
    // what makes this a tree and not a heap of branches.
    for node in 0..tree.nodes().len() {
        let mut walk = node;
        let mut steps = 0;
        while let Some(Some(next)) = toward[walk] {
            walk = next;
            steps += 1;
            assert!(steps <= tree.nodes().len(), "walking round in a circle");
        }
        assert_eq!(walk, laid.centre, "node {node} does not reach the middle");
    }

    // The terminals are the leaves, and no more.
    let mut tips = laid.terminals.clone();
    tips.sort_unstable();
    let mut leaves = tree.leaves();
    leaves.sort_unstable();
    assert_eq!(tips, leaves, "the terminals are not the leaves");

    // And the middle really is toward the middle: no branch off it leads to
    // more than half the tips.
    let mut beyond = vec![0usize; tree.nodes().len()];
    for tip in &laid.terminals {
        let mut walk = *tip;
        while let Some(Some(next)) = toward[walk] {
            if next == laid.centre {
                beyond[walk] += 1;
                break;
            }
            walk = next;
        }
    }
    let most = beyond.iter().copied().max().unwrap_or(0);
    assert!(
        most * 2 <= laid.terminals.len() + 1,
        "one branch off the middle holds {most} of {} tips",
        laid.terminals.len()
    );
}

/// A cladogram counts branches where a phylogram measures them, and without a
/// root that shows up in how long each branch is drawn.
#[test]
fn an_unrooted_cladogram_counts_rather_than_measures() {
    let tree =
        Tree::parse_newick("(((a:0.01,b:5.0):0.2,c:1.0):0.5,(d:3.0,e:0.001):2.0,f:0.7);").unwrap();

    let edges = |laid: &Unrooted| {
        let mut at = vec![(0.0f64, 0.0f64); tree.nodes().len()];
        for spot in &laid.spots {
            at[spot.node] = (spot.x, spot.y);
        }
        let mut out = Vec::new();
        for spot in &laid.spots {
            if let Some(toward) = spot.toward {
                let (ax, ay) = at[spot.node];
                let (bx, by) = at[toward];
                // The branch this node hangs by, whichever end the walk came
                // from: unrooted, one of the two owns the length.
                let length = tree.nodes()[spot.node]
                    .branch_length
                    .or(tree.nodes()[toward].branch_length)
                    .unwrap_or(1.0);
                out.push(((ax - bx).hypot(ay - by), length));
            }
        }
        out
    };

    for (drawn, _) in edges(&tree.unrooted(true)) {
        assert!(
            (drawn - 1.0).abs() < 1e-9,
            "counting branches drew one of them {drawn} long"
        );
    }
    let measured = edges(&tree.unrooted(false));
    assert!(
        measured.iter().any(|(drawn, _)| (drawn - 1.0).abs() > 0.1),
        "measuring branches drew them all the same length"
    );
    for (drawn, length) in &measured {
        assert!(
            (drawn - length).abs() < 1e-9,
            "a branch of {length} was drawn {drawn} long"
        );
    }
}
