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
fn a_row_cap_collapses_until_the_tree_fits_and_keeps_every_tip() {
    // The point of collapsing rather than cutting the list: a pileup that
    // meets its cap counts the reads it left out, and a tree cannot, because
    // a tip is not interchangeable with the tip below it.
    let tree = balanced(64);
    let rows = |cap: Option<usize>| {
        let svg = Figure::new(Region::new("phylo", 0, 1).unwrap())
            .push(TreeTrack::new(tree.clone()).max_rows(cap))
            .to_svg();
        let mut drawn = 0usize;
        let mut held = 0usize;
        for piece in svg.split("<text").skip(1) {
            let body = piece
                .split('>')
                .nth(1)
                .unwrap_or("")
                .split('<')
                .next()
                .unwrap_or("");
            // A collapsed clade says how many tips it stands for, either as
            // "NAME (12 tips)" when the clade is named or "t0 +11 more" when
            // it is not and the first tip has to speak for it.
            let folded = body
                .rfind(" (")
                .map(|open| (open + 2, 0usize))
                .or_else(|| body.rfind(" +").map(|open| (open + 2, 1usize)));
            if let Some((start, extra)) = folded {
                if let Some(count) = body[start..]
                    .split(' ')
                    .next()
                    .and_then(|word| word.replace(',', "").parse::<usize>().ok())
                {
                    drawn += 1;
                    held += count + extra;
                    continue;
                }
            }
            if body.starts_with('t') && body[1..].chars().all(|c| c.is_ascii_digit()) {
                drawn += 1;
                held += 1;
            }
        }
        (drawn, held)
    };

    assert_eq!(
        rows(None),
        (64, 64),
        "no cap draws every tip on its own row"
    );
    for cap in [8usize, 16, 32] {
        let (drawn, held) = rows(Some(cap));
        assert_eq!(drawn, cap, "a cap of {cap} should draw {cap} rows");
        assert_eq!(held, 64, "a cap of {cap} lost tips: {held} of 64");
    }
}

#[test]
fn a_row_cap_bounds_a_height_nothing_else_could_bound() {
    let tree = balanced(512);
    let height = |cap: Option<usize>| {
        TreeTrack::new(tree.clone())
            .max_rows(cap)
            .height(&Scale::new(
                &Region::new("phylo", 0, 1).unwrap(),
                0.0,
                900.0,
            ))
    };
    // row_height floors at 2.0, so without a cap this is the least tall the
    // tree can be drawn, and it is still several screens.
    assert!(
        height(None) > 1000.0,
        "512 tips is a tall figure: {}",
        height(None)
    );
    assert!(height(Some(40)) < height(None) / 4.0);
    assert!(height(Some(40)) < height(Some(200)));
}

/// A balanced tree of `tips` leaves named t0 upwards, for the cap tests.
fn balanced(tips: usize) -> Tree {
    let mut level: Vec<String> = (0..tips).map(|i| format!("t{i}:0.1")).collect();
    while level.len() > 1 {
        level = level
            .chunks(2)
            .map(|pair| {
                if pair.len() == 2 {
                    format!("({},{}):0.1", pair[0], pair[1])
                } else {
                    pair[0].clone()
                }
            })
            .collect();
    }
    let root = level[0].rsplit_once(':').map(|(head, _)| head).unwrap();
    Tree::parse_newick(&format!("{root};")).unwrap()
}

#[test]
fn a_clade_of_one_is_one_tip() {
    // Four places in this module counted tips and all four wrote "1 tips".
    // Every assertion in this file used a clade of two, so none of them
    // ever saw it.
    let tree = Tree::parse_newick("((ONLY:0.1):0.2,(A:0.1,B:0.1,C:0.1):0.2,OUT:0.4);").unwrap();
    let of_size = |want: usize| {
        (0..tree.nodes().len())
            .find(|node| !tree.nodes()[*node].is_leaf() && tree.clade_size(*node) == want)
            .unwrap_or_else(|| panic!("no clade of {want}"))
    };
    let svg = |node: usize| {
        Figure::new(Region::new("phylo", 0, 1).unwrap())
            .push(TreeTrack::new(tree.clone()).collapse(node))
            .to_svg()
    };
    assert!(
        svg(of_size(1)).contains("(1 tip)"),
        "a clade of one is one tip"
    );
    assert!(!svg(of_size(1)).contains("1 tips"));
    assert!(svg(of_size(3)).contains("(3 tips)"));
}

#[test]
fn a_tip_name_shrinks_with_the_row_it_sits_on() {
    // The help for --row-height promises this in as many words, and both
    // the gutter measurement and the drawing have to agree on the size or
    // the gutter is held open for text that is no longer that big.
    let tree = Tree::parse_newick("((one:0.1,two:0.1):0.2,(three:0.1,four:0.1):0.2);").unwrap();
    let drawn = |row_height: f64| {
        Figure::new(Region::new("phylo", 0, 1).unwrap())
            .push(TreeTrack::new(tree.clone()).row_height(row_height))
            .to_svg()
    };
    assert!(
        drawn(2.0).contains(r#"font-size="2""#),
        "a two pixel row wants a two pixel name"
    );
    assert!(drawn(6.0).contains(r#"font-size="6""#));
    // And it stops shrinking upwards: a tall row keeps the theme's size.
    assert!(!drawn(40.0).contains(r#"font-size="40""#));
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
fn support_can_be_encoded_visibly_and_filtered_without_losing_tooltips() {
    for track in [
        TreeTrack::new(tree()),
        TreeTrack::new(tree()).circular(),
        TreeTrack::new(tree()).unrooted(),
    ] {
        let visible = Figure::new(region())
            .width(540.0)
            .show_region_label(false)
            .push(
                track
                    .clone()
                    .support_style(SupportStyle::SymbolsAndLabels)
                    .support_threshold(80.0),
            )
            .to_svg();
        assert!(visible.contains(">0.9</text>"), "{visible}");
        assert!(visible.contains("clade support 0.9"), "{visible}");

        let filtered = Figure::new(region())
            .width(540.0)
            .show_region_label(false)
            .push(
                track
                    .support_style(SupportStyle::SymbolsAndLabels)
                    .support_threshold(0.95),
            )
            .to_svg();
        assert!(!filtered.contains(">0.9</text>"), "{filtered}");
        assert!(filtered.contains("clade support 0.9"), "{filtered}");
    }
}

#[test]
fn branch_event_labels_are_direct_exact_and_projection_independent() {
    let event_tree = Tree::parse_annotated_newick(
        "((A[&event=S_D614G]:0.8,B:0.8)0.95:0.8,C[&event=N_R203K]:1.6);",
    )
    .unwrap();
    for track in [
        TreeTrack::new(event_tree.clone()),
        TreeTrack::new(event_tree.clone()).circular(),
        TreeTrack::new(event_tree.clone()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(640.0)
            .show_region_label(false)
            .push(track.branch_labels("event").branch_label_size(6.0))
            .to_svg();
        for event in ["S_D614G", "N_R203K"] {
            assert!(svg.contains(&format!(">{event}</text>")), "{svg}");
            assert!(svg.contains(&format!("event {event}")), "{svg}");
        }
        assert_eq!(svg.matches(">S_D614G</text>").count(), 1, "{svg}");
    }
}

#[test]
fn dnds_is_direct_diverging_and_projection_independent() {
    let source = concat!(
        "((A[&omega=0.2,p=0.01]:0.8,B:0.8)",
        "AB[&omega=5,p=0.03]:0.6,",
        "C[&omega=1,p=0.4]:1.4);"
    );
    for track in [
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).circular(),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(720.0)
            .show_region_label(false)
            .push(
                track
                    .show_tips(false)
                    .dnds("omega")
                    .dnds_significance("p", 0.05),
            )
            .to_svg();
        assert!(
            svg.contains("dN/dS ω 0.2 (purifying); p 0.01 (≤ 0.05)"),
            "{svg}"
        );
        assert!(
            svg.contains("dN/dS ω 5 (diversifying); p 0.03 (≤ 0.05)"),
            "{svg}"
        );
        assert!(
            svg.contains("dN/dS ω 1 (approximately neutral); p 0.4 (&gt; 0.05)"),
            "{svg}"
        );
        assert!(
            svg.contains("<title>B; dN/dS missing</title>"),
            "an ancestor's omega must not be copied onto B: {svg}"
        );
        assert!(svg.contains("stroke=\"#0072b2\""), "{svg}");
        assert!(svg.contains("stroke=\"#d55e00\""), "{svg}");
        assert!(svg.contains("stroke-width=\"2.22\""), "{svg}");
        assert!(svg.contains("stroke-dasharray=\"1.5 3\""), "{svg}");
        for label in ["purifying", "near neutral", "diversifying", "missing"] {
            assert!(svg.contains(&format!(">{label}</text>")), "{svg}");
        }
        assert!(!svg.contains("NaN"), "{svg}");
    }
}

#[test]
fn the_last_branch_colour_encoding_wins() {
    let tree = Tree::parse_newick("(A:1,B:1);").unwrap();
    let dnds = TreeTrack::new(tree.clone())
        .color_by("country")
        .dnds("omega");
    assert!(dnds.color_by.is_none());
    assert_eq!(
        dnds.dnds.as_ref().map(|layer| layer.key.as_str()),
        Some("omega")
    );

    let categorical = TreeTrack::new(tree).dnds("omega").color_by("country");
    assert_eq!(categorical.color_by.as_deref(), Some("country"));
    assert!(categorical.dnds.is_none());
}

#[test]
fn weighted_rate_classes_remain_visible_in_every_projection() {
    let source = concat!(
        "((A[&omega1=0.15,w1=0.72,omega2=4.8,w2=0.28]:0.8,B:0.8):0.6,",
        "C[&omega1=0.7,w1=2,omega2=1.4,w2=1]:1.4);"
    );
    let mixture =
        BranchRateMixture::new(["omega1", "omega2"], ["w1", "w2"]).label("aBSREL classes");
    for track in [
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).circular(),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(700.0)
            .show_region_label(false)
            .push(track.branch_rate_mixture(mixture.clone()))
            .to_svg();
        assert!(
            svg.contains(
                "aBSREL classes | class 1 omega 0.15 weight 0.72; class 2 omega 4.8 weight 0.28"
            ),
            "{svg}"
        );
        assert!(
            svg.contains("class 1 omega 0.7 weight 2; class 2 omega 1.4 weight 1"),
            "source weights must remain exact even when geometry is normalised: {svg}"
        );
        assert!(svg.contains("stroke=\"#0072b2\""), "{svg}");
        assert!(svg.contains("stroke=\"#d55e00\""), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }
}

#[test]
fn recurrent_events_connect_branches_without_inheriting_singletons() {
    let source = concat!(
        "((A[&event=S45N]:0.8,B[&event=private]:0.8):0.6,",
        "C[&event=S45N]:1.4);"
    );
    for track in [
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).circular(),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(700.0)
            .show_region_label(false)
            .push(track.homoplasy_layer(HomoplasyLayer::new("event").label("homoplasy candidates")))
            .to_svg();
        assert!(
            svg.contains("recurrent event event = S45N; 2 branches"),
            "{svg}"
        );
        assert!(!svg.contains("recurrent event event = private"), "{svg}");
        assert!(svg.contains("stroke-dasharray=\"6 4\""), "{svg}");
        assert!(
            svg.contains("homoplasy candidates; dashed curves connect"),
            "{svg}"
        );
        assert!(!svg.contains("NaN"), "{svg}");
    }
}

#[test]
fn branch_geometry_changes_connections_without_changing_the_owned_tree() {
    let source = "((A:1,B:1)AB:1,C:2)root;";
    let orthogonal = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(Tree::parse_newick(source).unwrap())
                .branch_geometry(BranchGeometry::Orthogonal),
        )
        .to_svg();
    let diagonal = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(Tree::parse_newick(source).unwrap())
                .branch_geometry(BranchGeometry::Diagonal),
        )
        .to_svg();
    let curved = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(Tree::parse_newick(source).unwrap())
                .branch_geometry(BranchGeometry::Curved),
        )
        .to_svg();
    assert!(orthogonal.matches("<line").count() > diagonal.matches("<line").count());
    assert!(curved.contains(" C "), "{curved}");
    for svg in [orthogonal, diagonal, curved] {
        for tip in ["A", "B", "C"] {
            assert!(svg.contains(&format!(">{tip}</text>")), "{svg}");
        }
        assert!(!svg.contains("NaN"), "{svg}");
    }
}

#[test]
fn ancestral_events_and_intervals_project_together_without_inheritance() {
    let source = concat!(
        "((A[&p_A=0.1,p_B=0.9,events={S45N,E88K},cf=0.82,cf_lo=0.71,cf_hi=0.91]:1,",
        "B[&p_A=0.85,p_B=0.15,events={private},cf=0.44,cf_lo=0.31,cf_hi=0.58]:1)",
        "AB[&p_A=0.9,p_B=0.1,cf=0.91,cf_lo=0.84,cf_hi=0.96]:1,",
        "C[&p_A=0.08,p_B=0.92,events={S45N},cf=0.77,cf_lo=0.62,cf_hi=0.86]:2)",
        "root[&p_A=0.95,p_B=0.05];"
    );
    let states = AncestralStateLayer::new(["p_A", "p_B"])
        .label("ancestral host")
        .confidence(0.70);
    let events = BranchEventLayer::new("events").label("ancestral mutations");
    let concordance = BranchIntervalLayer::new("cf", "cf_lo", "cf_hi")
        .label("gene concordance")
        .threshold(0.70);
    for track in [
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).circular(),
        TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(740.0)
            .show_region_label(false)
            .push(
                track
                    .ancestral_states(states.clone())
                    .branch_event_layer(events.clone())
                    .branch_interval(concordance.clone()),
            )
            .to_svg();
        assert!(svg.contains("ancestral mutations | events = S45N"), "{svg}");
        assert!(svg.contains("ancestral mutations | events = E88K"), "{svg}");
        assert!(
            svg.contains("gene concordance | estimate 0.82 | interval 0.71 to 0.91"),
            "{svg}"
        );
        assert!(
            svg.contains("ancestral host transition p_A (0.9) to p_B (0.9)"),
            "{svg}"
        );
        assert!(svg.contains("AB; p_A 0.9; p_B 0.1"), "{svg}");
        assert!(!svg.contains("NaN"), "{svg}");
    }
}

#[test]
fn negative_ancestral_probabilities_do_not_create_transition_claims() {
    let source = "(A[&p_A=-0.2,p_B=1.2]:1)root[&p_A=0.9,p_B=0.1];";
    let svg = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(Tree::parse_annotated_newick(source).unwrap())
                .ancestral_states(AncestralStateLayer::new(["p_A", "p_B"])),
        )
        .to_svg();
    assert!(!svg.contains("transition p_A (0.9) to p_B"), "{svg}");
    assert!(!svg.contains("A; p_A -0.2; p_B 1.2"), "{svg}");
    assert!(!svg.contains("NaN"), "{svg}");
}

#[test]
fn branch_length_scale_bars_are_exact_across_phylogram_projections() {
    for track in [
        TreeTrack::new(tree()),
        TreeTrack::new(tree()).circular(),
        TreeTrack::new(tree()).unrooted(),
    ] {
        let svg = Figure::new(region())
            .width(640.0)
            .show_region_label(false)
            .push(
                track
                    .scale_bar()
                    .scale_bar_length(0.1)
                    .scale_bar_unit("substitutions/site"),
            )
            .to_svg();
        assert!(
            svg.contains("<title>branch length scale 0.1 substitutions/site</title>"),
            "{svg}"
        );
        assert!(svg.contains(">0.1 substitutions/site</text>"), "{svg}");
    }

    let cladogram = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(tree())
                .shape(TreeShape::Cladogram)
                .scale_bar(),
        )
        .to_svg();
    assert!(!cladogram.contains("branch length scale"), "{cladogram}");
}

#[test]
fn track_builders_reroot_by_node_name_outgroup_and_midpoint() {
    let named_tree = Tree::parse_newick("((A:1,B:1)AB:2,(C:1,D:1)CD:2);").unwrap();
    let ab = named_tree.node_named("AB").unwrap();
    let by_node = TreeTrack::new(named_tree.clone()).reroot(ab);
    assert_eq!(by_node.tree().root(), ab);
    assert!(by_node.show_root);

    let by_name = TreeTrack::new(named_tree).reroot_named("CD");
    assert_eq!(
        by_name.tree().nodes()[by_name.tree().root()]
            .name
            .as_deref(),
        Some("CD")
    );

    let outgroup_tree = Tree::parse_newick("(((A:1,B:1)AB:2,C:3)ING:4,(O1:2,O2:2)OUT:5);").unwrap();
    let by_outgroup = TreeTrack::new(outgroup_tree).reroot_outgroup(["O1", "O2"]);
    assert!(by_outgroup.show_root);
    assert_eq!(by_outgroup.tree().leaf_count(), 5);

    let midpoint_tree = Tree::parse_newick("((A:1,B:1)AB:1,C:4);").unwrap();
    let old_nodes = midpoint_tree.nodes().len();
    let by_midpoint = TreeTrack::new(midpoint_tree).reroot_midpoint();
    assert_eq!(by_midpoint.tree().root(), old_nodes);
    assert!(by_midpoint.show_root);
}

#[test]
fn selected_root_markers_are_explicit_and_only_belong_to_rooted_projections() {
    let tree = Tree::parse_newick("((A:1,B:1)AB:2,(C:1,D:1)CD:2);").unwrap();
    let rooted = TreeTrack::new(tree.clone()).reroot_named("AB");
    for track in [rooted.clone(), rooted.clone().circular()] {
        let svg = Figure::new(region())
            .width(560.0)
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("<title>selected root</title>"), "{svg}");
    }

    let unrooted = Figure::new(region())
        .width(560.0)
        .show_region_label(false)
        .push(rooted.unrooted())
        .to_svg();
    assert!(!unrooted.contains("selected root"), "{unrooted}");

    let unchanged = TreeTrack::new(tree).reroot_named("missing");
    assert!(!unchanged.show_root);
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
    let tree =
        Tree::parse_annotated_newick("((A[&date=2024]:1,B[&date=2025]:2)AB:1,C[&date=2023]:3);")
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
        3.0 * 15.0 + 22.0
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
    let bar = TraitColumn::bar("coverage");
    let binary = TraitColumn::binary("resistant");
    let symbol = TraitColumn::symbol("host");
    assert_eq!(categorical.key(), "lineage");
    assert_eq!(categorical.scale(), TraitScale::Categorical);
    assert_eq!(categorical.trait_style(), TraitStyle::Strip);
    assert_eq!(continuous.key(), "clock_rate");
    assert_eq!(continuous.scale(), TraitScale::Continuous);
    assert_eq!(bar.trait_style(), TraitStyle::Bar);
    assert_eq!(binary.trait_style(), TraitStyle::Binary);
    assert_eq!(symbol.trait_style(), TraitStyle::Symbol);
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
    let tree =
        Tree::parse_annotated_newick("((A[&date=2024]:1,B[&date=2025]:2)AB:1,C[&date=2023]:3);")
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
fn itol_style_bars_binary_marks_and_symbols_keep_exact_values() {
    let tree = Tree::parse_annotated_newick(
            "(A[&coverage=18,resistant=true,host=human]:1,B[&coverage=30,resistant=false,host=animal]:1,C[&coverage=42,resistant=true,host=water]:1);",
        )
        .unwrap();
    let svg = Figure::new(region())
        .width(620.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .circular()
                .trait_column(TraitColumn::bar("coverage").label("Depth"))
                .trait_column(TraitColumn::binary("resistant").label("AMR"))
                .trait_column(TraitColumn::symbol("host").label("Host")),
        )
        .to_svg();
    for title in ["A; coverage 18", "B; resistant false", "C; host water"] {
        assert!(svg.contains(&format!("<title>{title}</title>")), "{svg}");
    }
    for heading in [">Depth</text>", ">AMR</text>", ">Host</text>"] {
        assert!(svg.contains(heading), "{svg}");
    }
    assert!(
        svg.contains("fill-opacity=\"0.92\""),
        "numeric radial bars: {svg}"
    );
    assert!(
        svg.contains("<polygon"),
        "shape must reinforce colour: {svg}"
    );
    assert!(!svg.contains("NaN"));
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
    // The triangle names the clade it stands for, ends included, so a reader
    // sees which part of the tree it is and not only how big.
    assert!(
        svg.contains("<title>outbreak (2 tips), A to B</title>"),
        "{svg}"
    );
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
fn unrooted_height_is_explicit_and_the_source_root_is_not_the_centre() {
    let scale = Scale::new(&region(), 0.0, 100.0);
    let tree =
        Tree::parse_annotated_newick("[&U] (((((A:1,B:1):1,C:1):1,D:1):1,E:1):1,F:1);").unwrap();
    let scene = UnrootedScene::new(&tree, TreeShape::Phylogram, &BTreeSet::new(), -90.0);
    let centre = scene
        .visible
        .iter()
        .copied()
        .find(|node| scene.parents[*node].is_none())
        .unwrap();
    assert_ne!(
        centre,
        tree.root(),
        "the Newick root must not anchor the view"
    );
    assert_eq!(
        TreeTrack::new(tree)
            .unrooted()
            .unrooted_size(360.0)
            .height(&scale),
        360.0
    );
}

#[test]
fn an_unrooted_tree_keeps_branches_labels_support_and_annotation_rings() {
    let tree = Tree::parse_annotated_newick(
            "[&U] ((A[&country=Peru]:0.1,B[&country=Chile]:0.2)0.9:0.3,(C[&country=Peru]:0.15,D[&country=Chile]:0.05):0.2);",
        )
        .unwrap();
    let branches = tree.nodes().len() - 1;
    let svg = Figure::new(region())
        .width(560.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .unrooted()
                .show_nodes(true)
                .color_by("country")
                .trait_categorical("country"),
        )
        .to_svg();
    assert!(!svg.contains("NaN"), "{svg}");
    assert!(svg.matches("<line").count() >= branches + 4, "{svg}");
    for label in [">A</text>", ">B</text>", ">C</text>", ">D</text>"] {
        assert!(svg.contains(label), "{svg}");
    }
    assert!(svg.contains("clade support 0.9"), "{svg}");
    assert!(svg.contains("<title>A; country Peru</title>"), "{svg}");
    assert!(svg.contains(">country</text>"), "{svg}");
    let branch = svg.find("<title>country Peru</title>").unwrap();
    assert!(
        svg[branch..(branch + 180).min(svg.len())].contains("stroke=\"#0072b2\""),
        "{svg}"
    );
    let ring = svg.find("<title>A; country Peru</title>").unwrap();
    assert!(
        svg[ring..(ring + 260).min(svg.len())].contains("fill=\"#0072b2\""),
        "branch and ring must share one category mapping: {svg}"
    );
}

#[test]
fn a_single_unrooted_tip_is_finite() {
    let svg = Figure::new(region())
        .show_region_label(false)
        .push(TreeTrack::new(Tree::parse_newick("A;").unwrap()).unrooted())
        .to_svg();
    assert!(!svg.contains("NaN"));
    assert!(svg.contains(">A</text>"));
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

#[test]
fn node_glyphs_keep_exact_values_and_clade_fields_keep_tip_counts() {
    let tree = Tree::parse_annotated_newick(concat!(
        "((A[&load=4,human=3,animal=1]:1,B[&load=9,human=2,animal=2]:1)",
        "outbreak[&load=16,human=6,animal=2]:1,",
        "C[&load=1,human=1,animal=0]:2)root[&load=25,human=7,animal=3];"
    ))
    .unwrap();
    let outbreak = tree.node_named("outbreak").unwrap();
    let svg = Figure::new(region())
        .width(620.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .row_height(28.0)
                .node_glyph(
                    NodeGlyph::bubble("load")
                        .label("Isolates")
                        .target(NodeGlyphTarget::Internal),
                )
                .node_glyph(
                    NodeGlyph::donut(["human", "animal"])
                        .label("Host composition")
                        .target(NodeGlyphTarget::Leaves),
                )
                .clade_highlight(
                    CladeHighlight::new(outbreak)
                        .label("Transmission cluster")
                        .opacity(0.16),
                ),
        )
        .to_svg();
    assert!(svg.contains("outbreak; load 16"), "{svg}");
    assert!(svg.contains("A; human 3; animal 1"), "{svg}");
    assert!(svg.contains("Transmission cluster; 2 tips"), "{svg}");
    assert!(svg.contains("fill-opacity=\"0.16\""), "{svg}");
    assert!(svg.contains(">Isolates</text>"), "{svg}");
    assert!(svg.contains(">human</text>"), "{svg}");
    assert!(!svg.contains("NaN"));
}

#[test]
fn every_node_glyph_and_highlight_projects_without_losing_data() {
    let source = concat!(
        "((A[&load=4,x=3,y=1]:1,B[&load=9,x=2,y=2]:1)",
        "group[&load=16,x=6,y=2]:1,C[&load=1,x=1,y=0]:2)",
        "root[&load=25,x=7,y=3];"
    );
    for projection in [
        TreeProjection::Rectangular,
        TreeProjection::Circular,
        TreeProjection::Unrooted,
    ] {
        let tree = Tree::parse_annotated_newick(source).unwrap();
        let group = tree.node_named("group").unwrap();
        let mut track = TreeTrack::new(tree)
            .projection(projection)
            .radial_size(420.0)
            .node_glyph(NodeGlyph::bubble("load").target(NodeGlyphTarget::Internal))
            .node_glyph(NodeGlyph::pie(["x", "y"]).target(NodeGlyphTarget::Leaves))
            .node_glyph(NodeGlyph::stacked_bar(["x", "y"]).target(NodeGlyphTarget::Internal))
            .clade_highlight(CladeHighlight::new(group).label("group"));
        if projection == TreeProjection::Unrooted {
            track = track.unrooted();
        }
        let svg = Figure::new(region())
            .width(600.0)
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("group; load 16"), "{projection:?}: {svg}");
        assert!(svg.contains("A; x 3; y 1"), "{projection:?}: {svg}");
        assert!(svg.contains("group; 2 tips"), "{projection:?}: {svg}");
        assert!(!svg.contains("NaN"), "{projection:?}: {svg}");
    }
}

#[test]
fn the_tips_beyond_an_edge_are_counted_the_same_whether_walked_or_looked_up() {
    // The unrooted layout used to answer "how many tips lie that way" by
    // walking that way and counting, once per edge, which made a large tree
    // quadratic: 121 ms at 500 rows but 2836 at 8000, while the other two
    // projections stayed flat. ComponentTerminals answers the same question
    // from one pass. This is the check that the two answers agree, on every
    // directed edge, for shapes that break lazy reasoning: a star has one
    // node adjacent to everything, a caterpillar is as deep as it is wide,
    // and a hidden clade can leave the visible edges in disconnected pieces,
    // where a count taken across pieces would be wrong.
    let walked = |start: usize,
                  blocked: usize,
                  adjacency: &[Vec<(usize, f64)>],
                  terminals: &BTreeSet<usize>| {
        let mut count = 0usize;
        let mut stack = vec![(start, blocked)];
        while let Some((node, parent)) = stack.pop() {
            count += usize::from(terminals.contains(&node));
            for (next, _) in &adjacency[node] {
                if *next != parent {
                    stack.push((*next, node));
                }
            }
        }
        count
    };

    let shapes = [
        (
            "star",
            Tree::parse_newick("(a:0.1,b:0.2,c:0.3,d:0.4,e:0.5,f:0.6);").unwrap(),
        ),
        (
            "caterpillar",
            Tree::parse_newick(&format!("{}t8:0.1{};", "(t0:0.1,".repeat(8), ")".repeat(8)))
                .unwrap(),
        ),
        ("one tip", Tree::parse_newick("(only:0.5);").unwrap()),
        ("balanced", balanced(64)),
    ];

    let mut edges = 0usize;
    for (name, tree) in shapes {
        // Nothing collapsed, then every other internal node, then all of
        // them: collapsing is what hides nodes, and hidden nodes are what can
        // leave the visible edges in disconnected pieces.
        let internal: Vec<usize> = (0..tree.nodes().len())
            .filter(|node| !tree.nodes()[*node].is_leaf())
            .collect();
        let folds: [BTreeSet<usize>; 3] = [
            BTreeSet::new(),
            internal.iter().copied().step_by(2).collect(),
            internal.iter().copied().collect(),
        ];
        for (fold, collapsed) in folds.into_iter().enumerate() {
            let cap = fold;
            let visibility = visible_nodes(&tree, &collapsed);
            let terminals: BTreeSet<usize> =
                visible_terminals(&tree, &collapsed).into_iter().collect();
            let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); tree.nodes().len()];
            for (node, clade) in tree.nodes().iter().enumerate() {
                let Some(parent) = clade.parent else { continue };
                if !visibility[node] || !visibility[parent] {
                    continue;
                }
                adjacency[parent].push((node, 1.0));
                adjacency[node].push((parent, 1.0));
            }

            let counts = ComponentTerminals::new(&adjacency, &terminals);
            for node in 0..adjacency.len() {
                for (next, _) in &adjacency[node] {
                    let want = walked(*next, node, &adjacency, &terminals);
                    assert_eq!(
                        counts.beyond(*next, node),
                        want,
                        "{name} at cap {cap:?}: beyond({next}, {node})"
                    );
                    edges += 1;
                }
            }
        }
    }
    assert!(
        edges > 200,
        "only {edges} edges checked, too few to mean much"
    );
}

#[test]
fn an_unrooted_tree_is_centred_on_its_own_drawing_and_stays_inside_the_band() {
    // The layout's origin is whichever node the walk started from, and a tree
    // hangs off that node however its branches happen to fall. Fitting the
    // drawing by its distance from that origin left it off centre and small:
    // measured against the binary built before the fix, a two hundred tip tree
    // sat 47.8 px off the middle of its band and its branches spanned 454 px
    // of an 866 px band, and a tree with one long branch sat 82.3 px off and
    // spanned 74 px of 566. The same two are now 0.0 px off and 692 px, and
    // 5.8 px off and 174 px. Neither ran outside the clip; an earlier note
    // here said they did, which came from measuring the scale bar, a rule
    // drawn the full width of the figure, as though it were a branch.
    // One long branch, so the drawing is lopsided while the walk still starts
    // from the middle of the chain: the starting node is picked by how many
    // tips lie each way, which a long branch does not change, so the origin
    // and the middle of the picture are in different places.
    let caterpillar = Tree::parse_newick(&format!(
        "{}t24:40.0{};",
        "(t0:0.4,".repeat(24),
        ")".repeat(24)
    ))
    .unwrap();
    let svg = Figure::new(region())
        .width(600.0)
        .show_region_label(false)
        .push(TreeTrack::new(caterpillar).unrooted())
        .to_svg();

    let numbers = |tag: &str| {
        svg.split(tag)
            .skip(1)
            .filter_map(|piece| piece.split('"').next()?.parse::<f64>().ok())
            .collect::<Vec<_>>()
    };
    let xs: Vec<f64> = numbers("x1=\"")
        .into_iter()
        .chain(numbers("x2=\""))
        .collect();
    let ys: Vec<f64> = numbers("y1=\"")
        .into_iter()
        .chain(numbers("y2=\""))
        .collect();
    assert!(
        xs.len() > 40,
        "too few branches to mean anything: {}",
        xs.len()
    );

    let rect = svg
        .split("<clipPath")
        .nth(1)
        .and_then(|piece| piece.split("<rect ").nth(1))
        .expect("the track is clipped");
    let of = |name: &str| {
        rect.split(&format!("{name}=\""))
            .nth(1)
            .and_then(|piece| piece.split('"').next())
            .and_then(|value| value.parse::<f64>().ok())
            .expect("the clip carries this")
    };
    let (left, top) = (of("x"), of("y"));
    let (right, bottom) = (left + of("width"), top + of("height"));

    let (x0, x1) = (
        xs.iter().copied().fold(f64::MAX, f64::min),
        xs.iter().copied().fold(f64::MIN, f64::max),
    );
    let (y0, y1) = (
        ys.iter().copied().fold(f64::MAX, f64::min),
        ys.iter().copied().fold(f64::MIN, f64::max),
    );
    assert!(
        x0 >= left - 0.5 && x1 <= right + 0.5 && y0 >= top - 0.5 && y1 <= bottom + 0.5,
        "the drawing runs outside its band: x {x0}..{x1} and y {y0}..{y1} in {left}..{right} by {top}..{bottom}"
    );
    // Centred within a pixel, not merely inside.
    let slack = 1.0;
    assert!(
        ((x0 + x1) / 2.0 - (left + right) / 2.0).abs() < slack
            && ((y0 + y1) / 2.0 - (top + bottom) / 2.0).abs() < slack,
        "the drawing is off centre: middle ({}, {}) against band middle ({}, {})",
        (x0 + x1) / 2.0,
        (y0 + y1) / 2.0,
        (left + right) / 2.0,
        (top + bottom) / 2.0
    );
}

#[test]
fn an_unrooted_tree_moves_its_names_to_a_ring_only_once_they_would_collide() {
    // Measured on forty five trees of eight to forty tips: with the names at
    // the tips, every figure that really had two names touching had its
    // closest pair under eight tenths of a body, so that is where the drawing
    // gives up on tip names and gathers them onto a circle with a leader each.
    let leaders = |tips: usize| {
        let svg = Figure::new(region())
            .width(700.0)
            .show_region_label(false)
            .push(TreeTrack::new(balanced(tips)).unrooted())
            .to_svg();
        svg.matches("stroke-width=\"0.8\"").count()
    };
    assert_eq!(leaders(6), 0, "six names fit at their own tips");
    assert!(
        leaders(128) >= 128,
        "a hundred and twenty eight names cannot fit at their tips and must be gathered"
    );
}

#[test]
fn folding_ranks_clades_by_the_same_size_the_tree_would_report() {
    // Folding takes the smallest clade first, and it used to ask the tree for
    // each clade's size from inside the sort. That call collects the whole
    // subtree into a vector to count the leaves in it, so ranking n clades is
    // n traversals and a sort asks for a key more than once each: a hundred
    // thousand tip tree took 785 ms to reach sixty rows, against 69 ms to draw
    // all hundred thousand with no cap at all. The sizes now come off the walk
    // that is already being made, and these are the folds that ranking picks.
    // On a balanced sixty four they are exact: four rows is the four sixteens,
    // and nothing else adds up to it smallest first.
    let folds = |cap: usize| {
        let svg = Figure::new(region())
            .push(TreeTrack::new(balanced(64)).max_rows(Some(cap)))
            .to_svg();
        let mut folded: Vec<usize> = svg
            .split("<text")
            .skip(1)
            .filter_map(|piece| {
                let body = piece.split('>').nth(1)?.split('<').next()?;
                let at = body.rfind(" +")?;
                body[at + 2..]
                    .split(' ')
                    .next()?
                    .parse::<usize>()
                    .ok()
                    .map(|more| more + 1)
            })
            .collect();
        folded.sort_unstable();
        folded
    };
    assert_eq!(folds(4), vec![16, 16, 16, 16]);
    assert_eq!(folds(9), vec![4, 4, 8, 8, 8, 8, 8, 8, 8]);
    assert_eq!(folds(31), [vec![2; 30], vec![4]].concat());
    // Whatever it folds, no tip is lost and the cap is met exactly.
    for cap in [4usize, 9, 31] {
        let folded = folds(cap);
        let saved: usize = folded.iter().map(|size| size - 1).sum();
        assert_eq!(64 - saved, cap, "a cap of {cap} did not land on {cap} rows");
    }
}

#[test]
fn a_folded_triangle_names_the_clade_it_stands_for() {
    // A tooltip that said "clade (128 tips)" told a reader how big the
    // triangle is and nothing about which part of the tree it is. The two tips
    // at its ends name it instead: the tips under a node are a run, so the
    // first and the last of them pick out one clade and no other, which is
    // also an address anything driving the program can hand back to open it.
    let tree = balanced(64);
    let svg = Figure::new(region())
        .push(TreeTrack::new(tree.clone()).max_rows(Some(4)))
        .to_svg();
    let titles: Vec<&str> = svg
        .split("<title>")
        .skip(1)
        .filter_map(|piece| piece.split("</title>").next())
        .filter(|title| title.contains(" to "))
        .collect();
    assert_eq!(titles.len(), 4, "four folds, four addresses: {svg}");
    assert_eq!(
        titles[0], "clade (16 tips), t0 to t15",
        "the first triangle holds the first sixteen tips"
    );
    assert_eq!(titles[3], "clade (16 tips), t48 to t63");

    // And a one tip fold has nothing to span, so it says its size alone.
    let pair = Tree::parse_newick("((a:0.1,b:0.1):0.1,c:0.1);").unwrap();
    let one = Figure::new(region())
        .push(TreeTrack::new(pair).max_rows(Some(2)))
        .to_svg();
    assert!(
        one.contains("<title>clade (2 tips), a to b</title>"),
        "{one}"
    );
}

#[test]
fn a_triangle_prints_an_address_the_focus_flag_takes() {
    // The viewer used to read this off the page and hand it back, and does not
    // any more: it flies over coordinates the crate works out, so there are no
    // triangles in it to click. The address is still what a reader gets from a
    // figure and types at a shell, so the two halves are still checked against
    // each other here.
    let svg = Figure::new(region())
        .push(TreeTrack::new(balanced(32)).max_rows(Some(4)))
        .to_svg();
    let title = svg
        .split("<title>")
        .skip(1)
        .filter_map(|piece| piece.split("</title>").next())
        .find(|title| title.contains(" to "))
        .expect("a folded tree has a triangle");

    let (label, rest) = title.split_once(" (").expect("a bracket");
    let (held, span) = rest.split_once("), ").expect("a close and a comma");
    let (first, last) = span.split_once(" to ").expect("two names");
    assert_eq!(label, "clade");
    assert_eq!(held, "8 tips");
    assert_eq!((first, last), ("t0", "t7"));

    let line = format!("tree:1-1 --tree t.nwk --focus {first},{last}");
    let words: Vec<String> = line.split_whitespace().map(String::from).collect();
    match crate::cli::args::parse(&words).expect("the address parses as a flag") {
        crate::cli::args::Request::Draw(invocation) => assert_eq!(
            invocation.tracks[0].focus,
            Some(vec![first.to_string(), last.to_string()])
        ),
        other => panic!("expected a figure, got {other:?}"),
    }
}

#[test]
fn a_round_tree_is_as_tall_as_the_band_it_is_given_is_wide() {
    // The disc is drawn into the shorter of the band's two sides, so the height
    // reserved for it has to come from the same number the drawing gets. It
    // used to come from the figure's inner right edge instead, which never
    // moves: the height stayed put while the band narrowed under it, and the
    // disc shrank inside a reservation it no longer filled. Two hundred tips at
    // 900 px kept a height of 905.515 px while the disc went from radius 339.40
    // to 283.40 with nothing changed but the length of this track's own gutter
    // label, closing the gap between names from 10.66 px to 8.90 against an
    // 11 px body.
    let band = |label: &str| {
        let svg = Figure::new(region())
            .width(900.0)
            .show_region_label(false)
            .push(TreeTrack::new(balanced(400)).circular().label(label))
            .to_svg();
        let rect = svg
            .split("<clipPath")
            .nth(1)
            .and_then(|piece| piece.split("<rect ").nth(1))
            .expect("the track is clipped");
        let of = |name: &str| {
            rect.split(&format!("{name}=\""))
                .nth(1)
                .and_then(|piece| piece.split('"').next())
                .and_then(|value| value.parse::<f64>().ok())
                .expect("the clip carries this")
        };
        (of("width"), of("height"))
    };

    // Square, because the disc wants more room than the band has and takes all
    // of it in both directions.
    for label in ["t", "a rather longer name for this track"] {
        let (wide, tall) = band(label);
        assert!(
            (wide - tall).abs() < 0.01,
            "{label:?}: the band is {wide} by {tall}, so the disc does not fill what was reserved"
        );
    }

    // And the band really does move with the label, or the check above would
    // be measuring one case twice.
    let (short, _) = band("t");
    let (long, _) = band("a rather longer name for this track");
    assert!(
        long < short - 50.0,
        "a longer gutter label should leave a narrower band: {short} then {long}"
    );
}

#[test]
fn a_folded_clade_is_drawn_in_the_colour_its_value_earns() {
    // A cell with a value and no colour to draw it in is an empty outline, and
    // an empty outline is how this figure says it knows nothing. The colours
    // come from a domain built by walking the drawn nodes, and a folded row's
    // value belongs to none of them: it is what its tips agree on. So the
    // domain is fed the row values too, or forty rows of lineage come out as
    // forty empty outlines with forty tooltips that each name a lineage.
    let tree = Tree::parse_annotated_newick(concat!(
        "((a[&group=one]:0.1,b[&group=one]:0.1):0.1,",
        "(c[&group=two]:0.1,d[&group=two]:0.1):0.1);"
    ))
    .unwrap();
    let svg = Figure::new(region())
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .max_rows(Some(2))
                .trait_categorical("group"),
        )
        .to_svg();

    // Both rows are folded, both clades agree, and the two values are
    // different, so two cells are filled and they are not the same colour.
    let filled: Vec<&str> = svg
        .split("<rect")
        .skip(1)
        .filter_map(|piece| piece.split("fill=\"").nth(1)?.split('"').next())
        .filter(|fill| fill.starts_with('#') && *fill != "#ffffff")
        .collect();
    assert_eq!(filled.len(), 2, "one filled cell per folded row: {svg}");
    assert_ne!(
        filled[0], filled[1],
        "two clades that agree on different values are two colours: {svg}"
    );
    assert!(svg.contains("a +1 more; group one"), "{svg}");
    assert!(svg.contains("c +1 more; group two"), "{svg}");
}

#[test]
fn a_clade_whose_tips_agree_is_coloured_as_that_clade() {
    // Colouring by an annotation keyed on sample names reaches the tips and
    // stops there, because a branch above them carries no sample's value: a
    // two hundred tip tree coloured by lineage had its two hundred terminal
    // branches in six colours and all three hundred and ninety seven internal
    // ones left plain, which says the clades are of unknown lineage when every
    // tip in them says otherwise. A branch now takes the value its descendants
    // agree on, and a branch above two lineages stays plain, which is true.
    let tree = Tree::parse_annotated_newick(concat!(
        "((a[&group=one]:0.1,b[&group=one]:0.1):0.1,",
        "(c[&group=two]:0.1,d[&group=two]:0.1):0.1);"
    ))
    .unwrap();
    let svg = Figure::new(region())
        .show_region_label(false)
        .push(TreeTrack::new(tree).color_by("group"))
        .to_svg();

    let strokes: Vec<&str> = svg
        .split("<line")
        .skip(1)
        .filter_map(|piece| piece.split("stroke=\"").nth(1)?.split('"').next())
        .filter(|stroke| stroke.starts_with('#'))
        .collect();
    let coloured: Vec<&&str> = strokes
        .iter()
        .filter(|stroke| **stroke != "#1b1f23" && **stroke != "#d7dce2")
        .collect();
    assert!(
        coloured.len() >= 6,
        "each pair and the branch holding it takes a colour: {strokes:?}"
    );
    let hues: std::collections::BTreeSet<&&str> = coloured.into_iter().collect();
    assert_eq!(hues.len(), 2, "two groups, two colours: {svg}");
    // And the root, which holds both groups, is not either of them.
    assert!(
        strokes.contains(&"#1b1f23"),
        "the branch above two groups stays plain: {svg}"
    );
}
