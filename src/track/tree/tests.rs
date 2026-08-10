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
    assert!(svg.contains("<title>outbreak (2 tips)</title>"), "{svg}");
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
