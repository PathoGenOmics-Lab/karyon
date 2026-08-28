use super::*;

fn phylo_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "((A[&date=2023,country=Peru]:1,B[&date=2024,country=Peru]:1)",
        "[&country=Peru]:1,(C[&date=2023.5,country=Spain]:1,",
        "D[&date=2024.2,country=Kenya]:1):1)[&date=2022];"
    ))
    .unwrap()
}

fn phylo_coordinates() -> [GeoLocation; 3] {
    [
        GeoLocation::new("Peru", -9.19, -75.0152),
        GeoLocation::new("Spain", 40.4637, -3.7492),
        GeoLocation::new("Kenya", -0.0236, 37.9062),
    ]
}

#[test]
fn equirectangular_places_the_origin_at_the_centre() {
    let area = MapRect {
        x: 10.0,
        y: 20.0,
        w: 360.0,
        h: 180.0,
    };
    assert_eq!(
        project(
            GeoPosition::new(0.0, 0.0),
            GeoProjection::Equirectangular,
            area
        ),
        Some((190.0, 110.0))
    );
}

#[test]
fn orthographic_projection_hides_the_back_hemisphere() {
    let area = MapRect {
        x: 0.0,
        y: 0.0,
        w: 200.0,
        h: 200.0,
    };
    let projection = GeoProjection::orthographic(0.0, 0.0);
    assert_eq!(
        project(GeoPosition::new(0.0, 0.0), projection, area),
        Some((100.0, 100.0))
    );
    assert!(project(GeoPosition::new(0.0, 179.0), projection, area).is_none());
}

#[test]
fn invalid_coordinates_are_counted_and_named_on_the_page() {
    let map = Map::new()
        .push(GeoLocation::new("valid", 0.0, 0.0))
        .push(GeoLocation::new("invalid", 120.0, 0.0));
    assert_eq!(map.invalid_location_count(), 1);
    let svg = map.to_svg();
    assert!(svg.contains(">1 invalid location</text>"), "{svg}");
    assert!(svg.contains("1 invalid locations and 0 unresolved flows"));
}

#[test]
fn location_tooltips_keep_exact_coordinates_and_values() {
    let svg = Map::new()
        .push(
            GeoLocation::new("Lima", -12.046, -77.043)
                .category("Peru")
                .value(42.5)
                .count(7),
        )
        .to_svg();
    assert!(svg.contains(
            "<title>Lima; latitude -12.046; longitude -77.043; 7 observations; category Peru; value 42.5</title>"
        ));
}

#[test]
fn categories_use_shape_as_well_as_colour() {
    let svg = Map::new()
        .push(GeoLocation::new("A", 0.0, 0.0).category("one"))
        .push(GeoLocation::new("B", 10.0, 10.0).category("two"))
        .to_svg();
    assert!(svg.contains("<circle"), "{svg}");
    assert!(svg.contains("<rect"), "{svg}");
    assert!(svg.contains(">one</text>"), "{svg}");
    assert!(svg.contains(">two</text>"), "{svg}");
}

#[test]
fn flows_are_named_and_direction_is_visible() {
    let svg = Map::new()
        .push(GeoLocation::new("Lima", -12.046, -77.043))
        .push(GeoLocation::new("Madrid", 40.417, -3.704))
        .push_flow(GeoFlow::new("Lima", "Madrid").weight(4.0))
        .to_svg();
    assert!(svg.contains("<title>flow, source Lima, target Madrid, weight 4</title>"));
    assert!(svg.contains(" Q "), "{svg}");
    assert!(svg.contains("<polygon"), "arrowhead: {svg}");
}

#[test]
fn a_missing_or_ambiguous_flow_endpoint_is_not_guessed() {
    let missing = Map::new()
        .push(GeoLocation::new("A", 0.0, 0.0))
        .push_flow(GeoFlow::new("A", "B"));
    assert_eq!(missing.unresolved_flow_count(), 1);
    let ambiguous = Map::new()
        .push(GeoLocation::new("A", 0.0, 0.0))
        .push(GeoLocation::new("A", 1.0, 1.0))
        .push(GeoLocation::new("B", 2.0, 2.0))
        .push_flow(GeoFlow::new("A", "B"));
    assert_eq!(ambiguous.unresolved_flow_count(), 1);
}

#[test]
fn a_tooltip_value_does_not_run_off_into_three_hundred_digits() {
    // `f64::MAX` printed in full is three hundred and nine digits, of
    // which sixteen were measured and the rest are the formatter filling
    // in the gap between what an f64 holds and what a decimal can write.
    for extreme in [f64::MAX, f64::MIN, f64::MIN_POSITIVE, 1e308, -1e-308] {
        let text = data_number(extreme);
        let longest = text
            .split(|c: char| !c.is_ascii_digit())
            .map(str::len)
            .max()
            .unwrap_or(0);
        assert!(longest <= 20, "{extreme} came out as {text}");
    }
    // Everything a map is actually likely to carry stays plain.
    assert_eq!(data_number(40.417), "40.417");
    assert_eq!(data_number(-3.704), "-3.704");
    assert_eq!(data_number(0.0), "0");
    assert_eq!(data_number(1234.5), "1234.5");
    assert_eq!(data_number(f64::NAN), "not finite");
}

#[test]
fn base_land_is_compiled_into_every_projection_without_non_finite_output() {
    for projection in [
        GeoProjection::Equirectangular,
        GeoProjection::Mercator,
        GeoProjection::orthographic(20.0, 0.0),
    ] {
        let svg = Map::new().projection(projection).to_svg();
        assert!(svg.matches("<path").count() > 100, "{projection:?}");
        assert!(!svg.contains("NaN"), "{projection:?}");
        assert!(!svg.contains("inf"), "{projection:?}");
    }
}

#[test]
fn rendering_is_deterministic_and_accessible() {
    let map = Map::new()
        .title("Sites")
        .description("Three surveillance sites.")
        .push(GeoLocation::new("A", 0.0, 0.0));
    assert_eq!(map.to_svg(), map.to_svg());
    let svg = map.to_svg();
    assert!(svg.contains("role=\"img\""));
    assert!(svg.contains("<title id=\"karyon-title\">Sites</title>"));
    assert!(svg.contains("<desc id=\"karyon-desc\">Three surveillance sites.</desc>"));
}

#[test]
fn phylo_map_preserves_the_tree_and_aggregates_exact_location_counts() {
    let tree = phylo_tree();
    let leaves = tree.leaf_names();
    let map = PhyloMap::new(tree)
        .location_by("country")
        .coordinates(phylo_coordinates());
    assert_eq!(map.tree().leaf_names(), leaves);
    let svg = map.to_svg();
    assert!(svg.contains("<title>Peru; 2 mapped tips</title>"), "{svg}");
    assert!(svg.contains("<title>Spain; 1 mapped tip</title>"), "{svg}");
    assert!(svg.contains("<title>Kenya; 1 mapped tip</title>"), "{svg}");
    assert!(svg.contains("<title>Peru; latitude -9.19; longitude -75.0152; 2 mapped tips</title>"));
}

#[test]
fn individual_phylo_connectors_name_every_sample() {
    let svg = PhyloMap::new(phylo_tree())
        .location_by("country")
        .coordinates(phylo_coordinates())
        .connector(PhyloConnector::Individual)
        .to_svg();
    for title in [
        "A; location Peru",
        "B; location Peru",
        "C; location Spain",
        "D; location Kenya",
    ] {
        assert!(svg.contains(&format!("<title>{title}</title>")), "{title}");
    }
}

#[test]
fn unlocated_phylo_tips_are_counted_instead_of_dropped() {
    let map = PhyloMap::new(phylo_tree())
        .location_by("country")
        .coordinates(phylo_coordinates().into_iter().take(2));
    assert_eq!(map.unlocated_tip_count(), 1);
    let svg = map.to_svg();
    assert!(svg.contains(">1 unlocated tip</text>"), "{svg}");
    assert!(svg.contains("<title>D; country Kenya</title>"), "{svg}");
}

#[test]
fn ambiguous_coordinate_names_are_not_guessed() {
    let map = PhyloMap::new(phylo_tree())
        .location_by("country")
        .coordinates([
            GeoLocation::new("Peru", -9.19, -75.0152),
            GeoLocation::new("Peru", -12.046, -77.043),
            GeoLocation::new("Spain", 40.4637, -3.7492),
            GeoLocation::new("Kenya", -0.0236, 37.9062),
        ]);
    assert_eq!(map.unlocated_tip_count(), 2);
    assert!(map.to_svg().contains(">2 unlocated tips</text>"));
}

#[test]
fn hidden_phylo_locations_are_reported() {
    let tree = Tree::parse_annotated_newick("(A[&place=Far]:1);").unwrap();
    let svg = PhyloMap::new(tree)
        .location_by("place")
        .projection(GeoProjection::orthographic(0.0, 0.0))
        .coordinate(GeoLocation::new("Far", 0.0, 179.0))
        .to_svg();
    assert!(
        svg.contains(">1 mapped location outside projection</text>"),
        "{svg}"
    );
}

#[test]
fn phylo_time_guides_keep_exact_values_and_units() {
    let tree = Tree::parse_annotated_newick(
        "(A[&date=2023,place=Near]:1,B[&date=2025,place=Near]:1)[&date=2022];",
    )
    .unwrap();
    let svg = PhyloMap::new(tree)
        .location_by("place")
        .coordinate(GeoLocation::new("Near", 0.0, 0.0))
        .time("date")
        .time_unit("year")
        .to_svg();
    for value in ["2022 year", "2023.5 year", "2025 year"] {
        assert!(svg.contains(&format!(">{value}</text>")), "{value}: {svg}");
    }
}

#[test]
fn an_incomplete_phylo_time_layout_is_explicit() {
    let tree =
        Tree::parse_annotated_newick("(A[&date=2023,place=Near]:1,B[&place=Near]:1);").unwrap();
    let svg = PhyloMap::new(tree)
        .location_by("place")
        .coordinate(GeoLocation::new("Near", 0.0, 0.0))
        .time("date")
        .to_svg();
    assert!(
        svg.contains(">requested time layout unavailable</text>"),
        "{svg}"
    );
}

#[test]
fn phylo_map_rendering_is_deterministic_accessible_and_finite() {
    let map = PhyloMap::new(phylo_tree())
        .title("Circular surveillance")
        .description("Four synthetic samples at three supplied locations.")
        .location_by("country")
        .coordinates(phylo_coordinates());
    assert_eq!(map.to_svg(), map.to_svg());
    let svg = map.to_svg();
    assert!(svg.contains("role=\"img\""));
    assert!(svg.contains("<title id=\"karyon-title\">Circular surveillance</title>"));
    assert!(svg.contains(
        "<desc id=\"karyon-desc\">Four synthetic samples at three supplied locations.</desc>"
    ));
    assert!(!svg.contains("NaN"));
    assert!(!svg.contains("inf"));
}

#[test]
fn profiles_scale_map_themes_once() {
    let map = Map::new().profile(RenderProfile::Presentation);
    let phylo = PhyloMap::new(phylo_tree()).profile(RenderProfile::Presentation);
    assert_eq!(map.theme.title_font_size, Theme::light().title_font_size);
    assert_eq!(phylo.theme.title_font_size, Theme::light().title_font_size);
    assert_eq!(map.visual_scale, 1.35);
    assert_eq!(phylo.visual_scale, 1.35);
}

#[test]
fn a_coastline_outside_the_window_is_not_written_into_the_document() {
    // Two thirds of the world's coastline falls outside the window a map draws,
    // and it used to be written out in full and then thrown away by the clip:
    // the committed map figure carried 996 of 1,536 paths whose bounding box
    // missed the clip entirely, which is 646 kilobytes of geometry that drew
    // nothing. `touches` is what keeps them out, so this pins both directions.
    let area = MapRect {
        x: 100.0,
        y: 100.0,
        w: 200.0,
        h: 200.0,
    };
    let hairline = 1.0;

    // Inside, overlapping an edge, and touching a corner: all of them draw.
    assert!(draw::touches("M150 150L250 250Z", area, hairline));
    assert!(draw::touches("M50 150L150 150Z", area, hairline));
    assert!(draw::touches("M300 300L400 400Z", area, hairline));

    // Clear of every edge in each direction: none of them can.
    assert!(!draw::touches("M0 0L50 50Z", area, hairline));
    assert!(!draw::touches("M400 150L500 250Z", area, hairline));
    assert!(!draw::touches("M150 0L250 50Z", area, hairline));
    assert!(!draw::touches("M150 400L250 500Z", area, hairline));

    // A hairline's width outside is still drawn, because the stroke reaches in.
    assert!(draw::touches("M150 99.5L160 99.5Z", area, hairline));

    // A path with no coordinates in it cannot be shown to reach anything.
    assert!(!draw::touches("Z", area, hairline));
}
