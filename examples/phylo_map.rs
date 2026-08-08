//! Renders circular phylogeographic compositions from one synthetic tree.
//!
//! ```text
//! cargo run --example phylo_map -- assets
//! ```

use std::env;
use std::path::PathBuf;

use karyon::{GeoLocation, GeoProjection, Panels, PhyloConnector, PhyloMap, Tree, TreeShape};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));
    let tree = outbreak_tree();

    let aggregated = PhyloMap::new(tree.clone())
        .title("Calendar tree around geographic context")
        .subtitle("12 synthetic samples; one connector per reported location")
        .description(
            "A dated synthetic phylogeny around an orthographic map, with three exact supplied coordinates and one aggregated connector per location.",
        )
        .location_by("country")
        .coordinates(locations())
        .projection(GeoProjection::orthographic(10.0, -18.0))
        .time("date")
        .time_unit("year")
        .diameter(700.0)
        .radial_start(-96.0)
        .radial_sweep(350.0);

    let individual = PhyloMap::new(tree)
        .title("A fan with individual sample links")
        .subtitle("The same topology as a cladogram; no migration is inferred")
        .description(
            "The same synthetic samples shown as a circular cladogram, with one named connector per terminal sample and a partial angular sweep.",
        )
        .location_by("country")
        .coordinates(locations())
        .projection(GeoProjection::orthographic(10.0, -18.0))
        .connector(PhyloConnector::Individual)
        .shape(TreeShape::Cladogram)
        .diameter(700.0)
        .radial_start(-220.0)
        .radial_sweep(280.0)
        .show_tip_labels(true);

    let gallery = Panels::new()
        .title("Circular phylogeography (synthetic)")
        .columns(2)
        .gap(26.0)
        .push_captioned(
            &aggregated,
            "A",
            "Calendar radius, supplied place annotations and aggregated links",
        )
        .push_captioned(
            &individual,
            "B",
            "A partial cladogram with every sample-to-location link visible",
        );

    gallery.save_svg(out.join("example-phylo-map.svg"))?;
    let (width, height) = gallery.dimensions();
    println!("example-phylo-map.svg {width:.0} x {height:.0}, 2 circular map views");
    Ok(())
}

fn locations() -> [GeoLocation; 3] {
    [
        GeoLocation::new("Peru", -9.19, -75.0152),
        GeoLocation::new("Spain", 40.4637, -3.7492),
        GeoLocation::new("Kenya", -0.0236, 37.9062),
    ]
}

fn outbreak_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "[&R] (",
        "((PER_001[&date=2023.10,country=Peru]:0.18,",
        "PER_002[&date=2023.24,country=Peru]:0.32)",
        "PER_A[&date=2022.92,country=Peru]:0.42,",
        "(PER_003[&date=2023.68,country=Peru]:0.38,",
        "PER_004[&date=2023.91,country=Peru]:0.61)",
        "PER_B[&date=2023.30,country=Peru]:0.38)",
        "PER_outbreak[&date=2022.50,country=Peru]:0.65,",
        "((ESP_001[&date=2023.02,country=Spain]:0.30,",
        "ESP_002[&date=2023.44,country=Spain]:0.72,",
        "ESP_003[&date=2023.75,country=Spain]:1.03)",
        "ESP_outbreak[&date=2022.72,country=Spain]:0.52,",
        "(KEN_001[&date=2023.20,country=Kenya]:0.26,",
        "(KEN_002[&date=2023.63,country=Kenya]:0.33,",
        "KEN_003[&date=2024.08,country=Kenya]:0.78)",
        "KEN_B[&date=2023.30,country=Kenya]:0.36)",
        "KEN_outbreak[&date=2022.94,country=Kenya]:0.74)",
        "regional[&date=2022.20]:0.35)",
        "origin[&date=2021.85];"
    ))
    .expect("the annotated tree in this example is well formed")
}
