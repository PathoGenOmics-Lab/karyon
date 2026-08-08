//! Annotated, dated phylogenies for genomic surveillance.
//!
//! ```text
//! cargo run --example phylogenetics -- assets
//! ```
//!
//! The outbreak is synthetic. Its labels and metadata are fixed so the figure
//! is a deterministic visual regression target rather than an epidemiological
//! claim.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{Figure, Panels, RadialDirection, Region, TraitColumn, TreeShape, TreeTrack};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let tree = outbreak_tree();
    let sample_count = tree.leaf_count();
    let full = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Dated transmission context")
        .width(720.0)
        .show_region_label(false)
        .push(annotated_track(tree.clone()));

    let mut summary_tree = tree;
    summary_tree.ladderize(true);
    let peru = summary_tree.node_named("PER_outbreak").unwrap();
    let kenya = summary_tree.node_named("KEN_outbreak").unwrap();
    let summary = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Collapsed surveillance view")
        .width(720.0)
        .show_region_label(false)
        .push(annotated_track(summary_tree).collapse(peru).collapse(kenya));

    let sheet = Panels::new()
        .title("Annotated outbreak phylogeny (synthetic)")
        .columns(2)
        .gap(22.0)
        .push_captioned(
            &full,
            "A",
            "Calendar time, branch location and aligned sample metadata",
        )
        .push_captioned(
            &summary,
            "B",
            "Named clades collapse without changing the source topology",
        );
    sheet.save_svg(out.join("example-phylogenetics.svg"))?;
    let (width, height) = sheet.dimensions();
    println!(
        "example-phylogenetics.svg {width:.0} x {height:.0}, {} samples",
        sample_count
    );

    let radial = radial_layouts(outbreak_tree());
    radial.save_svg(out.join("example-phylo-layouts.svg"))?;
    let (width, height) = radial.dimensions();
    println!("example-phylo-layouts.svg {width:.0} x {height:.0}, 4 projections");
    Ok(())
}

fn annotated_track(tree: Tree) -> TreeTrack {
    base_track(tree)
        .trait_column(
            TraitColumn::categorical("country")
                .label("Country")
                .width(62.0)
                .ring_width(12.0),
        )
        .trait_column(
            TraitColumn::continuous("coverage")
                .label("Depth")
                .width(46.0)
                .ring_width(12.0),
        )
}

fn base_track(tree: Tree) -> TreeTrack {
    TreeTrack::new(tree)
        .time("date")
        .time_unit("year")
        .color_by("country")
        .show_nodes(true)
        .row_height(18.0)
}

fn radial_layouts(tree: Tree) -> Panels {
    let full = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Circular time tree")
        .width(650.0)
        .show_region_label(false)
        .push(annotated_track(tree.clone()).circular().radial_size(440.0));

    let mut fan_tree = tree.clone();
    fan_tree.ladderize(true);
    let peru = fan_tree.node_named("PER_outbreak").unwrap();
    let fan = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Collapsed 250° fan")
        .width(650.0)
        .show_region_label(false)
        .push(
            annotated_track(fan_tree)
                .collapse(peru)
                .fan(250.0)
                .radial_start(-215.0)
                .radial_size(520.0)
                .show_time_axis(false),
        );

    let inward = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Time radiating inwards")
        .width(650.0)
        .show_region_label(false)
        .push(
            base_track(tree.clone())
                .circular()
                .radial_direction(RadialDirection::Inward)
                .inner_radius(0.34)
                .radial_size(440.0)
                .show_tips(false),
        );

    let cladogram = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Circular cladogram")
        .width(650.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .shape(TreeShape::Cladogram)
                .circular()
                .radial_start(-72.0)
                .radial_size(440.0)
                .color_by("country")
                .show_nodes(true)
                .trait_column(
                    TraitColumn::categorical("country")
                        .label("Country")
                        .ring_width(12.0),
                ),
        );

    Panels::new()
        .title("Circular and radial phylogenies (synthetic)")
        .columns(2)
        .gap(22.0)
        .push_captioned(
            &full,
            "A",
            "Outward calendar radii with categorical and continuous trait rings",
        )
        .push_captioned(
            &fan,
            "B",
            "A partial sweep keeps room for annotation and a collapsed clade",
        )
        .push_captioned(
            &inward,
            "C",
            "The same dated topology points towards a controlled central gap",
        )
        .push_captioned(
            &cladogram,
            "D",
            "Branch counts replace lengths and every tip reaches one circumference",
        )
}

fn outbreak_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "[&R] (",
        "((PER_001[&date=2023.10,country=Peru,coverage=48]:0.18,",
        "PER_002[&date=2023.24,country=Peru,coverage=73]:0.32)",
        "PER_A[&date=2022.92,country=Peru]:0.42,",
        "(PER_003[&date=2023.68,country=Peru,coverage=31]:0.38,",
        "PER_004[&date=2023.91,country=Peru,coverage=59]:0.61)",
        "PER_B[&date=2023.30,country=Peru]:0.38)",
        "PER_outbreak[&date=2022.50,country=Peru]:0.65,",
        "((ESP_001[&date=2023.02,country=Spain,coverage=66]:0.30,",
        "ESP_002[&date=2023.44,country=Spain,coverage=41]:0.72,",
        "ESP_003[&date=2023.75,country=Spain,coverage=84]:1.03)",
        "ESP_outbreak[&date=2022.72,country=Spain]:0.52,",
        "(KEN_001[&date=2023.20,country=Kenya,coverage=28]:0.26,",
        "(KEN_002[&date=2023.63,country=Kenya,coverage=52]:0.33,",
        "KEN_003[&date=2024.08,country=Kenya,coverage=77]:0.78)",
        "KEN_B[&date=2023.30,country=Kenya]:0.36)",
        "KEN_outbreak[&date=2022.94,country=Kenya]:0.74)",
        "regional[&date=2022.20]:0.35)",
        "origin[&date=2021.85];"
    ))
    .expect("the annotated tree in this example is well formed")
}
