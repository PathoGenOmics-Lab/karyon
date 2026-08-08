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
use karyon::{
    Figure, Panels, RadialDirection, Region, SupportStyle, TraitColumn, TreeShape, TreeTrack,
};

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

    let advanced = advanced_layouts(outbreak_tree());
    advanced.save_svg(out.join("example-phylo-annotations.svg"))?;
    let (width, height) = advanced.dimensions();
    println!(
        "example-phylo-annotations.svg {width:.0} x {height:.0}, unrooted plus 4 annotation datasets"
    );

    let evidence = evidence_layouts(evidence_tree());
    evidence.save_svg(out.join("example-phylo-evidence.svg"))?;
    let (width, height) = evidence.dimensions();
    println!(
        "example-phylo-evidence.svg {width:.0} x {height:.0}, support plus branch events and distance scales"
    );
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

fn advanced_layouts(tree: Tree) -> Panels {
    let annotated_unrooted = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Unrooted topology with metadata halo")
        .width(690.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree.clone())
                .unrooted()
                .unrooted_start(-104.0)
                .unrooted_size(560.0)
                .color_by("country")
                .show_nodes(true)
                .trait_column(
                    TraitColumn::categorical("country")
                        .label("Country strip")
                        .ring_width(11.0)
                        .show_values(false),
                )
                .trait_column(
                    TraitColumn::bar("coverage")
                        .label("Depth bars")
                        .ring_width(22.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(13.0),
                )
                .trait_column(TraitColumn::symbol("host").label("Host").ring_width(15.0)),
        );

    let circular_annotations = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("iTOL-style annotation rings")
        .width(690.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(tree)
                .shape(TreeShape::Cladogram)
                .circular()
                .radial_start(-104.0)
                .radial_size(560.0)
                .color_by("country")
                .show_nodes(true)
                .trait_column(
                    TraitColumn::categorical("country")
                        .label("Country strip")
                        .ring_width(11.0)
                        .show_values(false),
                )
                .trait_column(
                    TraitColumn::bar("coverage")
                        .label("Depth bars")
                        .ring_width(22.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(13.0),
                )
                .trait_column(TraitColumn::symbol("host").label("Host").ring_width(15.0)),
        );

    Panels::new()
        .title("Unrooted trees and layered metadata (synthetic)")
        .columns(2)
        .gap(22.0)
        .push_captioned(
            &annotated_unrooted,
            "A",
            "A topology-balanced centre removes the arbitrary Newick root while leaders preserve unequal branch lengths",
        )
        .push_captioned(
            &circular_annotations,
            "B",
            "Colour strips, radial bars, binary marks and shaped categories retain exact values in tooltips",
        )
}

fn evidence_track(tree: Tree) -> TreeTrack {
    TreeTrack::new(tree)
        .color_by("lineage")
        .support_style(SupportStyle::SymbolsAndLabels)
        .support_threshold(0.70)
        .branch_labels("event")
        .branch_label_size(7.0)
        .scale_bar()
        .scale_bar_length(0.1)
        .scale_bar_unit("substitutions/site")
}

fn evidence_layouts(tree: Tree) -> Panels {
    let rectangular = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Evidence-rich phylogram")
        .width(530.0)
        .show_region_label(false)
        .push(evidence_track(tree.clone()).row_height(28.0));

    let circular = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Circular evidence view")
        .width(530.0)
        .show_region_label(false)
        .push(
            evidence_track(tree.clone())
                .circular()
                .radial_start(-112.0)
                .radial_size(500.0),
        );

    let unrooted = Figure::new(Region::new("phylogeny", 0, 1).unwrap())
        .title("Unrooted evidence view")
        .width(530.0)
        .show_region_label(false)
        .push(
            evidence_track(tree)
                .unrooted()
                .unrooted_start(-112.0)
                .unrooted_size(500.0),
        );

    Panels::new()
        .title("Branch evidence across phylogenetic projections (synthetic)")
        .columns(3)
        .gap(18.0)
        .push_captioned(
            &rectangular,
            "A",
            "Support, branch events and evolutionary distance remain independently readable",
        )
        .push_captioned(
            &circular,
            "B",
            "Event labels rotate with their branch while exact text remains in tooltips",
        )
        .push_captioned(
            &unrooted,
            "C",
            "The scale still measures edge length after removing the arbitrary source root",
        )
}

fn evidence_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "[&R] (",
        "((A01[&lineage=L1,event=rpoB-S450L]:0.18,A02[&lineage=L1]:0.14)",
        "0.98[&lineage=L1,event=katG-S315T]:0.22,",
        "(A03[&lineage=L2,event=gyrA-D94G]:0.21,A04[&lineage=L2]:0.16)",
        "0.84[&lineage=L2]:0.19)0.93:0.17,",
        "((B01[&lineage=L3,event=del-pks15]:0.13,B02[&lineage=L3]:0.20)",
        "0.76[&lineage=L3]:0.24,",
        "(B03[&lineage=L4,event=embB-M306V]:0.17,B04[&lineage=L4]:0.23)",
        "0.64[&lineage=L4]:0.16)0.88:0.21);"
    ))
    .expect("the evidence tree in this example is well formed")
}

fn outbreak_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "[&R] (",
        "((PER_001[&date=2023.10,country=Peru,coverage=48,resistant=true,host=human]:0.18,",
        "PER_002[&date=2023.24,country=Peru,coverage=73,resistant=true,host=human]:0.32)",
        "PER_A[&date=2022.92,country=Peru]:0.42,",
        "(PER_003[&date=2023.68,country=Peru,coverage=31,resistant=false,host=human]:0.38,",
        "PER_004[&date=2023.91,country=Peru,coverage=59,resistant=true,host=animal]:0.61)",
        "PER_B[&date=2023.30,country=Peru]:0.38)",
        "PER_outbreak[&date=2022.50,country=Peru]:0.65,",
        "((ESP_001[&date=2023.02,country=Spain,coverage=66,resistant=false,host=human]:0.30,",
        "ESP_002[&date=2023.44,country=Spain,coverage=41,resistant=false,host=water]:0.72,",
        "ESP_003[&date=2023.75,country=Spain,coverage=84,resistant=true,host=human]:1.03)",
        "ESP_outbreak[&date=2022.72,country=Spain]:0.52,",
        "(KEN_001[&date=2023.20,country=Kenya,coverage=28,resistant=false,host=animal]:0.26,",
        "(KEN_002[&date=2023.63,country=Kenya,coverage=52,resistant=true,host=human]:0.33,",
        "KEN_003[&date=2024.08,country=Kenya,coverage=77,resistant=true,host=water]:0.78)",
        "KEN_B[&date=2023.30,country=Kenya]:0.36)",
        "KEN_outbreak[&date=2022.94,country=Kenya]:0.74)",
        "regional[&date=2022.20]:0.35)",
        "origin[&date=2021.85];"
    ))
    .expect("the annotated tree in this example is well formed")
}
