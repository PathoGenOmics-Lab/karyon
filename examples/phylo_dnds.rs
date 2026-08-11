//! Branch-wise dN/dS across three phylogenetic projections.
//!
//! ```text
//! cargo run --example phylo_dnds -- assets
//! ```
//!
//! The tree and every annotation are synthetic. They are designed to exercise
//! the visual grammar, not to make a biological claim about a real lineage.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{Figure, Panels, Region, TraitColumn, TreeTrack};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));
    let tree = selection_tree();

    let rectangular = Figure::new(Region::new("codon model", 0, 1).unwrap())
        .title("Where selection changes along the tree")
        .width(700.0)
        .show_region_label(false)
        .push(
            selection_track(tree.clone())
                .row_height(43.0)
                .branch_labels("event")
                .branch_label_size(7.0)
                .scale_bar()
                .scale_bar_unit("substitutions / codon")
                .trait_column(TraitColumn::categorical("host").label("Host").width(56.0))
                .trait_column(TraitColumn::binary("resistant").label("AMR").width(38.0)),
        );

    let circular = Figure::new(Region::new("codon model", 0, 1).unwrap())
        .title("Selection with a genomic metadata halo")
        .width(700.0)
        .show_region_label(false)
        .push(
            selection_track(tree.clone())
                .circular()
                .radial_start(-105.0)
                .radial_size(560.0)
                .inner_radius(0.12)
                .trait_column(
                    TraitColumn::categorical("host")
                        .label("Host")
                        .ring_width(12.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(10.0),
                )
                .trait_column(
                    TraitColumn::continuous("identity")
                        .label("Protein identity")
                        .ring_width(12.0)
                        .show_values(false),
                ),
        );

    let unrooted = Figure::new(Region::new("codon model", 0, 1).unwrap())
        .title("The same branch estimates without a privileged root")
        .width(700.0)
        .show_region_label(false)
        .push(
            selection_track(tree)
                .unrooted()
                .unrooted_start(-112.0)
                .unrooted_size(560.0)
                .trait_column(
                    TraitColumn::categorical("host")
                        .label("Host")
                        .ring_width(12.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(10.0),
                ),
        );

    let cladogram = Figure::new(Region::new("codon model", 0, 1).unwrap())
        .title("Selection on topology alone")
        .width(640.0)
        .show_region_label(false)
        .push(
            selection_track(selection_tree())
                .shape(karyon::TreeShape::Cladogram)
                .row_height(43.0)
                .branch_labels("event")
                .branch_label_size(7.0),
        );

    let sheet = Panels::new()
        .title("Branch-wise molecular evolution (synthetic)")
        .columns(2)
        .row_major()
        .align_plot_areas(false)
        .gap(24.0)
        .column_gap(34.0)
        .push_captioned(
            &rectangular,
            "A",
            "A cool-to-warm logarithmic scale centres the branch model on ω = 1; thick branches pass q ≤ 0.05.",
        )
        .push_captioned(
            &circular,
            "B",
            "Host, resistance and protein identity remain separate terminal datasets around the selection tree.",
        )
        .push_captioned(
            &unrooted,
            "C",
            "Direct branch estimates survive the unrooted projection; dotted edges make missing fits visible.",
        )
        .push_captioned(
            &cladogram,
            "D",
            "A cladogram removes branch-length magnitude without changing the ω estimate owned by each edge.",
        );

    let path = out.join("example-phylo-dnds.svg");
    sheet.save_svg(&path)?;
    let (width, height) = sheet.dimensions();
    println!(
        "{} {width:.0} x {height:.0}, {} terminal genomes",
        path.display(),
        selection_tree().leaf_count()
    );
    Ok(())
}

fn selection_track(tree: Tree) -> TreeTrack {
    TreeTrack::new(tree)
        .dnds("omega")
        .dnds_label("Branch dN/dS (ω)")
        .dnds_neutral_band(0.9, 1.1)
        .dnds_saturation(4.0)
        .dnds_significance("q", 0.05)
        .line_width(1.5)
}

fn selection_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "(((H01[&omega=0.18,q=0.42,event=S12A,host=Human,resistant=false,identity=99.8]:0.10,",
        "H02[&omega=0.31,q=0.21,event=V88I,host=Human,resistant=false,identity=99.6]:0.12)",
        "[&omega=2.85,q=0.008,event=G45D]:0.18,",
        "(H03[&omega=1.42,q=0.038,event=A117T,host=Human,resistant=true,identity=98.9]:0.16,",
        "H04[&omega=0.72,q=0.16,event=T201I,host=Wastewater,resistant=true,identity=98.7]:0.15)",
        "[&omega=1.18,q=0.08,event=P33S]:0.20)",
        "[&omega=3.60,q=0.002,event=R164K]:0.24,",
        "((H05[&omega=0.14,q=0.62,event=L19F,host=Livestock,resistant=false,identity=97.9]:0.13,",
        "H06[&q=0.91,event=Q71H,host=Livestock,resistant=false,identity=97.7]:0.14)",
        "[&omega=0.24,q=0.44,event=I52V]:0.19,",
        "(H07[&omega=4.80,q=0.004,event=D96N,host=Wildlife,resistant=true,identity=96.8]:0.18,",
        "(H08[&omega=0.41,q=0.32,event=E143D,host=Wildlife,resistant=false,identity=96.5]:0.11,",
        "H09[&omega=1.02,q=0.77,event=A150V,host=Wastewater,resistant=false,identity=96.2]:0.12)",
        "[&omega=0.58,q=0.19,event=V132I]:0.13)",
        "[&omega=1.64,q=0.028,event=K102R]:0.22)",
        "[&omega=0.36,q=0.18,event=N7S]:0.30);"
    ))
    .expect("the synthetic branch-model tree is valid")
}
