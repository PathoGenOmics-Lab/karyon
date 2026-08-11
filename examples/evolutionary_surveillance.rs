//! Geometry, molecular evolution and genomic surveillance on one sheet.
//!
//! ```text
//! cargo run --example evolutionary_surveillance -- assets
//! ```
//!
//! Every value is synthetic. The example is a visual integration test for
//! upstream results from codon models, ancestral reconstruction, time-tree
//! workflows, tree comparison and lineage-frequency estimation.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{
    AncestralStateLayer, AxisTrack, BranchEventLayer, BranchGeometry, BranchIntervalLayer,
    BranchRateMixture, Feature, FeatureTrack, Figure, HomoplasyLayer, Panels, PhylodynamicPoint,
    PhylodynamicScale, PhylodynamicTrack, Region, SelectionEvidence, SelectionSite, SelectionTrack,
    SurveillanceObservation, SurveillanceStyle, SurveillanceTrack, TangleTieStyle, TanglegramTrack,
    TraitColumn, TreeTrack, Variant, VariantTrack,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let orthogonal = geometry_panel(BranchGeometry::Orthogonal, "Orthogonal phylogram");
    let diagonal = geometry_panel(BranchGeometry::Diagonal, "Diagonal phylogram");
    let curved = ancestral_panel();
    let circular = circular_evolution_panel();
    let unrooted = unrooted_evolution_panel();
    let comparison = comparison_panel();
    let genomic = genomic_selection_panel();
    let surveillance = surveillance_panel();

    let sheet = Panels::new()
        .title("Evolutionary genomics and surveillance atlas (synthetic)")
        .columns(2)
        .row_major()
        .align_plot_areas(false)
        .gap(27.0)
        .column_gap(34.0)
        .push_captioned(
            &orthogonal,
            "A",
            "Orthogonal branches preserve row alignment and make terminal comparison fastest.",
        )
        .push_captioned(
            &diagonal,
            "B",
            "Diagonal edges emphasise topology while retaining branch length and exact annotations.",
        )
        .push_captioned(
            &curved,
            "C",
            "Curves carry posterior ancestral states, inferred transitions, mutation lists and concordance intervals.",
        )
        .push_captioned(
            &circular,
            "D",
            "Circular selection combines mean omega, fitted rate mixtures, recurrence and terminal metadata rings.",
        )
        .push_captioned(
            &unrooted,
            "E",
            "The same direct branch evidence survives a topology-balanced unrooted projection.",
        )
        .push_captioned(
            &comparison,
            "F",
            "A tanglegram exposes core-versus-accessory discordance as crossings without altering either tree.",
        )
        .push_captioned(
            &genomic,
            "G",
            "Protein architecture, selected codons and observed variants meet on one genomic coordinate axis.",
        )
        .push_captioned(
            &surveillance,
            "H",
            "A coalescent trajectory, lineage composition and growth alerts share the same time pivots.",
        );

    let path = out.join("example-evolutionary-surveillance.svg");
    sheet.save_svg(&path)?;
    let (width, height) = sheet.dimensions();
    println!("{} {width:.0} x {height:.0}", path.display());
    Ok(())
}

fn tree_region() -> Region {
    Region::new("phylogeny", 0, 1).unwrap()
}

fn geometry_panel(geometry: BranchGeometry, title: &str) -> Figure {
    Figure::new(tree_region())
        .title(title)
        .width(650.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(evolution_tree())
                .branch_geometry(geometry)
                .row_height(31.0)
                .show_nodes(true)
                .support_style(karyon::SupportStyle::Symbols)
                .trait_categorical("host")
                .scale_bar()
                .scale_bar_unit("substitutions / codon"),
        )
}

fn ancestral_panel() -> Figure {
    Figure::new(tree_region())
        .title("Curved ancestral-state reconstruction")
        .width(650.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(evolution_tree())
                .branch_geometry(BranchGeometry::Curved)
                .row_height(39.0)
                .ancestral_states(
                    AncestralStateLayer::new(["state_human", "state_animal", "state_water"])
                        .label("ancestral host posterior")
                        .confidence(0.72)
                        .size(8.5),
                )
                .branch_event_layer(
                    BranchEventLayer::new("mutations")
                        .label("ancestral mutations")
                        .maximum_events(5),
                )
                .branch_interval(
                    BranchIntervalLayer::new("gcf", "gcf_low", "gcf_high")
                        .label("gene concordance")
                        .threshold(0.70),
                ),
        )
}

fn circular_evolution_panel() -> Figure {
    let mixture = BranchRateMixture::new(
        ["omega_1", "omega_2", "omega_3"],
        ["weight_1", "weight_2", "weight_3"],
    )
    .label("branch-site omega classes")
    .width(27.0)
    .saturation(7.0);
    Figure::new(tree_region())
        .title("Circular molecular evolution")
        .width(650.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(evolution_tree())
                .circular()
                .radial_start(-103.0)
                .radial_size(510.0)
                .inner_radius(0.14)
                .dnds("omega")
                .dnds_label("mean branch omega")
                .dnds_significance("q", 0.05)
                .branch_rate_mixture(mixture)
                .homoplasy_layer(HomoplasyLayer::new("change").label("recurrent amino-acid change"))
                .trait_column(
                    TraitColumn::categorical("host")
                        .label("Host")
                        .ring_width(11.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(9.0),
                ),
        )
}

fn unrooted_evolution_panel() -> Figure {
    Figure::new(tree_region())
        .title("Unrooted mutation and concordance map")
        .width(650.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(evolution_tree())
                .unrooted()
                .unrooted_start(-108.0)
                .unrooted_size(510.0)
                .branch_event_layer(BranchEventLayer::new("mutations").label("ancestral mutations"))
                .branch_interval(
                    BranchIntervalLayer::new("gcf", "gcf_low", "gcf_high")
                        .label("gene concordance")
                        .threshold(0.70),
                )
                .homoplasy("change")
                .trait_column(
                    TraitColumn::categorical("host")
                        .label("Host")
                        .ring_width(11.0),
                ),
        )
}

fn comparison_panel() -> Figure {
    let core = Tree::parse_annotated_newick(
        "(((S1[&lineage=L1]:1,S2[&lineage=L1]:1):1,(S3[&lineage=L2]:1,S4[&lineage=L2]:1):1):1,(S5[&lineage=L3]:1,S6[&lineage=L3]:1):2);",
    )
    .unwrap();
    let accessory = Tree::parse_annotated_newick(
        "(((S1[&lineage=L1]:1,S4[&lineage=L2]:1):1,(S3[&lineage=L2]:1,S6[&lineage=L3]:1):1):1,(S5[&lineage=L3]:1,S2[&lineage=L1]:1):2);",
    )
    .unwrap();
    Figure::new(Region::new("taxa", 0, 6).unwrap())
        .title("Core and accessory tree comparison")
        .width(650.0)
        .show_region_label(false)
        .push(
            TanglegramTrack::new(core, accessory)
                .names("core genome", "accessory genome")
                .color_by("lineage")
                .tie_style(TangleTieStyle::Ribbon)
                .untangle()
                .row_height(31.0),
        )
}

fn genomic_selection_panel() -> Figure {
    let region = Region::new("surface protein", 0, 300).unwrap();
    Figure::new(region)
        .title("Genomic layers around selected codons")
        .width(650.0)
        .show_region_label(false)
        .push(
            FeatureTrack::new(vec![
                Feature::new(4, 68).name("signal peptide"),
                Feature::new(79, 176).name("receptor-binding domain"),
                Feature::new(187, 248).name("stalk"),
                Feature::new(260, 293).name("TM helix"),
            ])
            .label("Protein"),
        )
        .push(
            SelectionTrack::new(site_selection())
                .label("FEL / MEME")
                .evidence(SelectionEvidence::PValue)
                .p_threshold(0.05)
                .neutral_band(0.85, 1.15)
                .saturation(8.0)
                .height(158.0),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(45).value(0.31).category("S45N recurrent"),
                Variant::new(88).value(0.18).category("escape candidate"),
                Variant::new(226).value(0.44).category("lineage marker"),
            ])
            .label("Observed"),
        )
        .push(AxisTrack::new().center_on_bases(true))
}

fn surveillance_panel() -> Figure {
    Figure::new(Region::new("month", 0, 13).unwrap())
        .title("Phylodynamics and lineage surveillance")
        .width(650.0)
        .show_region_label(false)
        .push(
            PhylodynamicTrack::new(phylodynamic_points())
                .label("Effective size")
                .scale(PhylodynamicScale::Log10)
                .unit("Ne")
                .height(128.0),
        )
        .push(
            SurveillanceTrack::new(surveillance_observations())
                .label("Lineages")
                .style(SurveillanceStyle::Stacked)
                .minimum_total(100)
                .frequency_alert(0.25)
                .growth_alert(0.08)
                .height(150.0),
        )
        .push(AxisTrack::new().center_on_bases(true).label("Month"))
}

fn site_selection() -> Vec<SelectionSite> {
    vec![
        SelectionSite::new(18).rates(0.44, 0.08).p_value(0.41),
        SelectionSite::new(45)
            .rates(0.18, 1.52)
            .p_value(0.0014)
            .episodic_rates(0.05, 3.8, 0.18)
            .label("recurrent S45N"),
        SelectionSite::new(72).rates(0.38, 0.14).p_value(0.22),
        SelectionSite::new(88)
            .rates(0.24, 1.13)
            .p_value(0.018)
            .episodic_rates(0.09, 2.7, 0.11),
        SelectionSite::new(104).rates(0.50, 0.07).p_value(0.008),
        SelectionSite::new(139)
            .rates(0.16, 0.91)
            .p_value(0.031)
            .episodic_rates(0.04, 4.2, 0.08),
        SelectionSite::new(181).rates(0.38, 0.42).p_value(0.71),
        SelectionSite::new(203).rates(0.29, 0.06).p_value(0.003),
        SelectionSite::new(226)
            .rates(0.22, 1.65)
            .p_value(0.006)
            .episodic_rates(0.07, 5.1, 0.14),
        SelectionSite::new(271).rates(0.20, 0.24).p_value(0.61),
    ]
}

fn phylodynamic_points() -> Vec<PhylodynamicPoint> {
    [
        (0, 180.0, 110.0, 280.0),
        (1, 220.0, 130.0, 360.0),
        (2, 310.0, 180.0, 520.0),
        (3, 480.0, 260.0, 850.0),
        (4, 760.0, 390.0, 1_300.0),
        (5, 1_200.0, 620.0, 2_100.0),
        (6, 2_100.0, 1_100.0, 3_800.0),
        (7, 3_600.0, 1_900.0, 6_400.0),
        (8, 4_800.0, 2_500.0, 8_100.0),
        (9, 4_100.0, 2_200.0, 7_000.0),
        (10, 3_200.0, 1_700.0, 5_600.0),
        (11, 2_700.0, 1_400.0, 4_800.0),
        (12, 2_300.0, 1_200.0, 4_100.0),
    ]
    .into_iter()
    .map(|(time, estimate, lower, upper)| {
        PhylodynamicPoint::new(time, estimate).interval(lower, upper)
    })
    .collect()
}

fn surveillance_observations() -> Vec<SurveillanceObservation> {
    let frequencies = [
        (0, 170, 25, 5),
        (1, 163, 30, 7),
        (2, 154, 37, 9),
        (3, 143, 45, 12),
        (4, 129, 56, 15),
        (5, 112, 69, 19),
        (6, 94, 82, 24),
        (7, 76, 94, 30),
        (8, 61, 101, 38),
        (9, 48, 103, 49),
        (10, 37, 99, 64),
        (11, 29, 91, 80),
        (12, 23, 80, 97),
    ];
    let mut observations = Vec::new();
    for (time, l1, l2, l3) in frequencies {
        observations.push(SurveillanceObservation::new(time, "L1", l1, 200));
        observations.push(SurveillanceObservation::new(time, "L2", l2, 200));
        observations.push(SurveillanceObservation::new(time, "L3", l3, 200));
    }
    observations
}

fn evolution_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "(((S1[&host=Human,resistant=false,omega=0.2,q=0.4,change=S45N,mutations={S_S45N,nuc_A135G},gcf=0.88,gcf_low=0.78,gcf_high=0.95,state_human=0.96,state_animal=0.02,state_water=0.02,omega_1=0.08,weight_1=0.78,omega_2=1,weight_2=0.18,omega_3=8,weight_3=0.04]:0.12,",
        "S2[&host=Human,resistant=false,omega=0.4,q=0.3,change=V88I,mutations={S_V88I},gcf=0.74,gcf_low=0.61,gcf_high=0.85,state_human=0.91,state_animal=0.04,state_water=0.05,omega_1=0.12,weight_1=0.66,omega_2=1,weight_2=0.28,omega_3=5,weight_3=0.06]:0.13)",
        "H[&omega=3.1,q=0.004,change=G45D,mutations={S_G45D},gcf=0.91,gcf_low=0.84,gcf_high=0.96,state_human=0.90,state_animal=0.04,state_water=0.06,omega_1=0.16,weight_1=0.52,omega_2=1,weight_2=0.24,omega_3=11,weight_3=0.24]:0.20,",
        "(S3[&host=Wastewater,resistant=true,omega=1.7,q=0.025,change=S45N,mutations={S_S45N,AMR_gain},gcf=0.69,gcf_low=0.55,gcf_high=0.81,state_human=0.08,state_animal=0.04,state_water=0.88,omega_1=0.2,weight_1=0.55,omega_2=1.2,weight_2=0.25,omega_3=7.8,weight_3=0.20]:0.16,",
        "S4[&host=Wastewater,resistant=true,omega=0.6,q=0.14,change=T201I,mutations={S_T201I},gcf=0.62,gcf_low=0.48,gcf_high=0.76,state_human=0.04,state_animal=0.03,state_water=0.93,omega_1=0.1,weight_1=0.70,omega_2=0.8,weight_2=0.25,omega_3=4.4,weight_3=0.05]:0.15)",
        "W[&state_human=0.06,state_animal=0.04,state_water=0.90,gcf=0.73,gcf_low=0.60,gcf_high=0.84]:0.18)",
        "HW[&state_human=0.61,state_animal=0.05,state_water=0.34,gcf=0.79,gcf_low=0.69,gcf_high=0.88]:0.22,",
        "(S5[&host=Livestock,resistant=false,omega=0.14,q=0.62,change=L19F,mutations={S_L19F},gcf=0.83,gcf_low=0.72,gcf_high=0.91,state_human=0.03,state_animal=0.94,state_water=0.03,omega_1=0.06,weight_1=0.83,omega_2=0.95,weight_2=0.14,omega_3=6,weight_3=0.03]:0.13,",
        "S6[&host=Livestock,resistant=false,omega=5.8,q=0.002,change=S45N,mutations={S_S45N,S_E226K},gcf=0.57,gcf_low=0.42,gcf_high=0.71,state_human=0.02,state_animal=0.95,state_water=0.03,omega_1=0.18,weight_1=0.41,omega_2=1,weight_2=0.23,omega_3=13,weight_3=0.36]:0.18)",
        "A[&state_human=0.03,state_animal=0.93,state_water=0.04,gcf=0.67,gcf_low=0.54,gcf_high=0.79]:0.28)",
        "root[&state_human=0.42,state_animal=0.34,state_water=0.24];"
    ))
    .expect("the synthetic evolutionary tree is valid")
}
