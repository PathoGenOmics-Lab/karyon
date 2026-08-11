//! A selection atlas joining branch models, recurrent changes and site scans.
//!
//! ```text
//! cargo run --example selection_atlas -- assets
//! ```
//!
//! All data are synthetic and intentionally contain missing fits, recurrent
//! events and heterogeneous rate classes so the visual semantics are visible.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{
    AxisTrack, BranchRateMixture, Feature, FeatureTrack, Figure, HomoplasyLayer, Panels, Region,
    SelectionEvidence, SelectionSite, SelectionTrack, TraitColumn, TreeTrack,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let mixture = BranchRateMixture::new(
        ["omega_slow", "omega_neutral", "omega_burst"],
        ["weight_slow", "weight_neutral", "weight_burst"],
    )
    .label("aBSREL omega classes")
    .width(29.0)
    .thickness(5.6)
    .saturation(6.0);
    let recurrences = HomoplasyLayer::new("change")
        .label("recurrent amino-acid change")
        .width(1.25);

    let rectangular = Figure::new(tree_region())
        .title("Branch mixtures preserve episodic classes")
        .width(700.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(selection_tree())
                .row_height(42.0)
                .line_width(1.25)
                .branch_rate_mixture(mixture.clone())
                .homoplasy_layer(recurrences.clone())
                .trait_column(TraitColumn::categorical("host").label("Host").width(52.0))
                .scale_bar()
                .scale_bar_unit("substitutions / codon"),
        );

    let circular = Figure::new(tree_region())
        .title("Selection and recurrence in circular context")
        .width(700.0)
        .show_region_label(false)
        .push(
            TreeTrack::new(selection_tree())
                .circular()
                .radial_start(-105.0)
                .radial_size(570.0)
                .inner_radius(0.16)
                .dnds("omega")
                .dnds_label("mean branch omega")
                .dnds_significance("q", 0.05)
                .homoplasy_layer(recurrences)
                .trait_column(
                    TraitColumn::categorical("host")
                        .label("Host")
                        .ring_width(12.0),
                )
                .trait_column(
                    TraitColumn::binary("resistant")
                        .label("AMR")
                        .ring_width(9.0),
                ),
        );

    let fel = Figure::new(codon_region())
        .title("Site-wise FEL and episodic MEME evidence")
        .width(700.0)
        .show_region_label(false)
        .push(protein_domains())
        .push(
            SelectionTrack::new(frequentist_sites())
                .label("FEL / MEME")
                .evidence(SelectionEvidence::PValue)
                .p_threshold(0.05)
                .neutral_band(0.85, 1.15)
                .saturation(8.0)
                .height(164.0),
        )
        .push(AxisTrack::new().center_on_bases(true));

    let fubar = Figure::new(codon_region())
        .title("Posterior selection scan over protein architecture")
        .width(700.0)
        .show_region_label(false)
        .push(protein_domains())
        .push(
            SelectionTrack::new(posterior_sites())
                .label("FUBAR")
                .evidence(SelectionEvidence::Posterior)
                .posterior_threshold(0.90)
                .neutral_band(0.85, 1.15)
                .saturation(8.0)
                .height(164.0),
        )
        .push(AxisTrack::new().center_on_bases(true));

    let sheet = Panels::new()
        .title("Molecular selection atlas (synthetic)")
        .columns(2)
        .row_major()
        .align_plot_areas(false)
        .gap(27.0)
        .column_gap(34.0)
        .push_captioned(
            &rectangular,
            "A",
            "Capsule width is fitted class weight; colour is omega. Dashed links mark a recurrent branch event.",
        )
        .push_captioned(
            &circular,
            "B",
            "Mean branch omega, statistical emphasis, terminal metadata and recurrent events use separate channels.",
        )
        .push_captioned(
            &fel,
            "C",
            "Evidence height and signed log2(omega) are separated; diamonds cross the p-value threshold.",
        )
        .push_captioned(
            &fubar,
            "D",
            "The same effect grammar works with posterior probability without pretending it is a p-value.",
        );

    let path = out.join("example-selection-atlas.svg");
    sheet.save_svg(&path)?;
    let (width, height) = sheet.dimensions();
    println!("{} {width:.0} x {height:.0}", path.display());
    Ok(())
}

fn tree_region() -> Region {
    Region::new("codon model", 0, 1).unwrap()
}

fn codon_region() -> Region {
    Region::new("surface protein", 0, 300).unwrap()
}

fn protein_domains() -> FeatureTrack {
    FeatureTrack::new(vec![
        Feature::new(5, 72).name("signal peptide"),
        Feature::new(82, 177).name("receptor-binding domain"),
        Feature::new(188, 248).name("stalk"),
        Feature::new(258, 293).name("TM helix"),
    ])
    .label("Protein")
    .row_height(18.0)
}

fn frequentist_sites() -> Vec<SelectionSite> {
    vec![
        SelectionSite::new(12).rates(0.42, 0.08).p_value(0.54),
        SelectionSite::new(28).rates(0.31, 0.11).p_value(0.22),
        SelectionSite::new(45)
            .rates(0.18, 1.52)
            .p_value(0.0014)
            .episodic_rates(0.05, 3.8, 0.18)
            .label("S45N recurrent site"),
        SelectionSite::new(63).rates(0.37, 0.21).p_value(0.17),
        SelectionSite::new(88)
            .rates(0.24, 1.13)
            .p_value(0.018)
            .episodic_rates(0.09, 2.7, 0.11)
            .label("surface loop"),
        SelectionSite::new(104).rates(0.50, 0.07).p_value(0.008),
        SelectionSite::new(121).rates(0.33, 0.31).p_value(0.46),
        SelectionSite::new(139)
            .rates(0.16, 0.91)
            .p_value(0.031)
            .episodic_rates(0.04, 4.2, 0.08),
        SelectionSite::new(162).rates(0.45, 0.14).p_value(0.12),
        SelectionSite::new(181).rates(0.38, 0.42).p_value(0.71),
        SelectionSite::new(203).rates(0.29, 0.06).p_value(0.003),
        SelectionSite::new(226)
            .rates(0.22, 1.65)
            .p_value(0.006)
            .episodic_rates(0.07, 5.1, 0.14)
            .label("escape candidate"),
        SelectionSite::new(247).rates(0.47, 0.18).p_value(0.29),
        SelectionSite::new(271).rates(0.20, 0.24).p_value(0.61),
        SelectionSite::new(286).rates(0.36, 0.04).p_value(0.044),
    ]
}

fn posterior_sites() -> Vec<SelectionSite> {
    let rates = [
        (12, 0.42, 0.08, 0.88),
        (28, 0.31, 0.11, 0.71),
        (45, 0.18, 1.52, 0.995),
        (63, 0.37, 0.21, 0.77),
        (88, 0.24, 1.13, 0.962),
        (104, 0.50, 0.07, 0.984),
        (121, 0.33, 0.31, 0.56),
        (139, 0.16, 0.91, 0.944),
        (162, 0.45, 0.14, 0.91),
        (181, 0.38, 0.42, 0.42),
        (203, 0.29, 0.06, 0.997),
        (226, 0.22, 1.65, 0.991),
        (247, 0.47, 0.18, 0.74),
        (271, 0.20, 0.24, 0.38),
        (286, 0.36, 0.04, 0.971),
    ];
    rates
        .into_iter()
        .map(|(pos, ds, dn, posterior)| SelectionSite::new(pos).rates(ds, dn).posterior(posterior))
        .collect()
}

fn selection_tree() -> Tree {
    Tree::parse_annotated_newick(concat!(
        "(((H01[&omega=0.20,q=0.42,change=S45N,host=Human,resistant=false,omega_slow=0.08,weight_slow=0.78,omega_neutral=1.1,weight_neutral=0.18,omega_burst=8.2,weight_burst=0.04]:0.12,",
        "H02[&omega=0.44,q=0.31,change=V88I,host=Human,resistant=false,omega_slow=0.12,weight_slow=0.66,omega_neutral=0.9,weight_neutral=0.28,omega_burst=5.4,weight_burst=0.06]:0.13)",
        "AB[&omega=3.10,q=0.004,change=G45D,omega_slow=0.16,weight_slow=0.52,omega_neutral=1.0,weight_neutral=0.24,omega_burst=11.0,weight_burst=0.24]:0.20,",
        "(H03[&omega=1.72,q=0.025,change=S45N,host=Wastewater,resistant=true,omega_slow=0.2,weight_slow=0.55,omega_neutral=1.2,weight_neutral=0.25,omega_burst=7.8,weight_burst=0.20]:0.16,",
        "H04[&omega=0.61,q=0.14,change=T201I,host=Wastewater,resistant=true,omega_slow=0.1,weight_slow=0.70,omega_neutral=0.8,weight_neutral=0.25,omega_burst=4.4,weight_burst=0.05]:0.15):0.18):0.22,",
        "((H05[&omega=0.14,q=0.62,change=L19F,host=Livestock,resistant=false,omega_slow=0.06,weight_slow=0.83,omega_neutral=0.95,weight_neutral=0.14,omega_burst=6.0,weight_burst=0.03]:0.13,",
        "H06[&q=0.91,change=Q71H,host=Livestock,resistant=false]:0.14):0.19,",
        "(H07[&omega=5.80,q=0.002,change=S45N,host=Wildlife,resistant=true,omega_slow=0.18,weight_slow=0.41,omega_neutral=1.0,weight_neutral=0.23,omega_burst=13.0,weight_burst=0.36]:0.18,",
        "H08[&omega=0.38,q=0.28,change=E143D,host=Wildlife,resistant=false,omega_slow=0.09,weight_slow=0.76,omega_neutral=0.85,weight_neutral=0.20,omega_burst=5.9,weight_burst=0.04]:0.12):0.24):0.28);"
    ))
    .expect("the synthetic selection tree is valid")
}
