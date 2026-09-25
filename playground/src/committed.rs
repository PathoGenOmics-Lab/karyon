//! Every figure the repository commits, as the function that draws it.
//!
//! The examples write the figures under `assets/`, and the functions they
//! write them with live in `examples/figures/`, one file per example. Those
//! files are compiled in here as they stand rather than copied, so a figure a
//! page draws for itself comes out of the code that drew the file it stands in
//! for: the same data, the same seed, the same calls. Only what the page asks
//! for, its theme, its width and its window, is different.
//!
//! The list below is a second copy of what `assets/` holds, and a second copy
//! is a thing that falls behind, so a test compares the two in both
//! directions: a figure committed and not listed here fails it by name, and so
//! does one listed and not committed.

use karyon::{Drawing, Region, Theme};

#[path = "../../examples/figures/association.rs"]
mod association;
#[path = "../../examples/figures/circular.rs"]
mod circular;
#[path = "../../examples/figures/copy_number.rs"]
mod copy_number;
#[path = "../../examples/figures/evolutionary_surveillance.rs"]
mod evolutionary_surveillance;
#[path = "../../examples/figures/gallery.rs"]
mod gallery;
#[path = "../../examples/figures/genomewide.rs"]
mod genomewide;
#[path = "../../examples/figures/ideogram.rs"]
mod ideogram;
#[path = "../../examples/figures/locus.rs"]
mod locus;
#[path = "../../examples/figures/logo.rs"]
mod logo;
#[path = "../../examples/figures/maps.rs"]
mod maps;
#[path = "../../examples/figures/msa.rs"]
mod msa;
#[path = "../../examples/figures/niche.rs"]
mod niche;
#[path = "../../examples/figures/pangenome.rs"]
mod pangenome;
#[path = "../../examples/figures/phylo_dnds.rs"]
mod phylo_dnds;
#[path = "../../examples/figures/phylo_map.rs"]
mod phylo_map;
#[path = "../../examples/figures/phylogenetics.rs"]
mod phylogenetics;
#[path = "../../examples/figures/pileup.rs"]
mod pileup;
#[path = "../../examples/figures/reading.rs"]
mod reading;
#[path = "../../examples/figures/regulation.rs"]
mod regulation;
#[path = "../../examples/figures/selection.rs"]
mod selection;
#[path = "../../examples/figures/selection_atlas.rs"]
mod selection_atlas;
#[path = "../../examples/figures/shapes.rs"]
mod shapes;
#[path = "../../examples/figures/snps.rs"]
mod snps;
#[path = "../../examples/figures/synteny.rs"]
mod synteny;
#[path = "../../examples/figures/visual_system.rs"]
mod visual_system;

/// What draws one figure: in a theme, at a width or its own, and over a window
/// or its own.
pub(crate) type Builder = fn(&Theme, Option<f64>, Option<&Region>) -> Box<dyn Drawing>;

/// Every committed figure, by the name of its file under `assets/` without the
/// `.svg`, sorted by that name.
pub(crate) const FIGURES: &[(&str, Builder)] = &[
    ("example", locus::example),
    ("example-association", association::example_association),
    ("example-bisulfite", shapes::example_bisulfite),
    ("example-circular", circular::example_circular),
    ("example-clades", reading::example_clades),
    ("example-cluster", niche::example_cluster),
    ("example-codons", reading::example_codons),
    ("example-copy-number", copy_number::example_copy_number),
    ("example-dark", locus::example_dark),
    (
        "example-evolutionary-surveillance",
        evolutionary_surveillance::example_evolutionary_surveillance,
    ),
    ("example-frames", shapes::example_frames),
    ("example-genomewide", genomewide::example_genomewide),
    ("example-ideogram", ideogram::example_ideogram),
    (
        "example-ideogram-bacterial",
        ideogram::example_ideogram_bacterial,
    ),
    ("example-logo", logo::example_logo),
    ("example-logo-protein", logo::example_logo_protein),
    ("example-logo-scores", logo::example_logo_scores),
    ("example-logo-stability", logo::example_logo_stability),
    ("example-maps", maps::example_maps),
    ("example-methylation", niche::example_methylation),
    ("example-msa", msa::example_msa),
    ("example-msa-protein", msa::example_msa_protein),
    ("example-pangenome", pangenome::example_pangenome),
    (
        "example-phylo-annotations",
        phylogenetics::example_phylo_annotations,
    ),
    ("example-phylo-dnds", phylo_dnds::example_phylo_dnds),
    (
        "example-phylo-evidence",
        phylogenetics::example_phylo_evidence,
    ),
    ("example-phylo-faces", phylogenetics::example_phylo_faces),
    (
        "example-phylo-layouts",
        phylogenetics::example_phylo_layouts,
    ),
    ("example-phylo-map", phylo_map::example_phylo_map),
    ("example-phylo-reroot", phylogenetics::example_phylo_reroot),
    (
        "example-phylogenetics",
        phylogenetics::example_phylogenetics,
    ),
    ("example-pileup", pileup::example_pileup),
    ("example-regulation", regulation::example_regulation),
    ("example-selection", selection::example_selection),
    (
        "example-selection-atlas",
        selection_atlas::example_selection_atlas,
    ),
    ("example-snps", snps::example_snps),
    ("example-split", reading::example_split),
    ("example-squiggle", niche::example_squiggle),
    ("example-structural", shapes::example_structural),
    ("example-synteny", synteny::example_synteny),
    (
        "example-synteny-inversion",
        synteny::example_synteny_inversion,
    ),
    ("example-tanglegram", shapes::example_tanglegram),
    ("example-transcripts", reading::example_transcripts),
    (
        "example-visual-system",
        visual_system::example_visual_system,
    ),
    ("example-zoom", locus::example_zoom),
    ("gallery", gallery::gallery),
];

/// The function that draws the figure committed as `assets/{stem}.svg`.
pub(crate) fn builder(stem: &str) -> Option<Builder> {
    FIGURES
        .iter()
        .find(|(name, _)| *name == stem)
        .map(|(_, build)| *build)
}
