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
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheets are built in a file of their own, which the documentation site's
// playground includes as well, so the committed SVGs and the ones drawn live
// in the page come out of the same code.
#[path = "figures/phylogenetics.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, sheet) in [
        (
            "example-phylogenetics.svg",
            figures::example_phylogenetics(&light, None, None),
        ),
        (
            "example-phylo-layouts.svg",
            figures::example_phylo_layouts(&light, None, None),
        ),
        (
            "example-phylo-annotations.svg",
            figures::example_phylo_annotations(&light, None, None),
        ),
        (
            "example-phylo-evidence.svg",
            figures::example_phylo_evidence(&light, None, None),
        ),
        (
            "example-phylo-reroot.svg",
            figures::example_phylo_reroot(&light, None, None),
        ),
        (
            "example-phylo-faces.svg",
            figures::example_phylo_faces(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), sheet.to_svg_with_id_prefix(""))?;
        let (width, height) = sheet.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
