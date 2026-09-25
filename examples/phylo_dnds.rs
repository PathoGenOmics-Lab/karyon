//! Branch-wise dN/dS across three phylogenetic projections.
//!
//! ```text
//! cargo run --example phylo_dnds -- assets
//! ```
//!
//! The tree and every annotation are synthetic. They are designed to exercise
//! the visual grammar, not to make a biological claim about a real lineage.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/phylo_dnds.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::example_phylo_dnds(&Theme::light(), None, None);
    let path = out.join("example-phylo-dnds.svg");
    fs::write(&path, sheet.to_svg_with_id_prefix(""))?;
    let (width, height) = sheet.dimensions();
    println!("{} {width:.0} x {height:.0}", path.display());
    Ok(())
}
