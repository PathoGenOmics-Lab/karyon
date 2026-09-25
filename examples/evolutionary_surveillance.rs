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
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/evolutionary_surveillance.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::example_evolutionary_surveillance(&Theme::light(), None, None);
    let path = out.join("example-evolutionary-surveillance.svg");
    fs::write(&path, sheet.to_svg_with_id_prefix(""))?;
    let (width, height) = sheet.dimensions();
    println!("{} {width:.0} x {height:.0}", path.display());
    Ok(())
}
