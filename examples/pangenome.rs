//! Renders the pangenome figure used in the README.
//!
//! ```text
//! cargo run --example pangenome -- assets
//! ```
//!
//! A presence and absence matrix of accessory genes, with the isolates sorted
//! by the phylogeny beside them. Sorted by sample name the same data is a
//! speckle; sorted by descent the accessory regions become rectangles, and a
//! rectangle is a claim about the biology.
//!
//! The isolates are *Klebsiella pneumoniae*, which has an open pangenome and
//! moves it around. The organism is load bearing rather than decorative: a
//! clonal species with no horizontal transfer has almost no accessory genome,
//! so this figure would have nothing in it.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/pangenome.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_pangenome(&Theme::light(), None, None);
    fs::write(
        out.join("example-pangenome.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-pangenome.svg {width:.0} x {height:.0}");
    Ok(())
}
