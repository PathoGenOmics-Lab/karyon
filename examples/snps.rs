//! Renders the variable site panel used in the README.
//!
//! ```text
//! cargo run --example snps -- assets
//! ```
//!
//! Twelve isolates against a reference over thirty kilobases. Almost all of it
//! agrees, so drawing all thirty thousand columns would spend every pixel on
//! the agreement. The panel keeps only the columns that vary and spaces them
//! evenly, which is why each one carries its own position underneath: the
//! spacing deliberately says nothing about distance.
//!
//! What is known about the isolates is drawn between their names and the panel,
//! out of the sample sheet in `figures/snps.rs`. Lineage runs in three blocks
//! because the phylogeny put the rows in that order, and resistance does not,
//! which is the thing the two strips side by side are for: the same phenotype
//! in two clades that did not inherit it from each other.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/snps.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_snps(&Theme::light(), None, None);
    fs::write(
        out.join("example-snps.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-snps.svg {width:.0} x {height:.0}");
    Ok(())
}
