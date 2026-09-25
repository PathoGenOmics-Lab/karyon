//! Renders the genome-wide figure used in the README.
//!
//! ```text
//! cargo run --example genomewide -- assets
//! ```
//!
//! What a figure over one region cannot show: a whole draft assembly, every
//! contig of it on one axis, with the association scan and the depth profile
//! laid over all of them at once.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/genomewide.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_genomewide(&Theme::light(), None, None);
    fs::write(
        out.join("example-genomewide.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-genomewide.svg {width:.0} x {height:.0}");
    Ok(())
}
