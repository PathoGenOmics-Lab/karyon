//! Renders the read pileup figure used in the README.
//!
//! ```text
//! cargo run --example pileup -- assets
//! ```
//!
//! A synthetic locus built to look like the thing you open a genome browser
//! for: a variant carried by about half the reads, a deletion carried by a
//! few, an insertion in one, and a patch of low mapping quality. The coverage
//! track is computed from the reads themselves, so the two tracks cannot
//! disagree.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/pileup.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_pileup(&Theme::light(), None, None);
    fs::write(
        out.join("example-pileup.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-pileup.svg {width:.0} x {height:.0}");
    Ok(())
}
