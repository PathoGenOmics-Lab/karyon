//! Renders the association figures used in the README.
//!
//! ```text
//! cargo run --example association -- assets
//! ```
//!
//! A locus that came out of a genome-wide association scan, and the isolates
//! behind it. The Manhattan panel says where the signal is; the matrix says who
//! carries it, which is the question the Manhattan panel always provokes.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/association.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_association(&Theme::light(), None, None);
    fs::write(
        out.join("example-association.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-association.svg {width:.0} x {height:.0}");
    Ok(())
}
