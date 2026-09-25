//! Renders the regulation figure used in the documentation.
//!
//! ```text
//! cargo run --example regulation -- assets
//! ```
//!
//! What an RNA-seq experiment and a model of it say about the same stretch:
//! how deep the reads lay, which introns they stepped over and how many took
//! each one, and which bases a sequence model leaned on when it predicted the
//! signal.
//!
//! The three bands answer three different questions and are drawn three
//! different ways on purpose. Depth is a height above a floor. A junction has
//! no height at all, so its arcs are put in lanes and the count is printed on
//! them. And an attribution is a signed number carried by a base, so the base
//! itself is the mark and it hangs below the line where the model pulled the
//! other way.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/regulation.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::example_regulation(&Theme::light(), None, None);
    fs::write(
        out.join("example-regulation.svg"),
        sheet.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = sheet.dimensions();
    println!("example-regulation.svg {width:.0} x {height:.0}");
    Ok(())
}
