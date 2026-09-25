//! Renders the circular genome figure used in the README.
//!
//! ```text
//! cargo run --example circular -- assets
//! ```
//!
//! A bacterial chromosome the way it actually is. Drawn as a line, a circular
//! genome gets an edge the biology does not have, straight through whatever
//! happens to sit at coordinate zero. Drawn as a circle it does not, and the
//! middle becomes usable: the chords join the two ends of the rearrangements.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/circular.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_circular(&Theme::light(), None, None);
    fs::write(
        out.join("example-circular.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-circular.svg {width:.0} x {height:.0}");
    Ok(())
}
