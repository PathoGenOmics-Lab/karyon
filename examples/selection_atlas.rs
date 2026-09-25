//! A selection atlas joining branch models, recurrent changes and site scans.
//!
//! ```text
//! cargo run --example selection_atlas -- assets
//! ```
//!
//! All data are synthetic and intentionally contain missing fits, recurrent
//! events and heterogeneous rate classes so the visual semantics are visible.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/selection_atlas.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::example_selection_atlas(&Theme::light(), None, None);
    let path = out.join("example-selection-atlas.svg");
    fs::write(&path, sheet.to_svg_with_id_prefix(""))?;
    let (width, height) = sheet.dimensions();
    println!("{} {width:.0} x {height:.0}", path.display());
    Ok(())
}
