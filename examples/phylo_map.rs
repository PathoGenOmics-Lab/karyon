//! Renders circular phylogeographic compositions from one synthetic tree.
//!
//! ```text
//! cargo run --example phylo_map -- assets
//! ```

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/phylo_map.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let gallery = figures::example_phylo_map(&Theme::light(), None, None);
    fs::write(
        out.join("example-phylo-map.svg"),
        gallery.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = gallery.dimensions();
    println!("example-phylo-map.svg {width:.0} x {height:.0}");
    Ok(())
}
