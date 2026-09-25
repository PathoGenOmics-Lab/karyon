//! Renders the deterministic geographic gallery used in the documentation.
//!
//! ```text
//! cargo run --example maps -- assets
//! ```

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/maps.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let gallery = figures::example_maps(&Theme::light(), None, None);
    fs::write(
        out.join("example-maps.svg"),
        gallery.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = gallery.dimensions();
    println!("example-maps.svg {width:.0} x {height:.0}");
    Ok(())
}
