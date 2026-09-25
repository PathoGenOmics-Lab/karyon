//! Renders the visual-system example used by the guide.
//!
//! ```text
//! cargo run --example visual_system -- assets
//! ```

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/visual_system.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::example_visual_system(&Theme::light(), None, None);
    fs::write(
        out.join("example-visual-system.svg"),
        sheet.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = sheet.dimensions();
    println!("example-visual-system.svg {width:.0} x {height:.0}");
    Ok(())
}
