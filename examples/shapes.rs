//! Renders the four figures whose shape belongs to genomics and nowhere else.
//!
//! ```text
//! cargo run --example shapes -- assets
//! ```
//!
//! Structural variants as arcs between their breakpoints, the six reading
//! frames with their stops, two trees face to face, and methylation one
//! molecule at a time. Four organisms, because the forms are not about any of
//! them.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/shapes.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-structural.svg",
            figures::example_structural(&light, None, None),
        ),
        (
            "example-frames.svg",
            figures::example_frames(&light, None, None),
        ),
        (
            "example-tanglegram.svg",
            figures::example_tanglegram(&light, None, None),
        ),
        (
            "example-bisulfite.svg",
            figures::example_bisulfite(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
