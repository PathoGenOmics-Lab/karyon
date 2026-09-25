//! Renders the sequence logo figures used in the README.
//!
//! ```text
//! cargo run --example logo -- assets
//! ```
//!
//! The point of the first figure is the comparison. The same eight columns are
//! drawn three ways, and only the third one can say anything about the symbols
//! that are missing.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/logo.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-logo.svg",
            figures::example_logo(&light, None, None),
        ),
        (
            "example-logo-protein.svg",
            figures::example_logo_protein(&light, None, None),
        ),
        (
            "example-logo-scores.svg",
            figures::example_logo_scores(&light, None, None),
        ),
        (
            "example-logo-stability.svg",
            figures::example_logo_stability(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
