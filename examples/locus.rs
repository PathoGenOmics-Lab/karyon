//! Renders the two figures used in the README.
//!
//! Run it with the output directory as an optional argument:
//!
//! ```text
//! cargo run --example locus -- assets
//! ```
//!
//! The depth and the reference bases are synthetic but shaped like the real
//! thing: a coverage dropout, and a GC content near the 65% this genome sits at.
//! The annotation and the variant positions are not synthetic. `rpoB` spans the
//! whole window and runs off both edges, the RRDR box is the real 81 base
//! hotspot, and every variant sits on the codon it is named for, which is worth
//! the trouble because a figure with an invented gene in it teaches the reader
//! the wrong thing about the tool as well as about the locus.
//!
//! Everything generated is generated from a fixed seed, so re-running the
//! example produces byte-identical files and a diff only appears when the
//! rendering actually changed.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/locus.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        ("example.svg", figures::example(&light, None, None)),
        (
            "example-zoom.svg",
            figures::example_zoom(&light, None, None),
        ),
        (
            "example-dark.svg",
            figures::example_dark(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
