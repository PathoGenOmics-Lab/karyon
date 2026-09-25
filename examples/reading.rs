//! Four figures about reading a coordinate as something other than a base.
//!
//! ```text
//! cargo run --example reading -- assets
//! ```
//!
//! As a residue, as a place one molecule visited, as a stretch a whole clade
//! shares, and as part of one RNA. Four organisms, because none of the four
//! forms is about any particular one.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/reading.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-codons.svg",
            figures::example_codons(&light, None, None),
        ),
        (
            "example-split.svg",
            figures::example_split(&light, None, None),
        ),
        (
            "example-clades.svg",
            figures::example_clades(&light, None, None),
        ),
        (
            "example-transcripts.svg",
            figures::example_transcripts(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
