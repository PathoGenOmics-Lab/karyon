//! Renders the ideogram figures used in the README.
//!
//! ```text
//! cargo run --example ideogram -- assets
//! ```
//!
//! An ideogram is a context strip: it sits above the detail tracks and says
//! where in the chromosome they are. The banding below is illustrative rather
//! than a real cytogenetic map; a real one comes straight from the UCSC
//! `cytoBand` table, whose `gieStain` column
//! [`Stain::from_name`](karyon::Stain::from_name) reads.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/ideogram.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-ideogram.svg",
            figures::example_ideogram(&light, None, None),
        ),
        (
            "example-ideogram-bacterial.svg",
            figures::example_ideogram_bacterial(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
