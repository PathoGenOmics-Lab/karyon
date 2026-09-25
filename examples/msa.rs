//! Renders the alignment figures used in the README.
//!
//! ```text
//! cargo run --example msa -- assets
//! ```
//!
//! The coordinates here are alignment columns, not genomic positions, so the
//! region spans the width of the alignment and the ruler counts columns.
//!
//! The conservation logo above the alignment is built from the same sequences
//! the alignment track draws, which is the point of both taking a list of rows.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/msa.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        ("example-msa.svg", figures::example_msa(&light, None, None)),
        (
            "example-msa-protein.svg",
            figures::example_msa_protein(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
