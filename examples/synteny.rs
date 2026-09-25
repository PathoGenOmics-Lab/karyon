//! Renders the sequence comparison figures used in the README.
//!
//! ```text
//! cargo run --example synteny -- assets
//! ```
//!
//! Two bacterial chromosomes that agree about most of themselves and disagree
//! about three places: an inversion, a translocated block, and a stretch one of
//! them simply does not have. The same blocks are drawn twice, once as a
//! dotplot and once as ribbons, because the two forms answer different halves
//! of the question.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/synteny.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-synteny.svg",
            figures::example_synteny(&light, None, None),
        ),
        (
            "example-synteny-inversion.svg",
            figures::example_synteny_inversion(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
