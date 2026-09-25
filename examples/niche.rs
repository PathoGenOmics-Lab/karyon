//! Renders the three specialised figures used in the README.
//!
//! ```text
//! cargo run --example niche -- assets
//! ```
//!
//! Three plots that belong to one sub-field each and that a general purpose
//! plotting library will never have: raw nanopore current, a gene cluster
//! compared across genomes, and per-strand methylation.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figures are built in a file of their own, which the documentation
// site's playground includes as well, so the committed SVGs and the ones drawn
// live in the page come out of the same code.
#[path = "figures/niche.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let light = Theme::light();
    for (file, figure) in [
        (
            "example-squiggle.svg",
            figures::example_squiggle(&light, None, None),
        ),
        (
            "example-cluster.svg",
            figures::example_cluster(&light, None, None),
        ),
        (
            "example-methylation.svg",
            figures::example_methylation(&light, None, None),
        ),
    ] {
        fs::write(out.join(file), figure.to_svg_with_id_prefix(""))?;
        let (width, height) = figure.dimensions();
        println!("{file} {width:.0} x {height:.0}");
    }
    Ok(())
}
