//! Renders the copy number figure used in the documentation.
//!
//! ```text
//! cargo run --example copy_number -- assets
//! ```
//!
//! One tumour's segmentation over a chromosome arm, and above it how often the
//! same arm was gained and lost across the cohort it came from. The two bands
//! are the two questions a copy number analysis is asked: what happened in this
//! sample, and does it happen in general.
//!
//! The cohort landscape is a [`WindowTrack`](karyon::WindowTrack) with two rows
//! over every span, one for the gains and one for the losses, because a locus
//! can be gained in a third of a cohort and lost in a fifth of it at the same
//! time and one signed number cannot say that. It is not a track type of its
//! own, and the reason it is not is that this already draws it.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/copy_number.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_copy_number(&Theme::light(), None, None);
    fs::write(
        out.join("example-copy-number.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-copy-number.svg {width:.0} x {height:.0}");
    Ok(())
}
