//! Renders the selection figure used in the README.
//!
//! ```text
//! cargo run --example selection -- assets
//! ```
//!
//! Two statistics that only mean anything relative to a line. pN/pS says
//! whether a stretch of gene is under purifying or diversifying selection, and
//! the answer is which side of one it fell on; GC skew changes sign at the
//! origin and the terminus of replication, and the sign is the whole signal.
//! Both are plotted as departures from a baseline rather than as heights above
//! a floor, which is what a coverage track would make of them.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The figure is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/selection.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let figure = figures::example_selection(&Theme::light(), None, None);
    fs::write(
        out.join("example-selection.svg"),
        figure.to_svg_with_id_prefix(""),
    )?;
    let (width, height) = figure.dimensions();
    println!("example-selection.svg {width:.0} x {height:.0}");
    Ok(())
}
