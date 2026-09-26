//! Most kinds of plot this crate draws, on one sheet.
//!
//! ```text
//! cargo run --example gallery -- assets
//! ```
//!
//! The twenty-three panels draw thirty of the thirty-seven track types
//! between them, and a circular chromosome. The seven with no panel are drawn
//! by other examples: `CopyNumberTrack` by `copy_number`, `DomainTrack` by
//! `phylogenetics`, `DynseqTrack` and `JunctionTrack` by `regulation`, and
//! `PhylodynamicTrack`, `SelectionTrack` and `SurveillanceTrack` by
//! `evolutionary_surveillance`.
//!
//! **When a new track type is added, give it a panel here**, in
//! `figures/gallery.rs`, where the panels are built. An overview that quietly
//! stops covering everything is worse than no overview, because it looks
//! complete.
//!
//! And the test a new track has to pass to earn a panel: **does it live on the
//! genomic coordinate axis?** That is the whole reason this crate exists rather
//! than a general plotting library. A track whose `draw` never reads
//! `ctx.scale` is a bar chart, a line chart or a heatmap that happens to have
//! been handed genomic data, and three of those were removed from this sheet
//! for exactly that reason.
//!
//! The panels do not share a coordinate system with each other, which is the
//! honest reason they are separate figures rather than one tall stack. Five of
//! them are not in genomic coordinates at all: the alignment counts columns,
//! the variable site panel counts sites, the squiggle counts samples of
//! current, and the two tree panels measure evolutionary distance.
//!
//! One is in genomic coordinates twice over. The codon ruler counts residues
//! along the same axis the bases are on, which is the only way to point at a
//! figure and say S450L.

use std::env;
use std::fs;
use std::path::PathBuf;

use karyon::Theme;

// The sheet is built in a file of its own, which the documentation site's
// playground includes as well, so the committed SVG and the one drawn live in
// the page come out of the same code.
#[path = "figures/gallery.rs"]
mod figures;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = figures::gallery(&Theme::light(), None, None);
    fs::write(out.join("gallery.svg"), sheet.to_svg_with_id_prefix(""))?;
    let (width, height) = sheet.dimensions();
    println!("gallery.svg {width:.0} x {height:.0}");
    Ok(())
}
