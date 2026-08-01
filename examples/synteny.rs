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
use std::path::PathBuf;

use karyon::{AlignmentBlock, Feature, Plot, Region, Strand};

/// Length of the query chromosome.
const QUERY: u64 = 4_400_000;
/// Length of the target chromosome.
const TARGET: u64 = 4_380_000;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    // Colinear for the first third, then an inversion, then a block that moved,
    // then colinear again. The gap between the third and fourth blocks on the
    // query is sequence the target does not have.
    let blocks = vec![
        AlignmentBlock::new(0, 1_500_000, 0, 1_500_000).identity(0.99),
        AlignmentBlock::new(1_520_000, 2_100_000, 1_520_000, 2_100_000)
            .reversed(true)
            .identity(0.97),
        AlignmentBlock::new(2_150_000, 2_600_000, 3_400_000, 3_850_000).identity(0.98),
        AlignmentBlock::new(2_900_000, 4_400_000, 2_150_000, 3_650_000).identity(0.99),
    ];

    let region = Region::new("H37Rv", 0, QUERY).unwrap();
    let figure = Plot::over(region)
        .title("Two chromosomes, three disagreements")
        .width(900.0)
        .add_dotplot(blocks.clone())
        .label("CDC1551")
        .adjust(|track| track.target_length(TARGET).height(210.0))
        .add_synteny(blocks.clone())
        .adjust(|track| {
            track
                .target_length(TARGET)
                .names("H37Rv", "CDC1551")
                .height(120.0)
        })
        .into_figure();
    figure.save_svg(out.join("example-synteny.svg"))?;

    // Zoomed onto the inversion, where the ribbon crossing is the whole story.
    let inversion = Region::new("H37Rv", 1_400_000, 2_250_000).unwrap();
    let detail = Plot::over(inversion)
        .title("The inversion, close up")
        .width(760.0)
        .add_synteny(blocks)
        .adjust(|track| {
            track
                .target_range(1_400_000, 2_250_000)
                .names("H37Rv", "CDC1551")
                .height(130.0)
        })
        .add_features(vec![Feature::new(1_520_000, 2_100_000)
            .name("inverted segment")
            .strand(Strand::Reverse)])
        .label("segment")
        .into_figure();
    detail.save_svg(out.join("example-synteny-inversion.svg"))?;

    let (w1, h1) = figure.dimensions();
    let (w2, h2) = detail.dimensions();
    println!("example-synteny.svg           {w1:.0} x {h1:.0}");
    println!("example-synteny-inversion.svg {w2:.0} x {h2:.0}");
    println!("{QUERY} bases of query against a {TARGET} base target");
    Ok(())
}
