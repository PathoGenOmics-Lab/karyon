//! The figures `cargo run --example synteny` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever window the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::{AlignmentBlock, Drawing, Feature, Plot, Region, Strand, Theme};

/// Length of the query chromosome.
const QUERY: u64 = 4_400_000;
/// Length of the target chromosome.
const TARGET: u64 = 4_380_000;

/// `example-synteny.svg`: two chromosomes, three disagreements.
///
/// `theme`, `width` and `region` replace the light theme, the 900 pixels and
/// the whole query chromosome the committed figure is drawn over. The region
/// moves along the query; the target is drawn whole whatever it says.
pub fn example_synteny(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let blocks = blocks();
    let own = Region::new("H37Rv", 0, QUERY).unwrap();
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("Two chromosomes, three disagreements")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .add_dotplot(blocks.clone())
        .label("CDC1551")
        .adjust(|track| track.target_length(TARGET).height(210.0))
        .add_synteny(blocks)
        .adjust(|track| {
            track
                .target_length(TARGET)
                .names("H37Rv", "CDC1551")
                .height(120.0)
        })
        .into_figure();
    Box::new(figure)
}

/// `example-synteny-inversion.svg`: the inversion, close up.
///
/// `theme`, `width` and `region` replace the light theme, the 760 pixels and
/// the stretch of the query the committed figure is drawn over. The target
/// keeps its own stretch, the one the inversion is in.
pub fn example_synteny_inversion(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // Zoomed onto the inversion, where the ribbon crossing is the whole story.
    let inversion = Region::new("H37Rv", 1_400_000, 2_250_000).unwrap();
    let detail = Plot::over(region.cloned().unwrap_or(inversion))
        .title("The inversion, close up")
        .theme(theme.clone())
        .width(width.unwrap_or(760.0))
        .add_synteny(blocks())
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
    Box::new(detail)
}

/// Colinear for the first third, then an inversion, then a block that moved,
/// then colinear again. The gap between the third and fourth blocks on the
/// query is sequence the target does not have.
fn blocks() -> Vec<AlignmentBlock> {
    vec![
        AlignmentBlock::new(0, 1_500_000, 0, 1_500_000).identity(0.99),
        AlignmentBlock::new(1_520_000, 2_100_000, 1_520_000, 2_100_000)
            .reversed(true)
            .identity(0.97),
        AlignmentBlock::new(2_150_000, 2_600_000, 3_400_000, 3_850_000).identity(0.98),
        AlignmentBlock::new(2_900_000, 4_400_000, 2_150_000, 3_650_000).identity(0.99),
    ]
}
