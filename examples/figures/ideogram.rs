//! The figures `cargo run --example ideogram` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever window the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::{Band, Drawing, Feature, IdeogramTrack, Plot, Region, Stain, Strand, Theme, Variant};

/// Length of the illustrative chromosome.
const CHROMOSOME: u64 = 48_000_000;
/// The window the detail tracks show.
const WINDOW: (u64, u64) = (31_200_000, 31_260_000);

/// `example-ideogram.svg`: sixty kilobases of a chromosome, and where they are.
///
/// `theme`, `width` and `region` replace the light theme, the 900 pixels and
/// the sixty kilobase window the committed figure is drawn with. The ideogram
/// marks whichever window that is.
pub fn example_ideogram(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let bands = vec![
        Band::new(0, 3_100_000, Stain::Gneg).name("p13"),
        Band::new(3_100_000, 6_400_000, Stain::Gpos50).name("p12"),
        Band::new(6_400_000, 9_800_000, Stain::Gneg).name("p11.2"),
        Band::new(9_800_000, 12_200_000, Stain::Gpos75).name("p11.1"),
        Band::new(12_200_000, 13_900_000, Stain::Acen),
        Band::new(13_900_000, 17_500_000, Stain::Gneg).name("q11"),
        Band::new(17_500_000, 22_100_000, Stain::Gpos100).name("q12"),
        Band::new(22_100_000, 26_000_000, Stain::Gneg).name("q13"),
        Band::new(26_000_000, 29_400_000, Stain::Gpos25).name("q21.1"),
        Band::new(29_400_000, 33_000_000, Stain::Gpos75).name("q21.2"),
        Band::new(33_000_000, 36_800_000, Stain::Gneg).name("q22"),
        Band::new(36_800_000, 41_500_000, Stain::Gpos100).name("q23"),
        Band::new(41_500_000, 44_200_000, Stain::Gvar).name("q24"),
        Band::new(44_200_000, 48_000_000, Stain::Gneg).name("q25"),
    ];

    let span = (WINDOW.1 - WINDOW.0) as usize;
    let depth: Vec<f64> = (0..span)
        .map(|i| 42.0 + 11.0 * ((i as f64) / 900.0).sin())
        .collect();

    let own = Region::new("chr7", WINDOW.0, WINDOW.1).unwrap();
    let context = Plot::over(region.cloned().unwrap_or(own))
        .title("Sixty kilobases of a chromosome, and where they are")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .add_ideogram(CHROMOSOME, bands)
        .label("chr7")
        .adjust(|track| track.show_band_names(true).height(36.0))
        // From the first base the depth was measured at rather than from the
        // left edge of the window, which is only the same place until the
        // window moves.
        .add_coverage_at(WINDOW.0, depth)
        .label("depth")
        .adjust(|track| track.height(52.0))
        .add_features(vec![
            Feature::new(31_210_000, 31_236_000)
                .name("GENEA")
                .strand(Strand::Forward),
            Feature::new(31_241_000, 31_255_000)
                .name("GENEB")
                .strand(Strand::Reverse),
        ])
        .label("genes")
        .add_variants(vec![
            Variant::new(31_218_400).value(0.91).category("missense"),
            Variant::new(31_247_900).value(0.44).category("splice"),
        ])
        .label("variants")
        .adjust(|track| track.height(42.0).axis_title("AF"))
        .into_figure();
    Box::new(context)
}

/// `example-ideogram-bacterial.svg`: one gene, and where it is on a
/// chromosome with no cytogenetics to speak of.
///
/// `theme`, `width` and `region` replace the light theme, the 760 pixels and
/// the gene the committed figure is drawn over.
pub fn example_ideogram_bacterial(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // A bacterial chromosome has no cytogenetics to speak of, so the ideogram
    // is a bare outline. It still answers the only question it was ever asked.
    let h37rv_length = 4_411_532;
    let rpo_b = Region::new("NC_000962.3", 759_806, 763_325).unwrap();
    let bacterial = Plot::over(region.cloned().unwrap_or(rpo_b))
        .title("Mycobacterium tuberculosis H37Rv, rpoB")
        .theme(theme.clone())
        .width(width.unwrap_or(760.0))
        .add_track(
            IdeogramTrack::bare(h37rv_length)
                .label("H37Rv")
                .height(20.0),
        )
        .add_features(vec![Feature::new(759_806, 763_325)
            .name("rpoB")
            .strand(Strand::Forward)])
        .label("gene")
        .into_figure();
    Box::new(bacterial)
}
