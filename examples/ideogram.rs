//! Renders the ideogram figures used in the README.
//!
//! ```text
//! cargo run --example ideogram -- assets
//! ```
//!
//! An ideogram is a context strip: it sits above the detail tracks and says
//! where in the chromosome they are. The banding below is illustrative rather
//! than a real cytogenetic map; a real one comes straight from the UCSC
//! `cytoBand` table, whose `gieStain` column [`Stain::from_name`] reads.

use std::env;
use std::path::PathBuf;

use karyon::{
    AxisTrack, Band, CoverageTrack, Feature, FeatureTrack, Figure, IdeogramTrack, Region, Stain,
    Strand, Variant, VariantTrack,
};

/// Length of the illustrative chromosome.
const CHROMOSOME: u64 = 48_000_000;
/// The window the detail tracks show.
const WINDOW: (u64, u64) = (31_200_000, 31_260_000);

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

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

    let region = Region::new("chr7", WINDOW.0, WINDOW.1).unwrap();
    let context = Figure::new(region)
        .title("Sixty kilobases of a chromosome, and where they are")
        .push(
            IdeogramTrack::new(CHROMOSOME, bands)
                .label("chr7")
                .show_band_names(true)
                .height(36.0),
        )
        .push(
            CoverageTrack::new(WINDOW.0, depth)
                .label("depth")
                .height(52.0),
        )
        .push(
            FeatureTrack::new(vec![
                Feature::new(31_210_000, 31_236_000)
                    .name("GENEA")
                    .strand(Strand::Forward),
                Feature::new(31_241_000, 31_255_000)
                    .name("GENEB")
                    .strand(Strand::Reverse),
            ])
            .label("genes"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(31_218_400).value(0.91).category("missense"),
                Variant::new(31_247_900).value(0.44).category("splice"),
            ])
            .label("variants")
            .height(42.0),
        )
        .push(AxisTrack::new());
    context.save_svg(out.join("example-ideogram.svg"))?;

    // A bacterial chromosome has no cytogenetics to speak of, so the ideogram
    // is a bare outline. It still answers the only question it was ever asked.
    let h37rv_length = 4_411_532;
    let rpo_b = Region::new("NC_000962.3", 759_806, 763_325).unwrap();
    let bacterial = Figure::new(rpo_b)
        .title("Mycobacterium tuberculosis H37Rv, rpoB")
        .width(760.0)
        .push(
            IdeogramTrack::bare(h37rv_length)
                .label("H37Rv")
                .height(20.0),
        )
        .push(
            FeatureTrack::new(vec![Feature::new(759_806, 763_325)
                .name("rpoB")
                .strand(Strand::Forward)])
            .label("gene"),
        )
        .push(AxisTrack::new());
    bacterial.save_svg(out.join("example-ideogram-bacterial.svg"))?;

    let (w1, h1) = context.dimensions();
    let (w2, h2) = bacterial.dimensions();
    println!("example-ideogram.svg           {w1:.0} x {h1:.0}");
    println!("example-ideogram-bacterial.svg {w2:.0} x {h2:.0}");
    Ok(())
}
