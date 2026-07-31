//! Renders the sequence logo figures used in the README.
//!
//! ```text
//! cargo run --example logo -- assets
//! ```
//!
//! The point of the first figure is the comparison. The same eight columns are
//! drawn three ways, and only the third one can say anything about the symbols
//! that are missing.

use std::env;
use std::path::PathBuf;

use karyon::{AxisTrack, Figure, LogoColumn, LogoScore, LogoTrack, Region, StackOrder};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    // Counts out of 100 aligned sequences, the shape real data arrives in:
    // mostly non-zero, with one position where a base genuinely never appears.
    // Column 4 is the one to look at. It is almost uniform, so the classic
    // logo says there is nothing here, and it is missing a T entirely, which
    // only the third panel can tell you.
    let motif = vec![
        LogoColumn::acgt(97.0, 1.0, 1.0, 1.0),    // all but fixed
        LogoColumn::acgt(48.0, 46.0, 3.0, 3.0),   // a two base split
        LogoColumn::acgt(26.0, 25.0, 24.0, 25.0), // uniform, nothing at all
        LogoColumn::acgt(34.0, 33.0, 33.0, 0.0),  // near uniform, but no T ever
        LogoColumn::acgt(8.0, 9.0, 74.0, 9.0),    // G rich
        LogoColumn::acgt(2.0, 44.0, 10.0, 44.0),  // A nearly gone
        LogoColumn::acgt(40.0, 30.0, 20.0, 10.0), // a gradient
        LogoColumn::acgt(32.0, 4.0, 32.0, 32.0),  // C rare
    ];

    let region = Region::new("motif", 0, motif.len() as u64).unwrap();
    let comparison = Figure::new(region.clone())
        .title("One motif, three scores")
        .show_region_label(false)
        .label_width(110.0)
        .push(
            LogoTrack::new(0, motif.clone())
                .alphabet_size(4)
                .score(LogoScore::Probability)
                .label("probability")
                .height(60.0),
        )
        .push(
            LogoTrack::new(0, motif.clone())
                .alphabet_size(4)
                .score(LogoScore::InformationContent)
                .label("bits")
                .height(70.0),
        )
        .push(
            LogoTrack::new(0, motif)
                .alphabet_size(4)
                .edlogo()
                .label("enrich / deplete")
                .height(130.0),
        )
        .push(AxisTrack::new().center_on_bases(true));
    comparison.save_svg(out.join("example-logo.svg"))?;

    // Symbols are arbitrary strings, so an alphabet does not have to be four
    // letters wide or one character long.
    let residues = vec![
        LogoColumn::new([("Trp", 62.0), ("Tyr", 25.0), ("Phe", 13.0)]),
        LogoColumn::new([("Gly", 88.0), ("Ala", 12.0)]),
        LogoColumn::new([("Asp", 40.0), ("Glu", 38.0), ("Asn", 12.0), ("Gln", 10.0)]),
        LogoColumn::new([("Cys", 96.0), ("Ser", 4.0)]),
        LogoColumn::new([("Leu", 30.0), ("Ile", 28.0), ("Val", 26.0), ("Met", 16.0)]),
    ];
    let residue_region = Region::new("active site", 0, residues.len() as u64).unwrap();
    let protein = Figure::new(residue_region)
        .title("An arbitrary alphabet: three letter residue codes")
        .show_region_label(false)
        .width(640.0)
        .push(
            LogoTrack::new(0, residues)
                .alphabet_size(20)
                .score(LogoScore::InformationContent)
                .order(StackOrder::LargestOutside)
                .label("residues")
                .height(110.0),
        )
        .push(AxisTrack::new().center_on_bases(true));
    protein.save_svg(out.join("example-logo-protein.svg"))?;

    // The scores are not interchangeable, and the difference is loudest where
    // a symbol is absent. Column 1 has no T and nothing else going on; column 2
    // carries a real gradient. Log odds is dominated by the first, the
    // divergence by the second, and the truth is that both are worth seeing.
    let contrast = vec![
        LogoColumn::acgt(34.0, 33.0, 33.0, 0.0),
        LogoColumn::acgt(50.0, 25.0, 15.0, 10.0),
        LogoColumn::acgt(25.0, 25.0, 25.0, 25.0),
        LogoColumn::acgt(70.0, 12.0, 12.0, 6.0),
    ];
    let contrast_region = Region::new("motif", 0, contrast.len() as u64).unwrap();
    let scores = Figure::new(contrast_region)
        .title("Five ways to score against a background")
        .show_region_label(false)
        .label_width(120.0)
        .width(760.0)
        .push(logo_panel(&contrast, LogoScore::LogOdds, "log odds"))
        .push(logo_panel(
            &contrast,
            LogoScore::KullbackLeibler,
            "KL divergence",
        ))
        .push(logo_panel(&contrast, LogoScore::Difference, "difference"))
        .push(logo_panel(&contrast, LogoScore::Ratio, "ratio"))
        .push(logo_panel(&contrast, LogoScore::OddsRatio, "odds ratio"))
        .push(AxisTrack::new().center_on_bases(true));
    scores.save_svg(out.join("example-logo-scores.svg"))?;

    let (w1, h1) = comparison.dimensions();
    let (w2, h2) = protein.dimensions();
    let (w3, h3) = scores.dimensions();
    println!("example-logo.svg         {w1:.0} x {h1:.0}");
    println!("example-logo-protein.svg {w2:.0} x {h2:.0}");
    println!("example-logo-scores.svg  {w3:.0} x {h3:.0}");
    Ok(())
}

/// One panel of the scoring comparison, all else held equal.
fn logo_panel(columns: &[LogoColumn], score: LogoScore, label: &str) -> LogoTrack {
    LogoTrack::new(0, columns.to_vec())
        .alphabet_size(4)
        .score(score)
        .label(label)
        .height(80.0)
}
