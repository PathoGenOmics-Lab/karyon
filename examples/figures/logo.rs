//! The figures `cargo run --example logo` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever columns the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::{Drawing, LogoColumn, LogoScore, LogoTrack, Plot, Region, StackOrder, Theme};

/// `example-logo.svg`: one motif, three scores.
///
/// The point of the figure is the comparison. The same eight columns are drawn
/// three ways, and only the third one can say anything about the symbols that
/// are missing. `theme`, `width` and `region` replace the light theme, the 900
/// pixels and the eight columns the committed figure is drawn with.
pub fn example_logo(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
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

    let own = columns("motif", motif.len());
    let comparison = Plot::over(region.cloned().unwrap_or(own))
        .title("One motif, three scores")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .remove_region_label()
        .label_width(110.0)
        // From the first column wherever the window is, since the columns
        // are where they are and the window is what moves.
        .add_logo_at(0, motif.clone())
        .label("probability")
        .adjust(|track| {
            track
                .alphabet_size(4)
                .score(LogoScore::Probability)
                .height(60.0)
        })
        .add_logo_at(0, motif.clone())
        .label("bits")
        .adjust(|track| {
            track
                .alphabet_size(4)
                .score(LogoScore::InformationContent)
                .height(70.0)
        })
        .add_logo_at(0, motif)
        .label("enrich / deplete")
        .adjust(|track| track.alphabet_size(4).edlogo().height(130.0))
        .add_axis()
        .adjust(|axis| axis.center_on_bases(true))
        .into_figure();
    Box::new(comparison)
}

/// `example-logo-protein.svg`: an arbitrary alphabet of three letter codes.
///
/// `theme`, `width` and `region` replace the light theme, the 640 pixels and
/// the five columns the committed figure is drawn with.
pub fn example_logo_protein(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // Symbols are arbitrary strings, so an alphabet does not have to be four
    // letters wide or one character long.
    let residues = vec![
        LogoColumn::new([("Trp", 62.0), ("Tyr", 25.0), ("Phe", 13.0)]),
        LogoColumn::new([("Gly", 88.0), ("Ala", 12.0)]),
        LogoColumn::new([("Asp", 40.0), ("Glu", 38.0), ("Asn", 12.0), ("Gln", 10.0)]),
        LogoColumn::new([("Cys", 96.0), ("Ser", 4.0)]),
        LogoColumn::new([("Leu", 30.0), ("Ile", 28.0), ("Val", 26.0), ("Met", 16.0)]),
    ];
    let own = columns("active site", residues.len());
    let protein = Plot::over(region.cloned().unwrap_or(own))
        .title("An arbitrary alphabet: three letter residue codes")
        .theme(theme.clone())
        .remove_region_label()
        .width(width.unwrap_or(640.0))
        .add_logo_at(0, residues)
        .label("residues")
        .adjust(|track| {
            track
                .alphabet_size(20)
                .score(LogoScore::InformationContent)
                .order(StackOrder::LargestOutside)
                .height(110.0)
        })
        .add_axis()
        .adjust(|axis| axis.center_on_bases(true))
        .into_figure();
    Box::new(protein)
}

/// `example-logo-scores.svg`: five ways to score against a background.
///
/// `theme`, `width` and `region` replace the light theme, the 760 pixels and
/// the four columns the committed figure is drawn with.
pub fn example_logo_scores(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
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
    let own = columns("motif", contrast.len());
    let scores = Plot::over(region.cloned().unwrap_or(own))
        .title("Five ways to score against a background")
        .theme(theme.clone())
        .remove_region_label()
        .label_width(120.0)
        .width(width.unwrap_or(760.0))
        .add_track(logo_panel(&contrast, LogoScore::LogOdds, "log odds"))
        .add_track(logo_panel(
            &contrast,
            LogoScore::KullbackLeibler,
            "KL divergence",
        ))
        .add_track(logo_panel(&contrast, LogoScore::Difference, "difference"))
        .add_track(logo_panel(&contrast, LogoScore::Ratio, "ratio"))
        .add_track(logo_panel(&contrast, LogoScore::OddsRatio, "odds ratio"))
        .add_axis()
        .adjust(|axis| axis.center_on_bases(true))
        .into_figure();
    Box::new(scores)
}

/// `example-logo-stability.svg`: the same proportions, five sequences to five
/// hundred.
///
/// `theme`, `width` and `region` replace the light theme, the 700 pixels and
/// the six columns the committed figure is drawn with.
pub fn example_logo_stability(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // How much of a logo should you believe? The proportions below are the
    // same at every sample size, so the raw panels are identical whether they
    // rest on five sequences or five hundred. Shrinkage is what puts the
    // difference back on the page.
    let shape = [
        [5.0, 0.0, 0.0, 0.0],
        [3.0, 1.0, 1.0, 0.0],
        [2.0, 1.0, 1.0, 1.0],
        [4.0, 1.0, 0.0, 0.0],
        [2.0, 2.0, 1.0, 0.0],
        [3.0, 0.0, 2.0, 0.0],
    ];
    let at_depth = |sequences: f64| -> Vec<LogoColumn> {
        let scale = sequences / 5.0;
        shape
            .iter()
            .map(|c| LogoColumn::acgt(c[0] * scale, c[1] * scale, c[2] * scale, c[3] * scale))
            .collect()
    };

    let own = columns("motif", shape.len());
    let mut stability = Plot::over(region.cloned().unwrap_or(own))
        .title("The same proportions, five sequences to five hundred")
        .theme(theme.clone())
        .remove_region_label()
        .label_width(150.0)
        .width(width.unwrap_or(700.0));
    for sequences in [5.0, 50.0, 500.0] {
        stability = stability
            .add_track(depth_panel(&at_depth(sequences), sequences, false))
            .add_track(depth_panel(&at_depth(sequences), sequences, true));
    }
    let stability = stability
        .add_axis()
        .adjust(|axis| axis.center_on_bases(true))
        .into_figure();
    Box::new(stability)
}

/// One column per position, the coordinates a logo is drawn in, named the way
/// the committed figures name them.
fn columns(name: &str, count: usize) -> Region {
    Region::parse(&format!("{name}:1-{count}")).expect("a logo has at least one column")
}

/// One panel of the scoring comparison, all else held equal.
fn logo_panel(columns: &[LogoColumn], score: LogoScore, label: &str) -> LogoTrack {
    LogoTrack::new(0, columns.to_vec())
        .alphabet_size(4)
        .score(score)
        .label(label)
        .height(80.0)
}

/// One panel of the sample size comparison, in bits so the axis is fixed at
/// two and the panels can be read against each other.
fn depth_panel(columns: &[LogoColumn], sequences: f64, stabilized: bool) -> LogoTrack {
    let track = LogoTrack::new(0, columns.to_vec())
        .alphabet_size(4)
        .score(LogoScore::InformationContent)
        .label(if stabilized {
            format!("n = {sequences:.0}, shrunk")
        } else {
            format!("n = {sequences:.0}, raw")
        })
        .height(55.0);
    if stabilized {
        track.stabilize()
    } else {
        track
    }
}
