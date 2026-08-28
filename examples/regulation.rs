//! Renders the regulation figure used in the documentation.
//!
//! ```text
//! cargo run --example regulation -- assets
//! ```
//!
//! What an RNA-seq experiment and a model of it say about the same stretch:
//! how deep the reads lay, which introns they stepped over and how many took
//! each one, and which bases a sequence model leaned on when it predicted the
//! signal.
//!
//! The three bands answer three different questions and are drawn three
//! different ways on purpose. Depth is a height above a floor. A junction has
//! no height at all, so its arcs are put in lanes and the count is printed on
//! them. And an attribution is a signed number carried by a base, so the base
//! itself is the mark and it hangs below the line where the model pulled the
//! other way.

use std::env;
use std::path::PathBuf;

use karyon::{
    CoverageTrack, DynseqTrack, Junction, JunctionTrack, Motif, Panels, Plot, Region, Strand,
};

/// Length of the window the figure is drawn over.
const SPAN: u64 = 3_600;
/// The three exons, 0-based half-open.
const EXONS: [(u64, u64); 3] = [(0, 520), (1_360, 1_610), (2_760, 3_600)];
/// The motif the model leaned on, sitting in the first exon's promoter.
const MOTIF: &[u8] = b"TTCCTCTTTCCTCC";
/// Where that motif starts.
const MOTIF_AT: u64 = 232;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let region = Region::new("chr17", 0, SPAN).unwrap();
    let mut rng = Lcg::new(9_017);

    // Depth. The middle exon is the skipped one, so it lies about a tenth as
    // deep as its neighbours.
    let depth: Vec<f64> = (0..SPAN)
        .map(
            |pos| match EXONS.iter().position(|(a, b)| (*a..*b).contains(&pos)) {
                Some(1) => 34.0 + (rng.next() % 9) as f64,
                Some(_) => 380.0 + (rng.next() % 70) as f64,
                None => 0.0,
            },
        )
        .collect();

    // The junctions. Two carry the dominant isoform and one skips the middle
    // exon, which is what makes that exon shallow. The last was crossed by
    // nobody, so it is counted and not drawn.
    let junctions = vec![
        Junction::new(520, 1_360, 412)
            .multi(18)
            .motif(Motif::GtAg)
            .strand(Strand::Forward)
            .annotated(true),
        Junction::new(1_610, 2_760, 396)
            .multi(11)
            .motif(Motif::GtAg)
            .strand(Strand::Forward)
            .annotated(true),
        Junction::new(520, 2_760, 31)
            .multi(4)
            .motif(Motif::GtAg)
            .strand(Strand::Forward)
            .annotated(false),
        Junction::new(820, 2_760, 2)
            .motif(Motif::Noncanonical)
            .strand(Strand::Forward)
            .annotated(false),
        Junction::new(520, 3_100, 0).annotated(false),
    ];

    // The model. It leaned on the motif, pulled away from a stretch downstream
    // of it, and was never run over the intron at all, which is why the rule
    // under that stretch is missing rather than flat.
    let filler = b"ACGTTGCAAGCTTAGCCATGGATCCGTAACGGCATTAGC";
    let mut seq: Vec<u8> = Vec::new();
    let mut scores: Vec<f64> = Vec::new();
    for pos in 0..SPAN {
        let inside_motif = (MOTIF_AT..MOTIF_AT + MOTIF.len() as u64).contains(&pos);
        seq.push(if inside_motif {
            MOTIF[(pos - MOTIF_AT) as usize]
        } else {
            filler[(pos as usize) % filler.len()]
        });
        scores.push(if inside_motif {
            0.30 + 0.60 * ((pos - MOTIF_AT) as f64 / 4.0).sin().abs()
        } else if (300..340).contains(&pos) {
            -0.20 - 0.45 * ((pos % 7) as f64) / 7.0
        } else if !EXONS.iter().any(|(a, b)| (*a..*b).contains(&pos)) {
            // Never scored: the model was run over the exons and nowhere else.
            f64::NAN
        } else {
            0.03 * (((pos % 11) as f64) - 5.0)
        });
    }

    let track = JunctionTrack::new(junctions)
        .label("junctions")
        .height(96.0);
    let held = track.discarded();
    let busiest = track.busiest().unwrap_or(0);

    let whole = Plot::over(region)
        .width(760.0)
        .add_track(track)
        .add_track(CoverageTrack::new(0, depth.clone()).label("depth"))
        .add_track(
            DynseqTrack::new(0, seq.clone(), scores.clone())
                .label("attribution")
                .height(96.0),
        )
        .into_figure();

    // The same attribution track over the promoter alone. Nothing about the
    // track changed: below about five pixels a base it draws an envelope of the
    // extremes, and above it the letters, and the zoom is what picks.
    let close = Plot::over(Region::new("chr17", 214, 268).unwrap())
        .width(760.0)
        .add_track(
            DynseqTrack::new(0, seq, scores)
                .label("attribution")
                .height(150.0),
        )
        .add_track(CoverageTrack::new(0, depth).label("depth").height(46.0))
        .into_figure();

    let sheet = Panels::new()
        .title("One promoter, three ways of asking what is happening at it")
        .columns(1)
        .gap(26.0)
        .push_captioned(
            &whole,
            "A",
            "Depth, the introns the reads stepped over, and what the model leaned on",
        )
        .push_captioned(
            &close,
            "B",
            "The same attribution over the motif, where a base is wide enough to be a letter",
        );

    sheet.save_svg(out.join("example-regulation.svg"))?;
    let (width, height) = sheet.dimensions();
    println!("example-regulation.svg {width:.0} x {height:.0}");
    println!("busiest junction {busiest} reads, {held} not drawn");
    Ok(())
}

/// A linear congruential generator, so the figure is reproducible without a
/// dependency.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
}
