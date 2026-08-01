//! Renders the two figures used in the README.
//!
//! Run it with the output directory as an optional argument:
//!
//! ```text
//! cargo run --example locus -- assets
//! ```
//!
//! The depth and the reference bases are synthetic but shaped like the real
//! thing: a coverage dropout, and a GC content near the 65% this genome sits at.
//! The annotation and the variant positions are not synthetic. `rpoB` spans the
//! whole window and runs off both edges, the RRDR box is the real 81 base
//! hotspot, and every variant sits on the codon it is named for, which is worth
//! the trouble because a figure with an invented gene in it teaches the reader
//! the wrong thing about the tool as well as about the locus.
//!
//! Everything generated is generated from a fixed seed, so re-running the
//! example produces byte-identical files and a diff only appears when the
//! rendering actually changed.

use std::env;
use std::path::PathBuf;

use karyon::{
    Aggregate, AxisTrack, CoverageTrack, Feature, FeatureTrack, Figure, Region, SequenceTrack,
    Strand, Theme, Variant, VariantTrack,
};

/// Start of the window, 0-based. H37Rv coordinates around rpoB.
const WINDOW_START: u64 = 760_999;
/// Length of the overview window in bases.
const WINDOW_LEN: usize = 2_000;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let bases = synthetic_sequence(WINDOW_LEN);
    let depth = synthetic_depth(WINDOW_LEN);

    // rpoB is Rv0667, 1-based 759,807 to 763,325, forward, 3,519 bp. It is far
    // wider than this window and draws with neither end in frame, which is what
    // the annotation actually looks like here. There is no second gene: the next
    // one along is rpoC, 371 bases past the right edge.
    let genes = vec![
        Feature::new(759_806, 763_325)
            .name("rpoB")
            .strand(Strand::Forward),
        // The rifampicin resistance determining region: codons 426 to 452, the
        // 81 base hotspot. CodonTrack pins the same frame, so this box and
        // `CodonTrack::span_of(450)` agree to the base.
        Feature::new(761_081, 761_162)
            .name("RRDR")
            .strand(Strand::Forward)
            .color("#d55e00"),
    ];

    // Codon N of rpoB spans 0-based 759,806 + 3(N-1) to +3, so each of these
    // lands on the base its name says it does.
    let variants = vec![
        Variant::new(761_108).value(0.98).category("missense"), // D435V
        Variant::new(761_138).value(0.55).category("missense"), // H445Y
        Variant::new(761_154).value(1.00).category("missense"), // S450L
        Variant::new(761_155).value(0.21).category("synonymous"), // S450 wobble
        Variant::new(761_051).value(0.12).category("synonymous"), // L426 wobble
    ];

    let overview = Figure::new(Region::new("NC_000962.3", WINDOW_START, 762_999).unwrap())
        .title("rpoB locus, resistance determining region")
        .push(
            // Min, not the default Max: at two and a half bases per pixel a
            // dropout is the thing worth not smoothing away.
            CoverageTrack::new(WINDOW_START, depth.clone())
                .label("depth")
                .aggregate(Aggregate::Min)
                .height(70.0),
        )
        .push(SequenceTrack::new(WINDOW_START, bases.clone()).label("reference"))
        .push(FeatureTrack::new(genes).label("annotation"))
        .push(VariantTrack::new(variants).label("variants"))
        .push(AxisTrack::new());
    overview.save_svg(out.join("example.svg"))?;

    // The same tracks, zoomed until individual bases are legible. Nothing about
    // the tracks changes; only the region does.
    // Onto the RRDR itself, 60 of its 81 bases, so the S450 codon is in frame.
    let zoom_start = 761_120;
    let zoom_len = 60usize;
    let offset = (zoom_start - WINDOW_START) as usize;
    let zoom_region = Region::new("NC_000962.3", zoom_start, zoom_start + zoom_len as u64).unwrap();
    let zoom = Figure::new(zoom_region)
        .title("The same locus at base resolution")
        .push(
            CoverageTrack::new(zoom_start, depth[offset..offset + zoom_len].to_vec())
                .label("depth")
                .height(45.0),
        )
        .push(
            SequenceTrack::new(zoom_start, bases[offset..offset + zoom_len].to_vec())
                .label("reference"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(761_138).value(0.55).category("missense"), // H445Y
                Variant::new(761_154).value(1.00).category("missense"), // S450L
                Variant::new(761_155).value(0.21).category("synonymous"),
            ])
            .label("variants")
            .height(40.0),
        )
        .push(AxisTrack::new());
    zoom.save_svg(out.join("example-zoom.svg"))?;

    // The dark theme is a selected set of colours rather than an inversion of
    // the light one, so it is worth rendering and looking at.
    let dark = Figure::new(Region::new("NC_000962.3", WINDOW_START, 762_999).unwrap())
        .title("The same locus, dark theme")
        .theme(Theme::dark())
        .push(
            CoverageTrack::new(WINDOW_START, depth)
                .label("depth")
                .aggregate(Aggregate::Min)
                .height(70.0),
        )
        .push(SequenceTrack::new(WINDOW_START, bases).label("reference"))
        .push(
            FeatureTrack::new(vec![
                Feature::new(759_806, 763_325)
                    .name("rpoB")
                    .strand(Strand::Forward),
                Feature::new(761_081, 761_162)
                    .name("RRDR")
                    .strand(Strand::Forward)
                    .color("#d55e00"),
            ])
            .label("annotation"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(761_108).value(0.98).category("missense"),
                Variant::new(761_154).value(1.00).category("missense"),
                Variant::new(761_155).value(0.21).category("synonymous"),
            ])
            .label("variants"),
        )
        .push(AxisTrack::new());
    dark.save_svg(out.join("example-dark.svg"))?;

    let (w1, h1) = overview.dimensions();
    let (w2, h2) = zoom.dimensions();
    let (w3, h3) = dark.dimensions();
    println!("example.svg      {w1:.0} x {h1:.0}");
    println!("example-zoom.svg {w2:.0} x {h2:.0}");
    println!("example-dark.svg {w3:.0} x {h3:.0}");
    Ok(())
}

/// A reproducible string at roughly the GC content of this genome, 65%.
///
/// Real sequence would come from a FASTA reader. An equiprobable ACGT draw would
/// be 50% GC, which in a coloured sequence track is visibly the wrong organism.
fn synthetic_sequence(len: usize) -> Vec<u8> {
    let mut rng = Lcg::new(20_260_731);
    (0..len)
        .map(|_| match rng.next() % 1_000 {
            r if r < 328 => b'G',
            r if r < 656 => b'C',
            r if r < 828 => b'A',
            _ => b'T',
        })
        .collect()
}

/// A depth profile with a dropout, the shape a deletion or a GC-poor stretch
/// leaves behind.
fn synthetic_depth(len: usize) -> Vec<f64> {
    let mut rng = Lcg::new(7);
    (0..len)
        .map(|i| {
            let wave = 58.0 + 9.0 * ((i as f64) / 130.0).sin();
            let noise = (rng.next() % 11) as f64 - 5.0;
            let dropout = if (900..1_030).contains(&i) { 0.06 } else { 1.0 };
            ((wave + noise) * dropout).max(0.0)
        })
        .collect()
}

/// A linear congruential generator, so the example needs no dependency to be
/// reproducible.
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
