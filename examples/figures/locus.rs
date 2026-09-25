//! The figures `cargo run --example locus` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever window the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::{Aggregate, Drawing, Feature, Plot, Region, Strand, Theme, Variant};

/// The overview window, as a 1-based inclusive locus string.
const WINDOW: &str = "NC_000962.3:761000-762999";
/// Start of the window, 0-based. H37Rv coordinates around rpoB.
const WINDOW_START: u64 = 760_999;
/// Length of the overview window in bases.
const WINDOW_LEN: usize = 2_000;

/// `example.svg`: the rpoB locus and its resistance determining region.
///
/// `theme`, `width` and `region` replace the light theme, the 900 pixels and
/// the two kilobase window the committed figure is drawn with.
pub fn example(theme: &Theme, width: Option<f64>, region: Option<&Region>) -> Box<dyn Drawing> {
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
            .color(theme.color(1)),
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

    let overview = Plot::over(region.cloned().unwrap_or_else(overview_window))
        .title("rpoB locus, resistance determining region")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        // Both from the first base of the window the example was built for,
        // rather than from the left edge of the one on display, which is only
        // the same place until the window moves.
        .add_coverage_at(WINDOW_START, depth)
        .label("depth")
        // Min, not the default Max: at two and a half bases per pixel a dropout
        // is the thing worth not smoothing away.
        .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
        .add_sequence_at(WINDOW_START, bases)
        .label("reference")
        .add_features(genes)
        .label("annotation")
        .add_variants(variants)
        .label("variants")
        .into_figure();
    Box::new(overview)
}

/// `example-zoom.svg`: the same locus at base resolution.
///
/// `theme`, `width` and `region` replace the light theme, the 900 pixels and
/// the sixty bases the committed figure is drawn over.
pub fn example_zoom(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let bases = synthetic_sequence(WINDOW_LEN);
    let depth = synthetic_depth(WINDOW_LEN);

    // The same tracks, zoomed until individual bases are legible. Nothing about
    // the tracks changes; only the region does.
    // Onto the RRDR itself, 60 of its 81 bases, so the S450 codon is in frame.
    let zoom_start = 761_120;
    let zoom_len = 60usize;
    let offset = (zoom_start - WINDOW_START) as usize;
    let own = Region::parse("NC_000962.3:761121-761180").expect("the locus in this example parses");
    let zoom = Plot::over(region.cloned().unwrap_or(own))
        .title("The same locus at base resolution")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .add_coverage_at(zoom_start, depth[offset..offset + zoom_len].to_vec())
        .label("depth")
        .adjust(|track| track.height(45.0))
        .add_sequence_at(zoom_start, bases[offset..offset + zoom_len].to_vec())
        .label("reference")
        .add_variants(vec![
            Variant::new(761_138).value(0.55).category("missense"), // H445Y
            Variant::new(761_154).value(1.00).category("missense"), // S450L
            Variant::new(761_155).value(0.21).category("synonymous"),
        ])
        .label("variants")
        .adjust(|track| track.height(40.0))
        .into_figure();
    Box::new(zoom)
}

/// `example-dark.svg`: the overview in the dark theme.
///
/// The dark theme is what this figure is for, so `theme` is ignored. `width`
/// and `region` replace the 900 pixels and the two kilobase window.
pub fn example_dark(
    _theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let bases = synthetic_sequence(WINDOW_LEN);
    let depth = synthetic_depth(WINDOW_LEN);

    // The dark theme is a selected set of colours rather than an inversion of
    // the light one, so it is worth rendering and looking at.
    let dark = Plot::over(region.cloned().unwrap_or_else(overview_window))
        .title("rpoB locus, resistance determining region")
        .theme(Theme::dark())
        .width(width.unwrap_or(900.0))
        .add_coverage_at(WINDOW_START, depth)
        .label("depth")
        .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
        .add_sequence_at(WINDOW_START, bases)
        .label("reference")
        .add_features(vec![
            Feature::new(759_806, 763_325)
                .name("rpoB")
                .strand(Strand::Forward),
            Feature::new(761_081, 761_162)
                .name("RRDR")
                .strand(Strand::Forward)
                .color(Theme::dark().color(1)),
        ])
        .label("annotation")
        .add_variants(vec![
            Variant::new(761_108).value(0.98).category("missense"),
            Variant::new(761_154).value(1.00).category("missense"),
            Variant::new(761_155).value(0.21).category("synonymous"),
        ])
        .label("variants")
        .into_figure();
    Box::new(dark)
}

/// The window the overview and the dark figure are drawn over.
fn overview_window() -> Region {
    Region::parse(WINDOW).expect("the locus in this example parses")
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
