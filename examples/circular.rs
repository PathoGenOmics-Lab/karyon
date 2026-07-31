//! Renders the circular genome figure used in the README.
//!
//! ```text
//! cargo run --example circular -- assets
//! ```
//!
//! A bacterial chromosome the way it actually is. Drawn as a line, a circular
//! genome gets an edge the biology does not have, straight through whatever
//! happens to sit at coordinate zero. Drawn as a circle it does not, and the
//! middle becomes usable: the chords join the two ends of the rearrangements.

use std::env;
use std::path::PathBuf;

use karyon::{
    AxisRing, Feature, FeatureRing, MarkerRing, Rings, SignalRing, Strand, Window, WindowTrack,
};

/// Length of the H37Rv chromosome.
const LENGTH: u64 = 4_411_532;
/// Size of one composition window.
const STEP: u64 = 10_000;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let mut rng = Lcg::new(1_906);

    // A sequence whose skew turns over at the origin and again at the terminus,
    // which is what puts the two halves of the ring on opposite sides.
    let bases: Vec<u8> = (0..LENGTH)
        .map(|index| {
            let leading = index < LENGTH / 2;
            match rng.next() % 100 {
                0..=32 if leading => b'G',
                0..=32 => b'C',
                33..=64 if leading => b'C',
                33..=64 => b'G',
                65..=82 => b'A',
                _ => b'T',
            }
        })
        .collect();

    let skew: Vec<Window> = WindowTrack::gc_skew(0, &bases, STEP).windows().to_vec();
    let gc = WindowTrack::gc_content(0, &bases, STEP);
    let content: Vec<Window> = gc.windows().to_vec();

    // A scattering of genes on both strands, and a handful of named loci.
    let mut genes: Vec<Feature> = (0..900)
        .map(|index| {
            let start = index * (LENGTH / 900) + rng.next() % 900;
            Feature::new(start, start + 800 + rng.next() % 1_800).strand(if rng.next() % 2 == 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            })
        })
        .collect();
    genes.push(
        Feature::new(759_807, 763_325)
            .name("rpoB")
            .strand(Strand::Forward),
    );
    genes.push(
        Feature::new(2_153_889, 2_156_111)
            .name("katG")
            .strand(Strand::Reverse),
    );
    genes.push(
        Feature::new(1_673_440, 1_674_183)
            .name("inhA")
            .strand(Strand::Forward),
    );

    // Resistance mutations, coloured by the drug they belong to.
    let mutations: Vec<(u64, usize)> = vec![
        (761_155, 0),
        (761_110, 0),
        (761_140, 0),
        (2_155_168, 1),
        (2_154_724, 1),
        (1_673_425, 2),
        (4_247_431, 3),
        (7_570, 4),
        (1_472_359, 5),
    ];

    let plot = Rings::new(LENGTH)
        .diameter(660.0)
        .title("H37Rv")
        .subtitle("4.41 Mb")
        .push(AxisRing::new())
        .push(FeatureRing::new(genes).thickness(20.0).show_names(false))
        .push(
            FeatureRing::new(vec![
                Feature::new(759_807, 763_325).name("rpoB"),
                Feature::new(2_153_889, 2_156_111).name("katG"),
                Feature::new(1_673_440, 1_674_183).name("inhA"),
            ])
            .thickness(2.0)
            .show_names(true)
            .split_strands(false),
        )
        .push(
            MarkerRing::categorised(mutations)
                .thickness(12.0)
                .width(1.8),
        )
        .push(
            SignalRing::new(content)
                .baseline(gc.baseline_value())
                .thickness(44.0)
                .colors("#7b3294", "#e69f00"),
        )
        .push(
            SignalRing::new(skew)
                .thickness(50.0)
                .colors("#0072b2", "#d55e00"),
        )
        // The two ends of an inversion and of a duplication, which is the one
        // thing the middle of a circle is good for.
        .link_colored((1_100_000, 1_180_000), (3_240_000, 3_320_000), None, 0.3)
        .link_colored(
            (420_000, 470_000),
            (2_760_000, 2_810_000),
            Some("#d55e00".to_string()),
            0.3,
        );

    plot.save_svg(out.join("example-circular.svg"))?;
    let (width, height) = plot.dimensions();
    println!(
        "example-circular.svg {width:.0} x {height:.0}, {} rings",
        plot.ring_count()
    );
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
