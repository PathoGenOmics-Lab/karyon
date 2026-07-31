//! Renders the pangenome accumulation figure used in the README.
//!
//! ```text
//! cargo run --example accumulation -- assets
//! ```
//!
//! Whether the pangenome closes. The total gene repertoire climbs as genomes
//! are added and the core falls, and both depend on the order they were added
//! in, so the order is shuffled two hundred times and the spread across those
//! shuffles is drawn with the curves. Without the spread the figure is one
//! arbitrary ordering presented as a result.

use std::env;
use std::path::PathBuf;

use karyon::{AccumulationTrack, AxisTrack, Figure, Region};

/// How many genomes are in the collection.
const GENOMES: usize = 40;
/// How many gene families every genome carries.
const CORE: usize = 3_100;
/// How many families are shared by some genomes but not all.
const ACCESSORY: usize = 1_400;
/// How many families turn up in a single genome and nowhere else.
const UNIQUE: usize = 40;
/// How many orderings the spread is measured over.
const PERMUTATIONS: usize = 200;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let mut rng = Lcg::new(8_675_309);
    let families = CORE + ACCESSORY + UNIQUE;

    let genomes: Vec<Vec<bool>> = (0..GENOMES)
        .map(|genome| {
            (0..families)
                .map(|family| {
                    if family < CORE {
                        true
                    } else if family < CORE + ACCESSORY {
                        // An accessory family sits in a minority of genomes,
                        // which is what keeps the core falling for a while.
                        rng.next() % 100 < 35
                    } else {
                        // A family of its own for every genome, which is what
                        // keeps the pangenome open.
                        family - (CORE + ACCESSORY) == genome
                    }
                })
                .collect()
        })
        .collect();

    let figure = Figure::new(Region::new("genomes", 0, GENOMES as u64).unwrap())
        .title("Does the pangenome close?")
        .width(880.0)
        .show_region_label(false)
        .push(
            AccumulationTrack::from_presence(&genomes, PERMUTATIONS, 4_242)
                .label("gene families")
                .height(190.0),
        )
        .push(AxisTrack::new().center_on_bases(true).label("genomes"));

    figure.save_svg(out.join("example-accumulation.svg"))?;

    let track = AccumulationTrack::from_presence(&genomes, PERMUTATIONS, 4_242);
    let pan = track.curves()[0].quantile(GENOMES, 0.5).unwrap();
    let core = track.curves()[1].quantile(GENOMES, 0.5).unwrap();
    let (width, height) = figure.dimensions();
    println!(
        "example-accumulation.svg {width:.0} x {height:.0}, \
         {pan:.0} families across {GENOMES} genomes and {core:.0} in all of them"
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
