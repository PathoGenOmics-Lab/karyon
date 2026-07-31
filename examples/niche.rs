//! Renders the three specialised figures used in the README.
//!
//! ```text
//! cargo run --example niche -- assets
//! ```
//!
//! Three plots that belong to one sub-field each and that a general purpose
//! plotting library will never have: raw nanopore current, a gene cluster
//! compared across genomes, and per-strand methylation.

use std::env;
use std::path::{Path, PathBuf};

use karyon::{
    AxisTrack, Feature, Figure, Homology, Locus, LocusTrack, MethylSite, MethylationTrack, Move,
    Region, SquiggleTrack, Strand,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    squiggle(&out)?;
    cluster(&out)?;
    methylation(&out)?;
    Ok(())
}

/// Raw current, and the bases a basecaller made of it.
fn squiggle(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(4_000);
    let bases = b"GATCAGGCTAGCTTGAAACGT";

    // Each base holds the current at a level of its own for a dwell of its own,
    // which is what a nanopore read is: a staircase with noise on it.
    let mut signal: Vec<f64> = Vec::new();
    let mut moves: Vec<Move> = Vec::new();
    for base in bases {
        moves.push(Move::new(signal.len(), *base));
        let level = match base {
            b'A' => 86.0,
            b'C' => 98.0,
            b'G' => 110.0,
            _ => 122.0,
        };
        let dwell = 8 + rng.next() % 14;
        for _ in 0..dwell {
            signal.push(level + (rng.next() % 100) as f64 / 25.0 - 2.0);
        }
    }

    let samples = signal.len() as u64;
    let figure = Figure::new(Region::new("read", 0, samples).unwrap())
        .title("One read, before it was a sequence")
        .width(880.0)
        .show_region_label(false)
        .push(SquiggleTrack::new(0, signal).moves(moves).label("current"))
        .push(AxisTrack::new().label("sample"));

    figure.save_svg(out.join("example-squiggle.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-squiggle.svg {width:.0} x {height:.0}, {samples} samples");
    Ok(())
}

/// The same locus in three genomes, and what each one is missing.
fn cluster(out: &Path) -> std::io::Result<()> {
    // One colour per gene family, so a rearrangement is visible without
    // following a single ribbon.
    let family = [
        "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#7b3294",
    ];
    let gene = |start: u64, len: u64, name: &str, group: usize, forward: bool| {
        Feature::new(start, start + len)
            .name(name)
            .color(family[group])
            .strand(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
    };

    let loci = vec![
        Locus::new(
            "H37Rv",
            vec![
                gene(200, 1_500, "esxB", 0, true),
                gene(1_800, 1_400, "esxA", 1, true),
                gene(3_300, 2_600, "espI", 2, true),
                gene(6_100, 1_900, "eccA1", 3, false),
                gene(8_200, 3_100, "eccB1", 4, true),
            ],
        ),
        Locus::new(
            "CDC1551",
            vec![
                gene(200, 1_500, "esxB", 0, true),
                gene(1_800, 1_400, "esxA", 1, true),
                // espI is gone from this one.
                gene(3_400, 1_900, "eccA1", 3, false),
                gene(5_500, 3_100, "eccB1", 4, true),
            ],
        ),
        Locus::new(
            "BCG",
            vec![
                gene(200, 1_500, "esxB", 0, true),
                // esxA and espI both gone: the RD1 deletion.
                gene(1_900, 1_900, "eccA1", 3, true),
                gene(4_000, 3_100, "eccB1", 4, true),
                gene(7_300, 1_200, "IS6110", 5, false),
            ],
        ),
    ];

    let links = vec![
        Homology::new(0, 0, 0, 0.99),
        Homology::new(0, 1, 1, 0.98),
        Homology::new(0, 3, 2, 0.97),
        Homology::new(0, 4, 3, 0.99),
        Homology::new(1, 0, 0, 0.99),
        Homology::new(1, 2, 1, 0.94),
        Homology::new(1, 3, 2, 0.98),
    ];

    let track = LocusTrack::new(loci).links(links).label("ESX-1");
    let missing = track.unmatched(2).len();

    let figure = Figure::new(Region::new("ESX-1", 0, 11_600).unwrap())
        .title("One locus, three genomes, and what each one is missing")
        .width(880.0)
        .show_region_label(false)
        .push(track)
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-cluster.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-cluster.svg {width:.0} x {height:.0}, {missing} unmatched in BCG");
    Ok(())
}

/// Methylation, one lane per strand.
fn methylation(out: &Path) -> std::io::Result<()> {
    let start = 1_460_000u64;
    let span = 3_000u64;
    let mut rng = Lcg::new(6_000);

    // GATC sites, mostly methylated on both strands, with a stretch that was
    // caught between replication and maintenance and is modified on one only.
    let mut sites = Vec::new();
    let mut pos = start + 40;
    while pos < start + span {
        let hemi = (start + 1_100..start + 1_700).contains(&pos);
        let coverage = 8 + rng.next() % 55;
        let forward = 0.88 + (rng.next() % 12) as f64 / 100.0;
        let reverse = if hemi {
            (rng.next() % 12) as f64 / 100.0
        } else {
            0.86 + (rng.next() % 14) as f64 / 100.0
        };
        sites.push(MethylSite::new(
            pos,
            Strand::Forward,
            forward,
            coverage as u32,
        ));
        sites.push(MethylSite::new(
            pos,
            Strand::Reverse,
            reverse,
            (coverage - 1) as u32,
        ));
        pos += 40 + rng.next() % 130;
    }

    let track = MethylationTrack::new(sites).label("6mA");
    let hemi = track.hemimethylated(0.5).len();

    let figure = Figure::new(Region::new("NC_000962.3", start, start + span).unwrap())
        .title("6mA, one lane per strand, faded by how many reads called it")
        .width(880.0)
        .push(track)
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-methylation.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-methylation.svg {width:.0} x {height:.0}, {hemi} hemimethylated sites");
    Ok(())
}

/// A linear congruential generator, so the figures are reproducible without a
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
