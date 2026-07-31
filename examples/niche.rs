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

use karyon::theme::mix;
use karyon::{
    AxisTrack, Feature, Figure, Homology, Legend, LegendTrack, Locus, LocusTrack, MethylSite,
    MethylationTrack, Move, Region, SquiggleTrack, Strand, Theme,
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
    // Colour is spent only on what differs. The conserved backbone of the
    // locus is one quiet slate, the block RD1 removed is one hue and the
    // element that replaced it is another: three colours where eight would
    // have made the reader learn a key before reading anything.
    let theme = Theme::light();
    let kept = mix(theme.surface(), &theme.rule, 0.85);
    let lost = theme.color(0).to_string();
    let gained = theme.color(1).to_string();

    let gene = |start: u64, len: u64, name: &str, color: &str, forward: bool| {
        Feature::new(start, start.saturating_add(len))
            .name(name)
            .color(color)
            .strand(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
    };

    // Rv3868 onwards. RD1 takes out Rv3872 to Rv3875, which is PE35, PPE68,
    // esxB and esxA: one contiguous block, which is why the figure shows one
    // contiguous hole rather than four scattered ones. The three genomes
    // number the locus differently, as three assemblies of one region do.
    let core = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(300), 1_800, "eccA1", &kept, true),
            gene(at(2_400), 1_500, "eccB1", &kept, true),
            gene(at(4_100), 3_500, "eccCa1", &kept, true),
        ]
    };
    let rd1 = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(7_800), 300, "PE35", &lost, true),
            gene(at(8_300), 1_100, "PPE68", &lost, true),
            gene(at(9_600), 300, "esxB", &lost, true),
            gene(at(10_100), 300, "esxA", &lost, true),
        ]
    };

    let mut h37rv = core(0);
    h37rv.extend(rd1(0));
    let mut bovis = core(-350);
    bovis.extend(rd1(-350));

    let mut bcg = core(150);
    // The scar RD1 left, with an insertion sequence sitting in it.
    bcg.push(gene(8_050, 1_350, "IS6110", &gained, false));

    let loci = vec![
        Locus::new("H37Rv", h37rv),
        Locus::new("M. bovis", bovis),
        Locus::new("BCG", bcg),
    ];

    let links = vec![
        Homology::new(0, 0, 0, 0.999),
        Homology::new(0, 1, 1, 0.997),
        Homology::new(0, 2, 2, 0.996),
        Homology::new(0, 3, 3, 0.95),
        // PPE68 is the one that has drifted between the two.
        Homology::new(0, 4, 4, 0.76),
        Homology::new(0, 5, 5, 0.999),
        Homology::new(0, 6, 6, 0.998),
        Homology::new(1, 0, 0, 0.998),
        Homology::new(1, 1, 1, 0.999),
        Homology::new(1, 2, 2, 0.997),
    ];

    let track = LocusTrack::new(loci).links(links).label("ESX-1");
    // Genes of the middle row that nothing in BCG matches: the deletion.
    let removed = (0..track.loci()[1].genes.len())
        .filter(|index| {
            !track
                .homologies()
                .iter()
                .any(|link| link.row == 1 && link.from == *index)
        })
        .count();

    // The ramp comes off the track rather than being named again here, so the
    // key cannot drift away from the ribbons it explains.
    let (pale, dark) = track.ramp_ends(&theme);
    let legend = Legend::new()
        .key("deleted in BCG", &lost)
        .key("insertion sequence", &gained)
        .ramp("identity", pale, dark, "70%", "100%")
        .outline("in no neighbouring locus", theme.foreground.clone());

    let figure = Figure::new(Region::new("ESX-1", 0, 12_000).unwrap())
        .title("ESX-1, and the deletion that made BCG a vaccine")
        .width(880.0)
        .show_region_label(false)
        .push(track)
        .push(LegendTrack::new(legend))
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-cluster.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-cluster.svg {width:.0} x {height:.0}, RD1 removed {removed} genes");
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
