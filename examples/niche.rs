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
    AxisTrack, Feature, FeatureTrack, Figure, Homology, Legend, LegendTrack, Locus, LocusTrack,
    MethylSite, MethylationTrack, Move, Region, SquiggleTrack, Strand, Theme,
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
    // locus is one quiet slate and the block RD1 removed is one hue: two
    // colours where eight would have made the reader learn a key before
    // reading anything.
    let theme = Theme::light();
    let kept = mix(&theme.muted, theme.surface(), 0.35);
    let lost = theme.color(0).to_string();

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

    // Rv3868 onwards, with the real gene lengths. RD1 removes Rv3871 through
    // Rv3879c, so the block that goes is eccCb1, PE35, PPE68, esxB and esxA and
    // what follows: one contiguous stretch, which is why the figure shows one
    // contiguous hole rather than several scattered ones. The three genomes
    // number the locus differently, as three assemblies of one region do.
    let core = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(300), 1_722, "eccA1", &kept, true),
            gene(at(2_200), 1_455, "eccB1", &kept, true),
            gene(at(3_800), 2_463, "eccCa1", &kept, true),
        ]
    };
    let rd1 = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(6_400), 1_776, "eccCb1", &lost, true),
            gene(at(8_300), 300, "PE35", &lost, true),
            gene(at(8_700), 1_107, "PPE68", &lost, true),
            gene(at(9_900), 303, "esxB", &lost, true),
            gene(at(10_300), 288, "esxA", &lost, true),
        ]
    };

    let mut h37rv = core(0);
    h37rv.extend(rd1(0));
    let mut bovis = core(-350);
    bovis.extend(rd1(-350));

    // BCG keeps the backbone and nothing else: RD1 left a clean junction, not
    // a scar with something in it.
    let bcg = core(150);

    let loci = vec![
        Locus::new("H37Rv", h37rv),
        Locus::new("M. bovis", bovis),
        Locus::new("BCG", bcg),
    ];

    let links = vec![
        Homology::new(0, 0, 0, 0.999),
        Homology::new(0, 1, 1, 0.997),
        Homology::new(0, 2, 2, 0.996),
        Homology::new(0, 3, 3, 0.99),
        Homology::new(0, 4, 4, 0.95),
        // PPE68 is the one that has drifted between the two.
        Homology::new(0, 5, 5, 0.76),
        Homology::new(0, 6, 6, 0.999),
        Homology::new(0, 7, 7, 0.998),
        // BCG matches only the backbone, so only three ribbons leave this row.
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
        .area("deleted in BCG", &lost)
        .ramp("identity", pale, dark, "70%", "100%")
        .outline("in no neighbouring locus", theme.foreground.clone());

    let figure = Figure::new(Region::new("ESX-1", 0, 13_000).unwrap())
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
    // Escherichia coli K-12 MG1655 at oriC. Dam writes 6mA at GATC on both
    // strands, and a few seconds behind the fork the new strand has not been
    // methylated yet, so every site is briefly hemimethylated. At oriC SeqA
    // holds that state for about a third of the cell cycle, which is the
    // sequestration that stops the origin firing twice, and it is the reason
    // the two strands need their own lanes rather than an average.
    //
    // This is the assay's own organism. A genome with no dam has no GATC
    // methylation to plot.
    let start = 3_924_500u64;
    let span = 3_000u64;
    // oriC itself, 1-based 3,925,744 to 3,925,975.
    let oric = (3_925_743u64, 3_925_975u64);
    let mut rng = Lcg::new(6_000);

    // The eleven GATC sites inside oriC, 0-based at the G. Hard-coded rather
    // than generated: eleven is the textbook count and their spacing is the
    // point.
    let inside: [u64; 11] = [
        3_925_743, 3_925_759, 3_925_774, 3_925_787, 3_925_813, 3_925_827, 3_925_845, 3_925_862,
        3_925_869, 3_925_891, 3_925_959,
    ];

    let mut positions: Vec<u64> = inside.to_vec();
    // Genome-wide a GATC turns up about every 240 bases, so the flanks are
    // sparser than the origin by design.
    let mut pos = start + 60;
    while pos < start + span {
        if pos < oric.0 || pos >= oric.1 {
            positions.push(pos);
        }
        pos += 120 + rng.next() % 260;
    }
    positions.sort_unstable();

    let mut sites = Vec::new();
    for pos in positions {
        let sequestered = pos >= oric.0 && pos < oric.1;
        let coverage = 8 + rng.next() % 55;
        // GATC is a palindrome, so the two adenines are one base apart: the
        // forward one at G+1 and the reverse one at G+2.
        let forward = 0.88 + (rng.next() % 12) as f64 / 100.0;
        let reverse = if sequestered {
            (rng.next() % 12) as f64 / 100.0
        } else {
            0.86 + (rng.next() % 14) as f64 / 100.0
        };
        sites.push(MethylSite::new(
            pos + 1,
            Strand::Forward,
            forward,
            coverage as u32,
        ));
        sites.push(MethylSite::new(
            pos + 2,
            Strand::Reverse,
            reverse,
            (coverage - 1) as u32,
        ));
    }

    let track = MethylationTrack::new(sites).label("6mA at GATC");
    let hemi = track.hemimethylated(0.5).len();

    let figure = Figure::new(Region::new("NC_000913.3", start, start + span).unwrap())
        .title("Dam methylation across oriC, one lane per strand")
        .width(880.0)
        .push(track)
        .push(
            FeatureTrack::new(vec![Feature::new(oric.0, oric.1).name("oriC")])
                .label("origin")
                .row_height(14.0),
        )
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
