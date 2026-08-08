//! Renders the four figures whose shape belongs to genomics and nowhere else.
//!
//! ```text
//! cargo run --example shapes -- assets
//! ```
//!
//! Structural variants as arcs between their breakpoints, the six reading
//! frames with their stops, two trees face to face, and methylation one
//! molecule at a time. Four organisms, because the forms are not about any of
//! them.

use std::env;
use std::path::{Path, PathBuf};

use karyon::tree::Tree;
use karyon::{
    plot, BisulfiteTrack, Legend, Molecule, OrfTrack, Plot, Region, StructuralTrack,
    StructuralVariant, SvKind, TanglegramTrack, Theme,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    structural(&out)?;
    frames(&out)?;
    tanglegram(&out)?;
    bisulfite(&out)?;
    Ok(())
}

/// Arcs between breakpoints, over the depth that should agree with them.
fn structural(out: &Path) -> std::io::Result<()> {
    // The RD1 neighbourhood of M. tuberculosis H37Rv, so the named call is where
    // the call actually is. Drawing RD1 at 8 kb would put it on top of gyrA and
    // gyrB, which is a different and much less survivable event.
    let from = 4_340_000u64;
    let span = 60_000u64;
    let mut rng = Lcg::new(7_919);

    let calls = vec![
        // RD1: the 9.5 kb block missing from every BCG, Rv3871 to Rv3879c,
        // carrying esxB and esxA.
        StructuralVariant::new(4_350_000, 4_359_600, SvKind::Deletion)
            .support(34)
            .name("RD1"),
        StructuralVariant::new(4_366_000, 4_373_000, SvKind::Duplication).support(21),
        StructuralVariant::new(4_381_000, 4_387_000, SvKind::Inversion).support(9),
        // One circular chromosome and no plasmid, so the far breakpoint of a
        // rearrangement is elsewhere on the same molecule and off this window.
        // The arc leaving the frame is the honest picture of that.
        StructuralVariant::new(4_345_000, 3_120_000, SvKind::Translocation).support(6),
        StructuralVariant::new(4_362_000, 4_362_000, SvKind::Insertion).support(17),
    ];

    // Depth that agrees: a hole under the deletion, a step under the
    // duplication, and nothing under the inversion, which is how an inversion
    // looks and why the arc has to say what the depth cannot.
    let depth: Vec<f64> = (from..from + span)
        .map(|at| {
            let base = if (4_350_000..4_359_600).contains(&at) {
                2.0
            } else if (4_366_000..4_373_000).contains(&at) {
                84.0
            } else {
                42.0
            };
            base + (rng.next() % 70) as f64 / 10.0
        })
        .collect();

    let track = StructuralTrack::new(calls).label("SV").height(96.0);
    let theme = Theme::light();
    let legend = Legend::new()
        .line("deletion", track.color_of(SvKind::Deletion, &theme))
        .line("duplication", track.color_of(SvKind::Duplication, &theme))
        .line("inversion", track.color_of(SvKind::Inversion, &theme))
        .line(
            "translocation",
            track.color_of(SvKind::Translocation, &theme),
        )
        .line("insertion", track.color_of(SvKind::Insertion, &theme));

    let figure = Plot::over(Region::new("NC_000962.3", from, from + span).unwrap())
        .title("Structural variants, and whether the depth agrees")
        .width(880.0)
        .remove_region_label()
        .add_track(track)
        .add_coverage(depth)
        .label("depth")
        .adjust(|track| track.height(56.0))
        .add_legend(legend)
        .into_figure();

    figure.save_svg(out.join("example-structural.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-structural.svg {width:.0} x {height:.0}");
    Ok(())
}

/// Six lanes, their stops, and what is open between them.
fn frames(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(2_718);
    let span = 3_600usize;

    // Random sequence, then a real gene written into frame +1: random codons
    // drawn from the sixty-one that are not stops. Writing one repeated codon
    // instead would leave the other two forward frames open as well, and the
    // figure would show three genes where the biology has one.
    let mut seq: Vec<u8> = (0..span)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    let gene_at = 900usize;
    seq[gene_at..gene_at + 3].copy_from_slice(b"ATG");
    let mut at = gene_at + 3;
    while at + 3 <= gene_at + 1_500 {
        let codon = loop {
            let picked = [
                b"ACGT"[(rng.next() % 4) as usize],
                b"ACGT"[(rng.next() % 4) as usize],
                b"ACGT"[(rng.next() % 4) as usize],
            ];
            if !matches!(
                picked,
                [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
            ) {
                break picked;
            }
        };
        seq[at..at + 3].copy_from_slice(&codon);
        at += 3;
    }
    seq[at..at + 3].copy_from_slice(b"TAA");

    let track = OrfTrack::new(0, seq).min_codons(60).label("frames");
    let orfs = track.orfs().len();

    let figure = Plot::over(Region::new("plasmid", 0, span as u64).unwrap())
        .title("Six reading frames: the stops, and what is open between them")
        .width(880.0)
        .remove_region_label()
        .add_track(track.lane_height(11.0))
        .into_figure();

    figure.save_svg(out.join("example-frames.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-frames.svg {width:.0} x {height:.0}, {orfs} open frames");
    Ok(())
}

/// Two trees, and where they disagree.
fn tanglegram(out: &Path) -> std::io::Result<()> {
    // Klebsiella pneumoniae, which has a large accessory genome and moves it
    // around. The species matters here: in a clonal organism with a closed
    // pangenome the two trees would agree by construction and the figure would
    // have nothing to show.
    let core = Tree::parse_annotated_newick(
        "(((ERR5001[&ward=ICU]:0.004,ERR5002[&ward=ICU]:0.004):0.010,\
          (ERR5003[&ward=Surgery]:0.003,ERR5004[&ward=Surgery]:0.005):0.009):0.020,\
          ((ERR5005[&ward=Medicine]:0.004,ERR5006[&ward=Medicine]:0.003):0.011,\
          (ERR5007[&ward=Community]:0.005,ERR5008[&ward=Community]:0.004):0.008):0.018);",
    )
    .expect("the tree in this example is well formed");
    // The accessory genome tells a different story: two isolates have swapped
    // sides, which is what horizontal transfer looks like in a tanglegram.
    let accessory = Tree::parse_annotated_newick(
        "(((ERR5001[&ward=ICU]:0.030,ERR5006[&ward=Medicine]:0.028):0.040,\
          (ERR5003[&ward=Surgery]:0.031,ERR5004[&ward=Surgery]:0.029):0.038):0.050,\
          ((ERR5005[&ward=Medicine]:0.030,ERR5002[&ward=ICU]:0.032):0.041,\
          (ERR5007[&ward=Community]:0.028,ERR5008[&ward=Community]:0.030):0.039):0.048);",
    )
    .expect("the tree in this example is well formed");

    let track = TanglegramTrack::new(core, accessory)
        .names("core genome", "accessory genome")
        .color_by("ward")
        .untangle()
        .tie_widths(1.0, 1.7)
        .row_height(20.0)
        .label("8 isolates");
    let crossings = track.crossings();

    let figure = plot("taxa:1-8")?
        .title("K. pneumoniae: core and accessory genome trees over one collection")
        .width(760.0)
        .remove_region_label()
        .remove_axis()
        .add_track(track)
        .into_figure();

    figure.save_svg(out.join("example-tanglegram.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-tanglegram.svg {width:.0} x {height:.0}, {crossings} crossings");
    Ok(())
}

/// One row per molecule, one column per cytosine.
fn bisulfite(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(1_729);
    // The H19/IGF2 imprinting control region on human chromosome 11. This is
    // the assay's own ground: bisulfite conversion reads 5-methylcytosine, and
    // the ICR is methylated on the paternal allele and unmethylated on the
    // maternal one, so a single sample really does carry two populations of
    // molecule rather than one loosely modified set.
    let start = 2_002_400u64;

    // CpGs across the ICR, clustered around the CTCF sites and sparser after.
    let sites: Vec<u64> = (0..14)
        .map(|index| {
            if index < 9 {
                start + 40 + index * 22
            } else {
                start + 260 + (index - 9) * 40
            }
        })
        .collect();

    // Two populations, and the reason for them is imprinting rather than
    // chance: the paternal allele is methylated across the ICR and the maternal
    // one is not. The site fractions alone would sit near a half everywhere and
    // could not tell this apart from every molecule being half done.
    let molecules: Vec<Molecule> = (0..16)
        .map(|index| {
            let paternal = index % 2 == 0;
            let calls: Vec<Option<bool>> = sites
                .iter()
                .enumerate()
                .map(|(site, _)| {
                    // A read does not always reach the far sites.
                    if site > 10 && rng.next() % 100 < 30 {
                        return None;
                    }
                    Some(if paternal {
                        rng.next() % 100 > 8
                    } else {
                        rng.next() % 100 < 6
                    })
                })
                .collect();
            Molecule::new(format!("read_{:02}", index + 1), calls)
        })
        .collect();

    let track = BisulfiteTrack::new(sites, molecules).label("CpG");
    let discordance = track.discordance().unwrap_or(0.0);

    let figure = Plot::over(Region::new("NC_000011.10", start, start + 500).unwrap())
        .title("The H19/IGF2 ICR one molecule at a time: two alleles, not one loose pattern")
        .width(880.0)
        .add_track(track.row_height(13.0))
        .into_figure();

    figure.save_svg(out.join("example-bisulfite.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-bisulfite.svg {width:.0} x {height:.0}, discordance {discordance:.2}");
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
