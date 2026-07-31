//! Renders the variable site panel used in the README.
//!
//! ```text
//! cargo run --example snps -- assets
//! ```
//!
//! Twelve isolates against a reference over thirty kilobases. Almost all of it
//! agrees, so drawing all thirty thousand columns would spend every pixel on
//! the agreement. The panel keeps only the columns that vary and spaces them
//! evenly, which is why each one carries its own position underneath: the
//! spacing deliberately says nothing about distance.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{Figure, MsaSequence, Region, SnpTrack};

/// Length of the stretch the isolates were compared over.
const LOCUS: usize = 30_000;
/// Where that stretch starts, 0-based.
const START: u64 = 1_472_000;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let alignment = synthetic_alignment();

    // A phylogeny of the same isolates. Passing it to the panel sorts the rows
    // by descent, which is what turns a clade's shared substitutions from
    // scattered noise into a block.
    let tree = Tree::parse_newick(
        "(((ERR400113:0.0002,ERR400152:0.0002)0.99:0.0011,\
          (ERR400191:0.0003,ERR400230:0.0002)0.97:0.0009)0.88:0.0021,\
          ((ERR400100:0.0015,ERR400126:0.0014):0.0006,\
           (ERR400139:0.0012,ERR400165:0.0013):0.0007):0.0018,\
          (ERR400178:0.0016,(ERR400204:0.0011,ERR400217:0.0012):0.0005):0.0014);",
    )
    .expect("the tree in this example is well formed");

    let panel = SnpTrack::from_alignment(0, &alignment)
        .offset(START)
        .tree(tree)
        .label("isolates")
        .row_height(16.0);

    let sites = panel.sites().len();
    let region = Region::new("sites", 0, sites as u64).unwrap();
    let figure = Figure::new(region)
        .title(format!(
            "{sites} variable sites in {} kb, and nothing else",
            LOCUS / 1_000
        ))
        .width(900.0)
        .show_region_label(false)
        .push(panel);

    figure.save_svg(out.join("example-snps.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-snps.svg {width:.0} x {height:.0}");
    println!(
        "{sites} sites kept out of {LOCUS} columns, {:.3}% of them",
        sites as f64 / LOCUS as f64 * 100.0
    );
    Ok(())
}

/// A reference and twelve isolates, two of which form a clade.
fn synthetic_alignment() -> Vec<MsaSequence> {
    let mut rng = Lcg::new(1_848);
    let reference: Vec<u8> = (0..LOCUS)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();

    // Sites the whole clade shares, and sites that belong to one isolate.
    let clade_sites: Vec<usize> = (0..5).map(|_| (rng.next() as usize) % LOCUS).collect();

    let mut rows = vec![MsaSequence::new("H37Rv", reference.clone())];
    for isolate in 0..12 {
        let mut row = reference.clone();
        let in_clade = isolate % 4 == 1;
        if in_clade {
            for site in &clade_sites {
                row[*site] = mutate(reference[*site], &mut rng);
            }
        }
        for _ in 0..(1 + rng.next() % 3) {
            let site = (rng.next() as usize) % LOCUS;
            row[site] = mutate(reference[site], &mut rng);
        }
        // One isolate carries a small deletion.
        if isolate == 7 {
            let start = (rng.next() as usize) % (LOCUS - 3);
            for column in start..start + 3 {
                row[column] = b'-';
            }
        }
        rows.push(MsaSequence::new(
            format!("ERR{:06}", 400_100 + isolate * 13),
            row,
        ));
    }
    rows
}

/// A base that is not the one it replaces.
fn mutate(from: u8, rng: &mut Lcg) -> u8 {
    let options: Vec<u8> = b"ACGT".iter().copied().filter(|b| *b != from).collect();
    options[(rng.next() as usize) % options.len()]
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
