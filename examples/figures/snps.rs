//! The figure `cargo run --example snps` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page, in the page's colours and
//! at the width of the column. One copy of the figure serves both, so the two
//! cannot drift apart.

use karyon::read;
use karyon::track::traits::Traits;
use karyon::tree::Tree;
use karyon::{Drawing, Figure, MsaSequence, Region, SnpTrack, Theme};

/// Length of the stretch the isolates were compared over.
const LOCUS: usize = 30_000;
/// Where that stretch starts, 0-based.
const START: u64 = 1_472_000;

/// `example-snps.svg`: the variable sites of thirty kilobases, and nothing
/// else.
///
/// `theme` and `width` replace the light theme and the 900 pixels. The panel
/// spaces the sites evenly rather than by where they are, so there is no
/// window along the genome to move, and `region` is ignored.
pub fn example_snps(
    theme: &Theme,
    width: Option<f64>,
    _region: Option<&Region>,
) -> Box<dyn Drawing> {
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

    // Read rather than built, so the example exercises the same path a shell
    // does. The join is names: the reference and the one isolate the tree does
    // not carry are absent from the sheet, and their cells say so.
    let sheet = read::sheet::sheet(SAMPLES).expect("the sheet in this example is well formed");
    let columns = sheet.columns.clone();
    let traits = Traits::new(sheet.rows).spread(columns);

    let panel = SnpTrack::from_alignment(0, &alignment)
        .offset(START)
        .tree(tree)
        .traits(traits)
        .label("isolates")
        .row_height(16.0);

    let sites = panel.sites().len();
    let own = Region::new("sites", 0, sites as u64).unwrap();
    let figure = Figure::new(own)
        .title(format!(
            "{sites} variable sites in {} kb, and nothing else",
            LOCUS / 1_000
        ))
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .show_region_label(false)
        .push(panel);
    Box::new(figure)
}

/// What is known about the isolates, in the shape a sample sheet arrives in.
///
/// Lineage follows the phylogeny, because a lineage is inherited. Resistance
/// does not: it is in one member of the first clade and two of the third, which
/// is what convergence looks like in a strip. ERR400243 is in the panel and not
/// in the sheet, so its cells are drawn as absent rather than as anything.
const SAMPLES: &str = "\
sample\tlineage\tresistance\tyear
ERR400113\tL4\trifampicin\t2016
ERR400152\tL4\tsusceptible\t2017
ERR400191\tL4\tsusceptible\t2015
ERR400230\tL4\tsusceptible\t2019
ERR400100\tL2\tsusceptible\t2014
ERR400126\tL2\tsusceptible\t2018
ERR400139\tL2\tisoniazid\t2020
ERR400165\tL2\tsusceptible\t2016
ERR400178\tL1\trifampicin\t2021
ERR400204\tL1\trifampicin\t2019
ERR400217\tL1\tsusceptible\t
";

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
            for cell in row.iter_mut().skip(start).take(3) {
                *cell = b'-';
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
