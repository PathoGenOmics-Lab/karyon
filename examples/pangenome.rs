//! Renders the pangenome figure used in the README.
//!
//! ```text
//! cargo run --example pangenome -- assets
//! ```
//!
//! A presence and absence matrix of accessory genes, with the isolates sorted
//! by the phylogeny beside them. Sorted by sample name the same data is a
//! speckle; sorted by descent the accessory regions become rectangles, and a
//! rectangle is a claim about the biology.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{
    AxisTrack, CellScale, Feature, FeatureTrack, Figure, MatrixRow, MatrixTrack, Region, Strand,
};

/// Start of the window, 0-based.
const START: u64 = 2_100_000;
/// Length of the window in bases.
const WINDOW: u64 = 60_000;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    // Two clades, and an accessory island each of them has or lacks.
    let tree = Tree::parse_newick(
        "(((ERR5001:0.0004,ERR5002:0.0004)0.99:0.0016,\
           (ERR5003:0.0005,ERR5004:0.0003)0.96:0.0014)0.91:0.0030,\
          ((ERR5005:0.0011,ERR5006:0.0010)0.93:0.0008,\
           (ERR5007:0.0009,(ERR5008:0.0006,ERR5009:0.0007)0.88:0.0004)0.90:0.0011):0.0026);",
    )
    .expect("the tree in this example is well formed");

    let mut rng = Lcg::new(52_814);
    let genes: Vec<u64> = (0..26).map(|i| START + 900 + i * 2_200).collect();

    // The first clade carries the island in the middle of the window; the
    // second carries a different one near the end. Everything else is core.
    let rows: Vec<MatrixRow> = (1..=9)
        .map(|isolate| {
            let first_clade = isolate <= 4;
            let values: Vec<f64> = genes
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    let core = !(6..=11).contains(&index) && !(18..=22).contains(&index);
                    if core {
                        // Core genes are present, bar the occasional dropout.
                        if rng.next() % 100 < 3 {
                            0.0
                        } else {
                            1.0
                        }
                    } else if rng.next() % 100 < 4 {
                        f64::NAN
                    } else {
                        let island_one = (6..=11).contains(&index);
                        let carried = if island_one {
                            first_clade
                        } else {
                            !first_clade
                        };
                        if carried && rng.next() % 100 < 92 {
                            1.0
                        } else {
                            0.0
                        }
                    }
                })
                .collect();
            MatrixRow::new(format!("ERR500{isolate}"), values)
        })
        .collect();

    let region = Region::new("NC_000962.3", START, START + WINDOW).unwrap();
    let figure = Figure::new(region)
        .title("An accessory genome, sorted by the tree beside it")
        .width(940.0)
        .push(
            FeatureTrack::new(vec![
                Feature::new(START + 13_000, START + 26_000)
                    .name("island A")
                    .strand(Strand::Forward),
                Feature::new(START + 39_000, START + 50_000)
                    .name("island B")
                    .strand(Strand::Reverse),
            ])
            .label("islands"),
        )
        .push(
            MatrixTrack::new(genes, rows)
                .tree(tree)
                .scale(CellScale::Sequential {
                    max: Some(1.0),
                    hue: None,
                })
                .row_height(15.0)
                .min_cell_width(11.0)
                .label("presence"),
        )
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-pangenome.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-pangenome.svg {width:.0} x {height:.0}");
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
