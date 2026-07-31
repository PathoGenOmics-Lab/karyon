//! Renders the distance matrix figure used in the README.
//!
//! ```text
//! cargo run --example distance -- assets
//! ```
//!
//! Pairwise SNP distances between isolates, ordered by the tree beside them and
//! with every pair inside a twelve SNP threshold outlined. Sorted by accession
//! the two clusters are confetti; sorted by descent they are two blocks on the
//! diagonal, and the outlines say which pairs the threshold actually admits.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{DistanceTrack, Figure, Region};

/// SNP distance under which a pair is called recent transmission.
///
/// A convention for *M. tuberculosis*, not a fact: it belongs to a mutation
/// rate and a sampling window, and it goes in the caption.
const THRESHOLD: f64 = 12.0;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    // Two clusters of four, and two isolates that belong to neither.
    let names: Vec<String> = [
        "ERR5001", "ERR5002", "ERR5003", "ERR5004", "ERR5005", "ERR5006", "ERR5007", "ERR5008",
        "ERR5009", "ERR5010",
    ]
    .iter()
    .map(|name| name.to_string())
    .collect();

    let cluster = |index: usize| match index {
        0..=3 => Some(0),
        4..=7 => Some(1),
        _ => None,
    };

    let mut rng = Lcg::new(31_415);
    let size = names.len();
    let mut distances = vec![0.0f64; size * size];
    for from in 0..size {
        for to in (from + 1)..size {
            let together = cluster(from).is_some() && cluster(from) == cluster(to);
            let value = if together {
                (rng.next() % 11) as f64
            } else {
                380.0 + (rng.next() % 240) as f64
            };
            distances[from * size + to] = value;
            distances[to * size + from] = value;
        }
    }

    let tree = Tree::parse_newick(
        "(((ERR5001:0.0002,ERR5003:0.0002)0.99:0.0009,\
           (ERR5002:0.0003,ERR5004:0.0001)0.97:0.0008)0.92:0.0140,\
          ((ERR5009:0.0180,ERR5010:0.0165):0.0060,\
           ((ERR5005:0.0002,ERR5007:0.0003)0.98:0.0011,\
            (ERR5006:0.0004,ERR5008:0.0002)0.95:0.0009):0.0130):0.0040);",
    )
    .expect("the tree in this example is well formed");

    let figure = Figure::new(Region::new("samples", 0, size as u64).unwrap())
        .title("Pairwise SNP distances, ordered by descent")
        // Cells are square and capped, so the figure is only as wide as the
        // matrix needs: ten samples at thirty pixels plus the tree and the
        // names beside them.
        .width(560.0)
        .show_region_label(false)
        .push(
            DistanceTrack::new(names, distances)
                .tree(tree)
                .cluster_threshold(THRESHOLD)
                .cell_size(30.0)
                .label("SNP distance"),
        );

    figure.save_svg(out.join("example-distance.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-distance.svg {width:.0} x {height:.0}, threshold {THRESHOLD:.0} SNPs");
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
