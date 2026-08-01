//! Renders the association figures used in the README.
//!
//! ```text
//! cargo run --example association -- assets
//! ```
//!
//! A locus that came out of a genome-wide association scan, and the isolates
//! behind it. The Manhattan panel says where the signal is; the matrix says who
//! carries it, which is the question the Manhattan panel always provokes.

use std::env;
use std::path::PathBuf;

use karyon::{Association, CellScale, Feature, MatrixRow, Plot, Region, Strand};

/// Window start, 0-based, around the rpoB locus of H37Rv.
const START: u64 = 759_000;
/// Window length in bases.
const WINDOW: u64 = 6_000;
/// Where the association peaks.
const PEAK: u64 = 761_155;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let mut rng = Lcg::new(4_242);

    // A scan over the window: background noise everywhere, and a tower of
    // linked markers around the peak, which is what a real signal looks like.
    let mut points = Vec::new();
    let mut sites = Vec::new();
    let mut position = START + 40;
    while position < START + WINDOW {
        let distance = position.abs_diff(PEAK) as f64;
        let linkage = (-distance / 260.0).exp();
        let noise = (rng.next() % 1000) as f64 / 1000.0;
        // Background is a uniform p-value; the peak lifts what is near it.
        let value = 0.3 + 2.4 * noise + 9.5 * linkage * (0.55 + 0.45 * noise);
        points.push(Association::new(position, value));
        if value > 3.0 {
            sites.push(position);
        }
        position += 30 + rng.next() % 40;
    }
    sites.truncate(14);

    // Isolates typed at the sites that passed. The ones carrying the haplotype
    // are the reason the tower is there.
    let rows: Vec<MatrixRow> = (0..10)
        .map(|isolate| {
            let carrier = isolate % 3 != 2;
            let values: Vec<f64> = sites
                .iter()
                .map(|_| {
                    if rng.next() % 100 < 4 {
                        f64::NAN
                    } else if carrier && rng.next() % 100 < 88 {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect();
            MatrixRow::new(format!("ERR{:04}", 3_100 + isolate * 7), values)
        })
        .collect();

    let carriers = rows
        .iter()
        .filter(|row| row.values.iter().filter(|v| **v == 1.0).count() > sites.len() / 2)
        .count();

    let region = Region::new("NC_000962.3", START, START + WINDOW).unwrap();
    let figure = Plot::over(region)
        .title("An association, and the isolates behind it")
        .width(940.0)
        .add_manhattan(points)
        .label("association")
        .adjust(|track| track.genome_wide_threshold().unit(" -log10 p").height(96.0))
        .add_features(vec![Feature::new(759_806, 763_325)
            .name("rpoB")
            .strand(Strand::Forward)])
        .label("gene")
        .add_matrix(sites.clone(), rows)
        .label("genotypes")
        .adjust(|track| {
            track
                .scale(CellScale::Sequential {
                    max: Some(1.0),
                    hue: None,
                })
                .min_cell_width(9.0)
        })
        .into_figure();

    figure.save_svg(out.join("example-association.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-association.svg {width:.0} x {height:.0}");
    println!(
        "{} sites past the threshold, {carriers} of 10 isolates carrying the haplotype",
        sites.len()
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
