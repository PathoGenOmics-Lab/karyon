//! Renders the genome-wide figures used in the README.
//!
//! ```text
//! cargo run --example genomewide -- assets
//! ```
//!
//! Two things a figure over one region cannot show. A whole draft assembly,
//! every contig of it on one axis, with the association scan and the depth
//! profile laid over all of them at once. And the gene frequency spectrum, the
//! U that says what a pangenome is made of rather than whether it closes.

use std::env;
use std::path::{Path, PathBuf};

use karyon::{
    Association, AxisTrack, CoverageTrack, Figure, Frequency, FrequencyTrack, Genome, GenomeTrack,
    Legend, LegendTrack, ManhattanTrack, Region,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    genome_wide(&out)?;
    spectrum(&out)?;
    Ok(())
}

/// A scan across a whole draft assembly.
fn genome_wide(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(11_235);

    // Twelve contigs, longest first, as an assembly comes out of a scaffolder.
    let genome = Genome::new(
        (0..12)
            .map(|index| {
                let length = 900_000 / (index + 1) + 40_000;
                (format!("contig_{:02}", index + 1), length)
            })
            .collect::<Vec<_>>(),
    );

    // A scan reported the way one arrives: per contig, in its own coordinates.
    let mut hits: Vec<(String, u64, f64)> = Vec::new();
    for contig in genome.sequences() {
        let mut at = 200u64;
        while at < contig.length {
            // One peak, on the fourth contig.
            let peak = contig.name == "contig_04" && at.abs_diff(120_000) < 9_000;
            let noise = (rng.next() % 1000) as f64 / 1000.0;
            let value = if peak {
                5.0 + 4.5 * noise
            } else {
                0.2 + 2.6 * noise
            };
            hits.push((contig.name.clone(), at, value));
            at += 900 + rng.next() % 2_600;
        }
    }
    let (mapped, dropped) = genome.map(hits);
    let points: Vec<Association> = mapped
        .iter()
        .map(|(at, value)| Association::new(*at, *value))
        .collect();

    // Depth over the same axis, dropping off at the ends of every contig the
    // way real coverage does.
    let total = genome.total();
    let depth: Vec<f64> = (0..total)
        .step_by(500)
        .map(|at| {
            let edge = genome
                .locate(at)
                .map(|(name, offset)| {
                    let length = genome
                        .sequences()
                        .iter()
                        .find(|seq| seq.name == name)
                        .map_or(1, |seq| seq.length);
                    offset.min(length - offset) as f64 / 12_000.0
                })
                .unwrap_or(0.0);
            42.0 * edge.min(1.0) + (rng.next() % 60) as f64 / 10.0
        })
        .collect();
    // One value per five hundred bases, spread back over the axis.
    let dense: Vec<f64> = depth
        .iter()
        .flat_map(|value| std::iter::repeat(*value).take(500))
        .take(total as usize)
        .collect();

    let figure = Figure::new(genome.region())
        .title("An association scan across a whole draft assembly")
        .width(940.0)
        .show_region_label(false)
        .push(
            ManhattanTrack::new(points)
                .bands(genome.boundaries())
                .genome_wide_threshold()
                .unit(" -log10 p")
                .label("association")
                .height(96.0),
        )
        .push(CoverageTrack::new(0, dense).label("depth").height(52.0))
        .push(GenomeTrack::new(genome.clone()).label("contigs"))
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-genomewide.svg"))?;
    let (width, height) = figure.dimensions();
    println!(
        "example-genomewide.svg {width:.0} x {height:.0}, {} contigs over {:.2} Mb, {dropped} unplaced",
        genome.len(),
        genome.total() as f64 / 1e6
    );
    Ok(())
}

/// What the pangenome is made of.
fn spectrum(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(90_210);
    let count = 60usize;
    let core = 3_100usize;
    let shell = 900usize;
    let cloud = 1_500usize;

    let genomes: Vec<Vec<bool>> = (0..count)
        .map(|genome| {
            let mut carried = Vec::with_capacity(core + shell + cloud);
            carried.extend(std::iter::repeat(true).take(core));
            // A shell family sits in a middling share of the collection, which
            // is what fills the trough of the U.
            for family in 0..shell {
                let share = 20 + (family * 60 / shell);
                carried.push((rng.next() % 100) < share as u64);
            }
            // A cloud family is in one or two, which is what builds the tower
            // on the left.
            for family in 0..cloud {
                carried.push(family % count == genome || family % (count * 3) == genome);
            }
            carried
        })
        .collect();

    let track = FrequencyTrack::from_presence(&genomes).label("gene families");
    let (core_total, shell_total, cloud_total) = (
        track.total(Frequency::Core),
        track.total(Frequency::Shell),
        track.total(Frequency::Cloud),
    );

    let legend = Legend::new()
        .area(
            format!("core, {core_total} families"),
            track.color(Frequency::Core, &Default::default()),
        )
        .area(
            format!("shell, {shell_total}"),
            track.color(Frequency::Shell, &Default::default()),
        )
        .area(
            format!("cloud, {cloud_total}"),
            track.color(Frequency::Cloud, &Default::default()),
        );

    let figure = Figure::new(Region::new("genomes", 0, count as u64).unwrap())
        .title("What the pangenome is made of")
        .width(880.0)
        .show_region_label(false)
        .push(track.height(170.0))
        .push(LegendTrack::new(legend))
        .push(AxisTrack::new().center_on_bases(true).label("genomes"));

    figure.save_svg(out.join("example-spectrum.svg"))?;
    let (width, height) = figure.dimensions();
    println!(
        "example-spectrum.svg {width:.0} x {height:.0}, \
         core {core_total}, shell {shell_total}, cloud {cloud_total}"
    );
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
