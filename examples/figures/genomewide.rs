//! The figure `cargo run --example genomewide` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page, in the page's colours, at
//! the width of the column and over whatever window the reader has moved to.
//! One copy of the figure serves both, so the two cannot drift apart.

use karyon::{Association, Drawing, Genome, Plot, Region, Theme};

/// `example-genomewide.svg`: a scan across a whole draft assembly.
///
/// `theme`, `width` and `region` replace the light theme, the 940 pixels and
/// the whole assembly the committed figure is drawn with.
pub fn example_genomewide(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
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
    let (mapped, _) = genome.map(hits);
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

    let figure = Plot::over(region.cloned().unwrap_or_else(|| genome.region()))
        .title("An association scan across a whole draft assembly")
        .theme(theme.clone())
        .width(width.unwrap_or(940.0))
        .remove_region_label()
        .add_manhattan(points)
        .label("association")
        .adjust(|track| {
            track
                .bands(genome.boundaries())
                .genome_wide_threshold()
                .axis_title("-log10 p")
                .height(96.0)
        })
        // From the first base of the assembly rather than from the left edge
        // of the window, which is only the same place until the window moves.
        .add_coverage_at(0, dense)
        .label("depth")
        .adjust(|track| track.height(52.0))
        .add_genome(genome)
        .label("contigs")
        .into_figure();
    Box::new(figure)
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
