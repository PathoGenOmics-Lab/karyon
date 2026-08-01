//! Renders the selection figure used in the README.
//!
//! ```text
//! cargo run --example selection -- assets
//! ```
//!
//! Two statistics that only mean anything relative to a line. pN/pS says
//! whether a stretch of gene is under purifying or diversifying selection, and
//! the answer is which side of one it fell on; GC skew changes sign at the
//! origin and the terminus of replication, and the sign is the whole signal.
//! Both are plotted as departures from a baseline rather than as heights above
//! a floor, which is what a coverage track would make of them.

use std::env;
use std::path::PathBuf;

use karyon::{Plot, Region, Window, WindowStyle, WindowTrack};

/// Start of the window, 0-based.
const START: u64 = 2_150_000;
/// Length of the window in bases.
const SPAN: u64 = 40_000;
/// Size of one pN/pS window.
const STEP: u64 = 500;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let mut rng = Lcg::new(90_210);

    // Mostly purifying, as a bacterial chromosome is, with one stretch under
    // diversifying selection in the middle of it.
    let windows: Vec<Window> = (0..(SPAN / STEP))
        .map(|index| {
            let start = START + index * STEP;
            let diversifying = (28..36).contains(&index);
            let noise = (rng.next() % 1000) as f64 / 1000.0;
            let ratio = if diversifying {
                1.4 + 1.6 * noise
            } else {
                0.08 + 0.42 * noise
            };
            Window::new(start, start + STEP, ratio)
        })
        .collect();

    // A sequence whose skew turns over a third of the way along, as it does
    // either side of a replication origin.
    let bases: Vec<u8> = (0..SPAN)
        .map(|index| {
            let leading = index < SPAN / 3;
            match rng.next() % 100 {
                0..=32 if leading => b'G',
                0..=32 => b'C',
                33..=64 if leading => b'C',
                33..=64 => b'G',
                65..=82 => b'A',
                _ => b'T',
            }
        })
        .collect();

    let figure = Plot::over(Region::new("NC_000962.3", START, START + SPAN).unwrap())
        .title("Selection and strand composition, read against their baselines")
        .width(880.0)
        .add_track(
            WindowTrack::ratios(windows)
                .label("pN/pS")
                .height(70.0)
                .colors("#d55e00", "#0072b2"),
        )
        .add_track(
            WindowTrack::gc_skew(START, &bases, 1_000)
                .style(WindowStyle::Line)
                .label("GC skew")
                .height(56.0),
        )
        .into_figure();

    figure.save_svg(out.join("example-selection.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-selection.svg {width:.0} x {height:.0}");
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
