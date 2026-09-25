//! The figure `cargo run --example selection` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page, in the page's colours, at
//! the width of the column and over whatever window the reader has moved to.
//! One copy of the figure serves both, so the two cannot drift apart.

use karyon::{Drawing, Plot, Region, Theme, Window, WindowStyle, WindowTrack};

/// Start of the window, 0-based.
const START: u64 = 2_150_000;
/// Length of the window in bases.
const SPAN: u64 = 40_000;
/// Size of one pN/pS window.
const STEP: u64 = 500;

/// `example-selection.svg`: selection and strand composition, read against
/// their baselines.
///
/// `theme`, `width` and `region` replace the light theme, the 880 pixels and
/// the forty kilobase window the committed figure is drawn with. The colours
/// above and below the pN/pS baseline come from `theme` as well.
pub fn example_selection(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
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

    let own = Region::new("NC_000962.3", START, START + SPAN).unwrap();
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("Selection and strand composition, read against their baselines")
        .theme(theme.clone())
        .width(width.unwrap_or(880.0))
        .add_track(
            WindowTrack::ratios(windows)
                .label("pN/pS")
                .height(70.0)
                .colors(theme.color(1), theme.color(0)),
        )
        .add_track(
            WindowTrack::gc_skew(START, &bases, 1_000)
                .style(WindowStyle::Line)
                .label("GC skew")
                .height(56.0),
        )
        .into_figure();
    Box::new(figure)
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
