//! The figure `cargo run --example copy_number` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page, in the page's colours, at
//! the width of the column and over whatever window the reader has moved to.
//! One copy of the figure serves both, so the two cannot drift apart.

use karyon::{
    CopyNumberSegment, CopyNumberTrack, Drawing, Plot, QuantitativeAxis, Region, Theme, Window,
    WindowTrack,
};

/// Length of the arm the figure is drawn over.
const SPAN: u64 = 46_000_000;
/// How many samples the cohort holds, which is the denominator of the landscape.
const COHORT: f64 = 120.0;

/// `example-copy-number.svg`: a cohort's copy number, and one sample under it.
///
/// `theme`, `width` and `region` replace the light theme, the 900 pixels and
/// the whole arm the committed figure is drawn with.
pub fn example_copy_number(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let own = Region::new("chr8", 0, SPAN).unwrap();

    // The cohort. Two rows over each span, so a place that went both ways is
    // drawn going both ways rather than netted into a number nobody measured.
    let mut landscape: Vec<Window> = Vec::new();
    for i in 0..230u64 {
        let from = i * 200_000;
        let to = from + 200_000;
        let at = i as f64 / 230.0;
        // A recurrent amplification over the distal arm, and a recurrent loss
        // proximal to it, with one stretch that does both.
        let gained = (0.06 + 0.62 * ((at - 0.78) / 0.10).powi(2).exp().recip()).min(0.74);
        let lost = (0.05 + 0.40 * ((at - 0.30) / 0.13).powi(2).exp().recip()).min(0.48);
        landscape.push(Window::new(from, to, gained));
        landscape.push(Window::new(from, to, -lost));
    }

    let cohort = WindowTrack::new(landscape)
        // Warm for gained and cool for lost, the field's convention, which is
        // the other way round from the default: the default is about a signed
        // statistic, and here the sign is which of two things happened.
        .colors(theme.color(1), theme.color(0))
        .axis(QuantitativeAxis::new().range(-1.0, 1.0).ticks(3))
        .label(format!("{COHORT:.0} samples"))
        .height(84.0);

    // One tumour, allele-specific. The arm it sits on is diploid, and what the
    // ladder cannot show is the stretch that kept two copies and lost one of
    // its two alleles: only the lane along the foot says so.
    let segments = vec![
        CopyNumberSegment::allelic(0, 8_400_000, 1.0, 1.0),
        CopyNumberSegment::allelic(8_400_000, 13_100_000, 1.0, 0.0),
        // Copy neutral loss of heterozygosity: two copies, both from one
        // allele, sitting exactly on the rule that means unchanged.
        CopyNumberSegment::allelic(13_100_000, 19_700_000, 2.0, 0.0),
        CopyNumberSegment::allelic(19_700_000, 27_300_000, 1.0, 1.0),
        // A homozygous deletion, which is a real state and not a missing value.
        CopyNumberSegment::allelic(27_300_000, 27_900_000, 0.0, 0.0),
        CopyNumberSegment::allelic(27_900_000, 35_200_000, 1.0, 1.0),
        // A focal amplification a fortieth of a pixel wide at this zoom, drawn
        // at its own level with a riser down to the arm rather than averaged
        // into it.
        CopyNumberSegment::allelic(35_260_000, 35_290_000, 9.0, 1.0),
        CopyNumberSegment::allelic(35_290_000, 38_600_000, 3.0, 1.0),
        // A stretch the caller could not resolve the alleles over. The total is
        // known, the split is not, and the lane says which of the two it is.
        CopyNumberSegment::total(38_600_000, 42_000_000, 3.0),
        CopyNumberSegment::allelic(42_000_000, SPAN, 4.0, 2.0),
    ];

    let tumour = CopyNumberTrack::diploid(segments)
        .label("one tumour")
        .height(150.0);

    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("Recurrent copy number across a cohort, and one sample under it")
        .theme(theme.clone())
        .width(width.unwrap_or(900.0))
        .add_track(cohort)
        .add_track(tumour)
        .into_figure();
    Box::new(figure)
}
