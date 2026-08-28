//! Renders the copy number figure used in the documentation.
//!
//! ```text
//! cargo run --example copy_number -- assets
//! ```
//!
//! One tumour's segmentation over a chromosome arm, and above it how often the
//! same arm was gained and lost across the cohort it came from. The two bands
//! are the two questions a copy number analysis is asked: what happened in this
//! sample, and does it happen in general.
//!
//! The cohort landscape is a [`WindowTrack`] with two rows over every span, one
//! for the gains and one for the losses, because a locus can be gained in a
//! third of a cohort and lost in a fifth of it at the same time and one signed
//! number cannot say that. It is not a track type of its own, and the reason it
//! is not is that this already draws it.

use std::env;
use std::path::PathBuf;

use karyon::{
    CopyNumberSegment, CopyNumberTrack, Plot, QuantitativeAxis, Region, Window, WindowTrack,
};

/// Length of the arm the figure is drawn over.
const SPAN: u64 = 46_000_000;
/// How many samples the cohort holds, which is the denominator of the landscape.
const COHORT: f64 = 120.0;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let region = Region::new("chr8", 0, SPAN).unwrap();

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
        .colors("#d55e00", "#0072b2")
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
    let loh = tumour.loh_spans().len();
    let deletions = tumour.homozygous_deletions().len();
    let unresolved = tumour.without_allele_call();

    let figure = Plot::over(region.clone())
        .title("Recurrent copy number across a cohort, and one sample under it")
        .width(900.0)
        .add_track(cohort)
        .add_track(tumour)
        .into_figure();

    figure.save_svg(out.join("example-copy-number.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-copy-number.svg {width:.0} x {height:.0}");
    println!(
        "{loh} spans lost heterozygosity, {deletions} lost both copies, \
         {unresolved} segment had no allele split"
    );
    Ok(())
}
