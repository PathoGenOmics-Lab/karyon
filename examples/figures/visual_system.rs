//! The figure `cargo run --example visual_system` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page. One copy of the figure
//! serves both, so the two cannot drift apart.

use karyon::{
    AxisFormat, AxisTrack, CoverageTrack, Drawing, Emphasis, Feature, FeatureTrack, Figure,
    LinePattern, Panels, QuantitativeAxis, ReferenceLine, Region, RenderProfile, Strand, Theme,
    Variant, VariantTrack, Window, WindowStyle, WindowTrack,
};

/// `example-visual-system.svg`: one visual grammar across genomic tracks.
///
/// The sheet shows what the web profile does, palette included, so it keeps
/// the profile it sets and `theme` is ignored. A sheet is as wide as its
/// panels and has no window of its own, so `width` and `region` are ignored
/// too.
pub fn example_visual_system(
    _theme: &Theme,
    _width: Option<f64>,
    _region: Option<&Region>,
) -> Box<dyn Drawing> {
    let region = Region::new("NC_000962.3", 760_999, 762_999).unwrap();

    let depth: Vec<f64> = (0..2_000)
        .map(|index| {
            let wave = 62.0 + 11.0 * (index as f64 / 95.0).sin();
            if (880..1_060).contains(&index) {
                wave * 0.18
            } else {
                wave
            }
        })
        .collect();
    let depth_axis = QuantitativeAxis::new()
        .range(0.0, 100.0)
        .ticks(3)
        .unit("x")
        .reference(
            ReferenceLine::new(30.0)
                .label("QC 30x")
                .emphasis(Emphasis::Alert)
                .pattern(LinePattern::Dashed),
        );

    let overview = Figure::new(region.clone())
        .profile(RenderProfile::Web)
        .width(780.0)
        .title("One visual grammar across genomic tracks")
        .push(
            CoverageTrack::new(region.start(), depth)
                .axis(depth_axis)
                .label("read depth after filtering"),
        )
        .push(
            FeatureTrack::new(vec![
                Feature::new(759_806, 763_325)
                    .name("rpoB")
                    .strand(Strand::Forward),
                Feature::new(761_081, 761_162)
                    .name("RRDR")
                    .strand(Strand::Forward),
            ])
            .label("annotation"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(761_108).value(0.98).category("missense"),
                Variant::new(761_138).value(0.55).category("indel"),
                Variant::new(761_154).value(1.00).category("synonymous"),
                Variant::new(761_155).value(0.21).category("promoter"),
            ])
            .axis(
                QuantitativeAxis::new()
                    .range(0.0, 1.0)
                    .format(AxisFormat::Percent(0)),
            )
            .label("allele frequency"),
        )
        .push(AxisTrack::new().label("position"));

    let windows: Vec<Window> = (0..40)
        .map(|index| {
            let start = region.start() + index * 50;
            Window::new(start, start + 50, (index as f64 / 4.5).sin())
        })
        .collect();
    let centred = Figure::new(region)
        .profile(RenderProfile::Web)
        .width(780.0)
        .title("Shared axes, reference lines and non-colour encodings")
        .push(
            WindowTrack::new(windows)
                .style(WindowStyle::Line)
                .axis(
                    QuantitativeAxis::new()
                        .range(-1.0, 1.0)
                        .ticks(5)
                        .format(AxisFormat::Fixed(1))
                        .reference(ReferenceLine::new(0.5).label("upper guide")),
                )
                .label("strand-composition score"),
        )
        .push(AxisTrack::new().label("position"));

    let sheet = Panels::new()
        .title("Karyon visual system")
        .push_captioned(&overview, "A", "Profile, annotation and categorical events")
        .push_captioned(
            &centred,
            "B",
            "Centred statistic with a shared quantitative axis",
        );
    Box::new(sheet)
}
