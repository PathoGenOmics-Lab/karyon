//! End to end checks on the rendered document.
//!
//! The unit tests inside the crate cover the arithmetic. These cover the thing
//! a user actually gets: a string of SVG that has to be well formed, stable
//! between runs, free of non-finite numbers, and correct about where a base
//! lands on the page.

use std::fs;

use karyon::{
    AxisTrack, CoverageStyle, CoverageTrack, Feature, FeatureTrack, Figure, LogoColumn,
    LogoScaling, LogoTrack, Region, SequenceTrack, Strand, Theme, Variant, VariantStyle,
    VariantTrack,
};

fn demo_figure() -> Figure {
    let region = Region::parse("NC_000962.3:761001-763000").unwrap();
    let depth: Vec<f64> = (0..2000).map(|i| 40.0 + (i % 31) as f64).collect();
    let bases: Vec<u8> = b"ACGTN".iter().cycle().take(2000).copied().collect();

    Figure::new(region)
        .title("integration figure")
        .push(CoverageTrack::new(761_000, depth).label("depth"))
        .push(SequenceTrack::new(761_000, bases).label("reference"))
        .push(
            FeatureTrack::new(vec![
                Feature::new(761_100, 761_900)
                    .name("geneA")
                    .strand(Strand::Forward),
                Feature::new(761_500, 762_400)
                    .name("geneB")
                    .strand(Strand::Reverse),
            ])
            .label("genes"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(761_250).value(0.9).category("missense"),
                Variant::new(762_100).value(0.4).category("synonymous"),
            ])
            .label("variants"),
        )
        .push(AxisTrack::new())
}

/// Counts opening and closing tags, ignoring self-closing ones.
fn tag_balance(svg: &str, tag: &str) -> (usize, usize) {
    let open = svg.matches(&format!("<{tag} ")).count() + svg.matches(&format!("<{tag}>")).count();
    let close = svg.matches(&format!("</{tag}>")).count();
    (open, close)
}

#[test]
fn the_document_is_well_formed() {
    let svg = demo_figure().to_svg();

    assert!(svg.starts_with("<svg "));
    assert!(svg.ends_with("</svg>"));
    for tag in ["svg", "g", "text", "defs", "clipPath"] {
        let (open, close) = tag_balance(&svg, tag);
        assert_eq!(
            open, close,
            "unbalanced <{tag}>: {open} open, {close} close"
        );
    }
    assert_eq!(svg.matches("<svg ").count(), 1);
}

#[test]
fn no_non_finite_number_reaches_the_output() {
    let svg = demo_figure().to_svg();
    for poison in ["NaN", "nan", "inf", "Infinity"] {
        assert!(!svg.contains(poison), "found {poison} in the output");
    }
}

#[test]
fn rendering_is_deterministic() {
    let first = demo_figure().to_svg();
    let second = demo_figure().to_svg();
    assert_eq!(first, second);
}

#[test]
fn every_track_is_clipped_and_the_clip_ids_are_unique() {
    let svg = demo_figure().to_svg();
    let ids: Vec<&str> = svg
        .match_indices("<clipPath id=\"")
        .map(|(index, prefix)| {
            let rest = &svg[index + prefix.len()..];
            &rest[..rest.find('"').unwrap()]
        })
        .collect();
    assert_eq!(ids.len(), 5, "one clip path per track");
    let mut unique = ids.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(unique.len(), ids.len(), "duplicate clip path id");
}

#[test]
fn a_variant_lands_on_the_centre_of_its_base() {
    // A 1000 bp region across the default 900 px figure, with no track labels
    // so the plotting area is the full width minus the margins.
    let region = Region::parse("chr1:1-1000").unwrap();
    let svg = Figure::new(region)
        .show_region_label(false)
        .push(VariantTrack::new(vec![Variant::new(499)]).show_legend(false))
        .to_svg();

    let plot_x = 12.0;
    let plot_width = 900.0 - 12.0 - 16.0;
    let expected = plot_x + (499.5 / 1000.0) * plot_width;

    let start = svg.find("<circle cx=\"").expect("no variant head drawn");
    let rest = &svg[start + "<circle cx=\"".len()..];
    let cx: f64 = rest[..rest.find('"').unwrap()].parse().unwrap();
    assert!(
        (cx - expected).abs() < 0.01,
        "variant drawn at {cx}, expected {expected}"
    );
}

#[test]
fn user_supplied_names_cannot_break_the_document() {
    let region = Region::new("chr<1> & \"friends\"", 0, 100).unwrap();
    let svg = Figure::new(region)
        .title("a & b")
        .push(FeatureTrack::new(vec![
            Feature::new(10, 90).name("gene<script>alert(1)</script>")
        ]))
        .to_svg();

    assert!(!svg.contains("<script>"));
    assert!(svg.contains("&lt;script&gt;"));
    assert!(svg.contains("&amp;"));
    let (open, close) = tag_balance(&svg, "text");
    assert_eq!(open, close);
}

#[test]
fn a_genome_wide_figure_stays_small() {
    // Four megabases of per-base depth, the size of a bacterial genome.
    let span = 4_000_000usize;
    let depth: Vec<f64> = (0..span).map(|i| ((i % 97) as f64) + 10.0).collect();
    let region = Region::new("NC_000962.3", 0, span as u64).unwrap();

    let svg = Figure::new(region)
        .push(CoverageTrack::new(0, depth).label("depth"))
        .push(SequenceTrack::new(0, vec![b'A'; span]).label("reference"))
        .push(AxisTrack::new())
        .to_svg();

    // One point per pixel column, not one per base.
    assert!(
        svg.len() < 100_000,
        "genome wide figure grew to {} bytes",
        svg.len()
    );
    assert!(svg.contains("Mb"));
}

#[test]
fn styles_and_themes_change_the_output_without_breaking_it() {
    let region = Region::parse("chr1:1-500").unwrap();
    let depth: Vec<f64> = (0..500).map(|i| (i % 40) as f64).collect();

    for style in [
        CoverageStyle::Area,
        CoverageStyle::Line,
        CoverageStyle::Bars,
    ] {
        for theme in [Theme::light(), Theme::dark()] {
            let svg = Figure::new(region.clone())
                .theme(theme)
                .push(CoverageTrack::new(0, depth.clone()).style(style))
                .push(VariantTrack::new(vec![Variant::new(100)]).style(VariantStyle::Tick))
                .push(AxisTrack::new())
                .to_svg();
            assert!(svg.starts_with("<svg "));
            assert!(svg.ends_with("</svg>"));
            assert!(!svg.contains("NaN"));
        }
    }
}

#[test]
fn save_svg_writes_exactly_what_to_svg_returns() {
    let figure = demo_figure();
    let path = std::env::temp_dir().join(format!("karyon-{}.svg", std::process::id()));
    figure.save_svg(&path).unwrap();
    let written = fs::read_to_string(&path).unwrap();
    fs::remove_file(&path).unwrap();
    assert_eq!(written, figure.to_svg());
}

#[test]
fn a_logo_figure_is_well_formed_and_deterministic() {
    let alignment = [
        "ACGTACGT", "ACGTACGA", "ACGAACGT", "ACGTTCGT", "ACGTACGT", "ACCTACGA",
    ];
    let figure = || {
        Figure::new(Region::new("motif", 0, 8).unwrap())
            .title("motif")
            .push(
                LogoTrack::from_sequences(0, &alignment)
                    .alphabet_size(4)
                    .scaling(LogoScaling::InformationContent)
                    .label("bits"),
            )
            .push(
                LogoTrack::from_sequences(0, &alignment)
                    .alphabet_size(4)
                    .scaling(LogoScaling::EnrichmentDepletion)
                    .label("enrich"),
            )
            .push(AxisTrack::new().center_on_bases(true))
            .to_svg()
    };

    let svg = figure();
    assert!(svg.starts_with("<svg "));
    assert!(svg.ends_with("</svg>"));
    assert!(!svg.contains("NaN"));
    // Glyphs are stretched to their boxes rather than set at a font size.
    assert!(svg.contains(r#"lengthAdjust="spacingAndGlyphs""#));
    for tag in ["svg", "g", "text", "clipPath"] {
        let (open, close) = tag_balance(&svg, tag);
        assert_eq!(open, close, "unbalanced <{tag}>");
    }
    assert_eq!(svg, figure());
}

#[test]
fn only_the_enrichment_logo_can_show_an_absent_base() {
    // Near uniform over three bases, with the fourth never observed.
    let column = vec![LogoColumn::acgt(34.0, 33.0, 33.0, 0.0)];
    let classic = LogoTrack::new(0, column.clone()).alphabet_size(4);
    let edlogo = LogoTrack::new(0, column)
        .alphabet_size(4)
        .scaling(LogoScaling::EnrichmentDepletion);

    assert!(classic.stacks()[0].down_total() == 0.0);
    assert!(classic.stacks()[0].up_total() < 0.45, "should look empty");
    assert!(edlogo.stacks()[0].down_total() > 3.0, "should be loud");
    assert_eq!(edlogo.stacks()[0].down[0].0, "T");
}

#[test]
fn an_empty_data_set_still_produces_a_figure() {
    let region = Region::parse("chr1:1-1000").unwrap();
    let svg = Figure::new(region)
        .push(CoverageTrack::new(0, Vec::new()).label("depth"))
        .push(FeatureTrack::new(Vec::new()).label("genes"))
        .push(VariantTrack::new(Vec::new()).label("variants"))
        .push(AxisTrack::new())
        .to_svg();
    assert!(svg.starts_with("<svg "));
    assert!(svg.ends_with("</svg>"));
    assert!(svg.contains("chr1:1-1000"));
}
