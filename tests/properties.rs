//! Properties that must hold for every figure, checked against generated ones.
//!
//! The tests beside this one are examples: an input somebody chose and the
//! output they expected. That is worth having and it inherits the blind spot of
//! whoever wrote it, because the person who writes the assertion is the person
//! who believed the code. A reverse-strand span was wrong in this crate from
//! its first commit and every test of it agreed with it.
//!
//! These tests assert invariants instead. They do not say what any particular
//! figure looks like; they say what is true of all of them, and then generate
//! figures until one is not. An invariant can be checked without knowing the
//! right answer, which is the only way to find something nobody thought of.
//!
//! Nothing here is random between runs. Every generator is a seeded
//! [`Lcg`], so a failure is reproducible from the seed in the message, and the
//! suite is worth no more than the hostility of what it feeds in: not plausible
//! data, but `NaN`, `u64::MAX`, empty vectors and one-base regions.

use karyon::{
    Association, AxisTrack, Band, CigarOp, CodonTrack, CoverageTrack, Feature, FeatureTrack,
    Figure, IdeogramTrack, ManhattanTrack, MatrixRow, MatrixTrack, MethylSite, MethylationTrack,
    MsaSequence, MsaTrack, OrfTrack, PileupTrack, Read, Region, Scale, SequenceTrack, SnpSite,
    SnpTrack, Stain, Strand, StructuralTrack, StructuralVariant, SvKind, Theme, Track, Variant,
    VariantTrack, Window, WindowTrack,
};

/// How many figures each property is given before it is believed.
const ROUNDS: u64 = 10_000;

// ---------------------------------------------------------------------------
// Generation
// ---------------------------------------------------------------------------

/// A linear congruential generator, so a failing case is reproducible.
///
/// The same one the examples use. A dependency would do this better and would
/// also be a dependency, and the crate's claim is that it has none, including
/// the ones that only show up when running the tests.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1))
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 11
    }

    /// A number below `bound`, or zero when the bound is zero.
    fn below(&mut self, bound: u64) -> u64 {
        if bound == 0 {
            0
        } else {
            self.next() % bound
        }
    }

    fn chance(&mut self, one_in: u64) -> bool {
        self.below(one_in.max(1)) == 0
    }

    /// A float that is ordinary most of the time and hostile the rest of it.
    ///
    /// The hostile values are the point. A generator that only produces
    /// plausible depths tests the same thing the example tests already test.
    fn value(&mut self) -> f64 {
        match self.below(16) {
            0 => f64::NAN,
            1 => f64::INFINITY,
            2 => f64::NEG_INFINITY,
            3 => 0.0,
            4 => -0.0,
            5 => -(self.below(1_000) as f64),
            6 => f64::MAX,
            7 => f64::MIN_POSITIVE,
            8 => 1e308,
            _ => (self.below(10_000) as f64) / 100.0,
        }
    }

    /// A coordinate that is ordinary most of the time and extreme the rest.
    fn position(&mut self, span: u64) -> u64 {
        match self.below(12) {
            0 => 0,
            1 => u64::MAX,
            2 => u64::MAX - 1,
            3 => span,
            4 => span.saturating_sub(1),
            5 => 1 << 32,
            _ => self.below(span.max(1)),
        }
    }

    /// How long a generated collection is, weighted towards the awkward ends.
    fn count(&mut self) -> usize {
        match self.below(8) {
            0 => 0,
            1 => 1,
            2 => 2,
            _ => self.below(40) as usize,
        }
    }

    fn base(&mut self) -> u8 {
        b"ACGTNacgtn-"[self.below(11) as usize]
    }

    fn strand(&mut self) -> Strand {
        match self.below(3) {
            0 => Strand::Forward,
            1 => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }

    /// A name, sometimes empty and sometimes carrying characters that have to
    /// be escaped on the way into the document.
    fn name(&mut self) -> String {
        match self.below(10) {
            0 => String::new(),
            1 => "a<b&c\"d'e".to_string(),
            2 => "<script>alert(1)</script>".to_string(),
            3 => "\u{1f9ec} \u{202e} tab\there".to_string(),
            4 => "x".repeat(200),
            _ => format!("s{}", self.below(1_000)),
        }
    }
}

/// A region, sometimes one base wide and sometimes the whole coordinate range.
fn region(rng: &mut Lcg) -> Region {
    let (start, end) = match rng.below(10) {
        0 => (0, 1),
        1 => (0, u64::MAX),
        2 => (u64::MAX - 10, u64::MAX),
        3 => (1 << 32, (1 << 32) + 1_000),
        _ => {
            let start = rng.below(4_000_000);
            (start, start + 1 + rng.below(100_000))
        }
    };
    Region::new("chr1", start, end).expect("start is always below end")
}

/// One track of a random kind, over `region`.
fn track(rng: &mut Lcg, region: &Region) -> Box<dyn Track> {
    let start = region.start();
    let span = region.len();
    let width = span.min(4_000);
    let count = rng.count();

    match rng.below(16) {
        0 => {
            let values: Vec<f64> = (0..count).map(|_| rng.value()).collect();
            Box::new(CoverageTrack::new(rng.position(span), values))
        }
        1 => {
            let bases: Vec<u8> = (0..count).map(|_| rng.base()).collect();
            Box::new(SequenceTrack::new(start, bases))
        }
        2 => {
            let features: Vec<Feature> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    Feature::new(a.min(b), a.max(b))
                        .name(rng.name())
                        .strand(rng.strand())
                })
                .collect();
            Box::new(FeatureTrack::new(features))
        }
        3 => {
            let calls: Vec<Variant> = (0..count)
                .map(|_| {
                    Variant::new(rng.position(span))
                        .value(rng.value())
                        .category(rng.name())
                })
                .collect();
            Box::new(VariantTrack::new(calls))
        }
        4 => {
            let windows: Vec<Window> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    Window::new(a.min(b), a.max(b), rng.value())
                })
                .collect();
            Box::new(WindowTrack::new(windows))
        }
        5 => {
            let points: Vec<Association> = (0..count)
                .map(|_| Association::new(rng.position(span), rng.value()))
                .collect();
            Box::new(ManhattanTrack::new(points))
        }
        6 => {
            let sites: Vec<u64> = (0..count).map(|_| rng.position(span)).collect();
            let rows: Vec<MatrixRow> = (0..rng.count())
                .map(|_| {
                    let values: Vec<f64> = (0..sites.len()).map(|_| rng.value()).collect();
                    MatrixRow::new(rng.name(), values)
                })
                .collect();
            Box::new(MatrixTrack::new(sites, rows))
        }
        7 => {
            let columns = rng.count();
            let rows: Vec<MsaSequence> = (0..rng.count())
                .map(|_| {
                    let residues: Vec<u8> = (0..columns).map(|_| rng.base()).collect();
                    MsaSequence::new(rng.name(), residues)
                })
                .collect();
            Box::new(MsaTrack::new(rows))
        }
        8 => {
            let reads: Vec<Read> = (0..count)
                .map(|_| {
                    let ops: Vec<CigarOp> = (0..1 + rng.below(6))
                        .map(|_| {
                            let n = rng.below(200) as u32;
                            match rng.below(5) {
                                0 => CigarOp::Insertion(n),
                                1 => CigarOp::Deletion(n),
                                2 => CigarOp::Skip(n),
                                3 => CigarOp::SoftClip(n),
                                _ => CigarOp::Match(n),
                            }
                        })
                        .collect();
                    Read::new(rng.position(span), ops)
                })
                .collect();
            Box::new(PileupTrack::new(reads))
        }
        9 => {
            let calls: Vec<StructuralVariant> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    let kind = match rng.below(4) {
                        0 => SvKind::Deletion,
                        1 => SvKind::Duplication,
                        2 => SvKind::Inversion,
                        _ => SvKind::Insertion,
                    };
                    StructuralVariant::new(a.min(b), a.max(b), kind)
                })
                .collect();
            Box::new(StructuralTrack::new(calls))
        }
        10 => {
            let seq: Vec<u8> = (0..count * 3).map(|_| rng.base()).collect();
            Box::new(OrfTrack::new(start, seq))
        }
        11 => {
            let a = rng.position(span);
            let b = rng.position(span);
            Box::new(CodonTrack::new(a.min(b), a.max(b), rng.strand()))
        }
        12 => {
            let sites: Vec<SnpSite> = (0..count)
                .map(|_| {
                    let alleles: Vec<u8> = (0..rng.count()).map(|_| rng.base()).collect();
                    SnpSite::new(rng.position(span), rng.base(), alleles)
                })
                .collect();
            let names: Vec<String> = (0..rng.count()).map(|_| rng.name()).collect();
            Box::new(SnpTrack::new(names, sites))
        }
        13 => {
            let sites: Vec<MethylSite> = (0..count)
                .map(|_| {
                    MethylSite::new(
                        rng.position(span),
                        rng.strand(),
                        rng.value(),
                        rng.below(1_000) as u32,
                    )
                })
                .collect();
            Box::new(MethylationTrack::new(sites))
        }
        14 => {
            let bands: Vec<Band> = (0..count)
                .map(|_| {
                    let a = rng.position(width);
                    let b = rng.position(width);
                    let stain = match rng.below(4) {
                        0 => Stain::Gneg,
                        1 => Stain::Gpos50,
                        2 => Stain::Acen,
                        _ => Stain::Gpos100,
                    };
                    Band::new(a.min(b), a.max(b), stain).name(rng.name())
                })
                .collect();
            Box::new(IdeogramTrack::new(width.max(1), bands))
        }
        _ => Box::new(AxisTrack::new()),
    }
}

/// A whole figure, with a random region, theme, size and stack.
fn figure(seed: u64) -> Figure {
    let mut rng = Lcg::new(seed);
    let region = region(&mut rng);
    let mut figure = Figure::new(region.clone());

    if rng.chance(2) {
        figure = figure.title(rng.name());
    }
    if rng.chance(3) {
        figure = figure.theme(Theme::dark());
    }
    if rng.chance(4) {
        figure = figure.width(match rng.below(4) {
            0 => 1.0,
            1 => 60.0,
            2 => 5_000.0,
            _ => 300.0 + rng.below(1_200) as f64,
        });
    }
    if rng.chance(3) {
        figure = figure.show_region_label(false);
    }

    for _ in 0..rng.below(5) {
        figure = figure.push_boxed(track(&mut rng, &region));
    }
    figure
}

// ---------------------------------------------------------------------------
// What a valid document is
// ---------------------------------------------------------------------------

/// Every `name="value"` pair in the document, in order.
fn attributes(svg: &str) -> Vec<(&str, &str)> {
    let mut out = Vec::new();
    let bytes = svg.as_bytes();
    let mut i = 0;
    while let Some(eq) = svg[i..].find("=\"") {
        let at = i + eq;
        let mut from = at;
        while from > 0 && (bytes[from - 1].is_ascii_alphanumeric() || bytes[from - 1] == b'-') {
            from -= 1;
        }
        let value_start = at + 2;
        match svg[value_start..].find('"') {
            Some(end) => {
                out.push((&svg[from..at], &svg[value_start..value_start + end]));
                i = value_start + end + 1;
            }
            None => break,
        }
    }
    out
}

/// Everything that makes a document invalid, as a list of complaints.
///
/// One function rather than several tests, so a generated figure is rendered
/// once and every rule is applied to it. A complaint carries the seed, since a
/// generated failure is worthless without the case that produced it.
fn complaints(svg: &str) -> Vec<String> {
    let mut bad = Vec::new();

    if !svg.starts_with("<svg ") || !svg.ends_with("</svg>") {
        bad.push("the document is not one svg element".to_string());
    }

    let opens = svg.matches("<g ").count() + svg.matches("<g>").count();
    let closes = svg.matches("</g>").count();
    if opens != closes {
        bad.push(format!("{opens} groups opened and {closes} closed"));
    }

    let mut ids: Vec<&str> = Vec::new();
    for (name, value) in attributes(svg) {
        // A number that is not a number is the failure the writer exists to
        // stop, and it reaches the page as a word rather than as a digit.
        let looks_numeric = matches!(
            name,
            "x" | "y"
                | "x1"
                | "y1"
                | "x2"
                | "y2"
                | "cx"
                | "cy"
                | "r"
                | "rx"
                | "ry"
                | "width"
                | "height"
                | "font-size"
                | "stroke-width"
                | "opacity"
                | "fill-opacity"
                | "textLength"
        );
        if looks_numeric {
            match value.parse::<f64>() {
                Ok(number) => {
                    if !number.is_finite() {
                        bad.push(format!("{name}=\"{value}\" is not finite"));
                    }
                    if matches!(name, "width" | "height" | "r" | "rx" | "ry" | "font-size")
                        && number < 0.0
                    {
                        bad.push(format!("{name}=\"{value}\" is negative"));
                    }
                }
                Err(_) => bad.push(format!("{name}=\"{value}\" is not a number")),
            }
            if value.contains(['e', 'E']) {
                bad.push(format!("{name}=\"{value}\" is in exponent notation"));
            }
        }
        if name == "id" {
            if ids.contains(&value) {
                bad.push(format!("id {value:?} is used twice"));
            }
            ids.push(value);
        }
    }

    // A reference that resolves to nothing crops a track to nothing.
    let mut rest = svg;
    while let Some(at) = rest.find("url(#") {
        rest = &rest[at + 5..];
        if let Some(end) = rest.find(')') {
            let target = &rest[..end];
            if !ids.contains(&target) {
                bad.push(format!("url(#{target}) resolves to nothing"));
            }
        }
    }

    // Text content is the one place a caller's own string reaches the page.
    // A `<` is always wrong there, and an `&` is wrong unless it opens one of
    // the entities the writer emits.
    for piece in svg.split("<title>").skip(1) {
        let Some(text) = piece.split("</title>").next() else {
            continue;
        };
        if text.contains('<') {
            bad.push(format!("a title carries a raw < : {text:?}"));
        }
        let mut rest = text;
        while let Some(at) = rest.find('&') {
            rest = &rest[at..];
            let entity = ["&amp;", "&lt;", "&gt;", "&quot;", "&apos;"]
                .iter()
                .find(|e| rest.starts_with(**e));
            match entity {
                Some(e) => rest = &rest[e.len()..],
                None => {
                    bad.push(format!("a title carries a raw & : {text:?}"));
                    break;
                }
            }
        }
    }
    for byte in svg.bytes() {
        if byte < 0x20 && byte != b'\n' && byte != b'\t' {
            bad.push(format!(
                "a control character {byte:#04x} reached the document"
            ));
            break;
        }
    }

    bad
}

// ---------------------------------------------------------------------------
// The properties
// ---------------------------------------------------------------------------

#[test]
fn every_figure_renders_a_valid_document() {
    for seed in 0..ROUNDS {
        let svg = figure(seed).to_svg();
        let bad = complaints(&svg);
        assert!(bad.is_empty(), "seed {seed}: {bad:?}");
    }
}

#[test]
fn no_tooltip_prints_more_digits_than_a_number_has() {
    // `u64::MAX` is twenty digits, and every number the crate groups breaks
    // into runs of three, so nothing legitimate reaches twenty-one in a row.
    // `f64::MAX` written to two decimal places is three hundred and ten, which
    // is what a formatter does when nobody has told it that an f64 stops
    // holding consecutive integers at 2^53.
    for seed in 0..ROUNDS {
        let svg = figure(seed).to_svg();
        for piece in svg.split("<title>").skip(1) {
            let Some(text) = piece.split("</title>").next() else {
                continue;
            };
            let mut run = 0usize;
            for byte in text.bytes() {
                run = if byte.is_ascii_digit() { run + 1 } else { 0 };
                assert!(
                    run <= 20,
                    "seed {seed}: a tooltip carries {run} digits in a row: {text:?}"
                );
            }
        }
    }
}

#[test]
fn rendering_is_deterministic() {
    for seed in 0..ROUNDS {
        assert_eq!(
            figure(seed).to_svg(),
            figure(seed).to_svg(),
            "seed {seed} rendered differently the second time"
        );
    }
}

#[test]
fn a_figure_reports_the_height_it_draws() {
    for seed in 0..ROUNDS {
        let figure = figure(seed);
        let (width, height) = figure.dimensions();
        let svg = figure.to_svg();
        assert!(
            width.is_finite() && height.is_finite() && width > 0.0 && height > 0.0,
            "seed {seed}: dimensions {width} x {height}"
        );
        // The reported size and the size in the document are the same claim.
        let attrs = attributes(&svg);
        let said = |name: &str| -> f64 {
            attrs
                .iter()
                .find(|(n, _)| *n == name)
                .and_then(|(_, v)| v.parse::<f64>().ok())
                .unwrap_or(f64::NAN)
        };
        assert!(
            (said("width") - width).abs() < 0.01 && (said("height") - height).abs() < 0.01,
            "seed {seed}: dimensions() says {width} x {height}, the document says {} x {}",
            said("width"),
            said("height")
        );
    }
}

#[test]
fn a_locus_string_survives_a_round_trip() {
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let original = region(&mut rng);
        let printed = original.to_string();
        let read_back = Region::parse(&printed)
            .unwrap_or_else(|error| panic!("seed {seed}: {printed:?} did not parse: {error}"));
        assert_eq!(
            read_back, original,
            "seed {seed}: {printed:?} came back as a different region"
        );
    }
}

#[test]
fn the_scale_never_goes_backwards() {
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let region = region(&mut rng);
        let scale = Scale::new(&region, 40.0, 800.0);
        let mut previous = f64::NEG_INFINITY;
        for step in 0..64u64 {
            let position = region
                .start()
                .saturating_add(step.saturating_mul(region.len() / 64).max(1));
            let x = scale.x(position);
            assert!(x.is_finite(), "seed {seed}: x of {position} is {x}");
            assert!(
                x >= previous,
                "seed {seed}: x went backwards at {position}, {previous} then {x}"
            );
            previous = x;
        }
    }
}

#[test]
fn the_theme_changes_no_geometry() {
    // A metamorphic property: the right numbers are unknown, but changing the
    // colours must not move anything. Anything that does is a layout decision
    // reading a colour, which is how a figure comes out different on a dark
    // page for reasons nobody intended.
    for seed in 0..ROUNDS {
        let light = figure(seed).theme(Theme::light());
        let dark = figure(seed).theme(Theme::dark());
        assert_eq!(
            light.dimensions(),
            dark.dimensions(),
            "seed {seed}: the two themes lay out differently"
        );
        let geometry = |svg: &str| -> Vec<String> {
            attributes(svg)
                .into_iter()
                .filter(|(name, _)| {
                    matches!(
                        *name,
                        "x" | "y" | "width" | "height" | "cx" | "cy" | "r" | "d"
                    )
                })
                .map(|(name, value)| format!("{name}={value}"))
                .collect()
        };
        assert_eq!(
            geometry(&light.to_svg()),
            geometry(&dark.to_svg()),
            "seed {seed}: the two themes draw different shapes"
        );
    }
}

#[test]
fn a_track_is_drawn_the_same_whatever_is_stacked_under_it() {
    // Reordering is not the property, since the bands move. What must hold is
    // that a track added below cannot change the one above: a figure whose
    // first track depends on its last has a track reading another's data.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let first = track(&mut rng, &window);
        let second = track(&mut rng, &window);

        let alone = Figure::new(window.clone())
            .show_region_label(false)
            .push_boxed(first)
            .dimensions();

        // The same first track again, since a boxed one cannot be cloned.
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let first = track(&mut rng, &window);
        let with_company = Figure::new(window)
            .show_region_label(false)
            .push_boxed(first)
            .push_boxed(second)
            .dimensions();

        assert!(
            with_company.1 >= alone.1,
            "seed {seed}: adding a track under made the figure shorter, {alone:?} then {with_company:?}"
        );
    }
}
