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

use std::collections::BTreeMap;

use karyon::{
    read, AnnotationValue, Association, AxisTrack, Band, CigarOp, CodonTrack, CoverageTrack,
    Feature, FeatureTrack, Figure, Format, GeoFlow, GeoLocation, GeoProjection, IdeogramTrack,
    ManhattanTrack, Map, MatrixRow, MatrixTrack, MethylSite, MethylationTrack, MsaSequence,
    MsaTrack, OrfTrack, PhyloConnector, PhyloMap, PileupTrack, Read, Region, RenderProfile, Scale,
    SequenceTrack, SnpSite, SnpTrack, Stain, Strand, StructuralTrack, StructuralVariant, SvKind,
    Theme, TimeDirection, Track, Tree, TreeShape, Variant, VariantTrack, Window, WindowTrack,
};

/// How many figures each property is given before it is believed.
const ROUNDS: u64 = 10_000;

/// How many drawings each map property is given.
///
/// Fewer, because a map carries the world outline and costs more to render
/// than a stack of tracks does. Enough that the hostile coordinates all come
/// up: every projection, every kind of unplaceable latitude.
const MAP_ROUNDS: u64 = 600;

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

/// The longest run of digits in any tooltip in the document.
fn longest_digit_run(svg: &str) -> (usize, String) {
    let mut worst = (0usize, String::new());
    for piece in svg.split("<title>").skip(1) {
        let Some(text) = piece.split("</title>").next() else {
            continue;
        };
        let mut run = 0usize;
        for byte in text.bytes() {
            run = if byte.is_ascii_digit() { run + 1 } else { 0 };
            if run > worst.0 {
                worst = (run, text.to_string());
            }
        }
    }
    worst
}

#[test]
fn no_tooltip_prints_more_digits_than_a_number_has() {
    // `u64::MAX` is twenty digits, and every number the crate groups breaks
    // into runs of three, so nothing legitimate reaches twenty-one in a row.
    // `f64::MAX` written out in full is three hundred and nine, which is what
    // a formatter does when nobody has told it that an f64 stops holding
    // consecutive integers at 2^53.
    //
    // Both kinds of drawing, because they have separate writers over separate
    // arithmetic and the second one had this wrong while the first was fixed.
    for seed in 0..ROUNDS {
        let (run, text) = longest_digit_run(&figure(seed).to_svg());
        assert!(
            run <= 20,
            "seed {seed}: a tooltip carries {run} digits in a row: {text:?}"
        );
    }
    for seed in 0..MAP_ROUNDS {
        let (run, text) = longest_digit_run(&map(seed).to_svg());
        assert!(
            run <= 20,
            "map seed {seed}: a tooltip carries {run} digits in a row: {text:?}"
        );
        let (run, text) = longest_digit_run(&phylo_map(seed).to_svg());
        assert!(
            run <= 20,
            "phylo seed {seed}: a tooltip carries {run} digits in a row: {text:?}"
        );
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

// ---------------------------------------------------------------------------
// Trees
// ---------------------------------------------------------------------------

/// A branch length, tame when the property needs arithmetic to mean something.
///
/// The hostile ones are not decoration. A negative length is what a rate
/// smoothing program writes when it overshoots, and a zero length is what a
/// polytomy resolved into a bifurcation looks like.
fn branch_length(rng: &mut Lcg, tame: bool) -> String {
    if tame {
        match rng.below(8) {
            0 => "0".to_string(),
            1 => "0.0".to_string(),
            2 => "0.000000001".to_string(),
            3 => format!("{}", rng.below(1_000_000)),
            _ => format!("{}", rng.below(100_000) as f64 / 1_000.0),
        }
    } else {
        // A length the parser rejects throws the whole tree away, so the two
        // that do are rare on purpose. Left common, they cost the layout
        // properties most of their cases and nothing says so.
        match rng.below(60) {
            0 => "NaN".to_string(),
            1 => "inf".to_string(),
            2..=10 => "-1".to_string(),
            11..=19 => "1e308".to_string(),
            20..=28 => "-0.0".to_string(),
            29..=37 => "1e-308".to_string(),
            38..=46 => "0".to_string(),
            _ => format!("{}", rng.below(100_000) as f64 / 1_000.0),
        }
    }
}

/// A Newick string with a random topology and tip names that are all different.
///
/// All different because every property below keys a tip by its name. Rerooting
/// renumbers the nodes, so an index means nothing across it, and a duplicate
/// name would make two tips indistinguishable in exactly the comparison the
/// property exists to make.
fn newick(rng: &mut Lcg, tame: bool) -> String {
    let mut tips = 0usize;
    format!("{};", newick_clade(rng, tame, 0, &mut tips))
}

fn newick_clade(rng: &mut Lcg, tame: bool, depth: usize, tips: &mut usize) -> String {
    let leaf = depth > 0 && (depth >= 4 || rng.chance(3));
    let mut out = if leaf {
        *tips += 1;
        format!("t{tips}")
    } else {
        let children = 2 + rng.below(2) as usize;
        let parts: Vec<String> = (0..children)
            .map(|_| newick_clade(rng, tame, depth + 1, tips))
            .collect();
        let mut inner = format!("({})", parts.join(","));
        // The format writes support values and internal names in the same
        // place, so both readings have to be generated.
        match rng.below(6) {
            0 => inner.push_str(&format!("{}", rng.below(101))),
            1 => inner.push_str(&format!("node{}", rng.below(100))),
            _ => {}
        }
        inner
    };
    // The root has no incoming branch. A tame tree gives every other node one,
    // because a missing length is an undefined distance.
    if depth > 0 && (tame || !rng.chance(4)) {
        out.push_str(&format!(":{}", branch_length(rng, tame)));
    }
    out
}

/// Every tip-to-tip distance, worked out from the undirected edges alone.
///
/// Rerooting changes which end of an edge is the parent, and nothing else, so
/// this is the number that has to survive it. It is computed here by walking
/// the graph rather than by asking the tree, because a helper shared with the
/// code under test would agree with a mistake in it. `None` when some branch
/// has no length, since then there is no distance to preserve.
fn tip_distances(tree: &Tree) -> Option<BTreeMap<(String, String), f64>> {
    let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); tree.nodes().len()];
    for (child, clade) in tree.nodes().iter().enumerate() {
        let Some(parent) = clade.parent else {
            continue;
        };
        let length = clade.branch_length?;
        if !length.is_finite() || length < 0.0 {
            return None;
        }
        adjacency[parent].push((child, length));
        adjacency[child].push((parent, length));
    }

    let mut out = BTreeMap::new();
    for start in tree.leaves() {
        let from = tree.nodes()[start].name.clone()?;
        let mut distance = vec![0.0f64; adjacency.len()];
        let mut seen = vec![false; adjacency.len()];
        seen[start] = true;
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            for (next, length) in &adjacency[node] {
                if !seen[*next] {
                    seen[*next] = true;
                    distance[*next] = distance[node] + length;
                    stack.push(*next);
                }
            }
        }
        for end in tree.leaves() {
            let to = tree.nodes()[end].name.clone()?;
            if from < to {
                out.insert((from.clone(), to), distance[end]);
            }
        }
    }
    Some(out)
}

/// Whether two distances agree, allowing for the arithmetic that got them here.
///
/// Rooting on an edge splits it in two, and the halves do not have to add back
/// up to the original bit for bit.
fn same_distance(before: f64, after: f64) -> bool {
    if before == after {
        return true;
    }
    (before - after).abs() <= 1e-9 * before.abs().max(after.abs()).max(1.0)
}

// ---------------------------------------------------------------------------
// Maps
// ---------------------------------------------------------------------------

/// A latitude and longitude, on Earth most of the time and off it the rest.
///
/// The exact poles and antimeridian are generated deliberately: they are the
/// values a projection divides by zero at, and the ones a bounds check written
/// with `<` instead of `<=` rejects.
/// Both coordinates have to be ordinary for a location to be drawn at all, so
/// the odds here are per coordinate and the drawable share is their square.
/// Hostile more often than this and the properties about drawing something
/// stop having much of anything to draw.
fn geo_position(rng: &mut Lcg) -> (f64, f64) {
    let latitude = match rng.below(24) {
        0 => f64::NAN,
        1 => f64::INFINITY,
        2 => f64::NEG_INFINITY,
        3 => 90.0,
        4 => -90.0,
        5 => 90.000_001,
        6 => -1e9,
        _ => (rng.below(18_001) as f64) / 100.0 - 90.0,
    };
    let longitude = match rng.below(24) {
        0 => f64::NAN,
        1 => f64::NEG_INFINITY,
        2 => f64::INFINITY,
        3 => 180.0,
        4 => -180.0,
        5 => 180.000_001,
        6 => 1e9,
        _ => (rng.below(36_001) as f64) / 100.0 - 180.0,
    };
    (latitude, longitude)
}

/// A place name, never one that could be mistaken for part of a tooltip.
///
/// The properties below count locations by their tooltips, so a name is
/// allowed to be hostile to the escaping but not to the counting.
fn place(rng: &mut Lcg) -> String {
    match rng.below(8) {
        0 => String::new(),
        1 => "a<b&c\"d'e".to_string(),
        2 => "\u{1f30d} \u{202e}".to_string(),
        3 => "x".repeat(120),
        _ => format!("p{}", rng.below(30)),
    }
}

/// A whole map, with a random projection, locations and flows between them.
fn map(seed: u64) -> Map {
    let mut rng = Lcg::new(seed);
    let mut drawing = Map::new();

    drawing = drawing.projection(match rng.below(5) {
        0 => GeoProjection::Equirectangular,
        1 => GeoProjection::Mercator,
        2 => GeoProjection::orthographic(0.0, 0.0),
        3 => GeoProjection::orthographic(90.0, 180.0),
        _ => {
            let (latitude, longitude) = geo_position(&mut rng);
            GeoProjection::orthographic(latitude, longitude)
        }
    });

    if rng.chance(3) {
        drawing = drawing.width(match rng.below(4) {
            0 => 1.0,
            1 => 40.0,
            2 => 4_000.0,
            _ => 300.0 + rng.below(900) as f64,
        });
    }
    if rng.chance(3) {
        drawing = drawing.height(60.0 + rng.below(900) as f64);
    }
    if rng.chance(4) {
        drawing = drawing.margin(rng.below(200) as f64);
    }
    if rng.chance(3) {
        drawing = drawing.theme(Theme::dark());
    }
    if rng.chance(3) {
        drawing = drawing.profile(match rng.below(5) {
            0 => RenderProfile::Compact,
            1 => RenderProfile::Presentation,
            2 => RenderProfile::Web,
            3 => RenderProfile::Dark,
            _ => RenderProfile::Manuscript,
        });
    }
    if rng.chance(3) {
        drawing = drawing.title(place(&mut rng));
    }
    if rng.chance(4) {
        drawing = drawing.subtitle(place(&mut rng));
    }
    drawing = drawing
        .show_graticule(!rng.chance(3))
        .show_labels(!rng.chance(3))
        .show_legend(!rng.chance(3));

    let mut names = Vec::new();
    for _ in 0..rng.count() {
        let (latitude, longitude) = geo_position(&mut rng);
        let name = place(&mut rng);
        let mut location = GeoLocation::new(name.clone(), latitude, longitude);
        if rng.chance(2) {
            location = location.category(place(&mut rng));
        }
        if rng.chance(2) {
            location = location.value(rng.value());
        }
        if rng.chance(3) {
            location = location.count(rng.below(1_000));
        }
        names.push(name);
        drawing = drawing.push(location);
    }

    for _ in 0..rng.below(6) {
        // An endpoint that names nothing, and one that names two places at
        // once, are both things a real table of flows contains.
        let pick = |rng: &mut Lcg, names: &[String]| -> String {
            if names.is_empty() || rng.chance(4) {
                place(rng)
            } else {
                names[rng.below(names.len() as u64) as usize].clone()
            }
        };
        let from = pick(&mut rng, &names);
        let to = pick(&mut rng, &names);
        let mut flow = GeoFlow::new(from, to).weight(rng.value());
        if rng.chance(2) {
            flow = flow.undirected();
        }
        if rng.chance(3) {
            flow = flow.category(place(&mut rng));
        }
        drawing = drawing.push_flow(flow);
    }

    drawing
}

/// A whole phylogeographic composition: a tree, an annotation naming where
/// each tip came from, and a table of coordinates for those names.
///
/// The two drawings share a module and a writer but not a code path, and a
/// tooltip bug lived in one of them while the other was clean.
fn phylo_map(seed: u64) -> PhyloMap {
    let mut rng = Lcg::new(seed);
    let tame = rng.chance(2);
    let text = newick(&mut rng, tame);
    let mut tree = Tree::parse_newick(&text).unwrap_or_else(|_| {
        Tree::parse_newick("((a:1,b:1):1,(c:1,d:1):1);").expect("a fixed tree parses")
    });

    // Where each tip came from, and how far the table of coordinates agrees.
    // A tip whose place is not in the table is the ordinary case, not an edge
    // one: a sample sheet is never complete.
    let places: Vec<String> = (0..1 + rng.below(5)).map(|_| place(&mut rng)).collect();
    for leaf in tree.leaves() {
        if rng.chance(6) {
            continue;
        }
        let where_from = places[rng.below(places.len() as u64) as usize].clone();
        if let Some(annotations) = tree.annotations_mut(leaf) {
            annotations.insert("location".to_string(), AnnotationValue::Text(where_from));
            if rng.chance(2) {
                annotations.insert("date".to_string(), AnnotationValue::Number(rng.value()));
            }
        }
    }

    let mut coordinates: Vec<GeoLocation> = Vec::new();
    for name in &places {
        // A place with no row in the coordinate table is the ordinary case.
        if rng.chance(4) {
            continue;
        }
        let (latitude, longitude) = geo_position(&mut rng);
        let mut location = GeoLocation::new(name.clone(), latitude, longitude);
        if rng.chance(2) {
            location = location.value(rng.value());
        }
        if rng.chance(3) {
            location = location.count(rng.below(500));
        }
        coordinates.push(location);
    }

    let mut drawing = PhyloMap::new(tree)
        .location_by("location")
        .coordinates(coordinates)
        .connector(match rng.below(3) {
            0 => PhyloConnector::Aggregated,
            1 => PhyloConnector::Individual,
            _ => PhyloConnector::None,
        })
        .shape(if rng.chance(2) {
            TreeShape::Phylogram
        } else {
            TreeShape::Cladogram
        })
        .projection(match rng.below(3) {
            0 => GeoProjection::Equirectangular,
            1 => GeoProjection::Mercator,
            _ => {
                let (latitude, longitude) = geo_position(&mut rng);
                GeoProjection::orthographic(latitude, longitude)
            }
        });

    if rng.chance(2) {
        drawing = drawing.time("date").time_direction(if rng.chance(2) {
            TimeDirection::Increasing
        } else {
            TimeDirection::Decreasing
        });
    }
    if rng.chance(3) {
        drawing = drawing.diameter(match rng.below(4) {
            0 => 1.0,
            1 => 50.0,
            2 => 3_000.0,
            _ => 300.0 + rng.below(700) as f64,
        });
    }
    if rng.chance(3) {
        drawing = drawing.margin(rng.below(150) as f64);
    }
    if rng.chance(3) {
        drawing = drawing.radial_start(rng.below(1_000) as f64 - 500.0);
    }
    if rng.chance(3) {
        drawing = drawing.radial_sweep(rng.below(800) as f64 - 100.0);
    }
    if rng.chance(3) {
        drawing = drawing.theme(Theme::dark());
    }
    if rng.chance(3) {
        drawing = drawing.title(place(&mut rng));
    }
    drawing
}

/// Every number the document positions something with, path data included.
///
/// The attribute walk on its own would miss the projection entirely: the land
/// outlines are a `d` string, so that is where a latitude that came out as ten
/// to the fifteenth would land.
fn geometry_numbers(svg: &str) -> Vec<f64> {
    let mut out = Vec::new();
    for (name, value) in attributes(svg) {
        if matches!(
            name,
            "x" | "y" | "x1" | "y1" | "x2" | "y2" | "cx" | "cy" | "r" | "width" | "height"
        ) {
            if let Ok(number) = value.parse::<f64>() {
                out.push(number);
            }
        }
        if name == "d" || name == "points" {
            for piece in value.split([' ', ',']) {
                if let Ok(number) = piece.parse::<f64>() {
                    out.push(number);
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The properties: trees
// ---------------------------------------------------------------------------

#[test]
fn a_layout_places_every_node_exactly_once() {
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(2);
        let text = newick(&mut rng, tame);
        let Ok(tree) = Tree::parse_newick(&text) else {
            continue;
        };
        for cladogram in [false, true] {
            let places = tree.layout(cladogram);
            assert_eq!(
                places.len(),
                tree.nodes().len(),
                "seed {seed}: {} nodes placed of {}",
                places.len(),
                tree.nodes().len()
            );
            let mut seen = vec![false; tree.nodes().len()];
            for place in &places {
                assert!(
                    !std::mem::replace(&mut seen[place.node], true),
                    "seed {seed}: node {} placed twice",
                    place.node
                );
                assert!(
                    place.row.is_finite(),
                    "seed {seed}: node {} sits at row {}",
                    place.node,
                    place.row
                );
            }
            // A leaf takes its own row, so the tips fill nought upwards with
            // nothing shared and nothing skipped. Two tips on one row is two
            // samples drawn on top of each other.
            let mut rows: Vec<f64> = tree
                .leaves()
                .into_iter()
                .map(|leaf| places[leaf].row)
                .collect();
            rows.sort_by(|a, b| a.partial_cmp(b).expect("leaf rows are finite"));
            for (expected, row) in rows.iter().enumerate() {
                assert_eq!(
                    *row, expected as f64,
                    "seed {seed}: the tips do not fill the rows, {rows:?}"
                );
            }
        }
    }
}

#[test]
fn a_parent_sits_between_its_children() {
    // The module says so, and a panel sorted by the leaf order lines up only
    // if it holds. A parent outside its children's rows draws a branch that
    // crosses the clade it belongs to.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(2);
        let text = newick(&mut rng, tame);
        let Ok(tree) = Tree::parse_newick(&text) else {
            continue;
        };
        for cladogram in [false, true] {
            let places = tree.layout(cladogram);
            for (index, clade) in tree.nodes().iter().enumerate() {
                if clade.children.is_empty() {
                    continue;
                }
                let rows: Vec<f64> = clade
                    .children
                    .iter()
                    .map(|child| places[*child].row)
                    .collect();
                let low = rows.iter().copied().fold(f64::INFINITY, f64::min);
                let high = rows.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let row = places[index].row;
                assert!(
                    row >= low && row <= high,
                    "seed {seed}: node {index} sits at {row}, outside its children {rows:?}"
                );
            }
        }
    }
}

#[test]
fn counting_branches_and_adding_them_up_agree_on_the_rows() {
    // A metamorphic property. What the depths should be is a matter of which
    // measure was asked for, but the rows are not: the module promises they do
    // not move between a phylogram and a cladogram, which is what lets a panel
    // beside the tree be sorted once and stay aligned in both.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(2);
        let text = newick(&mut rng, tame);
        let Ok(tree) = Tree::parse_newick(&text) else {
            continue;
        };
        let lengths = tree.layout(false);
        let branches = tree.layout(true);
        for (with, without) in lengths.iter().zip(branches.iter()) {
            assert_eq!(with.node, without.node, "seed {seed}: the node order moved");
            assert_eq!(
                with.row, without.row,
                "seed {seed}: node {} sits at row {} counted and {} added up",
                with.node, without.row, with.row
            );
        }
    }
}

#[test]
fn a_cladogram_depth_is_the_number_of_branches() {
    // Countable without the tree's help, which is the point: the oracle here
    // is the definition, not another call into the code being checked.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(2);
        let text = newick(&mut rng, tame);
        let Ok(tree) = Tree::parse_newick(&text) else {
            continue;
        };
        let places = tree.layout(true);
        for (index, clade) in tree.nodes().iter().enumerate() {
            let expected = tree.ancestors(index).len() as f64;
            assert_eq!(
                places[index].depth, expected,
                "seed {seed}: node {index} is {} branches from the root and drawn at {}",
                expected, places[index].depth
            );
            if let Some(parent) = clade.parent {
                assert_eq!(
                    places[index].depth,
                    places[parent].depth + 1.0,
                    "seed {seed}: node {index} is not one branch past its parent"
                );
            }
        }
    }
}

/// One way of asking for a new root, applied and reporting whether it took.
type Reroot = Box<dyn Fn(&mut Tree) -> bool>;

#[test]
fn rerooting_moves_the_root_and_nothing_else() {
    // The one thing rerooting must not do. Every tip-to-tip distance is a
    // property of the undirected edges, and reorienting them cannot change a
    // distance, so a changed one means a branch length was dropped, doubled or
    // left behind on the old root. The tip set must survive too: rooting on a
    // leaf would quietly turn a sample into an internal node.
    let mut checked = 0u64;
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let text = newick(&mut rng, true);
        let Ok(source) = Tree::parse_newick(&text) else {
            continue;
        };
        let Some(before) = tip_distances(&source) else {
            continue;
        };
        if source.leaf_count() < 3 {
            continue;
        }
        let mut names_before = source.leaf_names();
        names_before.sort();

        // Every way in: an internal node, a monophyletic outgroup, and the
        // midpoint of the longest path.
        let internal: Vec<usize> = (0..source.nodes().len())
            .filter(|node| !source.nodes()[*node].children.is_empty())
            .collect();
        let mut ways: Vec<(&str, Reroot)> = vec![(
            "midpoint",
            Box::new(|tree: &mut Tree| tree.reroot_midpoint().is_some()),
        )];
        let pick = internal[rng.below(internal.len() as u64) as usize];
        ways.push((
            "internal",
            Box::new(move |tree: &mut Tree| tree.reroot(pick)),
        ));
        let leaf = source.leaves()[rng.below(source.leaf_count() as u64) as usize];
        ways.push((
            "outgroup",
            Box::new(move |tree: &mut Tree| tree.reroot_outgroup(&[leaf]).is_some()),
        ));

        for (how, reroot) in ways {
            let mut tree = source.clone();
            if !reroot(&mut tree) {
                continue;
            }
            checked += 1;

            let mut names_after = tree.leaf_names();
            names_after.sort();
            assert_eq!(
                names_before, names_after,
                "seed {seed}: rerooting on the {how} changed which tips there are"
            );
            assert!(
                tree.nodes()[tree.root()].parent.is_none(),
                "seed {seed}: the {how} root still has a parent"
            );

            let after = tip_distances(&tree)
                .unwrap_or_else(|| panic!("seed {seed}: the {how} root lost a branch length"));
            assert_eq!(
                before.len(),
                after.len(),
                "seed {seed}: the {how} root left {} tip pairs of {}",
                after.len(),
                before.len()
            );
            for (pair, distance) in &before {
                let now = after.get(pair).unwrap_or_else(|| {
                    panic!("seed {seed}: {pair:?} is gone after the {how} root")
                });
                assert!(
                    same_distance(*distance, *now),
                    "seed {seed}: the {how} root moved {pair:?} from {distance} to {now}"
                );
            }
        }
    }
    assert!(
        checked > 1_000,
        "only {checked} reroots actually happened, so this proves little"
    );
}

#[test]
fn a_refused_reroot_leaves_the_tree_alone() {
    // Returning false is a promise about what did not happen. A half-applied
    // rerooting is worse than a refused one, because nothing downstream is
    // told about it.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(3);
        let text = newick(&mut rng, tame);
        let Ok(source) = Tree::parse_newick(&text) else {
            continue;
        };
        let before = source.nodes().to_vec();

        let mut tree = source.clone();
        // A leaf, an index that is not a node, and an outgroup that is neither
        // monophyletic nor made only of tips.
        let leaf = source.leaves()[rng.below(source.leaf_count() as u64) as usize];
        if !tree.reroot(leaf) {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: a refused leaf root moved something"
            );
        }
        let mut tree = source.clone();
        if !tree.reroot(source.nodes().len() + 7) {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: a root on nothing moved something"
            );
        }
        let mut tree = source.clone();
        if tree.reroot_outgroup(&[]).is_none() {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: an empty outgroup moved something"
            );
        }
        let mut tree = source.clone();
        let everything: Vec<usize> = source.leaves();
        if tree.reroot_outgroup(&everything).is_none() {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: a whole-tree outgroup moved something"
            );
        }
        let mut tree = source.clone();
        let internal: Vec<usize> = (0..source.nodes().len())
            .filter(|node| !source.nodes()[*node].children.is_empty())
            .collect();
        if tree.reroot_outgroup(&internal).is_none() {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: an internal outgroup moved something"
            );
        }
        let mut tree = source.clone();
        if tree.reroot_midpoint().is_none() {
            assert_eq!(
                tree.nodes(),
                before,
                "seed {seed}: a refused midpoint moved something"
            );
        }
    }
}

#[test]
fn a_string_that_is_not_a_tree_comes_back_as_an_error() {
    // Not that it rejects them: that it decides. A parser that panics on a
    // truncated file takes the caller down with it, and every one of these is
    // something a real Newick file has been found to contain.
    let fixed = [
        "",
        ";",
        "(",
        ")",
        "()",
        "(,)",
        "(();",
        "((a,b),c",
        "a:b;",
        "(a:1,b:2):;",
        "[&R]",
        "(a,b);extra",
        "(:1,:2);",
        "((((((((((a))))))))));",
        "(a,(b,(c,(d,(e)))));",
    ];
    for text in fixed {
        let _ = Tree::parse_newick(text);
        let _ = Tree::parse_annotated_newick(text);
        let _ = Tree::parse_nexus(text);
    }
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let tame = rng.chance(2);
        let good = newick(&mut rng, tame);
        let bytes = good.as_bytes();
        let text = match rng.below(6) {
            0 => good[..rng.below(good.len() as u64 + 1) as usize].to_string(),
            1 => good.replace('(', ""),
            2 => good.replace(')', "("),
            3 => good.replace(':', ":x"),
            4 => {
                let at = rng.below(bytes.len() as u64) as usize;
                format!("{}{}{}", &good[..at], "\u{0}\u{1f9ec},();:", &good[at..])
            }
            _ => good.repeat(2),
        };
        // No assertion on the verdict, only on getting one. A panic here fails
        // the test on its own.
        let _ = Tree::parse_newick(&text);
        let _ = Tree::parse_annotated_newick(&text);
        let _ = Tree::parse_nexus(&text);
    }
}

// ---------------------------------------------------------------------------
// The properties: maps
// ---------------------------------------------------------------------------

#[test]
fn every_map_renders_a_valid_document() {
    // The same rules the figures are held to. A map is a separate writer over
    // separate arithmetic, so it can fail all of them separately: a latitude
    // that came through as NaN reaches the page as `cy="NaN"`.
    for seed in 0..MAP_ROUNDS {
        let svg = map(seed).to_svg();
        let bad = complaints(&svg);
        assert!(bad.is_empty(), "seed {seed}: {bad:?}");
    }
    for seed in 0..MAP_ROUNDS {
        let svg = phylo_map(seed).to_svg();
        let bad = complaints(&svg);
        assert!(bad.is_empty(), "phylo seed {seed}: {bad:?}");
    }
}

#[test]
fn a_map_renders_the_same_way_twice() {
    for seed in 0..MAP_ROUNDS {
        assert_eq!(
            map(seed).to_svg(),
            map(seed).to_svg(),
            "seed {seed} rendered differently the second time"
        );
        assert_eq!(
            phylo_map(seed).to_svg(),
            phylo_map(seed).to_svg(),
            "phylo seed {seed} rendered differently the second time"
        );
    }
}

#[test]
fn no_location_is_dropped_without_saying_so() {
    // The failure this exists for is silent: a sample given to the map, not
    // drawn, and not mentioned. Nobody counts the dots on a map against the
    // rows of their table. Every location is drawn, counted as invalid, or
    // named in the notice, and the notice says how many.
    for seed in 0..MAP_ROUNDS {
        let drawing = map(seed);
        let total = drawing.locations().len();
        let invalid = drawing.invalid_location_count();
        let svg = drawing.to_svg();

        // One tooltip per drawn location, and its shape is the contract: the
        // name, then the latitude, then the longitude.
        let drawn = svg.matches("; latitude ").count();
        assert!(
            drawn + invalid <= total,
            "seed {seed}: {drawn} drawn plus {invalid} invalid, of {total} given"
        );
        let hidden = total - drawn - invalid;
        if hidden > 0 {
            assert!(
                svg.contains("outside projection"),
                "seed {seed}: {hidden} of {total} locations vanished without a word"
            );
            assert!(
                svg.contains(&format!("{hidden} location")),
                "seed {seed}: {hidden} locations went missing and the notice does not say {hidden}"
            );
        }
    }
}

#[test]
fn a_map_keeps_the_coordinates_it_was_given() {
    // Clamping an out-of-range coordinate would put a sample somewhere it was
    // never sampled, and the map would look right. The module's answer is to
    // keep the number and refuse to draw it, which only holds if the number
    // comes back unchanged.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let (latitude, longitude) = geo_position(&mut rng);
        let location = GeoLocation::new(place(&mut rng), latitude, longitude);
        let position = location.position();
        let same = |given: f64, kept: f64| given == kept || (given.is_nan() && kept.is_nan());
        assert!(
            same(latitude, position.latitude()) && same(longitude, position.longitude()),
            "seed {seed}: {latitude},{longitude} came back as {},{}",
            position.latitude(),
            position.longitude()
        );
        // And a coordinate off the Earth is not drawable, whatever else it is.
        let plausible = (-90.0..=90.0).contains(&latitude)
            && (-180.0..=180.0).contains(&longitude)
            && latitude.is_finite()
            && longitude.is_finite();
        assert_eq!(
            position.is_valid(),
            plausible,
            "seed {seed}: {latitude},{longitude} is called {}",
            position.is_valid()
        );
    }
}

#[test]
fn a_projection_stays_on_the_page() {
    // A finite number can still be nonsense. Mercator sends a latitude of
    // ninety to infinity, and the near miss is worse than the hit: a tangent
    // taken at 89.999 degrees is finite, passes every check for a number, and
    // puts a coastline several million pixels off the canvas, where it drags
    // the bounding box of the file with it.
    for seed in 0..MAP_ROUNDS {
        let drawing = map(seed);
        let (width, height) = drawing.dimensions();
        let svg = drawing.to_svg();
        let bound = 10.0 * width.max(height).max(1.0);
        for number in geometry_numbers(&svg) {
            assert!(
                number.abs() <= bound,
                "seed {seed}: something is drawn at {number} on a {width} by {height} map"
            );
        }
        let drawing = phylo_map(seed);
        let (width, height) = drawing.dimensions();
        let svg = drawing.to_svg();
        let bound = 10.0 * width.max(height).max(1.0);
        for number in geometry_numbers(&svg) {
            assert!(
                number.abs() <= bound,
                "phylo seed {seed}: something is drawn at {number} on a {width} by {height} map"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The properties: readers
// ---------------------------------------------------------------------------

#[test]
fn every_format_puts_the_same_interval_on_the_same_bases() {
    // The audit beside the readers aims every format at one known base. This
    // asks the same question without picking the base, which is the half that
    // a fixed fixture cannot answer: a reader can be right at position 100 and
    // wrong at position 0, or at the last base of the window, or wherever the
    // arithmetic happens to carry.
    //
    // Two claims. Each reader inverts its own specification, so a span written
    // out and read back is the span again. And every reader agrees with every
    // other, which needs no specification at all: if one of them is a base out,
    // the figure that stacks it on another is a base out, and nothing says so.
    let region = Region::new("chr1", 0, 200_000).expect("a fixed window");
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let start = match rng.below(8) {
            0 => 0,
            1 => 1,
            2 => 199_000,
            _ => rng.below(150_000),
        };
        let length = 1 + rng.below(40);
        let end = start + length;
        // 1-based inclusive is the same span counted from one, and the end
        // does not move because inclusive-of-the-end and half-open-past-it
        // land on the same number.
        let first = start + 1;

        // Whole spans, in every format that can carry one.
        let bed = format!("chr1\t{start}\t{end}\tg\t0\t+\n");
        let gff3 = format!("chr1\t.\tgene\t{first}\t{end}\t.\t+\t.\tID=g\n");
        let cytoband = format!("chr1\t{start}\t{end}\tp1\tgneg\n");
        let sam = format!(
            "r1\t0\tchr1\t{first}\t60\t{length}M\t*\t0\t0\t{}\t*\n",
            "A".repeat(length as usize)
        );

        let from_bed = read::interval::features(&bed, &region, None)
            .unwrap_or_else(|e| panic!("seed {seed}: bed {bed:?}: {e}"));
        let from_gff3 = read::interval::features(&gff3, &region, None)
            .unwrap_or_else(|e| panic!("seed {seed}: gff3 {gff3:?}: {e}"));
        let (_, bands) = read::interval::cytoband(&cytoband, "chr1")
            .unwrap_or_else(|e| panic!("seed {seed}: cytoband {cytoband:?}: {e}"));
        let reads = read::align::sam(&sam, &region)
            .unwrap_or_else(|e| panic!("seed {seed}: sam {sam:?}: {e}"));

        let spans = [
            ("bed", from_bed[0].start, from_bed[0].end),
            ("gff3", from_gff3[0].start, from_gff3[0].end),
            ("cytoband", bands[0].start, bands[0].end),
            ("sam", reads[0].start, reads[0].end()),
        ];
        for (format, read_start, read_end) in spans {
            assert_eq!(
                (read_start, read_end),
                (start, end),
                "seed {seed}: {format} put {start}..{end} at {read_start}..{read_end}"
            );
        }

        // Signals, where the same span is one row in one format and one row
        // per base in the other.
        let bedgraph = format!("chr1\t{start}\t{end}\t7\n");
        let depth: String = (first..=end).map(|at| format!("chr1\t{at}\t7\n")).collect();
        let covered = |pairs: Vec<(u64, f64)>| -> (u64, u64, usize) {
            let low = pairs.iter().map(|(at, _)| *at).min().unwrap_or(u64::MAX);
            let high = pairs.iter().map(|(at, _)| *at).max().unwrap_or(0);
            (low, high, pairs.len())
        };
        let from_bedgraph = covered(
            read::signal::dense(&bedgraph, &region, None)
                .unwrap_or_else(|e| panic!("seed {seed}: bedgraph: {e}")),
        );
        let from_depth = covered(
            read::signal::dense(&depth, &region, None)
                .unwrap_or_else(|e| panic!("seed {seed}: depth: {e}")),
        );
        assert_eq!(
            from_bedgraph, from_depth,
            "seed {seed}: bedGraph and samtools depth disagree about {start}..{end}"
        );
        assert_eq!(
            (from_bedgraph.0, from_bedgraph.1),
            (start, end - 1),
            "seed {seed}: a signal over {start}..{end} covers {}..={}",
            from_bedgraph.0,
            from_bedgraph.1
        );

        // Being told the format and being left to work it out are the same
        // answer. The module says format identification, not the arithmetic,
        // is where these readers went wrong in practice.
        for (text, format) in [(&bedgraph, Format::BedGraph), (&depth, Format::Depth)] {
            assert_eq!(
                read::signal::dense(text, &region, None).ok(),
                read::signal::dense(text, &region, Some(format)).ok(),
                "seed {seed}: {format:?} reads differently when it is named"
            );
        }
        assert_eq!(
            read::interval::features(&bed, &region, None).ok(),
            read::interval::features(&bed, &region, Some(Format::Bed)).ok(),
            "seed {seed}: bed reads differently when it is named"
        );
        assert_eq!(
            read::interval::features(&gff3, &region, None).ok(),
            read::interval::features(&gff3, &region, Some(Format::Gff3)).ok(),
            "seed {seed}: gff3 reads differently when it is named"
        );

        // Single points, where three more formats count from one.
        let vcf = format!("chr1\t{first}\t.\tC\tT\t.\t.\t.\n");
        let association = format!("{first}\t1e-9\n");
        let matrix = format!("sample\t{first}\nS1\t1\n");
        let calls = read::point::variants(&vcf, &region)
            .unwrap_or_else(|e| panic!("seed {seed}: vcf: {e}"));
        let points = read::point::associations(&association, &region)
            .unwrap_or_else(|e| panic!("seed {seed}: association: {e}"));
        let (sites, _) = read::table::matrix(&matrix, &region)
            .unwrap_or_else(|e| panic!("seed {seed}: matrix: {e}"));
        for (format, at) in [
            ("vcf", calls[0].pos),
            ("association", points[0].pos),
            ("matrix", sites[0]),
        ] {
            assert_eq!(
                at, start,
                "seed {seed}: {format} put 1-based {first} at 0-based {at}, not {start}"
            );
        }
    }
}
