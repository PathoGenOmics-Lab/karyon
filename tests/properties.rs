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

use karyon::track::traits::Traits;
use karyon::tree::Annotations;

use karyon::{
    read, AlignmentBlock, AncestralStateLayer, AnnotationValue, Association, AxisTrack, Band,
    BisulfiteTrack, BranchEventLayer, BranchGeometry, BranchIntervalLayer, BranchRateMixture,
    CigarOp, CladeBlock, CladeTrack, CodonTrack, CoverageTrack, DomainArchitecture, DomainFeature,
    DomainTrack, DotplotTrack, Feature, FeatureTrack, Figure, Format, Genome, GenomeTrack, GeoFlow,
    GeoLocation, GeoProjection, HomoplasyLayer, IdeogramTrack, Legend, LegendTrack, Locus,
    LocusTrack, LogoColumn, LogoTrack, ManhattanTrack, Map, MatrixRow, MatrixTrack, MethylSite,
    MethylationTrack, Molecule, MsaSequence, MsaTrack, OrfTrack, PhyloConnector, PhyloMap,
    PhylodynamicPoint, PhylodynamicScale, PhylodynamicTrack, PileupTrack, Read, Region,
    RenderProfile, Scale, SelectionEvidence, SelectionSite, SelectionTrack, SequenceTrack, SnpSite,
    SnpTrack, SplitRead, SplitReadTrack, SplitSegment, SquiggleTrack, Stain, Strand,
    StructuralTrack, StructuralVariant, SupportStyle, SurveillanceMetric, SurveillanceObservation,
    SurveillanceStyle, SurveillanceTrack, SvKind, SyntenyTrack, TanglegramTrack, Theme,
    TimeDirection, Track, TranscriptionUnit, TranscriptionUnitTrack, Tree, TreeProjection,
    TreeShape, TreeTrack, Variant, VariantTrack, Window, WindowTrack,
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

/// A sheet of metadata for `names`, sometimes for nobody at all.
///
/// The hostile cases are the point. A value that is not a number, a name the
/// rows have not got, a column every row leaves blank, a level per row so the
/// palette runs out, and names carrying the characters an SVG cannot hold: all
/// of those reach a tooltip and a swatch through the same two functions.
fn traits(rng: &mut Lcg, names: &[String]) -> Traits {
    let mut rows: BTreeMap<String, Annotations> = BTreeMap::new();
    // Half the time the sheet is about somebody else entirely, which is the
    // strip of nothing a caller who did not check the join would draw.
    let joined = rng.chance(2);
    for (index, name) in names.iter().enumerate() {
        let key = if joined {
            name.clone()
        } else {
            format!("{name}-elsewhere")
        };
        let mut held = Annotations::new();
        if !rng.chance(4) {
            held.insert(
                "lineage".to_string(),
                AnnotationValue::Text(format!("L{}", index % (1 + rng.below(9) as usize))),
            );
        }
        if !rng.chance(4) {
            held.insert("depth".to_string(), AnnotationValue::Number(rng.value()));
        }
        if !rng.chance(4) {
            held.insert("kept".to_string(), AnnotationValue::Boolean(rng.chance(2)));
        }
        if rng.chance(6) {
            held.insert("odd".to_string(), AnnotationValue::Text(rng.name()));
        }
        rows.insert(key, held);
    }

    let columns = ["lineage", "depth", "kept", "odd", "never"];
    let wanted: Vec<String> = columns
        .iter()
        .filter(|_| !rng.chance(3))
        .map(|key| key.to_string())
        .collect();
    let traits = Traits::new(rows).spread(wanted);
    if rng.chance(3) {
        traits.heading_room(0.0)
    } else {
        traits
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

    match rng.below(33) {
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
            let names: Vec<String> = rows.iter().map(|row| row.name.clone()).collect();
            let track = MatrixTrack::new(sites, rows);
            Box::new(if rng.chance(2) {
                track.traits(traits(rng, &names))
            } else {
                track
            })
        }
        7 => {
            let columns = rng.count();
            let rows: Vec<MsaSequence> = (0..rng.count())
                .map(|_| {
                    let residues: Vec<u8> = (0..columns).map(|_| rng.base()).collect();
                    MsaSequence::new(rng.name(), residues)
                })
                .collect();
            let names: Vec<String> = rows.iter().map(|row| row.name.clone()).collect();
            let track = MsaTrack::new(rows);
            Box::new(if rng.chance(2) {
                track.traits(traits(rng, &names))
            } else {
                track
            })
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
        15 => {
            let calls: Vec<Option<bool>> = (0..count)
                .map(|_| match rng.below(3) {
                    0 => None,
                    1 => Some(true),
                    _ => Some(false),
                })
                .collect();
            let sites: Vec<u64> = (0..calls.len()).map(|_| rng.position(span)).collect();
            let molecules: Vec<Molecule> = (0..rng.count())
                .map(|_| Molecule::new(rng.name(), calls.clone()))
                .collect();
            Box::new(BisulfiteTrack::new(sites, molecules))
        }
        16 => {
            let blocks: Vec<CladeBlock> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    let taxa: Vec<String> = (0..rng.count()).map(|_| rng.name()).collect();
                    CladeBlock::new(a.min(b), a.max(b), taxa).name(rng.name())
                })
                .collect();
            Box::new(CladeTrack::new(tree(rng), blocks))
        }
        17 => {
            let rows: Vec<DomainArchitecture> = (0..count)
                .map(|_| {
                    let mut row = DomainArchitecture::new(rng.name(), rng.position(span));
                    for _ in 0..rng.count() {
                        let a = rng.position(span);
                        let b = rng.position(span);
                        row = row.feature(DomainFeature::new(a.min(b), a.max(b)).label(rng.name()));
                    }
                    row
                })
                .collect();
            Box::new(DomainTrack::new(rows))
        }
        18 => {
            let sequences: Vec<(String, u64)> = (0..count)
                .map(|_| (rng.name(), rng.position(span)))
                .collect();
            Box::new(GenomeTrack::new(Genome::new(sequences)))
        }
        19 => {
            let mut legend = Legend::new();
            for _ in 0..count {
                legend = match rng.below(4) {
                    0 => legend.key(rng.name(), rng.name()),
                    1 => legend.dot(rng.name(), rng.name()),
                    2 => legend.line(rng.name(), rng.name()),
                    _ => legend.area(rng.name(), rng.name()),
                };
            }
            Box::new(LegendTrack::new(legend))
        }
        20 => {
            let loci: Vec<Locus> = (0..count)
                .map(|_| {
                    let genes: Vec<Feature> = (0..rng.count())
                        .map(|_| {
                            let a = rng.position(span);
                            let b = rng.position(span);
                            Feature::new(a.min(b), a.max(b)).strand(rng.strand())
                        })
                        .collect();
                    Locus::new(rng.name(), genes)
                })
                .collect();
            Box::new(LocusTrack::new(loci))
        }
        21 => {
            let columns: Vec<LogoColumn> = (0..count)
                .map(|_| LogoColumn::acgt(rng.value(), rng.value(), rng.value(), rng.value()))
                .collect();
            Box::new(LogoTrack::new(rng.position(span), columns))
        }
        22 => {
            let reads: Vec<SplitRead> = (0..count)
                .map(|_| {
                    let segments: Vec<SplitSegment> = (0..rng.count())
                        .map(|_| {
                            let a = rng.position(span);
                            let b = rng.position(span);
                            SplitSegment::new(a.min(b), a.max(b), rng.strand())
                        })
                        .collect();
                    SplitRead::new(segments)
                })
                .collect();
            Box::new(SplitReadTrack::new(reads))
        }
        23 => {
            let signal: Vec<f64> = (0..count).map(|_| rng.value()).collect();
            Box::new(SquiggleTrack::new(rng.count(), signal))
        }
        24 | 25 => {
            let blocks: Vec<AlignmentBlock> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    let c = rng.position(span);
                    let d = rng.position(span);
                    AlignmentBlock::new(a.min(b), a.max(b), c.min(d), c.max(d))
                })
                .collect();
            if rng.chance(2) {
                Box::new(SyntenyTrack::new(blocks))
            } else {
                Box::new(DotplotTrack::new(blocks))
            }
        }
        26 => Box::new(TanglegramTrack::new(tree(rng), tree(rng))),
        27 => {
            let units: Vec<TranscriptionUnit> = (0..count)
                .map(|_| {
                    let a = rng.position(span);
                    let b = rng.position(span);
                    TranscriptionUnit::new(a.min(b), a.max(b), rng.strand())
                })
                .collect();
            Box::new(TranscriptionUnitTrack::new(units))
        }
        28 => {
            // The builders, not just the constructor. Every one of these takes
            // a number nobody checks on the way in, and a figure has to survive
            // a caller who passes the wrong one.
            let mut track = TreeTrack::new(interacting_tree(rng));
            if rng.chance(2) {
                track = track.projection(match rng.below(3) {
                    0 => TreeProjection::Rectangular,
                    1 => TreeProjection::Circular,
                    _ => TreeProjection::Unrooted,
                });
            }
            if rng.chance(2) {
                track = track.shape(if rng.chance(2) {
                    TreeShape::Phylogram
                } else {
                    TreeShape::Cladogram
                });
            }
            if rng.chance(2) {
                track = track.row_height(rng.value());
            }
            if rng.chance(2) {
                track = track.fan(rng.value());
            }
            if rng.chance(2) {
                track = track.radial_start(rng.value()).radial_sweep(rng.value());
            }
            if rng.chance(2) {
                track = track.inner_radius(rng.value());
            }
            if rng.chance(2) {
                track = track.radial_size(rng.value());
            }
            if rng.chance(2) {
                track = track.unrooted_start(rng.value());
            }
            if rng.chance(2) {
                track = track.line_width(rng.value());
            }
            if rng.chance(2) {
                track = track.support_threshold(rng.value());
            }
            if rng.chance(2) {
                track = track.branch_label_size(rng.value());
            }
            if rng.chance(3) {
                track = track.scale_bar_length(rng.value());
            }
            if rng.chance(3) {
                track = track.reroot(rng.count());
            }
            if rng.chance(3) {
                track = track.collapse(rng.count());
            }
            if rng.chance(3) {
                track = track.reroot_midpoint();
            }
            if rng.chance(2) {
                track = track.support_style(match rng.below(3) {
                    0 => SupportStyle::Symbols,
                    1 => SupportStyle::Labels,
                    _ => SupportStyle::None,
                });
            }
            if rng.chance(2) {
                track = track.time(rng.name()).show_time_axis(true);
            }
            if rng.chance(2) {
                track = track.color_by(rng.name());
            }
            if rng.chance(2) {
                track = track.scale_bar();
            }
            // Last, because `dnds` clears `color_by` and `color_by` clears
            // `dnds`: putting the layers here exercises the rate colouring
            // whenever they reach for it and the categorical colouring the rest
            // of the time.
            let (track, _, _) = branch_layers(rng, track);
            Box::new(track)
        }
        29 => {
            let mut track = selection(rng, region).track;
            // A panel that is switched off still has to leave a valid document
            // behind, and the other half has to fill the band it is given.
            if rng.chance(4) {
                track = track.show_evidence(false);
            }
            if rng.chance(4) {
                track = track.show_effect(false);
            }
            Box::new(track)
        }
        30 => Box::new(phylodynamics(rng, region)),
        31 => {
            let observations = surveillance_observations(rng, region);
            Box::new(surveillance(rng, observations).track)
        }
        _ => Box::new(AxisTrack::new()),
    }
}

/// A tree for the tracks that take one, falling back when the text will not
/// parse, since a track wants a tree rather than an opinion about Newick.
fn tree(rng: &mut Lcg) -> Tree {
    let tame = rng.chance(2);
    let text = newick(rng, tame);
    Tree::parse_newick(&text)
        .unwrap_or_else(|_| Tree::parse_newick("((a:1,b:1):1,(c:1,d:1):1);").expect("a fixed tree"))
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
        // The writer's `num` guard covers attributes, and a tooltip is not an
        // attribute: it is a caller's number put through `format!` into text a
        // person reads. That is the gap `identity NaN` came through, having
        // been clamped by a `clamp` that propagates one.
        //
        // Not the infinities. A file can state a number too large for an f64
        // and it parses to one, so `text_exact` reports it verbatim on
        // purpose, and a reader seeing `inf` has been told what the file said.
        // A NaN is not a number a file states, it is the absence of one, and
        // printing it where a measurement goes is this crate's own bug class.
        for word in text.split(|c: char| !c.is_ascii_alphanumeric()) {
            if word.eq_ignore_ascii_case("nan") {
                bad.push(format!("a title presents a NaN as a measurement: {text:?}"));
            }
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

// ---------------------------------------------------------------------------
// Reading a drawing back
// ---------------------------------------------------------------------------

/// Every titled group in the document, as its title and what is drawn inside.
///
/// A tooltip is the only place a track writes a number a reader can hold the
/// picture to, and the group is what ties the two together. Depth is counted
/// rather than assumed, because a group can hold another one. The name of the
/// document never turns up here: it is a child of the root, not of a group.
fn titled_groups(svg: &str) -> Vec<(&str, &str)> {
    let mut out = Vec::new();
    let mut open: Vec<(Option<&str>, usize)> = Vec::new();
    let mut at = 0usize;
    while let Some(next) = svg[at..].find('<') {
        let here = at + next;
        let rest = &svg[here..];
        if rest.starts_with("</g>") {
            if let Some((Some(title), from)) = open.pop() {
                out.push((title, &svg[from..here]));
            }
            at = here + 4;
            continue;
        }
        if !rest.starts_with("<g>") && !rest.starts_with("<g ") {
            at = here + 1;
            continue;
        }
        let Some(end) = rest.find('>') else {
            break;
        };
        let inside = here + end + 1;
        open.push(match svg[inside..].strip_prefix("<title>") {
            Some(after) => match after.find("</title>") {
                Some(end) => (
                    Some(&after[..end]),
                    inside + "<title>".len() + end + "</title>".len(),
                ),
                None => (None, inside),
            },
            None => (None, inside),
        });
        at = inside;
    }
    out
}

/// One attribute of one element, by its exact name.
///
/// Through [`attributes`], so the name is matched whole: looking for `y` in an
/// element that carries a `stroke-dasharray` finds neither the dash pattern nor
/// the last letter of its name.
fn attribute<'a>(element: &'a str, name: &str) -> Option<&'a str> {
    attributes(element)
        .into_iter()
        .find(|(key, _)| *key == name)
        .map(|(_, value)| value)
}

/// One attribute of one element, as the number it has to be.
fn number(element: &str, name: &str) -> Option<f64> {
    attribute(element, name)?.parse().ok()
}

/// Every element of one kind in a fragment, as the text of its attributes.
fn elements<'a>(fragment: &'a str, tag: &str) -> Vec<&'a str> {
    fragment
        .split(tag)
        .skip(1)
        .filter_map(|piece| piece.split('>').next())
        .collect()
}

/// Every straight line in a fragment, in the order it was drawn.
fn lines(fragment: &str) -> Vec<((f64, f64), (f64, f64))> {
    elements(fragment, "<line ")
        .into_iter()
        .filter_map(|element| {
            Some((
                (number(element, "x1")?, number(element, "y1")?),
                (number(element, "x2")?, number(element, "y2")?),
            ))
        })
        .collect()
}

/// Every circle in a fragment, as its centre.
fn circles(fragment: &str) -> Vec<(f64, f64)> {
    elements(fragment, "<circle ")
        .into_iter()
        .filter_map(|element| Some((number(element, "cx")?, number(element, "cy")?)))
        .collect()
}

/// Every coordinate pair in a path or a point list, in the order it is visited.
///
/// The letters of a path drop out by failing to parse, which leaves the numbers
/// in the order they were written.
fn path_points(list: &str) -> Vec<(f64, f64)> {
    let numbers: Vec<f64> = list
        .split([' ', ','])
        .filter_map(|piece| piece.parse::<f64>().ok())
        .collect();
    numbers
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect()
}

/// The five entities the writer emits, put back.
///
/// The exact inverse of the escaping for everything these generators produce.
/// It cannot restore a character below `U+0020`, which is dropped rather than
/// escaped, and no name here carries one.
fn unescaped(text: &str) -> String {
    text.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
}

/// One field of a tooltip, by the words it opens with.
fn field<'a>(title: &'a str, prefix: &str) -> Option<&'a str> {
    title
        .split(" | ")
        .find_map(|part| part.strip_prefix(prefix))
}

/// Where a point sits along a segment, as a fraction of its length.
fn along(axis: ((f64, f64), (f64, f64)), point: (f64, f64)) -> Option<f64> {
    let (dx, dy) = (axis.1 .0 - axis.0 .0, axis.1 .1 - axis.0 .1);
    let length = dx.hypot(dy);
    (length > 0.0)
        .then(|| ((point.0 - axis.0 .0) * dx + (point.1 - axis.0 .1) * dy) / (length * length))
}

/// How far a number written for a person may sit from the number it came from.
///
/// Past the point where an f64 stops holding consecutive integers the writer
/// switches to an exponent and loses nothing.
fn rounding(printed: f64, decimals: i32) -> f64 {
    if printed.abs() >= (1u64 << 53) as f64 {
        0.0
    } else {
        0.5 * 10f64.powi(-decimals)
    }
}

// ---------------------------------------------------------------------------
// Selection scans
// ---------------------------------------------------------------------------

/// A rate estimate, hostile at the value that decides which of omega's three
/// answers a site gets.
///
/// Zero is not an edge case here, it is a third of the point: a dS of nought
/// with a dN above it is the site a scan cares most about, and it is the one
/// with no ratio a division can produce.
fn rate(rng: &mut Lcg) -> f64 {
    match rng.below(8) {
        0 => 0.0,
        1 => -0.0,
        2 => rng.value(),
        _ => (rng.below(1_000) as f64) / 100.0,
    }
}

/// A p-value or a posterior, hostile at both ends of the unit interval.
///
/// The thresholds themselves are in here on purpose. What separates a rule that
/// reads `<=` from one that reads `<` is one site at exactly 0.05, and a
/// generator that only rolls uniform numbers never produces one. So is a
/// genome-wide threshold, which is where a label written to three decimals
/// stops being a rounding of a number and becomes a different number.
fn probability(rng: &mut Lcg) -> f64 {
    match rng.below(16) {
        0 => 0.0,
        1 => 1.0,
        2 => 0.05,
        3 => 0.9,
        4 => -0.0,
        5 => 2.0,
        6 => -1.0,
        7 => f64::NAN,
        8 => f64::MIN_POSITIVE,
        9 => 5e-8,
        10 => 1e-300,
        11 => f64::INFINITY,
        _ => (rng.below(1_000) as f64) / 1_000.0,
    }
}

/// A selection track and the three things its drawing keeps to itself.
struct Selection {
    track: SelectionTrack,
    evidence: SelectionEvidence,
    /// The threshold the track settled on, worked out here rather than read
    /// back out of it. Both builders filter what they are handed, so this
    /// repeats the filter: a threshold taken from the track would agree with a
    /// mistake in it, and the disagreement is the case worth catching.
    threshold: f64,
    /// The two nonsynonymous rate classes of each site, in site order, as the
    /// builder accepted them. Nothing public reports them.
    classes: Vec<Option<(f64, f64)>>,
}

/// A selection scan over `region`, with rates, evidence and rate classes that
/// go missing as often as a real one's do.
///
/// No two sites share a position, for the reason the tip names do not: the
/// properties below key a marker by the position its tooltip opens with, and
/// two sites at one position are indistinguishable in exactly the comparison
/// they exist to make. Positions are offsets from the start of the window,
/// since a coordinate drawn absolutely falls outside every region that does not
/// begin at nought, and a scan with nothing in the window draws nothing.
fn selection(rng: &mut Lcg, region: &Region) -> Selection {
    let span = region.len();
    let count = rng.count();
    let mut sites = Vec::new();
    let mut classes = Vec::new();
    let mut used = std::collections::BTreeSet::new();
    for _ in 0..count {
        let pos = region.start().saturating_add(rng.position(span));
        // The tooltip counts from one and saturates, so the top two positions
        // print as one number.
        if !used.insert(pos.saturating_add(1)) {
            continue;
        }
        let mut site = SelectionSite::new(pos);
        if !rng.chance(8) {
            site = site.rates(rate(rng), rate(rng));
        }
        match rng.below(4) {
            0 => {}
            1 => site = site.p_value(probability(rng)),
            2 => site = site.posterior(probability(rng)),
            _ => site = site.p_value(probability(rng)).posterior(probability(rng)),
        }
        let mut accepted = None;
        if rng.chance(3) {
            let (minus, plus, weight) = (rate(rng), rate(rng), probability(rng));
            site = site.episodic_rates(minus, plus, weight);
            if minus.is_finite()
                && minus >= 0.0
                && plus.is_finite()
                && plus >= 0.0
                && weight.is_finite()
            {
                accepted = Some((minus, plus));
            }
        }
        if rng.chance(2) {
            site = site.label(rng.name());
        }
        sites.push(site);
        classes.push(accepted);
    }

    let evidence = if rng.chance(2) {
        SelectionEvidence::PValue
    } else {
        SelectionEvidence::Posterior
    };
    let mut track = SelectionTrack::new(sites).evidence(evidence);
    let mut threshold = match evidence {
        SelectionEvidence::PValue => 0.05,
        SelectionEvidence::Posterior => 0.90,
    };
    if rng.chance(2) {
        let asked = probability(rng);
        track = track.p_threshold(asked);
        if matches!(evidence, SelectionEvidence::PValue)
            && asked.is_finite()
            && asked > 0.0
            && asked <= 1.0
        {
            threshold = asked;
        }
    }
    if rng.chance(2) {
        let asked = probability(rng);
        track = track.posterior_threshold(asked);
        if matches!(evidence, SelectionEvidence::Posterior) && asked.is_finite() {
            threshold = asked.clamp(0.0, 1.0);
        }
    }
    if rng.chance(2) {
        track = track.neutral_band(rng.value(), rng.value());
    }
    if rng.chance(2) {
        track = track.saturation(rng.value());
    }
    if rng.chance(3) {
        track = track.height(match rng.below(6) {
            0 => f64::NAN,
            1 => -1.0,
            2 => 0.0,
            3 => 54.0,
            4 => 55.0,
            _ => 60.0 + rng.below(400) as f64,
        });
    }
    if rng.chance(2) {
        track = track.label(rng.name());
    }
    Selection {
        track,
        evidence,
        threshold,
        classes,
    }
}

/// One marker: the site it belongs to, how high it reaches, its shape and the
/// elements drawn inside it.
struct Marker {
    site: usize,
    title: String,
    body: String,
    y: f64,
    diamond: bool,
}

/// One selection track drawn on its own, taken apart into its two panels.
struct Scan {
    region: Region,
    selection: Selection,
    svg: String,
    /// The dashed significance rule.
    rule: f64,
    evidence: Vec<Marker>,
    effect: Vec<Marker>,
}

/// Draws one generated track on its own and files every marker under its site.
///
/// Both panels are turned on whatever the generator chose, since a panel that
/// is switched off has no markers to check. Nothing here is told where the
/// panels are: a marker is a stem and a symbol, and the baseline the stem rises
/// from is what says which panel drew it, matched against the rules the track
/// lays down before any marker. That way a property compares the marks against
/// the lines the same drawing put beside them.
fn scan(seed: u64) -> Scan {
    let mut rng = Lcg::new(seed);
    let region = region(&mut rng);
    let selection = selection(&mut rng, &region);
    let svg = Figure::new(region.clone())
        .show_region_label(false)
        .push(
            selection
                .track
                .clone()
                .show_evidence(true)
                .show_effect(true),
        )
        .to_svg();

    let (plain, dashed) = rules(&svg);
    assert_eq!(
        (plain.len(), dashed.len()),
        (2, 1),
        "seed {seed}: a track with both panels on drew {} plain rules and {} dashed",
        plain.len(),
        dashed.len()
    );
    let (baseline, mid) = (plain[0], plain[1]);

    let mut at: BTreeMap<String, usize> = BTreeMap::new();
    for (index, site) in selection.track.sites().iter().enumerate() {
        at.insert(format!("position {}", site.pos.saturating_add(1)), index);
    }

    let (mut evidence, mut effect) = (Vec::new(), Vec::new());
    for (title, body) in titled_groups(&svg) {
        let opening = title.split(" | ").next().unwrap_or_default();
        if !opening.starts_with("position ") {
            continue;
        }
        let site = *at
            .get(opening)
            .unwrap_or_else(|| panic!("seed {seed}: a marker for no site: {title:?}"));
        let stem = *lines(body)
            .first()
            .unwrap_or_else(|| panic!("seed {seed}: a marker with no stem: {title:?}"));
        let marker = Marker {
            site,
            title: title.to_string(),
            body: body.to_string(),
            y: stem.1 .1,
            diamond: body.contains("<polygon"),
        };
        if stem.0 .1 == baseline {
            evidence.push(marker);
        } else if stem.0 .1 == mid {
            effect.push(marker);
        } else {
            panic!(
                "seed {seed}: a marker rises from {}, not {baseline} or {mid}",
                stem.0 .1
            );
        }
    }
    Scan {
        region,
        selection,
        svg,
        rule: dashed[0],
        evidence,
        effect,
    }
}

impl Scan {
    /// Whether the track's threshold picks this site out, worked out from the
    /// numbers the site kept rather than from the track's own answer.
    fn crossed(&self, marker: &Marker) -> bool {
        let site = &self.selection.track.sites()[marker.site];
        match self.selection.evidence {
            SelectionEvidence::PValue => site.p().is_some_and(|p| p <= self.selection.threshold),
            SelectionEvidence::Posterior => site
                .probability()
                .is_some_and(|p| p >= self.selection.threshold),
        }
    }
}

/// The horizontal rules a track lays down before any marker, plain and dashed:
/// the skyline baseline, the line at omega one, and the significance rule
/// between them.
fn rules(svg: &str) -> (Vec<f64>, Vec<f64>) {
    let (mut plain, mut dashed) = (Vec::new(), Vec::new());
    let mut outside = String::new();
    let mut rest = svg;
    while let Some(at) = rest.find("<g><title>") {
        outside.push_str(&rest[..at]);
        match rest[at..].find("</g>") {
            Some(end) => rest = &rest[at + end + 4..],
            None => {
                rest = "";
                break;
            }
        }
    }
    outside.push_str(rest);
    for element in elements(&outside, "<line ") {
        let (Some(y1), Some(y2)) = (number(element, "y1"), number(element, "y2")) else {
            continue;
        };
        if y1 != y2 {
            continue;
        }
        if element.contains("stroke-dasharray") {
            dashed.push(y1);
        } else {
            plain.push(y1);
        }
    }
    (plain, dashed)
}

#[test]
fn every_estimated_site_is_drawn_in_every_panel_that_can_show_it() {
    // The module's promise is that a missing estimate is left out rather than
    // drawn as zero, which is worth nothing unless the ones that are there all
    // arrive. The two panels decide separately, and one of the decisions is the
    // point of the track: a dS of nought with a dN above it has no ratio a
    // division can produce and is still the strongest site on the page. A
    // drawing that quietly drops it loses the finding, and nobody counts the
    // marks on a plot against the rows of their table.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        let mut expected: BTreeMap<usize, usize> = BTreeMap::new();
        for (index, site) in scan.selection.track.sites().iter().enumerate() {
            if !scan.region.contains(site.pos) {
                continue;
            }
            // A p-value that reached the site is between nought and one, so its
            // -log10 is a finite number the skyline has room for.
            let evidence = match scan.selection.evidence {
                SelectionEvidence::PValue => site.p().is_some(),
                SelectionEvidence::Posterior => site.probability().is_some(),
            };
            let count = usize::from(evidence) + usize::from(site.omega().is_some());
            if count > 0 {
                expected.insert(index, count);
            }
        }
        let mut drawn: BTreeMap<usize, usize> = BTreeMap::new();
        for marker in scan.evidence.iter().chain(scan.effect.iter()) {
            *drawn.entry(marker.site).or_default() += 1;
        }
        assert_eq!(
            expected, drawn,
            "seed {seed}: the wrong sites reached the page"
        );
    }
}

#[test]
fn the_effect_strip_never_draws_a_larger_omega_below_a_smaller_one() {
    // omega has three answers and the strip has to honour all three. Its axis
    // says omega at the top and one over omega at the bottom, so height is the
    // whole encoding, and the answer that breaks it is the division by zero:
    // drawn at the neutral line an infinite omega says the site is under no
    // selection at all, and drawn at the floor it says the opposite of what it
    // is. Where any particular omega belongs is a matter of the cap and the
    // saturation point. The order is not.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        let sites = scan.selection.track.sites();
        for above in &scan.effect {
            for below in &scan.effect {
                let (high, low) = (
                    sites[above.site].omega().expect("a marker means an omega"),
                    sites[below.site].omega().expect("a marker means an omega"),
                );
                assert!(
                    high <= low || above.y <= below.y + 0.01,
                    "seed {seed}: omega {high} is drawn at {} and omega {low} at {}",
                    above.y,
                    below.y
                );
            }
        }
    }
}

#[test]
fn the_skyline_says_twice_whether_a_site_crossed_the_threshold() {
    // Two channels carry one claim, and they are computed apart. The shape is
    // the only cue that survives a greyscale print, since the colour is spent
    // on the direction of selection rather than on the evidence, and it is the
    // one the header explains: it prints `p<=0.05`, and a site at exactly 0.05
    // is what that promise is about, because the difference between a rule
    // written `<=` and one written `<` is invisible in every example anyone
    // writes and decides which sites a reader reports.
    //
    // The dashed rule is the second channel, and it is the reading instruction
    // for the whole skyline: a marker on the wrong side of it is a finding
    // thrown away or invented. That half needs no threshold at all, only the
    // line the same drawing put beside the mark. A marker sitting on the line
    // makes no claim either way, which is where a site at exactly the threshold
    // belongs; the shape carries that case instead.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        for marker in scan.evidence.iter().chain(scan.effect.iter()) {
            assert_eq!(
                marker.diamond,
                scan.crossed(marker),
                "seed {seed}: {:?} is drawn as a {}",
                marker.title,
                if marker.diamond { "diamond" } else { "circle" }
            );
        }
        for marker in &scan.evidence {
            if (marker.y - scan.rule).abs() <= 0.01 {
                continue;
            }
            assert_eq!(
                marker.y < scan.rule,
                marker.diamond,
                "seed {seed}: {:?} is drawn at {} and its rule is at {}",
                marker.title,
                marker.y,
                scan.rule
            );
        }
    }
}

#[test]
fn a_tooltip_divides_the_two_numbers_it_prints() {
    // The tooltip is the only place the rates survive, and it prints all three,
    // so a reader can check the division. One that reads `dS 0 | dN 1 | omega
    // 10000000` asks them to believe that nought divides into one ten million
    // times. Rounding is allowed to lose a digit of a rate; it is not allowed to
    // turn a rate into a nought that contradicts the ratio printed beside it,
    // and it is not allowed to swap the two, which nothing else here would
    // notice.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        for (title, _) in titled_groups(&scan.svg) {
            let (Some(ds), Some(dn), Some(omega)) = (
                field(title, "dS "),
                field(title, "dN "),
                field(title, "omega "),
            ) else {
                continue;
            };
            let (Ok(ds), Ok(dn)) = (ds.parse::<f64>(), dn.parse::<f64>()) else {
                panic!("seed {seed}: unreadable rates in {title:?}")
            };
            if ds == 0.0 {
                let expected = if dn > 0.0 {
                    "infinity"
                } else {
                    "not estimable"
                };
                assert_eq!(omega, expected, "seed {seed}: {title:?}");
                continue;
            }
            let printed: f64 = omega
                .parse()
                .unwrap_or_else(|_| panic!("seed {seed}: {title:?} has rates and no ratio"));
            let low = (dn - rounding(dn, 5)).max(0.0) / (ds + rounding(ds, 5));
            let floor = (ds - rounding(ds, 5)).max(0.0);
            let high = if floor > 0.0 {
                (dn + rounding(dn, 5)) / floor
            } else {
                f64::INFINITY
            };
            let slack = rounding(printed, 4);
            assert!(
                printed + slack >= low && printed - slack <= high,
                "seed {seed}: {title:?} says {printed}, its own numbers say {low} to {high}"
            );
        }
    }
}

#[test]
fn two_markers_a_tooltip_cannot_tell_apart_are_drawn_at_one_height() {
    // The skyline is drawn on a log scale and the tooltip is rounded on a
    // linear one, which is a bug waiting for a small enough p-value: a reader
    // who hovers two marks at different heights to find out why is told they
    // are the same number. No oracle here, and none needed, since this is the
    // document contradicting itself. The tolerance is a twentieth of a pixel,
    // which is smaller than any difference the drawing can mean and larger than
    // any two numbers that print alike can produce.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        let prefix = match scan.selection.evidence {
            SelectionEvidence::PValue => "p ",
            SelectionEvidence::Posterior => "posterior ",
        };
        let mut seen: BTreeMap<&str, &Marker> = BTreeMap::new();
        for marker in &scan.evidence {
            let Some(printed) = field(&marker.title, prefix) else {
                continue;
            };
            if let Some(other) = seen.insert(printed, marker) {
                assert!(
                    (other.y - marker.y).abs() <= 0.05,
                    "seed {seed}: {prefix}{printed} is drawn at {} and at {}",
                    other.y,
                    marker.y
                );
            }
        }
    }
}

#[test]
fn two_equal_rate_classes_are_painted_one_colour() {
    // The capsule divides both nonsynonymous classes by the synonymous rate,
    // and a synonymous rate of nought is the case `omega` answers three ways on
    // purpose. Whatever the right colour for it is, two equal rates have to get
    // the same one: painting them the two ends of the scale tells the reader
    // that one site is under the strongest purifying and the strongest positive
    // selection at once. Colours need no oracle for that.
    for seed in 0..ROUNDS {
        let scan = scan(seed);
        for marker in &scan.evidence {
            let Some((minus, plus)) = scan.selection.classes[marker.site] else {
                continue;
            };
            if minus != plus {
                continue;
            }
            // The surround, then the two halves. A class with no width is not
            // drawn at all, and then there is nothing to compare.
            let fills: Vec<&str> = elements(&marker.body, "<rect ")
                .into_iter()
                .filter_map(|element| attribute(element, "fill"))
                .collect();
            if fills.len() < 3 {
                continue;
            }
            assert_eq!(
                fills[1], fills[2],
                "seed {seed}: two rate classes of {minus} are painted {} and {}: {:?}",
                fills[1], fills[2], marker.title
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Trajectories
// ---------------------------------------------------------------------------

/// A trajectory, hostile in the two ways a real one is.
///
/// The estimate and the bounds are drawn from the same pool and never from each
/// other, because nothing ties them together on the way in either: `interval`
/// checks that the bounds are finite and the right way round and says nothing
/// about the estimate, so a point sitting outside its own confidence interval is
/// not an edge case here. It is what two columns of a table that were never
/// checked against each other look like. The log axis gets the same values as
/// the linear one, since an estimate at or below nought is what a growth rate
/// crossing zero is.
fn phylodynamics(rng: &mut Lcg, region: &Region) -> PhylodynamicTrack {
    let span = region.len();
    let points: Vec<PhylodynamicPoint> = (0..rng.count())
        .map(|_| {
            let time = region.start().saturating_add(rng.position(span));
            let mut point = PhylodynamicPoint::new(time, rng.value());
            if rng.chance(2) {
                let (a, b) = (rng.value(), rng.value());
                point = point.interval(a.min(b), a.max(b));
            }
            if rng.chance(3) {
                point = point.label(rng.name());
            }
            point
        })
        .collect();
    let mut track = PhylodynamicTrack::new(points);
    if rng.chance(2) {
        track = track.scale(PhylodynamicScale::Log10);
    }
    if rng.chance(2) {
        track = track.label(rng.name());
    }
    if rng.chance(3) {
        track = track.unit(rng.name());
    }
    if rng.chance(3) {
        track = track.height(rng.value());
    }
    if rng.chance(3) {
        track = track.reference(rng.value(), rng.name());
    }
    if rng.chance(3) {
        track = track.show_points(false);
    }
    if rng.chance(4) {
        track = track.show_interval(false);
    }
    track
}

/// A trajectory kept as what it was drawn from.
///
/// Separate from the [`track`] arm because the properties below read numbers
/// back out of the document and hold the marks to them. A unit in the way would
/// make an axis label the caller's string rather than the crate's arithmetic, a
/// point label would make a tooltip one, and a second track would put another
/// track's shapes in the same file. Everything that changes what gets drawn is
/// still generated: the scale, the reference, the two switches.
struct Trajectory {
    region: Region,
    points: Vec<PhylodynamicPoint>,
    scale: PhylodynamicScale,
    reference: Option<f64>,
    markers: bool,
    ribbon: bool,
}

impl Trajectory {
    fn track(&self) -> PhylodynamicTrack {
        let mut track = PhylodynamicTrack::new(self.points.clone())
            .scale(self.scale)
            .show_points(self.markers)
            .show_interval(self.ribbon);
        if let Some(reference) = self.reference {
            track = track.reference(reference, "ref");
        }
        track
    }

    fn to_svg(&self) -> String {
        Figure::new(self.region.clone())
            .show_region_label(false)
            .push(self.track())
            .to_svg()
    }
}

fn trajectory(seed: u64) -> Trajectory {
    let mut rng = Lcg::new(seed);
    let region = region(&mut rng);
    let span = region.len();
    let scale = if rng.chance(2) {
        PhylodynamicScale::Log10
    } else {
        PhylodynamicScale::Linear
    };
    let points = (0..rng.count())
        .map(|_| {
            let time = region.start().saturating_add(rng.position(span));
            let mut point = PhylodynamicPoint::new(time, rng.value());
            if rng.chance(2) {
                let (a, b) = (rng.value(), rng.value());
                point = point.interval(a.min(b), a.max(b));
            }
            point
        })
        .collect();
    let reference = if rng.chance(3) {
        Some(rng.value())
    } else {
        None
    };
    Trajectory {
        region,
        points,
        scale,
        reference,
        markers: !rng.chance(3),
        ribbon: !rng.chance(4),
    }
}

/// The time and the estimate a point's tooltip claims.
///
/// Read from the front, so that a unit and a label, neither of which is parsed
/// here, cannot move the two numbers.
fn point_estimate(title: &str) -> Option<(u64, f64)> {
    let rest = title.strip_prefix("time ")?;
    let (time, rest) = rest.split_once(" | estimate ")?;
    Some((time.parse().ok()?, rest.split(' ').next()?.parse().ok()?))
}

/// Every `<text>` in the document, as its anchor and what it says.
fn labels(svg: &str) -> Vec<(&str, &str)> {
    let mut out = Vec::new();
    let mut rest = svg;
    while let Some(at) = rest.find("<text ") {
        rest = &rest[at..];
        let Some(head) = rest.find('>') else {
            break;
        };
        let anchor = attribute(&rest[..head], "text-anchor").unwrap_or("");
        rest = &rest[head + 1..];
        let Some(end) = rest.find("</text>") else {
            break;
        };
        out.push((anchor, &rest[..end]));
        rest = &rest[end..];
    }
    out
}

#[test]
fn every_point_a_trajectory_names_is_drawn() {
    // The track filters its points, then opens a titled group for each survivor
    // and draws the mark inside it. A titled group with nothing in it is a
    // tooltip over blank paper: the reader is told a value is there, cannot
    // point at it, and has no way to tell it apart from a value that was never
    // in the file. The writer drops geometry it cannot write rather than writing
    // a number that is not one, which is the right answer for the document and a
    // silent one for the figure, so this is where an arithmetic failure upstream
    // of the writer goes quiet.
    for seed in 0..ROUNDS {
        let svg = trajectory(seed).to_svg();
        for (title, body) in titled_groups(&svg) {
            if point_estimate(title).is_none() {
                continue;
            }
            assert!(
                body.contains("<circle"),
                "seed {seed}: {title:?} names a point with nothing drawn for it"
            );
        }
    }
}

#[test]
fn every_number_up_a_value_axis_is_a_number() {
    // The rules up the side are the only thing in the panel that says what
    // height means; without them the line is a shape with no units. Worse than a
    // missing label is a wrong one: an axis tick reading `inf` is a figure that
    // will be read and believed. Nothing the caller supplied is in the way here,
    // because the figure carries one track, no unit and no region label, so
    // every string counted is the crate's own arithmetic written out.
    for seed in 0..ROUNDS {
        let svg = trajectory(seed).to_svg();
        if !titled_groups(&svg)
            .iter()
            .any(|(title, _)| point_estimate(title).is_some())
        {
            continue;
        }
        let ticks: Vec<&str> = labels(&svg)
            .into_iter()
            .filter(|(anchor, text)| *anchor == "end" && *text != "ref")
            .map(|(_, text)| text)
            .collect();
        assert!(
            !ticks.is_empty(),
            "seed {seed}: a trajectory was drawn against no axis at all"
        );
        for tick in &ticks {
            let number: f64 = tick
                .parse()
                .unwrap_or_else(|_| panic!("seed {seed}: the axis is labelled {tick:?}"));
            assert!(
                number.is_finite(),
                "seed {seed}: the axis is labelled {tick:?}"
            );
        }
    }
}

#[test]
fn a_bigger_estimate_is_never_drawn_lower_down() {
    // One vertical map for the whole track, and it points up. The ribbon, the
    // reference guide and the axis are all measured against it, so a mark placed
    // by anything other than its own value stops the panel being readable as a
    // trajectory. The case this is really for is the honest one: `interval` does
    // not check that the estimate lies between the bounds, and a file where it
    // does not is a file, so the standing temptation is to pull the mark into
    // its own ribbon on the way past. That draws an agreement the data does not
    // have, and it arrives here as a bigger estimate sitting lower down. The
    // comparison is on the tooltip's rounded value, which is safe: rounding
    // cannot swap two values over, only make them equal, and equal ones are
    // skipped.
    for seed in 0..ROUNDS {
        let svg = trajectory(seed).to_svg();
        let mut marks: Vec<(f64, f64)> = Vec::new();
        for (title, body) in titled_groups(&svg) {
            let Some((_, estimate)) = point_estimate(title) else {
                continue;
            };
            let Some(centre) = circles(body).first().copied() else {
                continue;
            };
            marks.push((estimate, centre.1));
        }
        for above in &marks {
            for below in &marks {
                if above.0 <= below.0 {
                    continue;
                }
                assert!(
                    above.1 <= below.1 + 0.002,
                    "seed {seed}: {} is drawn at {} and the smaller {} at {}",
                    above.0,
                    above.1,
                    below.0,
                    below.1
                );
            }
        }
    }
}

#[test]
fn a_point_outside_the_window_changes_nothing_inside_it() {
    // The vertical axis is a claim about what is on screen. A phylodynamic table
    // covers years and a figure shows a month of it, so most of the table is
    // always off the window, and an axis fitted to all of it would draw the
    // visible part as a flat line against a scale nothing on the page explains.
    // The same track would then come out different in two panels that differ
    // only in where they are looking, which is the comparison such a figure is
    // made for.
    for seed in 0..ROUNDS {
        let mut trajectory = trajectory(seed);
        let before = trajectory.to_svg();
        // One past the last base of the window, which saturates to a coordinate
        // still outside it when the window ends at the top of the range.
        let outside = trajectory.region.end();
        trajectory
            .points
            .push(PhylodynamicPoint::new(outside, 1e6).interval(-1e6, 1e6));
        assert_eq!(
            before,
            trajectory.to_svg(),
            "seed {seed}: a point at {outside}, outside {}, changed the drawing",
            trajectory.region
        );
    }
}

#[test]
fn the_same_trajectory_drawn_twice_draws_the_same_shape_twice() {
    // A track is a function of the band it is handed, so the same data in two
    // bands is the same shape at two heights: the same x, one constant drop. A
    // coordinate that does not move between the copies was never computed from
    // the band. That matters because a number that is not a number is written as
    // a nought rather than dropped when it reaches a path, which has no other
    // check: a ribbon vertex whose arithmetic came out NaN is not left out the
    // way a circle is. It is drawn, at the top of the sheet, over whatever track
    // is up there.
    for seed in 0..ROUNDS {
        let trajectory = trajectory(seed);
        let svg = Figure::new(trajectory.region.clone())
            .show_region_label(false)
            .push(trajectory.track())
            .push(trajectory.track())
            .to_svg();
        let mut shapes: Vec<(&str, Vec<(f64, f64)>)> = Vec::new();
        for (title, body) in titled_groups(&svg) {
            let mut points = Vec::new();
            for element in elements(body, "<path ") {
                if let Some(d) = attribute(element, "d") {
                    points.extend(path_points(d));
                }
            }
            points.extend(circles(body));
            if !points.is_empty() {
                shapes.push((title, points));
            }
        }
        assert!(
            shapes.len() % 2 == 0,
            "seed {seed}: {} shapes, which is not two of each",
            shapes.len()
        );
        let half = shapes.len() / 2;
        let mut drop: Option<f64> = None;
        for index in 0..half {
            let (first_title, first) = &shapes[index];
            let (second_title, second) = &shapes[index + half];
            assert_eq!(
                first_title, second_title,
                "seed {seed}: the two copies drew different things"
            );
            assert_eq!(
                first.len(),
                second.len(),
                "seed {seed}: {first_title:?} has {} vertices in one copy and {} in the other",
                first.len(),
                second.len()
            );
            for (one, two) in first.iter().zip(second.iter()) {
                assert!(
                    (one.0 - two.0).abs() < 0.002,
                    "seed {seed}: {first_title:?} moved sideways between the copies, {} then {}",
                    one.0,
                    two.0
                );
                let moved = two.1 - one.1;
                let known = *drop.get_or_insert(moved);
                assert!(
                    (moved - known).abs() < 0.01,
                    "seed {seed}: {first_title:?} dropped {moved} between the copies where \
                     everything else dropped {known}"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Surveillance panels
// ---------------------------------------------------------------------------

/// Observations over few enough lineages and times that they collide.
///
/// Drawn freely across the coordinate range, two of them would never share a
/// time, and the panel's opinions about a duplicated lineage and an incomplete
/// column would never be reached. The count and the total are drawn separately
/// on purpose: a count above its total is not a typo, it is a numerator counted
/// over a day against a denominator counted over a week, and it is the case the
/// division has no answer for.
fn surveillance_observations(rng: &mut Lcg, region: &Region) -> Vec<SurveillanceObservation> {
    let span = region.len();
    let pool = if rng.chance(3) { 40 } else { 4 };
    let lineages: Vec<String> = (0..1 + rng.below(pool)).map(|_| rng.name()).collect();
    // Offsets from the start of the window, since a time is a coordinate and a
    // window rarely starts at nought: drawn absolutely, nearly every observation
    // would fall outside the region and the panel would have almost nothing to
    // draw. They are still hostile, because nought is the first base, the span
    // is one past the last, and the top of the range saturates.
    let times: Vec<u64> = (0..1 + rng.below(5))
        .map(|_| region.start().saturating_add(rng.position(span)))
        .collect();
    let counted = |rng: &mut Lcg, time: u64, lineage: &str| -> SurveillanceObservation {
        let total = match rng.below(8) {
            0 => 0,
            1 => 1,
            2 => u64::MAX,
            _ => rng.below(1_000),
        };
        let count = match rng.below(8) {
            0 => 0,
            1 => total,
            2 => total.saturating_add(1),
            3 => u64::MAX,
            _ => rng.below(1_000),
        };
        let observation = SurveillanceObservation::new(time, lineage, count, total);
        if rng.chance(3) {
            observation.label(rng.name())
        } else {
            observation
        }
    };

    // A composition table is a grid: every lineage counted at every time. Row by
    // row at random it practically never comes out complete, and a complete
    // column is the case the stacked drawing exists for. The holes and the
    // repeats are here because a sample sheet is never complete either, and both
    // are things the panel has an opinion about.
    if rng.chance(3) {
        let mut out = Vec::new();
        for lineage in &lineages {
            for time in &times {
                if rng.chance(8) {
                    continue;
                }
                out.push(counted(rng, *time, lineage));
                if rng.chance(12) {
                    out.push(counted(rng, *time, lineage));
                }
            }
        }
        return out;
    }

    (0..rng.count())
        .map(|_| {
            let time = times[rng.below(times.len() as u64) as usize];
            let lineage = lineages[rng.below(lineages.len() as u64) as usize].clone();
            counted(rng, time, &lineage)
        })
        .collect()
}

/// One surveillance panel and the two settings it keeps to itself.
///
/// The metric and the sampling floor come back with the track because it does
/// not report them, and both are about what a setting is allowed to leave out of
/// a document.
struct Panel {
    track: SurveillanceTrack,
    observations: Vec<SurveillanceObservation>,
    metric: SurveillanceMetric,
    floor: u64,
}

/// One surveillance panel over `observations`, with every switch it has.
///
/// The observations are passed in rather than generated here so that a property
/// can build the same panel twice over the same numbers in a different order:
/// replaying the generator from the seed gives back the same switches, and only
/// the order changes.
fn surveillance(rng: &mut Lcg, observations: Vec<SurveillanceObservation>) -> Panel {
    let mut track = SurveillanceTrack::new(observations.clone());
    let metric = if rng.chance(2) {
        SurveillanceMetric::Count
    } else {
        SurveillanceMetric::Frequency
    };
    track = track.metric(metric);
    if rng.chance(2) {
        track = track.style(if rng.chance(2) {
            SurveillanceStyle::Lines
        } else {
            SurveillanceStyle::Stacked
        });
    }
    let floor = match rng.below(5) {
        0 => 0,
        1 => 1,
        2 => u64::MAX,
        _ => rng.below(100),
    };
    track = track.minimum_total(floor);
    if rng.chance(2) {
        track = track.frequency_alert(rng.value());
    }
    if rng.chance(2) {
        track = track.growth_alert(rng.value());
    }
    if rng.chance(3) {
        // Bounded, unlike the other builders here: the guard on this one takes
        // any finite number, and a band ten to the three hundred and eighth tall
        // makes every row below it arithmetic about nothing.
        track = track.height(match rng.below(6) {
            0 => f64::NAN,
            1 => f64::NEG_INFINITY,
            2 => -1.0,
            3 => 0.0,
            4 => 48.0,
            _ => 60.0 + rng.below(600) as f64,
        });
    }
    if rng.chance(3) {
        track = track.show_points(false);
    }
    if rng.chance(3) {
        track = track.label(rng.name());
    }
    Panel {
        track,
        observations,
        metric,
        floor,
    }
}

/// The observation a tooltip names, as the three numbers that identify it.
fn named_observation(tooltip: &str) -> Option<(u64, u64, u64)> {
    let time = field(tooltip, "time ")?.parse().ok()?;
    let (count, total) = field(tooltip, "count ")?.split_once(" of ")?;
    Some((time, count.parse().ok()?, total.parse().ok()?))
}

/// A figure holding one surveillance panel and nothing else.
fn one_panel(panel: &Panel, region: &Region) -> Figure {
    Figure::new(region.clone())
        .show_region_label(false)
        .push(panel.track.clone())
}

/// The band a one-track figure gave its track, read from the clip it wrote.
///
/// The clip is what a renderer will actually crop to, so it is the honest
/// boundary to measure against rather than a rectangle worked out here.
fn clipped_band(svg: &str) -> (f64, f64) {
    let inside = svg
        .split("<clipPath")
        .nth(1)
        .expect("one clip for one track");
    let rectangle = inside.split("/>").next().expect("a clipping rectangle");
    let read = |name: &str| number(rectangle, name).expect("a number");
    let top = read("y");
    (top, top + read("height"))
}

/// Every tooltip and the rows the shapes under it are drawn on.
fn mark_rows(svg: &str) -> Vec<(&str, Vec<f64>)> {
    let mut out = Vec::new();
    for (tooltip, body) in titled_groups(svg) {
        let mut rows = Vec::new();
        for element in elements(body, "<circle ") {
            rows.extend(number(element, "cy"));
        }
        for element in elements(body, "<line ") {
            rows.extend(number(element, "y1"));
            rows.extend(number(element, "y2"));
        }
        // A polygon and a filled area carry their coordinates in pairs, and the
        // second of each pair is the one that can leave the band.
        for (tag, name) in [("<polygon ", "points"), ("<path ", "d")] {
            for element in elements(body, tag) {
                if let Some(list) = attribute(element, name) {
                    rows.extend(path_points(list).into_iter().map(|(_, row)| row));
                }
            }
        }
        out.push((tooltip, rows));
    }
    out
}

/// Whether two numbers can be read as a frequency, by the definition rather than
/// by asking the track, which would agree with a mistake in it.
fn is_a_frequency(observation: &SurveillanceObservation) -> bool {
    observation.total > 0 && observation.count <= observation.total
}

/// Two rows for one lineage at one time, which the panel says it will not join.
///
/// Generated lists produce these on their own but not often enough to trust, so
/// one is added, along with the later row whose change is measured against them.
/// A table joined on a date produces nothing else.
fn twinned(
    rng: &mut Lcg,
    observations: &[SurveillanceObservation],
) -> Option<Vec<SurveillanceObservation>> {
    let readable: Vec<&SurveillanceObservation> =
        observations.iter().filter(|o| is_a_frequency(o)).collect();
    if readable.is_empty() {
        return None;
    }
    let pick = readable[rng.below(readable.len() as u64) as usize];
    // The complement, so the two rows disagree by as much as two rows can. Two
    // counts a few apart out of a denominator of ten to the nineteenth print the
    // same frequency, and then the order has nothing to change.
    let twin = SurveillanceObservation::new(
        pick.time,
        pick.lineage.clone(),
        pick.total - pick.count,
        pick.total,
    );
    if twin == *pick {
        return None;
    }
    let mut out = observations.to_vec();
    out.push(twin);
    if let Some(later) = pick.time.checked_add(1) {
        out.push(SurveillanceObservation::new(
            later,
            pick.lineage.clone(),
            pick.count,
            pick.total,
        ));
    }
    Some(out)
}

/// The first two observations that share a lineage and a time and differ.
fn twins(observations: &[SurveillanceObservation]) -> Option<(usize, usize)> {
    for first in 0..observations.len() {
        for second in (first + 1)..observations.len() {
            if observations[first].time == observations[second].time
                && observations[first].lineage == observations[second].lineage
                && observations[first] != observations[second]
            {
                return Some((first, second));
            }
        }
    }
    None
}

#[test]
fn every_observation_is_named_once_or_left_out_on_purpose() {
    // The failure this exists for is silent: a row given to the panel, not
    // drawn, and not mentioned. Nobody counts the marks on a plot against the
    // rows of their table. There are exactly two reasons a row may be left out
    // without a word, and both are the caller's own instructions: it sits
    // outside the drawn region, or its denominator is below the sampling floor
    // the caller set. Anything else is named, once, with its own numbers,
    // whether or not the panel could place it.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let observations = surveillance_observations(&mut rng, &window);
        let panel = surveillance(&mut rng, observations);
        let svg = one_panel(&panel, &window).to_svg();

        let mut named: Vec<(u64, u64, u64)> = titled_groups(&svg)
            .into_iter()
            .filter_map(|(tooltip, _)| named_observation(tooltip))
            .collect();
        let mut given: Vec<(u64, u64, u64)> = panel
            .observations
            .iter()
            .filter(|o| window.contains(o.time))
            .filter(|o| {
                o.total >= panel.floor
                    || (panel.metric == SurveillanceMetric::Frequency && !is_a_frequency(o))
            })
            .map(|o| (o.time, o.count, o.total))
            .collect();
        named.sort();
        given.sort();
        assert_eq!(
            named,
            given,
            "seed {seed}: the panel names {} of the {} observations it has to account for",
            named.len(),
            given.len()
        );
    }
}

#[test]
fn no_mark_is_drawn_outside_the_band_it_belongs_to() {
    // A track is clipped to its own band, so a mark placed outside one is not a
    // mark drawn in the wrong place: it is a mark drawn nowhere, with its
    // tooltip attached to nothing. That is how a frequency above one used to
    // leave the figure, and it is the same loss whichever end it leaves by. The
    // band comes from the clip the document itself carries, so the check is
    // against what the reader's renderer will do rather than against a number
    // this test worked out for itself.
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let observations = surveillance_observations(&mut rng, &window);
        let panel = surveillance(&mut rng, observations);
        let svg = one_panel(&panel, &window).to_svg();
        let (top, bottom) = clipped_band(&svg);
        for (tooltip, rows) in mark_rows(&svg) {
            for row in rows {
                assert!(
                    row >= top - 0.001 && row <= bottom + 0.001,
                    "seed {seed}: a mark at {row} is outside the band {top} to {bottom}: {tooltip}"
                );
            }
        }
    }
}

#[test]
fn reordering_two_rows_the_panel_will_not_join_changes_nothing() {
    // A metamorphic property, and the panel's own words are the specification:
    // it says duplicate lineage/time observations are not joined. Then the order
    // the two arrived in cannot matter, because the only thing that could read
    // it is a join. The lineages keep their order, so nothing else about the
    // figure is allowed to move either.
    let mut checked = 0u64;
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let observations = surveillance_observations(&mut rng, &window);
        let Some(observations) = twinned(&mut rng, &observations) else {
            continue;
        };
        let Some((first, second)) = twins(&observations) else {
            continue;
        };
        let mut swapped = observations.clone();
        swapped.swap(first, second);

        // The same panel twice, since the switches are chosen from the seed and a
        // track cannot be rebuilt from itself.
        let drawn = |list: Vec<SurveillanceObservation>| -> Vec<String> {
            let mut rng = Lcg::new(seed);
            let window = region(&mut rng);
            let _ = surveillance_observations(&mut rng, &window);
            let panel = surveillance(&mut rng, list);
            let svg = one_panel(&panel, &window).to_svg();
            let mut out: Vec<String> = titled_groups(&svg)
                .into_iter()
                .map(|(tooltip, _)| tooltip.to_string())
                .collect();
            out.sort();
            out
        };
        checked += 1;
        let before = drawn(observations);
        let after = drawn(swapped);
        let differences: Vec<(&String, &String)> = before
            .iter()
            .zip(after.iter())
            .filter(|(one, two)| one != two)
            .collect();
        assert_eq!(
            before, after,
            "seed {seed}: swapping two rows the panel refuses to join changed it: {differences:?}"
        );
    }
    assert!(
        checked > 1_000,
        "only {checked} panels had two rows to swap, so this proves little"
    );
}

#[test]
fn a_change_is_reported_when_there_is_something_to_measure_against() {
    // Growth is a difference between two observations, so a panel reporting one
    // is claiming there was an earlier observation to subtract, and the claim is
    // checkable in the document alone. Both directions matter. A change with
    // nothing behind it is a number invented out of one measurement, and a
    // missing change after two measurements the panel drew is an alert nobody
    // gets. What it does across an observation it refused to draw is left open
    // here, because that is a decision and not an invariant.
    let mut reported = 0u64;
    let mut required = 0u64;
    for seed in 0..ROUNDS {
        let mut rng = Lcg::new(seed);
        let window = region(&mut rng);
        let observations = surveillance_observations(&mut rng, &window);
        let panel = surveillance(&mut rng, observations);
        let svg = one_panel(&panel, &window).to_svg();

        // Per lineage: whether any earlier observation of it had a frequency, and
        // whether the one named just before was drawn and had one.
        let mut ever: BTreeMap<String, bool> = BTreeMap::new();
        let mut last: BTreeMap<String, bool> = BTreeMap::new();
        for (tooltip, _) in titled_groups(&svg) {
            let parts: Vec<&str> = tooltip.split(" | ").collect();
            let Some(at) = parts.iter().position(|part| part.starts_with("time ")) else {
                continue;
            };
            // The lineage is the part before the time, and anything before that
            // is the reason the panel gives for not drawing this one.
            let Some(lineage) = at.checked_sub(1).map(|before| parts[before].to_string()) else {
                continue;
            };
            let refused = at > 1;
            // The alert opens with the same word as the number it is about, and
            // only the number says the panel could read the two counts.
            let frequency = parts
                .iter()
                .any(|part| part.starts_with("frequency ") && !part.starts_with("frequency alert"));
            let change = parts.iter().any(|part| part.starts_with("change "));

            if change {
                reported += 1;
                assert!(
                    !refused && frequency,
                    "seed {seed}: a change on an observation the panel did not draw \
                     or could not read: {tooltip:?}"
                );
                assert_eq!(
                    ever.get(&lineage),
                    Some(&true),
                    "seed {seed}: a change measured against nothing: {tooltip:?}"
                );
            }
            if !refused && frequency && last.get(&lineage) == Some(&true) {
                required += 1;
                assert!(
                    change,
                    "seed {seed}: two drawn observations of one lineage and no change \
                     between them: {tooltip:?}"
                );
            }
            *ever.entry(lineage.clone()).or_default() |= frequency;
            last.insert(lineage, !refused && frequency);
        }
    }
    assert!(
        reported > 500 && required > 500,
        "only {reported} changes were reported and {required} were owed, \
         so this proves little"
    );
}

// ---------------------------------------------------------------------------
// Branch layers
// ---------------------------------------------------------------------------

/// The annotation keys the branch layers read.
const RATE_KEYS: [&str; 3] = ["omega1", "omega2", "omega3"];
const WEIGHT_KEYS: [&str; 3] = ["w1", "w2", "w3"];
const STATE_KEYS: [&str; 3] = ["p_A", "p_B", "p_C"];
const INTERVAL_KEYS: [&str; 3] = ["cf", "cf_lo", "cf_hi"];

/// A number that is a plausible reading two thirds of the time.
///
/// Hostile the rest of it, and not more, because every one of these layers
/// refuses to draw a value that is not finite: a generator that only produced
/// [`Lcg::value`] would spend most of its rounds checking that nothing is drawn.
fn measurement(rng: &mut Lcg, scale: f64) -> f64 {
    if rng.chance(3) {
        rng.value()
    } else {
        scale * (rng.below(1_001) as f64) / 1_000.0
    }
}

/// A short vocabulary of events, so the same one lands on several branches.
///
/// A homoplasy layer draws nothing until two branches agree, so the names have
/// to collide on purpose. The empty and the XML-hostile ones are in the
/// vocabulary because an ancestral-state reconstruction writes whatever it was
/// given.
fn event(rng: &mut Lcg) -> String {
    match rng.below(6) {
        0 => String::new(),
        1 => "a<b&c\"d'e".to_string(),
        2 => "\u{1f9ec}".to_string(),
        _ => format!("e{}", rng.below(3)),
    }
}

/// A tree whose branches carry what the interaction layers read.
///
/// Half the nodes get nothing. That is the case that matters: these layers
/// promise never to inherit, and a tree where every branch has its own value
/// cannot tell an inherited mark from a direct one.
fn interacting_tree(rng: &mut Lcg) -> Tree {
    let mut tree = tree(rng);
    for node in 0..tree.nodes().len() {
        if rng.chance(2) {
            continue;
        }
        let mut writes: Vec<(String, AnnotationValue)> = Vec::new();
        for index in 0..3 {
            if rng.chance(4) {
                continue;
            }
            writes.push((
                RATE_KEYS[index].to_string(),
                AnnotationValue::Number(measurement(rng, 5.0)),
            ));
            // A rate whose weight is missing, which is what a table one column
            // short looks like from inside the drawing.
            if rng.chance(6) {
                continue;
            }
            writes.push((
                WEIGHT_KEYS[index].to_string(),
                AnnotationValue::Number(measurement(rng, 1.0)),
            ));
        }
        for key in STATE_KEYS {
            if rng.chance(4) {
                continue;
            }
            writes.push((
                key.to_string(),
                AnnotationValue::Number(measurement(rng, 1.0)),
            ));
        }
        for key in INTERVAL_KEYS {
            if rng.chance(5) {
                continue;
            }
            writes.push((
                key.to_string(),
                AnnotationValue::Number(measurement(rng, 1.0)),
            ));
        }
        if !rng.chance(3) {
            writes.push((
                "omega".to_string(),
                AnnotationValue::Number(measurement(rng, 5.0)),
            ));
        }
        if !rng.chance(3) {
            writes.push((
                "p_value".to_string(),
                AnnotationValue::Number(measurement(rng, 1.0)),
            ));
        }
        // Every reading of an event annotation: one value of each kind, and the
        // brace-delimited list that becomes several ordered marks.
        if !rng.chance(3) {
            let value = match rng.below(4) {
                0 => AnnotationValue::Text(event(rng)),
                1 => AnnotationValue::Number(measurement(rng, 1.0)),
                2 => AnnotationValue::Boolean(rng.chance(2)),
                _ => AnnotationValue::List(
                    (0..rng.count())
                        .map(|_| match rng.below(3) {
                            0 => AnnotationValue::Number(measurement(rng, 1.0)),
                            1 => AnnotationValue::Boolean(rng.chance(2)),
                            _ => AnnotationValue::Text(event(rng)),
                        })
                        .collect(),
                ),
            };
            writes.push(("events".to_string(), value));
        }
        if !rng.chance(3) {
            writes.push(("event".to_string(), AnnotationValue::Text(event(rng))));
        }
        if let Some(annotations) = tree.annotations_mut(node) {
            for (key, value) in writes {
                annotations.insert(key, value);
            }
        }
    }
    tree
}

/// Every layer that hangs off a branch rather than standing on its own.
///
/// The extras come back with the track because they are the two things a
/// property cannot read out of the document: the interval domain is a caller's
/// choice that never reaches the page, and two mixtures on one branch are told
/// apart by their offsets rather than by their tooltips.
fn branch_layers(rng: &mut Lcg, mut track: TreeTrack) -> (TreeTrack, (f64, f64), usize) {
    let mut mixtures = 0usize;
    // The two key lists are the same length most of the time, because a
    // mismatched pair is refused whole and then draws nothing at all.
    if rng.chance(2) {
        let rates = 1 + rng.below(3) as usize;
        let weights = if rng.chance(4) {
            1 + rng.below(3) as usize
        } else {
            rates
        };
        let mut mixture = BranchRateMixture::new(
            RATE_KEYS[..rates].iter().copied(),
            WEIGHT_KEYS[..weights].iter().copied(),
        );
        if rng.chance(2) {
            mixture = mixture.label(rng.name());
        }
        if rng.chance(2) {
            mixture = mixture.width(rng.value());
        }
        if rng.chance(2) {
            mixture = mixture.thickness(rng.value());
        }
        if rng.chance(2) {
            mixture = mixture.neutral_band(rng.value(), rng.value());
        }
        if rng.chance(2) {
            mixture = mixture.saturation(rng.value());
        }
        track = track.branch_rate_mixture(mixture);
        mixtures += usize::from(rates == weights);
    }
    // A second one, since the layers stack by index and the offsets are the only
    // thing keeping two capsules off each other.
    if rng.chance(3) {
        track = track.branch_rate_mixture(BranchRateMixture::new(["omega1"], ["w1"]));
        mixtures += 1;
    }
    if rng.chance(2) {
        let mut layer = HomoplasyLayer::new(if rng.chance(6) {
            rng.name()
        } else {
            "event".to_string()
        });
        if rng.chance(2) {
            layer = layer.label(rng.name());
        }
        if rng.chance(2) {
            layer = layer.minimum_occurrences(rng.count());
        }
        if rng.chance(2) {
            layer = layer.maximum_connections(rng.count());
        }
        if rng.chance(2) {
            layer = layer.width(rng.value());
        }
        track = track.homoplasy_layer(layer);
    }
    if rng.chance(2) {
        let mut layer = BranchEventLayer::new(if rng.chance(6) {
            rng.name()
        } else {
            "events".to_string()
        });
        if rng.chance(2) {
            layer = layer.label(rng.name());
        }
        if rng.chance(2) {
            layer = layer.maximum_events(rng.count());
        }
        if rng.chance(2) {
            layer = layer.size(rng.value());
        }
        track = track.branch_event_layer(layer);
    }
    // One interval layer only, so the domain below belongs to every mark.
    let mut domain = (0.0, 1.0);
    if rng.chance(2) {
        let mut layer =
            BranchIntervalLayer::new(INTERVAL_KEYS[0], INTERVAL_KEYS[1], INTERVAL_KEYS[2]);
        if rng.chance(2) {
            layer = layer.label(rng.name());
        }
        if rng.chance(2) {
            // A domain whose own width overflows is what a caller writes when
            // they mean "any number at all".
            let (low, high) = (-rng.value(), rng.value());
            layer = layer.range(low, high);
            if low.is_finite() && high.is_finite() && high > low {
                domain = (low, high);
            }
        }
        if rng.chance(2) {
            layer = layer.threshold(rng.value());
        }
        if rng.chance(2) {
            layer = layer.width(rng.value());
        }
        track = track.branch_interval(layer);
    }
    if rng.chance(2) {
        let states = 1 + rng.below(3) as usize;
        let mut layer = AncestralStateLayer::new(STATE_KEYS[..states].iter().copied());
        if rng.chance(2) {
            layer = layer.label(rng.name());
        }
        if rng.chance(2) {
            layer = layer.confidence(rng.value());
        }
        if rng.chance(2) {
            layer = layer.size(rng.value());
        }
        if rng.chance(3) {
            layer = layer.show_transitions(false);
        }
        track = track.ancestral_states(layer);
    }
    if rng.chance(2) {
        track = track.dnds("omega");
        if rng.chance(2) {
            track = track.dnds_label(rng.name());
        }
        if rng.chance(2) {
            track = track.dnds_neutral_band(rng.value(), rng.value());
        }
        if rng.chance(2) {
            track = track.dnds_saturation(rng.value());
        }
        if rng.chance(2) {
            track = track.dnds_significance("p_value", rng.value());
        }
    }
    (track, domain, mixtures)
}

/// One drawing of one annotated tree, with what the document cannot say.
struct Interacting {
    figure: Figure,
    tree: Tree,
    interval_domain: (f64, f64),
    mixtures: usize,
}

fn interacting(seed: u64) -> Interacting {
    let mut rng = Lcg::new(seed);
    let window = region(&mut rng);
    let source = interacting_tree(&mut rng);
    let mut track = TreeTrack::new(source.clone());
    if rng.chance(2) {
        track = track.projection(match rng.below(3) {
            0 => TreeProjection::Rectangular,
            1 => TreeProjection::Circular,
            _ => TreeProjection::Unrooted,
        });
    }
    if rng.chance(2) {
        track = track.shape(if rng.chance(2) {
            TreeShape::Phylogram
        } else {
            TreeShape::Cladogram
        });
    }
    if rng.chance(3) {
        track = track.branch_geometry(match rng.below(3) {
            0 => BranchGeometry::Orthogonal,
            1 => BranchGeometry::Diagonal,
            _ => BranchGeometry::Curved,
        });
    }
    if rng.chance(4) {
        track = track.support_style(SupportStyle::Symbols);
    }
    let (track, interval_domain, mixtures) = branch_layers(&mut rng, track);
    // The row height and the branch width are left alone here, unlike in the
    // stack the whole-figure properties draw: those two decide the layout rather
    // than a layer, and a tree laid out a hundred million pixels tall has no
    // marks left on the page to hold to anything.
    let mut figure = Figure::new(window).show_region_label(false);
    if rng.chance(4) {
        figure = figure.width(match rng.below(4) {
            0 => 60.0,
            1 => 5_000.0,
            _ => 300.0 + rng.below(1_200) as f64,
        });
    }
    Interacting {
        figure: figure.push(track),
        tree: source,
        interval_domain,
        mixtures,
    }
}

#[test]
fn an_interval_mark_sits_where_its_number_says() {
    // The position of the dot is the claim the reader takes away, and the
    // tooltip beside it is the same claim in digits. A concordance factor, a
    // rate, a transition support: if the two disagree, every branch in the
    // figure is saying something the numbers do not, and nothing in the document
    // looks wrong.
    //
    // The fraction is worked out here from the domain the caller set, in
    // quarters, because the library's own `maximum - minimum` can overflow and
    // this must not overflow with it. A domain that wide is what somebody writes
    // when they mean "any number at all", and it is where a mark stops being
    // drawn at all: the writer refuses geometry it cannot write, which turns an
    // arithmetic mistake into a caption over blank paper. So the axis is what
    // says the mark was meant to be there.
    for seed in 0..ROUNDS {
        let drawing = interacting(seed);
        let (minimum, maximum) = drawing.interval_domain;
        let svg = drawing.figure.to_svg();
        for (title, content) in titled_groups(&svg) {
            let Some((_, rest)) = title.split_once(" | estimate ") else {
                continue;
            };
            let Some((estimate, _)) = rest.split_once(" | interval ") else {
                continue;
            };
            let Ok(estimate) = estimate.parse::<f64>() else {
                continue;
            };
            let Some(axis) = lines(content).first().copied() else {
                continue;
            };
            let dot = circles(content).first().copied().unwrap_or_else(|| {
                panic!("seed {seed}: {title:?} draws an axis and no estimate on it")
            });
            let wanted = ((estimate * 0.25 - minimum * 0.25) / (maximum * 0.25 - minimum * 0.25))
                .clamp(0.0, 1.0);
            let drawn = along(axis, dot)
                .unwrap_or_else(|| panic!("seed {seed}: {title:?} drew an axis of no length"));
            assert!(
                (drawn - wanted).abs() <= 0.02,
                "seed {seed}: an estimate of {estimate} on a {minimum:e} to {maximum:e} axis \
                 belongs {wanted} of the way along it and is drawn at {drawn}"
            );
        }
    }
}

#[test]
fn a_rate_mixture_draws_the_weights_it_prints() {
    // The capsule and its tooltip are two statements of one composition, and
    // this puts them side by side. The capsule normalises, so weights of 0.1 and
    // 0.1 draw the same picture as 0.5 and 0.5, and the tooltip is the only
    // place a reader can see that a mixture does not sum to one. That makes it
    // load-bearing rather than decorative.
    //
    // First it has to be a composition at all: every class the library drew had
    // a weight above zero, so a printed weight of nought is a number the reader
    // can neither use nor check. Then the segments, in order, against the shares
    // those weights imply, allowing for a class too small to reach the page.
    for seed in 0..ROUNDS {
        let svg = interacting(seed).figure.to_svg();
        for (title, content) in titled_groups(&svg) {
            if !title.contains(" | class 1 omega ") {
                continue;
            }
            let weights: Vec<f64> = title
                .split(" weight ")
                .skip(1)
                .filter_map(|piece| {
                    piece
                        .split(';')
                        .next()
                        .and_then(|number| number.trim().parse::<f64>().ok())
                })
                .collect();
            let total: f64 = weights.iter().sum();
            assert!(
                !weights.is_empty()
                    && weights.iter().all(|weight| *weight > 0.0)
                    && total.is_finite(),
                "seed {seed}: {title:?} prints weights {weights:?}, which cannot be read back"
            );
            let drawn = lines(content);
            let capsule = drawn[0];
            let span = (capsule.1 .0 - capsule.0 .0).hypot(capsule.1 .1 - capsule.0 .1);
            let segments: Vec<f64> = drawn[1..]
                .iter()
                .map(|(from, to)| (to.0 - from.0).hypot(to.1 - from.1))
                .collect();
            let tolerance = 0.02 * span.max(1.0) + 0.01;
            let covered: f64 = segments.iter().sum();
            assert!(
                (covered - span).abs() <= tolerance,
                "seed {seed}: {title:?} covers {covered} of a {span} capsule"
            );
            let mut segments = segments.iter().copied().peekable();
            for weight in &weights {
                let wanted = span * (weight / total);
                if wanted <= 1e-9 {
                    // Narrower than the page can hold: drawn as nothing, or
                    // swallowed by the running total, and both are honest.
                    if segments.peek().is_some_and(|drawn| *drawn <= 1e-6) {
                        segments.next();
                    }
                    continue;
                }
                let got = segments.next();
                assert!(
                    got.is_some_and(|got| (got - wanted).abs() <= tolerance),
                    "seed {seed}: {title:?} gives a weight of {weight} of {total} a segment of \
                     {got:?} rather than {wanted} of {span}"
                );
            }
            assert!(
                segments.next().is_none(),
                "seed {seed}: {title:?} draws more segments than it prints weights"
            );
        }
    }
}

#[test]
fn a_branch_layer_marks_only_the_branches_that_carry_its_data() {
    // The module says four times over that these values are direct and are never
    // inherited, and the reason is that a rate fitted on one branch is a
    // statement about that branch. An inherited value puts a fitted omega on
    // every descendant of the branch it was measured on, which is a claim nobody
    // made and which reads, on the page, exactly like a clade under selection.
    //
    // A count, not a match: a branch that carries the data can still go unmarked,
    // because it is too short to hold a capsule or because the bounds came back
    // the wrong way round. What cannot happen is a mark on a branch with nothing
    // on it, and counting is enough to catch that, since inheritance can only
    // ever add marks to branches that have none.
    for seed in 0..ROUNDS {
        let drawing = interacting(seed);
        let carrying = |keys: &[&str]| {
            (0..drawing.tree.nodes().len())
                .filter(|node| {
                    keys.iter()
                        .any(|key| drawing.tree.annotation(*node, key).is_some())
                })
                .count()
        };
        let intervals = carrying(&INTERVAL_KEYS[..1]);
        let rated = carrying(&RATE_KEYS);
        let svg = drawing.figure.to_svg();
        let groups = titled_groups(&svg);
        let count = |needle: &str| {
            groups
                .iter()
                .filter(|(title, _)| title.contains(needle))
                .count()
        };
        let drawn = count(" | estimate ");
        let capsules = count(" | class 1 omega ");
        assert!(
            drawn <= intervals,
            "seed {seed}: {drawn} interval marks over {intervals} branches that carry one"
        );
        assert!(
            capsules <= drawing.mixtures * rated,
            "seed {seed}: {capsules} capsules over {rated} branches that carry a rate, \
             from {} mixtures",
            drawing.mixtures
        );
    }
}

#[test]
fn every_event_drawn_on_a_branch_is_one_the_branch_carried() {
    // An event mark is a mutation, and its name is the whole of what it says. A
    // mark naming something no branch carries is a mutation the analysis never
    // reported, drawn where somebody will read it off and go looking for it. The
    // count catches the other direction: one branch's events copied onto its
    // descendants.
    //
    // The set is built from the tree rather than from the drawing, and the titles
    // are unescaped rather than the tree's strings escaped, so the two sides
    // never share a line of code.
    for seed in 0..ROUNDS {
        let drawing = interacting(seed);
        let mut carried: Vec<String> = Vec::new();
        for node in 0..drawing.tree.nodes().len() {
            match drawing.tree.annotation(node, "events") {
                Some(AnnotationValue::List(values)) => {
                    carried.extend(values.iter().map(ToString::to_string));
                }
                Some(other) => carried.push(other.to_string()),
                None => {}
            }
        }
        let svg = drawing.figure.to_svg();
        let mut drawn = 0usize;
        for (title, _) in titled_groups(&svg) {
            let Some((_, rest)) = title.split_once(" | events = ") else {
                continue;
            };
            drawn += 1;
            let event = unescaped(rest.split(" | ").next().unwrap_or(rest));
            assert!(
                carried.contains(&event),
                "seed {seed}: an event mark says {event:?}, which is on no branch of the tree"
            );
        }
        assert!(
            drawn <= carried.len(),
            "seed {seed}: {drawn} event marks for the {} events the tree carries",
            carried.len()
        );
    }
}

#[test]
fn a_transition_names_two_states_and_prints_two_fractions() {
    // The mark exists to say that the reconstructed state changed along this
    // branch. If the two names are the same there was no change and the figure
    // has invented one. If a posterior is above one it is not a posterior: it is
    // the raw number from before it was divided by the total, and a reader who
    // takes it for a probability will read an arbitrary weight as certainty.
    // Both halves come off the tooltip alone, with no knowledge of the tree, the
    // layer or the confidence threshold.
    for seed in 0..ROUNDS {
        let svg = interacting(seed).figure.to_svg();
        for (title, _) in titled_groups(&svg) {
            let Some((_, rest)) = title.split_once(" transition ") else {
                continue;
            };
            let Some((from, to)) = rest.split_once(" to ") else {
                continue;
            };
            let read = |piece: &str| -> Option<(String, f64)> {
                let (state, probability) = piece.rsplit_once(" (")?;
                let probability = probability.strip_suffix(')')?;
                Some((state.to_string(), probability.parse::<f64>().ok()?))
            };
            let (Some((was, before)), Some((now, after))) = (read(from), read(to)) else {
                continue;
            };
            assert!(
                was != now,
                "seed {seed}: {title:?} claims a change from a state to itself"
            );
            for probability in [before, after] {
                assert!(
                    (0.0..=1.0).contains(&probability),
                    "seed {seed}: {title:?} prints {probability} as a posterior"
                );
            }
        }
    }
}

#[test]
fn the_word_beside_an_omega_is_the_one_the_number_earns() {
    // This is the sentence a reader quotes. Everything else in these layers is a
    // position or a colour that somebody might squint at; the regime is a word,
    // in English, saying which way selection went on this branch, and it is the
    // one output of this crate that could be pasted into a paper unchanged.
    //
    // The band is the caller's, so the exact edges of "approximately neutral"
    // are not asserted here. Which side of one a named regime falls on is not
    // negotiable: the band is required to contain one at both ends, so nothing
    // purifying can sit above it and nothing diversifying below it. An unknown
    // word fails too, since a regime the reader has to guess at is worse than
    // none.
    for seed in 0..ROUNDS {
        let svg = interacting(seed).figure.to_svg();
        for (title, _) in titled_groups(&svg) {
            let Some((_, rest)) = title.split_once("dN/dS \u{3c9} ") else {
                continue;
            };
            let Some((number, regime)) = rest.split_once(" (") else {
                continue;
            };
            let Ok(omega) = number.parse::<f64>() else {
                continue;
            };
            let regime = regime.split(')').next().unwrap_or("");
            let wrong = match regime {
                "purifying" => omega >= 1.0,
                "diversifying" => omega <= 1.0,
                "approximately neutral" => false,
                _ => true,
            };
            assert!(
                !wrong,
                "seed {seed}: an omega of {omega} is called {regime:?}"
            );
        }
    }
}

#[test]
fn a_significant_branch_says_so_twice() {
    // Evidence and effect size are deliberately kept in separate channels here:
    // colour carries omega, weight carries the test. The weight is what a reader
    // sees from across the room and the tooltip is what they check, and the two
    // are computed by different expressions in different files over the same
    // threshold. Nothing but a property compares them, so a change to one of the
    // two comparisons would emphasise the branches that failed the test and
    // nobody would find out.
    //
    // A tie is skipped, and only a tie: a p-value of 0.0001 is printed as `0`,
    // so two numbers that print the same cannot be compared from the page.
    for seed in 0..ROUNDS {
        let svg = interacting(seed).figure.to_svg();
        let mut met: Vec<f64> = Vec::new();
        let mut missed: Vec<f64> = Vec::new();
        for (title, content) in titled_groups(&svg) {
            if !title.contains("dN/dS \u{3c9} ") {
                continue;
            }
            let Some((_, test)) = title.split_once("p_value ") else {
                continue;
            };
            let Some((value, rest)) = test.split_once(" (") else {
                continue;
            };
            let Some((relation, bound)) = rest.split_once(' ') else {
                continue;
            };
            let (Ok(value), Ok(bound)) = (
                value.parse::<f64>(),
                bound.trim_end_matches(')').parse::<f64>(),
            ) else {
                continue;
            };
            let passed = unescaped(relation) == "\u{2264}";
            assert!(
                value == bound || passed == (value <= bound),
                "seed {seed}: {title:?} puts {relation} between {value} and {bound}"
            );
            let Some(width) = number(content, "stroke-width") else {
                continue;
            };
            if passed {
                met.push(width);
            } else {
                missed.push(width);
            }
        }
        let thinnest = met.iter().copied().fold(f64::INFINITY, f64::min);
        let thickest = missed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            met.is_empty() || missed.is_empty() || thinnest > thickest,
            "seed {seed}: a branch that met the test is drawn {thinnest} wide and one that \
             missed it {thickest}"
        );
    }
}
