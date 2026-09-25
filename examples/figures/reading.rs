//! The figures `cargo run --example reading` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever window the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::tree::Tree;
use karyon::{
    CladeBlock, CladeTrack, CodonTrack, Drawing, Feature, Plot, Region, SplitRead, SplitReadTrack,
    SplitSegment, Strand, Terminator, Theme, TranscriptionUnit, TranscriptionUnitTrack, Variant,
};

/// `example-codons.svg`: a coding sequence in the coordinates its variants are
/// named in.
///
/// `theme`, `width` and `region` replace the light theme, the 880 pixels and
/// the eighty bases the committed figure is drawn over.
pub fn example_codons(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let mut rng = Lcg::new(4_493);
    // rpoB runs 759,807 to 763,325 in H37Rv, 1-based inclusive, forward strand.
    let cds_start = 759_806u64;
    let cds_end = 763_325u64;
    let view_from = 761_120u64;
    let view_to = 761_200u64;

    // A stretch of coding sequence with no stop in frame, since a stop inside a
    // gene would be a nonsense the figure never has to draw.
    let bases: Vec<u8> = (0..(view_to - view_from + 6))
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    let mut bases = bases;
    let ruler = CodonTrack::new(cds_start, cds_end, Strand::Forward);
    for codon in 1..=ruler.codons() {
        let Some((from, to)) = ruler.span_of(codon) else {
            continue;
        };
        if to <= view_from || from >= view_to + 6 {
            continue;
        }
        let at = (from - view_from) as usize;
        if at + 3 > bases.len() {
            continue;
        }
        // Serine at 450, and something that is not a stop everywhere else.
        let triplet: &[u8; 3] = match codon {
            450 => b"TCG",
            445 => b"CAC",
            435 => b"GAC",
            _ => b"CTG",
        };
        bases[at..at + 3].copy_from_slice(triplet);
    }

    let ruler = CodonTrack::new(cds_start, cds_end, Strand::Forward)
        .sequence(view_from, bases)
        .label("rpoB");
    let s450 = ruler.span_of(450).expect("rpoB has a codon 450");

    let own = Region::new("NC_000962.3", view_from, view_to).unwrap();
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("A coding sequence read as protein: rpoB codons 439 to 465")
        .theme(theme.clone())
        .width(width.unwrap_or(880.0))
        .remove_region_label()
        .add_variants(vec![
            Variant::new(s450.0 + 1).category("S450L"),
            Variant::new(ruler.span_of(445).unwrap().0 + 1).category("H445Y"),
        ])
        .label("variants")
        .adjust(|track| track.height(34.0))
        .add_track(ruler)
        .into_figure();
    Box::new(figure)
}

/// `example-split.svg`: an IS6110 transposition, which is a three-segment read
/// and not an arc.
///
/// `theme`, `width` and `region` replace the light theme, the 880 pixels and
/// the forty kilobases the committed figure is drawn over.
pub fn example_split(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    let mut rng = Lcg::new(8_101);
    // IS6110-1 of M. tuberculosis H37Rv, 1-based 889,021 to 890,375, 1,355 bp,
    // forward strand. The reference carries this copy; the sample carries it
    // and one more.
    let from = 880_000u64;
    let span = 40_000u64;
    let donor = (889_020u64, 890_375u64);
    let landing = 906_000u64;

    // Eight molecules crossing the new junction. Every one of them shows the
    // element in the same orientation, because one insertion event puts it in
    // one orientation: reads that disagreed about that would be two events or a
    // mapping artefact, not evidence for this one. Here the new copy went in
    // backwards relative to the reference copy, so the middle segment aligns to
    // the donor on the reverse strand and both connectors run backwards.
    let reads: Vec<SplitRead> = (0..8)
        .map(|index| {
            let left = 900 + rng.next() % 900;
            let right = 800 + rng.next() % 1_000;
            let element = donor.1 - donor.0;
            SplitRead::new(vec![
                SplitSegment::new(landing - left, landing, Strand::Forward).read_span(0, left),
                SplitSegment::new(donor.0, donor.1, Strand::Reverse)
                    .read_span(left, left + element),
                SplitSegment::new(landing, landing + right, Strand::Forward)
                    .read_span(left + element, left + element + right),
            ])
            .name(format!("read_{:02}", index + 1))
        })
        .collect();

    let track = SplitReadTrack::new(reads).label("reads").row_height(11.0);

    // Depth that agrees, and this is the half of the evidence the arcs cannot
    // carry: the sample has two copies of an element the reference has once, so
    // reads from both pile onto the one reference copy and the donor reads at
    // twice the background. A figure that also showed reads crossing an empty
    // donor site would be claiming the element had left, which is the opposite
    // of what this depth says.
    let depth: Vec<f64> = (from..from + span)
        .map(|at| {
            let base = if at >= donor.0 && at < donor.1 {
                96.0
            } else {
                48.0
            };
            base + (rng.next() % 60) as f64 / 10.0
        })
        .collect();

    let own = Region::new("NC_000962.3", from, from + span).unwrap();
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("An IS6110 copy in a new place: three segments of one molecule")
        .theme(theme.clone())
        .width(width.unwrap_or(880.0))
        .remove_region_label()
        .add_track(track)
        // From the first base the depth was measured at rather than from the
        // left edge of the window, which is only the same place until the
        // window moves.
        .add_coverage_at(from, depth)
        .label("depth")
        .adjust(|track| track.height(48.0))
        .add_features(vec![Feature::new(donor.0, donor.1)
            .name("IS6110")
            .strand(Strand::Forward)])
        .label("donor")
        .into_figure();
    Box::new(figure)
}

/// `example-clades.svg`: deletions painted onto the branch that lost them.
///
/// `theme`, `width` and `region` replace the light theme, the 880 pixels and
/// the whole genome the committed figure is drawn over.
pub fn example_clades(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // SARS-CoV-2 rather than a clonal bacterium, because the point of this
    // panel is a block that was lost more than once, and a recurrent deletion
    // has to come from a lineage that actually produced one.
    // The variants of concern are separate descendants of B.1 with no
    // well resolved order among them, so they are drawn as the polytomy they
    // are rather than given an invented branching order.
    let tree = Tree::parse_newick(
        "(A_Wuhan:0.002,(B_1_1_7_Alpha:0.005,\
          (B_1_351_Beta:0.004,P_1_Gamma:0.004):0.002,\
          B_1_617_2_Delta:0.005,\
          (BA_1_Omicron:0.008,BA_2_Omicron:0.008):0.003):0.003);",
    )
    .expect("the tree in this example is well formed");

    // 1-based coordinates on NC_045512.2, all lineage-defining.
    let blocks = vec![
        // The nsp6 SGF deletion: lost four times over, and Delta sits between
        // the carriers without carrying it. That is why the block is drawn with
        // Delta's row cut out rather than as one ancestral event.
        CladeBlock::new(
            11_287,
            11_296,
            [
                "B_1_1_7_Alpha",
                "B_1_351_Beta",
                "P_1_Gamma",
                "BA_1_Omicron",
                "BA_2_Omicron",
            ],
        )
        .name("nsp6 SGF del, recurrent"),
        // HV69-70, the deletion behind S gene target failure on the TaqPath
        // assay. Also recurrent: Alpha and BA.1 but not BA.2 between them.
        CladeBlock::new(21_764, 21_770, ["B_1_1_7_Alpha", "BA_1_Omicron"]).name("S HV69-70, SGTF"),
        // A clean single-branch block, the kind that is one event.
        CladeBlock::new(22_280, 22_289, ["B_1_351_Beta"]).name("S LAL242-244"),
    ];

    let track = CladeTrack::new(tree, blocks)
        .label("lineages")
        .row_height(14.0)
        .tree_width(150.0);

    let own = Region::parse("NC_045512.2:1-29,903").expect("the locus in this example parses");
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("Lineage-defining deletions, painted onto the branch that lost them")
        .theme(theme.clone())
        .width(width.unwrap_or(880.0))
        .remove_region_label()
        .add_track(track)
        .into_figure();
    Box::new(figure)
}

/// `example-transcripts.svg`: one RNA molecule at a time, and where
/// translation starts on it.
///
/// `theme`, `width` and `region` replace the light theme, the 880 pixels and
/// the ESX-1 locus the committed figure is drawn over.
pub fn example_transcripts(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // The ESX-1 locus of M. tuberculosis H37Rv, at its real coordinates. Locus
    // tags carry a trailing c only when the gene is on the complementary
    // strand, so Rv3879c is the reverse strand gene here and everything from
    // Rv3871 to Rv3878 runs forward.
    let units = vec![
        // esxB and esxA are one transcript, and it is leaderless: transcription
        // starts on the start codon itself. That is common enough in this
        // organism to be the point of the figure rather than a footnote.
        TranscriptionUnit::new(4_352_272, 4_352_945, Strand::Forward)
            .cds_start(4_352_272)
            .terminator(Terminator::Intrinsic)
            .name("esxB-esxA"),
        TranscriptionUnit::new(4_352_950, 4_355_120, Strand::Forward)
            .cds_start(4_353_008)
            .terminator(Terminator::Intrinsic)
            .name("espI"),
        // Rv3879c, read the other way, so its start site is the high coordinate.
        TranscriptionUnit::new(4_359_850, 4_357_591, Strand::Reverse)
            .cds_start(4_359_781)
            .terminator(Terminator::RhoDependent)
            .name("espK"),
    ];

    let track = TranscriptionUnitTrack::new(units)
        .label("transcripts")
        .row_height(28.0);

    // 1-based inclusive in the comments, 0-based half-open in the call.
    let genes = vec![
        Feature::new(4_348_825, 4_350_602) // Rv3871, eccCb1
            .name("eccCb1")
            .strand(Strand::Forward),
        Feature::new(4_350_743, 4_351_044) // Rv3872
            .name("PE35")
            .strand(Strand::Forward),
        Feature::new(4_351_073, 4_352_181) // Rv3873
            .name("PPE68")
            .strand(Strand::Forward),
        Feature::new(4_352_272, 4_352_576) // Rv3874
            .name("esxB")
            .strand(Strand::Forward),
        Feature::new(4_352_607, 4_352_896) // Rv3875
            .name("esxA")
            .strand(Strand::Forward),
        Feature::new(4_353_008, 4_355_010) // Rv3876
            .name("espI")
            .strand(Strand::Forward),
        Feature::new(4_355_005, 4_356_542) // Rv3877
            .name("eccD1")
            .strand(Strand::Forward),
        Feature::new(4_356_691, 4_357_535) // Rv3878
            .name("espJ")
            .strand(Strand::Forward),
        Feature::new(4_357_591, 4_359_782) // Rv3879c, the reverse one
            .name("espK")
            .strand(Strand::Reverse),
    ];

    let own =
        Region::parse("NC_000962.3:4,348,701-4,360,100").expect("the locus in this example parses");
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("Transcription units: one arrow, one molecule, one hairpin")
        .theme(theme.clone())
        .width(width.unwrap_or(880.0))
        .remove_region_label()
        .add_track(track)
        .add_features(genes)
        .label("genes")
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
