//! Four figures about reading a coordinate as something other than a base.
//!
//! ```text
//! cargo run --example reading -- assets
//! ```
//!
//! As a residue, as a place one molecule visited, as a stretch a whole clade
//! shares, and as part of one RNA.

use std::env;
use std::path::{Path, PathBuf};

use karyon::tree::Tree;
use karyon::{
    AxisTrack, CladeBlock, CladeTrack, CodonTrack, CoverageTrack, Feature, FeatureTrack, Figure,
    Region, SplitRead, SplitReadTrack, SplitSegment, Strand, Terminator, TranscriptionUnit,
    TranscriptionUnitTrack, Variant, VariantTrack,
};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    codons(&out)?;
    split_reads(&out)?;
    clades(&out)?;
    transcripts(&out)?;
    Ok(())
}

/// The rpoB resistance determining region, in the coordinates it is named in.
fn codons(out: &Path) -> std::io::Result<()> {
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
        .sequence(view_from, bases.clone())
        .label("rpoB");
    let s450 = ruler.span_of(450).expect("rpoB has a codon 450");

    let figure = Figure::new(Region::new("NC_000962.3", view_from, view_to).unwrap())
        .title("The rpoB RRDR, in the coordinates the clinic names it in")
        .width(880.0)
        .show_region_label(false)
        .push(
            VariantTrack::new(vec![
                Variant::new(s450.0 + 1).category("S450L"),
                Variant::new(ruler.span_of(445).unwrap().0 + 1).category("H445Y"),
            ])
            .label("variants")
            .height(34.0),
        )
        .push(ruler)
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-codons.svg"))?;
    let (width, height) = figure.dimensions();
    println!(
        "example-codons.svg {width:.0} x {height:.0}, S450 at {}..{}",
        s450.0, s450.1
    );
    Ok(())
}

/// An IS6110 transposition, which is a three-segment read and not an arc.
fn split_reads(out: &Path) -> std::io::Result<()> {
    let mut rng = Lcg::new(8_101);
    let span = 40_000u64;
    let donor = (12_000u64, 13_355u64);
    let landing = 29_400u64;

    // Six molecules crossing the new junction, and two crossing the empty donor
    // site, which is the pair of observations that makes it a transposition
    // rather than a duplication.
    let mut reads: Vec<SplitRead> = (0..6)
        .map(|index| {
            let jitter = rng.next() % 600;
            let flip = index % 3 == 0;
            SplitRead::new(vec![
                SplitSegment::new(landing - 1_400 - jitter, landing, Strand::Forward)
                    .read_span(0, 1_400 + jitter),
                SplitSegment::new(
                    donor.0,
                    donor.1,
                    if flip {
                        Strand::Reverse
                    } else {
                        Strand::Forward
                    },
                )
                .read_span(1_400 + jitter, 2_755 + jitter),
                SplitSegment::new(landing, landing + 900 + jitter, Strand::Forward)
                    .read_span(2_755 + jitter, 3_655 + 2 * jitter),
            ])
            .name(format!("read_{:02}", index + 1))
        })
        .collect();
    for index in 0..2 {
        let jitter = rng.next() % 400;
        reads.push(
            SplitRead::new(vec![
                SplitSegment::new(donor.0 - 1_100 - jitter, donor.0, Strand::Forward)
                    .read_span(0, 1_100 + jitter),
                SplitSegment::new(donor.1, donor.1 + 1_000, Strand::Forward)
                    .read_span(1_100 + jitter, 2_100 + jitter),
            ])
            .name(format!("read_{:02}", index + 7)),
        );
    }

    let track = SplitReadTrack::new(reads).label("reads").row_height(11.0);
    let backwards = track
        .reads()
        .iter()
        .filter(|read| read.goes_backwards())
        .count();

    // Depth that agrees: the donor site keeps its copy, so nothing drops out.
    let depth: Vec<f64> = (0..span)
        .map(|at| {
            let base = if at >= donor.0 && at < donor.1 {
                96.0
            } else {
                48.0
            };
            base + (rng.next() % 60) as f64 / 10.0
        })
        .collect();

    let figure = Figure::new(Region::new("NC_000962.3", 0, span).unwrap())
        .title("An IS6110 copy in a new place: three segments of one molecule")
        .width(880.0)
        .show_region_label(false)
        .push(track)
        .push(CoverageTrack::new(0, depth).label("depth").height(48.0))
        .push(
            FeatureTrack::new(vec![Feature::new(donor.0, donor.1)
                .name("IS6110")
                .strand(Strand::Forward)])
            .label("donor"),
        )
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-split.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-split.svg {width:.0} x {height:.0}, {backwards} reads run backwards");
    Ok(())
}

/// Regions of difference painted onto the lineage tree that carries them.
fn clades(out: &Path) -> std::io::Result<()> {
    let tree = Tree::parse_newick(
        "(((M_canettii:0.09,\
            (((L1_ERR1:0.011,L1_ERR2:0.010):0.014,\
              (L2_ERR3:0.009,L2_ERR4:0.011):0.013):0.008,\
             ((L3_ERR5:0.010,L4_ERR6:0.012):0.011,\
              (L4_ERR7:0.009,L4_ERR8:0.010):0.012):0.009):0.030):0.020,\
           M_bovis:0.055):0.010,BCG_Pasteur:0.062);",
    )
    .expect("the tree in this example is well formed");

    let blocks = vec![
        // TbD1 is missing from every modern lineage and present in the ancestral
        // ones, which is why it is a block on an internal branch.
        CladeBlock::new(
            1_779_000,
            1_788_000,
            [
                "L1_ERR1", "L1_ERR2", "L2_ERR3", "L2_ERR4", "L3_ERR5", "L4_ERR6", "L4_ERR7",
                "L4_ERR8",
            ],
        )
        .name("TbD1 deleted"),
        // RD1 is gone from BCG alone, which is what made it a vaccine.
        CladeBlock::new(4_350_000, 4_359_600, ["BCG_Pasteur"]).name("RD1"),
        // RD9 marks the M. bovis and BCG side.
        CladeBlock::new(2_330_000, 2_342_000, ["M_bovis", "BCG_Pasteur"]).name("RD9"),
        // A block two lineages carry and the one between them does not: not a
        // clade, so the figure cuts the rows out rather than claiming one event.
        CladeBlock::new(
            3_100_000,
            3_128_000,
            ["L1_ERR1", "L1_ERR2", "L4_ERR7", "L4_ERR8"],
        )
        .name("recurrent"),
    ];

    let track = CladeTrack::new(tree, blocks)
        .label("isolates")
        .row_height(12.0)
        .tree_width(120.0);
    let cut = track.cut_rows(3);
    let clades = (0..track.blocks().len())
        .filter(|index| track.is_clade(*index))
        .count();

    let figure = Figure::new(Region::new("NC_000962.3", 0, 4_411_532).unwrap())
        .title("Regions of difference, painted onto the branch that lost them")
        .width(880.0)
        .show_region_label(false)
        .push(track)
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-clades.svg"))?;
    let (width, height) = figure.dimensions();
    println!(
        "example-clades.svg {width:.0} x {height:.0}, {clades} of 4 are clades, {cut} rows cut out"
    );
    Ok(())
}

/// One RNA molecule at a time, and where translation starts on it.
fn transcripts(out: &Path) -> std::io::Result<()> {
    // The ESX-1 locus. esxB and esxA are one transcript, and leaderless
    // transcription is common enough in this organism to be the point of the
    // figure rather than a footnote.
    let units = vec![
        TranscriptionUnit::new(4_350_200, 4_352_400, Strand::Forward)
            .cds_start(4_350_200)
            .terminator(Terminator::Intrinsic)
            .name("esxB-esxA"),
        TranscriptionUnit::new(4_352_900, 4_355_600, Strand::Forward)
            .cds_start(4_353_064)
            .terminator(Terminator::Intrinsic)
            .name("espI"),
        TranscriptionUnit::new(4_358_900, 4_356_100, Strand::Reverse)
            .cds_start(4_358_820)
            .terminator(Terminator::RhoDependent)
            .name("eccCa1"),
        TranscriptionUnit::new(4_355_800, 4_356_400, Strand::Forward).name("unresolved 3' end"),
    ];

    let track = TranscriptionUnitTrack::new(units)
        .label("transcripts")
        .row_height(28.0);
    let leaderless = track.leaderless();

    let genes = FeatureTrack::new(vec![
        Feature::new(4_350_200, 4_350_500)
            .name("esxB")
            .strand(Strand::Forward),
        Feature::new(4_350_600, 4_350_900)
            .name("esxA")
            .strand(Strand::Forward),
        Feature::new(4_353_064, 4_355_400)
            .name("espI")
            .strand(Strand::Forward),
        Feature::new(4_356_100, 4_358_820)
            .name("eccCa1")
            .strand(Strand::Reverse),
    ])
    .label("genes");

    let figure = Figure::new(Region::new("NC_000962.3", 4_349_900, 4_359_200).unwrap())
        .title("Transcription units: one arrow, one molecule, one hairpin")
        .width(880.0)
        .show_region_label(false)
        .push(track)
        .push(genes)
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-transcripts.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-transcripts.svg {width:.0} x {height:.0}, {leaderless} leaderless");
    Ok(())
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
