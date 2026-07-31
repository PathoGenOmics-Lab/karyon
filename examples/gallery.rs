//! Every kind of plot this crate draws, on one sheet.
//!
//! ```text
//! cargo run --example gallery -- assets
//! ```
//!
//! Keep it that way: **when a new track type is added, give it a panel here.**
//! An overview that quietly stops covering everything is worse than no
//! overview, because it looks complete.
//!
//! The panels do not share a coordinate system with each other, which is the
//! honest reason they are separate figures rather than one tall stack. Three of
//! them are not even in genomic coordinates: the alignment counts columns, the
//! variable site panel counts sites, and the tree measures evolutionary
//! distance.

use std::env;
use std::path::PathBuf;

use karyon::tree::Tree;
use karyon::{
    Aggregate, AlignmentBlock, Association, AxisTrack, Band, CellScale, CigarOp, CoverageTrack,
    DotplotTrack, Feature, FeatureTrack, Figure, IdeogramTrack, LogoColumn, LogoScore, LogoTrack,
    ManhattanTrack, MatrixRow, MatrixTrack, MsaColoring, MsaDisplay, MsaSequence, MsaTrack, Panels,
    PileupTrack, Read, ReadColoring, Region, SequenceTrack, SnpTrack, Stain, Strand, SyntenyTrack,
    TreeTrack, Variant, VariantTrack,
};

/// Width every panel is drawn at.
const WIDTH: f64 = 860.0;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sheet = Panels::new()
        .title("karyon: every representation")
        .gap(20.0)
        .push_captioned(
            &locus(),
            "A",
            "Ideogram, coverage, sequence, genes, variants, ruler",
        )
        .push_captioned(
            &pileup(),
            "B",
            "Read pileup with mismatches against the reference",
        )
        .push_captioned(
            &logos(),
            "C",
            "Sequence logos: information content and enrichment/depletion",
        )
        .push_captioned(
            &association(),
            "D",
            "Association statistics and a genotype matrix",
        )
        .push_captioned(
            &comparison(),
            "E",
            "Dotplot and synteny ribbons between two sequences",
        )
        .push_captioned(
            &alignment(),
            "F",
            "Multiple sequence alignment, differences only",
        )
        .push_captioned(
            &variable_sites(),
            "G",
            "Variable sites with the phylogeny that orders them",
        )
        .push_captioned(&phylogeny(), "H", "A tree on its own");

    sheet.save_svg(out.join("gallery.svg"))?;
    let (width, height) = sheet.dimensions();
    println!(
        "gallery.svg {width:.0} x {height:.0}, {} panels",
        sheet.len()
    );
    Ok(())
}

/// A: the genomic stack.
fn locus() -> Figure {
    let start = 761_000u64;
    let span = 2_000usize;
    let depth: Vec<f64> = (0..span)
        .map(|i| 48.0 + 12.0 * ((i as f64) / 140.0).sin())
        .collect();
    let bases: Vec<u8> = b"ACGT".iter().cycle().take(span).copied().collect();

    Figure::new(Region::new("NC_000962.3", start, start + span as u64).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            IdeogramTrack::new(
                4_411_532,
                vec![
                    Band::new(0, 1_800_000, Stain::Gneg),
                    Band::new(1_800_000, 2_600_000, Stain::Gpos50),
                    Band::new(2_600_000, 4_411_532, Stain::Gneg),
                ],
            )
            .label("H37Rv")
            .height(18.0),
        )
        .push(
            CoverageTrack::new(start, depth)
                .aggregate(Aggregate::Min)
                .label("depth")
                .height(46.0),
        )
        .push(SequenceTrack::new(start, bases).label("reference"))
        .push(
            FeatureTrack::new(vec![
                Feature::new(761_050, 762_100)
                    .name("rpoB")
                    .strand(Strand::Forward),
                Feature::new(762_250, 762_800)
                    .name("Rv0668c")
                    .strand(Strand::Reverse),
            ])
            .label("genes"),
        )
        .push(
            VariantTrack::new(vec![
                Variant::new(761_109).value(0.98).category("missense"),
                Variant::new(761_452).value(0.44).category("synonymous"),
            ])
            .label("variants")
            .height(40.0),
        )
        .push(AxisTrack::new())
}

/// B: reads.
fn pileup() -> Figure {
    let start = 4_000u64;
    let span = 260usize;
    let reference = vec![b'A'; span];
    let reads: Vec<Read> = (0..24u64)
        .map(|i| {
            let at = start + (i * 7) % 150;
            let mut sequence = vec![b'A'; 100];
            // A variant at reference 4,120, in the reads that reach it.
            if i % 3 == 0 && (at..at + 100).contains(&4_120) {
                sequence[(4_120 - at) as usize] = b'G';
            }
            Read::new(at, vec![CigarOp::Match(100)])
                .sequence(sequence)
                .strand(if i % 2 == 0 {
                    Strand::Forward
                } else {
                    Strand::Reverse
                })
                .mapping_quality(60)
        })
        .collect();

    Figure::new(Region::new("NC_000962.3", start, start + span as u64).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            PileupTrack::new(reads)
                .reference(start, reference)
                .coloring(ReadColoring::Strand)
                .max_rows(Some(12))
                .label("reads"),
        )
        .push(AxisTrack::new())
}

/// C: logos.
fn logos() -> Figure {
    let motif = vec![
        LogoColumn::acgt(97.0, 1.0, 1.0, 1.0),
        LogoColumn::acgt(48.0, 46.0, 3.0, 3.0),
        LogoColumn::acgt(26.0, 25.0, 24.0, 25.0),
        LogoColumn::acgt(34.0, 33.0, 33.0, 0.0),
        LogoColumn::acgt(8.0, 9.0, 74.0, 9.0),
        LogoColumn::acgt(2.0, 44.0, 10.0, 44.0),
        LogoColumn::acgt(40.0, 30.0, 20.0, 10.0),
        LogoColumn::acgt(32.0, 4.0, 32.0, 32.0),
    ];
    Figure::new(Region::new("motif", 0, motif.len() as u64).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            LogoTrack::new(0, motif.clone())
                .alphabet_size(4)
                .score(LogoScore::InformationContent)
                .label("bits")
                .height(58.0),
        )
        .push(
            LogoTrack::new(0, motif)
                .alphabet_size(4)
                .edlogo()
                .label("enrich / deplete")
                .height(96.0),
        )
        .push(AxisTrack::new().center_on_bases(true))
}

/// D: association and genotypes.
fn association() -> Figure {
    let start = 759_000u64;
    let span = 6_000u64;
    let peak = 761_155u64;
    let mut rng = Lcg::new(4_242);

    let mut points = Vec::new();
    let mut sites = Vec::new();
    let mut position = start + 40;
    while position < start + span {
        let linkage = (-(position.abs_diff(peak) as f64) / 260.0).exp();
        let noise = (rng.next() % 1000) as f64 / 1000.0;
        let value = 0.3 + 2.4 * noise + 9.5 * linkage * (0.55 + 0.45 * noise);
        points.push(Association::new(position, value));
        if value > 3.0 && sites.len() < 12 {
            sites.push(position);
        }
        position += 30 + rng.next() % 40;
    }

    let rows: Vec<MatrixRow> = (0..6)
        .map(|isolate| {
            let carrier = isolate % 3 != 2;
            let values: Vec<f64> = sites
                .iter()
                .map(|_| {
                    if carrier && rng.next() % 100 < 88 {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect();
            MatrixRow::new(format!("ERR{:04}", 3_100 + isolate * 7), values)
        })
        .collect();

    Figure::new(Region::new("NC_000962.3", start, start + span).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            ManhattanTrack::new(points)
                .genome_wide_threshold()
                .unit(" -log10 p")
                .label("association")
                .height(74.0),
        )
        .push(
            MatrixTrack::new(sites, rows)
                .scale(CellScale::Sequential {
                    max: Some(1.0),
                    hue: None,
                })
                .min_cell_width(9.0)
                .label("genotypes"),
        )
        .push(AxisTrack::new())
}

/// E: two sequences.
fn comparison() -> Figure {
    let blocks = vec![
        AlignmentBlock::new(0, 1_500_000, 0, 1_500_000),
        AlignmentBlock::new(1_520_000, 2_100_000, 1_520_000, 2_100_000).reversed(true),
        AlignmentBlock::new(2_150_000, 2_600_000, 3_400_000, 3_850_000),
        AlignmentBlock::new(2_900_000, 4_400_000, 2_150_000, 3_650_000),
    ];
    Figure::new(Region::new("H37Rv", 0, 4_400_000).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            DotplotTrack::new(blocks.clone())
                .target_length(4_380_000)
                .label("CDC1551")
                .height(150.0),
        )
        .push(
            SyntenyTrack::new(blocks)
                .target_length(4_380_000)
                .names("H37Rv", "CDC1551")
                .height(92.0),
        )
        .push(AxisTrack::new())
}

/// F: an alignment.
fn alignment() -> Figure {
    let rows = vec![
        MsaSequence::new("KatG_H37Rv", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_CDC1551", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_S315T", b"MPEQHPPITETTTGATSNGCPV".to_vec()),
        MsaSequence::new("KatG_Beijing", b"MPEQHPPVTETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_Erdman", b"MPEQ-PPITETTTGAASNGCPV".to_vec()),
    ];
    Figure::new(Region::new("alignment", 0, 22).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            MsaTrack::new(rows)
                .display(MsaDisplay::Bases)
                .coloring(MsaColoring::Residue)
                .compare_to(0)
                .label("KatG")
                .row_height(14.0),
        )
        .push(AxisTrack::new().center_on_bases(true))
}

/// G: variable sites, ordered by a tree.
fn variable_sites() -> Figure {
    let mut rng = Lcg::new(1_848);
    let columns = 4_000usize;
    let reference: Vec<u8> = (0..columns)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();

    let mut rows = vec![MsaSequence::new("H37Rv", reference.clone())];
    for isolate in 0..8 {
        let mut row = reference.clone();
        if isolate % 3 == 1 {
            for site in [310usize, 1_020, 2_440] {
                row[site] = b"ACGT"[(rng.next() % 4) as usize];
            }
        }
        for _ in 0..2 {
            let site = (rng.next() as usize) % columns;
            row[site] = b"ACGT"[(rng.next() % 4) as usize];
        }
        rows.push(MsaSequence::new(format!("ERR{:04}", 5_001 + isolate), row));
    }

    let tree = Tree::parse_newick(
        "(((ERR5002:0.002,ERR5005:0.002)0.98:0.004,ERR5008:0.005):0.006,\
          ((ERR5001:0.004,ERR5003:0.003):0.002,(ERR5004:0.004,\
           (ERR5006:0.003,ERR5007:0.003):0.002):0.003):0.004);",
    )
    .expect("the tree in this example is well formed");

    let panel = SnpTrack::from_alignment(0, &rows)
        .tree(tree)
        .label("isolates");
    Figure::new(Region::new("sites", 0, panel.sites().len() as u64).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(panel)
}

/// H: a tree by itself.
fn phylogeny() -> Figure {
    let tree = Tree::parse_newick(
        "((lineage_1:0.012,lineage_2:0.011)0.99:0.006,\
          (lineage_3:0.009,(lineage_4:0.007,lineage_4.9:0.008)0.95:0.003):0.005);",
    )
    .expect("the tree in this example is well formed");

    Figure::new(Region::new("tree", 0, 1).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(TreeTrack::new(tree).label("phylogeny").row_height(17.0))
}

/// A linear congruential generator, so the sheet is reproducible without a
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
