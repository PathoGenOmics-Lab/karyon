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
//! And the test a new track has to pass to earn a panel: **does it live on the
//! genomic coordinate axis?** That is the whole reason this crate exists rather
//! than a general plotting library. A track whose `draw` never reads
//! `ctx.scale` is a bar chart, a line chart or a heatmap that happens to have
//! been handed genomic data, and three of those were removed from this sheet
//! for exactly that reason.
//!
//! The panels do not share a coordinate system with each other, which is the
//! honest reason they are separate figures rather than one tall stack. Four of
//! them are not in genomic coordinates at all: the alignment counts columns,
//! the variable site panel counts sites, and the two tree panels measure
//! evolutionary distance.

use std::env;
use std::path::PathBuf;

use karyon::theme::mix;
use karyon::tree::Tree;
use karyon::{
    Aggregate, AlignmentBlock, Association, AxisRing, AxisTrack, Band, BisulfiteTrack, CellScale,
    CigarOp, CoverageTrack, DotplotTrack, Feature, FeatureRing, FeatureTrack, Figure, Genome,
    GenomeTrack, Homology, IdeogramTrack, Legend, LegendTrack, Locus, LocusTrack, LogoColumn,
    LogoScore, LogoTrack, ManhattanTrack, MarkerRing, MatrixRow, MatrixTrack, MethylSite,
    MethylationTrack, Molecule, Move, MsaColoring, MsaDisplay, MsaSequence, MsaTrack, OrfTrack,
    Panels, PileupTrack, Read, ReadColoring, Region, Rings, SequenceTrack, SignalRing, SnpTrack,
    SquiggleTrack, Stain, Strand, StructuralTrack, StructuralVariant, SvKind, SyntenyTrack,
    TanglegramTrack, Theme, TreeTrack, Variant, VariantTrack, Window, WindowStyle, WindowTrack,
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
        // Fifteen panels in one column is a scroll rather than a figure.
        .columns(3)
        .gap(20.0)
        .push_captioned(
            &locus(),
            "A",
            "Ideogram, coverage, reference, genes, variants, ruler",
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
            "Multiple sequence alignment, coloured by residue class",
        )
        .push_captioned(
            &variable_sites(),
            "G",
            "Variable sites with the phylogeny that orders them",
        )
        .push_captioned(&phylogeny(), "H", "A tree on its own")
        .push_captioned(
            &selection(),
            "I",
            "Windowed statistics read against a baseline they can fall below",
        )
        .push_captioned(
            &circular(),
            "J",
            "A circular chromosome: annotation, composition, and chords across the middle",
        )
        .push_captioned(
            &squiggle(),
            "K",
            "Raw nanopore current, with the bases a basecaller made of it",
        )
        .push_captioned(
            &cluster(),
            "L",
            "One locus in three genomes: the deletion that made BCG a vaccine",
        )
        .push_captioned(
            &methylation(),
            "M",
            "Methylation per site, one lane per strand, faded by read depth",
        )
        .push_captioned(
            &genome_wide(),
            "N",
            "A whole draft assembly on one axis, every contig of it",
        )
        .push_captioned(
            &structural(),
            "O",
            "Structural variants as arcs between their breakpoints",
        )
        .push_captioned(
            &frames(),
            "P",
            "Six reading frames: the stops, and what is open between them",
        )
        .push_captioned(
            &tanglegram(),
            "Q",
            "Two trees face to face, and where they disagree",
        )
        .push_captioned(
            &bisulfite(),
            "R",
            "Methylation one molecule at a time, not one fraction per site",
        );

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
                Feature::new(759_807, 763_325)
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
        MsaSequence::new("KatG_lineage2", b"MPEQHPPITETTTGATSNGCPV".to_vec()),
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

/// I: statistics that only mean anything relative to a line.
fn selection() -> Figure {
    let start = 2_150_000u64;
    let span = 40_000u64;
    let step = 500u64;
    let mut rng = Lcg::new(90_210);

    let windows: Vec<Window> = (0..(span / step))
        .map(|index| {
            let at = start + index * step;
            let diversifying = (28..36).contains(&index);
            let noise = (rng.next() % 1000) as f64 / 1000.0;
            let ratio = if diversifying {
                1.4 + 1.6 * noise
            } else {
                0.08 + 0.42 * noise
            };
            Window::new(at, at + step, ratio)
        })
        .collect();

    let bases: Vec<u8> = (0..span)
        .map(|index| {
            let leading = index < span / 3;
            match rng.next() % 100 {
                0..=32 if leading => b'G',
                0..=32 => b'C',
                33..=64 if leading => b'C',
                33..=64 => b'G',
                65..=82 => b'A',
                _ => b'T',
            }
        })
        .collect();

    Figure::new(Region::new("NC_000962.3", start, start + span).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            // No colour override: above the line and below it mean the same
            // thing here as in the GC skew track under it.
            WindowTrack::ratios(windows).label("pN/pS").height(58.0),
        )
        .push(
            WindowTrack::gc_skew(start, &bases, 1_000)
                .style(WindowStyle::Line)
                .label("GC skew")
                .height(48.0),
        )
        .push(AxisTrack::new())
}

/// J: the chromosome as the circle it actually is.
///
/// Not a `Figure`: a circle maps position to an angle rather than to an x, so
/// it is a different container. Both go on the sheet because both know how to
/// render themselves to an SVG of a stated size, which is all `Panels` asks.
fn circular() -> Rings {
    let length = 4_411_532u64;
    let step = 10_000u64;
    let mut rng = Lcg::new(1_906);

    let bases: Vec<u8> = (0..length)
        .map(|index| {
            let leading = index < length / 2;
            match rng.next() % 100 {
                0..=32 if leading => b'G',
                0..=32 => b'C',
                33..=64 if leading => b'C',
                33..=64 => b'G',
                65..=82 => b'A',
                _ => b'T',
            }
        })
        .collect();
    let skew = WindowTrack::gc_skew(0, &bases, step).windows().to_vec();

    let genes: Vec<Feature> = (0..700)
        .map(|index| {
            let start = index * (length / 700) + rng.next() % 900;
            Feature::new(start, start + 800 + rng.next() % 1_600).strand(if rng.next() % 2 == 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            })
        })
        .collect();

    Rings::new(length)
        .diameter(520.0)
        .title("H37Rv")
        .subtitle("4.41 Mb")
        .push(AxisRing::new())
        .push(FeatureRing::new(genes).thickness(18.0))
        .push(
            FeatureRing::new(vec![
                Feature::new(759_807, 763_325).name("rpoB"),
                Feature::new(2_153_889, 2_156_111).name("katG"),
            ])
            .thickness(2.0)
            .show_names(true)
            .split_strands(false),
        )
        .push(
            MarkerRing::categorised(vec![(761_155, 0), (2_155_168, 1), (1_673_425, 2)])
                .thickness(11.0)
                .width(2.0),
        )
        .push(
            SignalRing::new(skew)
                .thickness(46.0)
                .colors("#0072b2", "#d55e00"),
        )
        .link_colored((1_100_000, 1_180_000), (3_240_000, 3_320_000), None, 0.3)
}

/// K: current, before it was a sequence.
fn squiggle() -> Figure {
    let mut rng = Lcg::new(4_000);
    let mut signal: Vec<f64> = Vec::new();
    let mut moves: Vec<Move> = Vec::new();
    for base in b"GATCAGGCTAGCTTGAAACGT" {
        moves.push(Move::new(signal.len(), *base));
        let level = match base {
            b'A' => 86.0,
            b'C' => 98.0,
            b'G' => 110.0,
            _ => 122.0,
        };
        for _ in 0..(8 + rng.next() % 14) {
            signal.push(level + (rng.next() % 100) as f64 / 25.0 - 2.0);
        }
    }
    let samples = signal.len() as u64;
    Figure::new(Region::new("read", 0, samples).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            SquiggleTrack::new(0, signal)
                .moves(moves)
                .label("current")
                .height(76.0),
        )
        .push(AxisTrack::new().label("sample"))
}

/// L: the same locus in three genomes.
fn cluster() -> Figure {
    // Colour is spent only on what differs: the conserved backbone is one
    // quiet slate, the block RD1 removed is one hue and what replaced it is
    // another. Eight colours would make the reader learn a key before reading
    // anything.
    let theme = Theme::light();
    let kept = mix(&theme.muted, theme.surface(), 0.35);
    let lost = theme.color(0).to_string();
    let gained = theme.color(1).to_string();

    let gene = |start: u64, len: u64, name: &str, color: &str, forward: bool| {
        Feature::new(start, start.saturating_add(len))
            .name(name)
            .color(color)
            .strand(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
    };

    let core = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(300), 1_800, "eccA1", &kept, true),
            gene(at(2_400), 1_500, "eccB1", &kept, true),
            gene(at(4_100), 3_500, "eccCa1", &kept, true),
        ]
    };
    let rd1 = |shift: i64| {
        let at = |base: i64| (base + shift).max(0) as u64;
        vec![
            gene(at(7_800), 300, "PE35", &lost, true),
            gene(at(8_300), 1_100, "PPE68", &lost, true),
            gene(at(9_600), 300, "esxB", &lost, true),
            gene(at(10_100), 300, "esxA", &lost, true),
        ]
    };

    let mut h37rv = core(0);
    h37rv.extend(rd1(0));
    let mut bovis = core(-350);
    bovis.extend(rd1(-350));
    let mut bcg = core(150);
    bcg.push(gene(8_050, 1_350, "IS6110", &gained, false));

    let track = LocusTrack::new(vec![
        Locus::new("H37Rv", h37rv),
        Locus::new("M. bovis", bovis),
        Locus::new("BCG", bcg),
    ])
    .links(vec![
        Homology::new(0, 0, 0, 0.999),
        Homology::new(0, 1, 1, 0.997),
        Homology::new(0, 2, 2, 0.996),
        Homology::new(0, 3, 3, 0.95),
        Homology::new(0, 4, 4, 0.76),
        Homology::new(0, 5, 5, 0.999),
        Homology::new(0, 6, 6, 0.998),
        Homology::new(1, 0, 0, 0.998),
        Homology::new(1, 1, 1, 0.999),
        Homology::new(1, 2, 2, 0.997),
    ])
    .label("ESX-1");

    let (pale, dark) = track.ramp_ends(&theme);
    let legend = Legend::new()
        .area("deleted in BCG", &lost)
        .area("insertion sequence", &gained)
        .ramp("identity", pale, dark, "70%", "100%")
        .outline("in no neighbouring locus", theme.foreground.clone());

    Figure::new(Region::new("ESX-1", 0, 12_000).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(track)
        .push(LegendTrack::new(legend))
        .push(AxisTrack::new())
}

/// M: modification, not mutation.
fn methylation() -> Figure {
    let start = 1_460_000u64;
    let span = 3_000u64;
    let mut rng = Lcg::new(6_000);

    let mut sites = Vec::new();
    let mut pos = start + 40;
    while pos < start + span {
        let hemi = (start + 1_100..start + 1_700).contains(&pos);
        let coverage = (8 + rng.next() % 55) as u32;
        let forward = 0.88 + (rng.next() % 12) as f64 / 100.0;
        let reverse = if hemi {
            (rng.next() % 12) as f64 / 100.0
        } else {
            0.86 + (rng.next() % 14) as f64 / 100.0
        };
        sites.push(MethylSite::new(pos, Strand::Forward, forward, coverage));
        sites.push(MethylSite::new(pos, Strand::Reverse, reverse, coverage - 1));
        pos += 40 + rng.next() % 130;
    }

    Figure::new(Region::new("NC_000962.3", start, start + span).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(MethylationTrack::new(sites).label("6mA").height(66.0))
        .push(AxisTrack::new())
}

/// N: a scan across every contig at once.
fn genome_wide() -> Figure {
    let mut rng = Lcg::new(11_235);
    let genome = Genome::new(
        (0..12)
            .map(|index| {
                let length = 900_000 / (index + 1) + 40_000;
                (format!("contig_{:02}", index + 1), length)
            })
            .collect::<Vec<_>>(),
    );

    let mut hits: Vec<(String, u64, f64)> = Vec::new();
    for contig in genome.sequences() {
        let mut at = 200u64;
        while at < contig.length {
            let peak = contig.name == "contig_04" && at.abs_diff(120_000) < 9_000;
            let noise = (rng.next() % 1000) as f64 / 1000.0;
            let value = if peak {
                5.0 + 4.5 * noise
            } else {
                0.2 + 2.6 * noise
            };
            hits.push((contig.name.clone(), at, value));
            at += 900 + rng.next() % 2_600;
        }
    }
    let (mapped, _) = genome.map(hits);

    Figure::new(genome.region())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            ManhattanTrack::new(
                mapped
                    .iter()
                    .map(|(at, value)| Association::new(*at, *value))
                    .collect::<Vec<_>>(),
            )
            .bands(genome.boundaries())
            .genome_wide_threshold()
            .unit(" -log10 p")
            .label("association")
            .height(92.0),
        )
        .push(GenomeTrack::new(genome).label("contigs"))
        .push(AxisTrack::new())
}

/// O: arcs between breakpoints.
fn structural() -> Figure {
    let span = 60_000u64;
    let mut rng = Lcg::new(7_919);
    let calls = vec![
        StructuralVariant::new(8_000, 14_000, SvKind::Deletion)
            .support(34)
            .name("RD1"),
        StructuralVariant::new(24_000, 31_000, SvKind::Duplication).support(21),
        StructuralVariant::new(38_000, 44_000, SvKind::Inversion).support(9),
        StructuralVariant::new(5_000, 52_000, SvKind::Translocation).support(6),
        StructuralVariant::new(20_000, 20_000, SvKind::Insertion).support(17),
    ];
    let depth: Vec<f64> = (0..span)
        .map(|at| {
            let base = if (8_000..14_000).contains(&at) {
                2.0
            } else if (24_000..31_000).contains(&at) {
                84.0
            } else {
                42.0
            };
            base + (rng.next() % 70) as f64 / 10.0
        })
        .collect();

    let track = StructuralTrack::new(calls).label("SV").height(88.0);
    let theme = Theme::light();
    let legend = Legend::new()
        .line("deletion", track.color_of(SvKind::Deletion, &theme))
        .line("duplication", track.color_of(SvKind::Duplication, &theme))
        .line("inversion", track.color_of(SvKind::Inversion, &theme))
        .line(
            "translocation",
            track.color_of(SvKind::Translocation, &theme),
        )
        .line("insertion", track.color_of(SvKind::Insertion, &theme));

    Figure::new(Region::new("NC_000962.3", 0, span).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(track)
        .push(CoverageTrack::new(0, depth).label("depth").height(50.0))
        .push(LegendTrack::new(legend))
        .push(AxisTrack::new())
}

/// P: the six frames.
fn frames() -> Figure {
    let mut rng = Lcg::new(2_718);
    let span = 3_600usize;
    let mut seq: Vec<u8> = (0..span)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    let gene_at = 900usize;
    seq[gene_at..gene_at + 3].copy_from_slice(b"ATG");
    let mut at = gene_at + 3;
    while at + 3 <= gene_at + 1_500 {
        let codon = loop {
            let picked = [
                b"ACGT"[(rng.next() % 4) as usize],
                b"ACGT"[(rng.next() % 4) as usize],
                b"ACGT"[(rng.next() % 4) as usize],
            ];
            if !matches!(
                picked,
                [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
            ) {
                break picked;
            }
        };
        seq[at..at + 3].copy_from_slice(&codon);
        at += 3;
    }
    seq[at..at + 3].copy_from_slice(b"TAA");

    Figure::new(Region::new("plasmid", 0, span as u64).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(OrfTrack::new(0, seq).min_codons(60).label("frames"))
        .push(AxisTrack::new())
}

/// Q: two trees that disagree.
fn tanglegram() -> Figure {
    let core = Tree::parse_newick(
        "(((ERR5001:0.004,ERR5002:0.004):0.010,(ERR5003:0.003,ERR5004:0.005):0.009):0.020,\
          ((ERR5005:0.004,ERR5006:0.003):0.011,(ERR5007:0.005,ERR5008:0.004):0.008):0.018);",
    )
    .expect("the tree in this example is well formed");
    let accessory = Tree::parse_newick(
        "(((ERR5001:0.030,ERR5006:0.028):0.040,(ERR5003:0.031,ERR5004:0.029):0.038):0.050,\
          ((ERR5005:0.030,ERR5002:0.032):0.041,(ERR5007:0.028,ERR5008:0.030):0.039):0.048);",
    )
    .expect("the tree in this example is well formed");

    Figure::new(Region::new("taxa", 0, 8).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            TanglegramTrack::new(core, accessory)
                .names("core genome", "accessory genome")
                .row_height(18.0)
                .label("8 isolates"),
        )
}

/// R: one row per molecule.
fn bisulfite() -> Figure {
    let mut rng = Lcg::new(1_729);
    let start = 1_460_000u64;
    let sites: Vec<u64> = (0..14)
        .map(|index| {
            if index < 9 {
                start + 40 + index * 22
            } else {
                start + 420 + (index - 9) * 70
            }
        })
        .collect();
    let molecules: Vec<Molecule> = (0..14)
        .map(|index| {
            let methylated = index % 2 == 0;
            let calls: Vec<Option<bool>> = (0..sites.len())
                .map(|site| {
                    if site > 10 && rng.next() % 100 < 30 {
                        return None;
                    }
                    Some(if methylated {
                        rng.next() % 100 > 8
                    } else {
                        rng.next() % 100 < 6
                    })
                })
                .collect();
            Molecule::new(format!("read_{:02}", index + 1), calls)
        })
        .collect();

    Figure::new(Region::new("NC_000962.3", start, start + 800).unwrap())
        .width(WIDTH)
        .show_region_label(false)
        .push(
            BisulfiteTrack::new(sites, molecules)
                .label("CpG")
                .row_height(12.0),
        )
        .push(AxisTrack::new())
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
