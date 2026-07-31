//! Renders the read pileup figure used in the README.
//!
//! ```text
//! cargo run --example pileup -- assets
//! ```
//!
//! A synthetic locus built to look like the thing you open a genome browser
//! for: a variant carried by about half the reads, a deletion carried by a
//! few, an insertion in one, and a patch of low mapping quality. The coverage
//! track is computed from the reads themselves, so the two tracks cannot
//! disagree.

use std::env;
use std::path::PathBuf;

use karyon::{
    AxisTrack, CigarOp, CoverageTrack, Figure, PileupTrack, Read, ReadColoring, Region,
    SequenceTrack, Strand, Variant, VariantTrack,
};

/// Window start, 0-based.
const START: u64 = 4_000;
/// Window length in bases.
const WINDOW: usize = 420;
/// Where the variant sits, 0-based.
const VARIANT: u64 = 4_210;
/// Where the deletion starts, 0-based.
const DELETION: u64 = 4_120;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let reference = synthetic_reference(WINDOW);
    let reads = synthetic_reads(&reference);

    // Depth straight off the reads, so the profile and the pileup agree by
    // construction rather than by hope.
    let mut depth = vec![0.0; WINDOW];
    for read in &reads {
        for segment in read.segments() {
            if let karyon::Segment::Aligned { start, len, .. } = segment {
                for position in start..start + len {
                    if let Some(slot) = depth.get_mut((position - START) as usize) {
                        *slot += 1.0;
                    }
                }
            }
        }
    }

    let carriers = reads
        .iter()
        .filter(|r| r.base_at(VARIANT) == Some(b'T'))
        .count();
    let covering = reads
        .iter()
        .filter(|r| r.base_at(VARIANT).is_some())
        .count();
    let fraction = carriers as f64 / covering.max(1) as f64;

    let region = Region::new("NC_000962.3", START, START + WINDOW as u64).unwrap();
    let figure = Figure::new(region)
        .title("A variant worth looking at twice")
        .width(920.0)
        .push(CoverageTrack::new(START, depth).label("depth").height(50.0))
        .push(
            VariantTrack::new(vec![Variant::new(VARIANT)
                .value(fraction)
                .category("candidate SNV")])
            .label("call")
            .height(34.0),
        )
        .push(SequenceTrack::new(START, reference.clone()).label("reference"))
        .push(
            PileupTrack::new(reads)
                .reference(START, reference)
                .coloring(ReadColoring::Strand)
                .fade_by_quality(true)
                .max_rows(Some(22))
                .label("reads"),
        )
        .push(AxisTrack::new());

    figure.save_svg(out.join("example-pileup.svg"))?;
    let (width, height) = figure.dimensions();
    println!("example-pileup.svg {width:.0} x {height:.0}");
    println!(
        "{carriers} of {covering} reads carry the variant, an allele fraction of {fraction:.2}"
    );
    Ok(())
}

fn synthetic_reference(len: usize) -> Vec<u8> {
    let mut rng = Lcg::new(20_260_731);
    (0..len)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect()
}

/// Reads tiled across the window, some of them carrying something.
fn synthetic_reads(reference: &[u8]) -> Vec<Read> {
    let mut rng = Lcg::new(31);
    let mut reads = Vec::new();

    for index in 0..46u64 {
        let start = START + (index * 9) % 300;
        let length = 110 + (rng.next() % 20) as u32;
        let end = (start + length as u64).min(START + WINDOW as u64);
        if end <= start {
            continue;
        }
        let length = (end - start) as u32;

        let carries_variant = rng.next() % 100 < 45;
        let carries_deletion = !carries_variant && rng.next() % 100 < 12;
        let carries_insertion = index == 17;
        let strand = if index % 2 == 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        // A patch of the locus maps poorly, the way a repeat does.
        let quality = if (4_060..4_090).contains(&start) {
            (rng.next() % 8) as u8
        } else {
            60
        };

        let spans = |position: u64| position >= start && position < end;

        let (cigar, sequence) = if carries_deletion && spans(DELETION) && spans(DELETION + 6) {
            let before = (DELETION - start) as u32;
            let after = length - before - 6;
            let mut bases = slice(reference, start, before as usize);
            bases.extend(slice(reference, DELETION + 6, after as usize));
            (
                vec![
                    CigarOp::Match(before),
                    CigarOp::Deletion(6),
                    CigarOp::Match(after),
                ],
                bases,
            )
        } else if carries_insertion && spans(4_180) {
            let before = (4_180 - start) as u32;
            let after = length - before;
            let mut bases = slice(reference, start, before as usize);
            bases.extend_from_slice(b"TTTT");
            bases.extend(slice(reference, 4_180, after as usize));
            (
                vec![
                    CigarOp::Match(before),
                    CigarOp::Insertion(4),
                    CigarOp::Match(after),
                ],
                bases,
            )
        } else {
            let mut bases = slice(reference, start, length as usize);
            if carries_variant && spans(VARIANT) {
                bases[(VARIANT - start) as usize] = b'T';
            }
            // A sprinkle of sequencing error, so the real signal has to be
            // told apart from the noise the way it is in life.
            if rng.next() % 100 < 30 {
                let at = (rng.next() as usize) % bases.len();
                bases[at] = b"ACGT"[(rng.next() % 4) as usize];
            }
            (vec![CigarOp::Match(length)], bases)
        };

        reads.push(
            Read::new(start, cigar)
                .sequence(sequence)
                .strand(strand)
                .mapping_quality(quality),
        );
    }

    reads
}

fn slice(reference: &[u8], start: u64, len: usize) -> Vec<u8> {
    let offset = (start - START) as usize;
    reference
        .iter()
        .skip(offset)
        .take(len)
        .copied()
        .collect::<Vec<u8>>()
}

/// A linear congruential generator, so the figure is reproducible without a
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
