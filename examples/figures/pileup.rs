//! The figure `cargo run --example pileup` writes, as a function.
//!
//! The example writes it to a file, and the documentation site's playground
//! includes this file too and draws it in the page, in the page's colours, at
//! the width of the column and over whatever window the reader has moved to.
//! One copy of the figure serves both, so the two cannot drift apart.

use karyon::{CigarOp, Drawing, Plot, Read, ReadColoring, Region, Strand, Theme, Variant};

/// Window start, 0-based.
const START: u64 = 4_000;
/// Window length in bases.
const WINDOW: usize = 420;
/// Where the variant sits, 0-based.
const VARIANT: u64 = 4_210;
/// Where the deletion starts, 0-based.
const DELETION: u64 = 4_120;

/// `example-pileup.svg`: a variant worth looking at twice.
///
/// `theme`, `width` and `region` replace the light theme, the 920 pixels and
/// the 420 base window the committed figure is drawn with.
pub fn example_pileup(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
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

    let own = Region::parse("NC_000962.3:4001-4420").expect("the locus in this example parses");
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("A variant worth looking at twice")
        .theme(theme.clone())
        .width(width.unwrap_or(920.0))
        // Both from the first base of the window the reads were made for,
        // rather than from the left edge of the one on display, which is only
        // the same place until the window moves.
        .add_coverage_at(START, depth)
        .label("depth")
        .adjust(|track| track.height(50.0))
        .add_variants(vec![Variant::new(VARIANT)
            .value(fraction)
            .category("candidate SNV")])
        .label("call")
        .adjust(|track| track.height(34.0))
        .add_sequence_at(START, reference.clone())
        .label("reference")
        .add_pileup(reads)
        .adjust(|track| {
            track
                .reference(START, reference)
                .coloring(ReadColoring::Strand)
                .fade_by_quality(true)
                .max_rows(Some(22))
        })
        .label("reads")
        .into_figure();
    Box::new(figure)
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
