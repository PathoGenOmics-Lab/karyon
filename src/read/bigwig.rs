//! bigWig, a signal along a genome, read a window at a time through the index
//! it holds.
//!
//! A bigWig is the values of a bedGraph or a wiggle file, packed in blocks
//! that an R-tree points into, and the same values summed up again at a few
//! coarser scales, the zoom levels, each a list of bins holding the least, the
//! most and the sum of the values in them. [`window`] reads the blocks over a
//! window and turns them into spans of `(start, end, value)`, 0-based and
//! half-open as the file holds them, which is what a
//! [`CoverageTrack`](crate::CoverageTrack) is painted from.
//!
//! # A window many bases to a pixel
//!
//! A chromosome drawn 900 pixels wide is a quarter of a million bases to a
//! pixel, and reading its values as written would read the whole of it to
//! draw 900 columns. So where a pixel holds two bins or more of a zoom level,
//! the coarsest such level is read instead, and each bin is painted over its
//! bases with the value the track's [`Aggregate`] takes of a pixel: its most,
//! its least, or its sum spread over every base of it. A base no value
//! covers counts as nought, as it does when a track is painted from spans,
//! so a bin with gaps in it paints at least nought for its most, at most
//! nought for its least, and its mean with the gaps counted in.
//!
//! Drawn that way a column is the one the values as written draw, but for a
//! bin that straddles two columns, which lends its most to both: a peak is
//! drawn at most a column wider than it is, as in UCSC's own browser, and no
//! column lower. Over a chromosome of 248,956,422 bases written as 4,684,581
//! spans, 638 of 792 columns came out the same as from the values as written
//! and 154 higher, read in 2.6 ms and 2.6 MB where the values as written took
//! 1.34 s and 298 MB.
//!
//! The levels are read from the file rather than assumed, since kent picks
//! them from the data: 3,184 bases for the finest of one file these tests
//! read, 1,904 for another, and 119 for that one's values written stored
//! rather than compressed, each level four times the one under it. A bin
//! starts at the first base it covers, not on a multiple of its size, and
//! runs a level's reduction from there.
//!
//! A bigWig written without zoom levels, as some writers allow, is read as
//! written over any window, which is all of it for a whole chromosome.
//!
//! # What is read
//!
//! Each of the three kinds of section a bigWig holds its values in: bedGraph
//! sections, a start, an end and a value per item; variable-step ones, a
//! start and a value, each item the section's span long; and fixed-step ones,
//! a value alone, each item a step after the one before. All three are what
//! `bigWigToBedGraph` prints, and a test holds every one to it. The values
//! are single precision in the file and widened here, exactly.
//!
//! # What is refused
//!
//! A file that is not a bigWig, one damaged or cut short, a section of a kind
//! the format does not have, a span that ends before it starts, and a
//! sequence the file does not have, named with the ones it has.
//!
//! ```
//! use std::io::Cursor;
//! use karyon::{read, Aggregate, Region};
//!
//! # let bytes = include_bytes!("fixtures/signal.bw").to_vec();
//! let region = Region::parse("chr1:1-200")?;
//! let signal = read::bigwig::window(Cursor::new(bytes), &region, 1.0, Aggregate::Max)?;
//! assert_eq!(signal.spans[0], (10, 20, 1.5));
//! assert_eq!(signal.zoom, None);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::io::{Read, Seek};

use super::bbi::{Bbi, Over, BIGWIG};
use super::ReadError;
use crate::{Aggregate, Region};

/// What a window of a bigWig holds.
#[derive(Debug, Clone, PartialEq)]
pub struct Signal {
    /// Each value over the bases it covers, as `(start, end, value)`,
    /// 0-based and half-open, in the order the file holds them.
    pub spans: Vec<(u64, u64, f64)>,
    /// The bin size of the zoom level the spans were read from, or `None`
    /// where they are the values as written.
    pub zoom: Option<u32>,
}

/// The sequences a bigWig names, each with its length, in the order of its
/// index, which is the order of their names.
///
/// # Errors
///
/// A file that is not a bigWig, or is damaged.
pub fn sequences<R: Read + Seek>(reader: R) -> Result<Vec<(String, u64)>, ReadError> {
    let mut file = Bbi::open(reader, BIGWIG, "bigWig")?;
    Ok(file
        .sequences()?
        .into_iter()
        .map(|named| (named.name, u64::from(named.length)))
        .collect())
}

/// The values over `region`, read from the zoom level a pixel of
/// `bases_per_pixel` bases holds two bins of where there is one, and as
/// written otherwise.
///
/// `aggregate` is what a pixel is drawn with, which picks what a zoomed bin
/// is painted with; the values as written are the same whichever it is. A
/// `bases_per_pixel` of one or less always reads them as written.
///
/// # Errors
///
/// A file that is not a bigWig or is damaged, and a region on a sequence the
/// file does not have, naming the ones it has.
pub fn window<R: Read + Seek>(
    reader: R,
    region: &Region,
    bases_per_pixel: f64,
    aggregate: Aggregate,
) -> Result<Signal, ReadError> {
    let mut file = Bbi::open(reader, BIGWIG, "bigWig")?;
    let named = file.named(region.seq())?;
    let over = Over {
        id: named.id,
        start: clamp(region.start()),
        end: clamp(region.end()),
    };
    let zoom = file
        .zooms
        .iter()
        .rev()
        .find(|zoom| zoom.reduction > 0 && f64::from(zoom.reduction) * 2.0 <= bases_per_pixel)
        .copied();
    let mut spans = Vec::new();
    match zoom {
        Some(zoom) => {
            for (offset, size) in file.blocks(zoom.index, Some(over))? {
                let block = file.block(offset, size)?;
                summaries(&file, &block, over, aggregate, &mut spans)?;
            }
        }
        None => {
            let full = file.full;
            for (offset, size) in file.blocks(full, Some(over))? {
                let block = file.block(offset, size)?;
                section(&file, &block, over, &mut spans)?;
            }
        }
    }
    Ok(Signal {
        spans,
        zoom: zoom.map(|zoom| zoom.reduction),
    })
}

/// Spans as bedGraph, on `sequence`: what the command line hands the readers
/// of text, so a bigWig is drawn by the same reader as the bedGraph it was
/// written from. Each value is written as the shortest decimal that reads
/// back as the same number.
pub fn bedgraph(sequence: &str, spans: &[(u64, u64, f64)]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(spans.len() * (sequence.len() + 24));
    for (start, end, value) in spans {
        // Writing to a string does not fail.
        let _ = writeln!(out, "{sequence}\t{start}\t{end}\t{value}");
    }
    out
}

/// A position as the file holds one, which is never past `u32::MAX`.
fn clamp(position: u64) -> u32 {
    u32::try_from(position).unwrap_or(u32::MAX)
}

/// The items of one block of the values as written that overlap `over`.
fn section<R: Read + Seek>(
    file: &Bbi<R>,
    block: &[u8],
    over: Over,
    spans: &mut Vec<(u64, u64, f64)>,
) -> Result<(), ReadError> {
    let mut bytes = file.bytes(block);
    // A block may hold several sections, one after another.
    while bytes.left() > 0 {
        let id = bytes.u32()?;
        let first = u64::from(bytes.u32()?);
        let _last = bytes.u32()?;
        let step = u64::from(bytes.u32()?);
        let span = u64::from(bytes.u32()?);
        let kind = bytes.u8()?;
        let _reserved = bytes.u8()?;
        let count = bytes.u16()?;
        for item in 0..u64::from(count) {
            let (start, end, value) = match kind {
                1 => {
                    let start = u64::from(bytes.u32()?);
                    let end = u64::from(bytes.u32()?);
                    (start, end, bytes.f32()?)
                }
                2 => {
                    let start = u64::from(bytes.u32()?);
                    (start, start + span, bytes.f32()?)
                }
                3 => {
                    let start = first + item * step;
                    (start, start + span, bytes.f32()?)
                }
                _ => {
                    return Err(ReadError::whole(format!(
                        "the bigWig is damaged: a section is of kind {kind}, and the format has \
                         kinds 1 to 3"
                    )))
                }
            };
            if end < start {
                return Err(ReadError::whole(format!(
                    "the bigWig is damaged: a span ends at {end}, before it starts at {start}"
                )));
            }
            if id == over.id && start < u64::from(over.end) && end > u64::from(over.start) {
                spans.push((start, end, f64::from(value)));
            }
        }
    }
    Ok(())
}

/// The bins of one block of a zoom level that overlap `over`, each painted
/// with what `aggregate` takes of a pixel.
fn summaries<R: Read + Seek>(
    file: &Bbi<R>,
    block: &[u8],
    over: Over,
    aggregate: Aggregate,
    spans: &mut Vec<(u64, u64, f64)>,
) -> Result<(), ReadError> {
    let mut bytes = file.bytes(block);
    while bytes.left() > 0 {
        let id = bytes.u32()?;
        let start = u64::from(bytes.u32()?);
        let end = u64::from(bytes.u32()?);
        let covered = u64::from(bytes.u32()?);
        let least = f64::from(bytes.f32()?);
        let most = f64::from(bytes.f32()?);
        let sum = f64::from(bytes.f32()?);
        let _squares = bytes.f32()?;
        if id != over.id
            || end <= start
            || start >= u64::from(over.end)
            || end <= u64::from(over.start)
        {
            continue;
        }
        // A base of the bin no value covers is nought, as it is in a track
        // painted from spans.
        let gaps = covered < end - start;
        let value = match aggregate {
            Aggregate::Max if gaps => most.max(0.0),
            Aggregate::Max => most,
            Aggregate::Min if gaps => least.min(0.0),
            Aggregate::Min => least,
            Aggregate::Mean => sum / (end - start) as f64,
        };
        if value.is_finite() {
            spans.push((start, end, value));
        }
    }
    Ok(())
}

/// A bigWig written by hand, for what kent's own tools do not write: the
/// variable-step and fixed-step sections only `wigToBigWig` writes, a file of
/// no zoom levels, and a file in the other byte order.
#[cfg(test)]
pub(crate) mod fixture {
    use crate::read::bbi::fixture::Numbers;

    /// The items of one section, as each kind holds them.
    pub(crate) enum Items {
        /// Each item's start, end and value.
        BedGraph(Vec<(u32, u32, f32)>),
        /// Each item's start and value, all `span` long.
        Variable { span: u32, items: Vec<(u32, f32)> },
        /// A value per item, the first at `start` and each `step` after.
        Fixed {
            start: u32,
            step: u32,
            span: u32,
            values: Vec<f32>,
        },
    }

    /// A bigWig of one sequence holding `sections`, each a block of its own
    /// and stored rather than compressed, with no zoom levels, written
    /// big-endian where `big` says so.
    pub(crate) fn written(big: bool, name: &str, length: u32, sections: &[Items]) -> Vec<u8> {
        let n = Numbers(big);
        let mut blocks = Vec::new();
        for items in sections {
            let mut block = Vec::new();
            let (first, last, step, span, kind, count) = match items {
                Items::BedGraph(items) => {
                    for (start, end, value) in items {
                        n.u32(&mut block, *start);
                        n.u32(&mut block, *end);
                        n.f32(&mut block, *value);
                    }
                    let first = items.iter().map(|item| item.0).min().unwrap_or(0);
                    let last = items.iter().map(|item| item.1).max().unwrap_or(0);
                    (first, last, 0, 0, 1u8, items.len())
                }
                Items::Variable { span, items } => {
                    for (start, value) in items {
                        n.u32(&mut block, *start);
                        n.f32(&mut block, *value);
                    }
                    let first = items.first().map_or(0, |item| item.0);
                    let last = items.last().map_or(0, |item| item.0 + span);
                    (first, last, 0, *span, 2, items.len())
                }
                Items::Fixed {
                    start,
                    step,
                    span,
                    values,
                } => {
                    for value in values {
                        n.f32(&mut block, *value);
                    }
                    let last = start + (values.len() as u32 - 1) * step + span;
                    (*start, last, *step, *span, 3, values.len())
                }
            };
            let mut section = Vec::new();
            for number in [0, first, last, step, span] {
                n.u32(&mut section, number);
            }
            section.extend([kind, 0]);
            n.u16(&mut section, count as u16);
            section.extend(block);
            blocks.push((first, last, section));
        }
        n.file(super::BIGWIG, 0, name, length, &blocks)
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::{written, Items};
    use super::*;
    use crate::read::bbi::Bbi;
    use std::io::Cursor;

    const SIGNAL: &[u8] = include_bytes!("fixtures/signal.bw");
    const STORED: &[u8] = include_bytes!("fixtures/signal.unc.bw");
    const STEPS: &[u8] = include_bytes!("fixtures/steps.bw");

    /// The rows `bigWigToBedGraph` printed for a file, on one sequence, with
    /// each value read as the single precision number the file holds.
    fn printed(text: &str, sequence: &str) -> Vec<(u64, u64, f64)> {
        text.lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|cols| cols[0] == sequence)
            .map(|cols| {
                let value: f32 = cols[3].parse().unwrap();
                (
                    cols[1].parse().unwrap(),
                    cols[2].parse().unwrap(),
                    f64::from(value),
                )
            })
            .collect()
    }

    fn read(bytes: &[u8], locus: &str, per_pixel: f64, aggregate: Aggregate) -> Signal {
        let region = Region::parse(locus).unwrap();
        window(Cursor::new(bytes), &region, per_pixel, aggregate).unwrap()
    }

    /// Compressed and stored, every sequence reads as the tool that undoes
    /// the file prints it, and the bedGraph written from it reads back as the
    /// same spans.
    #[test]
    fn a_window_holds_what_bigwigtobedgraph_prints() {
        let tool = include_str!("fixtures/signal.bw.bedgraph");
        for bytes in [SIGNAL, STORED] {
            for (locus, sequence) in [("chr1:1-1000", "chr1"), ("chr2", "chr2"), ("chr3", "chr3")] {
                let locus = if locus.contains(':') {
                    locus.to_string()
                } else {
                    format!("{locus}:1-1000")
                };
                let signal = read(bytes, &locus, 1.0, Aggregate::Max);
                assert_eq!(signal.zoom, None);
                assert_eq!(signal.spans, printed(tool, sequence), "{locus}");
                let text = bedgraph(sequence, &signal.spans);
                let region = Region::parse(&locus).unwrap();
                let again =
                    crate::read::signal::spans(&text, &region, Some(crate::Format::BedGraph))
                        .unwrap();
                assert_eq!(again, signal.spans, "{locus}");
            }
        }
        // A window takes the spans with a base in it, whole.
        let signal = read(SIGNAL, "chr1:16-100", 1.0, Aggregate::Max);
        assert_eq!(
            signal.spans,
            [(10, 20, 1.5), (20, 30, 3.0), (99, 100, 4.25)]
        );
        assert!(read(SIGNAL, "chr1:31-99", 1.0, Aggregate::Max)
            .spans
            .is_empty());
        // -0.1 is not a single precision number, and the one the file holds
        // is written out whole, so it reads back as itself.
        let signal = read(SIGNAL, "chr3:1-60", 1.0, Aggregate::Max);
        assert_eq!(
            bedgraph("chr3", &signal.spans),
            "chr3\t0\t1\t-1.25\nchr3\t40\t50\t-0.10000000149011612\n"
        );
    }

    #[test]
    fn three_hundred_spans_read_as_written_where_a_pixel_holds_few_bases() {
        let tool = include_str!("fixtures/steps.bw.bedgraph");
        let every = printed(tool, "chr1");
        assert_eq!(every.len(), 300);
        assert_eq!(
            read(STEPS, "chr1:1-125,000", 100.0, Aggregate::Max).spans,
            every
        );
        let inside: Vec<(u64, u64, f64)> = every
            .iter()
            .copied()
            .filter(|(start, end, _)| *start < 60_000 && *end > 40_000)
            .collect();
        assert_eq!(
            read(STEPS, "chr1:40,001-60,000", 3.0, Aggregate::Mean).spans,
            inside
        );
    }

    /// The coarsest level of which a pixel holds two bins, read off the
    /// levels the file has, which kent chose from the data.
    #[test]
    fn a_window_many_bases_to_a_pixel_reads_the_zoom_level() {
        let at = |per_pixel: f64| read(STEPS, "chr1:1-125,000", per_pixel, Aggregate::Max).zoom;
        assert_eq!(at(6_367.0), None);
        assert_eq!(at(6_368.0), Some(3_184));
        assert_eq!(at(25_471.0), Some(3_184));
        assert_eq!(at(25_472.0), Some(12_736));
        assert_eq!(at(1e9), Some(203_776));
        assert_eq!(at(f64::NAN), None);
        // A file with one level, coarser than its whole sequence.
        assert_eq!(
            read(SIGNAL, "chr1:1-1000", 1e6, Aggregate::Max).zoom,
            Some(1_904)
        );
    }

    /// Each bin paints what the values as written add up to over its bases,
    /// a base no value covers counting as nought: its most, its least, its
    /// mean. Held to the values themselves at every level, under all three.
    #[test]
    fn a_zoomed_bin_paints_what_the_values_under_it_add_up_to() {
        // Every level of the stepped file, and the one level of the small
        // file stored and compressed, whose chr3 holds only values under
        // nought, so a bin with gaps there has nought for its most.
        let mut cases = Vec::new();
        for per_pixel in [6_368.0, 25_472.0, 101_888.0, 407_552.0] {
            cases.push((STEPS, "chr1:1-125,000", per_pixel));
        }
        for bytes in [SIGNAL, STORED] {
            for locus in ["chr1:1-1,000", "chr2:1-500", "chr3:1-60"] {
                cases.push((bytes, locus, 1e7));
            }
        }
        for (bytes, locus, per_pixel) in cases {
            let every = read(bytes, locus, 1.0, Aggregate::Max).spans;
            let value_at = |base: u64| {
                every
                    .iter()
                    .find(|(start, end, _)| *start <= base && base < *end)
                    .map_or(0.0, |span| span.2)
            };
            for aggregate in [Aggregate::Max, Aggregate::Min, Aggregate::Mean] {
                let signal = read(bytes, locus, per_pixel, aggregate);
                assert!(signal.zoom.is_some());
                assert!(!signal.spans.is_empty());
                for (start, end, painted) in &signal.spans {
                    let values: Vec<f64> = (*start..*end).map(value_at).collect();
                    let expected = match aggregate {
                        Aggregate::Max => values.iter().copied().fold(f64::MIN, f64::max),
                        Aggregate::Min => values.iter().copied().fold(f64::MAX, f64::min),
                        Aggregate::Mean => values.iter().sum::<f64>() / values.len() as f64,
                    };
                    // The sum is kept in single precision.
                    assert!(
                        (painted - expected).abs() <= 1e-6 * expected.abs().max(1.0),
                        "{locus} {per_pixel} {aggregate:?} {start}..{end}: {painted} against \
                         {expected}"
                    );
                }
            }
        }
        // The bin over all of chr3 has values under nought and gaps between
        // them, and paints nought for its most.
        let chr3 = read(STORED, "chr3:1-60", 1e7, Aggregate::Max);
        assert_eq!(chr3.spans, [(0, 60, 0.0)]);
    }

    /// The sections only `wigToBigWig` writes, in a file written by hand,
    /// which `bigWigToBedGraph` 482 printed as these rows.
    #[test]
    fn fixed_and_variable_step_sections_read_as_bigwigtobedgraph_prints_them() {
        let sections = [
            Items::BedGraph(vec![(10, 20, 1.5), (20, 30, -2.0)]),
            Items::Variable {
                span: 5,
                items: vec![(100, 4.0), (120, 0.25)],
            },
            Items::Fixed {
                start: 200,
                step: 10,
                span: 3,
                values: vec![7.0, 8.5, 9.0],
            },
        ];
        let tool = "chrZ\t10\t20\t1.5\nchrZ\t20\t30\t-2\nchrZ\t100\t105\t4\nchrZ\t120\t125\t0.25\n\
                    chrZ\t200\t203\t7\nchrZ\t210\t213\t8.5\nchrZ\t220\t223\t9\n";
        for big in [false, true] {
            let bytes = written(big, "chrZ", 400, &sections);
            // No zoom levels at all, so a window of any scale reads the
            // values as written.
            for per_pixel in [1.0, 1e6] {
                let signal = read(&bytes, "chrZ:1-400", per_pixel, Aggregate::Max);
                assert_eq!(signal.zoom, None);
                assert_eq!(signal.spans, printed(tool, "chrZ"), "big-endian {big}");
            }
            assert_eq!(
                read(&bytes, "chrZ:205-212", 1.0, Aggregate::Max).spans,
                [(210, 213, 8.5)]
            );
            assert_eq!(
                sequences(Cursor::new(&bytes)).unwrap(),
                [("chrZ".to_string(), 400)]
            );
        }
    }

    /// A block may hold sections of more than one sequence, as a zoom
    /// level's blocks do, and a window keeps its own sequence's alone.
    #[test]
    fn a_block_over_two_sequences_gives_a_window_its_own_values() {
        let file = Bbi::open(Cursor::new(SIGNAL), BIGWIG, "bigWig").unwrap();
        let mut block = Vec::new();
        for (id, value) in [(1u32, 9.0f32), (0, 2.0)] {
            for number in [id, 10, 20, 0, 0] {
                block.extend(number.to_le_bytes());
            }
            block.extend([1, 0]);
            block.extend(1u16.to_le_bytes());
            for number in [10u32, 20] {
                block.extend(number.to_le_bytes());
            }
            block.extend(value.to_le_bytes());
        }
        let over = Over {
            id: 0,
            start: 0,
            end: 100,
        };
        let mut spans = Vec::new();
        section(&file, &block, over, &mut spans).unwrap();
        assert_eq!(spans, [(10, 20, 2.0)]);
        // Two bins on two sequences, likewise.
        let mut block = Vec::new();
        for (id, most) in [(1u32, 9.0f32), (0, 2.0)] {
            for number in [id, 10, 20, 10] {
                block.extend(number.to_le_bytes());
            }
            for number in [most, most, most * 10.0, 0.0] {
                block.extend(number.to_le_bytes());
            }
        }
        let mut spans = Vec::new();
        summaries(&file, &block, over, Aggregate::Max, &mut spans).unwrap();
        assert_eq!(spans, [(10, 20, 2.0)]);
    }

    #[test]
    fn the_sequences_are_named_with_their_lengths() {
        let held = sequences(Cursor::new(SIGNAL)).unwrap();
        assert_eq!(
            held,
            [
                ("chr1".to_string(), 1000),
                ("chr2".to_string(), 500),
                ("chr3".to_string(), 60)
            ]
        );
    }

    #[test]
    fn a_sequence_it_has_not_got_is_refused_with_the_ones_it_has() {
        let region = Region::parse("chrX:1-100").unwrap();
        let error = window(Cursor::new(SIGNAL), &region, 1.0, Aggregate::Max).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the bigWig has no sequence called chrX; it has chr1, chr2, chr3"
        );
        // A name longer than any key is no key.
        let region = Region::parse("chromosome_1:1-100").unwrap();
        assert!(window(Cursor::new(SIGNAL), &region, 1.0, Aggregate::Max).is_err());
        let error = window(
            Cursor::new(include_bytes!("fixtures/genes.bb")),
            &region,
            1.0,
            Aggregate::Max,
        )
        .unwrap_err();
        assert!(error.to_string().contains("not a bigWig"), "{error}");
    }

    /// Cut at every length, and every byte changed in turn, a bigWig is
    /// either read or refused, and nothing panics or asks for more memory
    /// than the file is long. Release builds check for overflow, so an offset
    /// added unchecked would panic here rather than wrap.
    #[test]
    fn a_damaged_file_is_an_error_and_never_a_panic() {
        let region = Region::parse("chr1:1-1000").unwrap();
        for bytes in [SIGNAL, STORED] {
            for cut in 0..bytes.len() {
                let _ = window(Cursor::new(&bytes[..cut]), &region, 1.0, Aggregate::Max);
                let _ = window(Cursor::new(&bytes[..cut]), &region, 1e6, Aggregate::Mean);
                let _ = sequences(Cursor::new(&bytes[..cut]));
            }
            for at in 0..bytes.len() {
                for change in [0xff, 0x55] {
                    let mut bent = bytes.to_vec();
                    bent[at] ^= change;
                    let _ = window(Cursor::new(&bent), &region, 1.0, Aggregate::Max);
                    let _ = window(Cursor::new(&bent), &region, 1e6, Aggregate::Min);
                    let _ = sequences(Cursor::new(&bent));
                }
            }
        }
        // The size every name is padded to, as large as it goes, is refused
        // before a name is padded to it: it asked for four gigabytes.
        let mut bent = SIGNAL.to_vec();
        let names = u64::from_le_bytes(bent[8..16].try_into().unwrap()) as usize;
        bent[names + 8..names + 12].copy_from_slice(&u32::MAX.to_le_bytes());
        let error = window(Cursor::new(&bent), &region, 1.0, Aggregate::Max).unwrap_err();
        assert!(
            error.to_string().contains("names are longer than the file"),
            "{error}"
        );
        // An index whose first child is the node it hangs from is walked
        // once, not for ever.
        let mut bent = STORED.to_vec();
        let index = u64::from_le_bytes(bent[24..32].try_into().unwrap());
        let root = index as usize + 48;
        assert_eq!(
            bent[root], 0,
            "the root of this file's index is an inner node"
        );
        bent[root + 4 + 16..root + 4 + 24].copy_from_slice(&(index + 48).to_le_bytes());
        let error = window(Cursor::new(&bent), &region, 1.0, Aggregate::Max).unwrap_err();
        assert!(
            error.to_string().contains("leads back into a node"),
            "{error}"
        );
        // And one whose child starts inside it, which is a node of no child
        // and no leaf, and walked unchecked leads into its neighbours.
        bent[root + 4 + 16..root + 4 + 24].copy_from_slice(&(index + 48 + 2).to_le_bytes());
        let error = window(Cursor::new(&bent), &region, 1.0, Aggregate::Max).unwrap_err();
        assert!(
            error.to_string().contains("leads back into a node"),
            "{error}"
        );
        // Each offset a window reads by, set as far as it goes: where the
        // names are, the blocks as written, and the first zoom level's.
        for (at, per_pixel) in [(8usize, 1.0), (24, 1.0), (80, 1e6)] {
            let mut bent = SIGNAL.to_vec();
            bent[at..at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
            let error = window(Cursor::new(&bent), &region, per_pixel, Aggregate::Max).unwrap_err();
            assert!(
                error.to_string().contains("past its own end"),
                "{at}: {error}"
            );
        }
    }
}
