//! bigBed, intervals along a genome, read a window at a time through the
//! index it holds.
//!
//! A bigBed is the rows of a BED file packed in blocks that an R-tree points
//! into: each row as its sequence's number, its start and its end, and the
//! rest of it as the text it was written as. [`bed`] reads the blocks over a
//! window and writes their rows back out as BED, which is what
//! [`interval::features`](super::interval::features) reads, so a bigBed is
//! drawn by the same reader as the BED it was made from.
//!
//! # Which columns are kept
//!
//! The first `definedFieldCount` columns, as the header says, and none after.
//! Those are the ones that mean what BED says they mean; the rest are the
//! file's own, named by the autoSql it carries. A peak caller's narrowPeak is
//! BED6 and four of its own, so its seventh column is a signal value, and kept
//! it would be read as a BED12's `thickStart`, a coding start at 5.2. A BED12
//! keeps all twelve, so a gene keeps its exons and the stretch that codes.
//!
//! # Coordinates
//!
//! 0-based and half-open, as BED is, and passed through.
//!
//! # What is refused
//!
//! A file that is not a bigBed, one damaged or cut short, and a sequence the
//! file does not have, named with the ones it has.
//!
//! ```
//! use std::io::Cursor;
//! use karyon::{read, Region};
//!
//! # let bytes = include_bytes!("fixtures/peaks.bb").to_vec();
//! let region = Region::parse("chr1:1-250")?;
//! let bed = read::bigbed::bed(Cursor::new(bytes), Some(&region))?;
//! assert_eq!(bed, "chr1\t100\t200\tpeak1\t500\t.\n");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::fmt::Write as _;
use std::io::{Read, Seek};

use super::bbi::{Bbi, Over, BIGBED};
use super::ReadError;
use crate::Region;

/// The sequences a bigBed names, each with its length, in the order of its
/// index, which is the order of their names.
///
/// # Errors
///
/// A file that is not a bigBed, or is damaged.
pub fn sequences<R: Read + Seek>(reader: R) -> Result<Vec<(String, u64)>, ReadError> {
    let mut file = Bbi::open(reader, BIGBED, "bigBed")?;
    Ok(file
        .sequences()?
        .into_iter()
        .map(|named| (named.name, u64::from(named.length)))
        .collect())
}

/// The rows over `region` as BED text, or every row of the file for `None`,
/// each cut to the columns the header says are BED's own.
///
/// A row is over the window when any base of it is, so a gene that starts
/// before the window and runs into it is kept whole.
///
/// # Errors
///
/// A file that is not a bigBed or is damaged, and a region on a sequence the
/// file does not have, naming the ones it has.
pub fn bed<R: Read + Seek>(reader: R, region: Option<&Region>) -> Result<String, ReadError> {
    let mut file = Bbi::open(reader, BIGBED, "bigBed")?;
    // The extra columns a row has, past its sequence, start and end.
    let extra = usize::from(file.defined).saturating_sub(3);
    let (over, names) = match region {
        Some(region) => {
            let named = file.named(region.seq())?;
            let over = Over {
                id: named.id,
                start: u32::try_from(region.start()).unwrap_or(u32::MAX),
                end: u32::try_from(region.end()).unwrap_or(u32::MAX),
            };
            (Some(over), vec![(named.id, named.name)])
        }
        None => (
            None,
            file.sequences()?
                .into_iter()
                .map(|named| (named.id, named.name))
                .collect(),
        ),
    };
    let full = file.full;
    let mut out = String::new();
    for (offset, size) in file.blocks(full, over)? {
        let block = file.block(offset, size)?;
        rows(&file, &block, over, &names, extra, &mut out)?;
    }
    Ok(out)
}

/// The rows of one block over `over`, or all of them for `None`, as BED on
/// `out`, each named by the sequence `names` gives its number and cut to
/// `extra` columns after its end.
///
/// A row on a sequence `names` does not give is left out: for a window,
/// `names` is the window's sequence alone, which leaves out the rows of a
/// block that runs on to the next sequence, and for the whole file, a row on
/// a sequence the index does not name is nowhere a figure can put it.
fn rows<R: Read + Seek>(
    file: &Bbi<R>,
    block: &[u8],
    over: Option<Over>,
    names: &[(u32, String)],
    extra: usize,
    out: &mut String,
) -> Result<(), ReadError> {
    let mut bytes = file.bytes(block);
    while bytes.left() > 0 {
        let id = bytes.u32()?;
        let start = bytes.u32()?;
        let end = bytes.u32()?;
        let rest = bytes.until_nought()?;
        let Some((_, name)) = names.iter().find(|(known, _)| *known == id) else {
            continue;
        };
        if over.is_some_and(|over| start >= over.end || end <= over.start) {
            continue;
        }
        // Writing to a string does not fail.
        let _ = write!(out, "{name}\t{start}\t{end}");
        let rest = String::from_utf8_lossy(rest);
        for column in rest.split('\t').filter(|_| !rest.is_empty()).take(extra) {
            out.push('\t');
            out.push_str(column);
        }
        out.push('\n');
    }
    Ok(())
}

/// A bigBed written by hand, for rows kent's own tools would not write.
#[cfg(test)]
pub(crate) mod fixture {
    use crate::read::bbi::fixture::Numbers;

    /// A bigBed of one sequence holding `rows`, each its start, its end and
    /// the rest of it as text, in one block, with `defined` columns of BED.
    pub(crate) fn written(
        name: &str,
        length: u32,
        rows: &[(u32, u32, &str)],
        defined: u16,
    ) -> Vec<u8> {
        let n = Numbers(false);
        let mut block = Vec::new();
        for (start, end, rest) in rows {
            for number in [0, *start, *end] {
                n.u32(&mut block, number);
            }
            block.extend(rest.as_bytes());
            block.push(0);
        }
        let first = rows.iter().map(|row| row.0).min().unwrap_or(0);
        let last = rows.iter().map(|row| row.1).max().unwrap_or(0);
        n.file(
            super::BIGBED,
            defined,
            name,
            length,
            &[(first, last, block)],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const GENES: &[u8] = include_bytes!("fixtures/genes.bb");
    const PEAKS: &[u8] = include_bytes!("fixtures/peaks.bb");

    fn window(bytes: &[u8], locus: &str) -> String {
        bed(Cursor::new(bytes), Some(&Region::parse(locus).unwrap())).unwrap()
    }

    /// Every row, and every window's rows, as `bigBedToBed` prints them: a
    /// BED12 keeps its twelve columns, exons and coding stretch with them.
    #[test]
    fn a_window_holds_the_rows_bigbedtobed_prints() {
        let tool = include_str!("fixtures/genes.bb.bed");
        assert_eq!(bed(Cursor::new(GENES), None).unwrap(), tool);
        let scores = include_bytes!("fixtures/scores.bb");
        assert_eq!(
            bed(Cursor::new(scores), None).unwrap(),
            include_str!("fixtures/scores.bb.bed")
        );
        let on = |sequence: &str| -> String {
            tool.lines()
                .filter(|line| line.starts_with(&format!("{sequence}\t")))
                .map(|line| format!("{line}\n"))
                .collect()
        };
        assert_eq!(window(GENES, "chr1:1-1000"), on("chr1"));
        assert_eq!(window(GENES, "chr2:1-500"), on("chr2"));
        // Read as BED, a gene keeps its three exons.
        let region = Region::parse("chr1:1-1000").unwrap();
        let features =
            crate::read::interval::features(&window(GENES, "chr1:1-1000"), &region, None).unwrap();
        let gene = features
            .iter()
            .find(|feature| feature.name.as_deref() == Some("geneA"))
            .unwrap();
        assert_eq!(gene.exons, [(200, 300), (500, 700), (750, 900)]);
        assert_eq!(gene.coding, [(250, 850)]);
    }

    /// A row that starts before the window and runs into it is kept whole,
    /// and a window between rows holds none.
    #[test]
    fn an_item_starting_before_the_window_is_kept() {
        assert_eq!(
            window(GENES, "chr1:801-860"),
            "chr1\t200\t900\tgeneA\t0\t+\t250\t850\t0\t3\t100,200,150,\t0,300,550,\n"
        );
        assert_eq!(window(GENES, "chr1:901-950"), "");
        assert_eq!(
            window(GENES, "chr1:100-100"),
            "chr1\t99\t100\tone\t0\t+\t99\t100\t0\t1\t1,\t0,\n"
        );
    }

    /// A peak caller's BED6 and four columns of its own keeps six: its
    /// signal value is not read as where a gene's coding starts.
    #[test]
    fn a_bed6_plus_4_keeps_six_columns() {
        let tool = include_str!("fixtures/peaks.bb.bed");
        assert_eq!(tool.lines().next().unwrap().split('\t').count(), 10);
        let rows = bed(Cursor::new(PEAKS), None).unwrap();
        let cut: String = tool
            .lines()
            .map(|line| line.split('\t').take(6).collect::<Vec<_>>().join("\t") + "\n")
            .collect();
        assert_eq!(rows, cut);
        let region = Region::parse("chr1:1-1000").unwrap();
        let features = crate::read::interval::features(&rows, &region, None).unwrap();
        assert_eq!(features.len(), 2);
        assert!(features
            .iter()
            .all(|feature| feature.exons.is_empty() && feature.coding.is_empty()));
    }

    #[test]
    fn a_sequence_it_has_not_got_is_refused_with_the_ones_it_has() {
        assert_eq!(
            sequences(Cursor::new(GENES)).unwrap(),
            [("chr1".to_string(), 1000), ("chr2".to_string(), 500)]
        );
        let region = Region::parse("chr9:1-100").unwrap();
        let error = bed(Cursor::new(GENES), Some(&region)).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the bigBed has no sequence called chr9; it has chr1, chr2"
        );
        let error = bed(Cursor::new(include_bytes!("fixtures/signal.bw")), None).unwrap_err();
        assert!(error.to_string().contains("not a bigBed"), "{error}");
    }

    /// A block may run from one sequence on to the next, and a window keeps
    /// the rows of its own sequence alone, however their positions fall.
    #[test]
    fn a_block_over_two_sequences_gives_a_window_its_own_rows() {
        let file = Bbi::open(Cursor::new(GENES), BIGBED, "bigBed").unwrap();
        let mut block = Vec::new();
        for (id, name) in [(0u32, "a"), (1, "b"), (0, "c")] {
            for number in [id, 5, 10] {
                block.extend(number.to_le_bytes());
            }
            block.extend(name.as_bytes());
            block.extend(b"\t0\t+\0");
        }
        let over = Over {
            id: 1,
            start: 0,
            end: 100,
        };
        let mut out = String::new();
        rows(
            &file,
            &block,
            Some(over),
            &[(1, "chr2".to_string())],
            3,
            &mut out,
        )
        .unwrap();
        assert_eq!(out, "chr2\t5\t10\tb\t0\t+\n");
        // Every row of the block, for the file read whole.
        let names = [(0, "chr1".to_string()), (1, "chr2".to_string())];
        let mut out = String::new();
        rows(&file, &block, None, &names, 1, &mut out).unwrap();
        assert_eq!(out, "chr1\t5\t10\ta\nchr2\t5\t10\tb\nchr1\t5\t10\tc\n");
    }

    #[test]
    fn a_damaged_file_is_an_error_and_never_a_panic() {
        let region = Region::parse("chr1:1-1000").unwrap();
        for bytes in [GENES, PEAKS] {
            for cut in 0..bytes.len() {
                let _ = bed(Cursor::new(&bytes[..cut]), Some(&region));
                let _ = bed(Cursor::new(&bytes[..cut]), None);
            }
            for at in 0..bytes.len() {
                for change in [0xff, 0x55] {
                    let mut bent = bytes.to_vec();
                    bent[at] ^= change;
                    let _ = bed(Cursor::new(&bent), Some(&region));
                    let _ = bed(Cursor::new(&bent), None);
                    let _ = sequences(Cursor::new(&bent));
                }
            }
        }
    }
}
