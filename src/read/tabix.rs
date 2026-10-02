//! The rows of a bgzipped text file over a window, found through the tabix
//! index beside it, as `tabix -h` prints them.
//!
//! A VCF, a BED, a bedGraph or a GFF3 compressed with bgzip and indexed with
//! `tabix` (a `.tbi`) or `tabix -C` and `bcftools index` (a `.csi`) can be read
//! a window at a time: [`index`](super::index) says which stretches of the file
//! hold the rows over a window, and [`Bgzf`] reads those stretches as lines.
//! A figure of one gene out of a whole genome's calls reads the few blocks
//! over that gene, where the whole file is otherwise inflated and parsed: on
//! a VCF of 200 samples and 825 MB of text, a window of two thousand bases
//! took 4.4 s and 869 MB read whole, and takes 5 ms and 4 MB through its
//! index.
//!
//! Nothing here opens a file: each reader takes anything that reads and
//! seeks, and an index already read by [`index::parse`](super::index::parse).
//!
//! ```
//! use std::io::Cursor;
//! use karyon::read::{index, point, tabix};
//! use karyon::Region;
//!
//! // The calls the documentation draws, and the index tabix wrote for them.
//! let vcf = include_bytes!("../../docs/data/calls.vcf.gz");
//! let tbi = include_bytes!("../../docs/data/calls.vcf.gz.tbi");
//! let index = index::parse(tbi)?;
//! let region = Region::parse("NC_000962.3:761,000-763,000")?;
//! let text = tabix::window(Cursor::new(&vcf[..]), &index, &region)?;
//! assert!(text.starts_with("##fileformat=VCF"));
//! let calls = point::variants(&text, &region)?;
//! assert!(!calls.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # What comes out
//!
//! Text the readers already take, as they take the whole file. [`head`] is
//! the lines before the first row, which is where a VCF names its samples on
//! its `#CHROM` line, a GFF3 says it is one and a table names its columns,
//! and the first row itself, which says what kind of file it is when a window
//! holds no row of it. [`rows`] is the rows the index holds over a window, in
//! the order of the file, and [`window`] is the two together.
//!
//! The header is what tabix itself takes for one: the lines it was told to
//! skip (`-S`), and every line at the top starting with the comment
//! character, `#` unless it was told otherwise (`-c`).
//!
//! A window's rows include some that lie near it and not over it: htslib
//! folds the bins of a small stretch of file into their parents, so a bin over
//! the window can hold rows a few kilobases to its left. Every reader drops
//! them, as it drops the rows outside the window of a whole file. A sequence
//! the index does not name holds no rows.
//!
//! # What is refused
//!
//! An index that does not describe the file beside it, which is what an index
//! written for another file, or for an earlier version of this one, is: the
//! rows not starting where the index says they start, a first row on another
//! sequence than the index names first, and a row the index puts on one
//! sequence that names another, that has no start where the index says, or
//! that starts before the row read before it, which tabix refuses to index.
//! An index of a BAM or a BCF, which says nothing of a row's columns, is
//! refused too, and so is any block [`Bgzf`] refuses. None of it panics.
//!
//! Each row is checked only as far as its sequence and its start, which is
//! all the index says about it; what the rest of it holds is its reader's to
//! refuse.

use std::io::{Read, Seek};

use crate::Region;

use super::bgzf::Bgzf;
use super::index::{Columns, Index};
use super::ReadError;

/// What a bgzipped text file holds before its rows, and its first row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Head {
    /// The header lines, each with its newline after it.
    pub text: String,
    /// The first row, without its newline, or `None` for a file of no rows.
    pub first_row: Option<String>,
}

/// The columns an index says a row has, or the refusal of an index that says
/// nothing of them.
fn columns_of(index: &Index) -> Result<Columns, ReadError> {
    index.columns().ok_or_else(|| {
        ReadError::whole(
            "the index is one of a BAM or a BCF, which says nothing of the columns a row \
             of text has",
        )
    })
}

/// The comment character as the byte a line starts with.
fn comment_byte(columns: &Columns) -> u8 {
    u8::try_from(u32::from(columns.comment)).unwrap_or(b'#')
}

/// Where a virtual offset is, in words: the compressed offset of a block and,
/// where it is not the block's first byte, how far into it.
fn at(offset: u64) -> String {
    match (offset >> 16, offset & 0xffff) {
        (block, 0) => format!("at byte {block}"),
        (block, within) => format!("{within} bytes into the block at byte {block}"),
    }
}

/// A line read out of a file as text.
fn text_of(line: &[u8]) -> Result<&str, ReadError> {
    std::str::from_utf8(line).map_err(|_| ReadError::whole("a line of the file is not UTF-8 text"))
}

/// The sequence a row names and where it starts, 0-based, by the columns the
/// index keeps. `None` where the row has no such columns or no number in the
/// start's, which a row tabix indexed always has.
fn placed<'a>(line: &'a str, columns: &Columns) -> Option<(&'a str, u64)> {
    let line = line.trim_end_matches('\r');
    let mut sequence = None;
    let mut start = None;
    for (at, field) in line.split('\t').enumerate() {
        if at + 1 == columns.sequence {
            sequence = Some(field);
        }
        if at + 1 == columns.start {
            start = Some(field.trim().parse::<u64>().ok()?);
        }
        if sequence.is_some() && start.is_some() {
            break;
        }
    }
    let start = start?;
    // Counted from one unless the index says from nought, as tabix counts
    // it: `-p bed` and `-0` are 0-based, VCF, SAM, GFF3 and plain columns
    // 1-based.
    let start = if columns.zero_based {
        start
    } else {
        start.saturating_sub(1)
    };
    Some((sequence?, start))
}

/// The header of a bgzipped text file and its first row, checked against the
/// index: the first row has to start where the index says the rows of its
/// first sequence start, and be on that sequence.
///
/// # Errors
///
/// An index of a BAM or a BCF, which has no columns; a file whose rows do not
/// start where the index says, or whose first row is on another sequence,
/// which is a file the index was not written for; and any block that is not
/// BGZF or is damaged.
pub fn head<R: Read + Seek>(reader: R, index: &Index) -> Result<Head, ReadError> {
    let columns = columns_of(index)?;
    let comment = comment_byte(&columns);
    let mut bgzf = Bgzf::new(reader);
    let mut text = String::new();
    let mut line = Vec::new();
    let mut number = 0usize;
    let mut first_row = None;
    loop {
        let offset = bgzf.tell();
        if !bgzf.line(&mut line)? {
            break;
        }
        number += 1;
        // As tabix reads a header: the lines it was told to skip, whatever
        // they hold, and then the lines that start with the comment
        // character, up to the first that does not.
        if number <= columns.skip || line.first() == Some(&comment) {
            text.push_str(text_of(&line)?);
            text.push('\n');
            continue;
        }
        let row = text_of(&line)?;
        if let Some(summary) = index.summary(0) {
            if summary.first != offset {
                return Err(ReadError::whole(format!(
                    "the index puts the first row {}, and the file's header ends {}",
                    at(summary.first),
                    at(offset)
                )));
            }
        }
        if let (Some(named), Some((sequence, _))) = (index.names().first(), placed(row, &columns)) {
            if named != sequence {
                return Err(ReadError::whole(format!(
                    "the index names {named} first, and the file's first row is on {sequence}"
                )));
            }
        }
        first_row = Some(row.to_string());
        break;
    }
    // A file of no rows has an index of no sequences.
    if first_row.is_none() && index.summary(0).is_some_and(|summary| summary.placed > 0) {
        return Err(ReadError::whole(
            "the index counts rows, and the file has none after its header",
        ));
    }
    Ok(Head { text, first_row })
}

/// The rows the index holds over `region`, each with its newline after it,
/// in the order of the file, without the header. Some lie near the window
/// rather than over it, since a bin over the window holds them, and a reader
/// drops them as it drops the rows outside the window of a whole file.
///
/// The read stops at the first row that starts at or past the window's end,
/// since tabix indexes only a file sorted by where each row starts and no
/// row after that one can reach back into the window. A sequence the index
/// does not name holds no rows.
///
/// # Errors
///
/// An index of a BAM or a BCF; a row where the index puts one of the window's
/// sequence that is on another, has no start in the column the index says,
/// or starts before the row before it, all of which say the index was written
/// for another file; and any block that is not BGZF or is damaged.
pub fn rows<R: Read + Seek>(
    reader: R,
    index: &Index,
    region: &Region,
) -> Result<String, ReadError> {
    let columns = columns_of(index)?;
    let comment = comment_byte(&columns);
    let Some(reference) = index.reference(region.seq()) else {
        return Ok(String::new());
    };
    let mut bgzf = Bgzf::new(reader);
    let mut out = String::new();
    let mut line = Vec::new();
    let mut last: Option<u64> = None;
    'chunks: for (begin, end) in index.chunks(reference, region.start(), region.end()) {
        bgzf.seek(begin)?;
        while bgzf.tell() < end {
            let offset = bgzf.tell();
            if !bgzf.line(&mut line)? {
                return Err(ReadError::whole(format!(
                    "the index has rows of {} past the end of the file",
                    region.seq()
                )));
            }
            // A comment between rows, which tabix passes over as it indexes.
            if line.first() == Some(&comment) {
                continue;
            }
            let row = text_of(&line)?;
            let Some((sequence, start)) = placed(row, &columns) else {
                return Err(ReadError::whole(format!(
                    "the index puts a row of {} {}, and it has no start in column {}",
                    region.seq(),
                    at(offset),
                    columns.start
                )));
            };
            if sequence != region.seq() {
                return Err(ReadError::whole(format!(
                    "the index puts a row of {} where the file has one of {sequence}",
                    region.seq()
                )));
            }
            if last.is_some_and(|last| start < last) {
                return Err(ReadError::whole(format!(
                    "the index puts rows of {} in order where the file has one starting at \
                     {} after one at {}",
                    region.seq(),
                    start + 1,
                    last.map_or(0, |last| last + 1)
                )));
            }
            last = Some(start);
            if start >= region.end() {
                break 'chunks;
            }
            out.push_str(row);
            out.push('\n');
        }
    }
    Ok(out)
}

/// The header and then the rows over `region`: what `tabix -h` prints, and
/// what a reader takes for the window as it takes the whole file for it.
///
/// # Errors
///
/// What [`head`] and [`rows`] refuse.
pub fn window<R: Read + Seek>(
    mut reader: R,
    index: &Index,
    region: &Region,
) -> Result<String, ReadError> {
    let head = head(&mut reader, index)?;
    let rows = rows(&mut reader, index, region)?;
    Ok(head.text + &rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::index::fixture::*;
    use crate::read::index::parse;
    use crate::read::{gzip, interval, point};
    use std::io::Cursor;

    fn parsed(bytes: &[u8]) -> Index {
        parse(bytes).unwrap_or_else(|error| panic!("{error}"))
    }

    /// The windows of a grid over both files' sequences: the edges of a
    /// leaf, a megabase, nothing at all, and a sequence the index does not
    /// name.
    fn grid() -> Vec<Region> {
        let mut regions = Vec::new();
        for sequence in ["chr1", "chr2", "chr3", "chr4"] {
            for start in [
                0u64, 8, 9, 10, 5_000, 16_379, 16_383, 16_384, 16_399, 49_999, 50_000, 99_999,
                100_000, 600_050, 1_048_575, 1_048_576, 1_499_999, 2_999_999, 5_000_000,
            ] {
                for width in [1u64, 2, 100, 16_384, 100_000, 3_000_000] {
                    regions.push(Region::new(sequence, start, start + width).unwrap());
                }
            }
        }
        regions
    }

    /// Every index of text the foundation holds, with the file and the text
    /// inside it.
    fn every() -> [(&'static str, Index, &'static [u8], &'static str); 4] {
        [
            ("calls, tbi", parsed(&CALLS_TBI), &CALLS[..], CALLS_TEXT),
            ("calls, csi", parsed(&CALLS_CSI), &CALLS[..], CALLS_TEXT),
            ("rows, tbi", parsed(&ROWS_TBI), &ROWS[..], ROWS_TEXT),
            ("rows, csi", parsed(&ROWS_CSI), &ROWS[..], ROWS_TEXT),
        ]
    }

    /// For every window of the grid, a VCF's calls and a BED's features read
    /// from the window are the ones read from the whole file, through a TBI
    /// and through a CSI. A row lost to a floor set too high or a stop made
    /// too early is a call or a gene missing here.
    #[test]
    fn a_window_draws_what_the_whole_file_draws() {
        let mut windows = 0;
        for (what, index, data, text) in every() {
            for region in grid() {
                let window = window(Cursor::new(data), &index, &region)
                    .unwrap_or_else(|error| panic!("{what}, {region}: {error}"));
                if what.starts_with("calls") {
                    assert_eq!(
                        point::variants(&window, &region).unwrap(),
                        point::variants(text, &region).unwrap(),
                        "{what}, {region}"
                    );
                    assert_eq!(
                        point::genotypes(&window, &region, None).unwrap().sites,
                        point::genotypes(text, &region, None).unwrap().sites,
                        "{what}, {region}"
                    );
                } else {
                    assert_eq!(
                        interval::features(&window, &region, None).unwrap(),
                        interval::features(text, &region, None).unwrap(),
                        "{what}, {region}"
                    );
                }
                windows += 1;
            }
        }
        assert_eq!(windows, 4 * 4 * 19 * 6);
    }

    /// The `#CHROM` line that names the samples comes with every window, even
    /// one that holds no row and one on a sequence the index does not name,
    /// since a reader of genotypes has nothing to name its rows by without it.
    #[test]
    fn the_header_comes_with_every_window_even_an_empty_one() {
        let header = &CALLS_TEXT[..CALLS_TEXT.find("chr1\t10").unwrap()];
        for index in [parsed(&CALLS_TBI), parsed(&CALLS_CSI)] {
            for locus in ["chr1:10-10", "chr1:11-12", "chr2:6-79999", "chr9:1-100"] {
                let region = Region::parse(locus).unwrap();
                let window = window(Cursor::new(&CALLS[..]), &index, &region).unwrap();
                assert!(window.starts_with(header), "{locus}: {window}");
                assert_eq!(point::samples(&window), ["A", "B"], "{locus}");
            }
            let head = head(Cursor::new(&CALLS[..]), &index).unwrap();
            assert_eq!(head.text, header);
            assert_eq!(
                head.first_row.as_deref(),
                Some("chr1\t10\tv1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1")
            );
        }
    }

    /// The rows of a window, as the lines of the file they are, in its order.
    /// Row v3 is cut across two blocks, and is one row.
    #[test]
    fn a_row_split_across_two_blocks_is_one_row() {
        let index = parsed(&CALLS_TBI);
        let region = Region::parse("chr1:16390-16390").unwrap();
        let rows = rows(Cursor::new(&CALLS[..]), &index, &region).unwrap();
        let v3 = CALLS_TEXT
            .lines()
            .find(|line| line.contains("\tv3\t"))
            .unwrap();
        assert!(rows.lines().any(|row| row == v3), "{rows}");
        assert!(
            rows.lines()
                .all(|row| CALLS_TEXT.lines().any(|line| line == row)),
            "{rows}"
        );
    }

    /// An index read against the file it was not written for is refused, by
    /// where it says the rows begin, by the sequence it names first, or by a
    /// row of the window that is not where it says.
    #[test]
    fn an_index_for_another_file_is_refused() {
        let error = head(Cursor::new(&ROWS[..]), &parsed(&CALLS_TBI)).unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("the index puts the first row"),
            "{error}"
        );
        // The same rows behind a header of another length: the rows begin
        // in another place.
        let longer = gzip::decompress(&CALLS).unwrap();
        let longer = format!("##another=line\n{}", String::from_utf8(longer).unwrap());
        let longer = crate::read::bgzf::fixture::blocks(longer.as_bytes(), &[]);
        let error = head(Cursor::new(&longer), &parsed(&CALLS_TBI)).unwrap_err();
        assert!(
            error.to_string().contains("the file's header ends"),
            "{error}"
        );
        // A header and no rows, beside an index that counts six on chr1.
        let header = &CALLS_TEXT[..CALLS_TEXT.find("chr1\t10").unwrap()];
        let bare = crate::read::bgzf::fixture::blocks(header.as_bytes(), &[]);
        let error = head(Cursor::new(&bare), &parsed(&CALLS_TBI)).unwrap_err();
        assert!(
            error.to_string().contains("none after its header"),
            "{error}"
        );
        // An index of a BAM says nothing of a row's columns.
        let bam = parsed(&crate::read::bam::fixture::BAI);
        let error = head(Cursor::new(&CALLS[..]), &bam).unwrap_err();
        assert!(error.to_string().contains("BAM"), "{error}");
        // Rows read where the index points that are not the window's: the
        // VCF's index read over the BED, whose first block it skips.
        let region = Region::parse("chr2:1-90000").unwrap();
        assert!(rows(Cursor::new(&ROWS[..]), &parsed(&CALLS_TBI), &region).is_err());
    }

    /// A TBI of a BED with `names` and one stretch for the first, from its
    /// root bin, which every window on it is in.
    fn root(names: &[&str], stretch: (u64, u64)) -> Index {
        parsed(&crate::read::index::fixture::rooted(
            names,
            [0x10000, 1, 2, 3],
            stretch,
        ))
    }

    /// The rows an index says are of one sequence, in order, that are not:
    /// a row of another sequence among them, one starting before the row
    /// before it, one with no start where the index says, and rows the index
    /// says go on past the end of the file. Each says the index does not
    /// describe the file, and none is handed over as a row.
    #[test]
    fn rows_that_are_not_where_the_index_says_are_refused() {
        let window = Region::parse("chr1:1-1,000").unwrap();
        let all = |text: &str| {
            let data = crate::read::bgzf::fixture::blocks(text.as_bytes(), &[]);
            // The last block, the empty one these fixtures end with, is 31
            // bytes: where the rows end.
            let end = ((data.len() as u64) - 31) << 16;
            rows(Cursor::new(&data), &root(&["chr1"], (0, end)), &window)
        };
        assert_eq!(
            all("chr1\t5\t6\nchr1\t9\t20\n").unwrap(),
            "chr1\t5\t6\nchr1\t9\t20\n"
        );
        // A comment between rows is passed over, as tabix passes over it.
        assert_eq!(
            all("chr1\t5\t6\n#a note\nchr1\t9\t20\n").unwrap(),
            "chr1\t5\t6\nchr1\t9\t20\n"
        );
        for (text, said) in [
            ("chr1\t5\t6\nchr2\t7\t8\n", "one of chr2"),
            (
                "chr1\t50\t60\nchr1\t10\t20\n",
                "starting at 11 after one at 51",
            ),
            ("chr1\tfive\t6\n", "no start in column 2"),
        ] {
            let error = all(text).unwrap_err().to_string();
            assert!(error.contains(said), "{text:?}: {error}");
        }
        let data = crate::read::bgzf::fixture::blocks(b"chr1\t5\t6\n", &[]);
        let error = rows(Cursor::new(&data), &root(&["chr1"], (0, 1 << 40)), &window)
            .unwrap_err()
            .to_string();
        assert!(error.contains("past the end of the file"), "{error}");
    }

    /// A first row on another sequence than the one the index names first is
    /// a file the index was not written for, though its rows begin where the
    /// index says: the same index with its two names swapped.
    #[test]
    fn the_first_row_is_on_the_sequence_the_index_names_first() {
        let mut swapped = gzip::decompress(&CALLS_TBI).unwrap();
        let at = swapped
            .windows(10)
            .position(|bytes| bytes == b"chr1\0chr2\0")
            .unwrap();
        swapped[at..at + 10].copy_from_slice(b"chr2\0chr1\0");
        let error = head(Cursor::new(&CALLS[..]), &parsed(&swapped)).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the index names chr2 first, and the file's first row is on chr1"
        );
    }

    /// The header ends where the index says the rows begin, in the files it
    /// was written for: each fixture's header fills a block of its own, so
    /// the first row starts the next block, which a reader that says the end
    /// of a block as that block at its length would put elsewhere.
    #[test]
    fn the_header_ends_where_the_index_says_the_rows_begin() {
        for (what, index, data, text) in every() {
            let head = head(Cursor::new(data), &index).unwrap_or_else(|e| panic!("{what}: {e}"));
            let first = text
                .lines()
                .find(|line| !line.starts_with('#'))
                .map(str::to_string);
            assert_eq!(head.first_row, first, "{what}");
        }
    }

    /// A reader that counts how many bytes are asked of it.
    struct Counting<'a> {
        inner: Cursor<&'a [u8]>,
        furthest: u64,
    }

    impl Read for Counting<'_> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.furthest = self.furthest.max(self.inner.position());
            Ok(n)
        }
    }

    impl Seek for Counting<'_> {
        fn seek(&mut self, to: std::io::SeekFrom) -> std::io::Result<u64> {
            self.inner.seek(to)
        }
    }

    /// A window stops at the first row past its end, and never reads the
    /// blocks after it. The stretch the index gives for the first base of
    /// chr1 runs on to v5, a megabase further, since htslib folds the bins
    /// of a small file into their parents: read to its end, it is three of
    /// the file's four blocks of rows. Stopped at v2, the first row past the
    /// window, it is the one block that holds v1 and v2.
    #[test]
    fn a_window_stops_at_the_first_row_past_its_end() {
        let index = parsed(&CALLS_TBI);
        let region = Region::parse("chr1:10-10").unwrap();
        let stretch = index.chunks(0, region.start(), region.end());
        // The block that holds v4, after the one that holds v1 and v2.
        let v4 = index.chunks(0, 49_999, 50_000)[0].0 >> 16;
        assert!(stretch[0].1 >> 16 >= v4, "the stretch reaches past v4");
        let mut counting = Counting {
            inner: Cursor::new(&CALLS[..]),
            furthest: 0,
        };
        let rows = rows(&mut counting, &index, &region).unwrap();
        let v1 = CALLS_TEXT
            .lines()
            .find(|line| line.contains("\tv1\t"))
            .unwrap();
        assert_eq!(rows, format!("{v1}\n"));
        assert!(
            counting.furthest <= v4,
            "read to {} of {v4}",
            counting.furthest
        );
    }

    /// A file or an index cut short or with a byte changed anywhere is read
    /// or refused, never a panic.
    #[test]
    fn a_damaged_index_or_block_is_an_error_and_never_a_panic() {
        let region = Region::parse("chr1:1-2,000,000").unwrap();
        for (_, index, data, _) in every() {
            for cut in 0..data.len() {
                let _ = window(Cursor::new(&data[..cut]), &index, &region);
            }
            for at in 0..data.len() {
                let mut bent = data.to_vec();
                bent[at] ^= 0x55;
                let _ = window(Cursor::new(&bent), &index, &region);
            }
        }
        for bytes in [&CALLS_TBI[..], &CALLS_CSI[..]] {
            let inflated = gzip::decompress(bytes).unwrap();
            for at in 0..inflated.len() {
                let mut bent = inflated.clone();
                bent[at] ^= 0x80;
                if let Ok(index) = parse(&bent) {
                    let _ = window(Cursor::new(&CALLS[..]), &index, &region);
                }
            }
        }
    }
}
