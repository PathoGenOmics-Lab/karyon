//! BGZF, the blocked gzip a BAM, a bgzipped text file and their indexes are
//! written in, read from any place an index points to.
//!
//! bgzip writes a file as a series of gzip members holding at most 64 KiB
//! each, and every member says its own compressed size in its header, so a
//! reader can start at any of them without inflating the ones before it. An
//! index names a place in such a file by a virtual offset: the compressed
//! offset of a block's first byte in the high 48 bits, and a position in that
//! block's inflated bytes in the low 16. [`Bgzf`] goes to one and reads on from
//! there, block after block, as the bytes of binary records ([`Bgzf::fill`])
//! or as lines of text ([`Bgzf::line`]).
//!
//! Nothing here opens a file: [`Bgzf`] takes anything that reads and seeks,
//! which is a file for the command line and a buffer for a test or a page.
//! [`bam`](super::bam) reads its records through it, and
//! [`tabix`](super::tabix) the lines of a bgzipped text file that a tabix
//! index points into; which blocks to go to is [`index`](super::index)'s
//! answer.
//!
//! # Where a block ends
//!
//! A record or a line that ends exactly where its block ends is followed by
//! the next block's first byte, and htslib writes that place as the next block
//! at offset nought, never as this block at its length. bgzip 1.24 ends a
//! block where a VCF's header ends, so the index beside `docs/data/calls.vcf.gz`
//! puts the first row at 210 << 16, the start of the second block.
//! [`Bgzf::tell`] answers the same way. Answering the other way, a reader that
//! compares where it is with where an index says a stretch ends reads one
//! record past every stretch that ends on a block boundary, and a check that
//! the header ends where the index says the rows begin fails on a sound file.
//!
//! # What is refused
//!
//! A block that does not start the way bgzip writes one, that does not say its
//! size, or that is cut short; a block whose data does not inflate or does not
//! match its CRC, as [`gzip::member`] refuses one; and a virtual offset past
//! the end of its block. None of it panics. Ending early is not refused here,
//! since only the format inside knows whether a file may end where it does:
//! [`Bgzf::fill`] says how much it filled, and [`bam`](super::bam) says what a
//! BAM cut short is.

use std::io::{Read, Seek, SeekFrom};

use super::{gzip, ReadError};

/// Whether `data` starts the way bgzip writes a block: a gzip member whose
/// extra field holds the subfield called BC that says the block's size.
///
/// That subfield is what makes a file one an index can point into. Plain gzip
/// has none, and a gzip file of one member has to be inflated from its start
/// to reach any byte of it, so an index beside one cannot be used.
pub fn is_bgzf(data: &[u8]) -> bool {
    if !gzip::is_gzip(data) || data.len() < 12 || data[3] & 0x04 == 0 {
        return false;
    }
    let length = usize::from(u16::from_le_bytes([data[10], data[11]]));
    data.get(12..12 + length)
        .is_some_and(|extra| block_size(extra).is_some())
}

/// The size of a block, from the subfield called BC of its extra field: the
/// size less one, as two bytes.
fn block_size(extra: &[u8]) -> Option<usize> {
    let mut at = 0;
    let mut size = None;
    while at + 4 <= extra.len() {
        let length = usize::from(u16::from_le_bytes([extra[at + 2], extra[at + 3]]));
        if extra[at] == b'B' && extra[at + 1] == b'C' && length == 2 && at + 6 <= extra.len() {
            size = Some(usize::from(u16::from_le_bytes([extra[at + 4], extra[at + 5]])) + 1);
        }
        at += 4 + length;
    }
    size
}

/// A BGZF file, read block by block from wherever an index points.
pub struct Bgzf<R> {
    inner: R,
    block: Vec<u8>,
    at: usize,
    offset: u64,
    next: u64,
}

impl<R: Read + Seek> Bgzf<R> {
    /// Starts at the first block.
    pub fn new(inner: R) -> Self {
        Bgzf {
            inner,
            block: Vec::new(),
            at: 0,
            // No block is loaded yet, and no block is at this offset, so the
            // first seek loads the one it names. At nought, a fresh reader
            // took a seek into its first block for a seek to the end of an
            // empty file there, and refused every place in that block as
            // past its end: the rows of a small bgzipped file share its
            // first block with its header.
            offset: u64::MAX,
            next: 0,
        }
    }

    /// Reads the block at compressed offset `offset`, and says whether there
    /// was one: false at the end of the file.
    fn load(&mut self, offset: u64) -> Result<bool, ReadError> {
        self.inner
            .seek(SeekFrom::Start(offset))
            .map_err(|error| ReadError::whole(error.to_string()))?;
        let mut head = [0u8; 12];
        let got = read_up_to(&mut self.inner, &mut head)?;
        self.offset = offset;
        self.block.clear();
        self.at = 0;
        if got == 0 {
            self.next = offset;
            return Ok(false);
        }
        if got < head.len() || !gzip::is_gzip(&head) || head[3] & 0x04 == 0 {
            return Err(ReadError::whole(
                "not BGZF: a block does not start the way bgzip writes one",
            ));
        }
        let extra_length = usize::from(u16::from_le_bytes([head[10], head[11]]));
        let mut extra = vec![0u8; extra_length];
        if read_up_to(&mut self.inner, &mut extra)? < extra_length {
            return Err(ReadError::whole("a BGZF block is cut short"));
        }
        let size = block_size(&extra)
            .ok_or_else(|| ReadError::whole("a BGZF block does not say its size"))?;
        let before = head.len() + extra.len();
        if size < before {
            return Err(ReadError::whole(
                "a BGZF block is smaller than its own header",
            ));
        }
        let mut whole = Vec::with_capacity(size);
        whole.extend_from_slice(&head);
        whole.extend_from_slice(&extra);
        whole.resize(size, 0);
        if read_up_to(&mut self.inner, &mut whole[before..])? < size - before {
            return Err(ReadError::whole("a BGZF block is cut short"));
        }
        let (data, _) = gzip::member(&whole)?;
        self.block = data;
        self.next = offset + size as u64;
        Ok(true)
    }

    /// Moves to a virtual offset, as an index gives one: the compressed
    /// offset of a block in the high 48 bits and a position within it in the
    /// low 16.
    ///
    /// # Errors
    ///
    /// A block there that is not BGZF or is damaged, and a position past the
    /// end of the block.
    pub fn seek(&mut self, virtual_offset: u64) -> Result<(), ReadError> {
        let (offset, within) = (virtual_offset >> 16, (virtual_offset & 0xffff) as usize);
        if offset != self.offset || (self.block.is_empty() && self.next != offset) {
            self.load(offset)?;
        }
        if within > self.block.len() {
            return Err(ReadError::whole("an index points past the end of a block"));
        }
        self.at = within;
        Ok(())
    }

    /// The virtual offset of the next byte.
    ///
    /// At the end of a block that is the next block at offset nought, which
    /// is how htslib writes it into an index; see the module documentation.
    pub fn tell(&self) -> u64 {
        if self.at == self.block.len() {
            self.next << 16
        } else {
            (self.offset << 16) | self.at as u64
        }
    }

    /// Fills `buf` from as many blocks as it takes, and says how much it
    /// filled: all of it, or less only where the file ends first.
    ///
    /// # Errors
    ///
    /// A block on the way that is not BGZF or is damaged.
    pub fn fill(&mut self, buf: &mut [u8]) -> Result<usize, ReadError> {
        let mut filled = 0;
        while filled < buf.len() {
            if self.at == self.block.len() {
                if !self.load(self.next)? {
                    break;
                }
                continue;
            }
            let n = (buf.len() - filled).min(self.block.len() - self.at);
            buf[filled..filled + n].copy_from_slice(&self.block[self.at..self.at + n]);
            self.at += n;
            filled += n;
        }
        Ok(filled)
    }

    /// The next line into `into`, without its newline and across as many
    /// blocks as it runs over, and whether there was one: false at the end of
    /// the file. A last line with no newline after it is a line.
    ///
    /// The bytes are handed over as they are, a carriage return included,
    /// because the readers that take the text already drop one.
    ///
    /// # Errors
    ///
    /// A block on the way that is not BGZF or is damaged.
    pub fn line(&mut self, into: &mut Vec<u8>) -> Result<bool, ReadError> {
        into.clear();
        let mut any = false;
        loop {
            if self.at == self.block.len() {
                if !self.load(self.next)? {
                    return Ok(any);
                }
                continue;
            }
            let rest = &self.block[self.at..];
            match rest.iter().position(|byte| *byte == b'\n') {
                Some(end) => {
                    into.extend_from_slice(&rest[..end]);
                    self.at += end + 1;
                    return Ok(true);
                }
                None => {
                    into.extend_from_slice(rest);
                    self.at = self.block.len();
                    any = true;
                }
            }
        }
    }
}

/// Reads until `buf` is full or the reader ends, and says how much it read.
fn read_up_to(reader: &mut impl Read, buf: &mut [u8]) -> Result<usize, ReadError> {
    let mut got = 0;
    while got < buf.len() {
        match reader.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(ReadError::whole(error.to_string())),
        }
    }
    Ok(got)
}

/// Blocks of text cut where a test needs them, as bgzip writes them.
#[cfg(test)]
pub(crate) mod fixture {
    /// `data` as BGZF with a block ending at each of `cuts`, and the empty
    /// block bgzip ends a file with. Each block is one stored DEFLATE block,
    /// which is all a test needs: the reader takes any kind of block.
    pub(crate) fn blocks(data: &[u8], cuts: &[usize]) -> Vec<u8> {
        let mut bounds = vec![0];
        bounds.extend_from_slice(cuts);
        bounds.push(data.len());
        let mut out = Vec::new();
        for pair in bounds.windows(2) {
            out.extend(block(&data[pair[0]..pair[1]]));
        }
        out.extend(block(&[]));
        out
    }

    fn block(data: &[u8]) -> Vec<u8> {
        // Header and BC subfield, then a stored block, the CRC and the length.
        let size = 18 + 5 + data.len() + 8;
        let mut out = vec![
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ];
        out.extend_from_slice(&((size - 1) as u16).to_le_bytes());
        out.push(1);
        let length = data.len() as u16;
        out.extend_from_slice(&length.to_le_bytes());
        out.extend_from_slice(&(!length).to_le_bytes());
        out.extend_from_slice(data);
        let mut crc = !0u32;
        for byte in data {
            crc ^= u32::from(*byte);
            for _ in 0..8 {
                crc = if crc & 1 != 0 {
                    0xedb8_8320 ^ (crc >> 1)
                } else {
                    crc >> 1
                };
            }
        }
        out.extend_from_slice(&(!crc).to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
        out
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::blocks;
    use super::*;
    use std::io::Cursor;

    fn lines_of(data: &[u8]) -> Vec<(u64, String)> {
        let mut bgzf = Bgzf::new(Cursor::new(data));
        let mut line = Vec::new();
        let mut out = Vec::new();
        loop {
            let at = bgzf.tell();
            if !bgzf.line(&mut line).unwrap() {
                return out;
            }
            out.push((at, String::from_utf8(line.clone()).unwrap()));
        }
    }

    /// A line that ends where its block ends leaves the reader at the next
    /// block's start, as htslib writes that place, and a line cut across two
    /// blocks is one line.
    #[test]
    fn a_line_that_ends_a_block_leaves_the_reader_at_the_next_block() {
        let text = b"#header\nchr1\t5\nchr1\t9\n";
        // A block of the header, then one ending part way through a row.
        // Each block here is its text and 31 bytes around it, so the second
        // starts at 39 and the third at 74.
        let data = blocks(text, &[8, 12]);
        let (second, third) = (39u64, 74u64);
        assert_eq!(
            lines_of(&data),
            [
                (0, "#header".to_string()),
                (second << 16, "chr1\t5".to_string()),
                ((third << 16) | 3, "chr1\t9".to_string()),
            ]
        );
        let mut bgzf = Bgzf::new(Cursor::new(&data));
        bgzf.seek(second << 16).unwrap();
        let mut line = Vec::new();
        assert!(bgzf.line(&mut line).unwrap());
        assert_eq!(line, b"chr1\t5");
        // And the end of the file is the empty block's end, once read past.
        while bgzf.line(&mut line).unwrap() {}
        assert_eq!(bgzf.tell(), (data.len() as u64) << 16);
    }

    /// A reader made to go straight to a place in the first block goes
    /// there, as a reader of a window does: the rows of a small bgzipped
    /// file begin in the block that holds its header.
    #[test]
    fn a_new_reader_goes_straight_to_a_place_in_the_first_block() {
        let data = blocks(b"#header\nchr1\t5\nchr1\t9\n", &[15]);
        let mut bgzf = Bgzf::new(Cursor::new(&data));
        bgzf.seek(8).unwrap();
        let mut line = Vec::new();
        assert!(bgzf.line(&mut line).unwrap());
        assert_eq!(line, b"chr1\t5");
        assert!(bgzf.line(&mut line).unwrap());
        assert_eq!(line, b"chr1\t9");
        // And a place past the end of that block is still refused.
        let mut bgzf = Bgzf::new(Cursor::new(&data));
        assert!(bgzf.seek(16).is_err());
    }

    #[test]
    fn a_last_line_with_no_newline_is_a_line() {
        let data = blocks(b"one\ntwo", &[5]);
        let read: Vec<String> = lines_of(&data).into_iter().map(|(_, line)| line).collect();
        assert_eq!(read, ["one", "two"]);
        assert!(lines_of(&blocks(b"", &[])).is_empty());
        let read: Vec<String> = lines_of(&blocks(b"\n\nx\n", &[1]))
            .into_iter()
            .map(|(_, line)| line)
            .collect();
        assert_eq!(read, ["", "", "x"]);
    }

    /// Bytes are filled across blocks, and less than asked for only at the
    /// end of the file.
    #[test]
    fn fill_reads_across_blocks_and_says_where_the_file_ends() {
        let data = blocks(b"abcdefgh", &[3, 5]);
        let mut bgzf = Bgzf::new(Cursor::new(&data));
        let mut buf = [0u8; 6];
        assert_eq!(bgzf.fill(&mut buf).unwrap(), 6);
        assert_eq!(&buf, b"abcdef");
        assert_eq!(bgzf.fill(&mut buf).unwrap(), 2);
        assert_eq!(&buf[..2], b"gh");
        assert_eq!(bgzf.fill(&mut buf).unwrap(), 0);
    }

    #[test]
    fn is_bgzf_tells_bgzip_from_gzip() {
        assert!(is_bgzf(&blocks(b"chr1\t1\n", &[])));
        // "karyon\n" as `gzip -n` writes it: no extra field at all.
        let gzip = [
            0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xcb, 0x4e, 0x2c, 0xaa,
            0xcc, 0xcf, 0xe3, 0x02, 0x00, 0x0d, 0x64, 0x70, 0x1a, 0x07, 0x00, 0x00, 0x00,
        ];
        assert!(gzip::is_gzip(&gzip) && !is_bgzf(&gzip));
        // An extra field of another subfield, and one cut short.
        let mut other = blocks(b"x", &[]);
        other[12] = b'Z';
        assert!(!is_bgzf(&other));
        assert!(!is_bgzf(&blocks(b"x", &[])[..14]));
        assert!(!is_bgzf(b""));
    }

    #[test]
    fn a_damaged_file_is_an_error_and_never_a_panic() {
        let data = blocks(b"#h\nchr1\t5\nchr1\t9\n", &[3, 7]);
        for cut in 0..data.len() {
            let mut bgzf = Bgzf::new(Cursor::new(&data[..cut]));
            let mut line = Vec::new();
            while let Ok(true) = bgzf.line(&mut line) {}
        }
        for at in 0..data.len() {
            let mut bent = data.clone();
            bent[at] ^= 0x55;
            let mut bgzf = Bgzf::new(Cursor::new(&bent));
            let mut line = Vec::new();
            while let Ok(true) = bgzf.line(&mut line) {}
        }
        let mut bgzf = Bgzf::new(Cursor::new(&data));
        let error = bgzf.seek(0xffff).unwrap_err();
        assert!(error.to_string().contains("past the end"), "{error}");
        let mut bgzf = Bgzf::new(Cursor::new(&b"plain text, not gzip"[..]));
        let error = bgzf.line(&mut Vec::new()).unwrap_err();
        assert!(error.to_string().contains("not BGZF"), "{error}");
    }
}
