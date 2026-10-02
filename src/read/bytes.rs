//! Numbers read out of a binary file in the order it was written in, and
//! stretches of it read by their offset, never past its end.
//!
//! The UCSC formats, bigWig, bigBed and 2bit, write every number in the order
//! of the machine that wrote them, and say which by how their magic number
//! reads, so [`Bytes`] reads either order. A file found anywhere is a file a
//! damaged copy of may be found as, and the offsets and counts in one are what
//! a damaged copy gets wrong, so nothing here believes them: a stretch is read
//! only where the file holds all of it, which caps whatever is allocated for
//! it at the bytes the file has left, and every offset is added checked. An
//! offset near the top of `u64` wraps unchecked, and with overflow checks on,
//! as this crate's release builds have them, the wrap is a panic rather than a
//! wrong number.

use std::io::{Read, Seek, SeekFrom};

use super::ReadError;

/// Bytes read from the front as numbers, in the byte order a file says.
pub(crate) struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
    big: bool,
    what: &'static str,
}

impl<'a> Bytes<'a> {
    /// `data`, read as `what` writes it: big-endian where `big` says so.
    pub(crate) fn new(data: &'a [u8], big: bool, what: &'static str) -> Self {
        Bytes {
            data,
            at: 0,
            big,
            what,
        }
    }

    /// The next `n` bytes.
    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8], ReadError> {
        let end = self
            .at
            .checked_add(n)
            .filter(|end| *end <= self.data.len())
            .ok_or_else(|| short(self.what))?;
        let taken = &self.data[self.at..end];
        self.at = end;
        Ok(taken)
    }

    /// The bytes up to the next nought, which is passed over, as a string
    /// that ends with one is stored.
    pub(crate) fn until_nought(&mut self) -> Result<&'a [u8], ReadError> {
        let length = self.data[self.at..]
            .iter()
            .position(|byte| *byte == 0)
            .ok_or_else(|| short(self.what))?;
        let taken = self.take(length)?;
        self.at += 1;
        Ok(taken)
    }

    /// How many bytes are left to read.
    pub(crate) fn left(&self) -> usize {
        self.data.len() - self.at
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], ReadError> {
        let mut out = [0u8; N];
        out.copy_from_slice(self.take(N)?);
        if self.big {
            out.reverse();
        }
        Ok(out)
    }

    pub(crate) fn u8(&mut self) -> Result<u8, ReadError> {
        Ok(self.take(1)?[0])
    }

    pub(crate) fn u16(&mut self) -> Result<u16, ReadError> {
        self.array().map(u16::from_le_bytes)
    }

    pub(crate) fn u32(&mut self) -> Result<u32, ReadError> {
        self.array().map(u32::from_le_bytes)
    }

    pub(crate) fn u64(&mut self) -> Result<u64, ReadError> {
        self.array().map(u64::from_le_bytes)
    }

    pub(crate) fn f32(&mut self) -> Result<f32, ReadError> {
        self.array().map(f32::from_le_bytes)
    }
}

/// A file cut short, as `what` says it.
pub(crate) fn short(what: &str) -> ReadError {
    ReadError::whole(format!("the {what} is cut short"))
}

/// A file read by offset, which knows how long it is.
pub(crate) struct File<R> {
    reader: R,
    /// How many bytes the file holds.
    pub(crate) size: u64,
    /// What the file is, for a message.
    pub(crate) what: &'static str,
}

impl<R: Read + Seek> File<R> {
    /// The file `reader` reads, measured once.
    pub(crate) fn new(mut reader: R, what: &'static str) -> Result<Self, ReadError> {
        let size = reader
            .seek(SeekFrom::End(0))
            .map_err(|error| ReadError::whole(error.to_string()))?;
        Ok(File { reader, size, what })
    }

    /// The `length` bytes at `offset`, refused where the file does not hold
    /// all of them, before anything is allocated for them.
    pub(crate) fn read_at(&mut self, offset: u64, length: u64) -> Result<Vec<u8>, ReadError> {
        let inside = offset
            .checked_add(length)
            .is_some_and(|end| end <= self.size);
        let length = usize::try_from(length)
            .ok()
            .filter(|_| inside)
            .ok_or_else(|| {
                ReadError::whole(format!(
                    "the {} points past its own end, so it is damaged or cut short",
                    self.what
                ))
            })?;
        self.reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| ReadError::whole(error.to_string()))?;
        let mut out = vec![0u8; length];
        self.reader
            .read_exact(&mut out)
            .map_err(|_| short(self.what))?;
        Ok(out)
    }

    /// Up to `length` bytes at `offset`, as many as the file holds: for
    /// reading a header whose length is not known until it is read.
    pub(crate) fn read_up_to(&mut self, offset: u64, length: u64) -> Result<Vec<u8>, ReadError> {
        let length = length.min(self.size.saturating_sub(offset));
        self.read_at(offset.min(self.size), length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn numbers_read_in_either_order() {
        let data = [1, 0, 0, 0, 0, 0, 0x80, 0x3f];
        let mut little = Bytes::new(&data, false, "file");
        assert_eq!(little.u32().unwrap(), 1);
        assert_eq!(little.f32().unwrap(), 1.0);
        let mut big = Bytes::new(&data, true, "file");
        assert_eq!(big.u32().unwrap(), 1 << 24);
        assert_eq!(big.u16().unwrap(), 0);
        assert_eq!(big.left(), 2);
        let error = big.u32().unwrap_err();
        assert_eq!(error.to_string(), "the file is cut short");
    }

    /// An offset or a length near the top of `u64` is refused, never added
    /// past it, and nothing is allocated for a stretch the file has not got.
    #[test]
    fn a_stretch_past_the_end_is_refused_before_it_is_allocated() {
        let mut file = File::new(Cursor::new(vec![7u8; 16]), "bigWig").unwrap();
        assert_eq!(file.read_at(8, 8).unwrap(), vec![7; 8]);
        for (offset, length) in [(8, 9), (u64::MAX, 2), (2, u64::MAX), (17, 0)] {
            let error = file.read_at(offset, length).unwrap_err();
            assert!(error.to_string().contains("past its own end"), "{error}");
        }
        assert_eq!(file.read_up_to(12, 100).unwrap(), vec![7; 4]);
        assert!(file.read_up_to(u64::MAX, 100).unwrap().is_empty());
    }
}
