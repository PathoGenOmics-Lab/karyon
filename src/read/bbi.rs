//! What a bigWig and a bigBed share, which is everything but what their
//! blocks hold: the header, the index of sequence names, the R-tree that says
//! which blocks hold a window, and the blocks themselves.
//!
//! Both are UCSC's BBI layout, written by kent's library and by every tool
//! built on it. The header is 64 bytes in the byte order of the machine that
//! wrote the file, which the magic number says; then a header per zoom level;
//! then, wherever the header points, a B+ tree from each sequence's name to a
//! number and a length, and an R-tree from each stretch of the file to the
//! block holding it. A block is compressed with zlib when the header gives the
//! size the largest one inflates to, and stored as it is when that is nought.
//!
//! Reading a window is a walk down the R-tree to the blocks that overlap it,
//! and nothing else of the file is read: a window of ten thousand bases of a
//! whole genome's bigWig reads a few kilobytes of it. kent pads every node of
//! both trees to the block size the writer was given, 256 entries by default,
//! which is why four spans written that way are 13 KB, and why a reader walks
//! nodes by the count each one says rather than by their length.
//!
//! # What is refused
//!
//! An offset or a size that reaches past the end of the file, a block that
//! does not inflate or inflates to more than the header allows, and a tree
//! that leads back into a node it has been through, round to itself or into
//! the middle of another, which a damaged offset can do and which would
//! otherwise be walked for ever. Nothing in a damaged file
//! panics, and nothing it says is allocated before the file is seen to hold
//! it.

use std::collections::BTreeMap;
use std::io::{Read, Seek};

use super::bytes::{Bytes, File};
use super::{gzip, ReadError};

/// The magic number a bigWig starts with, read in its own byte order.
pub(crate) const BIGWIG: u32 = 0x888F_FC26;
/// The magic number a bigBed starts with.
pub(crate) const BIGBED: u32 = 0x8789_F2EB;
/// The B+ tree of sequence names.
const NAMES: u32 = 0x78CA_8C91;
/// The R-tree of blocks.
const BLOCKS: u32 = 0x2468_ACE0;

/// One zoom level: the data summarised in bins of `reduction` bases.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Zoom {
    pub(crate) reduction: u32,
    /// Where its R-tree is.
    pub(crate) index: u64,
}

/// A bigWig or a bigBed, its header read.
pub(crate) struct Bbi<R> {
    file: File<R>,
    big: bool,
    /// The zoom levels, finest first.
    pub(crate) zooms: Vec<Zoom>,
    names: u64,
    /// Where the R-tree over the data as written is.
    pub(crate) full: u64,
    /// How many of a bigBed's columns its autoSql says are BED's own.
    pub(crate) defined: u16,
    /// The most a block inflates to, or nought where blocks are stored.
    inflated: u32,
}

/// A sequence, as the B+ tree names it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Named {
    pub(crate) name: String,
    pub(crate) id: u32,
    pub(crate) length: u32,
}

/// The stretch a query asks for, by sequence number and position.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Over {
    pub(crate) id: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

/// The nodes a walk down a tree has been through, each as where it starts
/// and where it ends.
#[derive(Default)]
struct Walk {
    nodes: BTreeMap<u64, u64>,
}

impl Walk {
    /// Takes in the node from `start` to `end`, unless it overlaps one taken
    /// in already.
    fn visit(&mut self, start: u64, end: u64) -> bool {
        let before = self.nodes.range(..=start).next_back();
        let after = self.nodes.range(start..).next();
        if before.is_some_and(|(_, last)| *last > start)
            || after.is_some_and(|(first, _)| *first < end)
        {
            return false;
        }
        self.nodes.insert(start, end);
        true
    }
}

/// How deep a tree may go. A tree of 2^32 sequences two to a node is 32
/// levels, so anything deeper is a damaged offset going round.
const DEEPEST: usize = 64;

impl<R: Read + Seek> Bbi<R> {
    /// Reads the header of a file that starts with `magic`, either way round.
    pub(crate) fn open(reader: R, magic: u32, what: &'static str) -> Result<Self, ReadError> {
        let mut file = File::new(reader, what)?;
        let head = file.read_up_to(0, 64)?;
        if head.len() < 64 {
            return Err(ReadError::whole(format!("not a {what}: it is too short")));
        }
        let first = u32::from_le_bytes([head[0], head[1], head[2], head[3]]);
        let big = if first == magic {
            false
        } else if first.swap_bytes() == magic {
            true
        } else {
            return Err(ReadError::whole(format!(
                "not a {what}: it does not start with a {what}'s magic"
            )));
        };
        let mut bytes = Bytes::new(&head[4..], big, what);
        let _version = bytes.u16()?;
        let levels = bytes.u16()?;
        let names = bytes.u64()?;
        let _data = bytes.u64()?;
        let full = bytes.u64()?;
        let _fields = bytes.u16()?;
        let defined = bytes.u16()?;
        let _auto_sql = bytes.u64()?;
        let _summary = bytes.u64()?;
        let inflated = bytes.u32()?;
        let levels = file.read_at(64, u64::from(levels) * 24)?;
        let mut bytes = Bytes::new(&levels, big, what);
        let mut zooms = Vec::with_capacity(levels.len() / 24);
        while bytes.left() > 0 {
            let reduction = bytes.u32()?;
            let _reserved = bytes.u32()?;
            let _data = bytes.u64()?;
            let index = bytes.u64()?;
            zooms.push(Zoom { reduction, index });
        }
        // kent writes them finest first, and a reader that picks one by its
        // reduction should not depend on it.
        zooms.sort_by_key(|zoom| zoom.reduction);
        Ok(Bbi {
            file,
            big,
            zooms,
            names,
            full,
            defined,
            inflated,
        })
    }

    /// What the file is, `bigWig` or `bigBed`, for a message.
    pub(crate) fn what(&self) -> &'static str {
        self.file.what
    }

    /// The sequence called `name`, or a refusal naming the ones there are,
    /// the first dozen of them, as a BAM's is.
    pub(crate) fn named(&mut self, name: &str) -> Result<Named, ReadError> {
        if let Some(named) = self.find(name)? {
            return Ok(named);
        }
        let held = self.sequences()?;
        let mut names: Vec<&str> = held
            .iter()
            .take(12)
            .map(|named| named.name.as_str())
            .collect();
        if held.len() > names.len() {
            names.push("and more");
        }
        Err(ReadError::whole(format!(
            "the {} has no sequence called {name}; it has {}",
            self.what(),
            if names.is_empty() {
                "none".to_string()
            } else {
                names.join(", ")
            }
        )))
    }

    /// A node of either tree: whether it is a leaf, and the bytes of its
    /// items, `leaf` bytes each in a leaf and `inner` bytes each otherwise.
    fn node(
        &mut self,
        at: u64,
        leaf: u64,
        inner: u64,
        walk: &mut Walk,
    ) -> Result<(bool, Vec<u8>), ReadError> {
        let head = self.file.read_at(at, 4)?;
        let mut bytes = Bytes::new(&head, self.big, self.what());
        let is_leaf = bytes.u8()? != 0;
        let _reserved = bytes.u8()?;
        let count = u64::from(bytes.u16()?);
        let length = count.saturating_mul(if is_leaf { leaf } else { inner });
        // The nodes of a tree do not overlap, so a node that does was reached
        // by a path no sound file has: back round to itself, or into another,
        // and with a few children each, either one is walked for ever or for
        // as long as there are ways to go.
        if !walk.visit(at, at.saturating_add(4).saturating_add(length)) {
            return Err(self.damaged("an index leads back into a node it has been through"));
        }
        let items = self.file.read_at(at.saturating_add(4), length)?;
        Ok((is_leaf, items))
    }

    fn damaged(&self, reason: &str) -> ReadError {
        ReadError::whole(format!("the {} is damaged: {reason}", self.what()))
    }

    /// The header of the B+ tree of names: its key size, and where its root
    /// node is.
    fn names_root(&mut self) -> Result<(u64, u64), ReadError> {
        let head = self.file.read_at(self.names, 32)?;
        let mut bytes = Bytes::new(&head, self.big, self.what());
        if bytes.u32()? != NAMES {
            return Err(self.damaged("its index of sequence names is not where it says"));
        }
        let _block = bytes.u32()?;
        let key = u64::from(bytes.u32()?);
        let value = u64::from(bytes.u32()?);
        if value != 8 {
            return Err(self.damaged("a sequence's entry is not a number and a length"));
        }
        // The size every name is padded to, which a name looked up is padded
        // to as well: one byte of the size flipped asked for four gigabytes
        // of noughts, and took half a minute to be refused.
        if key == 0 || key > self.file.size || usize::try_from(key).is_err() {
            return Err(self.damaged("its names are longer than the file"));
        }
        Ok((key, self.names.saturating_add(32)))
    }

    /// Every sequence the file names, in the order of its index, which is
    /// the order of their names.
    pub(crate) fn sequences(&mut self) -> Result<Vec<Named>, ReadError> {
        let (key, root) = self.names_root()?;
        let mut found = Vec::new();
        let mut walk = Walk::default();
        let mut stack = vec![(root, 0usize)];
        while let Some((at, depth)) = stack.pop() {
            if depth > DEEPEST {
                return Err(self.damaged("its index of sequence names runs too deep"));
            }
            let (leaf, items) = self.node(at, key + 8, key + 8, &mut walk)?;
            let mut bytes = Bytes::new(&items, self.big, self.what());
            let mut children = Vec::new();
            while bytes.left() > 0 {
                let name = bytes.take(key as usize)?;
                if leaf {
                    let id = bytes.u32()?;
                    let length = bytes.u32()?;
                    found.push(Named {
                        name: unpadded(name),
                        id,
                        length,
                    });
                } else {
                    children.push((bytes.u64()?, depth + 1));
                }
            }
            // Walked in order: the stack takes the last child first.
            stack.extend(children.into_iter().rev());
        }
        Ok(found)
    }

    /// The sequence called `name`, or `None` where the file names none so.
    pub(crate) fn find(&mut self, name: &str) -> Result<Option<Named>, ReadError> {
        let (key, mut at) = self.names_root()?;
        let wanted = name.as_bytes();
        if wanted.len() as u64 > key {
            return Ok(None);
        }
        let mut padded = wanted.to_vec();
        padded.resize(key as usize, 0);
        let mut walk = Walk::default();
        for _ in 0..=DEEPEST {
            let (leaf, items) = self.node(at, key + 8, key + 8, &mut walk)?;
            let mut bytes = Bytes::new(&items, self.big, self.what());
            let mut child = None;
            while bytes.left() > 0 {
                let here = bytes.take(key as usize)?;
                if leaf {
                    let id = bytes.u32()?;
                    let length = bytes.u32()?;
                    if here == padded.as_slice() {
                        return Ok(Some(Named {
                            name: name.to_string(),
                            id,
                            length,
                        }));
                    }
                } else {
                    let offset = bytes.u64()?;
                    // Down the last child whose first key is not past the
                    // name, as kent's own reader goes.
                    if child.is_some() && here > padded.as_slice() {
                        break;
                    }
                    child = Some(offset);
                }
            }
            match child {
                Some(offset) if !leaf => at = offset,
                _ => return Ok(None),
            }
        }
        Err(self.damaged("its index of sequence names runs too deep"))
    }

    /// The blocks the R-tree at `index` says hold data over `over`, or every
    /// block for `None`, each as its offset and its size, in the file's order.
    pub(crate) fn blocks(
        &mut self,
        index: u64,
        over: Option<Over>,
    ) -> Result<Vec<(u64, u64)>, ReadError> {
        let head = self.file.read_at(index, 48)?;
        if Bytes::new(&head, self.big, self.what()).u32()? != BLOCKS {
            return Err(self.damaged("its index of blocks is not where it says"));
        }
        let mut found = Vec::new();
        let mut walk = Walk::default();
        let mut stack = vec![(index.saturating_add(48), 0usize)];
        while let Some((at, depth)) = stack.pop() {
            if depth > DEEPEST {
                return Err(self.damaged("its index of blocks runs too deep"));
            }
            // A leaf's item is the stretch it covers and its block's offset
            // and size, 32 bytes; an inner node's is the stretch and where
            // the node under it is, 24.
            let (leaf, items) = self.node(at, 32, 24, &mut walk)?;
            let size = if leaf { 32 } else { 24 };
            let mut bytes = Bytes::new(&items, self.big, self.what());
            let mut children = Vec::new();
            for _ in 0..items.len() / size {
                let first = (bytes.u32()?, bytes.u32()?);
                let last = (bytes.u32()?, bytes.u32()?);
                let offset = bytes.u64()?;
                let length = if leaf { Some(bytes.u64()?) } else { None };
                let overlaps = over.map_or(true, |over| {
                    first < (over.id, over.end) && last > (over.id, over.start)
                });
                if !overlaps {
                    continue;
                }
                match length {
                    Some(length) => found.push((offset, length)),
                    None => children.push((offset, depth + 1)),
                }
            }
            stack.extend(children.into_iter().rev());
        }
        Ok(found)
    }

    /// A block's bytes, inflated where the file compresses its blocks.
    pub(crate) fn block(&mut self, offset: u64, size: u64) -> Result<Vec<u8>, ReadError> {
        let stored = self.file.read_at(offset, size)?;
        if self.inflated == 0 {
            return Ok(stored);
        }
        gzip::zlib(&stored, self.inflated as usize)
            .map_err(|error| self.damaged(&format!("a block of it does not inflate: {error}")))
    }

    /// Numbers out of a block, in the file's byte order.
    pub(crate) fn bytes<'a>(&self, data: &'a [u8]) -> Bytes<'a> {
        Bytes::new(data, self.big, self.what())
    }
}

/// A name out of the B+ tree, without the noughts that pad it to the key's
/// size.
fn unpadded(key: &[u8]) -> String {
    let end = key.iter().position(|byte| *byte == 0).unwrap_or(key.len());
    String::from_utf8_lossy(&key[..end]).into_owned()
}

/// A BBI file written by hand, a bigWig or a bigBed, for what kent's own
/// tools do not write.
#[cfg(test)]
pub(crate) mod fixture {
    /// Numbers written in one byte order, big-endian where it says so.
    pub(crate) struct Numbers(pub(crate) bool);

    /// A block of a file of several sequences: the sequence and base it
    /// starts at, the sequence and base it ends at, and its bytes.
    pub(crate) type Block = ((u32, u32), (u32, u32), Vec<u8>);

    impl Numbers {
        pub(crate) fn u16(&self, out: &mut Vec<u8>, value: u16) {
            out.extend(if self.0 {
                value.to_be_bytes()
            } else {
                value.to_le_bytes()
            });
        }

        pub(crate) fn u32(&self, out: &mut Vec<u8>, value: u32) {
            out.extend(if self.0 {
                value.to_be_bytes()
            } else {
                value.to_le_bytes()
            });
        }

        pub(crate) fn u64(&self, out: &mut Vec<u8>, value: u64) {
            out.extend(if self.0 {
                value.to_be_bytes()
            } else {
                value.to_le_bytes()
            });
        }

        pub(crate) fn f32(&self, out: &mut Vec<u8>, value: f32) {
            self.u32(out, value.to_bits());
        }

        /// A file of one sequence, `magic` its kind and `defined` the columns
        /// of a bigBed that are BED's own, holding `blocks`, each the first
        /// and last base it covers and its bytes, stored as they are, under
        /// one leaf of an index, and no zoom levels.
        pub(crate) fn file(
            &self,
            magic: u32,
            defined: u16,
            name: &str,
            length: u32,
            blocks: &[(u32, u32, Vec<u8>)],
        ) -> Vec<u8> {
            let blocks: Vec<Block> = blocks
                .iter()
                .map(|(first, last, block)| ((0, *first), (0, *last), block.clone()))
                .collect();
            self.genome(magic, defined, &[(name, length)], &blocks)
        }

        /// A file of `sequences`, each its name and its length, numbered
        /// from nought in the order given as kent numbers them in the order
        /// of their names, and all under one leaf of the index of names,
        /// holding `blocks`.
        pub(crate) fn genome(
            &self,
            magic: u32,
            defined: u16,
            sequences: &[(&str, u32)],
            blocks: &[Block],
        ) -> Vec<u8> {
            // The sequence index sits after the header, then the data, then
            // the block index.
            let names_at = 64u64;
            let key = sequences
                .iter()
                .map(|(name, _)| name.len())
                .max()
                .unwrap_or(1);
            let mut names = Vec::new();
            for number in [0x78CA_8C91, sequences.len() as u32, key as u32, 8] {
                self.u32(&mut names, number);
            }
            self.u64(&mut names, sequences.len() as u64);
            self.u64(&mut names, 0);
            names.extend([1, 0]);
            self.u16(&mut names, sequences.len() as u16);
            for (id, (name, length)) in sequences.iter().enumerate() {
                names.extend(name.as_bytes());
                names.extend(std::iter::repeat(0).take(key - name.len()));
                self.u32(&mut names, id as u32);
                self.u32(&mut names, *length);
            }
            let data_at = names_at + names.len() as u64;
            let mut data = Vec::new();
            self.u64(&mut data, blocks.len() as u64);
            let mut leaves = Vec::new();
            for (first, last, block) in blocks {
                leaves.push((
                    *first,
                    *last,
                    data_at + data.len() as u64,
                    block.len() as u64,
                ));
                data.extend(block);
            }
            let index_at = data_at + data.len() as u64;
            let mut index = Vec::new();
            self.u32(&mut index, 0x2468_ACE0);
            self.u32(&mut index, leaves.len() as u32);
            self.u64(&mut index, leaves.len() as u64);
            let (first, last) = (
                leaves.first().map_or((0, 0), |leaf| leaf.0),
                leaves.last().map_or((0, 0), |leaf| leaf.1),
            );
            for number in [first.0, first.1, last.0, last.1] {
                self.u32(&mut index, number);
            }
            self.u64(&mut index, index_at);
            self.u32(&mut index, 1);
            self.u32(&mut index, 0);
            index.extend([1, 0]);
            self.u16(&mut index, leaves.len() as u16);
            for (first, last, offset, size) in &leaves {
                for number in [first.0, first.1, last.0, last.1] {
                    self.u32(&mut index, number);
                }
                self.u64(&mut index, *offset);
                self.u64(&mut index, *size);
            }
            let mut out = Vec::new();
            self.u32(&mut out, magic);
            self.u16(&mut out, 4);
            self.u16(&mut out, 0);
            for offset in [names_at, data_at, index_at] {
                self.u64(&mut out, offset);
            }
            self.u16(&mut out, defined);
            self.u16(&mut out, defined);
            self.u64(&mut out, 0);
            self.u64(&mut out, 0);
            self.u32(&mut out, 0);
            self.u64(&mut out, 0);
            out.extend(names);
            out.extend(data);
            out.extend(index);
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// The four blocks of the small bigWig hold chr1 from 10 to 30 and from
    /// 99 to 900, chr2 from 5 to 10 and chr3 from 0 to 50, and a window is
    /// handed the blocks it overlaps, by sequence and by position, half-open.
    #[test]
    fn a_window_is_handed_the_blocks_it_overlaps() {
        let bytes = include_bytes!("fixtures/signal.bw");
        let mut file = Bbi::open(Cursor::new(&bytes[..]), BIGWIG, "bigWig").unwrap();
        let full = file.full;
        let mut count = |id: u32, start: u32, end: u32| {
            file.blocks(full, Some(Over { id, start, end }))
                .unwrap()
                .len()
        };
        assert_eq!(count(0, 30, 99), 0);
        assert_eq!(count(0, 29, 30), 1);
        assert_eq!(count(0, 29, 100), 2);
        assert_eq!(count(0, 899, 5000), 1);
        assert_eq!(count(1, 0, 5), 0);
        assert_eq!(count(1, 0, 500), 1);
        assert_eq!(count(2, 0, 60), 1);
        assert_eq!(count(3, 0, 60), 0);
        assert_eq!(file.blocks(full, None).unwrap().len(), 4);
    }

    #[test]
    fn a_sequence_is_found_by_its_name_through_inner_nodes() {
        let bytes = include_bytes!("fixtures/signal.bw");
        let mut file = Bbi::open(Cursor::new(&bytes[..]), BIGWIG, "bigWig").unwrap();
        for (name, id, length) in [("chr1", 0, 1000), ("chr2", 1, 500), ("chr3", 2, 60)] {
            let named = file.find(name).unwrap().unwrap();
            assert_eq!((named.id, named.length), (id, length), "{name}");
        }
        for absent in ["chr0", "chr10", "chr4", "chrX", "", "c"] {
            assert_eq!(file.find(absent).unwrap(), None, "{absent}");
        }
    }
}
