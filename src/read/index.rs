//! The indexes that say where in a BGZF file the records over a window lie:
//! BAI beside a BAM, TBI beside a bgzipped text file, and CSI beside either.
//!
//! All three file every record under the smallest bin of a fixed tree that
//! holds the whole of it, and keep for each bin the stretches of the file its
//! records lie in, as pairs of virtual offsets (see [`bgzf`](super::bgzf)). A
//! record over a window is in one of the bins over that window on some level
//! of the tree: a short one deep down, one starting far to the left and
//! reaching into the window near the root. So a window is read by asking every
//! level, which is a handful of bins, and reading their stretches in order.
//!
//! The three differ in the size of the tree and in how they say where the
//! first record over a place is:
//!
//! - BAI and TBI have leaves of 2^14 bases and five levels over them, so the
//!   root spans 2^29 bases. A linear index gives, for every 16 KiB window of
//!   a sequence, the first record over it.
//! - CSI writes its leaf size (`min_shift`) and its number of levels (`depth`)
//!   in the file, and gives every bin the first record over it in place of a
//!   linear index. That is what lets it index a sequence past 2^29 bases. The
//!   depth is whatever the writer needed and has to be read: samtools wrote 0
//!   for a BAM of a 1,000-base sequence, bcftools 6 and `tabix -C` 8.
//!
//! TBI, and CSI written for a text file, also say which columns of a row hold
//! its sequence, start and end ([`Columns`]), and name the sequences, which a
//! BAM's own header does for its BAI.
//!
//! Nothing here opens a file. [`parse`] takes an index's bytes, as samtools,
//! tabix and bcftools write them or already taken out of their wrapper, and
//! [`Index::chunks`] answers with offsets for [`Bgzf`](super::bgzf::Bgzf) to
//! go to.
//!
//! # What comes out
//!
//! For a window, 0-based and half-open as everywhere else, the stretches of the
//! file that can hold a record over it, in file order and merged where they
//! touch. They hold every such record and some that are not: a bin's stretch
//! holds every record filed under it, and htslib folds the bins of a small
//! file into their parents, so on a file of a few blocks a stretch can begin
//! at the sequence's first record. A reader drops what lies outside the
//! window, as it drops the rows outside the window of a whole file. Each
//! stretch begins no earlier than the first record over the window's 16 KiB
//! leaf, where the index says where that is, since nothing before it can
//! reach the window.
//!
//! # The pseudo-bin
//!
//! One bin past the last of the tree, 37,450 at BAI's depth and
//! `((1 << 3 * (depth + 1)) - 1) / 7 + 1` at any other, is not a bin. It holds
//! two pairs that are not stretches: where a sequence's first record is and
//! where its last ends, and how many records it has with a position and
//! without one. For text, the second count is nought and the first is the
//! sequence's rows. [`Index::summary`] gives them. The format leaves the bin
//! out where a writer chooses, so every use of it has an answer without it.
//!
//! # What is refused
//!
//! Bytes that start with none of the three magics, an index cut short, a count
//! below nought, a tree too deep for its bins to be numbered, a column below
//! nought, a comment character past a byte, a TBI whose names do not match its
//! count of sequences, and a CSI whose extra bytes are not tabix's columns.
//! None of it panics, and no count read from the file sets an allocation
//! larger than the bytes left could fill.
//!
//! # Who reads it
//!
//! [`bam`](super::bam) reads a BAM through its BAI, and through a CSI the same
//! way when a caller hands one. The command line looks for the `.bai` alone,
//! as `reads.bam.bai` or `reads.bai`. htslib looks for a `.csi` first, so
//! looking for one changes which index an existing BAM is read through, and
//! that comes with BCF, whose only index is a CSI, with htslib's order and a
//! test of a BAM that has both. The rows of a bgzipped text file are found
//! through the same chunks: [`tabix`](super::tabix) reads them from each with
//! [`Bgzf::line`](super::bgzf::Bgzf::line), and the [`Columns`] the index
//! keeps say where each row lies, so it stops at the first past the window.
//! The command line looks for a `.csi` before a `.tbi` beside such a file, as
//! htslib does.

use std::collections::BTreeMap;
use std::ops::RangeInclusive;

use super::{gzip, ReadError};

/// Which of the three an index is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A BAM's index, as `samtools index` writes it.
    Bai,
    /// A bgzipped text file's index, as `tabix` writes it.
    Tbi,
    /// Either, with a tree of the writer's size, as `samtools index -c`,
    /// `tabix -C` and `bcftools index` write it.
    Csi,
}

/// The kind of text file tabix was told it was indexing, which fixes how a
/// row's end is found where no column holds it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// Columns as tabix was told them, BED and GFF3 among them.
    Generic,
    /// SAM, whose end is worked out from its CIGAR.
    Sam,
    /// VCF, whose end is its position and the length of its reference
    /// allele, or what `END=` says.
    Vcf,
    /// Graph alignments, as minigraph and vg write them.
    Gaf,
}

/// Which columns of a row hold its sequence, start and end, as tabix was told
/// them, and which lines before the rows are not rows.
///
/// Column numbers count from one, as tabix counts them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Columns {
    /// What the file was indexed as.
    pub preset: Preset,
    /// Whether a start counts from nought and an interval is half-open, as
    /// `-p bed` and `-0` index it, rather than from one and inclusive.
    pub zero_based: bool,
    /// The column holding the sequence's name.
    pub sequence: usize,
    /// The column holding the start.
    pub start: usize,
    /// The column holding the end, or `None` where the end comes from the
    /// start, or for VCF and SAM from the row itself.
    pub end: Option<usize>,
    /// The character a header or comment line starts with.
    pub comment: char,
    /// How many lines at the top are a header whatever they start with.
    pub skip: usize,
}

/// What the pseudo-bin says about one sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Summary {
    /// The virtual offset of the sequence's first record.
    pub first: u64,
    /// The virtual offset just past its last.
    pub end: u64,
    /// How many of its records have a position: for text, its rows.
    pub placed: u64,
    /// How many do not, as an unmapped read placed by its mate. Nought for
    /// text.
    pub unplaced: u64,
}

/// A BAI, TBI or CSI index.
#[derive(Debug, Clone)]
pub struct Index {
    kind: Kind,
    min_shift: u32,
    depth: u32,
    columns: Option<Columns>,
    names: Vec<String>,
    references: Vec<Reference>,
    unplaced: Option<u64>,
}

/// An index of no sequences, as a BAI of an empty BAM would be.
impl Default for Index {
    fn default() -> Self {
        Index {
            kind: Kind::Bai,
            min_shift: 14,
            depth: 5,
            columns: None,
            names: Vec::new(),
            references: Vec::new(),
            unplaced: None,
        }
    }
}

/// One sequence's bins.
#[derive(Debug, Clone, Default)]
struct Reference {
    bins: BTreeMap<u32, Bin>,
    linear: Vec<u64>,
    summary: Option<Summary>,
}

/// One bin: the first record over it, which only a CSI says, and the
/// stretches its records lie in.
#[derive(Debug, Clone, Default)]
struct Bin {
    first: u64,
    chunks: Vec<(u64, u64)>,
}

/// The pseudo-bin of a tree `depth` levels deep: the first bin of the level
/// under its deepest.
fn pseudo_bin(depth: u32) -> u32 {
    (((1u64 << (3 * (depth + 1))) - 1) / 7 + 1) as u32
}

/// The bins on each level of a tree that can hold a record over
/// `[start, end)`, as the SAM specification's `reg2bins` lists them: one
/// range a level, from the root down. A window past the last position the
/// tree has room for holds nothing it could have filed.
fn levels(min_shift: u32, depth: u32, start: u64, end: u64) -> Vec<RangeInclusive<u32>> {
    let reach = 1u64 << (min_shift + 3 * depth);
    if start >= reach {
        return Vec::new();
    }
    // A window of no bases is asked about as the one base at its start.
    let last = end.saturating_sub(1).max(start).min(reach - 1);
    (0..=depth)
        .map(|level| {
            let shift = min_shift + 3 * (depth - level);
            let first = ((1u64 << (3 * level)) - 1) / 7;
            (first + (start >> shift)) as u32..=(first + (last >> shift)) as u32
        })
        .collect()
}

/// An index's bytes, read from the front.
struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
    /// What the index is called in a message, as "the BAM index".
    noun: &'static str,
}

impl<'a> Bytes<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], ReadError> {
        let bytes = self
            .at
            .checked_add(n)
            .and_then(|end| self.data.get(self.at..end))
            .ok_or_else(|| ReadError::whole(format!("{} is cut short", self.noun)))?;
        self.at += n;
        Ok(bytes)
    }

    fn left(&self) -> usize {
        self.data.len() - self.at
    }

    fn u32(&mut self) -> Result<u32, ReadError> {
        let bytes = self.take(4)?;
        Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn i32(&mut self) -> Result<i32, ReadError> {
        Ok(self.u32()? as i32)
    }

    fn u64(&mut self) -> Result<u64, ReadError> {
        let bytes = self.take(8)?;
        let mut eight = [0u8; 8];
        eight.copy_from_slice(bytes);
        Ok(u64::from_le_bytes(eight))
    }

    /// A count the format stores as an `int32`, which a sound index never
    /// writes below nought.
    fn count(&mut self, what: &str) -> Result<usize, ReadError> {
        let value = self.i32()?;
        usize::try_from(value).map_err(|_| {
            ReadError::whole(format!(
                "{} has an impossible count of {what}, {value}",
                self.noun
            ))
        })
    }
}

/// Reads an index of any of the three kinds, told by its magic: BAI as
/// samtools writes it, and TBI and CSI in the BGZF they are written in or
/// already taken out of it.
///
/// # Errors
///
/// Bytes that are none of the three, or one that is damaged; see the module
/// documentation.
pub fn parse(data: &[u8]) -> Result<Index, ReadError> {
    let inflated;
    let data = if gzip::is_gzip(data) {
        inflated = gzip::decompress(data)?;
        &inflated[..]
    } else {
        data
    };
    let (kind, noun) = match data.get(..4) {
        Some(b"BAI\x01") => (Kind::Bai, "the BAM index"),
        Some(b"TBI\x01") => (Kind::Tbi, "the tabix index"),
        Some(b"CSI\x01") => (Kind::Csi, "the CSI index"),
        Some(_) => {
            return Err(ReadError::whole(
                "not an index: it starts with none of BAI's, TBI's and CSI's magic",
            ))
        }
        None => return Err(ReadError::whole("the index is cut short")),
    };
    let mut bytes = Bytes { data, at: 4, noun };
    let mut index = Index {
        kind,
        ..Index::default()
    };
    let references = match kind {
        Kind::Bai => bytes.count("sequences")?,
        Kind::Tbi => {
            let references = bytes.count("sequences")?;
            let (columns, names) = columns(&mut bytes)?;
            index.columns = Some(columns);
            index.names = names;
            references
        }
        Kind::Csi => {
            let min_shift = bytes.count("bases a leaf holds")?;
            let depth = bytes.count("levels")?;
            // A bin is numbered in 32 bits and a position shifted in 64, which
            // a tree eleven levels deep or reaching past 2^62 bases outgrows.
            if depth > 10 || min_shift + 3 * depth > 62 {
                return Err(ReadError::whole(format!(
                    "{noun} has a tree of bins too deep to number, {min_shift} and {depth}"
                )));
            }
            index.min_shift = min_shift as u32;
            index.depth = depth as u32;
            let extra = bytes.count("extra bytes")?;
            let extra = bytes.take(extra)?;
            if !extra.is_empty() {
                let mut aux = Bytes {
                    data: extra,
                    at: 0,
                    noun,
                };
                let (columns, names) = columns(&mut aux).map_err(|_| {
                    ReadError::whole(format!("{noun} holds extra bytes that are not tabix's"))
                })?;
                index.columns = Some(columns);
                index.names = names;
            }
            bytes.count("sequences")?
        }
    };
    if !index.names.is_empty() && index.names.len() != references {
        return Err(ReadError::whole(format!(
            "{noun} names {} sequences and indexes {references}",
            index.names.len()
        )));
    }
    let pseudo = pseudo_bin(index.depth);
    for _ in 0..references {
        let mut reference = Reference::default();
        for _ in 0..bytes.count("bins")? {
            let number = bytes.u32()?;
            let first = if kind == Kind::Csi { bytes.u64()? } else { 0 };
            let n = bytes.count("chunks")?;
            let mut chunks = Vec::with_capacity(n.min(bytes.left() / 16));
            for _ in 0..n {
                chunks.push((bytes.u64()?, bytes.u64()?));
            }
            if number == pseudo {
                if let [(first, end), (placed, unplaced)] = chunks[..] {
                    reference.summary = Some(Summary {
                        first,
                        end,
                        placed,
                        unplaced,
                    });
                }
                continue;
            }
            // A bin written twice keeps both its lists, where the second
            // would otherwise hide the records the first one points to.
            let bin = reference.bins.entry(number).or_insert(Bin {
                first,
                chunks: Vec::new(),
            });
            bin.first = bin.first.min(first);
            bin.chunks.extend(chunks);
        }
        if kind != Kind::Csi {
            let n = bytes.count("linear offsets")?;
            reference.linear = Vec::with_capacity(n.min(bytes.left() / 8));
            for _ in 0..n {
                reference.linear.push(bytes.u64()?);
            }
        }
        index.references.push(reference);
    }
    // How many records have no sequence at all, which the format allows a
    // writer to leave off the end.
    if bytes.left() >= 8 {
        index.unplaced = Some(bytes.u64()?);
    }
    Ok(index)
}

/// TBI's description of a text file: its columns, then the names of its
/// sequences, which a CSI written for text carries in its extra bytes.
fn columns(bytes: &mut Bytes<'_>) -> Result<(Columns, Vec<String>), ReadError> {
    let noun = bytes.noun;
    let format = bytes.i32()?;
    let preset = match format & 0xffff {
        0 => Preset::Generic,
        1 => Preset::Sam,
        2 => Preset::Vcf,
        3 => Preset::Gaf,
        other => {
            return Err(ReadError::whole(format!(
                "{noun} is for a kind of file tabix does not write, {other}"
            )))
        }
    };
    let sequence = bytes.count("the sequence's column")?;
    let start = bytes.count("the start's column")?;
    let end = bytes.count("the end's column")?;
    let comment = bytes.count("the comment character")?;
    let comment = u8::try_from(comment).map(char::from).map_err(|_| {
        ReadError::whole(format!(
            "{noun} gives a comment character past a byte, {comment}"
        ))
    })?;
    let skip = bytes.count("lines to skip")?;
    let length = bytes.count("bytes of names")?;
    let names = bytes
        .take(length)?
        .split(|byte| *byte == 0)
        .filter(|name| !name.is_empty())
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();
    Ok((
        Columns {
            preset,
            zero_based: format & 0x10000 != 0,
            sequence,
            start,
            end: (end > 0).then_some(end),
            comment,
            skip,
        },
        names,
    ))
}

impl Index {
    /// Which of the three this is.
    pub fn kind(&self) -> Kind {
        self.kind
    }

    /// How many bases a leaf of the tree spans, as a power of two: 14 for BAI
    /// and TBI.
    pub fn min_shift(&self) -> u32 {
        self.min_shift
    }

    /// How many levels the tree has under its root: 5 for BAI and TBI.
    pub fn depth(&self) -> u32 {
        self.depth
    }

    /// The columns of a text file's rows, for TBI and for a CSI written for
    /// text; `None` for an index of a BAM or a BCF, whose records say where
    /// they are themselves.
    pub fn columns(&self) -> Option<Columns> {
        self.columns
    }

    /// The sequences, in the order the file first writes them, which is the
    /// order [`Index::chunks`] numbers them in. Empty for an index of a BAM
    /// or a BCF, whose header names them.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Where a sequence the index names is in [`Index::names`].
    pub fn reference(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|named| named == name)
    }

    /// How many sequences the index has bins for.
    pub fn references(&self) -> usize {
        self.references.len()
    }

    /// What the pseudo-bin says about a sequence, where the index has one.
    pub fn summary(&self, reference: usize) -> Option<Summary> {
        self.references.get(reference)?.summary
    }

    /// How many records the file holds with no sequence at all, where the
    /// index says.
    pub fn unplaced(&self) -> Option<u64> {
        self.unplaced
    }

    /// The stretches of the file that can hold a record over `[start, end)`
    /// of sequence `reference`, as pairs of virtual offsets, in order and
    /// merged where they touch. Empty for a sequence the index has no bins
    /// for.
    pub fn chunks(&self, reference: usize, start: u64, end: u64) -> Vec<(u64, u64)> {
        let Some(sequence) = self.references.get(reference) else {
            return Vec::new();
        };
        let floor = self.floor(sequence, start);
        let mut chunks: Vec<(u64, u64)> = levels(self.min_shift, self.depth, start, end)
            .into_iter()
            .flat_map(|range| sequence.bins.range(range))
            .flat_map(|(_, bin)| bin.chunks.iter().copied())
            .filter(|(_, stop)| *stop > floor)
            .map(|(begin, stop)| (begin.max(floor), stop))
            .collect();
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(chunks.len());
        for (begin, stop) in chunks {
            match merged.last_mut() {
                Some(last) if begin <= last.1 => last.1 = last.1.max(stop),
                _ => merged.push((begin, stop)),
            }
        }
        merged
    }

    /// The virtual offset no record over `start` lies before, or nought where
    /// the index cannot say.
    ///
    /// For BAI and TBI that is the linear index over `start`'s 16 KiB window.
    /// A CSI says where the first record over each bin is, and every bin at
    /// or to the left of `start` on a level, or over it on a higher one,
    /// starts no later than the first record over `start`, since the records
    /// are sorted by where they start. The latest of those, the nearest bin to
    /// the left on each level, is the floor.
    fn floor(&self, sequence: &Reference, start: u64) -> u64 {
        match self.kind {
            Kind::Bai | Kind::Tbi => sequence
                .linear
                .get((start >> self.min_shift) as usize)
                .copied()
                .unwrap_or(0),
            Kind::Csi => {
                if start >> (self.min_shift + 3 * self.depth) > 0 {
                    return 0;
                }
                (0..=self.depth)
                    .filter_map(|level| {
                        let shift = self.min_shift + 3 * (self.depth - level);
                        let first = (((1u64 << (3 * level)) - 1) / 7) as u32;
                        let over = first + (start >> shift) as u32;
                        sequence
                            .bins
                            .range(first..=over)
                            .next_back()
                            .map(|(_, bin)| bin.first)
                    })
                    .max()
                    .unwrap_or(0)
            }
        }
    }
}

/// Two bgzipped text files and the indexes htslib 1.24 wrote for them, for
/// the tests here and for the readers that go through an index.
///
/// Each file was cut into blocks where a test needs a block to end, by
/// bgzipping pieces of it and joining them, and indexed as it stands: the
/// header in a block of its own, so the first row starts the second block,
/// and one row cut across two blocks.
#[cfg(test)]
pub(crate) mod fixture {
    /// A TBI, already out of its BGZF wrapper, written here rather than by
    /// tabix, for what tabix does not write: one stretch for the first of
    /// `names`, `[from, to)` in virtual offsets, filed under the root bin that
    /// every window is in, and no pseudo-bin, which the format leaves to the
    /// writer. `format` is tabix's word for the kind of file, and the columns
    /// of the sequence, the start and the end, as `[0x10000, 1, 2, 3]` for
    /// `-p bed`.
    pub(crate) fn rooted(names: &[&str], format: [i32; 4], stretch: (u64, u64)) -> Vec<u8> {
        let mut out = b"TBI\x01".to_vec();
        let int = |out: &mut Vec<u8>, value: i32| out.extend_from_slice(&value.to_le_bytes());
        int(&mut out, names.len() as i32);
        // The comment character, `#`, and no lines skipped.
        for value in format.into_iter().chain([35, 0]) {
            int(&mut out, value);
        }
        let joined: Vec<u8> = names
            .iter()
            .flat_map(|name| name.bytes().chain([0]))
            .collect();
        int(&mut out, joined.len() as i32);
        out.extend(joined);
        for at in 0..names.len() {
            int(&mut out, i32::from(at == 0));
            if at == 0 {
                out.extend_from_slice(&0u32.to_le_bytes());
                int(&mut out, 1);
                out.extend_from_slice(&stretch.0.to_le_bytes());
                out.extend_from_slice(&stretch.1.to_le_bytes());
            }
            // No linear index.
            int(&mut out, 0);
        }
        out
    }

    /// Rows over three sequences as BED: one over the edge of the first leaf,
    /// one of 4.88 Mb that only a bin near the root holds, and one at the start
    /// of the leaf at 2^20.
    pub(crate) const ROWS_TEXT: &str = "\
#rows over three sequences
chr1\t0\t10\ta
chr1\t5\t20000\tb
chr1\t16383\t16385\tc
chr1\t100000\t100001\td
chr1\t120000\t5000000\te
chr1\t600000\t600100\tf
chr1\t1048576\t1048577\tg
chr1\t3000000\t3000500\th
chr2\t10\t20\ti
chr2\t70000\t70010\tj
chr3\t1\t2\tk
";

    /// `ROWS_TEXT` in four blocks, cut after the comment, inside row d and
    /// before row h.
    pub(crate) const ROWS: [u8; 339] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x39, 0x00, 0x01, 0x1b, 0x00, 0xe4, 0xff, 0x23, 0x72, 0x6f, 0x77, 0x73, 0x20, 0x6f,
        0x76, 0x65, 0x72, 0x20, 0x74, 0x68, 0x72, 0x65, 0x65, 0x20, 0x73, 0x65, 0x71, 0x75, 0x65,
        0x6e, 0x63, 0x65, 0x73, 0x0a, 0x0d, 0x76, 0xd1, 0x35, 0x1b, 0x00, 0x00, 0x00, 0x1f, 0x8b,
        0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x3d,
        0x00, 0x4b, 0xce, 0x28, 0x32, 0xe4, 0x34, 0xe0, 0x34, 0x34, 0xe0, 0x4c, 0xe4, 0x4a, 0x06,
        0xb1, 0x4d, 0x39, 0x8d, 0x0c, 0x80, 0x80, 0x33, 0x09, 0xc2, 0x35, 0x34, 0x33, 0xb6, 0x30,
        0x06, 0x93, 0xa6, 0x9c, 0xc9, 0x5c, 0x00, 0xee, 0xb6, 0x25, 0x9f, 0x2e, 0x00, 0x00, 0x00,
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x27, 0x00, 0x01, 0x09, 0x00, 0xf6, 0xff, 0x63, 0x68, 0x72, 0x31, 0x09, 0x31, 0x30,
        0x30, 0x30, 0xb1, 0x19, 0x4b, 0xac, 0x09, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00,
        0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x4c, 0x00, 0x33, 0x30,
        0xe0, 0x34, 0x34, 0x00, 0x02, 0x43, 0xce, 0x14, 0xae, 0xe4, 0x8c, 0x22, 0x43, 0x4e, 0x43,
        0x23, 0x10, 0x97, 0xd3, 0xd4, 0x00, 0x0c, 0x38, 0x53, 0x21, 0xa2, 0x66, 0x10, 0x1e, 0x90,
        0x32, 0x04, 0x52, 0x69, 0x50, 0xa5, 0x06, 0x26, 0x16, 0xa6, 0xe6, 0x66, 0x50, 0xda, 0x9c,
        0x33, 0x9d, 0x0b, 0x00, 0xce, 0x01, 0x7c, 0xd0, 0x4e, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08,
        0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x49, 0x00,
        0x4b, 0xce, 0x28, 0x32, 0xe4, 0x34, 0x36, 0x00, 0x03, 0x30, 0x6d, 0x0a, 0xa4, 0x33, 0xb8,
        0x92, 0x33, 0x8a, 0x8c, 0x38, 0x0d, 0x0d, 0x38, 0x8d, 0x0c, 0x38, 0x33, 0x21, 0x1c, 0x73,
        0xb0, 0x0a, 0x20, 0x09, 0x14, 0xcd, 0x02, 0x09, 0x19, 0x73, 0x1a, 0x72, 0x1a, 0x71, 0x66,
        0x73, 0x01, 0x00, 0x74, 0x64, 0x3f, 0xe0, 0x42, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04,
        0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `tabix -p bed rows.bed.gz`
    pub(crate) const ROWS_TBI: [u8; 237] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0xd0, 0x00, 0xed, 0x96, 0x31, 0x0a, 0xc2, 0x40, 0x10, 0x45, 0xdf, 0x26, 0x28, 0x0a,
        0x82, 0xa0, 0x85, 0xbd, 0xe6, 0x02, 0xc6, 0x22, 0x60, 0x69, 0x17, 0xc1, 0xce, 0x1b, 0xa4,
        0x49, 0x63, 0x93, 0xca, 0x2b, 0xd8, 0x7a, 0x02, 0x2f, 0xe8, 0x1d, 0x84, 0x9d, 0x19, 0x4c,
        0x16, 0x05, 0xc1, 0x76, 0xa6, 0xd8, 0x07, 0x7f, 0xe6, 0xcf, 0x4e, 0xf9, 0xcf, 0x87, 0x3a,
        0xe4, 0x00, 0x04, 0x02, 0x90, 0x01, 0x39, 0xb0, 0x41, 0x6a, 0x0e, 0x34, 0x6d, 0xb7, 0xa5,
        0x69, 0xbb, 0x92, 0xa6, 0xed, 0x76, 0x8c, 0x01, 0x96, 0xc4, 0x71, 0x78, 0xc6, 0x77, 0xa5,
        0xac, 0x33, 0xd3, 0xf7, 0xf1, 0x9d, 0xf1, 0x88, 0x3c, 0xde, 0x65, 0xb5, 0xe9, 0x36, 0x3f,
        0x61, 0x58, 0xb7, 0x85, 0xf8, 0x2b, 0xf5, 0xd9, 0xfe, 0x8b, 0xea, 0x6b, 0xd5, 0xad, 0x3f,
        0x45, 0x74, 0xfb, 0xc7, 0xfa, 0x65, 0xe8, 0xdf, 0x20, 0x84, 0xeb, 0x4f, 0xb4, 0x5d, 0x4e,
        0xa7, 0xd3, 0xe9, 0x74, 0x3a, 0x9d, 0xce, 0xff, 0x98, 0x03, 0xb5, 0xe6, 0x38, 0xcb, 0x7f,
        0x85, 0xd2, 0xf2, 0xa1, 0xe9, 0x95, 0x52, 0x32, 0xe3, 0xbb, 0x4e, 0xea, 0x2f, 0x92, 0xb9,
        0xd1, 0x87, 0x9d, 0xdf, 0x98, 0xf5, 0xee, 0x30, 0x3f, 0x54, 0xa1, 0x7f, 0x47, 0xaa, 0x87,
        0xe4, 0x8e, 0xa1, 0x57, 0xea, 0x05, 0x07, 0xa6, 0x43, 0x1a, 0x4b, 0x0b, 0x00, 0x00, 0x1f,
        0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `tabix -C -p bed rows.bed.gz`, which chose a depth of 8.
    pub(crate) const ROWS_CSI: [u8; 217] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0xbc, 0x00, 0x73, 0x0e, 0xf6, 0x64, 0xe4, 0x63, 0x60, 0x60, 0xe0, 0x00, 0x62, 0x6d,
        0x06, 0x10, 0x60, 0x04, 0x42, 0x06, 0x06, 0x26, 0x20, 0x66, 0x06, 0x62, 0x65, 0x06, 0x08,
        0xe0, 0x07, 0xe2, 0xe4, 0x8c, 0x22, 0x43, 0x10, 0x61, 0x04, 0x22, 0x8c, 0xc1, 0xd2, 0x6c,
        0x20, 0xb9, 0xc9, 0x2a, 0x0c, 0x3c, 0x0c, 0x0b, 0xa0, 0x9a, 0x41, 0xe0, 0x2d, 0x98, 0x14,
        0x87, 0xd2, 0x9e, 0x93, 0x58, 0x80, 0xa4, 0x15, 0x92, 0x3c, 0x84, 0x0d, 0xd3, 0xe3, 0x35,
        0x49, 0x05, 0x22, 0x0c, 0xb5, 0x16, 0x26, 0x0f, 0xd3, 0xcf, 0xc1, 0x80, 0x0a, 0x3a, 0x27,
        0xa1, 0xda, 0x67, 0x0e, 0x65, 0xc3, 0xec, 0xcd, 0x45, 0x93, 0x57, 0x82, 0xb2, 0x61, 0xea,
        0x3c, 0x85, 0x18, 0x50, 0xdc, 0x03, 0x53, 0x0b, 0x53, 0xc7, 0x0c, 0x76, 0xb3, 0x0a, 0xdc,
        0x7e, 0x46, 0x24, 0xb7, 0xa8, 0x40, 0x69, 0x74, 0x37, 0xc3, 0xe4, 0xcd, 0xa1, 0x34, 0x13,
        0x9a, 0x9b, 0x7d, 0x81, 0xe6, 0xa9, 0x20, 0x99, 0xa7, 0x82, 0x45, 0x3d, 0xc8, 0x4e, 0x73,
        0x24, 0x35, 0x30, 0x36, 0x90, 0xc5, 0x88, 0xcd, 0x4e, 0x74, 0x79, 0x46, 0x06, 0xec, 0x00,
        0x00, 0x05, 0x1f, 0xa1, 0x8e, 0xe3, 0x01, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
        0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// Calls of two samples over two sequences, one a deletion of twenty bases
    /// across the edge of the first leaf.
    pub(crate) const CALLS_TEXT: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=2000000>
##contig=<ID=chr2,length=90000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
chr1\t10\tv1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
chr1\t16380\tv2\tACGTACGTACGTACGTACGT\tA\t.\tPASS\t.\tGT\t0/1\t0/0
chr1\t16390\tv3\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1
chr1\t50000\tv4\tG\tGA\t.\tPASS\t.\tGT\t1/1\t0/1
chr1\t1048577\tv5\tT\tC\t.\tPASS\t.\tGT\t0/1\t0/1
chr1\t1500000\tv6\tA\tC,G\t.\tPASS\t.\tGT\t1/2\t0/1
chr2\t5\tv7\tC\tA\t.\tPASS\t.\tGT\t0/1\t./.
chr2\t80000\tv8\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0
";

    /// `CALLS_TEXT` in four blocks, cut after the header, inside row v3 and
    /// before row v6.
    pub(crate) const CALLS: [u8; 532] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0xb3, 0x00, 0x65, 0x8d, 0x41, 0x0b, 0x82, 0x30, 0x18, 0x86, 0xcf, 0x5f, 0x3f, 0x43,
        0xaf, 0x12, 0x29, 0x5d, 0x82, 0x16, 0x2c, 0x75, 0x36, 0x50, 0x57, 0x73, 0x75, 0x2f, 0x99,
        0x3a, 0xd0, 0x4d, 0xd6, 0x0a, 0xfa, 0xf7, 0x89, 0xe1, 0xa9, 0xf7, 0xfa, 0x3c, 0x2f, 0x8f,
        0xef, 0x37, 0xaa, 0x97, 0x8d, 0xb1, 0xc3, 0xdd, 0xa1, 0x5b, 0x4c, 0xde, 0xdb, 0x75, 0xb4,
        0xf2, 0xfd, 0xda, 0x68, 0xa7, 0x5a, 0xb4, 0xa7, 0x09, 0xaa, 0x3b, 0x1b, 0x06, 0xbd, 0xd4,
        0xad, 0xeb, 0x50, 0xb4, 0x99, 0x77, 0xf8, 0x33, 0xa2, 0xc5, 0xd8, 0x2d, 0x9c, 0x30, 0x5e,
        0x60, 0x31, 0xf3, 0x4c, 0x04, 0xe5, 0x6b, 0x78, 0x48, 0x8b, 0xc2, 0x40, 0x7c, 0x46, 0x89,
        0x2a, 0x67, 0x95, 0x6e, 0x83, 0x44, 0x3e, 0x6b, 0xab, 0x46, 0xa7, 0x8c, 0x46, 0x5e, 0x26,
        0xb5, 0x71, 0x13, 0xf4, 0xa6, 0x73, 0x7c, 0xe2, 0xac, 0x80, 0x33, 0xab, 0x80, 0x26, 0xc0,
        0x53, 0x02, 0x38, 0x17, 0x70, 0xb9, 0xe2, 0x1c, 0x08, 0xcd, 0x45, 0xca, 0x81, 0x96, 0x84,
        0xc1, 0xaf, 0x00, 0x18, 0x8e, 0xab, 0x2f, 0xd5, 0x9a, 0xa9, 0xa6, 0xc6, 0x00, 0x00, 0x00,
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x53, 0x00, 0x4b, 0xce, 0x28, 0x32, 0xe4, 0x34, 0x34, 0xe0, 0x2c, 0x33, 0xe4, 0x74,
        0xe4, 0x74, 0xe7, 0xd4, 0xe3, 0x0c, 0x70, 0x0c, 0x0e, 0x06, 0x52, 0xee, 0x21, 0x9c, 0x06,
        0xfa, 0x40, 0x19, 0x7d, 0x43, 0xae, 0x64, 0xb0, 0x12, 0x33, 0x63, 0x0b, 0xa0, 0x2a, 0x23,
        0x4e, 0x47, 0x67, 0xf7, 0x10, 0x74, 0x0c, 0xd4, 0x8a, 0xae, 0xd1, 0x40, 0xdf, 0x80, 0x0b,
        0x00, 0x0a, 0xa2, 0xf3, 0xae, 0x5c, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
        0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x26, 0x00, 0x01, 0x08, 0x00,
        0xf7, 0xff, 0x63, 0x68, 0x72, 0x31, 0x09, 0x31, 0x36, 0x33, 0x01, 0x77, 0xd9, 0x37, 0x08,
        0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
        0x42, 0x43, 0x02, 0x00, 0x5f, 0x00, 0xb3, 0x34, 0xe0, 0x2c, 0x33, 0xe6, 0x74, 0xe6, 0x0c,
        0xe1, 0xd4, 0xe3, 0x0c, 0x70, 0x0c, 0x0e, 0x06, 0x52, 0xee, 0x21, 0x9c, 0x06, 0xfa, 0x06,
        0x40, 0x6c, 0xc8, 0x95, 0x9c, 0x51, 0x64, 0xc8, 0x69, 0x6a, 0x00, 0x04, 0x9c, 0x65, 0x26,
        0x9c, 0xee, 0x9c, 0xee, 0x8e, 0x28, 0xca, 0x0c, 0xf5, 0x0d, 0x11, 0xca, 0x0c, 0x0d, 0x4c,
        0x2c, 0x4c, 0xcd, 0xcd, 0x39, 0xcb, 0x4c, 0x81, 0x86, 0x39, 0xa3, 0x19, 0x07, 0x51, 0x07,
        0x00, 0x76, 0x04, 0x8f, 0xb5, 0x6d, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
        0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x68, 0x00, 0x4b, 0xce, 0x28,
        0x32, 0xe4, 0x34, 0x34, 0x35, 0x00, 0x01, 0xce, 0x32, 0x33, 0x4e, 0x47, 0x4e, 0x67, 0x1d,
        0x77, 0x4e, 0x3d, 0xce, 0x00, 0xc7, 0xe0, 0x60, 0x20, 0xe5, 0x1e, 0xc2, 0x69, 0xa8, 0x6f,
        0xc4, 0x69, 0xa0, 0x6f, 0xc8, 0x95, 0x9c, 0x51, 0x64, 0xc4, 0x69, 0xca, 0x59, 0x66, 0xce,
        0xe9, 0x0c, 0x54, 0x85, 0xac, 0x02, 0x28, 0xcb, 0xa9, 0xa7, 0xaf, 0x07, 0x51, 0x61, 0x01,
        0x31, 0xc9, 0x82, 0xd3, 0x9d, 0x33, 0x04, 0x55, 0x55, 0x0d, 0xd0, 0xa2, 0x1a, 0x03, 0x2e,
        0x00, 0xff, 0xfc, 0xf4, 0xa6, 0x72, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
        0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `tabix -p vcf calls.vcf.gz`
    pub(crate) const CALLS_TBI: [u8; 195] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0xa6, 0x00, 0x0b, 0x71, 0xf2, 0x64, 0x64, 0x62, 0x60, 0x60, 0x60, 0x62, 0x60, 0x60,
        0x60, 0x84, 0xd2, 0x0c, 0x0c, 0x0c, 0x0c, 0xca, 0x50, 0x9a, 0x8b, 0x81, 0x81, 0x21, 0x39,
        0xa3, 0xc8, 0x90, 0x21, 0x39, 0xa3, 0xc8, 0x88, 0x81, 0x85, 0x81, 0x81, 0x61, 0x89, 0x10,
        0x44, 0x25, 0x03, 0x43, 0x3f, 0x98, 0xd2, 0x82, 0xd2, 0x9e, 0x4c, 0x30, 0xf1, 0x2d, 0x60,
        0xd2, 0x95, 0x41, 0x1f, 0xcc, 0xf5, 0x9a, 0x04, 0x33, 0x15, 0x22, 0x0e, 0x53, 0xcf, 0xc6,
        0x80, 0x0a, 0x3a, 0xa1, 0xe6, 0xc2, 0xf4, 0xc1, 0xcc, 0x8f, 0x41, 0xd2, 0xab, 0x0c, 0xa5,
        0xe5, 0xa0, 0x6a, 0x60, 0xb4, 0xeb, 0x28, 0x3d, 0x22, 0x68, 0x58, 0x9a, 0x18, 0x6e, 0x34,
        0x33, 0x03, 0x03, 0x83, 0x27, 0x34, 0xfd, 0xc3, 0xf2, 0x87, 0x0f, 0x94, 0x86, 0xe5, 0x1f,
        0x2d, 0xb8, 0x9e, 0x1f, 0x60, 0x1a, 0x96, 0x53, 0x61, 0xc0, 0x17, 0xaa, 0xdf, 0x07, 0x4d,
        0x1d, 0x2b, 0x16, 0x33, 0x71, 0xd1, 0x30, 0x00, 0x00, 0x4c, 0xd3, 0x70, 0x74, 0x16, 0x04,
        0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42,
        0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `bcftools index calls.vcf.gz`, which chose a depth of 6.
    pub(crate) const CALLS_CSI: [u8; 179] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x96, 0x00, 0x73, 0x0e, 0xf6, 0x64, 0xe4, 0x63, 0x60, 0x60, 0x60, 0x03, 0x62, 0x35,
        0x20, 0x66, 0x02, 0x62, 0x46, 0x28, 0x0d, 0x02, 0xca, 0x50, 0x9a, 0x0b, 0x88, 0x93, 0x33,
        0x8a, 0x0c, 0x41, 0x84, 0x11, 0x58, 0x96, 0x05, 0x88, 0x97, 0x4c, 0x02, 0xc9, 0xf5, 0x33,
        0x32, 0x40, 0x75, 0xc1, 0xd8, 0x5a, 0x50, 0xda, 0x53, 0x08, 0x44, 0x6e, 0x61, 0x40, 0xc8,
        0x43, 0xd8, 0xae, 0x0c, 0xfa, 0x60, 0xae, 0xd7, 0x24, 0x16, 0x06, 0x18, 0x60, 0x42, 0x92,
        0x87, 0xe9, 0x67, 0x63, 0x40, 0x05, 0x9d, 0x93, 0x10, 0x7a, 0x19, 0x91, 0xcc, 0x81, 0xd9,
        0xcb, 0x0c, 0xb2, 0x73, 0x12, 0x42, 0x3f, 0x23, 0x92, 0x59, 0x3e, 0x50, 0x1a, 0xdd, 0x4e,
        0x98, 0x3c, 0x03, 0xc3, 0x0f, 0x46, 0x06, 0x24, 0x9f, 0xc3, 0x80, 0xef, 0x24, 0x84, 0x5e,
        0x46, 0x24, 0x73, 0x60, 0xea, 0x61, 0x00, 0x00, 0xd1, 0xea, 0x48, 0xb1, 0x4a, 0x01, 0x00,
        0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
        0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// Two rows: one that runs past the first 128 kb, which only a bin of a
    /// megabase holds, and one in a leaf of its own after it, as `bgzip` wrote
    /// them in one block.
    pub(crate) const FLOOR: [u8; 83] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x36, 0x00, 0x4b, 0xce, 0x28, 0x32, 0xe4, 0x34, 0x34, 0xe0, 0x34, 0x34, 0x31, 0x00,
        0x02, 0xce, 0x02, 0xae, 0x64, 0x90, 0x80, 0x91, 0x01, 0x98, 0x07, 0xa2, 0x80, 0x72, 0x85,
        0x5c, 0x00, 0x30, 0xc9, 0xc1, 0xc9, 0x26, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00,
        0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `tabix -p bed floor.bed.gz`
    pub(crate) const FLOOR_TBI: [u8; 119] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x5a, 0x00, 0x0b, 0x71, 0xf2, 0x64, 0x64, 0x64, 0x00, 0x01, 0x46, 0x20, 0x64, 0x60,
        0x60, 0x02, 0x62, 0x66, 0x20, 0x56, 0x66, 0x80, 0x00, 0x56, 0x20, 0x4e, 0xce, 0x28, 0x32,
        0x04, 0x0b, 0x7a, 0x82, 0x95, 0x21, 0x80, 0x20, 0x94, 0xf6, 0x9a, 0x04, 0xd1, 0x88, 0x00,
        0xe6, 0x60, 0x12, 0x55, 0x8c, 0x81, 0x21, 0x54, 0x08, 0xa2, 0x5f, 0x10, 0x4d, 0x1d, 0x2f,
        0x03, 0xf5, 0x80, 0x20, 0x01, 0x1a, 0x06, 0x00, 0xb7, 0x7a, 0xef, 0x2d, 0xf9, 0x00, 0x00,
        0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
        0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// `tabix -C -p bed floor.bed.gz`
    pub(crate) const FLOOR_CSI: [u8; 129] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x64, 0x00, 0x73, 0x0e, 0xf6, 0x64, 0xe4, 0x63, 0x60, 0x60, 0xe0, 0x00, 0x62, 0x45,
        0x06, 0x10, 0x60, 0x04, 0x42, 0x06, 0x06, 0x26, 0x20, 0x66, 0x06, 0x62, 0x65, 0x06, 0x08,
        0x60, 0x05, 0xe2, 0xe4, 0x8c, 0x22, 0x43, 0xb0, 0x24, 0x48, 0xc2, 0x73, 0x12, 0x03, 0x1c,
        0x30, 0x32, 0x20, 0x80, 0x20, 0x94, 0xf6, 0x9a, 0xa4, 0x02, 0x17, 0x66, 0x62, 0x40, 0x06,
        0xe6, 0x58, 0xc4, 0x18, 0x18, 0x42, 0x27, 0xa9, 0xc0, 0xf5, 0x32, 0x22, 0x99, 0x03, 0x53,
        0x0f, 0x03, 0x00, 0x38, 0xc9, 0x4c, 0xf3, 0xb1, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04,
        0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
}

#[cfg(test)]
mod tests {
    use super::fixture::*;
    use super::*;
    use crate::read::bgzf::Bgzf;
    use std::io::Cursor;

    fn parsed(bytes: &[u8]) -> Index {
        parse(bytes).unwrap_or_else(|error| panic!("{error}"))
    }

    /// Where a row is, 0-based and half-open, as tabix places it with the
    /// columns its index keeps, and what the row is called.
    fn placed(line: &str, columns: Columns) -> (String, u64, u64, String) {
        let fields: Vec<&str> = line.split('\t').collect();
        let start: u64 = fields[columns.start - 1].parse().unwrap();
        let start = if columns.zero_based { start } else { start - 1 };
        let end = match (columns.end, columns.preset) {
            (Some(column), _) => fields[column - 1].parse().unwrap(),
            (None, Preset::Vcf) => start + fields[3].len() as u64,
            (None, _) => start + 1,
        };
        let name = if columns.preset == Preset::Vcf {
            fields[2]
        } else {
            fields[3]
        };
        (
            fields[columns.sequence - 1].to_string(),
            start,
            end,
            name.to_string(),
        )
    }

    /// The rows over a window, read through the index's chunks.
    fn through(index: &Index, data: &[u8], sequence: &str, start: u64, end: u64) -> Vec<String> {
        let columns = index.columns().expect("an index of text");
        let Some(reference) = index.reference(sequence) else {
            return Vec::new();
        };
        let mut bgzf = Bgzf::new(Cursor::new(data));
        let mut line = Vec::new();
        let mut found = Vec::new();
        for (begin, stop) in index.chunks(reference, start, end) {
            bgzf.seek(begin).unwrap();
            while bgzf.tell() < stop && bgzf.line(&mut line).unwrap() {
                let line = String::from_utf8(line.clone()).unwrap();
                if line.starts_with(columns.comment) {
                    continue;
                }
                let (on, from, to, name) = placed(&line, columns);
                if on == sequence && from < end && to > start {
                    found.push(name);
                }
            }
        }
        found
    }

    /// The rows over a window, found by reading every row.
    fn scanned(text: &str, columns: Columns, sequence: &str, start: u64, end: u64) -> Vec<String> {
        text.lines()
            .filter(|line| !line.starts_with(columns.comment))
            .map(|line| placed(line, columns))
            .filter(|(on, from, to, _)| on == sequence && *from < end && *to > start)
            .map(|(_, _, _, name)| name)
            .collect()
    }

    /// The four text indexes, each with the file and the text it indexes.
    fn indexed() -> [(&'static str, Index, &'static [u8], &'static str); 4] {
        [
            ("rows, tbi", parsed(&ROWS_TBI), &ROWS[..], ROWS_TEXT),
            ("rows, csi", parsed(&ROWS_CSI), &ROWS[..], ROWS_TEXT),
            ("calls, tbi", parsed(&CALLS_TBI), &CALLS[..], CALLS_TEXT),
            ("calls, csi", parsed(&CALLS_CSI), &CALLS[..], CALLS_TEXT),
        ]
    }

    #[test]
    fn a_tbi_says_its_columns_its_sequences_and_their_rows() {
        let calls = parsed(&CALLS_TBI);
        assert_eq!(calls.kind(), Kind::Tbi);
        assert_eq!((calls.min_shift(), calls.depth()), (14, 5));
        assert_eq!(
            calls.columns(),
            Some(Columns {
                preset: Preset::Vcf,
                zero_based: false,
                sequence: 1,
                start: 2,
                end: None,
                comment: '#',
                skip: 0,
            })
        );
        assert_eq!(calls.names(), ["chr1", "chr2"]);
        assert_eq!(calls.reference("chr2"), Some(1));
        assert_eq!(calls.reference("chr3"), None);
        let rows: Vec<Option<u64>> = (0..calls.references())
            .map(|at| calls.summary(at).map(|summary| summary.placed))
            .collect();
        assert_eq!(rows, [Some(6), Some(2)]);
        assert_eq!(calls.summary(0).map(|summary| summary.unplaced), Some(0));

        let rows = parsed(&ROWS_TBI);
        let columns = rows.columns().unwrap();
        assert!(columns.zero_based && columns.preset == Preset::Generic);
        assert_eq!(
            (columns.sequence, columns.start, columns.end),
            (1, 2, Some(3))
        );
        assert_eq!(rows.names(), ["chr1", "chr2", "chr3"]);
    }

    /// A CSI written for text carries the columns and names a TBI does, and
    /// a tree of the depth its writer chose, which has to be read.
    #[test]
    fn a_csi_carries_what_its_tbi_carries_with_a_tree_of_its_own() {
        for (tbi, csi, depth) in [
            (&ROWS_TBI[..], &ROWS_CSI[..], 8),
            (&CALLS_TBI[..], &CALLS_CSI[..], 6),
        ] {
            let (tbi, csi) = (parsed(tbi), parsed(csi));
            assert_eq!(csi.kind(), Kind::Csi);
            assert_eq!((csi.min_shift(), csi.depth()), (14, depth));
            assert_eq!(csi.columns(), tbi.columns());
            assert_eq!(csi.names(), tbi.names());
            for at in 0..tbi.references() {
                assert_eq!(csi.summary(at), tbi.summary(at), "sequence {at}");
            }
        }
        // samtools wrote a tree of one level for a BAM of 1,000 bases, and no
        // columns, since a BAM's records say where they are.
        let bam = parsed(&crate::read::bam::fixture::CSI);
        assert_eq!((bam.kind(), bam.depth()), (Kind::Csi, 0));
        assert!(bam.columns().is_none() && bam.names().is_empty());
        assert_eq!(bam.references(), 2);
    }

    /// The rows over each window are the rows tabix printed for it, through
    /// the TBI and through the CSI. The regions are written as tabix was
    /// given them, 1-based and inclusive.
    #[test]
    fn a_window_holds_the_rows_tabix_prints_for_it() {
        let rows = [
            ("chr1:1-1", "a"),
            ("chr1:11-16384", "b c"),
            ("chr1:16385-16385", "b c"),
            ("chr1:20001-100000", ""),
            ("chr1:100001-100001", "d"),
            ("chr1:4000000-4000001", "e"),
            ("chr1:5000001-6000000", ""),
            ("chr1:1048577-1048577", "e g"),
            ("chr2:1-70000", "i"),
            ("chr2:70011-80000", ""),
            ("chr3:1-2", "k"),
            ("chr4:1-10", ""),
        ];
        let calls = [
            ("chr1:10-10", "v1"),
            ("chr1:11-16379", ""),
            ("chr1:16399-16399", "v2"),
            ("chr1:16400-16400", ""),
            ("chr1:16390-50000", "v2 v3 v4"),
            ("chr1:50001-1048576", ""),
            ("chr1:1048577-2000000", "v5 v6"),
            ("chr2:1-5", "v7"),
            ("chr2:6-79999", ""),
            ("chr2:80000-90000", "v8"),
            ("chr3:1-10", ""),
        ];
        for (what, index, data, _) in indexed() {
            let printed = if what.starts_with("rows") {
                &rows[..]
            } else {
                &calls[..]
            };
            for (locus, names) in printed {
                let region = crate::Region::parse(locus).unwrap();
                let found = through(&index, data, region.seq(), region.start(), region.end());
                assert_eq!(found.join(" "), *names, "{what}, {locus}");
            }
        }
    }

    /// Every row over a window is in the chunks the index gives for it, over
    /// a grid of windows at the edges of the leaves and the levels above
    /// them. A floor set too high, a level left out or a bin misnumbered
    /// loses a row here.
    #[test]
    fn every_row_over_a_window_is_in_its_chunks() {
        let starts = [
            0u64, 1, 9, 10, 5_000, 16_379, 16_382, 16_383, 16_384, 16_385, 16_399, 20_000, 49_999,
            50_000, 99_999, 100_000, 100_001, 119_999, 120_000, 600_050, 1_048_575, 1_048_576,
            1_048_577, 1_499_999, 2_999_999, 3_000_499, 4_999_999, 5_000_000, 5_000_001, 8_000_000,
        ];
        let widths = [1u64, 2, 100, 16_384, 100_000, 3_000_000];
        let mut windows = 0;
        for (what, index, data, text) in indexed() {
            let columns = index.columns().unwrap();
            for sequence in ["chr1", "chr2", "chr3", "chr4"] {
                for start in starts {
                    for width in widths {
                        let end = start + width;
                        assert_eq!(
                            through(&index, data, sequence, start, end),
                            scanned(text, columns, sequence, start, end),
                            "{what}, {sequence}:{start}-{end}"
                        );
                        windows += 1;
                    }
                }
            }
        }
        assert_eq!(windows, 4 * 4 * 30 * 6);
    }

    /// The pseudo-bin says where the first row is, which is where the header
    /// ends. Each fixture's header fills a block of its own, so the first row
    /// starts the next block, and a reader that says the end of the header's
    /// block instead disagrees with the index on a sound file.
    #[test]
    fn the_header_ends_where_the_index_says_the_rows_begin() {
        for (what, index, data, _) in indexed() {
            let columns = index.columns().unwrap();
            let mut bgzf = Bgzf::new(Cursor::new(data));
            let mut line = Vec::new();
            let mut first = bgzf.tell();
            while bgzf.line(&mut line).unwrap() && line.starts_with(&[columns.comment as u8]) {
                first = bgzf.tell();
            }
            let summary = index.summary(0).unwrap();
            assert_eq!(first, summary.first, "{what}");
            assert_eq!(first & 0xffff, 0, "{what}: the rows start a block");
        }
        // And the file the documentation draws from, whose index puts its
        // first row at the start of its second block, 210 bytes in.
        let docs = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/data");
        let index = parsed(&std::fs::read(docs.join("calls.vcf.gz.tbi")).unwrap());
        assert_eq!(
            index.summary(0).map(|summary| summary.first),
            Some(210 << 16)
        );
        assert_eq!(index.names(), ["NC_000962.3"]);
    }

    /// The bins on each level are the ones the SAM specification's
    /// `reg2bins` lists, written out for BAI's tree as samtools has it, and
    /// the pseudo-bin is the first bin under the deepest level, at any depth.
    #[test]
    fn the_bins_of_a_window_are_the_specification_s() {
        fn reg2bins(start: u64, end: u64) -> Vec<u32> {
            let last = end.saturating_sub(1).max(start);
            let mut bins = vec![0u32];
            for (shift, first) in [(26u32, 1u64), (23, 9), (20, 73), (17, 585), (14, 4681)] {
                for bin in first + (start >> shift)..=first + (last >> shift) {
                    bins.push(bin as u32);
                }
            }
            bins
        }
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        for _ in 0..2_000 {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            let start = seed % (1 << 29);
            let end = (start + 1 + (seed >> 40) % (1 << 22)).min(1 << 29);
            let listed: Vec<u32> = levels(14, 5, start, end).into_iter().flatten().collect();
            assert_eq!(listed, reg2bins(start, end), "{start}-{end}");
        }
        // The edge of the first leaf, counted from nought.
        assert_eq!(levels(14, 5, 16_383, 16_384)[5], 4681..=4681);
        assert_eq!(levels(14, 5, 16_384, 16_385)[5], 4682..=4682);
        let pseudo: Vec<u32> = [5, 6, 8, 3, 0].into_iter().map(pseudo_bin).collect();
        assert_eq!(pseudo, [37_450, 299_594, 19_173_962, 586, 2]);
        // Past the last base a tree has room for, nothing could be filed.
        assert!(levels(14, 5, 1 << 29, (1 << 29) + 10).is_empty());
        assert!(levels(14, 0, 16_384, 20_000).is_empty());
    }

    /// The virtual offset each row of a bgzipped file starts at, by its name
    /// in the fourth column.
    fn starts(data: &[u8]) -> Vec<(String, u64)> {
        let mut bgzf = Bgzf::new(Cursor::new(data));
        let mut line = Vec::new();
        let mut at = bgzf.tell();
        let mut found = Vec::new();
        while bgzf.line(&mut line).unwrap() {
            let text = String::from_utf8(line.clone()).unwrap();
            if let Some(name) = text.split('\t').nth(3) {
                found.push((name.to_string(), at));
            }
            at = bgzf.tell();
        }
        found
    }

    /// A window is read from the first row that can reach it, not from the
    /// start of every stretch a bin over it holds. A bin near the root holds
    /// the long rows that start far to the left, and only some of them reach
    /// the window: read from the start of its stretch, a figure of one gene
    /// at the end of a chromosome would read every long row before it.
    #[test]
    fn a_window_starts_at_the_first_row_that_can_reach_it() {
        let at = |data: &[u8], name: &str| {
            starts(data)
                .into_iter()
                .find(|(named, _)| named == name)
                .map(|(_, offset)| offset)
                .unwrap()
        };
        // Row p runs from 10 to 140,000, so a bin of a megabase holds it, and
        // row q sits in a leaf at 200,000. The bins over q hold both.
        for (what, index) in [("tbi", parsed(&FLOOR_TBI)), ("csi", parsed(&FLOOR_CSI))] {
            let q = index.chunks(0, 200_000, 200_001);
            assert_eq!(
                q.first().map(|chunk| chunk.0),
                Some(at(&FLOOR, "q")),
                "{what}"
            );
            let p = index.chunks(0, 100_000, 100_001);
            assert_eq!(
                p.first().map(|chunk| chunk.0),
                Some(at(&FLOOR, "p")),
                "{what}"
            );
        }
        // htslib folds the leaves of rows a and d into the bin of 128 kb that
        // holds b, so one stretch holds all four; the linear index says d is
        // the first over 100,000, and a TBI reads from there. A CSI says only
        // where the first row over each bin is, which for that bin is a.
        let d = parsed(&ROWS_TBI).chunks(0, 100_000, 100_001);
        assert_eq!(d.first().map(|chunk| chunk.0), Some(at(&ROWS, "d")));
        let a = parsed(&ROWS_CSI).chunks(0, 100_000, 100_001);
        assert_eq!(a.first().map(|chunk| chunk.0), Some(at(&ROWS, "a")));
    }

    /// A window on a sequence the index has no bins for, or past its tree,
    /// is no chunks rather than a panic or the whole file.
    #[test]
    fn a_window_the_index_has_nothing_for_is_no_chunks() {
        let calls = parsed(&CALLS_TBI);
        assert!(calls.chunks(2, 0, 100).is_empty());
        assert!(calls.chunks(0, 1 << 29, (1 << 29) + 100).is_empty());
        assert!(!calls.chunks(0, 0, 100).is_empty());
        assert!(parsed(&ROWS_CSI)
            .chunks(0, 1 << 40, (1 << 40) + 1)
            .is_empty());
    }

    #[test]
    fn a_damaged_index_is_an_error_and_never_a_panic() {
        let all = [
            &ROWS_TBI[..],
            &ROWS_CSI[..],
            &CALLS_TBI[..],
            &CALLS_CSI[..],
            &FLOOR_TBI[..],
            &FLOOR_CSI[..],
            &crate::read::bam::fixture::BAI[..],
            &crate::read::bam::fixture::CSI[..],
        ];
        for bytes in all {
            // Out of their wrapper too, where a cut or a changed byte reaches
            // the counts rather than the CRC.
            let inflated = if gzip::is_gzip(bytes) {
                gzip::decompress(bytes).unwrap()
            } else {
                bytes.to_vec()
            };
            for data in [bytes.to_vec(), inflated] {
                for cut in 0..data.len() {
                    if let Ok(index) = parse(&data[..cut]) {
                        let _ = index.chunks(0, 0, 1 << 20);
                    }
                }
                for at in 0..data.len() {
                    for flip in [0x01, 0x80, 0xff] {
                        let mut bent = data.clone();
                        bent[at] ^= flip;
                        if let Ok(index) = parse(&bent) {
                            for reference in 0..index.references() {
                                let _ = index.chunks(reference, 16_000, 1 << 21);
                            }
                        }
                    }
                }
            }
        }
        let error = parse(b"not an index").unwrap_err();
        assert!(error.to_string().contains("none of BAI's"), "{error}");
        // A count below nought, which a sound index never writes.
        let mut bent = gzip::decompress(&CALLS_TBI).unwrap();
        bent[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        let error = parse(&bent).unwrap_err();
        assert!(error.to_string().contains("impossible count"), "{error}");
        let error = parse(&gzip::decompress(&CALLS_TBI).unwrap()[..50]).unwrap_err();
        assert_eq!(error.to_string(), "the tabix index is cut short");
    }
}
