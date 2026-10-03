//! Juicer's `.hic`, a contact map, read a window at a time through the
//! indexes it holds.
//!
//! A `.hic` is the contact map of a Hi-C experiment at several resolutions:
//! for each two sequences and each size of bin, how many read pairs joined
//! each two bins, kept in blocks compressed with zlib that an index at the end
//! of the file points into. [`contacts`] reads the map of one sequence with
//! itself over a window, at one resolution, and hands its cells over as
//! [`Pair`]s, two bins and a count, which a [`PairTrack`](crate::PairTrack)
//! draws as a triangle under the axis.
//!
//! # What is read
//!
//! The header, for the sequences, their lengths and the resolutions in
//! bases; the master index at the end of the file, for where the map of the
//! window's sequence with itself is kept; that map's list of blocks at the
//! resolution asked for; and of those, the blocks that can hold a cell of the
//! window, and no other. A block keeps its cells one of two ways, as rows of
//! the cells that hold something or as a dense square, and writes each bin
//! and each count short or long as its writer chose; every one of those is
//! read, and only the cells with both bins over the window are kept.
//!
//! The counts are the ones observed, raw, as `hictk dump` prints them with no
//! `--balance`. The normalisations a file may carry beside them, KR, VC or
//! SCALE, are not read: each file carries its own set, each a vector per
//! sequence and resolution, and a map drawn from one is a different map.
//!
//! # Which bins hold the window
//!
//! Version 9 files the cells of a sequence with itself along the diagonal: a
//! block is a stretch of `blockBinCount` bins of the diagonal, by where half
//! way between a cell's two bins falls, in one of a series of bands out from
//! it, each twice as wide as the one before, by how far apart its two bins
//! are, `log2(1 + |x - y| / √2 / blockBinCount)`. A block's number is its
//! band times the file's `blockColumnCount`, plus its stretch. A window is a
//! triangle on the diagonal, so it reaches the stretches under it and the
//! bands out to its own width, and those are the blocks read. In the map these
//! tests read, its 1,235 bins at 1 kb are four blocks, two stretches in each
//! of two bands, and a window of 120 bins reads one of them.
//!
//! # Coordinates
//!
//! A cell is two bins, each numbered from nought along its sequence, so bin
//! `n` at a resolution of `r` bases is the bases from `n × r` to
//! `(n + 1) × r`, 0-based and half-open, and the last bin of a sequence stops
//! where the sequence does. A window holds every bin it touches, as `hictk
//! dump -r` takes one.
//!
//! # Which version
//!
//! Version 9, which hictk writes. Version 8 and older lay their blocks on a
//! grid rather than along the diagonal, write a sequence's length in 32 bits
//! and the bins of a block in one width, where version 9 says which of two,
//! and no file of them could be made here to test against: they are refused by
//! their version, with the commands that write the same map as version 9. A
//! newer version is refused the same way.
//!
//! # What is refused
//!
//! A file that is not a `.hic`, one of another version, and one damaged or
//! cut short: an offset or a length past the end of the file, a count below
//! nought, a block that does not inflate or inflates past [`BLOCK_MOST`], a
//! block of a kind the format does not have, and a cell before the start of
//! its sequence. A sequence the file does not have, and a resolution it does
//! not hold, are refused naming the ones it does.
//!
//! ```
//! use std::io::Cursor;
//! use karyon::{read, Region};
//!
//! # let bytes = include_bytes!("fixtures/hic/contacts.hic").to_vec();
//! let region = Region::parse("chr1:1-1,234,567")?;
//! let header = read::hic::header_of(Cursor::new(&bytes))?;
//! let resolution = read::hic::resolution_for(&header.resolutions, &region, 250);
//! assert_eq!(resolution, Some(5_000));
//! let cells = read::hic::contacts(Cursor::new(&bytes), &region, 50_000)?;
//! // The last bin stops where the sequence does.
//! assert_eq!(cells.last().unwrap().second, (1_200_000, 1_234_567));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::io::{Read, Seek};

use super::bytes::{short, Bytes, File};
use super::gzip::zlib;
use super::ReadError;
use crate::region::Region;
use crate::track::axis::group_thousands;
use crate::track::pairs::Pair;

/// The most bins a window is cut into where the resolution is chosen for it.
///
/// A contact map of a window is a triangle of cells, up to half the square of
/// its bins, and each is a polygon of about 94 bytes of SVG. Over a simulated
/// chromosome of 249 Mb, its 249 bins at 1 Mb drew 31,125 cells in a figure of
/// 2.9 MB, and a window of it cut into 2,000 bins drew 173,516 in 16.2 MB.
pub const BINS: u64 = 250;

/// The one version read.
const VERSION: i32 = 9;

/// The most a block may inflate to.
///
/// A file does not say how large its blocks are once inflated, and DEFLATE
/// writes a run of 258 bytes in as little as a byte, so a damaged block could
/// ask for far more memory than the file is long. hictk writes a block per
/// stretch of about a thousand bins, and the largest of a simulated chromosome
/// of 249 Mb, 7.6 million cells at 5 kb, zoomed out to eleven resolutions,
/// inflated to 4.5 MB: this is over fifty times that.
pub const BLOCK_MOST: usize = 256 << 20;

/// The longest name or attribute read, past which a string is taken for
/// damage rather than read on into the rest of the file.
const LONGEST: usize = 16 << 20;

/// How much of the file is read at a time where how much is wanted is not
/// known until it is read: the header and the master index.
const CHUNK: u64 = 64 << 10;

/// What a `.hic` says about itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// The version of the format, which is 9 for every file read.
    pub version: i32,
    /// The assembly the file names, as its writer was told it.
    pub genome: String,
    /// Each sequence with its length, in the order of the file, less the
    /// `All` a file of several resolutions opens with, which is every
    /// sequence end to end in kilobases and no sequence of the genome.
    pub sequences: Vec<(String, u64)>,
    /// The resolutions the file holds in bases, in the order it lists them.
    pub resolutions: Vec<u32>,
    /// Where the master index starts.
    footer: u64,
    /// The number the file gives each of `sequences`, which its maps are
    /// named by in the master index.
    numbers: Vec<usize>,
}

/// Reads the header of a `.hic`.
///
/// # Errors
///
/// A file that is not a `.hic`, of another version than 9, or damaged.
pub fn header_of<R: Read + Seek>(reader: R) -> Result<Header, ReadError> {
    let mut file = File::new(reader, ".hic")?;
    read_header(&mut file)
}

/// The sequences a `.hic` holds, each with its length, in the order of the
/// file.
///
/// # Errors
///
/// As [`header_of`].
pub fn sequences<R: Read + Seek>(reader: R) -> Result<Vec<(String, u64)>, ReadError> {
    header_of(reader).map(|header| header.sequences)
}

/// The cells of the map of `region`'s sequence with itself, at `resolution`
/// bases a bin, that have both bins over the window: each the lower bin
/// first, its count as observed, in order of their bins.
///
/// A sequence the file holds no contacts on has no map, and a window on it
/// none.
///
/// # Errors
///
/// A file that is not a `.hic`, of another version than 9, or damaged; a
/// sequence the file does not have and a resolution it does not hold, each
/// naming the ones it does.
pub fn contacts<R: Read + Seek>(
    reader: R,
    region: &Region,
    resolution: u32,
) -> Result<Vec<Pair>, ReadError> {
    let mut file = File::new(reader, ".hic")?;
    let header = read_header(&mut file)?;
    header.read_contacts(&mut file, region, resolution)
}

/// The finest of `resolutions` that cuts `region` into at most `most` bins,
/// counting each bin the window touches, or the coarsest where none does.
/// `None` where there is no resolution, as in a file of restriction fragments
/// alone.
pub fn resolution_for(resolutions: &[u32], region: &Region, most: u64) -> Option<u32> {
    let bins = |resolution: u32| {
        let size = u64::from(resolution);
        region.end().saturating_sub(1) / size - region.start() / size + 1
    };
    let known = resolutions.iter().copied().filter(|size| *size > 0);
    known
        .clone()
        .filter(|size| bins(*size) <= most)
        .min()
        .or_else(|| known.max())
}

/// Cells as BEDPE, on `sequence`: what the command line hands the reader of
/// pairs, so a `.hic` is drawn by the same reader as the BEDPE `hictk dump
/// --join` writes for it. Each count is written as the shortest decimal that
/// reads back as the same number.
pub fn bedpe(sequence: &str, cells: &[Pair]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(cells.len() * (2 * sequence.len() + 48));
    for cell in cells {
        // Writing to a string does not fail.
        let _ = writeln!(
            out,
            "{sequence}\t{}\t{}\t{sequence}\t{}\t{}\t{}",
            cell.first.0, cell.first.1, cell.second.0, cell.second.1, cell.value
        );
    }
    out
}

impl Header {
    /// [`contacts`], from a file whose header this is, read again from
    /// `reader`.
    ///
    /// # Errors
    ///
    /// As [`contacts`].
    pub fn contacts<R: Read + Seek>(
        &self,
        reader: R,
        region: &Region,
        resolution: u32,
    ) -> Result<Vec<Pair>, ReadError> {
        let mut file = File::new(reader, ".hic")?;
        self.read_contacts(&mut file, region, resolution)
    }

    /// [`contacts`], from `file`, measured and with its header read.
    fn read_contacts<R: Read + Seek>(
        &self,
        file: &mut File<R>,
        region: &Region,
        resolution: u32,
    ) -> Result<Vec<Pair>, ReadError> {
        if !self.resolutions.contains(&resolution) {
            return Err(ReadError::whole(format!(
                "the .hic has no {}-base resolution; it holds {}",
                group_thousands(u64::from(resolution)),
                listed(&self.resolutions)
            )));
        }
        let Some(at) = self
            .sequences
            .iter()
            .position(|(name, _)| name == region.seq())
        else {
            let names: Vec<&str> = self
                .sequences
                .iter()
                .map(|(name, _)| name.as_str())
                .collect();
            return Err(ReadError::whole(format!(
                "the .hic has no sequence called {}; it has {}",
                region.seq(),
                some_of(&names)
            )));
        };
        let length = self.sequences[at].1;
        let number = self.numbers[at];
        let size = u64::from(resolution);
        // Every bin the window touches, and none past the sequence's last.
        let Some(last) = length.checked_sub(1).map(|end| end / size) else {
            return Ok(Vec::new());
        };
        let low = region.start() / size;
        let high = (region.end().saturating_sub(1) / size).min(last);
        if low > high {
            return Ok(Vec::new());
        }
        let Some((offset, bytes)) = master(file, self.footer, &format!("{number}_{number}"))?
        else {
            return Ok(Vec::new());
        };
        let matrix = file.read_at(offset, bytes)?;
        let Some(zoom) = zoom_of(&matrix, number, resolution)? else {
            return Ok(Vec::new());
        };
        let (first, last_stretch) = (low / zoom.bins, high / zoom.bins);
        let band = band_of(high - low, zoom.bins);
        let mut cells = Vec::new();
        let span = |bin: u64| (bin * size, (bin + 1).saturating_mul(size).min(length));
        for &(block, offset, bytes) in &zoom.blocks {
            let (out, along) = (block / zoom.columns, block % zoom.columns);
            if out > band || along < first || along > last_stretch {
                continue;
            }
            let packed = file.read_at(offset, bytes)?;
            let inflated = zlib(&packed, BLOCK_MOST).map_err(|error| {
                ReadError::whole(format!(
                    "the .hic's block at byte {offset} will not inflate, so it is damaged: {error}"
                ))
            })?;
            read_cells(&inflated, |x, y, count| {
                let (x, y) = match (u64::try_from(x), u64::try_from(y)) {
                    (Ok(x), Ok(y)) => (x.min(y), x.max(y)),
                    _ => return Err(damaged("holds a cell before the start of its sequence")),
                };
                if x >= low && y <= high {
                    cells.push(Pair::spans(span(x), span(y), f64::from(count)));
                }
                Ok(())
            })?;
        }
        cells.sort_by_key(|cell| (cell.first, cell.second));
        Ok(cells)
    }
}

/// The resolutions, finest first, as a sentence lists them.
pub(crate) fn listed(resolutions: &[u32]) -> String {
    let mut sorted: Vec<u32> = resolutions.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let words: Vec<String> = sorted
        .into_iter()
        .map(|size| group_thousands(u64::from(size)))
        .collect();
    match words.split_last() {
        None => "none in bases".to_string(),
        Some((only, [])) => only.clone(),
        Some((last, rest)) => format!("{} and {last}", rest.join(", ")),
    }
}

/// The first dozen names, and that there are more where there are, as a
/// bigWig's refusal lists them.
fn some_of(names: &[&str]) -> String {
    if names.is_empty() {
        return "none".to_string();
    }
    let mut shown: Vec<&str> = names.iter().copied().take(12).collect();
    if names.len() > shown.len() {
        shown.push("and more");
    }
    shown.join(", ")
}

/// A file damaged as `what` says.
fn damaged(what: &str) -> ReadError {
    ReadError::whole(format!("the .hic {what}, so it is damaged"))
}

/// The header, read from the start of the file.
fn read_header<R: Read + Seek>(file: &mut File<R>) -> Result<Header, ReadError> {
    let mut at = Reading::new(file, 0);
    if at.take(4).ok() != Some(&b"HIC\0"[..]) {
        return Err(ReadError::whole(
            "not a .hic: it does not start with HIC, as every Juicer contact map does",
        ));
    }
    let version = at.i32()?;
    if version != VERSION {
        return Err(ReadError::whole(format!(
            "the .hic is version {version}, and version {VERSION} is the one read; hictk \
             convert in.hic in.mcool, and then hictk convert in.mcool out.hic, write the same \
             map as version {VERSION}"
        )));
    }
    let footer = at.offset("its master index")?;
    let genome = at.string()?;
    // Where the normalisation vectors are indexed, and how long that index
    // is, which raw counts do not need.
    at.take(16)?;
    for _ in 0..at.count("attributes")? {
        at.string()?;
        at.string()?;
    }
    let mut sequences = Vec::new();
    let mut numbers = Vec::new();
    for number in 0..at.count("sequences")? {
        let name = at.string()?;
        let length = u64::try_from(at.i64()?)
            .map_err(|_| damaged(&format!("gives {name} a length below nought")))?;
        // A file of several resolutions opens with every sequence end to
        // end, in kilobases, under the name All.
        if number == 0 && name.eq_ignore_ascii_case("all") {
            continue;
        }
        sequences.push((name, length));
        numbers.push(number);
    }
    let mut resolutions = Vec::new();
    for _ in 0..at.count("resolutions")? {
        let size = at.i32()?;
        resolutions.push(
            u32::try_from(size)
                .ok()
                .filter(|size| *size > 0)
                .ok_or_else(|| damaged(&format!("lists a resolution of {size} bases")))?,
        );
    }
    // The resolutions in restriction fragments, and the sites of each
    // sequence they are counted in, come next, and are not read.
    Ok(Header {
        version,
        genome,
        sequences,
        resolutions,
        footer,
        numbers,
    })
}

/// Where the map named `key` is kept and how many bytes it is, from the
/// master index at `footer`, or `None` where the index does not name it,
/// which is a map with no contacts.
fn master<R: Read + Seek>(
    file: &mut File<R>,
    footer: u64,
    key: &str,
) -> Result<Option<(u64, u64)>, ReadError> {
    let mut at = Reading::new(file, footer);
    // How many bytes the index and the expected values after it take, which
    // the entries are read without.
    at.take(8)?;
    for _ in 0..at.count("maps")? {
        let name = at.string()?;
        let offset = at.offset("a map")?;
        let bytes = at.count("bytes of a map")?;
        if name == key {
            return Ok(Some((offset, bytes as u64)));
        }
    }
    Ok(None)
}

/// One resolution of a map: how its blocks are laid out, and where each is.
struct Zoom {
    /// `blockBinCount`, how many bins of the diagonal a block's stretch is.
    bins: u64,
    /// `blockColumnCount`, how many stretches a band holds, which a block's
    /// number is counted in.
    columns: u64,
    /// Each block, as its number, where it starts and how many bytes it is.
    blocks: Vec<(u64, u64, u64)>,
}

/// The resolution of `resolution` bases of a map of sequence `number` with
/// itself, from the map's bytes, or `None` where the map does not hold it.
fn zoom_of(matrix: &[u8], number: usize, resolution: u32) -> Result<Option<Zoom>, ReadError> {
    let mut at = Bytes::new(matrix, false, ".hic");
    let signed = |value: u32| value as i32;
    let (first, second) = (signed(at.u32()?), signed(at.u32()?));
    if usize::try_from(first).ok() != Some(number) || first != second {
        return Err(damaged(&format!(
            "points to the map of {first} with {second} for the map of {number} with itself"
        )));
    }
    for _ in 0..counted(signed(at.u32()?), "resolutions of a map")? {
        let unit = at.until_nought()?.to_vec();
        // Its number among the resolutions, then four statistics of it: the
        // sum of its counts, the cells that hold one, their spread and the
        // 95th percentile.
        at.take(20)?;
        let size = signed(at.u32()?);
        let bins = counted(signed(at.u32()?), "bins of a block")? as u64;
        let columns = counted(signed(at.u32()?), "columns of blocks")? as u64;
        let blocks = counted(signed(at.u32()?), "blocks")?;
        if unit != b"BP" || u32::try_from(size).ok() != Some(resolution) {
            // Sixteen bytes a block: its number, where it is, and its size.
            at.take(blocks.checked_mul(16).ok_or_else(|| short(".hic"))?)?;
            continue;
        }
        if bins == 0 || columns == 0 {
            return Err(damaged("lays its blocks out in stretches of no bins"));
        }
        let mut listed = Vec::with_capacity(blocks.min(at.left() / 16));
        for _ in 0..blocks {
            let block = counted(signed(at.u32()?), "a block's number")? as u64;
            let offset = u64::try_from(at.u64()? as i64)
                .map_err(|_| damaged("points to a block before its start"))?;
            let bytes = counted(signed(at.u32()?), "bytes of a block")? as u64;
            listed.push((block, offset, bytes));
        }
        return Ok(Some(Zoom {
            bins,
            columns,
            blocks: listed,
        }));
    }
    Ok(None)
}

/// A count as the file writes it, refused below nought.
fn counted(value: i32, what: &str) -> Result<usize, ReadError> {
    usize::try_from(value).map_err(|_| damaged(&format!("counts {value} {what}")))
}

/// The furthest band out from the diagonal a cell of two bins `distance`
/// apart is filed in, of blocks of `bins` bins.
///
/// A writer works this out in floating point too, Juicer's as a quotient of
/// two natural logarithms, and two ways of working it can fall either side of
/// a whole number. Within a hair of the next band, that band is read as well:
/// a block too many is read and its cells left out, where a block too few is
/// cells missing from the figure.
fn band_of(distance: u64, bins: u64) -> u64 {
    let depth = (1.0 + distance as f64 / std::f64::consts::SQRT_2 / bins as f64).log2();
    let whole = depth.floor();
    if depth - whole > 1.0 - 1e-9 {
        whole as u64 + 1
    } else {
        whole as u64
    }
}

/// Each cell of an inflated block, as its two bins and its count, handed to
/// `each`.
///
/// A block opens with its count of cells, which is not trusted for anything,
/// the bins its numbers are counted from, and four bytes saying how it is
/// written: whether its counts are 16-bit integers or 32-bit floats, whether
/// its column and its row bins are 16 or 32 bits, and whether it is rows of
/// the cells that hold something or a dense square. Each of the first three
/// is nought for the short one, the reverse of what a flag usually is, as
/// straw and hictk read it; hictk writes all three long.
fn read_cells(
    block: &[u8],
    mut each: impl FnMut(i64, i64, f32) -> Result<(), ReadError>,
) -> Result<(), ReadError> {
    let mut at = Bytes::new(block, false, ".hic block");
    at.u32()?;
    let x0 = i64::from(at.u32()? as i32);
    let y0 = i64::from(at.u32()? as i32);
    let short_counts = at.u8()? == 0;
    let short_x = at.u8()? == 0;
    let short_y = at.u8()? == 0;
    let kind = at.u8()?;
    let number = |at: &mut Bytes, short: bool| -> Result<i64, ReadError> {
        Ok(if short {
            i64::from(at.u16()? as i16)
        } else {
            i64::from(at.u32()? as i32)
        })
    };
    let count = |at: &mut Bytes| -> Result<f32, ReadError> {
        Ok(if short_counts {
            f32::from(at.u16()? as i16)
        } else {
            at.f32()?
        })
    };
    match kind {
        // Rows: how many, then each row's bin, how many cells it holds, and
        // each cell's column bin and count.
        1 => {
            let rows = number(&mut at, short_y)?;
            for _ in 0..rows.max(0) {
                let y = y0 + number(&mut at, short_y)?;
                let cells = number(&mut at, short_x)?;
                for _ in 0..cells.max(0) {
                    let x = x0 + number(&mut at, short_x)?;
                    each(x, y, count(&mut at)?)?;
                }
            }
        }
        // A dense square, `width` cells to a row, row after row, with a cell
        // that holds nothing written as the least 16-bit integer or as NaN.
        2 => {
            let cells = i64::from(at.u32()? as i32);
            let width = i64::from(at.u16()? as i16);
            if cells > 0 && width <= 0 {
                return Err(damaged(&format!("writes a dense block {width} cells wide")));
            }
            for cell in 0..cells.max(0) {
                let value = if short_counts {
                    let value = at.u16()? as i16;
                    (value != i16::MIN).then_some(f32::from(value))
                } else {
                    Some(at.f32()?).filter(|value| !value.is_nan())
                };
                if let Some(value) = value {
                    each(x0 + cell % width, y0 + cell / width, value)?;
                }
            }
        }
        other => {
            return Err(damaged(&format!(
                "holds a block of kind {other}, which the format does not have"
            )))
        }
    }
    Ok(())
}

/// A file read from an offset on, a stretch at a time, where how much of it
/// is wanted is only known as it is read. Nothing is read past the end of
/// the file, and a string runs at most [`LONGEST`] bytes.
struct Reading<'f, R> {
    file: &'f mut File<R>,
    /// The bytes read and not yet taken, from `at` on.
    held: Vec<u8>,
    /// Where in the file `held` starts.
    from: u64,
    /// How many of `held` are taken.
    at: usize,
}

impl<'f, R: Read + Seek> Reading<'f, R> {
    fn new(file: &'f mut File<R>, offset: u64) -> Self {
        Reading {
            file,
            held: Vec::new(),
            from: offset,
            at: 0,
        }
    }

    /// At least `n` bytes ahead where the file holds them, and as many as it
    /// does otherwise: how many there are.
    fn fill(&mut self, n: usize) -> Result<usize, ReadError> {
        let ahead = self.held.len() - self.at;
        if ahead >= n {
            return Ok(ahead);
        }
        let offset = self
            .from
            .checked_add(self.at as u64)
            .ok_or_else(|| short(self.file.what))?;
        self.held = self.file.read_up_to(offset, (n as u64).max(CHUNK))?;
        self.from = offset;
        self.at = 0;
        Ok(self.held.len())
    }

    fn take(&mut self, n: usize) -> Result<&[u8], ReadError> {
        if self.fill(n)? < n {
            return Err(short(self.file.what));
        }
        let taken = &self.held[self.at..self.at + n];
        self.at += n;
        Ok(taken)
    }

    fn i32(&mut self) -> Result<i32, ReadError> {
        let bytes = self.take(4)?;
        Ok(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn i64(&mut self) -> Result<i64, ReadError> {
        let bytes = self.take(8)?;
        let mut eight = [0u8; 8];
        eight.copy_from_slice(bytes);
        Ok(i64::from_le_bytes(eight))
    }

    /// A count, refused below nought.
    fn count(&mut self, what: &str) -> Result<usize, ReadError> {
        let value = self.i32()?;
        counted(value, what)
    }

    /// An offset into the file, refused below nought.
    fn offset(&mut self, what: &str) -> Result<u64, ReadError> {
        let value = self.i64()?;
        u64::try_from(value).map_err(|_| damaged(&format!("puts {what} before its start")))
    }

    /// A string ending in a nought, which is passed over.
    fn string(&mut self) -> Result<String, ReadError> {
        let mut looked = 0;
        loop {
            let ahead = &self.held[self.at..];
            if let Some(end) = ahead[looked..].iter().position(|byte| *byte == 0) {
                let text = String::from_utf8_lossy(&ahead[..looked + end]).into_owned();
                self.at += looked + end + 1;
                return Ok(text);
            }
            looked = ahead.len();
            if looked >= LONGEST {
                return Err(damaged(&format!(
                    "holds a name or an attribute longer than {} MB",
                    LONGEST >> 20
                )));
            }
            if self.fill(looked + CHUNK as usize)? == looked {
                return Err(short(self.file.what));
            }
        }
    }
}

#[cfg(test)]
pub(crate) mod fixture {
    //! A `.hic` written by hand, for what hictk does not write: short bins
    //! and counts, cells that hold nothing in a dense block, and damage.

    /// A zlib stream of `data` stored as it is, which needs no compressor.
    pub(crate) fn stored(data: &[u8]) -> Vec<u8> {
        let mut out = vec![0x78, 0x01];
        let mut runs = data.chunks(65_535).peekable();
        if runs.peek().is_none() {
            out.extend([1, 0, 0, 0xff, 0xff]);
        }
        while let Some(run) = runs.next() {
            out.push(u8::from(runs.peek().is_none()));
            let length = run.len() as u16;
            out.extend(length.to_le_bytes());
            out.extend((!length).to_le_bytes());
            out.extend(run);
        }
        out.extend(crate::read::gzip::adler32(data).to_be_bytes());
        out
    }

    fn string(out: &mut Vec<u8>, text: &str) {
        out.extend(text.as_bytes());
        out.push(0);
    }

    /// A version 9 file of one sequence, `name` of `length` bases, holding
    /// its map with itself at `resolution` bases a bin in blocks of `bins`
    /// bins and `columns` stretches, each block its number and its bytes
    /// before they are compressed.
    pub(crate) fn file(
        name: &str,
        length: i64,
        resolution: i32,
        (bins, columns): (i32, i32),
        blocks: &[(i32, Vec<u8>)],
    ) -> Vec<u8> {
        let mut out = b"HIC\0".to_vec();
        out.extend(9i32.to_le_bytes());
        let footer_at = out.len();
        out.extend(0i64.to_le_bytes());
        string(&mut out, "test");
        out.extend([0u8; 16]);
        out.extend(1i32.to_le_bytes());
        string(&mut out, "software");
        string(&mut out, "by hand");
        out.extend(1i32.to_le_bytes());
        string(&mut out, name);
        out.extend(length.to_le_bytes());
        out.extend(1i32.to_le_bytes());
        out.extend(resolution.to_le_bytes());
        out.extend(0i32.to_le_bytes());
        let mut index = Vec::new();
        for (number, block) in blocks {
            let packed = stored(block);
            index.push((*number, out.len() as i64, packed.len() as i32));
            out.extend(packed);
        }
        let matrix_at = out.len() as i64;
        out.extend(0i32.to_le_bytes());
        out.extend(0i32.to_le_bytes());
        out.extend(1i32.to_le_bytes());
        string(&mut out, "BP");
        out.extend([0u8; 20]);
        out.extend(resolution.to_le_bytes());
        out.extend(bins.to_le_bytes());
        out.extend(columns.to_le_bytes());
        out.extend((index.len() as i32).to_le_bytes());
        for (number, offset, size) in index {
            out.extend(number.to_le_bytes());
            out.extend(offset.to_le_bytes());
            out.extend(size.to_le_bytes());
        }
        let matrix_size = out.len() as i64 - matrix_at;
        let footer = out.len() as i64;
        out[footer_at..footer_at + 8].copy_from_slice(&footer.to_le_bytes());
        out.extend(0i64.to_le_bytes());
        out.extend(1i32.to_le_bytes());
        string(&mut out, "0_0");
        out.extend(matrix_at.to_le_bytes());
        out.extend((matrix_size as i32).to_le_bytes());
        out
    }

    /// A block's bytes before they are compressed: its offsets, its four
    /// bytes of how it is written, and then `body`.
    pub(crate) fn block(
        (x0, y0): (i32, i32),
        (short_counts, short_x, short_y): (bool, bool, bool),
        kind: u8,
        body: &[u8],
    ) -> Vec<u8> {
        let mut out = 0i32.to_le_bytes().to_vec();
        out.extend(x0.to_le_bytes());
        out.extend(y0.to_le_bytes());
        // Nought for the short one, as the format has it.
        out.extend([
            u8::from(!short_counts),
            u8::from(!short_x),
            u8::from(!short_y),
            kind,
        ]);
        out.extend(body);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::{block, file};
    use super::*;
    use std::io::{Cursor, SeekFrom};

    const CONTACTS: &[u8] = include_bytes!("fixtures/hic/contacts.hic");
    const DUMP: &str = include_str!("fixtures/hic/contacts.dump");
    const SQUARE: &[u8] = include_bytes!("fixtures/hic/square.hic");
    const SQUARE_DUMP: &str = include_str!("fixtures/hic/square.dump");

    /// Each window `make.sh` asked `hictk dump` for, in `dump`: its
    /// resolution, its region, 0-based and half-open as hictk takes one, and
    /// the cells it printed as BEDPE.
    fn dumped(dump: &str) -> Vec<(u32, Region, String)> {
        let mut out: Vec<(u32, Region, String)> = Vec::new();
        for line in dump.lines() {
            if let Some(window) = line.strip_prefix("# ") {
                let words: Vec<&str> = window.split(' ').collect();
                let region = Region::new(
                    words[1],
                    words[2].parse().unwrap(),
                    words[3].parse().unwrap(),
                )
                .unwrap();
                out.push((words[0].parse().unwrap(), region, String::new()));
            } else {
                let last = out.last_mut().unwrap();
                last.2.push_str(line);
                last.2.push('\n');
            }
        }
        out
    }

    /// Cells as `(start, end, start, end, count)`, for comparing.
    fn cells_of(bedpe: &str) -> Vec<(u64, u64, u64, u64, f64)> {
        let mut out: Vec<(u64, u64, u64, u64, f64)> = bedpe
            .lines()
            .map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                (
                    f[1].parse().unwrap(),
                    f[2].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                )
            })
            .collect();
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out
    }

    /// Every window hictk dumped is read cell for cell, count for count: a
    /// window inside one block and one across every block of the finest
    /// resolution, both bands out from the diagonal included; windows whose
    /// edges fall inside a bin; the last bin, cut where the sequence ends;
    /// a dense block; and a window with nothing in it.
    #[test]
    fn every_window_hictk_dumps_is_read_cell_for_cell() {
        let cases = dumped(DUMP);
        assert_eq!(cases.len(), 9);
        for (resolution, region, expected) in cases {
            let read = contacts(Cursor::new(CONTACTS), &region, resolution).unwrap();
            let written = bedpe(region.seq(), &read);
            assert_eq!(
                cells_of(&written),
                cells_of(&expected),
                "{resolution} {region}"
            );
            assert!(written.lines().all(|line| line.starts_with(region.seq())));
        }
    }

    #[test]
    fn the_header_names_its_sequences_and_resolutions_and_leaves_out_all() {
        let header = header_of(Cursor::new(CONTACTS)).unwrap();
        assert_eq!(header.version, 9);
        assert_eq!(header.genome, "unknown");
        assert_eq!(
            header.sequences,
            [
                ("chr1".to_string(), 1_234_567),
                ("chr2".to_string(), 345_678)
            ]
        );
        assert_eq!(
            header.resolutions,
            [1_000, 2_000, 5_000, 10_000, 50_000, 250_000]
        );
        // All is the file's sequence 0, so chr1's map is 1_1.
        assert_eq!(header.numbers, [1, 2]);
        assert_eq!(sequences(Cursor::new(CONTACTS)).unwrap(), header.sequences);
    }

    /// The last bin of a sequence stops where the sequence does, and a
    /// window past its end holds nothing.
    #[test]
    fn the_last_bin_is_cut_at_the_end_of_the_sequence() {
        let region = Region::new("chr1", 1_200_000, 1_300_000).unwrap();
        let read = contacts(Cursor::new(CONTACTS), &region, 250_000).unwrap();
        assert_eq!(
            read,
            [Pair::spans(
                (1_000_000, 1_234_567),
                (1_000_000, 1_234_567),
                2_038.5
            )]
        );
        let past = Region::new("chr1", 1_234_567, 1_300_000).unwrap();
        assert!(contacts(Cursor::new(CONTACTS), &past, 1_000)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn the_resolution_is_the_finest_that_keeps_the_window_to_250_bins() {
        let held = [250_000, 1_000, 5_000, 2_000, 50_000, 10_000];
        let whole = Region::new("chr1", 0, 1_234_567).unwrap();
        // 1,235 bins at 1 kb and 618 at 2 kb; 247 at 5 kb.
        assert_eq!(resolution_for(&held, &whole, BINS), Some(5_000));
        // 250 bins exactly is kept to; one more is not.
        let fits = Region::new("chr1", 0, 250_000).unwrap();
        assert_eq!(resolution_for(&held, &fits, BINS), Some(1_000));
        let over = Region::new("chr1", 0, 250_001).unwrap();
        assert_eq!(resolution_for(&held, &over, BINS), Some(2_000));
        // A window whose edges fall inside bins touches one bin more.
        let shifted = Region::new("chr1", 500, 250_000).unwrap();
        assert_eq!(resolution_for(&held, &shifted, BINS), Some(1_000));
        let across = Region::new("chr1", 500, 250_501).unwrap();
        assert_eq!(resolution_for(&held, &across, BINS), Some(2_000));
        // Too wide for any, the coarsest; and none where there is none.
        let genome = Region::new("chr1", 0, 3_000_000_000).unwrap();
        assert_eq!(resolution_for(&held, &genome, BINS), Some(250_000));
        assert_eq!(resolution_for(&[], &whole, BINS), None);
    }

    /// Reads and seeks a file, keeping where each read started.
    struct Watched<'a> {
        inner: Cursor<&'a [u8]>,
        reads: Vec<u64>,
    }

    impl Read for Watched<'_> {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            self.reads.push(self.inner.position());
            self.inner.read(out)
        }
    }

    impl Seek for Watched<'_> {
        fn seek(&mut self, to: SeekFrom) -> std::io::Result<u64> {
            self.inner.seek(to)
        }
    }

    /// The blocks of chr1's map at 1 kb, as hictk wrote them: where each
    /// starts, by its number.
    fn blocks_of_chr1() -> Vec<(u64, u64)> {
        let header = header_of(Cursor::new(CONTACTS)).unwrap();
        let mut file = File::new(Cursor::new(CONTACTS), ".hic").unwrap();
        let (offset, bytes) = master(&mut file, header.footer, "1_1").unwrap().unwrap();
        let matrix = file.read_at(offset, bytes).unwrap();
        let zoom = zoom_of(&matrix, 1, 1_000).unwrap().unwrap();
        assert_eq!((zoom.bins, zoom.columns), (618, 2));
        zoom.blocks
            .iter()
            .map(|(number, offset, _)| (*number, *offset))
            .collect()
    }

    /// A window reads the blocks that can hold its cells and no others: one
    /// near the start, the two stretches of the band on the diagonal for a
    /// window across them, and the band out from it only for a window wide
    /// enough to hold a cell there, as the whole sequence is.
    #[test]
    fn a_window_reads_only_the_blocks_over_it() {
        let blocks = blocks_of_chr1();
        assert_eq!(
            blocks.iter().map(|(number, _)| *number).collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        for (start, end, wanted) in [
            (300_000, 420_000, vec![0]),
            (900_000, 1_000_000, vec![1]),
            (500_000, 800_000, vec![0, 1]),
            (0, 1_234_567, vec![0, 1, 2, 3]),
            // Cells 849 bins apart at most, which the band on the diagonal
            // holds, and 899, which reach the next.
            (0, 850_000, vec![0, 1]),
            (0, 900_000, vec![0, 1, 2, 3]),
        ] {
            let mut watched = Watched {
                inner: Cursor::new(CONTACTS),
                reads: Vec::new(),
            };
            let region = Region::new("chr1", start, end).unwrap();
            contacts(&mut watched, &region, 1_000).unwrap();
            let read: Vec<u64> = blocks
                .iter()
                .filter(|(_, offset)| watched.reads.contains(offset))
                .map(|(number, _)| *number)
                .collect();
            assert_eq!(read, wanted, "{region}");
        }
    }

    #[test]
    fn a_version_other_than_9_is_refused_naming_it_and_how_to_write_9() {
        for version in [6u8, 7, 8, 10] {
            let mut bytes = CONTACTS.to_vec();
            bytes[4] = version;
            let error = header_of(Cursor::new(&bytes)).unwrap_err().to_string();
            assert!(
                error.starts_with(&format!("the .hic is version {version}, ")),
                "{error}"
            );
            assert!(error.contains("hictk convert in.mcool out.hic"), "{error}");
        }
        let error = header_of(Cursor::new(b"HIK\0\x09\0\0\0".to_vec())).unwrap_err();
        assert!(error.to_string().starts_with("not a .hic"), "{error}");
    }

    #[test]
    fn a_sequence_or_a_resolution_the_file_has_not_got_is_refused_naming_its_own() {
        let region = Region::new("chr3", 0, 1_000).unwrap();
        let error = contacts(Cursor::new(CONTACTS), &region, 1_000).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the .hic has no sequence called chr3; it has chr1, chr2"
        );
        let region = Region::new("chr1", 0, 1_000).unwrap();
        let error = contacts(Cursor::new(CONTACTS), &region, 7_000).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the .hic has no 7,000-base resolution; it holds 1,000, 2,000, 5,000, 10,000, \
             50,000 and 250,000"
        );
        // All is no sequence of the genome.
        let all = Region::new("All", 0, 1_000).unwrap();
        assert!(contacts(Cursor::new(CONTACTS), &all, 1_000).is_err());
    }

    /// Bins and counts written short read as they do written long, a dense
    /// block leaves out the cells that hold nothing, and the cells of a
    /// block are counted from its offsets.
    #[test]
    fn short_bins_short_counts_and_dense_blocks_read_as_written() {
        // Rows from bin 10: the first holds columns 0 and 2 after it, the
        // second column 1, each count after its column's bin.
        let written = |short_bins: bool, short_counts: bool| {
            let mut body = Vec::new();
            let bin = |body: &mut Vec<u8>, value: i32| {
                if short_bins {
                    body.extend((value as i16).to_le_bytes());
                } else {
                    body.extend(value.to_le_bytes());
                }
            };
            let count = |body: &mut Vec<u8>, value: i16| {
                if short_counts {
                    body.extend(value.to_le_bytes());
                } else {
                    body.extend(f32::from(value).to_le_bytes());
                }
            };
            bin(&mut body, 2);
            bin(&mut body, 0);
            bin(&mut body, 2);
            bin(&mut body, 0);
            count(&mut body, 7);
            bin(&mut body, 2);
            count(&mut body, 3);
            bin(&mut body, 1);
            bin(&mut body, 1);
            bin(&mut body, 1);
            count(&mut body, 9);
            block((10, 10), (short_counts, short_bins, short_bins), 1, &body)
        };
        let expected = [
            Pair::spans((100, 110), (100, 110), 7.0),
            Pair::spans((100, 110), (120, 125), 3.0),
            Pair::spans((110, 120), (110, 120), 9.0),
        ];
        let region = Region::new("s", 0, 125).unwrap();
        for short_bins in [false, true] {
            for short_counts in [false, true] {
                let bytes = file(
                    "s",
                    125,
                    10,
                    (20, 1),
                    &[(0, written(short_bins, short_counts))],
                );
                let read = contacts(Cursor::new(&bytes), &region, 10).unwrap();
                assert_eq!(read, expected, "{short_bins} {short_counts}");
            }
        }
        // Dense, two cells to a row from bin 10: the least 16-bit integer,
        // or NaN, is a cell that holds nothing.
        let dense = |short_counts: bool| {
            let mut body = 4i32.to_le_bytes().to_vec();
            body.extend(2i16.to_le_bytes());
            for value in [Some(7i16), None, None, Some(9)] {
                if short_counts {
                    body.extend(value.unwrap_or(i16::MIN).to_le_bytes());
                } else {
                    body.extend(value.map_or(f32::NAN, f32::from).to_le_bytes());
                }
            }
            block((10, 10), (short_counts, false, false), 2, &body)
        };
        for short_counts in [false, true] {
            let bytes = file("s", 125, 10, (20, 1), &[(0, dense(short_counts))]);
            let read = contacts(Cursor::new(&bytes), &region, 10).unwrap();
            assert_eq!(
                read,
                [
                    Pair::spans((100, 110), (100, 110), 7.0),
                    Pair::spans((110, 120), (110, 120), 9.0),
                ],
                "{short_counts}"
            );
        }
    }

    /// A dense block off the diagonal, as hictk writes a square of contacts
    /// far from it, reads cell for cell as hictk dumps it, whole and cut by
    /// a window. Its columns start at bin 0 and its rows at bin 40, and every
    /// count is a different number, so a cell read with its column and its
    /// row swapped lands on another cell's count: a dense block on the
    /// diagonal, which is all `contacts.hic` holds, is its own mirror and
    /// cannot show it.
    #[test]
    fn a_dense_block_off_the_diagonal_reads_as_hictk_dumps_it() {
        let header = header_of(Cursor::new(SQUARE)).unwrap();
        let mut file = File::new(Cursor::new(SQUARE), ".hic").unwrap();
        let (offset, bytes) = master(&mut file, header.footer, "0_0").unwrap().unwrap();
        let matrix = file.read_at(offset, bytes).unwrap();
        let zoom = zoom_of(&matrix, 0, 1_000).unwrap().unwrap();
        let [(_, offset, bytes)] = zoom.blocks[..] else {
            panic!("{} blocks", zoom.blocks.len());
        };
        let block = zlib(&file.read_at(offset, bytes).unwrap(), BLOCK_MOST).unwrap();
        // Its offsets, and dense.
        assert_eq!(block[4..12], [0, 0, 0, 0, 40, 0, 0, 0]);
        assert_eq!(block[15], 2);
        let cases = dumped(SQUARE_DUMP);
        assert_eq!(cases.len(), 2);
        for (resolution, region, expected) in cases {
            let read = contacts(Cursor::new(SQUARE), &region, resolution).unwrap();
            assert!(!read.is_empty(), "{region}");
            assert_eq!(
                cells_of(&bedpe(region.seq(), &read)),
                cells_of(&expected),
                "{region}"
            );
        }
    }

    /// A block off the diagonal reads its columns as columns and its rows as
    /// rows, written each way a block can be: dense, three cells to a row,
    /// where a cell's column is its place in its row; and in rows, with the
    /// column bins and the row bins each as wide as their own byte says,
    /// one short and the other long both ways round. Columns start at bin 10
    /// and rows at bin 30, and every count is a different number, so a cell
    /// read in its mirror's place, or bins read at the other's width, is a
    /// cell in the wrong place or a block refused.
    #[test]
    fn a_block_off_the_diagonal_reads_its_columns_and_its_rows_the_right_way_round() {
        let region = Region::new("s", 0, 400).unwrap();
        let read = |block: Vec<u8>| {
            let bytes = file("s", 400, 10, (40, 1), &[(0, block)]);
            contacts(Cursor::new(&bytes), &region, 10).unwrap()
        };
        let cell = |column: u64, row: u64, count: f64| {
            Pair::spans(
                (column * 10, column * 10 + 10),
                (row * 10, row * 10 + 10),
                count,
            )
        };
        // Two rows of three: 1, 2, 3 in row 30 and 4, 5, 6 in row 31.
        let mut body = 6i32.to_le_bytes().to_vec();
        body.extend(3i16.to_le_bytes());
        for count in 1..=6u8 {
            body.extend(f32::from(count).to_le_bytes());
        }
        assert_eq!(
            read(block((10, 30), (false, false, false), 2, &body)),
            [
                cell(10, 30, 1.0),
                cell(10, 31, 4.0),
                cell(11, 30, 2.0),
                cell(11, 31, 5.0),
                cell(12, 30, 3.0),
                cell(12, 31, 6.0),
            ]
        );
        // Row 30 holds columns 10 and 12, counting 1 and 3, and row 31
        // column 11, counting 5.
        for (short_x, short_y) in [(false, false), (true, false), (false, true), (true, true)] {
            let bin = |body: &mut Vec<u8>, short: bool, value: i32| {
                if short {
                    body.extend((value as i16).to_le_bytes());
                } else {
                    body.extend(value.to_le_bytes());
                }
            };
            let mut body = Vec::new();
            bin(&mut body, short_y, 2);
            bin(&mut body, short_y, 0);
            bin(&mut body, short_x, 2);
            bin(&mut body, short_x, 0);
            body.extend(1f32.to_le_bytes());
            bin(&mut body, short_x, 2);
            body.extend(3f32.to_le_bytes());
            bin(&mut body, short_y, 1);
            bin(&mut body, short_x, 1);
            bin(&mut body, short_x, 1);
            body.extend(5f32.to_le_bytes());
            assert_eq!(
                read(block((10, 30), (false, short_x, short_y), 1, &body)),
                [cell(10, 30, 1.0), cell(11, 31, 5.0), cell(12, 30, 3.0)],
                "{short_x} {short_y}"
            );
        }
    }

    /// Every way of cutting the file short is refused, or read as the whole
    /// file reads where the cut leaves every byte the window needs; every
    /// byte of it turned over is refused or read, never a panic; and a block
    /// that says it holds two thousand million rows in a few bytes is refused
    /// once its bytes run out, with nothing allocated for them.
    #[test]
    fn a_damaged_file_is_refused_not_a_panic() {
        let region = Region::new("chr1", 0, 1_234_567).unwrap();
        let whole = contacts(Cursor::new(CONTACTS), &region, 1_000).unwrap();
        let mut refused = 0;
        // Every seventh length and every seventh byte, which in a debug build
        // is a second rather than eleven, and still lands in every field of
        // the header, the indexes and the blocks.
        for length in (0..CONTACTS.len()).step_by(7) {
            match contacts(Cursor::new(&CONTACTS[..length]), &region, 1_000) {
                Ok(read) => assert_eq!(read, whole, "{length}"),
                Err(_) => refused += 1,
            }
        }
        // Past the master index's entry for chr1, the rest of the file is
        // what other windows and sequences need.
        assert!(refused > 11_000 / 7, "{refused}");
        for at in (0..CONTACTS.len()).step_by(7) {
            let mut turned = CONTACTS.to_vec();
            turned[at] = !turned[at];
            let _ = contacts(Cursor::new(&turned), &region, 1_000);
            let _ = header_of(Cursor::new(&turned));
        }
        let mut body = i32::MAX.to_le_bytes().to_vec();
        body.extend(0i32.to_le_bytes());
        body.extend(i32::MAX.to_le_bytes());
        let huge = file(
            "s",
            125,
            10,
            (20, 1),
            &[(0, block((0, 0), (false, false, false), 1, &body))],
        );
        let small = Region::new("s", 0, 125).unwrap();
        let error = contacts(Cursor::new(&huge), &small, 10).unwrap_err();
        assert!(error.to_string().contains("cut short"), "{error}");
        let strange = file(
            "s",
            125,
            10,
            (20, 1),
            &[(0, block((0, 0), (false, false, false), 3, &[]))],
        );
        let error = contacts(Cursor::new(&strange), &small, 10).unwrap_err();
        assert!(error.to_string().contains("block of kind 3"), "{error}");
        let before = file(
            "s",
            125,
            10,
            (20, 1),
            &[(
                0,
                block((-5, 0), (false, false, false), 2, &{
                    let mut body = 1i32.to_le_bytes().to_vec();
                    body.extend(1i16.to_le_bytes());
                    body.extend(1f32.to_le_bytes());
                    body
                }),
            )],
        );
        let error = contacts(Cursor::new(&before), &small, 10).unwrap_err();
        assert!(error.to_string().contains("before the start"), "{error}");
    }

    /// The BEDPE written for the cells is what the reader of pairs reads
    /// back as the same cells.
    #[test]
    fn the_bedpe_written_reads_back_as_the_cells() {
        let region = Region::new("chr1", 612_345, 987_654).unwrap();
        let read = contacts(Cursor::new(CONTACTS), &region, 5_000).unwrap();
        assert!(read.len() > 50);
        let (back, _) = crate::read::pairs::pairs(&bedpe("chr1", &read), &region).unwrap();
        assert_eq!(back, read);
    }
}
