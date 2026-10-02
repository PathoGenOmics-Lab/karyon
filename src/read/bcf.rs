//! BCF, the binary form of VCF, read a window at a time through the CSI
//! beside it and handed over as the VCF text the readers of calls take.
//!
//! A BCF holds what a VCF holds, a record to a site and a column of values to
//! a sample, with the numbers stored as numbers and every name a VCF spells
//! out on each row, a sequence, a filter, an INFO or a FORMAT key, stored as
//! its place in a dictionary the header keeps. It is written in BGZF, as a
//! bgzipped VCF is, and `bcftools index` writes a CSI beside it, so [`window`]
//! reads the few blocks of a whole genome's calls that hold the records over a
//! window, and the whole file only where no index is given.
//!
//! What comes out is the text `bcftools view` prints for the same records, the
//! header first and then a row per record over the window, so `--variants`,
//! `--genotypes` and `--structural` read a BCF through the readers they read a
//! VCF through and draw what they draw from its text. A test holds every
//! record of the files in `fixtures/bcf` to what `bcftools view` 1.24 prints
//! for it, byte for byte, numbers included: a float is printed as htslib
//! prints one, six significant digits by its own rounding, since a value
//! printed otherwise is a lollipop at another height.
//!
//! Most of a cohort's BCF is its samples' columns, and a track of calls reads
//! none of them, so [`Fields`] says how much of each record to write: the
//! eight columns of a site, as `bcftools view -G` prints them, which leaves
//! the samples' bytes undecoded; those and each sample's `GT` alone, which
//! is all a track of genotypes reads; or every field, as `bcftools view`
//! prints them.
//!
//! # Coordinates
//!
//! A record's position is stored 0-based, where the VCF text it stands for
//! writes it 1-based, so a record at 9 is written at 10 and read back by the
//! readers of VCF at 9 again.
//!
//! # What is read
//!
//! BCF 2.2, the one version htslib reads and bcftools has written since 2014,
//! in BGZF as `bcftools view -Ob` writes it, in BGZF blocks stored rather than
//! compressed as `-Ou` writes it, and as the bare stream `-Ou` wrote before
//! htslib 1.10 and gzip turns BGZF into. A file compressed with gzip rather
//! than bgzip is read whole once inflated, since an index has no blocks in it
//! to point to.
//!
//! # What is refused
//!
//! A file that is not BCF, one of another version, and one damaged or cut
//! short: a record whose parts do not fit the bytes it says it holds, a value
//! of a type BCF does not have, a sequence or a key the header's dictionaries
//! do not hold, and a record of more or fewer samples than the header names.
//! None of it panics, and no count read from the file allocates more than the
//! bytes it has left could fill.
//!
//! ```
//! use std::io::Cursor;
//! use karyon::{read, Region};
//!
//! # let bcf = include_bytes!("fixtures/bcf/tiny.bcf").to_vec();
//! let region = Region::parse("chr1:1-100")?;
//! let text = read::bcf::window(Cursor::new(bcf), None, &region, read::bcf::Fields::Sites)?;
//! let calls = read::point::variants(&text, &region)?;
//! assert_eq!(calls[0].pos, 9);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::io::{Cursor, Read, Seek, SeekFrom};

use super::bgzf::{self, Bgzf};
use super::bytes::Bytes;
use super::index::{Index, Kind};
use super::{gzip, ReadError};
use crate::Region;

/// How much of each record [`window`] writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fields {
    /// The eight columns of a site, and no column of the samples, as
    /// `bcftools view -G` prints them. The samples' bytes are passed over
    /// undecoded, which is most of a cohort's file.
    Sites,
    /// The site, and of each sample its `GT` alone, as `bcftools annotate -x
    /// ^FORMAT/GT` prints them: a record with no `GT` has `.` for its
    /// `FORMAT` and for every sample.
    Genotypes,
    /// Every field of every sample, as `bcftools view` prints them.
    All,
}

/// What a BCF says about itself before its records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// The header as `bcftools view` prints it: the text the file holds, less
    /// the `IDX` each line carries in a BCF to number its dictionaries.
    pub text: String,
    /// The sequences, each with the length its `##contig` line gives, in the
    /// order the records number them.
    pub contigs: Vec<(String, Option<u64>)>,
    /// The samples, as the `#CHROM` line names them, in the order of their
    /// columns.
    pub samples: Vec<String>,
    /// The sequence of each number a record can carry, where the header
    /// numbers one: [`Header::contigs`] by number, with the gaps a header
    /// with lines taken out of it leaves.
    sequences: Vec<Option<usize>>,
    /// The name of each number a FILTER, an INFO or a FORMAT key is stored
    /// as, `PASS` at nought.
    keys: Vec<Option<String>>,
    /// Whether the header says VCF 4.4 or later, in which a `GT` may say how
    /// its first allele is phased, and is printed saying so where the rest of
    /// the call does not imply it.
    prefixed: bool,
}

impl Header {
    /// The number a record and an index give the sequence called `name`, or
    /// `None` for one the header does not name.
    pub fn sequence(&self, name: &str) -> Option<usize> {
        let at = self.contigs.iter().position(|(held, _)| held == name)?;
        self.sequences
            .iter()
            .position(|sequence| *sequence == Some(at))
    }

    /// The sequence a record numbered `number` is on.
    fn named(&self, number: i32) -> Option<&str> {
        let at = (*self.sequences.get(usize::try_from(number).ok()?)?)?;
        Some(self.contigs[at].0.as_str())
    }

    /// The name a FILTER, INFO or FORMAT key is stored as `number` under.
    fn key(&self, number: i64) -> Result<&str, ReadError> {
        usize::try_from(number)
            .ok()
            .and_then(|at| self.keys.get(at)?.as_deref())
            .ok_or_else(|| {
                ReadError::whole(format!(
                    "a BCF record names key {number}, which the header's dictionary does not \
                     hold"
                ))
            })
    }

    /// The header as [`window`] writes it for `fields`: whole, or for sites
    /// alone as `bcftools view -G` prints it, with no FORMAT line and the
    /// `#CHROM` line cut after `INFO`.
    fn written(&self, fields: Fields) -> String {
        if fields != Fields::Sites {
            return self.text.clone();
        }
        let start = match self.text.rfind("\n#CHROM") {
            Some(at) => at + 1,
            None => 0,
        };
        let mut text: String = self.text[..start]
            .split_inclusive('\n')
            .filter(|line| !line.starts_with("##FORMAT="))
            .collect();
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        text
    }
}

/// The header of a BCF.
///
/// # Errors
///
/// A file that is not BCF 2.2, or whose header is damaged or cut short.
pub fn header_of<R: Read + Seek>(reader: R) -> Result<Header, ReadError> {
    let mut stream = open(reader)?;
    header(&mut stream)
}

/// The VCF text of the records of a BCF over `region`: the header, then every
/// record over the window, as `bcftools view` prints them, cut to `fields`.
///
/// A record is over the window where the bases it spans touch it: from its
/// position as far as its reference allele spells, or as its `END` says,
/// whichever is further, which is the readers' rule and tabix's together. The
/// readers drop the ones a rule of their own leaves out, as they drop the rows
/// of a whole VCF outside the window.
///
/// Read through `index` where one is given and the file is BGZF, which a CSI
/// points into: only the blocks that hold records over the window are read,
/// and reading stops at the first record past it. Without one, every record
/// is read, and kept where it is over the window, since a file with no index
/// need not be sorted. A window on a sequence the header does not name holds
/// no record, as a VCF's rows on another sequence go past without a word.
///
/// # Errors
///
/// A file that is not BCF or is damaged, and an index that is not a CSI,
/// the one index a BCF has.
pub fn window<R: Read + Seek>(
    reader: R,
    index: Option<&Index>,
    region: &Region,
    fields: Fields,
) -> Result<String, ReadError> {
    let mut stream = open(reader)?;
    let header = header(&mut stream)?;
    let mut out = header.written(fields).into_bytes();
    let Some(number) = header.sequence(region.seq()) else {
        return Ok(text(out));
    };
    let over = Over {
        number,
        start: region.start(),
        end: region.end(),
    };
    let mut record = Record::default();
    match index {
        Some(index) if stream.blocked() => {
            fits(index, stream.tell())?;
            'chunks: for (begin, stop) in index.chunks(number, over.start, over.end) {
                stream.seek(begin)?;
                while stream.tell() < stop {
                    if !record.shared(&mut stream)? {
                        break 'chunks;
                    }
                    match over.holds(&record)? {
                        Some(true) => record.write(&mut stream, &header, fields, &mut out)?,
                        Some(false) => record.pass(&mut stream)?,
                        None => break 'chunks,
                    }
                }
            }
        }
        _ => {
            while record.shared(&mut stream)? {
                match over.holds(&record)? {
                    Some(true) => record.write(&mut stream, &header, fields, &mut out)?,
                    _ => record.pass(&mut stream)?,
                }
            }
        }
    }
    Ok(text(out))
}

/// The VCF text of every record of a BCF, the header first, as `bcftools
/// view` prints them, cut to `fields`.
///
/// For a reader that draws a record away from where it lies, as a track of
/// structural calls draws the arc from a breakend's mate.
///
/// # Errors
///
/// A file that is not BCF or is damaged.
pub fn whole<R: Read + Seek>(reader: R, fields: Fields) -> Result<String, ReadError> {
    let mut stream = open(reader)?;
    let header = header(&mut stream)?;
    let mut out = header.written(fields).into_bytes();
    let mut record = Record::default();
    while record.shared(&mut stream)? {
        record.write(&mut stream, &header, fields, &mut out)?;
    }
    Ok(text(out))
}

/// How many records a BCF holds on each of its sequences, in the order the
/// header numbers them, leaving out the sequences with none: what an empty
/// window says the file does hold.
///
/// Counted from `index` where one is given and the file is BGZF, from the
/// count each sequence's pseudo-bin keeps, and no record is read. A sequence
/// the index has no bins for holds none, as `bcftools index` writes nothing
/// for a sequence the header names and no record is on, which most headers
/// of a whole genome have. Where the index keeps no count for a sequence it
/// has bins for, which the format leaves to the writer, and where no index
/// is given, every record is read as far as the number of its sequence.
///
/// # Errors
///
/// A file that is not BCF or is damaged, and an index [`window`] would not
/// read through: one that is not a CSI written for a BCF, or one that puts
/// the first record somewhere other than where the header ends, as the index
/// of another file does, whose counts are another file's.
pub fn counted<R: Read + Seek>(
    reader: R,
    index: Option<&Index>,
) -> Result<Vec<(String, usize)>, ReadError> {
    let mut stream = open(reader)?;
    let header = header(&mut stream)?;
    if let Some(index) = index.filter(|_| stream.blocked()) {
        fits(index, stream.tell())?;
        if let Some(counts) = indexed_counts(&header, index) {
            return Ok(counts);
        }
    }
    let mut counts = vec![0usize; header.sequences.len()];
    let mut record = Record::default();
    while record.shared(&mut stream)? {
        let number = record.sequence(&header)?;
        counts[number] += 1;
        record.pass(&mut stream)?;
    }
    Ok(header
        .sequences
        .iter()
        .zip(counts)
        .filter(|(_, count)| *count > 0)
        .filter_map(|(sequence, count)| Some((header.contigs[(*sequence)?].0.clone(), count)))
        .collect())
}

/// Whether `index` is one a BCF whose header ends at the virtual offset
/// `header_end` is read through, and why not where it is not: what [`window`]
/// reads through and what [`counted`] counts from are one index.
///
/// A CSI written for a BCF, which carries no columns, since a BCF's records
/// say where they are themselves. Its first record is where the header ends,
/// and an index of another file puts it somewhere else: read through, it
/// would hand over records of bytes that are not where it says, and counted
/// from, the records of the other file.
fn fits(index: &Index, header_end: u64) -> Result<(), ReadError> {
    if index.kind() != Kind::Csi || index.columns().is_some() {
        return Err(ReadError::whole(
            "the index is not a CSI written for a BCF, which is the one index a BCF has",
        ));
    }
    let first = (0..index.references())
        .filter_map(|reference| index.summary(reference))
        .map(|summary| summary.first)
        .min();
    if let Some(first) = first.filter(|first| *first != header_end) {
        return Err(ReadError::whole(format!(
            "the index puts the first record {}, and the file's header ends {}",
            super::tabix::at(first),
            super::tabix::at(header_end)
        )));
    }
    Ok(())
}

/// The records on each sequence the header names, as [`counted`] says them,
/// from the counts `index` keeps, in one pass over the header's numbers.
/// `None` where it keeps no count for a sequence it has bins for, or a count
/// past what a count here holds.
fn indexed_counts(header: &Header, index: &Index) -> Option<Vec<(String, usize)>> {
    let mut counts = Vec::new();
    for (number, sequence) in header.sequences.iter().enumerate() {
        let Some(at) = *sequence else {
            continue;
        };
        let rows = match index.summary(number) {
            Some(summary) => usize::try_from(summary.placed).ok()?,
            None if !index.binned(number) => 0,
            None => return None,
        };
        if rows > 0 {
            counts.push((header.contigs[at].0.clone(), rows));
        }
    }
    Some(counts)
}

/// Bytes written as text, which they are but where a value held bytes that
/// are not UTF-8, which are then replaced as a text file's would be refused.
fn text(out: Vec<u8>) -> String {
    String::from_utf8(out)
        .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned())
}

/// The bytes of a BCF in whichever of its three wrappings it came.
enum Stream<R> {
    /// BGZF, which an index points into.
    Blocked(Bgzf<R>),
    /// The bare stream.
    Plain(R),
    /// Plain gzip, inflated whole.
    Inflated(Cursor<Vec<u8>>),
}

impl<R: Read + Seek> Stream<R> {
    /// Whether the file is BGZF, which is the one wrapping an index points
    /// into.
    fn blocked(&self) -> bool {
        matches!(self, Stream::Blocked(_))
    }

    /// Goes to a virtual offset an index gives, in a file that is BGZF.
    fn seek(&mut self, virtual_offset: u64) -> Result<(), ReadError> {
        match self {
            Stream::Blocked(bgzf) => bgzf.seek(virtual_offset),
            _ => Ok(()),
        }
    }

    /// The virtual offset of the next byte, in a file that is BGZF.
    fn tell(&self) -> u64 {
        match self {
            Stream::Blocked(bgzf) => bgzf.tell(),
            _ => 0,
        }
    }

    /// Fills `buf` and says how much it filled: all of it, or less only where
    /// the file ends first.
    fn fill(&mut self, buf: &mut [u8]) -> Result<usize, ReadError> {
        match self {
            Stream::Blocked(bgzf) => bgzf.fill(buf),
            Stream::Plain(reader) => read_up_to(reader, buf),
            Stream::Inflated(cursor) => read_up_to(cursor, buf),
        }
    }

    /// Passes over `n` bytes without decoding them, read into `scratch` a
    /// piece at a time.
    ///
    /// They are read rather than sought past: a seek would cost a call to the
    /// system and the buffer a reader of a file on disk keeps, for each
    /// record, and BGZF has to be inflated to be passed through anyway. The
    /// pieces are at most 64 KiB, so a length from a damaged file asks for no
    /// more than that, and `scratch` is kept from one record to the next
    /// rather than made and zeroed afresh for each of a cohort's million.
    fn skip(&mut self, n: usize, scratch: &mut Vec<u8>) -> Result<(), ReadError> {
        let mut left = n;
        while left > 0 {
            let piece = left.min(1 << 16);
            if scratch.len() < piece {
                scratch.resize(piece, 0);
            }
            if self.fill(&mut scratch[..piece])? < piece {
                return Err(cut_short());
            }
            left -= piece;
        }
        Ok(())
    }

    /// Reads exactly `n` bytes into `into`, growing it as the bytes arrive,
    /// so a count from a damaged file asks for no more than the file holds.
    fn take(&mut self, n: usize, into: &mut Vec<u8>) -> Result<(), ReadError> {
        into.clear();
        while into.len() < n {
            let at = into.len();
            let piece = (n - at).min(1 << 16);
            into.resize(at + piece, 0);
            if self.fill(&mut into[at..])? < piece {
                return Err(cut_short());
            }
        }
        Ok(())
    }
}

/// A BCF cut short.
fn cut_short() -> ReadError {
    ReadError::whole("the BCF ends in the middle of a record")
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

/// The file at its start, unwrapped as its first bytes say: BGZF, plain gzip,
/// or the bare stream.
fn open<R: Read + Seek>(mut reader: R) -> Result<Stream<R>, ReadError> {
    let io = |error: std::io::Error| ReadError::whole(error.to_string());
    reader.seek(SeekFrom::Start(0)).map_err(io)?;
    // A gzip member's header and the extra field bgzip says a block's size in.
    let mut first = [0u8; 18];
    let got = read_up_to(&mut reader, &mut first)?;
    reader.seek(SeekFrom::Start(0)).map_err(io)?;
    let first = &first[..got];
    Ok(if bgzf::is_bgzf(first) {
        Stream::Blocked(Bgzf::new(reader))
    } else if gzip::is_gzip(first) {
        let mut all = Vec::new();
        reader.read_to_end(&mut all).map_err(io)?;
        Stream::Inflated(Cursor::new(gzip::decompress(&all)?))
    } else {
        Stream::Plain(reader)
    })
}

/// The header: the magic, the version, and the text that names everything a
/// record stores as a number.
fn header<R: Read + Seek>(stream: &mut Stream<R>) -> Result<Header, ReadError> {
    let mut magic = [0u8; 5];
    if stream.fill(&mut magic)? < magic.len() || &magic[..3] != b"BCF" {
        return Err(ReadError::whole(
            "not BCF: the file does not start with BCF's magic",
        ));
    }
    if magic[3..] != [2, 2] {
        return Err(ReadError::whole(format!(
            "the BCF is version {}.{}, and karyon reads 2.2, as htslib does",
            magic[3], magic[4]
        )));
    }
    let mut length = [0u8; 4];
    if stream.fill(&mut length)? < length.len() {
        return Err(ReadError::whole("the BCF's header is cut short"));
    }
    let mut text = Vec::new();
    stream
        .take(u32::from_le_bytes(length) as usize, &mut text)
        .map_err(|_| ReadError::whole("the BCF's header is cut short"))?;
    // The text ends with a nought, and may be padded with more.
    while text.last() == Some(&0) {
        text.pop();
    }
    parse_header(&String::from_utf8_lossy(&text))
}

/// The dictionaries a header's text numbers, and the text as `bcftools view`
/// prints it.
///
/// htslib writes `IDX=n` into every line that names a key in a BCF, and a
/// line that carries one is numbered by it: an INFO and a FORMAT key of one
/// name share a number, as `DP` usually is both. A line that carries none is
/// numbered next, after `PASS` at nought, as htslib numbers the lines of a
/// header that came from a VCF. A sequence is numbered among the `##contig`
/// lines alone, the same way.
fn parse_header(raw: &str) -> Result<Header, ReadError> {
    // A number past the length of the header is no number a line of it could
    // have been given, and allocating for it would let a damaged header ask
    // for gigabytes.
    let most = raw.len();
    let numbered = |idx: Option<&str>| -> Result<Option<usize>, ReadError> {
        idx.map(|idx| {
            idx.parse::<usize>()
                .ok()
                .filter(|idx| *idx <= most)
                .ok_or_else(|| {
                    ReadError::whole(format!(
                        "the BCF's header numbers a line {idx:?}, which no header numbers one"
                    ))
                })
        })
        .transpose()
    };
    let mut text = String::with_capacity(raw.len());
    let mut keys: Vec<Option<String>> = vec![Some("PASS".to_string())];
    let mut contigs: Vec<(String, Option<u64>)> = Vec::new();
    let mut sequences: Vec<Option<usize>> = Vec::new();
    let mut samples = None;
    let mut prefixed = false;
    for line in raw.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(rest) = line.strip_prefix("#CHROM") {
            samples = Some(
                rest.split('\t')
                    .skip(9)
                    .map(str::to_string)
                    .collect::<Vec<String>>(),
            );
            text.push_str(line);
            text.push('\n');
            continue;
        }
        if let Some(version) = line.strip_prefix("##fileformat=VCFv") {
            let mut parts = version.split('.').map(|part| part.parse::<u32>().ok());
            if let (Some(Some(major)), Some(Some(minor))) = (parts.next(), parts.next()) {
                prefixed = (major, minor) >= (4, 4);
            }
        }
        let Some((key, fields)) = structured(line) else {
            text.push_str(line);
            text.push('\n');
            continue;
        };
        let value = |name: &str| {
            fields.iter().find_map(|field| {
                field
                    .split_once('=')
                    .and_then(|(key, value)| (key == name).then_some(value))
            })
        };
        // The line as bcftools view prints it, its number taken out.
        text.push_str("##");
        text.push_str(key);
        text.push_str("=<");
        let kept: Vec<&str> = fields
            .iter()
            .copied()
            .filter(|field| !field.starts_with("IDX="))
            .collect();
        text.push_str(&kept.join(","));
        text.push_str(">\n");
        let Some(id) = value("ID") else {
            continue;
        };
        let idx = numbered(value("IDX"))?;
        match key {
            "FILTER" | "INFO" | "FORMAT" => {
                let at = match idx {
                    Some(at) => at,
                    None => match keys.iter().position(|held| held.as_deref() == Some(id)) {
                        Some(_) => continue,
                        None => keys.len(),
                    },
                };
                if keys.len() <= at {
                    keys.resize(at + 1, None);
                }
                keys[at] = Some(id.to_string());
            }
            "contig" => {
                let at = idx.unwrap_or(sequences.len());
                if sequences.len() <= at {
                    sequences.resize(at + 1, None);
                }
                let length = value("length").and_then(|length| length.parse().ok());
                sequences[at] = Some(contigs.len());
                contigs.push((id.to_string(), length));
            }
            _ => {}
        }
    }
    let Some(samples) = samples else {
        return Err(ReadError::whole(
            "the BCF's header has no #CHROM line, which every BCF's header ends with",
        ));
    };
    Ok(Header {
        text,
        contigs,
        samples,
        sequences,
        keys,
        prefixed,
    })
}

/// A header line of the form `##KEY=<ID=x,Number=1,...>`, as its key and its
/// fields in order, `ID=x`, a quoted value whole with its quotes, commas and
/// all.
fn structured(line: &str) -> Option<(&str, Vec<&str>)> {
    let rest = line.strip_prefix("##")?;
    let (key, body) = rest.split_once("=<")?;
    let body = body.strip_suffix('>')?;
    let mut fields = Vec::new();
    let mut start = 0;
    let mut quoted = false;
    let mut escaped = false;
    for (at, character) in body.char_indices() {
        match character {
            _ if escaped => escaped = false,
            '\\' if quoted => escaped = true,
            '"' => quoted = !quoted,
            ',' if !quoted => {
                fields.push(&body[start..at]);
                start = at + 1;
            }
            _ => {}
        }
    }
    fields.push(&body[start..]);
    Some((key, fields))
}

/// A window, as the number of its sequence and its 0-based half-open span.
struct Over {
    number: usize,
    start: u64,
    end: u64,
}

impl Over {
    /// Whether a record is over the window: `None` where it is past it on
    /// its sequence, where a sorted file has nothing more for the window.
    fn holds(&self, record: &Record) -> Result<Option<bool>, ReadError> {
        if usize::try_from(record.chrom).ok() != Some(self.number) {
            return Ok(Some(false));
        }
        let pos = u64::try_from(record.pos).unwrap_or(0);
        if pos >= self.end {
            return Ok(None);
        }
        let spans = u64::try_from(record.rlen)
            .unwrap_or(0)
            .max(record.reference()? as u64)
            .max(1);
        Ok(Some(pos.saturating_add(spans) > self.start))
    }
}

/// One record: its shared part, read, and its samples' part, read only where
/// it is written.
#[derive(Default)]
struct Record {
    shared: Vec<u8>,
    indiv: Vec<u8>,
    /// How many bytes the samples' part holds.
    l_indiv: usize,
    chrom: i32,
    pos: i32,
    rlen: i32,
}

/// What a record's fixed fields take: its sequence, position, length,
/// quality, and the two words its counts are packed in.
const FIXED: usize = 24;

impl Record {
    /// Reads the next record's shared part, and says whether there was one:
    /// false at the end of the file, which between records is where a BCF
    /// ends.
    fn shared<R: Read + Seek>(&mut self, stream: &mut Stream<R>) -> Result<bool, ReadError> {
        let mut lengths = [0u8; 8];
        match stream.fill(&mut lengths)? {
            0 => return Ok(false),
            8 => {}
            _ => return Err(cut_short()),
        }
        let l_shared = u32::from_le_bytes([lengths[0], lengths[1], lengths[2], lengths[3]]);
        let l_indiv = u32::from_le_bytes([lengths[4], lengths[5], lengths[6], lengths[7]]);
        let l_shared = l_shared as usize;
        if l_shared < FIXED {
            return Err(ReadError::whole(format!(
                "a BCF record says its sites take {l_shared} bytes, fewer than the {FIXED} its \
                 fixed fields take"
            )));
        }
        stream.take(l_shared, &mut self.shared)?;
        self.l_indiv = l_indiv as usize;
        let mut fixed = Bytes::new(&self.shared, false, "BCF record");
        self.chrom = fixed.u32()? as i32;
        self.pos = fixed.u32()? as i32;
        self.rlen = fixed.u32()? as i32;
        Ok(true)
    }

    /// Passes over the samples' part undecoded.
    fn pass<R: Read + Seek>(&mut self, stream: &mut Stream<R>) -> Result<(), ReadError> {
        stream.skip(self.l_indiv, &mut self.indiv)
    }

    /// The number of the record's sequence, which the header has to name.
    fn sequence(&self, header: &Header) -> Result<usize, ReadError> {
        header
            .named(self.chrom)
            .and_then(|_| usize::try_from(self.chrom).ok())
            .ok_or_else(|| unnamed(self.chrom))
    }

    /// How many bases the reference allele spells, as `bcftools view` would
    /// print it.
    fn reference(&self) -> Result<usize, ReadError> {
        let mut bytes = Bytes::new(&self.shared, false, "BCF record");
        bytes.take(16)?;
        let alleles = bytes.u32()? >> 16;
        bytes.take(4)?;
        // The ID, then the reference allele.
        Typed::read(&mut bytes)?;
        if alleles == 0 {
            return Ok(1);
        }
        let mut reference = Vec::new();
        Typed::read(&mut bytes)?.array(&mut reference)?;
        Ok(reference.len())
    }

    /// Writes the record as `bcftools view` prints its row, reading its
    /// samples' part where `fields` writes any of it.
    fn write<R: Read + Seek>(
        &mut self,
        stream: &mut Stream<R>,
        header: &Header,
        fields: Fields,
        out: &mut Vec<u8>,
    ) -> Result<(), ReadError> {
        let chrom = header
            .named(self.chrom)
            .ok_or_else(|| unnamed(self.chrom))?;
        let mut bytes = Bytes::new(&self.shared, false, "BCF record");
        bytes.take(12)?;
        let qual = bytes.u32()?;
        let counts = bytes.u32()?;
        let (infos, alleles) = (counts & 0xffff, counts >> 16);
        let packed = bytes.u32()?;
        let (samples, formats) = ((packed & 0x00ff_ffff) as usize, packed >> 24);
        if samples != header.samples.len() {
            return Err(ReadError::whole(format!(
                "a BCF record holds {samples} samples, and the header names {}",
                header.samples.len()
            )));
        }
        out.extend_from_slice(chrom.as_bytes());
        out.push(b'\t');
        integer(out, i64::from(self.pos) + 1);
        out.push(b'\t');
        Typed::read(&mut bytes)?.array(out)?;
        out.push(b'\t');
        if alleles == 0 {
            out.push(b'.');
        }
        for allele in 0..alleles {
            if allele == 1 {
                out.push(b'\t');
            } else if allele > 1 {
                out.push(b',');
            }
            Typed::read(&mut bytes)?.array(out)?;
        }
        if alleles <= 1 {
            out.extend_from_slice(b"\t.");
        }
        out.push(b'\t');
        if qual == FLOAT_MISSING {
            out.push(b'.');
        } else {
            float(out, f64::from(f32::from_bits(qual)));
        }
        out.push(b'\t');
        let filters = Typed::read(&mut bytes)?;
        if filters.count == 0 {
            out.push(b'.');
        }
        for at in 0..filters.count {
            if at > 0 {
                out.push(b';');
            }
            out.extend_from_slice(header.key(filters.integer(at)?)?.as_bytes());
        }
        out.push(b'\t');
        if infos == 0 {
            out.push(b'.');
        }
        for at in 0..infos {
            if at > 0 {
                out.push(b';');
            }
            let key = Typed::integer_of(&mut bytes)?;
            out.extend_from_slice(header.key(key)?.as_bytes());
            let value = Typed::read(&mut bytes)?;
            value.info(out)?;
        }
        if fields == Fields::Sites || samples == 0 {
            self.pass(stream)?;
            out.push(b'\n');
            return Ok(());
        }
        stream.take(self.l_indiv, &mut self.indiv)?;
        let mut indiv = Bytes::new(&self.indiv, false, "BCF record");
        let mut columns: Vec<(&str, Typed<'_>)> = Vec::with_capacity((formats as usize).min(64));
        for _ in 0..formats {
            let key = header.key(Typed::integer_of(&mut indiv)?)?;
            let (kind, count) = Typed::size(&mut indiv)?;
            let width = count
                .checked_mul(kind.width())
                .ok_or_else(|| short_of("a FORMAT field"))?;
            let all = width
                .checked_mul(samples)
                .ok_or_else(|| short_of("a FORMAT field"))?;
            let data = indiv.take(all)?;
            if fields == Fields::Genotypes && key != "GT" {
                continue;
            }
            columns.push((key, Typed { kind, count, data }));
        }
        if columns.is_empty() {
            for _ in 0..=samples {
                out.extend_from_slice(b"\t.");
            }
            out.push(b'\n');
            return Ok(());
        }
        for (at, (key, _)) in columns.iter().enumerate() {
            out.push(if at == 0 { b'\t' } else { b':' });
            out.extend_from_slice(key.as_bytes());
        }
        for sample in 0..samples {
            out.push(b'\t');
            for (at, (key, column)) in columns.iter().enumerate() {
                if at > 0 {
                    out.push(b':');
                }
                let width = column.count * column.kind.width();
                let one = Typed {
                    kind: column.kind,
                    count: column.count,
                    data: &column.data[sample * width..(sample + 1) * width],
                };
                if *key == "GT" {
                    one.genotype(out, header.prefixed)?;
                } else if one.count == 1 {
                    one.single(out)?;
                } else {
                    one.array(out)?;
                }
            }
        }
        out.push(b'\n');
        Ok(())
    }
}

/// A record on a sequence the header does not name.
fn unnamed(number: i32) -> ReadError {
    ReadError::whole(format!(
        "a BCF record is on sequence {number}, which the header's dictionary does not name"
    ))
}

/// A record's part that does not hold what it says.
fn short_of(what: &str) -> ReadError {
    ReadError::whole(format!("{what} of a BCF record does not fit the record"))
}

/// The bits a missing float is stored as.
const FLOAT_MISSING: u32 = 0x7F80_0001;
/// The bits that end a vector of floats shorter than its slot.
const FLOAT_END: u32 = 0x7F80_0002;

/// The type of a typed value, from the low four bits of its descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Type {
    /// No value, as a flag is stored.
    Missing,
    Int8,
    Int16,
    Int32,
    Float,
    Char,
}

impl Type {
    fn of(bits: u8) -> Result<Type, ReadError> {
        Ok(match bits {
            0 => Type::Missing,
            1 => Type::Int8,
            2 => Type::Int16,
            3 => Type::Int32,
            5 => Type::Float,
            7 => Type::Char,
            other => {
                return Err(ReadError::whole(format!(
                    "a BCF record holds a value of type {other}, which BCF does not have"
                )))
            }
        })
    }

    /// How many bytes one value of the type takes. One for no type, as
    /// htslib counts it.
    fn width(self) -> usize {
        match self {
            Type::Missing | Type::Int8 | Type::Char => 1,
            Type::Int16 => 2,
            Type::Int32 | Type::Float => 4,
        }
    }
}

/// A typed value: its type, how many values it holds, and their bytes.
struct Typed<'a> {
    kind: Type,
    count: usize,
    data: &'a [u8],
}

/// One value of a vector of numbers, as a number, a missing value or the end
/// of the vector.
enum Value {
    Int(i64),
    Float(f32),
    Missing,
    End,
}

impl<'a> Typed<'a> {
    /// A descriptor: the type in its low four bits and the count in its high
    /// four, where fifteen says the count follows as a typed integer.
    fn size(bytes: &mut Bytes<'a>) -> Result<(Type, usize), ReadError> {
        let descriptor = bytes.u8()?;
        let kind = Type::of(descriptor & 0x0f)?;
        let count = match descriptor >> 4 {
            15 => {
                let count = Typed::integer_of(bytes)?;
                usize::try_from(count).map_err(|_| {
                    ReadError::whole(format!("a BCF record holds a count of {count}"))
                })?
            }
            count => usize::from(count),
        };
        Ok((kind, count))
    }

    /// The next typed value, refused where its bytes run past the record's.
    fn read(bytes: &mut Bytes<'a>) -> Result<Typed<'a>, ReadError> {
        let (kind, count) = Typed::size(bytes)?;
        let length = count
            .checked_mul(kind.width())
            .filter(|length| *length <= bytes.left())
            .ok_or_else(|| short_of("a value"))?;
        Ok(Typed {
            kind,
            count,
            data: bytes.take(length)?,
        })
    }

    /// The next typed value that is one integer, as a key and an overflowing
    /// count are stored.
    fn integer_of(bytes: &mut Bytes<'a>) -> Result<i64, ReadError> {
        let typed = Typed::read(bytes)?;
        if typed.count != 1 || !matches!(typed.kind, Type::Int8 | Type::Int16 | Type::Int32) {
            return Err(ReadError::whole(
                "a BCF record stores a key or a count as something other than one integer",
            ));
        }
        typed.integer(0)
    }

    /// The `at`th value as an integer, as a FILTER is stored.
    fn integer(&self, at: usize) -> Result<i64, ReadError> {
        match self.value(at) {
            Value::Int(value) => Ok(value),
            _ => Err(ReadError::whole(
                "a BCF record stores a filter as something other than a number",
            )),
        }
    }

    /// The `at`th value, with the two values each type keeps for a missing
    /// value and for the end of a vector told apart, a float's by its bits.
    fn value(&self, at: usize) -> Value {
        let width = self.kind.width();
        let bytes = &self.data[at * width..(at + 1) * width];
        match self.kind {
            Type::Int8 => match bytes[0] as i8 {
                i8::MIN => Value::Missing,
                -127 => Value::End,
                value => Value::Int(i64::from(value)),
            },
            Type::Int16 => match i16::from_le_bytes([bytes[0], bytes[1]]) {
                i16::MIN => Value::Missing,
                -32767 => Value::End,
                value => Value::Int(i64::from(value)),
            },
            Type::Int32 => match i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) {
                i32::MIN => Value::Missing,
                -2147483647 => Value::End,
                value => Value::Int(i64::from(value)),
            },
            Type::Float => match u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) {
                FLOAT_MISSING => Value::Missing,
                FLOAT_END => Value::End,
                bits => Value::Float(f32::from_bits(bits)),
            },
            Type::Missing | Type::Char => Value::Int(i64::from(bytes[0])),
        }
    }

    /// The values as htslib's `bcf_fmt_array` prints them: `.` for none at
    /// all; a string up to its first nought; and numbers joined by commas, `.`
    /// for a missing one and stopping at the end of the vector.
    fn array(&self, out: &mut Vec<u8>) -> Result<(), ReadError> {
        if self.count == 0 {
            out.push(b'.');
            return Ok(());
        }
        match self.kind {
            Type::Missing => {
                return Err(ReadError::whole(
                    "a BCF record holds values of no type, which htslib refuses to print",
                ))
            }
            Type::Char => {
                let data = &self.data[..self.count];
                let end = data
                    .iter()
                    .position(|byte| *byte == 0)
                    .unwrap_or(data.len());
                out.extend_from_slice(&data[..end]);
                return Ok(());
            }
            _ => {}
        }
        for at in 0..self.count {
            let value = self.value(at);
            if matches!(value, Value::End) {
                break;
            }
            if at > 0 {
                out.push(b',');
            }
            match value {
                Value::Int(value) => integer(out, value),
                Value::Float(value) => float(out, f64::from(value)),
                Value::Missing | Value::End => out.push(b'.'),
            }
        }
        Ok(())
    }

    /// An INFO value after its key, as htslib's `vcf_format` prints one it
    /// has read from a BCF: nothing for a flag, which holds no value, and
    /// otherwise as [`Typed::array`] prints it, one value or several.
    fn info(&self, out: &mut Vec<u8>) -> Result<(), ReadError> {
        if self.count == 0 {
            return Ok(());
        }
        out.push(b'=');
        self.array(out)
    }

    /// One sample's value of a FORMAT field of one value, as htslib's
    /// `bcf_fmt_array1` prints it: a character, with the byte that stands for
    /// a missing one as `.`; a number, or `.` for a missing one; and nothing
    /// for the end of a vector.
    fn single(&self, out: &mut Vec<u8>) -> Result<(), ReadError> {
        match (self.kind, self.value(0)) {
            (Type::Missing, _) => {
                return Err(ReadError::whole(
                    "a BCF record holds values of no type, which htslib refuses to print",
                ))
            }
            (Type::Char, Value::Int(7)) => out.push(b'.'),
            (Type::Char, Value::Int(byte)) => out.push(byte as u8),
            (_, Value::End) => {}
            (_, Value::Missing) => out.push(b'.'),
            (_, Value::Int(value)) => integer(out, value),
            (_, Value::Float(value)) => float(out, f64::from(value)),
        }
        Ok(())
    }

    /// One sample's `GT` as htslib prints it: each allele as its number, or
    /// `.` for one nobody called, joined by `|` where the allele is phased
    /// with the one before it and `/` where not, and `.` for a call of none.
    /// An allele is stored as its number plus one, doubled, plus one where it
    /// is phased, so nought is an allele nobody called.
    ///
    /// Under VCF 4.4 the first allele says how it is phased too, and is
    /// printed with a `/` or `|` before it where the rest of the call does not
    /// imply it, as htslib's `bcf_format_gt_v2` decides: a call of several
    /// alleles is taken to be phased where every allele after the first is,
    /// so `/0|1` keeps its `/` and `|0|1` is printed `0|1`; a haploid one is
    /// taken to be phased, so `/0` keeps its `/` and `|0` is `0`, but for an
    /// allele nobody called, which is `.` unphased and `|.` phased.
    ///
    /// The missing value of the type is what htslib stores for a sample whose
    /// `GT` the text left out, as a sample written `4` under `DP:GT`, and
    /// before VCF 4.4 bcftools 1.24 prints it `.` where it is the whole call,
    /// alone or followed by the end of the vector. Anywhere else, and under
    /// VCF 4.4 everywhere, it is not looked for, and is printed as the number
    /// it is, `-65` for one byte, as bcftools prints it.
    fn genotype(&self, out: &mut Vec<u8>, prefixed: bool) -> Result<(), ReadError> {
        if self.kind == Type::Missing {
            out.push(b'.');
            return Ok(());
        }
        if !matches!(self.kind, Type::Int8 | Type::Int16 | Type::Int32) {
            return Err(ReadError::whole(
                "a BCF record stores a GT as something other than numbers",
            ));
        }
        let left_out = self.count > 0
            && matches!(self.value(0), Value::Missing)
            && (self.count == 1 || matches!(self.value(1), Value::End));
        if left_out && !prefixed {
            out.push(b'.');
            return Ok(());
        }
        let start = out.len();
        let mut first = 0i64;
        let mut unphased = false;
        let mut ploidy = 0;
        for at in 0..self.count {
            let value = match self.value(at) {
                Value::End => break,
                Value::Int(value) => value,
                // Anywhere but a whole call before VCF 4.4, the missing
                // value prints as the number it is, as htslib prints it.
                Value::Missing => match self.kind {
                    Type::Int8 => i64::from(i8::MIN),
                    Type::Int16 => i64::from(i16::MIN),
                    _ => i64::from(i32::MIN),
                },
                Value::Float(_) => 0,
            };
            if at == 0 {
                first = value;
            } else {
                out.push(if value & 1 == 1 { b'|' } else { b'/' });
                unphased |= value & 1 == 0;
            }
            match value >> 1 {
                0 => out.push(b'.'),
                allele => integer(out, allele - 1),
            }
            ploidy += 1;
        }
        if ploidy == 0 {
            out.push(b'.');
        }
        if prefixed {
            let prefix = if first & 1 == 1 {
                (ploidy > 1 && unphased || ploidy <= 1 && first >> 1 == 0).then_some(b'|')
            } else {
                (ploidy <= 1 && first != 0 || ploidy > 1 && !unphased).then_some(b'/')
            };
            if let Some(prefix) = prefix {
                out.insert(start, prefix);
            }
        }
        Ok(())
    }
}

/// An integer as `kputw` writes it.
fn integer(out: &mut Vec<u8>, value: i64) {
    let mut digits = [0u8; 20];
    let mut at = digits.len();
    let mut left = value.unsigned_abs();
    loop {
        at -= 1;
        digits[at] = b'0' + (left % 10) as u8;
        left /= 10;
        if left == 0 {
            break;
        }
    }
    if value < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(&digits[at..]);
}

/// A float as htslib's `kputd` writes it, which is `%g`, six significant
/// digits with the trailing noughts dropped, by a rounding of its own between
/// 0.0001 and 999,999 and by `printf` outside them.
///
/// Its own rounding scales the value to six digits and rounds those with
/// `rint`, a half to the even digit, where `printf` rounds the exact binary
/// value: 123,456.5 is printed 123456, not 123457. A value printed the other
/// way is another number to every reader of the text, so this is htslib 1.23's
/// arithmetic, step for step.
fn float(out: &mut Vec<u8>, value: f64) {
    if value == 0.0 {
        out.extend_from_slice(if value.is_sign_negative() {
            b"-0"
        } else {
            b"0"
        });
        return;
    }
    let mut d = value;
    if d < 0.0 {
        out.push(b'-');
        d = -d;
    }
    if !(0.0001..=999_999.0).contains(&d) {
        out.extend_from_slice(printf_g(d).as_bytes());
        return;
    }
    // The power of ten that makes six digits of it, and where the point goes,
    // counted as htslib counts it from the end of a buffer of twenty.
    let (scale, shift): (f64, usize) = if d < 0.001 {
        (1e9, 1)
    } else if d < 0.01 {
        (1e8, 2)
    } else if d < 0.1 {
        (1e7, 3)
    } else if d < 1.0 {
        (1e6, 4)
    } else if d < 10.0 {
        (1e5, 5)
    } else if d < 100.0 {
        (1e4, 6)
    } else if d < 1_000.0 {
        (1e3, 7)
    } else if d < 10_000.0 {
        (1e2, 8)
    } else if d < 100_000.0 {
        (1e1, 9)
    } else {
        (1.0, 10)
    };
    // Six digits, or seven where the rounding carried, as 0.0099999999 does.
    let digits = format!("{:06}", rint(d * scale) as u32).into_bytes();
    let placed = shift + digits.len();
    let (mut written, mut end) = if placed <= 10 {
        let mut written = b"0.".to_vec();
        written.extend(std::iter::repeat(b'0').take(10 - placed));
        let end = written.len() + 5;
        written.extend_from_slice(&digits);
        (written, end)
    } else {
        let whole = placed - 10;
        let mut written = digits[..whole].to_vec();
        written.push(b'.');
        written.extend_from_slice(&digits[whole..]);
        (written, 6)
    };
    // Back over the noughts, and past the last digit kept, or onto a point
    // left with nothing after it, which goes.
    while written[end] == b'0' && end > 0 {
        end -= 1;
    }
    if written[end] != b'.' {
        end += 1;
    }
    written.truncate(end);
    out.extend_from_slice(&written);
}

/// C's `rint` in its default mode, the nearest whole number with a half to
/// the even one, which `f64::round_ties_even` is from Rust 1.77 and this crate
/// builds with 1.74.
fn rint(value: f64) -> f64 {
    let floor = value.floor();
    match (value - floor).partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => floor,
        Some(std::cmp::Ordering::Greater) => floor + 1.0,
        _ if floor % 2.0 == 0.0 => floor,
        _ => floor + 1.0,
    }
}

/// `printf("%g")`, for the values `kputd` hands to it: six significant
/// digits, a half to the even digit, as an exponent where it is below -4 or
/// six and above, the trailing noughts dropped. macOS's libc keeps one where
/// it rounds a half down onto a nought, as `4.81980e+06` for 4,819,805, which
/// the C standard, glibc and this do not.
fn printf_g(d: f64) -> String {
    if d.is_nan() {
        return "nan".to_string();
    }
    if d.is_infinite() {
        return "inf".to_string();
    }
    let exponential = format!("{d:.5e}");
    let (mantissa, exponent) = exponential.split_once('e').unwrap_or((&exponential, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let trimmed = |digits: String| -> String {
        if digits.contains('.') {
            digits
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string()
        } else {
            digits
        }
    };
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            trimmed(mantissa.to_string()),
            exponent.unsigned_abs()
        )
    } else {
        let places = (5 - exponent).max(0) as usize;
        trimmed(format!("{d:.places$}"))
    }
}

/// The BCF files the tests here and the command line's read, written by
/// bcftools 1.24 from the VCF beside each, with what bcftools prints for them;
/// `fixtures/bcf/make.py` writes them all again.
#[cfg(test)]
pub(crate) mod fixture {
    /// Seven records of three samples over three sequences, one with no
    /// length, holding a value of every type and every way of being missing.
    pub(crate) const TINY: &[u8] = include_bytes!("fixtures/bcf/tiny.bcf");
    /// Its CSI, as `bcftools index` writes it.
    pub(crate) const TINY_CSI: &[u8] = include_bytes!("fixtures/bcf/tiny.bcf.csi");
    /// The same records as `bcftools view -Ou` writes them, BGZF of stored
    /// blocks.
    pub(crate) const TINY_U: &[u8] = include_bytes!("fixtures/bcf/tiny.u.bcf");
    /// The same records with no BGZF around them, as gzip leaves them.
    pub(crate) const TINY_RAW: &[u8] = include_bytes!("fixtures/bcf/tiny.raw.bcf");
    /// `bcftools view --no-version tiny.bcf`.
    pub(crate) const TINY_VIEW: &str = include_str!("fixtures/bcf/tiny.bcf.vcf");
    /// `bcftools view --no-version -G tiny.bcf`.
    pub(crate) const TINY_SITES: &str = include_str!("fixtures/bcf/tiny.bcf.sites.vcf");
    /// `bcftools annotate --no-version -x ^FORMAT/GT tiny.bcf`, which also
    /// takes the other FORMAT lines out of the header.
    pub(crate) const TINY_GT: &str = include_str!("fixtures/bcf/tiny.bcf.gt.vcf");
    /// The calls of `indexed/cohort.vcf.gz`, in blocks of a few kilobytes cut
    /// anywhere, so a record lies across two blocks.
    pub(crate) const COHORT: &[u8] = include_bytes!("fixtures/bcf/cohort.bcf");
    /// Its CSI.
    pub(crate) const COHORT_CSI: &[u8] = include_bytes!("fixtures/bcf/cohort.bcf.csi");
    /// `bcftools view --no-version cohort.bcf`.
    pub(crate) const COHORT_VIEW: &str = include_str!("fixtures/bcf/cohort.bcf.vcf");
    /// Floats of every size, each as bcftools prints it.
    pub(crate) const FLOATS: &[u8] = include_bytes!("fixtures/bcf/floats.bcf");
    /// `bcftools view --no-version floats.bcf`.
    pub(crate) const FLOATS_VIEW: &str = include_str!("fixtures/bcf/floats.bcf.vcf");
    /// Structural calls, breakends whose mates lie far off among them.
    pub(crate) const SV: &[u8] = include_bytes!("fixtures/bcf/sv.bcf");
    /// Its CSI.
    pub(crate) const SV_CSI: &[u8] = include_bytes!("fixtures/bcf/sv.bcf.csi");
    /// `bcftools view --no-version sv.bcf`.
    pub(crate) const SV_VIEW: &str = include_str!("fixtures/bcf/sv.bcf.vcf");
    /// Calls under VCF 4.4, whose first allele may say how it is phased.
    pub(crate) const PHASED: &[u8] = include_bytes!("fixtures/bcf/phased.bcf");
    /// `bcftools view --no-version phased.bcf`.
    pub(crate) const PHASED_VIEW: &str = include_str!("fixtures/bcf/phased.bcf.vcf");
}

#[cfg(test)]
mod tests {
    use super::fixture::*;
    use super::*;
    use crate::read::index;

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    fn all(bytes: &[u8], fields: Fields) -> String {
        whole(Cursor::new(bytes), fields).unwrap()
    }

    /// The rows of VCF text, its header left out.
    fn rows(text: &str) -> Vec<&str> {
        text.lines().filter(|line| !line.starts_with('#')).collect()
    }

    /// Every record of each file is the row bcftools prints for it, byte for
    /// byte: a value of every type, a string longer than a descriptor counts,
    /// a vector shorter than its slot, missing values of every width and the
    /// floats htslib rounds its own way. The header is bcftools', its `IDX`
    /// taken out, and the three wrappings of one file read alike.
    #[test]
    fn a_bcf_reads_as_bcftools_view_prints_it() {
        for (bytes, view) in [
            (TINY, TINY_VIEW),
            (TINY_U, TINY_VIEW),
            (TINY_RAW, TINY_VIEW),
            (COHORT, COHORT_VIEW),
            (FLOATS, FLOATS_VIEW),
            (PHASED, PHASED_VIEW),
        ] {
            assert_eq!(all(bytes, Fields::All), view);
        }
        assert!(!TINY_VIEW.contains("IDX="));
        // The same file with the gzip around its BGZF made plain gzip.
        let plain = gzip_of(TINY_RAW);
        assert_eq!(all(&plain, Fields::All), TINY_VIEW);
    }

    /// Sites alone are the eight columns `bcftools view -G` prints, and the
    /// samples' bytes are never decoded.
    #[test]
    fn sites_alone_are_what_bcftools_view_minus_g_prints() {
        assert_eq!(all(TINY, Fields::Sites), TINY_SITES);
        assert_eq!(all(TINY_RAW, Fields::Sites), TINY_SITES);
    }

    /// `GT` alone is what `bcftools annotate -x ^FORMAT/GT` keeps of each
    /// record, a record with none of it `.` for its FORMAT and each sample,
    /// and the header is left whole, which the readers do not mind.
    #[test]
    fn genotypes_alone_are_what_bcftools_keeps_of_gt() {
        let read = all(TINY, Fields::Genotypes);
        assert_eq!(rows(&read), rows(TINY_GT));
        assert_eq!(read.lines().last(), TINY_GT.lines().last());
        assert_eq!(
            read.lines().filter(|line| line.starts_with("##")).count(),
            TINY_VIEW
                .lines()
                .filter(|line| line.starts_with("##"))
                .count()
        );
        assert!(read.contains("chrM\t100\tv7\tA\tT\t.\tPASS\t.\t.\t.\t.\t."));
    }

    /// The header names what records store as numbers: the sequences with
    /// their lengths, where they have one, the samples, and a key INFO and
    /// FORMAT share, DP, under the one number its IDX gives it.
    #[test]
    fn the_header_names_what_the_records_number() {
        let header = header_of(Cursor::new(TINY)).unwrap();
        assert_eq!(
            header.contigs,
            [
                ("chr1".to_string(), Some(2_000_000)),
                ("chr2".to_string(), Some(90_000)),
                ("chrM".to_string(), None),
            ]
        );
        assert_eq!(header.samples, ["A", "B", "C"]);
        assert_eq!(header.sequence("chr2"), Some(1));
        assert_eq!(header.sequence("chr9"), None);
        let dp = header
            .keys
            .iter()
            .position(|key| key.as_deref() == Some("DP"));
        assert_eq!(dp, Some(3));
        assert_eq!(
            header
                .keys
                .iter()
                .filter(|key| key.as_deref() == Some("DP"))
                .count(),
            1
        );
        assert!(TINY_VIEW.contains("FORMAT\tA\tB\tC\n"));
        assert!(TINY_VIEW.contains("0/1:7:3,4:12.5:PASS"));
    }

    /// A header with no `IDX` is numbered in its order, `PASS` first and a
    /// key two lines share once, and one with `IDX` by it, gaps and all.
    #[test]
    fn the_dictionary_follows_idx_or_else_the_order_of_the_lines() {
        let plain = "##fileformat=VCFv4.2\n\
                     ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"a, b\">\n\
                     ##FILTER=<ID=q10,Description=\"x\">\n\
                     ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"y\">\n\
                     ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"z\">\n\
                     ##contig=<ID=b,length=5>\n##contig=<ID=a>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header = parse_header(plain).unwrap();
        let keys: Vec<Option<&str>> = header.keys.iter().map(Option::as_deref).collect();
        assert_eq!(keys, [Some("PASS"), Some("DP"), Some("q10"), Some("GT")]);
        assert_eq!(header.sequence("a"), Some(1));
        assert_eq!(header.text, plain);
        let numbered = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"q\",IDX=7>\n\
                        ##contig=<ID=chr9,IDX=4>\n\
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header = parse_header(numbered).unwrap();
        assert_eq!(header.keys[7].as_deref(), Some("AF"));
        assert_eq!(header.key(7).unwrap(), "AF");
        assert!(header.key(6).is_err() && header.key(-1).is_err());
        assert_eq!(header.sequence("chr9"), Some(4));
        assert_eq!(header.named(4), Some("chr9"));
        assert_eq!(header.named(0), None);
        assert!(!header.text.contains("IDX"));
        // A number no line of a header this long could be given.
        let far = "##contig=<ID=x,IDX=99999999>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        assert!(parse_header(far).is_err());
        assert!(parse_header("##fileformat=VCFv4.2\n").is_err());
    }

    /// Whether a record lies over a window, by the bases its reference
    /// allele spells or its END reaches, whichever is further.
    fn spans(row: &str) -> (String, u64, u64) {
        let fields: Vec<&str> = row.split('\t').collect();
        let pos: u64 = fields[1].parse::<u64>().unwrap() - 1;
        let end = fields[7]
            .split(';')
            .find_map(|pair| pair.strip_prefix("END="))
            .map(|end| end.parse::<u64>().unwrap())
            .unwrap_or(0);
        let reach = (fields[3].len() as u64).max(end.saturating_sub(pos)).max(1);
        (fields[0].to_string(), pos, pos + reach)
    }

    /// Through its CSI and read whole, a window holds the records over it and
    /// no others, every one as bcftools prints it, over a grid of windows of
    /// every size across the cohort's three sequences, blocks cut anywhere.
    #[test]
    fn a_window_holds_the_records_over_it_through_the_csi_or_without_one() {
        let csi = index::parse(COHORT_CSI).unwrap();
        assert_eq!(csi.kind(), Kind::Csi);
        let lengths = [
            ("chr1", 3_000_000u64),
            ("chr2", 1_000_000),
            ("chr3", 200_000),
        ];
        let mut over: Vec<Region> = Vec::new();
        for (sequence, length) in lengths {
            for width in [1u64, 150, 4_000, 60_000, 700_000] {
                let mut start = 0;
                while start < length {
                    over.push(Region::new(sequence, start, (start + width).min(length)).unwrap());
                    start += length / 29 + width / 3 + 1;
                }
            }
        }
        // And the base before every ninth record, its first base, and its
        // last, where a window's edge meets it.
        for row in rows(COHORT_VIEW).into_iter().step_by(9) {
            let (sequence, from, to) = spans(row);
            for (start, end) in [
                (from.saturating_sub(1), from),
                (from, from + 1),
                (to - 1, to),
            ] {
                if start < end {
                    over.push(Region::new(&sequence, start, end).unwrap());
                }
            }
        }
        let mut windows = 0;
        let mut held = 0;
        for over in over {
            let expected: Vec<&str> = rows(COHORT_VIEW)
                .into_iter()
                .filter(|row| {
                    let (name, from, to) = spans(row);
                    name == over.seq() && from < over.end() && to > over.start()
                })
                .collect();
            for index in [Some(&csi), None] {
                let read = window(Cursor::new(COHORT), index, &over, Fields::All).unwrap();
                assert_eq!(rows(&read), expected, "{over} through {}", index.is_some());
                assert!(read.starts_with("##fileformat=VCFv4.2\n"));
            }
            held += expected.len();
            windows += 1;
        }
        assert!(
            windows > 400 && held > 1_000,
            "{windows} windows, {held} records"
        );
        // The END of a deletion reaches past what its N spells.
        let read = window(
            Cursor::new(TINY),
            Some(&index::parse(TINY_CSI).unwrap()),
            &region("chr1:50,400-50,401"),
            Fields::Sites,
        )
        .unwrap();
        assert_eq!(
            rows(&read),
            ["chr1\t50000\tsv1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=50500"]
        );
    }

    /// A record reaches as far as its reference allele spells, whatever its
    /// stored length says, as the readers of VCF keep it: htslib stores the
    /// longer of the two, and a file another writer wrote with a length
    /// shorter than its allele still draws the allele.
    #[test]
    fn a_record_reaches_as_far_as_its_reference_allele_spells() {
        // chr1:16380, whose REF spells 20 bases, said to span one.
        let mut short = TINY_RAW.to_vec();
        let fixed: Vec<u8> = [0i32, 16379, 20]
            .iter()
            .flat_map(|word| word.to_le_bytes())
            .collect();
        let at = short.windows(12).position(|words| words == fixed).unwrap();
        short[at + 8..at + 12].copy_from_slice(&1i32.to_le_bytes());
        let read = window(
            Cursor::new(&short),
            None,
            &region("chr1:16,395-16,395"),
            Fields::Sites,
        )
        .unwrap();
        assert_eq!(rows(&read).len(), 1, "{read}");
        assert!(rows(&read)[0].starts_with("chr1\t16380\t"), "{read}");
    }

    /// The index is what is read: a block of the file the window's records
    /// are not in, damaged, stops a read of the whole file and not a read of
    /// the window through the index.
    #[test]
    fn a_window_through_the_csi_reads_only_the_blocks_over_it() {
        let csi = index::parse(COHORT_CSI).unwrap();
        // The last block that holds records, which are chr3's.
        let mut starts = Vec::new();
        let mut at = 0;
        while at + 18 <= COHORT.len() {
            starts.push(at);
            at += usize::from(u16::from_le_bytes([COHORT[at + 16], COHORT[at + 17]])) + 1;
        }
        let last = starts[starts.len() - 2];
        let mut bent = COHORT.to_vec();
        bent[last + 30] ^= 0xff;
        let early = region("chr1:1-200,000");
        let through = window(Cursor::new(&bent), Some(&csi), &early, Fields::Sites).unwrap();
        let sound = window(Cursor::new(COHORT), None, &early, Fields::Sites).unwrap();
        assert_eq!(through, sound);
        assert!(!rows(&through).is_empty());
        assert!(window(Cursor::new(&bent), None, &early, Fields::Sites).is_err());
    }

    /// The index of another file is refused by where it says the first
    /// record is, rather than read through into the middle of records.
    #[test]
    fn the_index_of_another_bcf_is_refused() {
        let other = index::parse(TINY_CSI).unwrap();
        let error = window(
            Cursor::new(COHORT),
            Some(&other),
            &region("chr1:1-100"),
            Fields::All,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("the index puts the first record"),
            "{error}"
        );
    }

    #[test]
    fn a_window_on_a_sequence_the_header_does_not_name_holds_no_record() {
        let read = window(Cursor::new(TINY), None, &region("chr9:1-100"), Fields::All).unwrap();
        assert_eq!(
            read,
            TINY_VIEW
                .lines()
                .take_while(|line| line.starts_with('#'))
                .map(|line| format!("{line}\n"))
                .collect::<String>()
        );
    }

    /// The records on each sequence are the ones bcftools prints, counted from
    /// the file and from its CSI alike. A sequence the header names and no
    /// record is on, which the CSI gives no bins, holds none, and the rest are
    /// counted from the index without a record read: a block of records
    /// damaged stops the count of every record and not the index's. An index
    /// of another file, or one that is not a BCF's, is refused, as a window
    /// refuses it, where its counts were the other file's.
    #[test]
    fn the_records_on_each_sequence_are_counted() {
        let csi = |bytes: &[u8]| index::parse(bytes).unwrap();
        for (bytes, index, expected) in [
            (
                COHORT,
                csi(COHORT_CSI),
                [("chr1", 300), ("chr2", 120), ("chr3", 30)],
            ),
            (TINY, csi(TINY_CSI), [("chr1", 5), ("chr2", 3), ("chrM", 1)]),
        ] {
            let expected: Vec<(String, usize)> = expected
                .iter()
                .map(|(name, count)| (name.to_string(), *count))
                .collect();
            assert_eq!(counted(Cursor::new(bytes), None).unwrap(), expected);
            assert_eq!(counted(Cursor::new(bytes), Some(&index)).unwrap(), expected);
        }
        let sv = csi(SV_CSI);
        // sv.bcf's header names scaffold_9 between chr1 and chr2.
        let header = header_of(Cursor::new(SV)).unwrap();
        let scaffold = header.sequence("scaffold_9").unwrap();
        assert!(scaffold < sv.references() && !sv.binned(scaffold));
        assert_eq!(sv.summary(scaffold), None);
        let expected = [("chr1".to_string(), 6), ("chr2".to_string(), 3)];
        assert_eq!(counted(Cursor::new(SV), None).unwrap(), expected);
        let mut starts = Vec::new();
        let mut at = 0;
        while at + 18 <= SV.len() {
            starts.push(at);
            at += usize::from(u16::from_le_bytes([SV[at + 16], SV[at + 17]])) + 1;
        }
        let mut bent = SV.to_vec();
        bent[starts[starts.len() - 2] + 30] ^= 0xff;
        assert!(counted(Cursor::new(&bent), None).is_err());
        assert_eq!(counted(Cursor::new(&bent), Some(&sv)).unwrap(), expected);
        let error = counted(Cursor::new(COHORT), Some(&csi(TINY_CSI))).unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("the index puts the first record"),
            "{error}"
        );
        let tbi = csi(&crate::read::index::fixture::ROWS_TBI);
        assert!(counted(Cursor::new(COHORT), Some(&tbi)).is_err());
    }

    /// Each width's missing value and end of a vector, a float's told by its
    /// bits, are told apart: `.` for one, the end of the values for the other.
    #[test]
    fn missing_and_end_of_vector_are_told_apart_in_every_width() {
        let typed = |kind: Type, data: &[u8]| -> String {
            let width = kind.width();
            let mut out = Vec::new();
            Typed {
                kind,
                count: data.len() / width,
                data,
            }
            .array(&mut out)
            .unwrap();
            String::from_utf8(out).unwrap()
        };
        assert_eq!(typed(Type::Int8, &[3, 0x80, 0x81, 9]), "3,.");
        let int16: Vec<u8> = [300i16, i16::MIN, -32767]
            .iter()
            .flat_map(|value| value.to_le_bytes())
            .collect();
        assert_eq!(typed(Type::Int16, &int16), "300,.");
        let int32: Vec<u8> = [i32::MIN, -70000, i32::MIN + 1]
            .iter()
            .flat_map(|value| value.to_le_bytes())
            .collect();
        assert_eq!(typed(Type::Int32, &int32), ".,-70000");
        let floats: Vec<u8> = [1.5f32.to_bits(), FLOAT_MISSING, FLOAT_END, 2.0f32.to_bits()]
            .iter()
            .flat_map(|bits| bits.to_le_bytes())
            .collect();
        assert_eq!(typed(Type::Float, &floats), "1.5,.");
        // A NaN that is neither is a number nobody wrote, and printed.
        let nan = 0x7FC0_0000u32.to_le_bytes();
        assert_eq!(typed(Type::Float, &nan), "nan");
        assert_eq!(typed(Type::Char, b"ab\0\0"), "ab");
        // One character, the byte for a missing one is `.`, and a number at
        // the end of its vector is nothing.
        let single = |kind: Type, data: &[u8]| -> String {
            let mut out = Vec::new();
            Typed {
                kind,
                count: 1,
                data,
            }
            .single(&mut out)
            .unwrap();
            String::from_utf8(out).unwrap()
        };
        assert_eq!(single(Type::Char, &[7]), ".");
        assert_eq!(single(Type::Char, b"x"), "x");
        assert_eq!(single(Type::Int8, &[0x81]), "");
        assert_eq!(single(Type::Int8, &[0x80]), ".");
        // A GT: a haploid call padded to a diploid's slot, and a call of none.
        let gt = |data: &[u8], prefixed: bool| -> String {
            let mut out = Vec::new();
            Typed {
                kind: Type::Int8,
                count: data.len(),
                data,
            }
            .genotype(&mut out, prefixed)
            .unwrap();
            String::from_utf8(out).unwrap()
        };
        assert_eq!(gt(&[0x02, 0x81], false), "0");
        assert_eq!(gt(&[0x00, 0x00], false), "./.");
        assert_eq!(gt(&[0x81, 0x81], false), ".");
        assert_eq!(gt(&[0x02, 0x05], false), "0|1");
        assert_eq!(gt(&[0x02, 0x05], true), "/0|1");
        assert_eq!(gt(&[0x03, 0x04], true), "|0/1");
        assert_eq!(gt(&[0x03, 0x81], true), "0");
        assert_eq!(gt(&[0x02, 0x81], true), "/0");
        assert_eq!(gt(&[0x01, 0x81], true), "|.");
        assert_eq!(gt(&[0x00, 0x81], true), ".");
        // A GT left out of a sample is stored as the missing value, which
        // bcftools 1.24 prints `.` as a whole call before VCF 4.4, and as the
        // number it is anywhere else, as measured on files patched to hold
        // each of these.
        for (data, before, under) in [
            (&[0x80][..], ".", "/-65"),
            (&[0x80, 0x81], ".", "/-65"),
            (&[0x80, 0x81, 0x81], ".", "/-65"),
            (&[0x80, 0x80], "-65/-65", "-65/-65"),
            (&[0x80, 0x80, 0x81], "-65/-65", "-65/-65"),
            (&[0x80, 0x04], "-65/1", "-65/1"),
            (&[0x02, 0x80], "0/-65", "0/-65"),
        ] {
            assert_eq!(gt(data, false), before, "{data:02x?}");
            assert_eq!(gt(data, true), under, "{data:02x?}");
        }
        let mut out = Vec::new();
        let wide = [i16::MIN.to_le_bytes(), (-32767i16).to_le_bytes()].concat();
        Typed {
            kind: Type::Int16,
            count: 2,
            data: &wide,
        }
        .genotype(&mut out, false)
        .unwrap();
        assert_eq!(out, b".");
        let mut out = Vec::new();
        Typed {
            kind: Type::Int8,
            count: 0,
            data: &[],
        }
        .genotype(&mut out, false)
        .unwrap();
        assert_eq!(out, b".");
    }

    /// Floats as htslib prints them, which is not always the nearest six
    /// digits: every one of the fixture's, chosen and drawn, is the text
    /// bcftools printed for it, and a few are checked here by hand.
    #[test]
    fn a_float_is_printed_as_htslib_prints_it() {
        let printed = |value: f32| {
            let mut out = Vec::new();
            float(&mut out, f64::from(value));
            String::from_utf8(out).unwrap()
        };
        for (value, text) in [
            (29.5, "29.5"),
            (1e6, "1e+06"),
            (0.333_333, "0.333333"),
            (0.0001, "0.0001"),
            (1.5e-7, "1.5e-07"),
            (123_456.5, "123456"),
            (238_148.5, "238148"),
            (0.009_999_999, "0.01"),
            (9.999_999, "10"),
            (999_999.0, "999999"),
            (-0.5, "-0.5"),
            (-0.0, "-0"),
            (3.0, "3"),
            (f32::INFINITY, "inf"),
        ] {
            assert_eq!(printed(value), text, "{value}");
        }
        let values: Vec<&str> = rows(FLOATS_VIEW)
            .iter()
            .flat_map(|row| row.split('\t').nth(7).unwrap()[2..].split(','))
            .collect();
        assert!(values.len() > 600, "{}", values.len());
    }

    #[test]
    fn a_file_that_is_not_bcf_2_2_is_refused_saying_what_it_is() {
        let error = header_of(Cursor::new(&b"##fileformat=VCFv4.2\n"[..])).unwrap_err();
        assert!(error.to_string().contains("not BCF"), "{error}");
        let mut older = TINY_RAW.to_vec();
        older[4] = 1;
        let error = header_of(Cursor::new(older)).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the BCF is version 2.1, and karyon reads 2.2, as htslib does"
        );
        let tbi = index::parse(&crate::read::index::fixture::ROWS_TBI).unwrap();
        let error = window(
            Cursor::new(TINY),
            Some(&tbi),
            &region("chr1:1-100"),
            Fields::Sites,
        )
        .unwrap_err();
        assert!(error.to_string().contains("not a CSI"), "{error}");
    }

    /// Cut anywhere, or with any byte of it changed, a BCF is an error or
    /// some other text and never a panic, and a count from a damaged file
    /// asks for no more than the file holds. The bare stream has no CRC to
    /// stop a changed byte before it is decoded, so every byte of it is bent
    /// three ways and every record decoded.
    #[test]
    fn a_damaged_bcf_or_csi_is_an_error_and_never_a_panic() {
        let over = region("chr1:1-2,000,000");
        for cut in 0..TINY_RAW.len() {
            let _ = whole(Cursor::new(&TINY_RAW[..cut]), Fields::All);
        }
        for at in 0..TINY_RAW.len() {
            for bend in [0xff, 0x00, 0x55] {
                let mut bent = TINY_RAW.to_vec();
                bent[at] ^= bend;
                for fields in [Fields::All, Fields::Genotypes, Fields::Sites] {
                    let _ = whole(Cursor::new(&bent), fields);
                    let _ = window(Cursor::new(&bent), None, &over, fields);
                }
                let _ = counted(Cursor::new(&bent), None);
            }
        }
        for cut in 0..TINY.len() {
            let _ = whole(Cursor::new(&TINY[..cut]), Fields::All);
        }
        let csi = index::parse(TINY_CSI).unwrap();
        for at in 0..TINY.len() {
            let mut bent = TINY.to_vec();
            bent[at] ^= 0x55;
            let _ = window(Cursor::new(&bent), Some(&csi), &over, Fields::All);
            let _ = counted(Cursor::new(&bent), Some(&csi));
        }
        for cut in 0..TINY_CSI.len() {
            if let Ok(index) = index::parse(&TINY_CSI[..cut]) {
                let _ = window(Cursor::new(TINY), Some(&index), &over, Fields::All);
                let _ = counted(Cursor::new(TINY), Some(&index));
            }
        }
        // Lengths at the top of their words: a header, a record's sites and
        // its samples, each said to be four gigabytes.
        let header_at = 5;
        let mut long = TINY_RAW.to_vec();
        long[header_at..header_at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(whole(Cursor::new(&long), Fields::All).is_err());
        let text = u32::from_le_bytes(TINY_RAW[5..9].try_into().unwrap()) as usize;
        for word in [0, 4] {
            let mut long = TINY_RAW.to_vec();
            let at = 9 + text + word;
            long[at..at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
            let error = whole(Cursor::new(&long), Fields::All).unwrap_err();
            assert!(error.to_string().contains("ends in the middle"), "{error}");
        }
    }

    fn gzip_of(data: &[u8]) -> Vec<u8> {
        // One stored block in a gzip member with no extra field.
        let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
        for chunk in data.chunks(0xffff) {
            let last = chunk.as_ptr() == data.chunks(0xffff).last().unwrap().as_ptr();
            out.push(u8::from(last));
            let length = chunk.len() as u16;
            out.extend_from_slice(&length.to_le_bytes());
            out.extend_from_slice(&(!length).to_le_bytes());
            out.extend_from_slice(chunk);
        }
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
