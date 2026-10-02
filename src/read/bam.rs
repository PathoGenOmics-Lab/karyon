//! BAM, read a window at a time through its index.
//!
//! A BAM is the reads of a SAM file in binary, compressed as BGZF, which is a
//! series of gzip members of at most 64 KiB each. Its `.bai` index says which
//! members hold the reads over any window, so a figure of one gene reads the
//! few blocks that gene is in rather than the whole file: a human genome's
//! BAM runs to tens of gigabytes, and the window is a few thousand bases of
//! it. Without an index the file is read from its start, and a file sorted by
//! position is left as soon as the reads pass the window.
//!
//! Nothing here opens a file. [`window`] takes anything that reads and seeks,
//! which is a file for the command line and a buffer for a test, and
//! [`index`] takes the index's bytes. A `.csi`, which `samtools index -c`
//! writes for a sequence longer than a BAI has room for, is read by
//! [`index::parse`](super::index::parse) into the same [`Index`], and
//! [`window`] reads a BAM through it as through a BAI. The command line
//! looks for the `.csi` first, as samtools does.
//!
//! # What comes out
//!
//! The records over the window, which [`sam`] writes as the SAM text the
//! pileup and split-read readers already take, and [`depth`] counts the way
//! `samtools depth -a` does by default: every base of the window, reads that
//! are unmapped, secondary, failing quality checks or duplicates left out,
//! and a deletion or a skipped intron not counted as covered.
//!
//! # What is refused
//!
//! A file that is not BAM, a block or a record cut short, a record whose
//! lengths do not add up to its size, and an index that is not a BAI. A
//! window on a sequence the BAM does not have is refused naming the ones it
//! has. None of it panics.

use std::io::{Read, Seek};

use crate::Region;

use super::ReadError;

/// The reader a BAM is read through, which lived here until the bgzipped text
/// files a tabix index points into came to need it too.
pub use super::bgzf::Bgzf;

/// The flag bits `samtools depth` leaves out by default: unmapped, secondary,
/// failing quality checks, and duplicate.
const NOT_COUNTED: u16 = 0x4 | 0x100 | 0x200 | 0x400;

/// A record larger than this is a damaged file rather than a read.
const LARGEST_RECORD: usize = 1 << 28;

/// Fills `buf` from the BAM, and says whether there was anything at all:
/// false when the file ends before the first byte, which between records is
/// where a BAM ends.
fn read_into<R: Read + Seek>(bgzf: &mut Bgzf<R>, buf: &mut [u8]) -> Result<bool, ReadError> {
    match bgzf.fill(buf)? {
        got if got == buf.len() => Ok(true),
        0 => Ok(false),
        _ => Err(ReadError::whole("the BAM ends in the middle of a record")),
    }
}

fn exact<R: Read + Seek>(bgzf: &mut Bgzf<R>, buf: &mut [u8]) -> Result<(), ReadError> {
    if read_into(bgzf, buf)? || buf.is_empty() {
        Ok(())
    } else {
        Err(ReadError::whole("the BAM ends early"))
    }
}

/// What a BAM says about itself before its reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// The SAM header, as text: `@HD`, `@SQ`, `@RG` and `@PG` lines.
    pub text: String,
    /// Each reference sequence's name and length, in the order the reads
    /// number them.
    pub references: Vec<(String, u64)>,
}

impl Header {
    /// Whether the reads are sorted by position, which lets a reader stop at
    /// the end of the window.
    fn sorted(&self) -> bool {
        self.text
            .lines()
            .find(|line| line.starts_with("@HD"))
            .is_some_and(|line| line.split('\t').any(|field| field == "SO:coordinate"))
    }
}

fn header<R: Read + Seek>(bgzf: &mut Bgzf<R>) -> Result<Header, ReadError> {
    let mut magic = [0u8; 4];
    exact(bgzf, &mut magic)?;
    if &magic != b"BAM\x01" {
        return Err(ReadError::whole(
            "not BAM: the file does not start with BAM's magic",
        ));
    }
    let text_length = count(bgzf, 1 << 30, "header")?;
    let mut text = vec![0u8; text_length];
    exact(bgzf, &mut text)?;
    // The text may be padded with noughts.
    let text = String::from_utf8_lossy(&text)
        .trim_end_matches('\0')
        .to_string();
    let references = count(bgzf, 1 << 24, "reference list")?;
    let mut named = Vec::with_capacity(references.min(4096));
    for _ in 0..references {
        let name_length = count(bgzf, 1 << 16, "reference name")?;
        let mut name = vec![0u8; name_length];
        exact(bgzf, &mut name)?;
        let name = String::from_utf8_lossy(&name)
            .trim_end_matches('\0')
            .to_string();
        let mut length = [0u8; 4];
        exact(bgzf, &mut length)?;
        named.push((name, u64::from(u32::from_le_bytes(length))));
    }
    Ok(Header {
        text,
        references: named,
    })
}

/// A count the format stores as an `int32`, held to what a sound file holds.
fn count<R: Read + Seek>(bgzf: &mut Bgzf<R>, most: usize, what: &str) -> Result<usize, ReadError> {
    let mut bytes = [0u8; 4];
    exact(bgzf, &mut bytes)?;
    let value = i32::from_le_bytes(bytes);
    usize::try_from(value)
        .ok()
        .filter(|value| *value <= most)
        .ok_or_else(|| {
            ReadError::whole(format!(
                "the BAM's {what} has an impossible length, {value}"
            ))
        })
}

/// One aligned read, as BAM stores it.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// Which reference, by its place in [`Header::references`], or -1.
    pub reference: i32,
    /// The leftmost aligned position, 0-based, or -1.
    pub position: i64,
    /// The mapping quality.
    pub mapping_quality: u8,
    /// The SAM flag.
    pub flag: u16,
    /// The read's name.
    pub name: String,
    /// The CIGAR, each operation as BAM packs it: the length shifted left four
    /// bits, and the operation's number in `MIDNSHP=X` order.
    pub cigar: Vec<u32>,
    /// The bases, as letters.
    pub sequence: Vec<u8>,
    /// The base qualities, as Phred scores, or empty where the file has none.
    pub qualities: Vec<u8>,
    /// The mate's reference, as `reference` is.
    pub mate_reference: i32,
    /// The mate's position, as `position` is.
    pub mate_position: i64,
    /// The observed template length.
    pub template_length: i32,
    /// The optional fields, in BAM's binary form.
    pub tags: Vec<u8>,
}

impl Record {
    /// One past the last reference base the read covers.
    pub fn end(&self) -> i64 {
        let covered: i64 = self
            .cigar
            .iter()
            .filter(|op| matches!(*op & 0xf, 0 | 2 | 3 | 7 | 8))
            .map(|op| i64::from(*op >> 4))
            .sum();
        self.position + covered.max(1)
    }
}

/// The next record, or `None` at the end of the file.
fn record<R: Read + Seek>(bgzf: &mut Bgzf<R>) -> Result<Option<Record>, ReadError> {
    let mut size = [0u8; 4];
    if !read_into(bgzf, &mut size)? {
        return Ok(None);
    }
    let size = i32::from_le_bytes(size);
    let size = usize::try_from(size)
        .ok()
        .filter(|size| (32..=LARGEST_RECORD).contains(size))
        .ok_or_else(|| {
            ReadError::whole(format!("a BAM record claims an impossible size, {size}"))
        })?;
    let mut body = vec![0u8; size];
    exact(bgzf, &mut body)?;
    let int = |at: usize| i32::from_le_bytes([body[at], body[at + 1], body[at + 2], body[at + 3]]);
    let short = |at: usize| u16::from_le_bytes([body[at], body[at + 1]]);
    let name_length = usize::from(body[8]);
    let operations = usize::from(short(12));
    let bases = usize::try_from(int(16))
        .map_err(|_| ReadError::whole("a BAM record has a negative sequence length"))?;
    let cigar_at = 32 + name_length;
    let sequence_at = cigar_at + operations * 4;
    let qualities_at = sequence_at + bases.div_ceil(2);
    let tags_at = qualities_at + bases;
    if tags_at > body.len() {
        return Err(ReadError::whole(
            "a BAM record's lengths add up to more than its size",
        ));
    }

    let name = String::from_utf8_lossy(&body[32..cigar_at])
        .trim_end_matches('\0')
        .to_string();
    let mut cigar: Vec<u32> = (0..operations)
        .map(|op| {
            let at = cigar_at + op * 4;
            u32::from_le_bytes([body[at], body[at + 1], body[at + 2], body[at + 3]])
        })
        .collect();
    const LETTERS: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let sequence: Vec<u8> = (0..bases)
        .map(|base| {
            let byte = body[sequence_at + base / 2];
            let code = if base % 2 == 0 { byte >> 4 } else { byte & 0xf };
            LETTERS[usize::from(code)]
        })
        .collect();
    let qualities = &body[qualities_at..tags_at];
    let qualities = if qualities.first() == Some(&0xff) {
        Vec::new()
    } else {
        qualities.to_vec()
    };
    let mut tags = body[tags_at..].to_vec();

    // A CIGAR of more than 65,535 operations does not fit the record, which
    // long reads reach. BAM keeps it in the CG tag and writes a placeholder
    // of the read's length soft clipped over a skip in its place. The tag has
    // said all it had to once the CIGAR is back, and samtools drops it then.
    if cigar.len() == 2
        && cigar[0] & 0xf == 4
        && (cigar[0] >> 4) as usize == bases
        && cigar[1] & 0xf == 3
    {
        if let Some(long) = long_cigar(&tags) {
            cigar = long;
            tags = without(&tags, *b"CG");
        }
    }

    Ok(Some(Record {
        reference: int(0),
        position: i64::from(int(4)),
        mapping_quality: body[9],
        flag: short(14),
        name,
        cigar,
        sequence,
        qualities,
        mate_reference: int(20),
        mate_position: i64::from(int(24)),
        template_length: int(28),
        tags,
    }))
}

/// The CIGAR a long read keeps in its CG tag.
fn long_cigar(tags: &[u8]) -> Option<Vec<u32>> {
    let mut found = None;
    walk_tags(tags, |name, kind, value| {
        if name == *b"CG" && kind == b'B' && value.first() == Some(&b'I') && value.len() >= 5 {
            let n = u32::from_le_bytes([value[1], value[2], value[3], value[4]]) as usize;
            let ops: Vec<u32> = value[5..]
                .chunks_exact(4)
                .take(n)
                .map(|op| u32::from_le_bytes([op[0], op[1], op[2], op[3]]))
                .collect();
            if ops.len() == n {
                found = Some(ops);
            }
        }
    });
    found
}

/// The optional fields with the one called `name` left out.
fn without(tags: &[u8], name: [u8; 2]) -> Vec<u8> {
    let mut kept = Vec::with_capacity(tags.len());
    walk_tags(tags, |field, kind, value| {
        if field != name {
            kept.extend_from_slice(&field);
            kept.push(kind);
            kept.extend_from_slice(value);
        }
    });
    kept
}

/// Each optional field: its name, its type, and its value's bytes.
///
/// Stops at the first field that does not fit, which is where a damaged
/// record's tags end.
fn walk_tags(tags: &[u8], mut each: impl FnMut([u8; 2], u8, &[u8])) {
    let mut at = 0;
    while at + 3 <= tags.len() {
        let name = [tags[at], tags[at + 1]];
        let kind = tags[at + 2];
        let start = at + 3;
        let length = match kind {
            b'A' | b'c' | b'C' => Some(1),
            b's' | b'S' => Some(2),
            b'i' | b'I' | b'f' => Some(4),
            b'Z' | b'H' => tags[start..]
                .iter()
                .position(|byte| *byte == 0)
                .map(|end| end + 1),
            b'B' => tags.get(start..start + 5).and_then(|head| {
                let width = match head[0] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                let n = u32::from_le_bytes([head[1], head[2], head[3], head[4]]) as usize;
                n.checked_mul(width).map(|bytes| 5 + bytes)
            }),
            _ => None,
        };
        let Some(length) = length.filter(|length| start + length <= tags.len()) else {
            return;
        };
        each(name, kind, &tags[start..start + length]);
        at = start + length;
    }
}

/// A BAM's index, which [`index`] reads. The same type holds a CSI, and
/// [`window`] reads a BAM through either.
pub use super::index::Index;

/// Reads a `.bai` index.
///
/// # Errors
///
/// Bytes that are not a BAI, or one cut short.
pub fn index(data: &[u8]) -> Result<Index, ReadError> {
    match data.get(..4) {
        None => Err(ReadError::whole("the BAM index is cut short")),
        Some(magic) if magic != b"BAI\x01" => Err(ReadError::whole(
            "not a BAI index: it does not start with BAI's magic",
        )),
        Some(_) => super::index::parse(data),
    }
}

/// The header of a BAM, and every record that overlaps `region`.
///
/// Read through `index` where one is given, and otherwise from the start of
/// the file, leaving it at the end of the window when the file says it is
/// sorted by position.
///
/// # Errors
///
/// A file that is not BAM or is damaged, and a region on a sequence the BAM
/// does not have, naming the ones it has.
pub fn window<R: Read + Seek>(
    reader: R,
    index: Option<&Index>,
    region: &Region,
) -> Result<(Header, Vec<Record>), ReadError> {
    let mut bgzf = Bgzf::new(reader);
    let header = header(&mut bgzf)?;
    let Some(target) = header
        .references
        .iter()
        .position(|(name, _)| name == region.seq())
    else {
        let mut named: Vec<&str> = header
            .references
            .iter()
            .take(12)
            .map(|(name, _)| name.as_str())
            .collect();
        if header.references.len() > named.len() {
            named.push("and more");
        }
        return Err(ReadError::whole(format!(
            "the BAM has no sequence called {}; it has {}",
            region.seq(),
            named.join(", ")
        )));
    };
    let sorted = header.sorted() || index.is_some();
    let (start, end) = (region.start() as i64, region.end() as i64);
    let mut records = Vec::new();

    let mut keep = |record: Record| -> bool {
        // False once the reads have passed the window, which only a sorted
        // file can say.
        if record.reference != target as i32 {
            return !(sorted && (record.reference > target as i32 || record.reference < 0));
        }
        if record.position >= end {
            return !sorted;
        }
        if record.end() > start {
            records.push(record);
        }
        true
    };

    match index {
        Some(index) => {
            'chunks: for (begin, stop) in index.chunks(target, region.start(), region.end()) {
                bgzf.seek(begin)?;
                while bgzf.tell() < stop {
                    let Some(record) = record(&mut bgzf)? else {
                        break 'chunks;
                    };
                    if !keep(record) {
                        break 'chunks;
                    }
                }
            }
        }
        None => {
            while let Some(record) = record(&mut bgzf)? {
                if !keep(record) {
                    break;
                }
            }
        }
    }
    Ok((header, records))
}

/// The records of one read, by its name, and of the reads Dorado split out
/// of it, which name it as their parent in `pi:Z`, read from the start of the
/// file to its end: a basecaller's BAM is neither sorted nor indexed.
///
/// # Errors
///
/// A file that is not BAM, or is damaged.
pub fn named<R: Read + Seek>(reader: R, name: &str) -> Result<(Header, Vec<Record>), ReadError> {
    let mut bgzf = Bgzf::new(reader);
    let header = header(&mut bgzf)?;
    let mut records = Vec::new();
    while let Some(record) = record(&mut bgzf)? {
        let mut parent = None;
        walk_tags(&record.tags, |tag, kind, value| {
            if tag == *b"pi" && kind == b'Z' {
                parent = value.split_last().map(|(_, text)| text.to_vec());
            }
        });
        if record.name == name || parent.as_deref() == Some(name.as_bytes()) {
            records.push(record);
        }
    }
    Ok((header, records))
}

/// Just the header of a BAM, for the lengths of its sequences.
///
/// # Errors
///
/// A file that is not BAM, or is damaged.
pub fn header_of<R: Read + Seek>(reader: R) -> Result<Header, ReadError> {
    header(&mut Bgzf::new(reader))
}

/// The records as SAM text, the header's lines first.
pub fn sam(header: &Header, records: &[Record]) -> String {
    let mut out = String::new();
    if header.text.trim().is_empty() {
        for (name, length) in &header.references {
            out.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
        }
    } else {
        out.push_str(header.text.trim_end_matches('\n'));
        out.push('\n');
    }
    let reference = |id: i32| -> &str {
        usize::try_from(id)
            .ok()
            .and_then(|id| header.references.get(id))
            .map_or("*", |(name, _)| name.as_str())
    };
    for record in records {
        const OPS: &[u8; 9] = b"MIDNSHP=X";
        let cigar: String = if record.cigar.is_empty() {
            "*".to_string()
        } else {
            record
                .cigar
                .iter()
                .map(|op| {
                    let letter = OPS.get((op & 0xf) as usize).copied().unwrap_or(b'?');
                    format!("{}{}", op >> 4, char::from(letter))
                })
                .collect()
        };
        let mate = if record.mate_reference < 0 {
            "*"
        } else if record.mate_reference == record.reference {
            "="
        } else {
            reference(record.mate_reference)
        };
        let sequence = if record.sequence.is_empty() {
            "*".to_string()
        } else {
            String::from_utf8_lossy(&record.sequence).into_owned()
        };
        let qualities = if record.qualities.is_empty() {
            "*".to_string()
        } else {
            record
                .qualities
                .iter()
                .map(|q| char::from(q.saturating_add(33).min(126)))
                .collect()
        };
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            if record.name.is_empty() {
                "*"
            } else {
                &record.name
            },
            record.flag,
            reference(record.reference),
            record.position + 1,
            record.mapping_quality,
            cigar,
            mate,
            record.mate_position + 1,
            record.template_length,
            sequence,
            qualities
        ));
        walk_tags(&record.tags, |name, kind, value| {
            out.push('\t');
            out.push_str(&tag_text(name, kind, value));
        });
        out.push('\n');
    }
    out
}

/// One optional field as SAM writes it.
fn tag_text(name: [u8; 2], kind: u8, value: &[u8]) -> String {
    let name = String::from_utf8_lossy(&name).into_owned();
    let number = |kind: u8, bytes: &[u8]| -> String {
        match kind {
            b'c' => (bytes[0] as i8).to_string(),
            b'C' => bytes[0].to_string(),
            b's' => i16::from_le_bytes([bytes[0], bytes[1]]).to_string(),
            b'S' => u16::from_le_bytes([bytes[0], bytes[1]]).to_string(),
            b'i' => i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]).to_string(),
            b'I' => u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]).to_string(),
            _ => f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]).to_string(),
        }
    };
    match kind {
        b'A' => format!("{name}:A:{}", char::from(value[0])),
        b'c' | b'C' | b's' | b'S' | b'i' | b'I' => format!("{name}:i:{}", number(kind, value)),
        b'f' => format!("{name}:f:{}", number(kind, value)),
        b'Z' | b'H' => format!(
            "{name}:{}:{}",
            char::from(kind),
            String::from_utf8_lossy(&value[..value.len() - 1])
        ),
        _ => {
            // B: a subtype, a count, and the values.
            let width = match value[0] {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                _ => 4,
            };
            let values: Vec<String> = value[5..]
                .chunks_exact(width)
                .map(|bytes| number(value[0], bytes))
                .collect();
            format!("{name}:B:{},{}", char::from(value[0]), values.join(","))
        }
    }
}

/// The depth over `region`, one count per base, as `samtools depth -a` counts
/// by default.
pub fn depth(records: &[Record], region: &Region) -> Vec<u32> {
    let (start, end) = (region.start() as i64, region.end() as i64);
    let mut depth = vec![0u32; (end - start).max(0) as usize];
    for record in records {
        if record.flag & NOT_COUNTED != 0 || record.position < 0 {
            continue;
        }
        let mut at = record.position;
        for op in &record.cigar {
            let length = i64::from(op >> 4);
            match op & 0xf {
                // Aligned bases: M, = and X.
                0 | 7 | 8 => {
                    for base in at.max(start)..(at + length).min(end) {
                        depth[(base - start) as usize] += 1;
                    }
                    at += length;
                }
                // A deletion and a skip cover the reference without a base.
                2 | 3 => at += length,
                _ => {}
            }
        }
    }
    depth
}

/// Depth as bedGraph, one line per run of the same count.
pub fn bedgraph(sequence: &str, region: &Region, depth: &[u32]) -> String {
    let mut out = String::new();
    let mut run = 0;
    while run < depth.len() {
        let mut stop = run + 1;
        while stop < depth.len() && depth[stop] == depth[run] {
            stop += 1;
        }
        let start = region.start() + run as u64;
        out.push_str(&format!(
            "{sequence}\t{start}\t{}\t{}\n",
            start + (stop - run) as u64,
            depth[run]
        ));
        run = stop;
    }
    out
}

/// A dozen reads written by samtools, and the index it wrote for them, for
/// the tests here and for the command line's.
#[cfg(test)]
pub(crate) mod fixture {
    pub(crate) const BAM: [u8; 467] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x9d, 0x00, 0x73, 0x72, 0xf4, 0x65, 0xec, 0x65, 0x60, 0x60, 0x70, 0xf0, 0x70, 0xe1,
        0x0c, 0xf3, 0xb3, 0x32, 0xd4, 0x33, 0xe3, 0x0c, 0xf6, 0xb7, 0x4a, 0xce, 0xcf, 0x2f, 0x4a,
        0xc9, 0xcc, 0x4b, 0x2c, 0x49, 0xe5, 0x72, 0x08, 0x0e, 0xe4, 0x0c, 0xf6, 0xb3, 0x4a, 0xce,
        0x28, 0x32, 0xe4, 0xf4, 0x01, 0x2a, 0x30, 0x30, 0x30, 0x40, 0x12, 0x33, 0x02, 0x89, 0x99,
        0x82, 0x84, 0x02, 0xdc, 0x39, 0x3d, 0x5d, 0xac, 0x8a, 0x13, 0x73, 0x4b, 0xf2, 0xf3, 0x73,
        0x8a, 0x39, 0x03, 0xfc, 0x10, 0x6c, 0xb0, 0xb9, 0x46, 0x26, 0x9c, 0xce, 0x3e, 0x70, 0x31,
        0x85, 0xb2, 0xcc, 0xd4, 0x72, 0x05, 0xdd, 0x24, 0x05, 0xdd, 0x7c, 0x85, 0x92, 0xcc, 0xbc,
        0x4a, 0xbd, 0xa4, 0xc4, 0x5c, 0x08, 0x03, 0xa8, 0x82, 0x8b, 0x09, 0xe8, 0x20, 0x56, 0x20,
        0x06, 0x59, 0xca, 0xf0, 0x82, 0x19, 0xce, 0x31, 0x62, 0xf8, 0xc2, 0xc8, 0xc0, 0x00, 0x00,
        0x86, 0x85, 0x2f, 0xf9, 0xb3, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00,
        0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x18, 0x01, 0x65, 0x91, 0x4f, 0x4b,
        0xc3, 0x30, 0x18, 0xc6, 0x7f, 0xab, 0x65, 0xa8, 0xa0, 0x6c, 0x11, 0x2f, 0x82, 0x1e, 0xa4,
        0xe8, 0x60, 0x01, 0xdb, 0x69, 0x05, 0xd9, 0x0e, 0x96, 0x79, 0xc8, 0x0e, 0x15, 0xa1, 0x08,
        0x22, 0x08, 0x6a, 0xfd, 0xb3, 0xa3, 0x08, 0xde, 0xd7, 0x6f, 0xe2, 0x87, 0xf4, 0x3e, 0xc9,
        0xf2, 0x76, 0x2b, 0xd9, 0x73, 0x09, 0x24, 0xf9, 0xbd, 0x79, 0x9e, 0x27, 0x57, 0x38, 0x6d,
        0x01, 0xc1, 0x68, 0xa2, 0x5a, 0xc0, 0x36, 0x30, 0x17, 0xd9, 0xb3, 0x17, 0x7e, 0x01, 0x65,
        0x94, 0x51, 0xbd, 0xa5, 0x6e, 0xf3, 0x31, 0x4f, 0x02, 0xef, 0x5a, 0xf8, 0x68, 0xa2, 0x42,
        0x3a, 0xec, 0x78, 0xf0, 0x2b, 0x17, 0xc0, 0x1d, 0x70, 0x2c, 0xeb, 0x4d, 0x16, 0x55, 0xca,
        0x74, 0x82, 0x5a, 0x76, 0x5a, 0x56, 0x94, 0x7f, 0x0f, 0xd9, 0x63, 0x39, 0xfd, 0x1e, 0xe8,
        0x7e, 0x12, 0xc7, 0x3a, 0x39, 0xcf, 0x75, 0xc2, 0x48, 0x1e, 0xe8, 0x8a, 0xbb, 0x0d, 0x71,
        0xda, 0x7c, 0xa0, 0xe4, 0x5a, 0x6e, 0xd8, 0x55, 0x99, 0x2a, 0x9a, 0xcd, 0x6b, 0x0d, 0x84,
        0xdf, 0xb3, 0x3c, 0x8b, 0x74, 0xad, 0x4d, 0x8f, 0x7f, 0x63, 0x26, 0xe9, 0xea, 0x68, 0x35,
        0xb5, 0xbf, 0xea, 0x24, 0xf4, 0xa9, 0xf7, 0x35, 0xea, 0x5e, 0xa8, 0xc3, 0x86, 0x57, 0xbf,
        0xc9, 0x8f, 0x45, 0x01, 0xe3, 0xb6, 0x2b, 0x42, 0x99, 0x6e, 0x54, 0xad, 0x1a, 0x2d, 0x24,
        0x7f, 0x1a, 0xeb, 0xbe, 0x4e, 0xf3, 0xb4, 0xd0, 0x97, 0xb1, 0x8e, 0x87, 0x9c, 0xc9, 0xe0,
        0xd3, 0xa5, 0x9d, 0x92, 0xb6, 0xec, 0x0d, 0x81, 0x03, 0xe0, 0x93, 0x67, 0x67, 0xc7, 0x8d,
        0x3a, 0x01, 0xec, 0x4f, 0x86, 0x8d, 0x5f, 0xc5, 0xf3, 0x32, 0xe5, 0x0b, 0xe8, 0x35, 0x76,
        0x03, 0x8c, 0xb2, 0x4c, 0xe8, 0xdd, 0xfc, 0xc1, 0xa5, 0xfc, 0x07, 0x8f, 0xbb, 0xc6, 0xde,
        0x2c, 0x02, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06,
        0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00,
    ];

    pub(crate) const BAI: [u8; 176] = [
        0x42, 0x41, 0x49, 0x01, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x49, 0x12, 0x00,
        0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0xd6, 0x01,
        0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4a, 0x92, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0xd6, 0x01, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00,
        0x00, 0x00, 0x49, 0x12, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0xd6, 0x01, 0x9e, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x02, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4a, 0x92, 0x00, 0x00,
        0x02, 0x00, 0x00, 0x00, 0xd6, 0x01, 0x9e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x9e,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0xd6, 0x01, 0x9e, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// The same reads indexed by `samtools index -c` 1.24, which chose a tree
    /// of one level for sequences of 1,000 and 500 bases.
    pub(crate) const CSI: [u8; 104] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x4b, 0x00, 0x73, 0x0e, 0xf6, 0x64, 0xe4, 0x63, 0x40, 0x00, 0x26, 0x28, 0x86, 0x80,
        0x79, 0x60, 0x92, 0x11, 0x89, 0x7d, 0x8d, 0x71, 0x1e, 0x5c, 0x1d, 0xb2, 0x1e, 0x74, 0x79,
        0x76, 0x06, 0x54, 0x00, 0x53, 0x0f, 0x93, 0x67, 0x44, 0x62, 0x33, 0x30, 0x61, 0x37, 0x13,
        0x5d, 0x9e, 0x11, 0xcd, 0x4c, 0x18, 0x1f, 0x00, 0xa4, 0x1c, 0x51, 0x22, 0xc4, 0x00, 0x00,
        0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
        0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
}

#[cfg(test)]
mod tests {
    use super::fixture::{BAI, BAM};
    use super::*;
    use std::io::Cursor;

    // A dozen reads written by samtools, and its index: every CIGAR operation,
    // a secondary, a duplicate and an unmapped read, a read with no qualities
    // and one with no bases, tags of several types, a pair, and a second
    // sequence. The SAM they came from:
    //
    //   a  0     chr1 10 10M              NM:i:0
    //   b  16    chr1 15 3S5M2D5M         AS:i:-12  XA:Z:chr2,+100,13M,1
    //   c  0     chr1 18 4M1I4M           no qualities
    //   d  256   chr1 20 8M               secondary
    //   e  1024  chr1 22 8M               duplicate
    //   f  0     chr1 30 5M100N5M         SA:Z:chr2,50,+,5M5S,60,0;
    //   g  99    chr1 40 6M               paired, mate at 60
    //   h  0     chr2 5  7M               no bases
    //   u  4     *                        unmapped
    /// What `samtools view tiny.bam chr1:15-20` prints.
    const VIEW_15_20: &str = "\
a\t0\tchr1\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n\
b\t16\tchr1\t15\t30\t3S5M2D5M\t*\t0\t0\tGGGACGTAACGTA\t#########IIII\tAS:i:-12\tXA:Z:chr2,+100,13M,1\n\
c\t0\tchr1\t18\t60\t4M1I4M\t*\t0\t0\tACGTTACGT\t*\n\
d\t256\tchr1\t20\t0\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
";

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    fn names(records: &[Record]) -> Vec<&str> {
        records.iter().map(|record| record.name.as_str()).collect()
    }

    #[test]
    fn the_header_names_each_sequence_and_its_length() {
        let header = header_of(Cursor::new(&BAM[..])).unwrap();
        assert_eq!(
            header.references,
            [("chr1".to_string(), 1000), ("chr2".to_string(), 500)]
        );
        assert!(header.sorted());
    }

    #[test]
    fn a_window_holds_what_samtools_view_prints_for_it() {
        let index = index(&BAI).unwrap();
        let (header, through) =
            window(Cursor::new(&BAM[..]), Some(&index), &region("chr1:15-20")).unwrap();
        let (_, scanned) = window(Cursor::new(&BAM[..]), None, &region("chr1:15-20")).unwrap();
        assert_eq!(names(&through), ["a", "b", "c", "d"]);
        assert_eq!(through, scanned, "the index and a scan disagree");
        let text = sam(&header, &through);
        let records: String = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(records, VIEW_15_20);
    }

    /// A CSI is read through the same index as a BAI, with the tree its
    /// writer chose in place of BAI's, and finds the same reads: over every
    /// window of every sequence, and the same as a scan of the file.
    #[test]
    fn a_bam_with_a_csi_beside_it_reads_what_its_bai_reads() {
        let bai = index(&BAI).unwrap();
        let csi = super::super::index::parse(&super::fixture::CSI).unwrap();
        assert_eq!(csi.depth(), 0, "samtools chose one level for 1,000 bases");
        for (sequence, length) in [("chr1", 1000u64), ("chr2", 500)] {
            for start in (0..length).step_by(7) {
                for width in [1, 5, 40, 1000] {
                    let end = (start + width).min(length);
                    let over = Region::new(sequence, start, end).unwrap();
                    let read = |index: Option<&Index>| {
                        window(Cursor::new(&BAM[..]), index, &over).unwrap().1
                    };
                    let through = read(Some(&csi));
                    assert_eq!(through, read(Some(&bai)), "{over}");
                    assert_eq!(through, read(None), "{over}");
                }
            }
        }
        assert_eq!(
            names(
                &window(Cursor::new(&BAM[..]), Some(&csi), &region("chr1:15-20"))
                    .unwrap()
                    .1
            ),
            ["a", "b", "c", "d"]
        );
        // The index of a BAM is a BAI where a caller asks for one.
        let error = index(&super::fixture::CSI).unwrap_err();
        assert!(error.to_string().contains("not a BAI"), "{error}");
    }

    /// A read is found by its name from one end of the file to the other,
    /// mapped or not, since a basecaller's BAM is neither sorted nor indexed.
    #[test]
    fn a_read_is_found_by_its_name_mapped_or_not() {
        let (header, found) = named(Cursor::new(&BAM[..]), "u").unwrap();
        assert_eq!(names(&found), ["u"]);
        assert!(sam(&header, &found).contains("u\t4\t*"));
        let (_, found) = named(Cursor::new(&BAM[..]), "b").unwrap();
        assert_eq!(names(&found), ["b"]);
        let (_, none) = named(Cursor::new(&BAM[..]), "zz").unwrap();
        assert!(none.is_empty());
    }

    #[test]
    fn depth_counts_as_samtools_depth_does() {
        // `samtools depth -a -r chr1:10-35`: the secondary and the duplicate
        // are not counted, nor the deletion in b, nor the skip in f.
        let (_, records) = window(Cursor::new(&BAM[..]), None, &region("chr1:10-35")).unwrap();
        let depth = depth(&records, &region("chr1:10-35"));
        assert_eq!(
            depth,
            [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 1, 1, 2, 2, 2, 2, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0]
        );
        assert_eq!(
            bedgraph("chr1", &region("chr1:10-14"), &depth[..5]),
            "chr1\t9\t14\t1\n"
        );
    }

    #[test]
    fn a_sequence_the_bam_has_not_got_is_refused_with_the_ones_it_has() {
        let error = window(Cursor::new(&BAM[..]), None, &region("chr3:1-10")).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("no sequence called chr3; it has chr1, chr2"),
            "{error}"
        );
    }

    #[test]
    fn a_damaged_file_or_index_is_an_error_and_never_a_panic() {
        for cut in 0..BAM.len() {
            let _ = window(Cursor::new(&BAM[..cut]), None, &region("chr1:1-1000"));
        }
        for cut in 0..BAI.len() {
            let _ = index(&BAI[..cut]);
        }
        let error = window(
            Cursor::new(&b"not a bam at all"[..]),
            None,
            &region("chr1:1-10"),
        )
        .unwrap_err();
        assert!(error.to_string().contains("BGZF"), "{error}");
    }

    #[test]
    fn a_long_cigar_comes_back_from_its_tag_and_the_tag_goes() {
        // CG:B:I with three operations, 5M 2I 5M, and an NM beside it.
        let mut tags = b"CGBI".to_vec();
        tags.extend_from_slice(&3u32.to_le_bytes());
        for op in [5u32 << 4, (2u32 << 4) | 1, 5u32 << 4] {
            tags.extend_from_slice(&op.to_le_bytes());
        }
        tags.extend_from_slice(b"NMC\x02");
        assert_eq!(long_cigar(&tags), Some(vec![5 << 4, (2 << 4) | 1, 5 << 4]));
        assert_eq!(without(&tags, *b"CG"), b"NMC\x02");
    }
}
