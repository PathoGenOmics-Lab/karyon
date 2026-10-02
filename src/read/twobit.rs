//! 2bit, a genome's sequences four bases to a byte, read a window at a time.
//!
//! UCSC's 2bit keeps each base in two bits, T, C, A and G as 0 to 3, the
//! first base of a byte in its top two, with the runs of N and the runs in
//! lower case listed apart, since two bits have no room for either. An index
//! at the front says where each sequence starts, so [`bases`] reads the index,
//! the lists of one sequence and the quarter of a byte per base of the window,
//! and none of the rest: a window of ten thousand bases is 2,500 bytes of
//! bases read, where the same window of a FASTA reads the whole FASTA first.
//!
//! The bases come back as `twoBitToFa` writes them, which a test holds them
//! to: N over a run of N, and lower case over a run the file marks, so a
//! soft-masked reference says what it says in a FASTA, and an N inside a
//! masked run is an `n`.
//!
//! # Coordinates
//!
//! None in the file: a sequence starts at its own first base, so byte n of it
//! is 0-based position n, as in a FASTA.
//!
//! # What is read
//!
//! Both versions of the format: 0, whose offsets are 32 bits and which every
//! 2bit under 4 GB is, and 1, whose offsets are 64, as `faToTwoBit -long`
//! writes; and either byte order, which the signature says.
//!
//! # What is refused
//!
//! A file that is not a 2bit, one damaged or cut short, and a sequence the
//! file does not have, or has twice, named with the ones it has. A file of
//! one sequence answers for it whatever the window calls it, as a FASTA of
//! one record does.
//!
//! ```
//! use std::io::Cursor;
//! use karyon::{read, Region};
//!
//! # let bytes = include_bytes!("fixtures/ref.2bit").to_vec();
//! let region = Region::parse("chr1:3-12")?;
//! let window = read::twobit::bases(Cursor::new(bytes), &region)?;
//! assert_eq!((window.start, window.bases.as_slice()), (2, &b"GTNNNNNTGA"[..]));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::io::{Read, Seek};

use super::bytes::{Bytes, File};
use super::ReadError;
use crate::Region;

/// The signature a 2bit starts with, read in its own byte order.
const MAGIC: u32 = 0x1A41_2743;

/// The bases of one sequence of a 2bit over a window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bases {
    /// The sequence's name, as the file gives it.
    pub name: String,
    /// How many bases the whole sequence holds.
    pub length: u64,
    /// The 0-based position of the first of `bases` on the sequence.
    pub start: u64,
    /// The bases the window holds, as far as the sequence runs, and none for
    /// a window that starts past its end.
    pub bases: Vec<u8>,
}

/// The sequences a 2bit holds, each with its length, in the order of its
/// index.
///
/// # Errors
///
/// A file that is not a 2bit, or is damaged.
pub fn sequences<R: Read + Seek>(reader: R) -> Result<Vec<(String, u64)>, ReadError> {
    let mut file = File::new(reader, "2bit")?;
    let index = index(&mut file)?;
    let mut out = Vec::with_capacity(index.entries.len());
    for (name, offset) in index.entries {
        let size = file.read_at(offset, 4)?;
        let length = Bytes::new(&size, index.big, "2bit").u32()?;
        out.push((name, u64::from(length)));
    }
    Ok(out)
}

/// The bases of the sequence `region` is on, over `region`.
///
/// # Errors
///
/// A file that is not a 2bit or is damaged, and a region on a sequence the
/// file does not have, or has more than once.
pub fn bases<R: Read + Seek>(reader: R, region: &Region) -> Result<Bases, ReadError> {
    let mut file = File::new(reader, "2bit")?;
    let index = index(&mut file)?;
    let (name, offset) = match index.entries.as_slice() {
        [one] => one.clone(),
        many => {
            let named: Vec<&(String, u64)> = many
                .iter()
                .filter(|(name, _)| name == region.seq())
                .collect();
            match named.as_slice() {
                [one] => (*one).clone(),
                [] => {
                    let mut held: Vec<&str> = many
                        .iter()
                        .take(12)
                        .map(|(name, _)| name.as_str())
                        .collect();
                    if many.len() > held.len() {
                        held.push("and more");
                    }
                    return Err(ReadError::whole(format!(
                        "the 2bit has no sequence called {}; it has {}",
                        region.seq(),
                        if held.is_empty() {
                            "none".to_string()
                        } else {
                            held.join(", ")
                        }
                    )));
                }
                twice => {
                    return Err(ReadError::whole(format!(
                        "the 2bit has {} sequences called {}, so the name does not pick one",
                        twice.len(),
                        region.seq()
                    )))
                }
            }
        }
    };
    let big = index.big;
    let numbers = |data: &[u8]| -> Result<Vec<u32>, ReadError> {
        let mut bytes = Bytes::new(data, big, "2bit");
        (0..data.len() / 4).map(|_| bytes.u32()).collect()
    };
    // The record: its length, its runs of N, its runs in lower case and a
    // word kept for later, then the bases. Each list is read only where the
    // file holds all of it, so a damaged count asks for no more memory than
    // the file is long.
    let after = |at: u64, by: u64| {
        at.checked_add(by)
            .ok_or_else(|| super::bytes::short("2bit"))
    };
    let head = file.read_at(offset, 8)?;
    let mut bytes = Bytes::new(&head, big, "2bit");
    let length = u64::from(bytes.u32()?);
    let unknown = u64::from(bytes.u32()?);
    let at = after(offset, 8)?;
    let unknown = numbers(&file.read_at(at, unknown * 8)?)?;
    let at = after(at, (unknown.len() * 4) as u64)?;
    let masked = file.read_at(at, 4)?;
    let masked = u64::from(Bytes::new(&masked, big, "2bit").u32()?);
    let at = after(at, 4)?;
    let masked = numbers(&file.read_at(at, masked * 8)?)?;
    let packed_at = after(after(at, (masked.len() * 4) as u64)?, 4)?;

    let from = region.start().min(length);
    let to = region.end().min(length);
    let mut bases = Vec::new();
    if from < to {
        // A quarter of a byte a base, from the byte the first is in to the
        // one the last is in: read before anything is allocated for them.
        let first = from / 4;
        let packed = file.read_at(after(packed_at, first)?, (to - 1) / 4 + 1 - first)?;
        bases.reserve_exact((to - from) as usize);
        for position in from..to {
            let byte = packed[(position / 4 - first) as usize];
            let code = (byte >> (6 - 2 * (position % 4))) & 3;
            bases.push(b"TCAG"[code as usize]);
        }
        let runs = |list: &[u32]| -> Vec<(u64, u64)> {
            let half = list.len() / 2;
            list[..half]
                .iter()
                .zip(&list[half..])
                .map(|(start, size)| (u64::from(*start), u64::from(*start) + u64::from(*size)))
                .filter(|(start, end)| *start < to && *end > from)
                .collect()
        };
        for (start, end) in runs(&unknown) {
            for base in &mut bases[(start.max(from) - from) as usize..(end.min(to) - from) as usize]
            {
                *base = b'N';
            }
        }
        for (start, end) in runs(&masked) {
            for base in &mut bases[(start.max(from) - from) as usize..(end.min(to) - from) as usize]
            {
                base.make_ascii_lowercase();
            }
        }
    }
    Ok(Bases {
        name,
        length,
        start: from,
        bases,
    })
}

/// A 2bit's index: each sequence's name and where its record starts.
struct Index {
    big: bool,
    entries: Vec<(String, u64)>,
}

/// Reads the header and the index, a stretch at a time: the index is read
/// once for every window, and an entry at a time it was two reads for each
/// of the hundred thousand contigs a draft assembly can have.
fn index<R: Read + Seek>(file: &mut File<R>) -> Result<Index, ReadError> {
    let head = file.read_up_to(0, 16)?;
    if head.len() < 16 {
        return Err(ReadError::whole("not a 2bit: it is too short"));
    }
    let first = u32::from_le_bytes([head[0], head[1], head[2], head[3]]);
    let big = if first == MAGIC {
        false
    } else if first.swap_bytes() == MAGIC {
        true
    } else {
        return Err(ReadError::whole(
            "not a 2bit: it does not start with a 2bit's signature",
        ));
    };
    let mut bytes = Bytes::new(&head[4..], big, "2bit");
    let version = bytes.u32()?;
    let long = match version {
        0 => false,
        1 => true,
        _ => {
            return Err(ReadError::whole(format!(
                "the 2bit is version {version}, and versions 0 and 1 are the ones there are"
            )))
        }
    };
    let count = bytes.u32()?;
    // An entry is a byte of length, a name and an offset: never more than
    // this, and never less than six bytes, which caps what the count asks
    // for at what the file has room for.
    const LONGEST: usize = 1 + 255 + 8;
    let mut entries = Vec::with_capacity((count as usize).min((file.size / 6) as usize));
    let mut chunk: Vec<u8> = Vec::new();
    let (mut chunk_at, mut within) = (16u64, 0usize);
    for _ in 0..count {
        if chunk.len() - within < LONGEST {
            chunk_at = chunk_at
                .checked_add(within as u64)
                .ok_or_else(|| super::bytes::short("2bit"))?;
            chunk = file.read_up_to(chunk_at, 1 << 16)?;
            within = 0;
        }
        let mut bytes = Bytes::new(&chunk[within..], big, "2bit");
        let length = usize::from(bytes.u8()?);
        let name = String::from_utf8_lossy(bytes.take(length)?).into_owned();
        let offset = if long {
            bytes.u64()?
        } else {
            u64::from(bytes.u32()?)
        };
        within = chunk.len() - bytes.left();
        entries.push((name, offset));
    }
    Ok(Index { big, entries })
}

/// A 2bit written by hand from FASTA records, as `faToTwoBit` writes one, in
/// either byte order and either version.
#[cfg(test)]
pub(crate) mod fixture {
    /// The runs of bases `inside` says are in, as starts then sizes.
    fn runs(bases: &[u8], inside: impl Fn(u8) -> bool) -> (Vec<u32>, Vec<u32>) {
        let (mut starts, mut sizes) = (Vec::new(), Vec::new());
        let mut at = 0;
        while at < bases.len() {
            if inside(bases[at]) {
                let start = at;
                while at < bases.len() && inside(bases[at]) {
                    at += 1;
                }
                starts.push(start as u32);
                sizes.push((at - start) as u32);
            } else {
                at += 1;
            }
        }
        (starts, sizes)
    }

    pub(crate) fn written(records: &[(&str, &[u8])], big: bool, long: bool) -> Vec<u8> {
        let u32b = |value: u32| {
            if big {
                value.to_be_bytes()
            } else {
                value.to_le_bytes()
            }
        };
        let u64b = |value: u64| {
            if big {
                value.to_be_bytes()
            } else {
                value.to_le_bytes()
            }
        };
        let mut out = Vec::new();
        out.extend(u32b(0x1A41_2743));
        out.extend(u32b(u32::from(long)));
        out.extend(u32b(records.len() as u32));
        out.extend(u32b(0));
        let index: usize = records
            .iter()
            .map(|(name, _)| 1 + name.len() + if long { 8 } else { 4 })
            .sum();
        let mut bodies = Vec::new();
        let mut at = (16 + index) as u64;
        for (name, bases) in records {
            out.push(name.len() as u8);
            out.extend(name.as_bytes());
            if long {
                out.extend(u64b(at));
            } else {
                out.extend(u32b(at as u32));
            }
            let mut body = Vec::new();
            body.extend(u32b(bases.len() as u32));
            for (starts, sizes) in [
                runs(bases, |base| base.eq_ignore_ascii_case(&b'N')),
                runs(bases, |base| base.is_ascii_lowercase()),
            ] {
                body.extend(u32b(starts.len() as u32));
                starts.iter().for_each(|start| body.extend(u32b(*start)));
                sizes.iter().for_each(|size| body.extend(u32b(*size)));
            }
            body.extend(u32b(0));
            for four in bases.chunks(4) {
                let mut byte = 0u8;
                for (place, base) in four.iter().enumerate() {
                    let code = match base.to_ascii_uppercase() {
                        b'C' => 1,
                        b'A' => 2,
                        b'G' => 3,
                        _ => 0,
                    };
                    byte |= code << (6 - 2 * place);
                }
                body.push(byte);
            }
            at += body.len() as u64;
            bodies.push(body);
        }
        for body in bodies {
            out.extend(body);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::written;
    use super::*;
    use std::io::Cursor;

    const TWOBIT: &[u8] = include_bytes!("fixtures/ref.2bit");
    const LONG: &[u8] = include_bytes!("fixtures/ref.long.2bit");

    /// The records of a FASTA, with its lines joined.
    fn records(text: &str) -> Vec<(String, Vec<u8>)> {
        let mut out: Vec<(String, Vec<u8>)> = Vec::new();
        for line in text.lines() {
            match line.strip_prefix('>') {
                Some(header) => out.push((
                    header.split_whitespace().next().unwrap().to_string(),
                    Vec::new(),
                )),
                None => out.last_mut().unwrap().1.extend(line.trim().bytes()),
            }
        }
        out
    }

    fn window(bytes: &[u8], locus: &str) -> Bases {
        bases(Cursor::new(bytes), &Region::parse(locus).unwrap()).unwrap()
    }

    /// Every sequence whole, and every window of every sequence, holds the
    /// bases `twoBitToFa` writes: N over the runs of N, lower case over the
    /// masked runs, and `n` where they meet.
    #[test]
    fn the_bases_are_what_twobittofa_writes() {
        let tool = records(include_str!("fixtures/ref.2bit.fa"));
        assert_eq!(tool.len(), 3);
        for bytes in [TWOBIT, LONG] {
            for (name, expected) in &tool {
                let whole = window(bytes, &format!("{name}:1-1000"));
                assert_eq!((whole.name.as_str(), whole.start), (name.as_str(), 0));
                assert_eq!(whole.length, expected.len() as u64);
                assert_eq!(whole.bases, *expected, "{name}");
                for from in 0..expected.len() {
                    for to in from + 1..=expected.len() {
                        let part = window(bytes, &format!("{name}:{}-{to}", from + 1));
                        assert_eq!(part.start, from as u64);
                        assert_eq!(part.bases, expected[from..to], "{name} {from}..{to}");
                    }
                }
            }
        }
        // And the window twoBitToFa -start=2 -end=50 wrote, under the header
        // it gives a 0-based start.
        let piece = records(include_str!("fixtures/ref.chr1.fa"));
        assert_eq!(window(TWOBIT, "chr1:3-50").bases, piece[0].1);
        assert!(piece[0].1.contains(&b'n') && piece[0].1.contains(&b'N'));
    }

    /// The writer here writes what faToTwoBit wrote, byte for byte, so the
    /// other byte order it also writes is the same file turned round.
    #[test]
    fn version_1_offsets_and_a_swapped_byte_order_read() {
        let fasta = records(include_str!("fixtures/ref.fa"));
        let borrowed: Vec<(&str, &[u8])> = fasta
            .iter()
            .map(|(name, bases)| (name.as_str(), bases.as_slice()))
            .collect();
        assert_eq!(written(&borrowed, false, false), TWOBIT);
        assert_eq!(written(&borrowed, false, true), LONG);
        for long in [false, true] {
            let swapped = written(&borrowed, true, long);
            assert_eq!(
                sequences(Cursor::new(&swapped)).unwrap(),
                sequences(Cursor::new(TWOBIT)).unwrap()
            );
            for locus in ["chr1:1-120", "chr1:30-47", "chr2:5-37", "chr3:1-12"] {
                assert_eq!(window(&swapped, locus), window(TWOBIT, locus), "{locus}");
            }
        }
    }

    /// An index longer than the stretch it is read in, as a draft assembly's
    /// is: three thousand names of thirty letters are 105 KB of index, read
    /// in stretches of 64 KB, and an entry cut by the end of one is read
    /// whole from the next.
    #[test]
    fn an_index_of_thousands_of_sequences_reads_across_its_stretches() {
        let records: Vec<(String, Vec<u8>)> = (0..3_000)
            .map(|at| {
                let bases = b"ACGTTGCA"[..1 + at % 8].to_vec();
                (format!("scaffold_{at:05}_of_a_draft_assembly"), bases)
            })
            .collect();
        let borrowed: Vec<(&str, &[u8])> = records
            .iter()
            .map(|(name, bases)| (name.as_str(), bases.as_slice()))
            .collect();
        let bytes = written(&borrowed, false, false);
        let held = sequences(Cursor::new(&bytes)).unwrap();
        assert_eq!(held.len(), 3_000);
        for ((name, length), (wanted, bases)) in held.iter().zip(&records) {
            assert_eq!((name, *length), (wanted, bases.len() as u64));
        }
        let last = window(&bytes, "scaffold_02999_of_a_draft_assembly:1-8");
        assert_eq!(last.bases, b"ACGTTGCA");
    }

    #[test]
    fn the_sequences_are_named_with_their_lengths() {
        let held = sequences(Cursor::new(TWOBIT)).unwrap();
        assert_eq!(
            held,
            [
                ("chr1".to_string(), 120),
                ("chr2".to_string(), 37),
                ("chr3".to_string(), 12)
            ]
        );
        assert_eq!(sequences(Cursor::new(LONG)).unwrap(), held);
    }

    /// A window that runs past the end holds the bases there are, and one
    /// that starts past it holds none, with the length to say why.
    #[test]
    fn a_window_past_the_end_is_clamped_and_one_beyond_it_is_empty() {
        let part = window(TWOBIT, "chr2:31-60");
        assert_eq!((part.start, part.bases.len(), part.length), (30, 7, 37));
        let none = window(TWOBIT, "chr2:38-60");
        assert_eq!((none.start, none.bases.len(), none.length), (37, 0, 37));
    }

    /// A 2bit of one sequence is that sequence whatever the window calls it,
    /// as a FASTA of one record is; of several, the name picks.
    #[test]
    fn one_sequence_is_that_sequence_whatever_the_region_calls_it() {
        let one = include_bytes!("fixtures/one.2bit");
        let fasta = records(include_str!("fixtures/one.fa"));
        let read = window(one, "chr1:1-10");
        assert_eq!(
            (read.name.as_str(), read.bases.as_slice()),
            ("contig_7", &fasta[0].1[..10])
        );
        let error = bases(Cursor::new(TWOBIT), &Region::parse("chr9:1-10").unwrap()).unwrap_err();
        assert_eq!(
            error.to_string(),
            "the 2bit has no sequence called chr9; it has chr1, chr2, chr3"
        );
        let twice = written(&[("chr1", b"ACGT"), ("chr1", b"TTTT")], false, false);
        let error = bases(Cursor::new(&twice), &Region::parse("chr1:1-4").unwrap()).unwrap_err();
        assert!(
            error.to_string().contains("2 sequences called chr1"),
            "{error}"
        );
    }

    #[test]
    fn a_damaged_file_is_an_error_and_never_a_panic() {
        let region = Region::parse("chr1:1-120").unwrap();
        for bytes in [TWOBIT, LONG] {
            for cut in 0..bytes.len() {
                // Cut after chr1's bases, chr1 still reads.
                let _ = bases(Cursor::new(&bytes[..cut]), &region);
                let _ = sequences(Cursor::new(&bytes[..cut]));
            }
            for at in 0..bytes.len() {
                for change in [0xff, 0x55] {
                    let mut bent = bytes.to_vec();
                    bent[at] ^= change;
                    let _ = bases(Cursor::new(&bent), &region);
                    let _ = sequences(Cursor::new(&bent));
                }
            }
        }
        // A count of runs as large as it goes asks for nothing it cannot see.
        let mut bent = TWOBIT.to_vec();
        let record = u32::from_le_bytes([bent[21], bent[22], bent[23], bent[24]]) as usize;
        bent[record + 4..record + 8].copy_from_slice(&u32::MAX.to_le_bytes());
        let error = bases(Cursor::new(&bent), &region).unwrap_err();
        assert!(error.to_string().contains("past its own end"), "{error}");
        let error = bases(Cursor::new(&b"not a 2bit at all"[..]), &region).unwrap_err();
        assert!(error.to_string().contains("not a 2bit"), "{error}");
    }
}
