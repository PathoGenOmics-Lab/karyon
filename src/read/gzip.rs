//! gzip and BGZF, decompressed with no dependency.
//!
//! Nearly every file a genomics pipeline leaves behind is compressed: a
//! `.vcf.gz`, a `.gff3.gz`, a `.bed.gz`, and a BAM, which is BGZF underneath.
//! Each of them is text, or a binary format, inside the same wrapper, so a
//! reader that takes the wrapper off turns every one of those files into one
//! karyon already reads. It is written here rather than taken from a crate
//! because the crate's first promise is that `cargo add karyon` brings nothing
//! else with it.
//!
//! # What is here
//!
//! A DEFLATE decoder (RFC 1951) in the canonical-Huffman style of zlib's
//! `puff`, with a lookup table for the short codes that make up most of a
//! stream, and the gzip wrapper around it (RFC 1952): the header, the CRC-32
//! and the length every member ends with, both checked. BGZF, which is what
//! bgzip, samtools and tabix write, is a series of ordinary gzip members, so a
//! reader of every member is a BGZF reader; [`member`] reads one, for a reader
//! that seeks to a block through an index.
//!
//! # What is refused
//!
//! Anything that is not what it says: a member that does not start with the
//! gzip magic, a stream that ends early, a code that does not decode, a
//! distance reaching before the start of the output, and a member whose CRC or
//! length does not match what it decompressed to. A damaged file is an error
//! naming what was wrong, never a panic and never a figure drawn from half of
//! it. The padding of zeros some writers leave after the last member is not
//! damage, and is read past.

use super::ReadError;

/// Whether `data` starts the way every gzip member does.
pub fn is_gzip(data: &[u8]) -> bool {
    data.len() >= 3 && data[0] == 0x1f && data[1] == 0x8b && data[2] == 8
}

/// Decompresses every member in `data`: plain gzip, a concatenation of gzip
/// files, or BGZF.
///
/// # Errors
///
/// A member that is not gzip, is cut short, does not decode, or does not match
/// its own CRC-32 or length.
pub fn decompress(data: &[u8]) -> Result<Vec<u8>, ReadError> {
    let mut out = Vec::with_capacity(data.len().saturating_mul(4));
    let mut at = 0;
    let mut members = 0usize;
    while at < data.len() {
        // Some writers pad the file with zeros after the last member.
        if data[at..].iter().all(|byte| *byte == 0) {
            break;
        }
        at += member_into(&data[at..], &mut out)
            .map_err(|reason| ReadError::whole(format!("gzip member {}: {reason}", members + 1)))?;
        members += 1;
    }
    Ok(out)
}

/// Decompresses the one member at the start of `data`, and says how many of
/// its bytes that member took.
///
/// For BGZF read through an index, which seeks to a block and wants that
/// block and no more.
///
/// # Errors
///
/// The same as [`decompress`], for this member.
pub fn member(data: &[u8]) -> Result<(Vec<u8>, usize), ReadError> {
    let mut out = Vec::new();
    let used = member_into(data, &mut out).map_err(ReadError::whole)?;
    Ok((out, used))
}

/// One member appended to `out`, and the bytes of `data` it took.
fn member_into(data: &[u8], out: &mut Vec<u8>) -> Result<usize, String> {
    if !is_gzip(data) {
        return Err("not gzip: it does not start with the gzip magic".to_string());
    }
    let flags = *data.get(3).ok_or("the header is cut short")?;
    // Magic, method, flags, time, extra flags and system: ten bytes.
    let mut at = 10;
    let byte = |at: usize| data.get(at).copied().ok_or("the header is cut short");
    if flags & 0x04 != 0 {
        let length = usize::from(byte(at)?) | usize::from(byte(at + 1)?) << 8;
        at += 2 + length;
    }
    for name in [0x08, 0x10] {
        // A file name, then a comment, each ended by a nought.
        if flags & name != 0 {
            while byte(at)? != 0 {
                at += 1;
            }
            at += 1;
        }
    }
    if flags & 0x02 != 0 {
        at += 2;
    }
    let body = data.get(at..).ok_or("the header is cut short")?;
    let before = out.len();
    let used = inflate(body, out)?;
    at += used;
    let trailer = data
        .get(at..at + 8)
        .ok_or("the stream ends before its CRC and length")?;
    let crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
    let written = out.len() - before;
    if crc32(&out[before..]) != crc {
        return Err("the CRC-32 does not match what it decompressed to".to_string());
    }
    // The length is kept modulo 2^32, which is what the format says.
    if written as u32 != length {
        return Err("the length does not match what it decompressed to".to_string());
    }
    Ok(at + 8)
}

/// The bits of a DEFLATE stream, least significant first.
struct Bits<'a> {
    data: &'a [u8],
    at: usize,
    held: u32,
    count: u32,
}

impl<'a> Bits<'a> {
    fn new(data: &'a [u8]) -> Self {
        Bits {
            data,
            at: 0,
            held: 0,
            count: 0,
        }
    }

    /// The next `n` bits, `n` at most sixteen.
    fn take(&mut self, n: u32) -> Result<u32, String> {
        while self.count < n {
            let byte = *self
                .data
                .get(self.at)
                .ok_or("the compressed stream ends early")?;
            self.at += 1;
            self.held |= u32::from(byte) << self.count;
            self.count += 8;
        }
        let value = self.held & ((1u32 << n) - 1);
        self.held >>= n;
        self.count -= n;
        Ok(value)
    }

    /// Drops what is left of the current byte, for a stored block, and gives
    /// back the whole bytes the lookahead took.
    fn align(&mut self) {
        self.at -= (self.count / 8) as usize;
        self.held = 0;
        self.count = 0;
    }

    /// Where the stream ended: the bytes read, less the whole bytes of
    /// lookahead still held.
    fn consumed(&self) -> usize {
        self.at - (self.count / 8) as usize
    }
}

/// How many bits the lookup table decodes at once.
const FAST: u32 = 9;

/// A canonical Huffman code, as DEFLATE describes one: a length per symbol.
struct Huffman {
    /// How many codes there are of each length.
    count: [u16; 16],
    /// The symbols in the order of their codes.
    symbol: Vec<u16>,
    /// The codes of `FAST` bits or fewer, looked up by their next `FAST` bits
    /// as the stream holds them: the symbol and the code's length, or a length
    /// of nought for a longer code, which is decoded a bit at a time.
    fast: Vec<(u16, u8)>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Result<Huffman, String> {
        let mut count = [0u16; 16];
        for length in lengths {
            count[usize::from(*length)] += 1;
        }
        // A code that promises more codes than its lengths can hold does not
        // decode. One that holds fewer is allowed: a stream may use one
        // distance, and DEFLATE writes that as a code of one.
        let mut left: i32 = 1;
        for codes in &count[1..] {
            left <<= 1;
            left -= i32::from(*codes);
            if left < 0 {
                return Err("a Huffman code is over-subscribed".to_string());
            }
        }
        let mut offsets = [0u16; 16];
        for length in 1..15 {
            offsets[length + 1] = offsets[length] + count[length];
        }
        let mut symbol = vec![0u16; lengths.len()];
        for (value, length) in lengths.iter().enumerate() {
            if *length != 0 {
                let slot = &mut offsets[usize::from(*length)];
                symbol[usize::from(*slot)] = value as u16;
                *slot += 1;
            }
        }

        let mut fast = vec![(0u16, 0u8); 1 << FAST];
        let mut code: u32 = 0;
        let mut index = 0usize;
        for length in 1..=15u32 {
            for _ in 0..count[length as usize] {
                if length <= FAST {
                    // The code is written most significant bit first and read
                    // least significant first, so the table is indexed by it
                    // reversed, once for every value of the bits after it.
                    let mut reversed = 0u32;
                    for bit in 0..length {
                        reversed |= ((code >> bit) & 1) << (length - 1 - bit);
                    }
                    let mut slot = reversed;
                    while slot < (1 << FAST) {
                        fast[slot as usize] = (symbol[index], length as u8);
                        slot += 1 << length;
                    }
                }
                code += 1;
                index += 1;
            }
            code <<= 1;
        }
        Ok(Huffman {
            count,
            symbol,
            fast,
        })
    }

    fn decode(&self, bits: &mut Bits<'_>) -> Result<u16, String> {
        while bits.count < FAST && bits.at < bits.data.len() {
            bits.held |= u32::from(bits.data[bits.at]) << bits.count;
            bits.at += 1;
            bits.count += 8;
        }
        if bits.count >= FAST {
            let (value, length) = self.fast[(bits.held & ((1 << FAST) - 1)) as usize];
            if length != 0 {
                bits.held >>= length;
                bits.count -= u32::from(length);
                return Ok(value);
            }
        }
        // A long code, or the last few bits of the stream: a bit at a time.
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for length in 1..16 {
            code |= bits.take(1)? as i32;
            let count = i32::from(self.count[length]);
            if code - count < first {
                return Ok(self.symbol[(index + (code - first)) as usize]);
            }
            index += count;
            first += count;
            first <<= 1;
            code <<= 1;
        }
        Err("a Huffman code does not decode".to_string())
    }
}

const LENGTH_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LENGTH_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DISTANCE_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DISTANCE_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

/// One block's literals and copies, until its end code.
fn codes(
    bits: &mut Bits<'_>,
    out: &mut Vec<u8>,
    start: usize,
    literal: &Huffman,
    distance: &Huffman,
) -> Result<(), String> {
    loop {
        let value = literal.decode(bits)?;
        match value {
            0..=255 => out.push(value as u8),
            256 => return Ok(()),
            _ => {
                let at = usize::from(value - 257);
                if at >= LENGTH_BASE.len() {
                    return Err("a length code is out of range".to_string());
                }
                let length =
                    usize::from(LENGTH_BASE[at]) + bits.take(u32::from(LENGTH_EXTRA[at]))? as usize;
                let at = usize::from(distance.decode(bits)?);
                if at >= DISTANCE_BASE.len() {
                    return Err("a distance code is out of range".to_string());
                }
                let back = usize::from(DISTANCE_BASE[at])
                    + bits.take(u32::from(DISTANCE_EXTRA[at]))? as usize;
                // Only as far back as this member's own output: a member does
                // not reach into the one before it.
                if back > out.len() - start {
                    return Err("a distance reaches before the start of the output".to_string());
                }
                // One byte at a time, since a copy may overlap what it writes,
                // which is how DEFLATE writes a run.
                let from = out.len() - back;
                for offset in 0..length {
                    let byte = out[from + offset];
                    out.push(byte);
                }
            }
        }
    }
}

/// Inflates one raw DEFLATE stream onto `out`, and says how many bytes of
/// `data` it took.
fn inflate(data: &[u8], out: &mut Vec<u8>) -> Result<usize, String> {
    let start = out.len();
    let mut bits = Bits::new(data);
    loop {
        let last = bits.take(1)?;
        match bits.take(2)? {
            0 => {
                bits.align();
                let at = bits.at;
                let header = data.get(at..at + 4).ok_or("a stored block is cut short")?;
                let length = usize::from(u16::from_le_bytes([header[0], header[1]]));
                let check = u16::from_le_bytes([header[2], header[3]]);
                if length != usize::from(!check) {
                    return Err("a stored block's length does not match its check".to_string());
                }
                let stored = data
                    .get(at + 4..at + 4 + length)
                    .ok_or("a stored block is cut short")?;
                out.extend_from_slice(stored);
                bits.at = at + 4 + length;
            }
            1 => {
                let mut lengths = [0u8; 288];
                lengths[..144].fill(8);
                lengths[144..256].fill(9);
                lengths[256..280].fill(7);
                lengths[280..].fill(8);
                let literal = Huffman::new(&lengths)?;
                let distance = Huffman::new(&[5u8; 30])?;
                codes(&mut bits, out, start, &literal, &distance)?;
            }
            2 => {
                let literals = bits.take(5)? as usize + 257;
                let distances = bits.take(5)? as usize + 1;
                let lengths_of_lengths = bits.take(4)? as usize + 4;
                const ORDER: [usize; 19] = [
                    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
                ];
                let mut code_lengths = [0u8; 19];
                for at in ORDER.iter().take(lengths_of_lengths) {
                    code_lengths[*at] = bits.take(3)? as u8;
                }
                let lengths_code = Huffman::new(&code_lengths)?;
                let mut lengths = vec![0u8; literals + distances];
                let mut at = 0;
                while at < lengths.len() {
                    let value = lengths_code.decode(&mut bits)?;
                    if value < 16 {
                        lengths[at] = value as u8;
                        at += 1;
                        continue;
                    }
                    let (repeated, times) = match value {
                        16 => {
                            let previous = *lengths
                                .get(at.wrapping_sub(1))
                                .ok_or("a repeat has no length before it")?;
                            (previous, 3 + bits.take(2)? as usize)
                        }
                        17 => (0, 3 + bits.take(3)? as usize),
                        _ => (0, 11 + bits.take(7)? as usize),
                    };
                    if at + times > lengths.len() {
                        return Err("the code lengths run past their count".to_string());
                    }
                    lengths[at..at + times].fill(repeated);
                    at += times;
                }
                if lengths[256] == 0 {
                    return Err("a block has no end code".to_string());
                }
                let literal = Huffman::new(&lengths[..literals])?;
                let distance = Huffman::new(&lengths[literals..])?;
                codes(&mut bits, out, start, &literal, &distance)?;
            }
            _ => return Err("a block has a type DEFLATE does not define".to_string()),
        }
        if last == 1 {
            return Ok(bits.consumed());
        }
    }
}

/// The CRC-32 table gzip uses, computed once at compile time.
const CRC_TABLE: [u32; 256] = {
    let mut table = [0u32; 256];
    let mut at = 0;
    while at < 256 {
        let mut crc = at as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = if crc & 1 != 0 {
                0xedb8_8320 ^ (crc >> 1)
            } else {
                crc >> 1
            };
            bit += 1;
        }
        table[at] = crc;
        at += 1;
    }
    table
};

fn crc32(data: &[u8]) -> u32 {
    let mut crc = !0u32;
    for byte in data {
        crc = CRC_TABLE[((crc ^ u32::from(*byte)) & 0xff) as usize] ^ (crc >> 8);
    }
    !crc
}

#[cfg(test)]
mod tests {
    use super::*;

    /// "karyon\n" as `gzip -n` writes it: fixed Huffman codes.
    const SMALL: [u8; 27] = [
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xcb, 0x4e, 0x2c, 0xaa, 0xcc,
        0xcf, 0xe3, 0x02, 0x00, 0x0d, 0x64, 0x70, 0x1a, 0x07, 0x00, 0x00, 0x00,
    ];

    const DYNAMIC: [u8; 451] = [
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x03, 0x55, 0xd4, 0x3b, 0x6e, 0xc3,
        0x30, 0x10, 0x04, 0xd0, 0xda, 0xb9, 0x4a, 0x10, 0x40, 0xfb, 0x93, 0xc4, 0xf3, 0x18, 0x46,
        0x52, 0xa5, 0xc8, 0xfd, 0x8b, 0x70, 0x68, 0x8e, 0x38, 0xea, 0x8c, 0x91, 0x4c, 0x60, 0x1e,
        0x77, 0xf5, 0xfc, 0xf9, 0xb3, 0xc7, 0xf6, 0x68, 0xdb, 0xe3, 0xfb, 0xf5, 0xfb, 0xda, 0xfa,
        0xcf, 0xaf, 0x8f, 0x27, 0x32, 0xdb, 0xb6, 0x87, 0xcd, 0x18, 0x6f, 0x7c, 0xbe, 0x63, 0xef,
        0xb1, 0xcf, 0xd8, 0x57, 0x1c, 0x3d, 0x8e, 0x19, 0xc7, 0x3a, 0x24, 0x7b, 0x9c, 0x33, 0xce,
        0xf5, 0x76, 0xf5, 0xb8, 0x66, 0x5c, 0x2b, 0xde, 0x7b, 0xbc, 0xcf, 0x78, 0x5f, 0x87, 0x1c,
        0x3d, 0x3e, 0x66, 0x7c, 0xac, 0xb7, 0xcf, 0x1e, 0x9f, 0x33, 0x3e, 0x57, 0xdc, 0x7a, 0xdc,
        0x66, 0xdc, 0x6e, 0x75, 0x7a, 0x9f, 0x8d, 0x85, 0xb6, 0xf5, 0x07, 0x1b, 0x4d, 0xaf, 0xaa,
        0xd2, 0xd5, 0x50, 0xd6, 0xd8, 0xd6, 0x5c, 0x4e, 0x43, 0x5f, 0x63, 0x61, 0x0b, 0xf9, 0x0f,
        0x2a, 0x1b, 0x3b, 0x9b, 0x94, 0x36, 0xb4, 0x36, 0xd6, 0xb6, 0x92, 0xd3, 0x50, 0xdc, 0xd8,
        0xdc, 0x76, 0xf9, 0x0f, 0xba, 0x1b, 0xcb, 0x9b, 0xb4, 0x37, 0xd4, 0x37, 0xf6, 0xb7, 0x53,
        0x4e, 0x6b, 0xe3, 0xe6, 0xf8, 0xa4, 0xdd, 0xee, 0xae, 0x5f, 0x1e, 0x0d, 0x5c, 0x0c, 0x1c,
        0x06, 0x4e, 0x03, 0xb7, 0x75, 0x9a, 0x8f, 0x0b, 0xbf, 0x6e, 0x5c, 0xae, 0xdc, 0x61, 0xe0,
        0x34, 0x70, 0x31, 0x70, 0x18, 0x38, 0x0d, 0x3c, 0xe5, 0x34, 0x18, 0x38, 0x0d, 0x5c, 0xee,
        0xde, 0x61, 0xe0, 0x34, 0x70, 0x31, 0x70, 0x18, 0x38, 0x0d, 0xfc, 0x90, 0xd3, 0x60, 0xe0,
        0x34, 0x70, 0x19, 0x02, 0x6f, 0x63, 0x4c, 0xf9, 0xa4, 0xdd, 0x06, 0xb5, 0x4f, 0x2a, 0x0d,
        0x42, 0x06, 0x3e, 0x60, 0x10, 0x34, 0x08, 0x99, 0x83, 0x80, 0x41, 0xd0, 0x20, 0x74, 0xec,
        0xc7, 0xdc, 0x5f, 0x83, 0x2f, 0x93, 0x1f, 0x30, 0x08, 0x1a, 0x84, 0xcc, 0x41, 0xc0, 0x20,
        0x68, 0x10, 0x62, 0x10, 0x30, 0x08, 0x1a, 0x84, 0xac, 0x40, 0xc0, 0x20, 0x68, 0x10, 0x32,
        0x07, 0x01, 0x83, 0xa0, 0x41, 0x88, 0x41, 0xb4, 0xb1, 0x93, 0x7c, 0xd2, 0x6e, 0x5b, 0xd9,
        0xd7, 0x92, 0x06, 0x29, 0x73, 0x90, 0x30, 0x48, 0x1a, 0xa4, 0x18, 0x24, 0x0c, 0x92, 0x06,
        0x29, 0xbb, 0x90, 0x30, 0x48, 0x1a, 0xa4, 0xcc, 0x41, 0x8e, 0xf5, 0xbf, 0xf6, 0x5f, 0x0c,
        0x12, 0x06, 0x49, 0x83, 0x94, 0x5d, 0x48, 0x18, 0x24, 0x0d, 0x52, 0xe6, 0x20, 0x61, 0x90,
        0x34, 0x48, 0x31, 0x48, 0x18, 0x24, 0x0d, 0x52, 0x76, 0x21, 0xdb, 0xf8, 0x00, 0xf1, 0x49,
        0xbb, 0x7d, 0x82, 0xfa, 0x37, 0x88, 0x06, 0x25, 0x06, 0x05, 0x83, 0xa2, 0x41, 0xc9, 0x2e,
        0x14, 0x0c, 0x8a, 0x06, 0x25, 0x73, 0x50, 0x30, 0x28, 0x1a, 0x94, 0x18, 0x14, 0x0c, 0x8a,
        0x06, 0x25, 0xbb, 0x50, 0xe3, 0x2b, 0x78, 0x7d, 0x06, 0x65, 0x0e, 0x0a, 0x06, 0x45, 0x83,
        0x12, 0x83, 0x82, 0x41, 0xd1, 0xa0, 0x64, 0x17, 0x0a, 0x06, 0x45, 0x83, 0x92, 0x39, 0xa8,
        0x36, 0xbe, 0xb6, 0x7c, 0xf2, 0x36, 0xf8, 0x07, 0xa4, 0x63, 0xc3, 0xed, 0xf7, 0x05, 0x00,
        0x00,
    ];

    const STORED: [u8; 423] = [
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x01, 0x90, 0x01, 0x6f, 0xfe,
        0x52, 0xf2, 0x26, 0x65, 0xa6, 0x0c, 0x12, 0xd2, 0x89, 0x18, 0x5d, 0x95, 0x0e, 0xe8, 0x81,
        0x36, 0x09, 0x16, 0x6f, 0x6b, 0x11, 0x3d, 0x17, 0x8d, 0x6c, 0x0f, 0xd3, 0x90, 0x1f, 0xf2,
        0x39, 0xa1, 0xa0, 0x95, 0xf2, 0x0f, 0x93, 0x95, 0x65, 0x0c, 0xf9, 0x38, 0x0b, 0x8e, 0xdb,
        0x22, 0x4a, 0x6b, 0x24, 0x8a, 0x1e, 0x92, 0x4e, 0x8f, 0xd0, 0xae, 0x2e, 0x1a, 0x94, 0x92,
        0xa3, 0x30, 0x5f, 0x18, 0x8c, 0xb6, 0x10, 0x90, 0x0f, 0x9e, 0x34, 0x7f, 0xae, 0x88, 0x6d,
        0xc6, 0x50, 0x77, 0x95, 0xec, 0x74, 0x5c, 0x4c, 0x3f, 0xcb, 0x2e, 0xb2, 0xc7, 0x3e, 0x14,
        0x93, 0x4c, 0x86, 0x7e, 0xe0, 0x57, 0xba, 0x72, 0x49, 0x9b, 0xfa, 0x12, 0x1e, 0x83, 0x6b,
        0x2a, 0xc1, 0x57, 0x26, 0xee, 0x7d, 0x6b, 0x0a, 0xf6, 0xab, 0x13, 0xc3, 0x8e, 0x92, 0xca,
        0xe0, 0xd1, 0x50, 0x57, 0xb1, 0x59, 0x98, 0x7f, 0x94, 0xcc, 0x74, 0x11, 0xd7, 0x17, 0xf1,
        0x45, 0x79, 0xb2, 0xaa, 0x10, 0x0f, 0xbb, 0xb3, 0x4f, 0xa5, 0x93, 0xfe, 0xae, 0xd2, 0x72,
        0x48, 0xb7, 0x62, 0xe3, 0xab, 0x58, 0x05, 0xf0, 0x76, 0x5a, 0x2b, 0x9c, 0x1d, 0x7e, 0x0f,
        0x37, 0xc4, 0x49, 0x21, 0xbd, 0x3f, 0x65, 0x64, 0xea, 0xdf, 0x7f, 0x14, 0x2a, 0x72, 0x66,
        0x8c, 0x47, 0xe2, 0x23, 0xd1, 0x6e, 0xdd, 0x8c, 0x47, 0xb4, 0x6a, 0xfc, 0x5b, 0xae, 0xe2,
        0x61, 0xf5, 0x3b, 0x26, 0x15, 0x2d, 0x26, 0x3b, 0xa8, 0x3b, 0x03, 0x7c, 0xd4, 0x96, 0x2e,
        0x43, 0x48, 0x01, 0x25, 0x6b, 0x88, 0x5e, 0x9c, 0x90, 0x51, 0xf3, 0x20, 0xb0, 0xdb, 0x83,
        0xf3, 0x9e, 0xa7, 0xad, 0xbd, 0x0d, 0x74, 0xe6, 0xde, 0xc7, 0xf3, 0xdf, 0xae, 0xcc, 0x8f,
        0x64, 0x65, 0x66, 0x64, 0x1a, 0x7b, 0xa2, 0x66, 0x0f, 0x30, 0x11, 0xfc, 0x35, 0x70, 0x29,
        0x1c, 0x57, 0x99, 0x0d, 0x1a, 0x00, 0x91, 0x26, 0x89, 0x19, 0xf2, 0x5d, 0x9d, 0x06, 0x12,
        0xdf, 0x35, 0x9d, 0x60, 0x26, 0xa2, 0x40, 0xf4, 0x58, 0x9a, 0x5d, 0x79, 0x1f, 0x1d, 0xd9,
        0x7c, 0xfe, 0xfa, 0x77, 0x7a, 0x7b, 0x4f, 0x15, 0x24, 0x1a, 0xbf, 0x57, 0xbd, 0x43, 0x7a,
        0xd4, 0xb1, 0x29, 0x84, 0x05, 0x34, 0xf3, 0xf3, 0x87, 0x5c, 0x25, 0xb0, 0x8b, 0xea, 0x06,
        0xc2, 0x87, 0x4c, 0xfa, 0xa4, 0xdd, 0x17, 0xb2, 0xd8, 0x42, 0x84, 0x5d, 0xe8, 0x2a, 0x5b,
        0xc5, 0x39, 0x88, 0x8a, 0xc7, 0x80, 0x54, 0xa2, 0x39, 0x9c, 0xcf, 0xc9, 0xfc, 0xc2, 0xda,
        0x31, 0xce, 0x3d, 0xd1, 0x66, 0xbd, 0xcd, 0x3a, 0x33, 0x84, 0x7e, 0x5b, 0xbb, 0x07, 0xfd,
        0x07, 0xca, 0x47, 0x78, 0x42, 0x31, 0xb1, 0x9a, 0xf4, 0x58, 0x72, 0xce, 0xef, 0xb9, 0xfc,
        0x59, 0xf4, 0xf9, 0x5d, 0x14, 0x38, 0x1a, 0x3a, 0x78, 0x32, 0x56, 0x34, 0x7b, 0x9f, 0xfc,
        0xe6, 0x9c, 0xd7, 0x00, 0x7a, 0xe8, 0xa7, 0x58, 0xcc, 0xa4, 0x21, 0x20, 0xe8, 0x7e, 0x90,
        0x01, 0x00, 0x00,
    ];

    /// The text `DYNAMIC` holds, which `gzip -9` wrote with a dynamic block.
    fn dynamic_text() -> Vec<u8> {
        let lines: Vec<String> = (0..60)
            .map(|i| {
                let strand = if i % 3 == 0 { '-' } else { '+' };
                format!("chr1\t{}\t{}\tgene{i}\t0\t{strand}", i * 100, i * 100 + 90)
            })
            .collect();
        (lines.join("\n") + "\n").into_bytes()
    }

    #[test]
    fn all_three_kinds_of_block_decode() {
        // Fixed codes in SMALL, a dynamic block here, and a stored block, which
        // is what gzip writes for bytes it cannot shrink.
        assert_eq!(decompress(&DYNAMIC).unwrap(), dynamic_text());
        let stored = decompress(&STORED).unwrap();
        assert_eq!(stored.len(), 400);
        assert_eq!(stored, STORED[15..415].to_vec());
    }

    #[test]
    fn a_small_member_decompresses_to_what_was_compressed() {
        assert_eq!(decompress(&SMALL).unwrap(), b"karyon\n");
        let (out, used) = member(&SMALL).unwrap();
        assert_eq!(
            (out.as_slice(), used),
            (b"karyon\n".as_slice(), SMALL.len())
        );
    }

    #[test]
    fn several_members_are_one_file_as_bgzf_writes_it() {
        let mut two = SMALL.to_vec();
        two.extend_from_slice(&SMALL);
        // And the zeros some writers pad the end with.
        two.extend_from_slice(&[0; 5]);
        assert_eq!(decompress(&two).unwrap(), b"karyon\nkaryon\n");
    }

    #[test]
    fn a_damaged_file_is_an_error_and_never_a_panic() {
        // Cut at every length, and every byte changed in turn: each is either
        // the right text or an error, and nothing panics.
        for cut in 0..SMALL.len() {
            assert!(
                decompress(&SMALL[..cut]).is_err() || cut == 0,
                "cut at {cut}"
            );
        }
        for at in 0..SMALL.len() {
            let mut bent = SMALL;
            bent[at] ^= 0x55;
            match decompress(&bent) {
                Ok(text) => assert_eq!(text, b"karyon\n", "byte {at} changed the text"),
                Err(error) => assert!(error.to_string().contains("gzip member 1"), "{error}"),
            }
        }
        let error = decompress(b"not compressed").unwrap_err();
        assert!(error.to_string().contains("gzip magic"), "{error}");
    }

    #[test]
    fn the_crc_catches_what_decodes_and_is_wrong() {
        let mut wrong = SMALL;
        wrong[19] ^= 1;
        let error = decompress(&wrong).unwrap_err();
        assert!(error.to_string().contains("CRC-32"), "{error}");
    }
}
