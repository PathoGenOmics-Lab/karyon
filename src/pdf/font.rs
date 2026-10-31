//! TrueType: reading the faces a PDF embeds, and cutting each down to the
//! glyphs one page uses.
//!
//! A face is three to four hundred kilobytes, and a figure uses fifty of its
//! two and a half thousand glyphs, so what goes into the file is a subset:
//! the same font with every other glyph's outline taken out. Glyph numbers
//! are kept as they are, so a glyph is the same number in the face, in the
//! subset and in the content stream, and the PDF maps one to the other with
//! `/CIDToGIDMap /Identity`; the tables indexed by glyph are cut short after
//! the highest one used, and a glyph below it that is not used is there with
//! no outline. A glyph built from others, as an accented letter is from its
//! letter and its accent, brings those others with it.
//!
//! What a PDF reader needs of a TrueType font is its outlines, `glyf` and
//! `loca`, the metrics that size them, `head`, `hhea`, `hmtx` and `maxp`, and
//! the hinting programs, `cvt `, `fpgm` and `prep`, that keep small text
//! crisp on a screen. The rest of the file is for a system that installs the
//! font and a layout engine that shapes text with it, and is left behind:
//! kerning, ligatures, glyph names and the full character map.

use std::collections::BTreeSet;

use crate::svg::num;

/// A TrueType face, read from the bytes of its file.
#[derive(Debug)]
pub(crate) struct Font {
    data: &'static [u8],
    /// Each table's tag, and where it is in `data`.
    tables: Vec<([u8; 4], std::ops::Range<usize>)>,
    /// The PostScript name, from the `name` table: `LiberationSans-Bold`.
    pub(crate) name: String,
    units_per_em: u32,
    glyph_count: u16,
    /// How many glyphs have an advance of their own in `hmtx`; the rest
    /// share the last one's.
    metric_count: u16,
    long_loca: bool,
    /// Where the Unicode character map, a format 4 subtable, starts; a
    /// subset has none, since a PDF reads its glyphs by number.
    cmap: Option<usize>,
}

impl Font {
    /// Reads the face in `data`, or `None` when it is not TrueType outlines;
    /// every face this crate carries is.
    pub(crate) fn read(data: &'static [u8]) -> Option<Font> {
        if u32_at(data, 0) != 0x0001_0000 {
            return None;
        }
        let count = usize::from(u16_at(data, 4));
        let mut tables = Vec::with_capacity(count);
        for at in 0..count {
            let record = 12 + 16 * at;
            let tag: [u8; 4] = data.get(record..record + 4)?.try_into().ok()?;
            let offset = u32_at(data, record + 8) as usize;
            let length = u32_at(data, record + 12) as usize;
            if offset.checked_add(length)? > data.len() {
                return None;
            }
            tables.push((tag, offset..offset + length));
        }
        let mut font = Font {
            data,
            tables,
            name: String::new(),
            units_per_em: 0,
            glyph_count: 0,
            metric_count: 0,
            long_loca: false,
            cmap: None,
        };
        for tag in [*b"head", *b"hhea", *b"hmtx", *b"maxp", *b"loca", *b"glyf"] {
            font.table(tag)?;
        }
        let head = font.table(*b"head")?.start;
        font.units_per_em = u32::from(u16_at(data, head + 18)).max(1);
        font.long_loca = u16_at(data, head + 50) == 1;
        font.glyph_count = u16_at(data, font.table(*b"maxp")?.start + 4);
        font.metric_count = u16_at(data, font.table(*b"hhea")?.start + 34).max(1);
        font.cmap = font.unicode_map();
        font.name = font.postscript_name().unwrap_or_default();
        Some(font)
    }

    fn table(&self, tag: [u8; 4]) -> Option<std::ops::Range<usize>> {
        self.tables
            .iter()
            .find(|(known, _)| *known == tag)
            .map(|(_, range)| range.clone())
    }

    /// The format 4 subtable that maps Unicode's first plane to glyphs:
    /// Windows' own, or Unicode's where Windows' is missing.
    fn unicode_map(&self) -> Option<usize> {
        let cmap = self.table(*b"cmap")?.start;
        let count = usize::from(u16_at(self.data, cmap + 2));
        let subtables: Vec<(u16, u16, usize)> = (0..count)
            .map(|at| {
                let record = cmap + 4 + 8 * at;
                let offset = cmap + u32_at(self.data, record + 4) as usize;
                (
                    u16_at(self.data, record),
                    u16_at(self.data, record + 2),
                    offset,
                )
            })
            .filter(|(_, _, offset)| u16_at(self.data, *offset) == 4)
            .collect();
        [(3, 1), (0, 3)].iter().find_map(|(platform, encoding)| {
            subtables
                .iter()
                .find(|(p, e, _)| p == platform && e == encoding)
                .map(|(_, _, offset)| *offset)
        })
    }

    /// Name 6, the PostScript name, which is ASCII in every record it has.
    fn postscript_name(&self) -> Option<String> {
        let name = self.table(*b"name")?.start;
        let count = usize::from(u16_at(self.data, name + 2));
        let strings = name + usize::from(u16_at(self.data, name + 4));
        (0..count).find_map(|at| {
            let record = name + 6 + 12 * at;
            if u16_at(self.data, record + 6) != 6 {
                return None;
            }
            let platform = u16_at(self.data, record);
            let length = usize::from(u16_at(self.data, record + 8));
            let start = strings + usize::from(u16_at(self.data, record + 10));
            let bytes = self.data.get(start..start + length)?;
            // Windows writes it in UTF-16, every other byte nought.
            let ascii: String = match platform {
                3 => bytes
                    .chunks(2)
                    .filter_map(|pair| pair.get(1))
                    .map(|byte| char::from(*byte))
                    .collect(),
                _ => bytes.iter().map(|byte| char::from(*byte)).collect(),
            };
            Some(ascii).filter(|ascii| !ascii.is_empty())
        })
    }

    /// The glyph that draws `c`, if the face has one.
    pub(crate) fn glyph(&self, c: char) -> Option<u16> {
        let code = u16::try_from(u32::from(c)).ok()?;
        let map = self.cmap?;
        let segments = usize::from(u16_at(self.data, map + 6) / 2);
        let ends = map + 14;
        let starts = ends + 2 * segments + 2;
        let deltas = starts + 2 * segments;
        let ranges = deltas + 2 * segments;
        // The segments are sorted by their last character, so the first
        // that ends at or after `code` is the only one that can hold it.
        let (mut low, mut high) = (0, segments);
        while low < high {
            let middle = (low + high) / 2;
            if u16_at(self.data, ends + 2 * middle) < code {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        if low == segments || u16_at(self.data, starts + 2 * low) > code {
            return None;
        }
        let delta = u16_at(self.data, deltas + 2 * low);
        let range = usize::from(u16_at(self.data, ranges + 2 * low));
        let glyph = if range == 0 {
            code.wrapping_add(delta)
        } else {
            let start = u16_at(self.data, starts + 2 * low);
            let at = ranges + 2 * low + range + 2 * usize::from(code - start);
            match u16_at(self.data, at) {
                0 => 0,
                glyph => glyph.wrapping_add(delta),
            }
        };
        Some(glyph).filter(|glyph| *glyph != 0 && *glyph < self.glyph_count)
    }

    /// How far `glyph` moves the pen, in the face's own units.
    fn advance(&self, glyph: u16) -> u16 {
        let hmtx = self.table(*b"hmtx").map_or(0, |range| range.start);
        let at = glyph.min(self.metric_count - 1);
        u16_at(self.data, hmtx + 4 * usize::from(at))
    }

    /// How far `glyph` moves the pen, in thousandths of an em to the nearest
    /// one, which is the unit of a PDF's `/W` and of Adobe's metrics.
    pub(crate) fn width(&self, glyph: u16) -> u16 {
        let units = u32::from(self.advance(glyph));
        let thousandths = (units * 1000 + self.units_per_em / 2) / self.units_per_em;
        u16::try_from(thousandths).unwrap_or(u16::MAX)
    }

    /// The space before `glyph`'s outline, in the face's own units, as the
    /// bits it is written in.
    fn left_bearing(&self, glyph: u16) -> u16 {
        let hmtx = self.table(*b"hmtx").map_or(0, |range| range.start);
        let at = if glyph < self.metric_count {
            4 * usize::from(glyph) + 2
        } else {
            4 * usize::from(self.metric_count) + 2 * usize::from(glyph - self.metric_count)
        };
        u16_at(self.data, hmtx + at)
    }

    /// How many glyphs the face has.
    #[cfg(test)]
    pub(crate) fn glyph_count(&self) -> u16 {
        self.glyph_count
    }

    /// The bytes of `glyph`'s outline in `glyf`, empty for one that draws
    /// nothing, such as the space.
    pub(crate) fn outline(&self, glyph: u16) -> &'static [u8] {
        let (Some(loca), Some(glyf)) = (self.table(*b"loca"), self.table(*b"glyf")) else {
            return &[];
        };
        let entry = |at: usize| -> usize {
            if self.long_loca {
                u32_at(self.data, loca.start + 4 * at) as usize
            } else {
                2 * usize::from(u16_at(self.data, loca.start + 2 * at))
            }
        };
        let at = usize::from(glyph);
        let (start, end) = (glyf.start + entry(at), glyf.start + entry(at + 1));
        if start >= end || end > glyf.end {
            return &[];
        }
        &self.data[start..end]
    }

    /// The glyphs a composite glyph is built from; none for a simple one.
    fn components(&self, glyph: u16) -> Vec<u16> {
        let outline = self.outline(glyph);
        let mut parts = Vec::new();
        // A composite is the one with a negative count of contours.
        if outline.len() < 10 || i16::from_be_bytes([outline[0], outline[1]]) >= 0 {
            return parts;
        }
        let mut at = 10;
        loop {
            let flags = u16_at(outline, at);
            parts.push(u16_at(outline, at + 2));
            at += 4;
            // Its offset, as two bytes or two words.
            at += if flags & 0x0001 != 0 { 4 } else { 2 };
            // And its scale: one number, one for each axis, or a matrix.
            at += if flags & 0x0008 != 0 {
                2
            } else if flags & 0x0040 != 0 {
                4
            } else if flags & 0x0080 != 0 {
                8
            } else {
                0
            };
            if flags & 0x0020 == 0 || at >= outline.len() {
                break;
            }
        }
        parts
    }

    /// `used`, the `.notdef` every font starts with, and every glyph any of
    /// them is built from, however deep.
    pub(crate) fn closure(&self, used: impl IntoIterator<Item = u16>) -> BTreeSet<u16> {
        let mut kept = BTreeSet::new();
        let mut pending: Vec<u16> = std::iter::once(0).chain(used).collect();
        while let Some(glyph) = pending.pop() {
            if glyph < self.glyph_count && kept.insert(glyph) {
                pending.extend(self.components(glyph));
            }
        }
        kept
    }

    /// The face with only `used` drawn, and what those glyphs are built from:
    /// a TrueType file a PDF embeds as its `/FontFile2`.
    pub(crate) fn subset(&self, used: &BTreeSet<u16>) -> Vec<u8> {
        let kept = self.closure(used.iter().copied());
        let count = kept.last().map_or(1, |last| last + 1);

        let mut glyf = Vec::new();
        let mut offsets = Vec::with_capacity(usize::from(count) + 1);
        for glyph in 0..count {
            offsets.push(glyf.len());
            if kept.contains(&glyph) {
                glyf.extend_from_slice(self.outline(glyph));
                // Every glyph starts on an even byte, which is what lets the
                // short form of `loca` halve its offsets.
                if glyf.len() % 2 == 1 {
                    glyf.push(0);
                }
            }
        }
        offsets.push(glyf.len());
        let long = glyf.len() > 2 * usize::from(u16::MAX);
        let mut loca = Vec::with_capacity(offsets.len() * 4);
        for offset in &offsets {
            if long {
                loca.extend_from_slice(&(*offset as u32).to_be_bytes());
            } else {
                loca.extend_from_slice(&((*offset / 2) as u16).to_be_bytes());
            }
        }

        // Every glyph up to the highest kept has an advance and a bearing of
        // its own, so no reader has to work out which share the last one's.
        let mut hmtx = Vec::with_capacity(4 * usize::from(count));
        for glyph in 0..count {
            hmtx.extend_from_slice(&self.advance(glyph).to_be_bytes());
            hmtx.extend_from_slice(&self.left_bearing(glyph).to_be_bytes());
        }

        let copy = |tag: [u8; 4]| -> Vec<u8> {
            self.table(tag)
                .map_or_else(Vec::new, |range| self.data[range].to_vec())
        };
        let mut head = copy(*b"head");
        // The whole-file checksum is set once the file is assembled.
        head[8..12].fill(0);
        head[50..52].copy_from_slice(&u16::from(long).to_be_bytes());
        let mut hhea = copy(*b"hhea");
        hhea[34..36].copy_from_slice(&count.to_be_bytes());
        let mut maxp = copy(*b"maxp");
        maxp[4..6].copy_from_slice(&count.to_be_bytes());

        let mut tables: Vec<([u8; 4], Vec<u8>)> = vec![
            (*b"glyf", glyf),
            (*b"head", head),
            (*b"hhea", hhea),
            (*b"hmtx", hmtx),
            (*b"loca", loca),
            (*b"maxp", maxp),
        ];
        for tag in [*b"cvt ", *b"fpgm", *b"prep"] {
            if self.table(tag).is_some() {
                tables.push((tag, copy(tag)));
            }
        }
        sfnt(tables)
    }

    /// The numbers a PDF's font descriptor gives for the face, in thousandths
    /// of an em: what a reader sizes a selection by and draws a face it could
    /// not read with, which for an embedded face is never.
    fn descriptor(&self) -> String {
        let at = |tag: [u8; 4], offset: usize| -> i16 {
            self.table(tag)
                .map_or(0, |range| u16_at(self.data, range.start + offset) as i16)
        };
        let scale = |units: i16| {
            let thousandths = f64::from(units) * 1000.0 / f64::from(self.units_per_em);
            num(thousandths.round())
        };
        let (ascent, descent) = (at(*b"hhea", 4), at(*b"hhea", 6));
        // The height of a capital, which OS/2 has from its second version.
        let capitals = if at(*b"OS/2", 0) >= 2 {
            at(*b"OS/2", 88)
        } else {
            ascent
        };
        let fixed = self
            .table(*b"post")
            .is_some_and(|range| u32_at(self.data, range.start + 12) != 0);
        let italic = self.table(*b"post").map_or(0.0, |range| {
            f64::from(u32_at(self.data, range.start + 4) as i32) / 65536.0
        });
        // TrueType has no stem width, and a reader uses it only to fake a
        // face it does not have, so it is estimated from the weight class.
        let weight = f64::from(at(*b"OS/2", 4));
        let stem = (10.0 + 220.0 * (weight - 50.0) / 900.0).round().max(0.0);
        // Symbolic, as a font read through its own glyph numbers is, and of
        // fixed pitch where every glyph is the same width.
        let flags = 4 + u8::from(fixed);
        format!(
            "/Flags {flags} /FontBBox [{} {} {} {}] /ItalicAngle {} /Ascent {} /Descent {} \
             /CapHeight {} /StemV {}",
            scale(at(*b"head", 36)),
            scale(at(*b"head", 38)),
            scale(at(*b"head", 40)),
            scale(at(*b"head", 42)),
            num(italic),
            scale(ascent),
            scale(descent),
            scale(capitals),
            num(stem),
        )
    }
}

/// A TrueType file of `tables`, each on a four-byte boundary and listed in
/// the order of its tag, with the checksums every table and the file carry.
fn sfnt(mut tables: Vec<([u8; 4], Vec<u8>)>) -> Vec<u8> {
    tables.sort_by_key(|(tag, _)| *tag);
    let count = tables.len() as u16;
    let mut power = 1u16;
    let mut log = 0u16;
    while power * 2 <= count {
        power *= 2;
        log += 1;
    }
    let mut out = Vec::new();
    out.extend_from_slice(&0x0001_0000u32.to_be_bytes());
    out.extend_from_slice(&count.to_be_bytes());
    out.extend_from_slice(&(power * 16).to_be_bytes());
    out.extend_from_slice(&log.to_be_bytes());
    out.extend_from_slice(&(count * 16 - power * 16).to_be_bytes());
    let mut offset = 12 + 16 * tables.len();
    let mut head = None;
    for (tag, bytes) in &tables {
        out.extend_from_slice(tag);
        out.extend_from_slice(&checksum(bytes).to_be_bytes());
        out.extend_from_slice(&(offset as u32).to_be_bytes());
        out.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
        if tag == b"head" {
            head = Some(offset);
        }
        offset += bytes.len().next_multiple_of(4);
    }
    for (_, bytes) in &tables {
        out.extend_from_slice(bytes);
        out.resize(out.len().next_multiple_of(4), 0);
    }
    if let Some(head) = head {
        let adjustment = 0xB1B0_AFBAu32.wrapping_sub(checksum(&out));
        out[head + 8..head + 12].copy_from_slice(&adjustment.to_be_bytes());
    }
    out
}

/// The sum of `bytes` as big-endian 32-bit words, the last padded with
/// noughts, as TrueType checks a table.
fn checksum(bytes: &[u8]) -> u32 {
    bytes.chunks(4).fold(0u32, |sum, word| {
        let mut padded = [0u8; 4];
        padded[..word.len()].copy_from_slice(word);
        sum.wrapping_add(u32::from_be_bytes(padded))
    })
}

fn u16_at(data: &[u8], at: usize) -> u16 {
    data.get(at..at + 2)
        .map_or(0, |bytes| u16::from_be_bytes([bytes[0], bytes[1]]))
}

fn u32_at(data: &[u8], at: usize) -> u32 {
    data.get(at..at + 4).map_or(0, |bytes| {
        u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three faces in `fonts/`.
    fn faces() -> [Font; 3] {
        [
            include_bytes!("../../fonts/LiberationSans-Regular.ttf").as_slice(),
            include_bytes!("../../fonts/LiberationSans-Bold.ttf"),
            include_bytes!("../../fonts/LiberationMono-Regular.ttf"),
        ]
        .map(|data| Font::read(data).expect("TrueType"))
    }

    #[test]
    fn each_face_reads_with_its_name_and_its_letters() {
        let [sans, bold, mono] = faces();
        for (font, name, fixed) in [
            (&sans, "LiberationSans", false),
            (&bold, "LiberationSans-Bold", false),
            (&mono, "LiberationMono", true),
        ] {
            assert_eq!(font.name, name);
            assert_eq!(font.units_per_em, 2048);
            assert_eq!(font.glyph('A'), Some(36));
            assert_eq!(font.glyph('?'), Some(34));
            // Two glyphs a sample name may hold that no face here draws: a
            // Chinese character, and a double helix past the first plane.
            assert_eq!(font.glyph('\u{4e2d}'), None);
            assert_eq!(font.glyph('\u{1f9ec}'), None);
            let flags = if fixed { "/Flags 5 " } else { "/Flags 4 " };
            assert!(font.descriptor().starts_with(flags), "{name}");
        }
        // Liberation Sans: the box every glyph fits in, the line's ascent and
        // descent, a capital's height and, from the weight, a stem.
        assert_eq!(
            sans.descriptor(),
            "/Flags 4 /FontBBox [-544 -303 1302 980] /ItalicAngle 0 /Ascent 905 \
             /Descent -212 /CapHeight 688 /StemV 96"
        );
    }

    /// A subset is a TrueType file any reader takes: its tables listed in the
    /// order of their tags, each on a four-byte boundary and carrying its own
    /// checksum, and the whole file summing to the number TrueType fixes.
    #[test]
    fn a_subset_is_a_sound_truetype_file() {
        let [font, _, _] = faces();
        let used: BTreeSet<u16> = "rpoB caf\u{e9} \u{3c9}"
            .chars()
            .filter_map(|c| font.glyph(c))
            .collect();
        let file = font.subset(&used);
        let count = usize::from(u16_at(&file, 4));
        let mut tags = Vec::new();
        for at in 0..count {
            let record = 12 + 16 * at;
            let tag = String::from_utf8(file[record..record + 4].to_vec()).unwrap();
            let sum = u32_at(&file, record + 4);
            let offset = u32_at(&file, record + 8) as usize;
            let length = u32_at(&file, record + 12) as usize;
            assert_eq!(offset % 4, 0, "{tag}");
            let mut table = file[offset..offset + length].to_vec();
            if tag == "head" {
                // Its own sum is taken with the file's left out.
                table[8..12].fill(0);
            }
            assert_eq!(checksum(&table), sum, "{tag}");
            tags.push(tag);
        }
        assert_eq!(
            tags,
            ["cvt ", "fpgm", "glyf", "head", "hhea", "hmtx", "loca", "maxp", "prep"]
        );
        assert_eq!(checksum(&file), 0xB1B0_AFBA);

        // It reads back with every glyph it keeps where it was, at its own
        // advance, and the highest glyph kept is the last there is.
        let file: &'static [u8] = Box::leak(file.into_boxed_slice());
        let subset = Font::read(file).expect("TrueType");
        let kept = font.closure(used.iter().copied());
        assert_eq!(subset.glyph_count, kept.last().unwrap() + 1);
        assert_eq!(subset.metric_count, subset.glyph_count);
        for glyph in 0..subset.glyph_count {
            assert_eq!(subset.advance(glyph), font.advance(glyph), "glyph {glyph}");
            assert_eq!(subset.left_bearing(glyph), font.left_bearing(glyph));
        }
        for glyph in &kept {
            assert_eq!(
                subset.outline(*glyph),
                font.outline(*glyph),
                "glyph {glyph}"
            );
        }
        // An e with an acute is drawn from an e and an acute, and those come
        // with it although nothing asked for them.
        let acute = font.glyph('\u{e9}').unwrap();
        assert_eq!(font.components(acute), [font.glyph('e').unwrap(), 118]);
        assert!(kept.contains(&118));
    }
}
