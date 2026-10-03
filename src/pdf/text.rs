//! Text: which of the base fonts sets a label, which bytes stand for its
//! characters, and how wide it comes out.
//!
//! No font is embedded. A face is 300 to 400 kilobytes here, and embedding one
//! means reading TrueType, cutting a subset and writing a composite font with a
//! map back to Unicode, which is most of a PDF library for a figure that is
//! mostly rectangles. Every PDF reader carries Helvetica, Times, Courier and
//! Symbol, or a face drawn to their exact widths, so a label is set in one of
//! those and stays text: it can be searched, copied and read aloud.
//!
//! What that costs is the SVG's own faces, Inter and JetBrains Mono, and
//! Helvetica is what a reader sees instead. The layout survives it.
//! [`text_width`](crate::svg::text_width) measures every printable ASCII
//! character, and every character beyond it that karyon writes, at least as
//! wide as Helvetica or Symbol draws it,
//! [`text_width_strong`](crate::svg::text_width_strong) at least as wide as
//! Helvetica-Bold, and [`mono_width`](crate::svg::mono_width) is Courier's
//! 600 exactly, so a label the SVG fitted in its room fits in the same room
//! here. Each label is anchored with Adobe's own widths whatever it was
//! measured with, so a label anchored at its end ends exactly where the SVG
//! ended it. A character karyon does not write itself, such as an accented
//! letter in a sample name, is measured at a flat 600 thousandths, and
//! Helvetica draws an accented capital wider than that, so such a label can
//! still run over at the end it was not anchored by.

use super::metrics;

/// One of the fonts every PDF reader has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Face {
    Helvetica,
    HelveticaBold,
    Times,
    TimesBold,
    Courier,
    CourierBold,
    Symbol,
}

impl Face {
    /// The name a PDF gives the font.
    pub(crate) fn base_font(self) -> &'static str {
        match self {
            Face::Helvetica => "Helvetica",
            Face::HelveticaBold => "Helvetica-Bold",
            Face::Times => "Times-Roman",
            Face::TimesBold => "Times-Bold",
            Face::Courier => "Courier",
            Face::CourierBold => "Courier-Bold",
            Face::Symbol => "Symbol",
        }
    }

    /// The advance of WinAnsi byte `code`, in thousandths of an em.
    fn width(self, code: u8) -> u16 {
        let index = usize::from(code.saturating_sub(32));
        match self {
            Face::Helvetica => metrics::HELVETICA[index],
            Face::HelveticaBold => metrics::HELVETICA_BOLD[index],
            Face::Times => metrics::TIMES_ROMAN[index],
            Face::TimesBold => metrics::TIMES_BOLD[index],
            Face::Courier | Face::CourierBold => 600,
            Face::Symbol => symbol_width(code),
        }
    }

    /// The face a font stack and a weight ask for.
    ///
    /// The stack decides by kind rather than by name, since no name in it is
    /// one a PDF reader has: a stack that names a monospaced face anywhere, or
    /// ends in `monospace`, is Courier; one that opens with a serif face, or
    /// ends in `serif`, is Times; anything else is Helvetica, which is what the
    /// default stack, Inter then Liberation Sans, Arial and Helvetica, is
    /// drawn to match. Bold is `bold`, `bolder`, or a weight of 600 or more,
    /// which is where CSS rounds a semibold to the bold face of a family that
    /// has only the two.
    pub(crate) fn choose(family: &str, weight: &str) -> Face {
        let weight = weight.trim();
        let bold = matches!(weight, "bold" | "bolder")
            || weight.parse::<f64>().is_ok_and(|value| value >= 600.0);
        let families: Vec<String> = family
            .split(',')
            .map(|name| {
                name.trim()
                    .trim_matches(|c: char| c == '"' || c == '\'')
                    .to_ascii_lowercase()
            })
            .collect();
        let mono = families.iter().any(|name| {
            name == "monospace"
                || ["mono", "courier", "consolas", "menlo", "monaco"]
                    .iter()
                    .any(|kind| name.contains(kind))
        });
        let serif = |name: &str| {
            name == "serif"
                || (name.contains("serif") && !name.contains("sans"))
                || ["times", "georgia", "garamond", "cambria", "palatino"]
                    .iter()
                    .any(|kind| name.contains(kind))
        };
        let serif = families.first().is_some_and(|name| serif(name))
            || families.iter().any(|name| name == "serif");
        match (mono, serif, bold) {
            (true, _, false) => Face::Courier,
            (true, _, true) => Face::CourierBold,
            (false, true, false) => Face::Times,
            (false, true, true) => Face::TimesBold,
            (false, false, false) => Face::Helvetica,
            (false, false, true) => Face::HelveticaBold,
        }
    }
}

/// A stretch of a label set in one face.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Run {
    pub(crate) face: Face,
    /// The bytes that draw it, in the face's encoding.
    pub(crate) bytes: Vec<u8>,
    /// Its advance, in thousandths of an em.
    pub(crate) width: f64,
}

/// `text` cut into runs by the face that can set each character.
///
/// A character goes to `base` when WinAnsi has it, to Symbol when Symbol has
/// it, and is otherwise drawn as a question mark in `base` and added to
/// `missing`, once. A question mark is a visible stand-in rather than a gap,
/// so the label keeps its length and a reader can see something was there;
/// the document's title and alt text keep the exact text either way.
pub(crate) fn runs(text: &str, base: Face, missing: &mut Vec<char>) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for c in text.chars().filter(|c| !invisible(*c)) {
        let (face, code) = if let Some(code) = winansi(c) {
            (base, code)
        } else if let Some(code) = symbol(c) {
            (Face::Symbol, code)
        } else {
            if !missing.contains(&c) {
                missing.push(c);
            }
            (base, b'?')
        };
        let width = f64::from(face.width(code));
        match runs.last_mut() {
            Some(run) if run.face == face => {
                run.bytes.push(code);
                run.width += width;
            }
            _ => runs.push(Run {
                face,
                bytes: vec![code],
                width,
            }),
        }
    }
    runs
}

/// Characters a browser draws as nothing: controls, the soft hyphen, the
/// zero-width spaces and joiners, the marks that steer bidirectional text, and
/// variation selectors. Drawn as a question mark they would add a letter the
/// SVG never showed.
fn invisible(c: char) -> bool {
    matches!(
        c,
        '\u{0}'..='\u{1f}'
            | '\u{7f}'..='\u{9f}'
            | '\u{ad}'
            | '\u{200b}'..='\u{200f}'
            | '\u{2028}'..='\u{202e}'
            | '\u{2060}'..='\u{206f}'
            | '\u{fe00}'..='\u{fe0f}'
            | '\u{feff}'
    )
}

/// The WinAnsi byte for `c`, if WinAnsi has one.
pub(crate) fn winansi(c: char) -> Option<u8> {
    let code = c as u32;
    if (0x20..0x7f).contains(&code) || (0xa0..=0xff).contains(&code) {
        // ASCII and Latin-1 are where Unicode put them.
        return u8::try_from(code).ok();
    }
    // The rest of Windows-1252, which puts typographic punctuation and a few
    // letters where Latin-1 has its second block of controls.
    const HIGH: [(char, u8); 27] = [
        ('\u{20ac}', 0x80),
        ('\u{201a}', 0x82),
        ('\u{0192}', 0x83),
        ('\u{201e}', 0x84),
        ('\u{2026}', 0x85),
        ('\u{2020}', 0x86),
        ('\u{2021}', 0x87),
        ('\u{02c6}', 0x88),
        ('\u{2030}', 0x89),
        ('\u{0160}', 0x8a),
        ('\u{2039}', 0x8b),
        ('\u{0152}', 0x8c),
        ('\u{017d}', 0x8e),
        ('\u{2018}', 0x91),
        ('\u{2019}', 0x92),
        ('\u{201c}', 0x93),
        ('\u{201d}', 0x94),
        ('\u{2022}', 0x95),
        ('\u{2013}', 0x96),
        ('\u{2014}', 0x97),
        ('\u{02dc}', 0x98),
        ('\u{2122}', 0x99),
        ('\u{0161}', 0x9a),
        ('\u{203a}', 0x9b),
        ('\u{0153}', 0x9c),
        ('\u{017e}', 0x9e),
        ('\u{0178}', 0x9f),
    ];
    HIGH.iter()
        .find(|(known, _)| *known == c)
        .map(|(_, code)| *code)
}

/// The Symbol byte for `c`, if Symbol has it and WinAnsi does not.
pub(crate) fn symbol(c: char) -> Option<u8> {
    metrics::SYMBOL
        .binary_search_by(|(known, _, _)| known.cmp(&c))
        .ok()
        .map(|at| metrics::SYMBOL[at].1)
}

fn symbol_width(code: u8) -> u16 {
    metrics::SYMBOL
        .iter()
        .find(|(_, known, _)| *known == code)
        .map_or(0, |(_, _, width)| *width)
}

/// Whitespace collapsed as a browser collapses it under `white-space:
/// normal`: every run of spaces, tabs and line breaks is one space, and none
/// is kept at either end.
pub(crate) fn collapse(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for word in text.split([' ', '\t', '\n', '\r', '\u{c}']) {
        if word.is_empty() {
            continue;
        }
        if !out.is_empty() {
            out.push(' ');
        }
        out.push_str(word);
    }
    out
}

/// Writes `bytes` as a PDF literal string that stays seven-bit: the three
/// characters a string has to escape are escaped, and every byte from 128 up
/// is an octal escape, so a content stream can be read and diffed as text.
pub(crate) fn write_string(out: &mut Vec<u8>, bytes: &[u8]) {
    out.push(b'(');
    for &byte in bytes {
        match byte {
            b'(' | b')' | b'\\' => {
                out.push(b'\\');
                out.push(byte);
            }
            0x20..=0x7e => out.push(byte),
            _ => out.extend_from_slice(format!("\\{byte:03o}").as_bytes()),
        }
    }
    out.push(b')');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_stacks_are_set_in_helvetica_and_courier() {
        let sans = "Inter, Liberation Sans, Arial, Helvetica, sans-serif";
        let mono = "JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace";
        assert_eq!(Face::choose(sans, "normal"), Face::Helvetica);
        assert_eq!(Face::choose(sans, "bold"), Face::HelveticaBold);
        assert_eq!(Face::choose(sans, "600"), Face::HelveticaBold);
        assert_eq!(Face::choose(sans, "500"), Face::Helvetica);
        assert_eq!(Face::choose(mono, "normal"), Face::Courier);
        assert_eq!(
            Face::choose("Fira Code, monospace", "bold"),
            Face::CourierBold
        );
        assert_eq!(Face::choose("Georgia, serif", "normal"), Face::Times);
        assert_eq!(
            Face::choose("'Source Serif 4', serif", "700"),
            Face::TimesBold
        );
        assert_eq!(
            Face::choose("Fira Sans, sans-serif", "normal"),
            Face::Helvetica
        );
    }

    #[test]
    fn omega_and_less_or_equal_are_set_in_symbol_and_the_rest_in_winansi() {
        let mut missing = Vec::new();
        let set = runs("\u{3c9} \u{2264} 0.05", Face::Helvetica, &mut missing);
        let faces: Vec<Face> = set.iter().map(|run| run.face).collect();
        assert_eq!(
            faces,
            [Face::Symbol, Face::Helvetica, Face::Symbol, Face::Helvetica]
        );
        assert_eq!(set[0].bytes, [0x77]);
        assert_eq!(set[2].bytes, [0xa3]);
        assert!(missing.is_empty());
        // Omega is 686 thousandths in Symbol, and a space 278 in Helvetica.
        assert_eq!(set[0].width, 686.0);
        assert_eq!(set[1].width, 278.0);
        // An ellipsis, a times sign, an en dash, a superscript two and a
        // middle dot are WinAnsi's.
        for (c, code) in [
            ('\u{2026}', 0x85),
            ('\u{d7}', 0xd7),
            ('\u{2013}', 0x96),
            ('\u{b2}', 0xb2),
            ('\u{b7}', 0xb7),
        ] {
            assert_eq!(winansi(c), Some(code), "{c}");
        }
        for (c, code) in [('\u{2265}', 0xb3), ('\u{2248}', 0xbb), ('\u{2192}', 0xae)] {
            assert_eq!(winansi(c), None, "{c}");
            assert_eq!(symbol(c), Some(code), "{c}");
        }
    }

    #[test]
    fn a_character_no_builtin_font_has_is_a_question_mark_and_said_once() {
        let mut missing = Vec::new();
        let set = runs("a\u{4e2d}b\u{4e2d}", Face::Helvetica, &mut missing);
        assert_eq!(set.len(), 1);
        assert_eq!(set[0].bytes, b"a?b?");
        assert_eq!(missing, ['\u{4e2d}']);
        // What a browser draws as nothing is left out rather than marked.
        let mut missing = Vec::new();
        let set = runs("rpoB\u{200f}\u{ad}", Face::Helvetica, &mut missing);
        assert_eq!(set[0].bytes, b"rpoB");
        assert!(missing.is_empty());
    }

    #[test]
    fn whitespace_in_a_label_collapses_as_a_browser_collapses_it() {
        assert_eq!(collapse("  gene\t A\n\nB  "), "gene A B");
        assert_eq!(collapse(" \n "), "");
        // A no-break space is not whitespace to collapse.
        assert_eq!(collapse("1\u{a0} kb"), "1\u{a0} kb");
    }

    #[test]
    fn a_string_stays_seven_bit_with_its_delimiters_escaped() {
        let mut out = Vec::new();
        write_string(&mut out, b"f(x) \\ \x85\xd7");
        assert_eq!(out, b"(f\\(x\\) \\\\ \\205\\327)");
    }
}
