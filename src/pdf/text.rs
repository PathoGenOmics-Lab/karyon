//! Text: which of the embedded faces sets a label, which glyphs draw its
//! characters, and how wide it comes out.
//!
//! A PDF carries its own faces, Liberation Sans, Liberation Sans Bold and
//! Liberation Mono, each as a subset of the glyphs the page uses, so every
//! reader draws the same letters and a journal's preflight finds every font
//! embedded. The SVG's own faces, Inter and JetBrains Mono, are what it gives
//! up. The layout survives it, because Liberation is drawn to the widths of
//! Helvetica and Arial, and Liberation Mono to Courier's:
//! [`text_width`](crate::svg::text_width) measures every printable ASCII
//! character, and every character beyond it that karyon writes, at least as
//! wide as Liberation Sans draws it,
//! [`text_width_strong`](crate::svg::text_width_strong) at least as wide as
//! Liberation Sans Bold, and [`mono_width`](crate::svg::mono_width) is
//! Liberation Mono's 600 exactly, so a label the SVG fitted in its room fits
//! in the same room here. Each label is anchored with the widths of the face
//! that draws it, the same ones the font's `/W` gives, so a label anchored at
//! its end ends exactly where the SVG ended it. A character karyon does not
//! write itself, such as an accented letter in a sample name, is measured at
//! a flat 600 thousandths, and Liberation Sans draws an accented capital
//! wider than that, so such a label can still run over at the end it was not
//! anchored by.

use std::sync::OnceLock;

use super::font::Font;

/// One of the faces a PDF embeds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Face {
    Sans,
    SansBold,
    Mono,
}

impl Face {
    /// The face, read once from the copy the crate carries in `fonts/`.
    pub(crate) fn font(self) -> &'static Font {
        static SANS: OnceLock<Font> = OnceLock::new();
        static SANS_BOLD: OnceLock<Font> = OnceLock::new();
        static MONO: OnceLock<Font> = OnceLock::new();
        let (cell, data): (&OnceLock<Font>, &'static [u8]) = match self {
            Face::Sans => (
                &SANS,
                include_bytes!("../../fonts/LiberationSans-Regular.ttf"),
            ),
            Face::SansBold => (
                &SANS_BOLD,
                include_bytes!("../../fonts/LiberationSans-Bold.ttf"),
            ),
            Face::Mono => (
                &MONO,
                include_bytes!("../../fonts/LiberationMono-Regular.ttf"),
            ),
        };
        cell.get_or_init(|| Font::read(data).expect("the faces in fonts/ are TrueType"))
    }

    /// The face a font stack and a weight ask for, and what the page says
    /// when it is not the face they asked for.
    ///
    /// The stack decides by kind rather than by name, since a PDF carries
    /// only the faces it embeds: a stack that names a monospaced face
    /// anywhere, or ends in `monospace`, is Liberation Mono, and anything
    /// else Liberation Sans, which is what the default stack, Inter then
    /// Liberation Sans, Arial and Helvetica, is drawn to match. Bold is
    /// `bold`, `bolder`, or a weight of 600 or more, which is where CSS
    /// rounds a semibold to the bold face of a family that has only the two.
    /// A serif stack, one that opens with a serif face or ends in `serif`,
    /// and a bold monospaced one ask for a face the PDF does not carry, and
    /// karyon writes neither, so each is set in the nearest it does carry
    /// and named in the notes.
    pub(crate) fn choose(family: &str, weight: &str) -> (Face, Option<&'static str>) {
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
        let sans = if bold { Face::SansBold } else { Face::Sans };
        match (mono, serif, bold) {
            (true, _, false) => (Face::Mono, None),
            (true, _, true) => (
                Face::Mono,
                Some(
                    "bold monospaced text is set in Liberation Mono, whose one weight the PDF \
                     embeds is regular",
                ),
            ),
            (false, true, _) => (
                sans,
                Some("a serif face is set in Liberation Sans, since the PDF embeds no serif face"),
            ),
            (false, false, _) => (sans, None),
        }
    }
}

/// A stretch of a label set in one face.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Run {
    pub(crate) face: Face,
    /// The glyphs that draw it, each with the character it stands for.
    pub(crate) glyphs: Vec<(u16, char)>,
    /// Its advance, in thousandths of an em.
    pub(crate) width: f64,
}

/// `text` cut into runs by the face that can set each character.
///
/// A character goes to `base` when it has it, and otherwise to Liberation
/// Sans, which has every glyph the bold and the monospaced face have and a
/// few more. One neither has, or one written right to left, is drawn as a
/// question mark in `base` and added to `missing`, once. A question mark is a
/// visible stand-in rather than a gap, so the label keeps its length and a
/// reader can see something was there; the document's title and alt text
/// keep the exact text either way.
pub(crate) fn runs(text: &str, base: Face, missing: &mut Vec<char>) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for c in text.chars().filter(|c| !invisible(*c)) {
        let found = if right_to_left(c) {
            None
        } else {
            base.font()
                .glyph(c)
                .map(|glyph| (base, glyph))
                .or_else(|| Face::Sans.font().glyph(c).map(|glyph| (Face::Sans, glyph)))
        };
        let (face, glyph, c) = match found {
            Some((face, glyph)) => (face, glyph, c),
            None => {
                if !missing.contains(&c) {
                    missing.push(c);
                }
                (base, base.font().glyph('?').unwrap_or(0), '?')
            }
        };
        let width = f64::from(face.font().width(glyph));
        match runs.last_mut() {
            Some(run) if run.face == face => {
                run.glyphs.push((glyph, c));
                run.width += width;
            }
            _ => runs.push(Run {
                face,
                glyphs: vec![(glyph, c)],
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

/// Letters written right to left: Hebrew, Arabic and the scripts beside
/// them. A browser lays a run of them out from the right, and a PDF draws
/// glyphs in the order they are given, so Hebrew, which Liberation has, would
/// come out backwards; it is marked as missing instead, as Arabic, which
/// Liberation does not have, would be anyway.
pub(crate) fn right_to_left(c: char) -> bool {
    matches!(
        c,
        '\u{590}'..='\u{8ff}' | '\u{fb1d}'..='\u{fdff}' | '\u{fe70}'..='\u{fefc}'
    )
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

/// Writes `glyphs` as the string a composite font reads them from: two bytes
/// each, the glyph's number, in hexadecimal.
pub(crate) fn write_glyphs(out: &mut Vec<u8>, glyphs: &[(u16, char)]) {
    out.push(b'<');
    for (glyph, _) in glyphs {
        out.extend_from_slice(format!("{glyph:04X}").as_bytes());
    }
    out.push(b'>');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_stacks_are_set_in_liberation_sans_and_mono() {
        let sans = "Inter, Liberation Sans, Arial, Helvetica, sans-serif";
        let mono = "JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace";
        assert_eq!(Face::choose(sans, "normal"), (Face::Sans, None));
        assert_eq!(Face::choose(sans, "bold"), (Face::SansBold, None));
        assert_eq!(Face::choose(sans, "600"), (Face::SansBold, None));
        assert_eq!(Face::choose(sans, "500"), (Face::Sans, None));
        assert_eq!(Face::choose(mono, "normal"), (Face::Mono, None));
        assert_eq!(
            Face::choose("Fira Sans, sans-serif", "normal"),
            (Face::Sans, None)
        );
        // The faces the PDF does not carry are set in the nearest it does,
        // and said.
        let (face, note) = Face::choose("Fira Code, monospace", "bold");
        assert_eq!(face, Face::Mono);
        assert!(note.is_some_and(|note| note.contains("Liberation Mono")));
        for (stack, weight, face) in [
            ("Georgia, serif", "normal", Face::Sans),
            ("'Source Serif 4', serif", "700", Face::SansBold),
        ] {
            let (chosen, note) = Face::choose(stack, weight);
            assert_eq!(chosen, face, "{stack}");
            assert!(note.is_some_and(|note| note.contains("serif")), "{stack}");
        }
    }

    #[test]
    fn omega_the_relations_and_the_ellipsis_are_set_in_the_label_s_own_face() {
        let mut missing = Vec::new();
        let set = runs("\u{3c9} \u{2264} 0.05\u{2026}", Face::Sans, &mut missing);
        assert_eq!(set.len(), 1);
        assert_eq!(set[0].face, Face::Sans);
        assert!(missing.is_empty());
        let text: String = set[0].glyphs.iter().map(|(_, c)| *c).collect();
        assert_eq!(text, "\u{3c9} \u{2264} 0.05\u{2026}");
        // Omega is 781 thousandths in Liberation Sans, a space 278, less
        // than or equal 549, the digits of 0.05 556 each and its full stop
        // 278, and the ellipsis an em.
        assert_eq!(
            set[0].width,
            781.0 + 278.0 + 549.0 + 278.0 + 1946.0 + 1000.0
        );
        let mono = runs("\u{3c9}\u{2192}", Face::Mono, &mut missing);
        assert_eq!(mono.len(), 1);
        assert_eq!(mono[0].width, 1200.0);
    }

    #[test]
    fn a_character_no_embedded_face_has_is_a_question_mark_and_said_once() {
        let mut missing = Vec::new();
        let set = runs("a\u{4e2d}b\u{4e2d}", Face::Sans, &mut missing);
        assert_eq!(set.len(), 1);
        let text: String = set[0].glyphs.iter().map(|(_, c)| *c).collect();
        assert_eq!(text, "a?b?");
        assert_eq!(missing, ['\u{4e2d}']);
        // An element-of sign was set in Symbol before the faces were
        // embedded, and Liberation has none.
        let mut missing = Vec::new();
        runs("x \u{2208} S", Face::SansBold, &mut missing);
        assert_eq!(missing, ['\u{2208}']);
        // Hebrew is in Liberation, and would be drawn backwards.
        let mut missing = Vec::new();
        let set = runs("\u{5e9}\u{5dc}\u{5d5}\u{5dd} 1", Face::Sans, &mut missing);
        let text: String = set[0].glyphs.iter().map(|(_, c)| *c).collect();
        assert_eq!(text, "???? 1");
        assert_eq!(missing.len(), 4);
        // What a browser draws as nothing is left out rather than marked.
        let mut missing = Vec::new();
        let set = runs("rpoB\u{200f}\u{ad}", Face::Sans, &mut missing);
        let text: String = set[0].glyphs.iter().map(|(_, c)| *c).collect();
        assert_eq!(text, "rpoB");
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
    fn glyphs_are_written_two_bytes_each_in_hexadecimal() {
        let mut out = Vec::new();
        write_glyphs(&mut out, &[(36, 'A'), (0x0b1c, '\u{3c9}'), (3, ' ')]);
        assert_eq!(out, b"<00240B1C0003>");
    }
}
