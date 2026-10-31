//! Advance widths of Helvetica and Helvetica-Bold, from Adobe's metrics, which
//! the faces a PDF embeds are held to.
//!
//! The layout measures text with [`text_width`](crate::svg::text_width) and
//! [`text_width_strong`](crate::svg::text_width_strong), whose tables were
//! made never to be narrower than these. Liberation Sans is drawn to the same
//! widths, and Liberation Mono to Courier's 600 thousandths of an em for every
//! glyph, which is what lets a PDF set its text in them without moving
//! anything; a test reads each face's own advances and compares them with
//! these, character by character. Nothing outside the tests reads them, since
//! the PDF anchors each label with the widths of the face that draws it.
//!
//! The numbers are the `WX` entries of the Adobe Font Metrics files for the
//! standard fourteen fonts, which carry this notice:
//!
//! > Copyright (c) 1985, 1987, 1989, 1990, 1997 Adobe Systems Incorporated.
//! > All Rights Reserved.
//!
//! They were read from the copy of those files that matplotlib ships, in
//! `mpl-data/fonts/pdfcorefonts`, and indexed by byte through the glyph names
//! of the Adobe Glyph List: each table runs from the space, byte 32, to byte
//! 255 of WinAnsiEncoding, the encoding a PDF sets those faces in. The six
//! bytes WinAnsi leaves undefined are nought.
//!
//! Five characters karyon writes are not in WinAnsi, and a PDF set them in
//! Symbol before it embedded its faces, at the widths in [`SYMBOL`].

/// Helvetica, by WinAnsi byte from 32.
pub(crate) const HELVETICA: [u16; 224] = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278,
    278, // 0x20 to 0x2f
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584,
    556, // 0x30 to 0x3f
    1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722,
    778, // 0x40 to 0x4f
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469,
    556, // 0x50 to 0x5f
    333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556,
    556, // 0x60 to 0x6f
    556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,
    0, // 0x70 to 0x7f
    556, 0, 222, 556, 333, 1000, 556, 556, 333, 1000, 667, 333, 1000, 0, 611,
    0, // 0x80 to 0x8f
    0, 222, 222, 333, 333, 350, 556, 1000, 333, 1000, 500, 333, 944, 0, 500,
    667, // 0x90 to 0x9f
    278, 333, 556, 556, 556, 556, 260, 556, 333, 737, 370, 556, 584, 333, 737,
    333, // 0xa0 to 0xaf
    400, 584, 333, 333, 333, 556, 537, 278, 333, 333, 365, 556, 834, 834, 834,
    611, // 0xb0 to 0xbf
    667, 667, 667, 667, 667, 667, 1000, 722, 667, 667, 667, 667, 278, 278, 278,
    278, // 0xc0 to 0xcf
    722, 722, 778, 778, 778, 778, 778, 584, 778, 722, 722, 722, 722, 667, 667,
    611, // 0xd0 to 0xdf
    556, 556, 556, 556, 556, 556, 889, 500, 556, 556, 556, 556, 278, 278, 278,
    278, // 0xe0 to 0xef
    556, 556, 556, 556, 556, 556, 556, 584, 611, 556, 556, 556, 556, 500, 556,
    500, // 0xf0 to 0xff
];

/// Helvetica-Bold, by WinAnsi byte from 32.
pub(crate) const HELVETICA_BOLD: [u16; 224] = [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278,
    278, // 0x20 to 0x2f
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584,
    611, // 0x30 to 0x3f
    975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722,
    778, // 0x40 to 0x4f
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584,
    556, // 0x50 to 0x5f
    333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611,
    611, // 0x60 to 0x6f
    611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584,
    0, // 0x70 to 0x7f
    556, 0, 278, 556, 500, 1000, 556, 556, 333, 1000, 667, 333, 1000, 0, 611,
    0, // 0x80 to 0x8f
    0, 278, 278, 500, 500, 350, 556, 1000, 333, 1000, 556, 333, 944, 0, 500,
    667, // 0x90 to 0x9f
    278, 333, 556, 556, 556, 556, 280, 556, 333, 737, 370, 556, 584, 333, 737,
    333, // 0xa0 to 0xaf
    400, 584, 333, 333, 333, 611, 556, 278, 333, 333, 365, 556, 834, 834, 834,
    611, // 0xb0 to 0xbf
    722, 722, 722, 722, 722, 722, 1000, 722, 667, 667, 667, 667, 278, 278, 278,
    278, // 0xc0 to 0xcf
    722, 722, 778, 778, 778, 778, 778, 584, 778, 722, 722, 722, 722, 667, 667,
    611, // 0xd0 to 0xdf
    556, 556, 556, 556, 556, 556, 889, 556, 556, 556, 556, 556, 278, 278, 278,
    278, // 0xe0 to 0xef
    611, 611, 611, 611, 611, 611, 611, 584, 611, 611, 611, 611, 611, 556, 611,
    556, // 0xf0 to 0xff
];

/// Symbol's widths for the five characters karyon writes beyond WinAnsi.
pub(crate) const SYMBOL: [(char, u16); 5] = [
    ('\u{03c9}', 686), // omega
    ('\u{2192}', 987), // arrowright
    ('\u{2248}', 549), // approxequal
    ('\u{2264}', 549), // lessequal
    ('\u{2265}', 549), // greaterequal
];

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
