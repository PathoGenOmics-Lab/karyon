//! Tests of the conversion as a whole: what the writer writes reaches the
//! page, and where on it.

use super::*;
use crate::style::{LinePattern, Symbol};
use crate::svg::{Anchor, SvgWriter, TextStyle};

fn convert(svg: &str) -> Pdf {
    Pdf::from_svg(svg).expect("a sized root")
}

/// The content stream, which is object 4 and the first in the file.
fn content(pdf: &Pdf) -> String {
    let text = String::from_utf8_lossy(&pdf.bytes);
    let start = text.find("4 0 obj\n").expect("object 4");
    let begin = start + text[start..].find("stream\n").expect("a stream") + "stream\n".len();
    let end = begin
        + text[begin..]
            .find("endstream")
            .expect("the end of the stream");
    text[begin..end].to_string()
}

fn text(pdf: &Pdf) -> String {
    String::from_utf8_lossy(&pdf.bytes).into_owned()
}

/// The bytes of the first stream after `marker`, found by byte rather than in
/// [`text`], whose replacement characters move every offset after the
/// binary comment on the second line.
fn stream_after<'a>(pdf: &'a Pdf, marker: &str) -> &'a [u8] {
    let find = |from: usize, needle: &[u8]| {
        from + pdf.bytes[from..]
            .windows(needle.len())
            .position(|window| window == needle)
            .unwrap_or_else(|| panic!("no {}", String::from_utf8_lossy(needle)))
    };
    let begin = find(find(0, marker.as_bytes()), b"stream\n") + b"stream\n".len();
    let end = find(begin, b"\nendstream");
    &pdf.bytes[begin..end]
}

/// Every public method of the writer, called at least once, with every
/// pattern, symbol and anchor it takes, and a sheet of two panels around it.
///
/// This is the lock that holds the writer and the converter together. The
/// writer is the only way out of the crate, so anything it can write has to
/// be something the converter reads; a new element or attribute that is not
/// comes back as a note, and this test prints the note. A new method is
/// caught by name before that, by
/// `every_method_of_the_svg_writer_is_called_by_the_lock`, which reads the
/// writer's source for its methods and this function's for the calls.
#[test]
fn every_element_the_svg_writer_writes_is_one_the_pdf_reads() {
    let mut svg = SvgWriter::with_id_prefix("lock-");
    svg.describe("every element", "one of each thing the writer writes");
    let fade = svg.fade_down("#0072b2", 0.6, 0.05);
    svg.path("M0 40 L20 10 L40 30 L40 40 Z", &fade, 1.0);
    svg.rect(1.0, 2.0, 3.0, 4.0, "#d55e00");
    svg.rect_opacity(1.0, 2.0, 3.0, 4.0, "#d55e00", 0.4);
    svg.rect_rounded(1.0, 2.0, 30.0, 4.0, 2.0, "#009e73");
    svg.rect_rounded_opacity(1.0, 2.0, 30.0, 4.0, 2.0, "#009e73", 0.5);
    svg.rect_rounded_edged(1.0, 2.0, 30.0, 8.0, 3.0, "#ffffff", "#000000", 1.0);
    svg.rect_outline(1.0, 2.0, 30.0, 8.0, "#000000", 0.5);
    svg.circle(10.0, 10.0, 3.0, "#cc79a7");
    svg.circle_ringed(10.0, 10.0, 3.0, "#cc79a7", "#ffffff", 1.0);
    for pattern in [LinePattern::Solid, LinePattern::Dashed, LinePattern::Dotted] {
        svg.line(0.0, 1.0, 20.0, 1.0, "#333333", 1.0);
        svg.line_pattern(0.0, 1.0, 20.0, 1.0, "#333333", 1.0, pattern);
        svg.polyline(&[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)], "#333333", 1.0);
        svg.polyline_pattern(&[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)], "#333", 1.0, pattern);
        svg.path_stroked("M0 0 Q5 10 10 0 C12 2 14 2 16 0", "#333333", 1.0);
        svg.path_stroked_pattern("M0 0 A5 5 0 0 1 10 0", "#333333", 1.0, pattern);
    }
    for symbol in [
        Symbol::Circle,
        Symbol::Square,
        Symbol::Diamond,
        Symbol::Triangle,
    ] {
        svg.symbol(10.0, 10.0, 3.0, symbol, "#56b4e9");
        svg.symbol_ringed(10.0, 10.0, 3.0, symbol, "#56b4e9", "#ffffff", 1.0);
    }
    svg.polygon(&[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)], "#e69f00");
    svg.polygon_edged(
        &[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)],
        "#e69f00",
        "#000",
        0.5,
    );
    svg.begin_clip(0.0, 0.0, 50.0, 50.0);
    svg.begin_clip_path("M0 0 L50 0 L25 40 Z");
    svg.begin_titled("a tooltip");
    svg.begin_titled("");
    svg.begin_titled_inert("an inert tooltip");
    for anchor in [Anchor::Start, Anchor::Middle, Anchor::End] {
        svg.text(
            10.0,
            10.0,
            "depth \u{2264} 5 \u{3c9} \u{2026}",
            "#111111",
            9.0,
            anchor,
        );
        svg.text_bold(10.0, 10.0, "title", "#111111", 12.0, anchor);
        svg.text_rotated((10.0, 10.0), -90.0, "column", "#111111", 9.0, anchor);
        for bold in [false, true] {
            svg.text_haloed(10.0, 10.0, "rpoB", "#111111", "#ffffff", 9.0, anchor, bold);
        }
        // A family and a weight, each as karyon writes it: the locus in the
        // monospaced stack, and a track's name set semibold. Both on one
        // label would ask for a bold monospaced face, which the PDF does not
        // carry and says so, as a test of its own holds.
        svg.text_styled(
            10.0,
            10.0,
            "761,000 r\u{b2}",
            "#111111",
            9.0,
            anchor,
            TextStyle {
                family: Some("JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace"),
                tracking: 0.03,
                ..TextStyle::default()
            },
        );
        svg.text_styled(
            10.0,
            10.0,
            "depth",
            "#111111",
            9.0,
            anchor,
            TextStyle {
                weight: Some(600),
                tracking: 0.03,
                ..TextStyle::default()
            },
        );
        svg.text_styled(
            10.0,
            10.0,
            "plain",
            "#111",
            9.0,
            anchor,
            TextStyle::default(),
        );
    }
    svg.glyph(10.0, 20.0, 6.0, 14.0, "A", "#109648");
    // Writes nothing, and is called so the lock below finds it.
    svg.keep_clear(10.0, 10.0, 20.0, 8.0);
    assert_eq!(svg.kept_clear(), [(10.0, 10.0, 20.0, 8.0)]);
    svg.end_group();
    svg.end_group();
    let document = svg.finish(200.0, 120.0, "#ffffff", "Inter, Arial, sans-serif");

    // A sheet: the panels nest whole documents, each moved by a group that
    // carries its name.
    let figure = || {
        crate::Figure::new(crate::Region::parse("chr1:1-1000").unwrap())
            .title("a panel")
            .push(crate::CoverageTrack::new(0, vec![30.0; 1000]).label("depth"))
    };
    let sheet = crate::Panels::new()
        .push(&figure(), "A")
        .push_captioned(&figure(), "B", "the same again")
        .to_svg();

    for (what, svg) in [("every method", document), ("a sheet", sheet)] {
        let pdf = convert(&svg);
        assert!(pdf.notes.is_empty(), "{what}: {:#?}", pdf.notes);
        assert!(content(&pdf).contains(" Tj\n"), "{what} drew no text");
    }
}

/// The source of the lock test above, as far as its closing brace.
fn lock_source() -> &'static str {
    const SOURCE: &str = include_str!("tests.rs");
    let start = SOURCE
        .find("fn every_element_the_svg_writer_writes_is_one_the_pdf_reads")
        .expect("the lock test");
    let end = start + SOURCE[start..].find("\n}\n").expect("its end");
    &SOURCE[start..end]
}

#[test]
fn every_method_of_the_svg_writer_is_called_by_the_lock() {
    const WRITER: &str = include_str!("../svg.rs");
    let start = WRITER
        .find("\nimpl SvgWriter {")
        .expect("the writer's methods");
    let end = start + WRITER[start..].find("\n}\n").expect("their end");
    let methods: Vec<&str> = WRITER[start..end]
        .lines()
        .filter_map(|line| line.strip_prefix("    pub fn "))
        .filter_map(|rest| rest.split(['(', '<']).next())
        .collect();
    assert!(
        methods.len() > 30,
        "only {} methods found; the writer has been rearranged and this no longer reads it",
        methods.len()
    );
    let lock = lock_source();
    for method in methods {
        assert!(
            lock.contains(&format!(".{method}(")) || lock.contains(&format!("::{method}(")),
            "SvgWriter::{method} is not called by every_element_the_svg_writer_writes_is_one_\
             the_pdf_reads, so nothing checks the PDF reads what it writes; call it there"
        );
    }
}

/// The other half of the lock: every element and attribute name the
/// writer's source spells in markup, including the branches the lock test
/// cannot reach with ordinary arguments, is one the converter reads.
/// panels.rs is read too, because it is the one other place that writes
/// markup of its own.
#[test]
fn every_name_the_writer_s_source_spells_in_markup_is_one_the_pdf_reads() {
    let mut elements: Vec<String> = Vec::new();
    let mut attributes: Vec<String> = Vec::new();
    for source in [include_str!("../svg.rs"), include_str!("../panels.rs")] {
        let code = &source[..source.find("#[cfg(test)]").unwrap_or(source.len())];
        for line in code.lines() {
            let line = line.trim_start();
            if line.starts_with("//") || !line.contains('"') {
                continue;
            }
            let bytes = line.as_bytes();
            for (at, _) in line.match_indices('<') {
                let before = at.checked_sub(1).map(|b| bytes[b]);
                let starts = bytes.get(at + 1).is_some_and(u8::is_ascii_alphabetic);
                if matches!(before, Some(b'"' | b'>')) && starts {
                    let name: String = line[at + 1..]
                        .chars()
                        .take_while(char::is_ascii_alphanumeric)
                        .collect();
                    if !elements.contains(&name) {
                        elements.push(name);
                    }
                }
            }
            for (at, _) in line.match_indices("=\"") {
                let name: String = line[..at]
                    .chars()
                    .rev()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '-' || *c == ':')
                    .collect::<Vec<char>>()
                    .into_iter()
                    .rev()
                    .collect();
                if !name.is_empty() && !attributes.contains(&name) {
                    attributes.push(name);
                }
            }
        }
    }
    assert!(elements.len() > 10, "{elements:?}");
    assert!(attributes.len() > 20, "{attributes:?}");
    for element in &elements {
        assert!(paint::known(element), "<{element}> is written and not read");
    }
    let readers = [
        "svg",
        "g",
        "rect",
        "circle",
        "line",
        "polyline",
        "polygon",
        "path",
        "text",
        "linearGradient",
        "stop",
        "clipPath",
        "title",
        "desc",
    ];
    for attribute in &attributes {
        assert!(
            readers
                .iter()
                .any(|element| paint::reads(element, attribute)),
            "the {attribute} attribute is written and not read"
        );
    }
}

/// [`drawn`]'s promise: every kind of document the crate makes, empty ones
/// included, has a size to make a page of.
#[test]
fn every_drawing_the_crate_makes_converts_with_nothing_left_out() {
    let region = crate::Region::parse("chr1:1-1000").unwrap();
    let tree = crate::Tree::parse_newick("((A:1,B:2):1,C:3);").unwrap();
    let pdfs = [
        crate::Figure::new(region.clone()).to_pdf(),
        crate::Figure::new(region.clone())
            .push(crate::CoverageTrack::new(0, vec![10.0; 1000]))
            .to_pdf(),
        crate::Panels::new().to_pdf(),
        crate::Rings::new(1_000).to_pdf(),
        crate::Map::new().to_pdf(),
        crate::Map::new()
            .push(crate::GeoLocation::new("Lima", -12.046, -77.043))
            .push(crate::GeoLocation::new("Valencia", 39.47, -0.376))
            .push_flow(crate::GeoFlow::new("Lima", "Valencia"))
            .to_pdf(),
        crate::PhyloMap::new(tree.clone()).to_pdf(),
        crate::plot_tree().add_tree(tree).to_pdf(),
        crate::plot("chr1:1-1000").unwrap().to_pdf(),
    ];
    for (at, pdf) in pdfs.iter().enumerate() {
        assert!(pdf.bytes.starts_with(b"%PDF-1.4\n"), "drawing {at}");
        assert!(pdf.notes.is_empty(), "drawing {at}: {:?}", pdf.notes);
    }
}

#[test]
fn the_same_figure_writes_the_same_bytes() {
    let figure = crate::plot("chr1:1-1000")
        .unwrap()
        .title("twice")
        .add_coverage(vec![20.0; 1000])
        .into_figure();
    let (one, two) = (figure.to_pdf(), figure.to_pdf());
    assert_eq!(one, two);
    // Nothing that changes from one run to the next: no date and no id.
    let written = text(&one);
    for unstable in ["/CreationDate", "/ModDate", "/ID"] {
        assert!(!written.contains(unstable), "{unstable}");
    }
}

/// A PDF that named its faces and embedded none met a journal's preflight
/// with `no` under `emb` for every face. Each face a page uses is now in the
/// file, once, as a subset: a composite font over one TrueType font with its
/// widths, its outlines and a map back to the text, and nothing is named that
/// a reader has to find for itself.
#[test]
fn each_face_a_page_uses_is_embedded_once_as_a_subset() {
    let pdf = crate::plot("chr1:1-1000")
        .unwrap()
        .title("rpoB")
        .add_coverage(vec![20.0; 1000])
        .label("depth")
        .into_figure()
        .to_pdf();
    let faces = embedded(&pdf);
    let mut named: Vec<&str> = faces
        .iter()
        .map(|face| {
            let (tag, name) = face.name.split_once('+').expect("a subset tag");
            assert!(
                tag.len() == 6 && tag.bytes().all(|byte| byte.is_ascii_uppercase()),
                "{}",
                face.name
            );
            name
        })
        .collect();
    named.sort_unstable();
    assert_eq!(
        named,
        ["LiberationMono", "LiberationSans", "LiberationSans-Bold"]
    );
    let written = text(&pdf);
    for (key, count) in [
        ("/Subtype /Type0", 3),
        ("/Encoding /Identity-H", 3),
        ("/Subtype /CIDFontType2", 3),
        ("/CIDToGIDMap /Identity", 3),
        ("/Type /FontDescriptor", 3),
        ("/FontFile2 ", 3),
        ("/ToUnicode ", 3),
        ("/Subtype /Type1", 0),
        ("/BaseFont /Helvetica", 0),
        ("/BaseFont /Courier", 0),
        ("/BaseFont /Symbol", 0),
    ] {
        assert_eq!(written.matches(key).count(), count, "{key}");
    }
    for face in &faces {
        assert!(
            face.file.starts_with(&[0, 1, 0, 0]),
            "{} is not TrueType",
            face.name
        );
    }

    // A face the page does not use is not written: a label in the regular
    // weight alone embeds the regular face alone.
    let mut svg = SvgWriter::new();
    svg.text(10.0, 20.0, "rpoB", "#000", 10.0, Anchor::Start);
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    let faces = embedded(&pdf);
    assert_eq!(faces.len(), 1);
    assert!(
        faces[0].name.ends_with("+LiberationSans"),
        "{}",
        faces[0].name
    );
    assert!(text(&pdf).contains("/Font << /F1 7 0 R >>"));
    // And a page with no text embeds nothing.
    let svg = SvgWriter::new().finish(100.0, 50.0, "#ffffff", "sans-serif");
    assert!(!text(&convert(&svg)).contains("/Font"));
}

/// The two pages that describe a PDF said, plainly, that it named its faces
/// and embedded none, with the Ghostscript line that embedded them, because
/// two readers had met a preflight over it. They are held to the other side
/// now: each names the faces the writer embeds, says they are embedded and
/// how `pdffonts` shows it, and offers no workaround for what is not there.
#[test]
fn the_pages_on_the_pdf_say_its_faces_are_embedded() {
    const GUIDE: &str = include_str!("../../docs/guide/cli.md");
    const HOW: &str = include_str!("../../docs/how-it-works/pdf.md");
    let flat = |text: &str| {
        text.split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .to_lowercase()
    };
    for (page, said) in [("cli.md", flat(GUIDE)), ("pdf.md", flat(HOW))] {
        for needed in [
            "liberation sans",
            "liberation sans bold",
            "liberation mono",
            "embedded",
            "pdffonts",
            "`yes` under `emb`",
        ] {
            assert!(said.contains(needed), "{page} does not say {needed}");
        }
        for gone in [
            "none of them is embedded",
            "no font is embedded",
            "`no` under `emb`",
            "-sdevice=pdfwrite",
            "helvetica-bold for",
        ] {
            assert!(!said.contains(gone), "{page} still says {gone}");
        }
    }
}

/// The subset is the face with every glyph the page does not draw emptied:
/// the glyphs of the text, `.notdef`, and what a composite glyph is built
/// from, here the e and the acute of an é, and no other.
#[test]
fn the_subset_holds_exactly_the_glyphs_the_page_draws() {
    let label = "caf\u{e9} 12";
    let mut svg = SvgWriter::new();
    svg.text(10.0, 20.0, label, "#000", 10.0, Anchor::Start);
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    let faces = embedded(&pdf);
    assert_eq!(faces.len(), 1);
    let source = text::Face::Sans.font();
    let drawn: std::collections::BTreeSet<u16> = label
        .chars()
        .map(|c| source.glyph(c).expect("Liberation Sans has it"))
        .collect();
    let (e, acute) = (source.glyph('e').unwrap(), 118);
    let mut kept = drawn.clone();
    kept.extend([0, e, acute]);
    assert_eq!(source.closure(drawn.iter().copied()), kept);

    let file: &'static [u8] = Box::leak(faces[0].file.clone().into_boxed_slice());
    let subset = font::Font::read(file).expect("the subset reads as TrueType");
    assert_eq!(subset.glyph_count(), kept.last().unwrap() + 1);
    for glyph in 0..subset.glyph_count() {
        let expected: &[u8] = if kept.contains(&glyph) {
            source.outline(glyph)
        } else {
            &[]
        };
        assert_eq!(subset.outline(glyph), expected, "glyph {glyph}");
    }
    // The widths are listed for the glyphs drawn and no others, and they are
    // the face's.
    let listed: Vec<(u16, u16)> = widths(&faces[0].widths);
    let expected: Vec<(u16, u16)> = drawn
        .iter()
        .map(|glyph| (*glyph, source.width(*glyph)))
        .collect();
    assert_eq!(listed, expected);
}

/// What a reader copies, searches and reads aloud is the text, through the
/// map from each glyph back to its character, in every face.
#[test]
fn the_text_copies_back_through_the_map_to_unicode() {
    let labels = [
        "\u{3c9} \u{2264} 0.05",
        "\u{c5}lesund (r\u{b2} \u{2265} 0.8)",
        "6 \u{2192} 3 crossings \u{b7} 6 linked\u{2026}",
    ];
    let mut svg = SvgWriter::new();
    svg.text(10.0, 20.0, labels[0], "#000", 10.0, Anchor::Start);
    svg.text_bold(10.0, 40.0, labels[1], "#000", 10.0, Anchor::Start);
    svg.text_styled(
        10.0,
        60.0,
        labels[2],
        "#000",
        10.0,
        Anchor::Start,
        TextStyle {
            family: Some("JetBrains Mono, monospace"),
            ..TextStyle::default()
        },
    );
    let pdf = convert(&svg.finish(300.0, 80.0, "none", "sans-serif"));
    assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
    let faces = embedded(&pdf);
    assert_eq!(faces.len(), 3);
    assert_eq!(copied(&pdf, &faces), labels);
}

#[test]
fn a_page_is_three_quarters_of_a_point_per_pixel() {
    let svg = SvgWriter::new().finish(900.0, 305.0, "#ffffff", "sans-serif");
    let pdf = convert(&svg);
    let written = text(&pdf);
    assert!(written.contains("/MediaBox [0 0 675 228.75]"), "{written}");
    assert!(written.starts_with("%PDF-1.4\n"));
    assert!(!written.contains("/UserUnit"));
    // Inside, coordinates are the SVG's own, under one matrix that scales
    // and flips them.
    assert!(content(&pdf).contains("0.75 0 0 -0.75 0 228.75 cm\n"));
    assert!(content(&pdf).contains("0 0 900 305 re\nf\n"));
}

#[test]
fn a_page_wider_than_two_hundred_inches_carries_a_user_unit() {
    let svg = SvgWriter::new().finish(100_000.0, 305.0, "#ffffff", "sans-serif");
    let written = text(&convert(&svg));
    // 75,000 points is six units of 12,500.
    assert!(written.starts_with("%PDF-1.6\n"), "{written}");
    assert!(
        written.contains("/MediaBox [0 0 12500 38.125] /Resources"),
        "{written}"
    );
    assert!(written.contains("/UserUnit 6"), "{written}");
    assert!(
        written.contains("0.125 0 0 -0.125 0 38.125 cm"),
        "{written}"
    );
}

#[test]
fn a_document_with_no_size_is_none_and_a_view_box_is_a_size() {
    assert!(Pdf::from_svg("").is_none());
    assert!(Pdf::from_svg("<g/>").is_none());
    assert!(Pdf::from_svg("<svg xmlns=\"http://www.w3.org/2000/svg\"/>").is_none());
    let pdf = convert(r#"<?xml version="1.0"?><svg viewBox="0 0 40 20"/>"#);
    assert!(text(&pdf).contains("/MediaBox [0 0 30 15]"));
    // A width of nought is still a page, as small as a reader allows.
    let pdf = convert(r#"<svg width="0" height="0"/>"#);
    assert!(text(&pdf).contains("/MediaBox [0 0 3 3]"));
}

fn one_label(anchor: Anchor, content_text: &str, size: f64) -> String {
    let mut svg = SvgWriter::new();
    svg.text(100.0, 50.0, content_text, "#000000", size, anchor);
    content(&convert(&svg.finish(200.0, 100.0, "none", "sans-serif")))
}

#[test]
fn text_anchored_at_its_end_ends_at_its_x() {
    // Liberation Sans draws d, e, p and h 556 thousandths of an em wide and
    // t 278, as Helvetica does,
    // so "depth" at 11 pixels is 27.522 wide and ends at 100.
    assert!(one_label(Anchor::End, "depth", 11.0).contains("1 0 0 -1 72.478 50 Tm"));
    // r 333, p 556, o 556, B 667: 21.12 wide at 10, so it starts 10.56 left
    // of its middle.
    assert!(one_label(Anchor::Middle, "rpoB", 10.0).contains("1 0 0 -1 89.44 50 Tm"));
    assert!(one_label(Anchor::Start, "rpoB", 10.0).contains("1 0 0 -1 100 50 Tm"));
}

#[test]
fn a_spaced_label_is_measured_with_its_spacing() {
    let mut svg = SvgWriter::new();
    let style = TextStyle {
        tracking: 0.03,
        ..TextStyle::default()
    };
    svg.text_styled(100.0, 50.0, "AB", "#000", 10.0, Anchor::End, style);
    let stream = content(&convert(&svg.finish(200.0, 100.0, "none", "sans-serif")));
    // A and B are 667 each, 13.34 at 10 pixels, and 0.3 after each letter.
    assert!(stream.contains("0.3 Tc\n"), "{stream}");
    assert!(stream.contains("1 0 0 -1 86.06 50 Tm"), "{stream}");
}

#[test]
fn a_logo_glyph_is_stretched_to_its_box() {
    let mut svg = SvgWriter::new();
    // A is 722 in Liberation Sans Bold: 14.44 at 20 pixels, asked to fill
    // 7.22.
    svg.glyph(10.0, 40.0, 7.22, 20.0, "A", "#109648");
    let stream = content(&convert(&svg.finish(50.0, 50.0, "none", "sans-serif")));
    assert!(stream.contains("50 Tz\n"), "{stream}");
    assert!(stream.contains("/F1 20 Tf\n<0024> Tj"), "{stream}");
}

#[test]
fn a_label_turned_minus_ninety_reads_upwards() {
    let mut svg = SvgWriter::new();
    svg.text_rotated((50.0, 60.0), -90.0, "column", "#000", 9.0, Anchor::Start);
    let stream = content(&convert(&svg.finish(100.0, 100.0, "none", "sans-serif")));
    // The turn is a matrix of its own, which takes the label's x axis to the
    // page's upward direction, and the label itself is only flipped upright.
    assert!(stream.contains("q\n0 -1 1 0 50 60 cm\n"), "{stream}");
    assert!(stream.contains("1 0 0 -1 0 0 Tm"), "{stream}");
    assert!(stream.contains("ET\nQ\n"), "{stream}");
}

#[test]
fn a_haloed_label_is_stroked_under_its_fill() {
    let mut svg = SvgWriter::new();
    svg.text_haloed(
        10.0,
        20.0,
        "rpoB",
        "#111111",
        "#ffffff",
        10.0,
        Anchor::Start,
        false,
    );
    let stream = content(&convert(&svg.finish(100.0, 50.0, "none", "sans-serif")));
    let halo = stream.find("2 Tr\n").expect("a halo");
    let fill = stream.find("0 Tr\n").expect("the letters");
    assert!(halo < fill, "{stream}");
    let before_halo = &stream[..halo];
    assert!(before_halo.contains("1 1 1 RG\n2.8 w\n1 j\n"), "{stream}");
}

#[test]
fn omega_less_or_equal_and_the_ellipsis_are_set_in_the_label_s_own_face() {
    let mut svg = SvgWriter::new();
    svg.text(
        10.0,
        20.0,
        "\u{3c9} \u{2264} 1\u{2026} 2\u{d7}",
        "#000",
        10.0,
        Anchor::End,
    );
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    let stream = content(&pdf);
    // One face, one run: omega, a space, less or equal, a space, one, the
    // ellipsis, a space, two and the multiplication sign.
    assert!(
        stream.contains("/F1 10 Tf\n<037B0003084C0003001407EE000300150099> Tj\n"),
        "{stream}"
    );
    // 781, 278, 549, 278, 556, 1000, 278, 556 and 584 thousandths: 48.6
    // pixels at 10, so the label starts that far left of where it ends.
    assert!(stream.contains("1 0 0 -1 -38.6 20 Tm"), "{stream}");
    assert_eq!(embedded(&pdf).len(), 1);
    assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
}

#[test]
fn a_character_no_embedded_face_has_is_a_question_mark_and_a_note() {
    let mut svg = SvgWriter::new();
    svg.describe("\u{4e2d} sample", "");
    svg.text(10.0, 20.0, "\u{4e2d}", "#000", 10.0, Anchor::Start);
    svg.text(10.0, 30.0, "a\u{4e2d}", "#000", 10.0, Anchor::Start);
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    assert!(content(&pdf).contains("<0022> Tj"));
    assert_eq!(pdf.notes.len(), 1, "{:?}", pdf.notes);
    assert!(pdf.notes[0].contains("U+4E2D"), "{:?}", pdf.notes);
    assert!(
        pdf.notes[0].contains("none of the faces"),
        "{:?}",
        pdf.notes
    );
    // The title keeps the character, since a PDF string can hold any.
    assert!(text(&pdf).contains("/Title <FEFF4E2D002000730061006D0070006C0065>"));
    // Hebrew is in the face, and is marked rather than drawn backwards.
    let mut svg = SvgWriter::new();
    svg.text(10.0, 20.0, "\u{5d0} \u{4e2d}", "#000", 10.0, Anchor::Start);
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    assert_eq!(pdf.notes.len(), 2, "{:?}", pdf.notes);
    assert!(
        pdf.notes[1].contains("U+05D0) is written right to left"),
        "{:?}",
        pdf.notes
    );
}

/// A serif stack and a bold monospaced one ask for faces the PDF does not
/// carry. karyon writes neither, and an SVG of someone else's that does is
/// set in the nearest face there is, with a note.
#[test]
fn a_face_the_pdf_does_not_carry_is_set_in_the_nearest_and_said() {
    let svg = r##"<svg width="100" height="40"><text x="1" y="10" font-family="Georgia, serif">a</text><text x="1" y="20" font-family="Georgia, serif">b</text><text x="1" y="30" font-family="monospace" font-weight="bold">c</text></svg>"##;
    let pdf = convert(svg);
    assert_eq!(pdf.notes.len(), 2, "{:?}", pdf.notes);
    assert!(pdf.notes[0].contains("serif"), "{:?}", pdf.notes);
    assert!(pdf.notes[1].contains("Liberation Mono"), "{:?}", pdf.notes);
    let names: Vec<String> = embedded(&pdf).into_iter().map(|face| face.name).collect();
    assert_eq!(names.len(), 2);
    assert!(names[0].ends_with("+LiberationSans") && names[1].ends_with("+LiberationMono"));
}

#[test]
fn a_fade_is_an_image_whose_alpha_runs_top_to_bottom_clipped_to_its_shape() {
    let mut svg = SvgWriter::new();
    let fade = svg.fade_down("#0072b2", 0.6, 0.05);
    svg.path("M0 10 L100 10 L100 90 L0 90 Z", &fade, 1.0);
    // The same paint again is the same image.
    svg.path("M0 20 L50 20 L50 90 L0 90 Z", &fade, 1.0);
    let pdf = convert(&svg.finish(100.0, 100.0, "none", "sans-serif"));
    let stream = content(&pdf);
    assert_eq!(stream.matches("/I1 Do").count(), 2, "{stream}");
    assert!(!stream.contains("/I2 Do"));
    // Clipped to the path, then stretched over its box, 100 by 80 at 0,10.
    assert!(
        stream.contains("0 10 m\n100 10 l\n100 90 l\n0 90 l\nh\nW n\n100 0 0 80 0 10 cm\n"),
        "{stream}"
    );
    // Which way up: an image's first row is drawn along the top of its unit
    // square, and the two matrices before it take that row to the top of the
    // box, y = 10, where the fade is strongest, and its last row to the foot.
    let placed = image_placement(&stream, "/I1 Do");
    for (x, y, lands) in [(0.0, 1.0, 10.0), (1.0, 1.0, 10.0), (0.0, 0.0, 90.0)] {
        let (_, at) = placed.apply(x, y);
        assert!(
            (at - lands).abs() < 1e-9,
            "({x}, {y}) lands at {at}: {stream}"
        );
    }
    let alpha = stream_after(&pdf, "/ColorSpace /DeviceGray");
    assert_eq!(alpha.len(), 256);
    // The top row is the top's opacity, the last the foot's, falling all
    // the way.
    assert_eq!(
        alpha[0],
        (255.0 * (0.6 - 0.55 * 0.5 / 256.0_f64)).round() as u8
    );
    assert_eq!(
        alpha[255],
        (255.0 * (0.05 + 0.55 * 0.5 / 256.0_f64)).round() as u8
    );
    assert!(alpha.windows(2).all(|pair| pair[0] >= pair[1]));
    // Its colour is the fade's throughout.
    let rgb = stream_after(&pdf, "/ColorSpace /DeviceRGB");
    assert_eq!(rgb.len(), 3 * 256);
    assert!(rgb.chunks(3).all(|pixel| pixel == [0x00, 0x72, 0xb2]));
    assert!(text(&pdf).contains("/Group << /S /Transparency /CS /DeviceRGB >>"));
}

/// The map from an image's unit square to the page's SVG coordinates: the
/// `cm` matrices written between the last `q` before `draw` and `draw`
/// itself, the later one applied first, as a PDF reader applies them.
fn image_placement(stream: &str, draw: &str) -> path::Matrix {
    let end = stream
        .find(draw)
        .unwrap_or_else(|| panic!("no {draw}: {stream}"));
    let start = stream[..end].rfind("q\n").expect("the image's own q");
    let mut placed = path::Matrix::IDENTITY;
    for line in stream[start..end].lines() {
        if let Some(numbers) = line.strip_suffix(" cm") {
            let values: Vec<f64> = numbers
                .split(' ')
                .map(|number| number.parse().expect("a number"))
                .collect();
            let matrix = path::Matrix(values.try_into().expect("six numbers"));
            placed = placed.then_inner(&matrix);
        }
    }
    placed
}

#[test]
fn a_fade_with_no_height_draws_nothing() {
    let mut svg = SvgWriter::new();
    let fade = svg.fade_down("#0072b2", 0.6, 0.05);
    svg.path("M0 10 L100 10 Z", &fade, 1.0);
    let pdf = convert(&svg.finish(100.0, 100.0, "none", "sans-serif"));
    assert!(!content(&pdf).contains(" Do"));
    assert!(pdf.notes.is_empty());
}

#[test]
fn a_clip_ends_where_its_group_ends() {
    let mut svg = SvgWriter::new();
    svg.begin_clip(0.0, 0.0, 10.0, 10.0);
    svg.rect(0.0, 0.0, 20.0, 20.0, "#ff0000");
    svg.end_group();
    svg.rect(0.0, 0.0, 30.0, 30.0, "#0000ff");
    let stream = content(&convert(&svg.finish(40.0, 40.0, "none", "sans-serif")));
    let clip = stream.find("0 0 10 10 re\nW n\n").expect("the clip");
    let inside = stream
        .find("0 0 20 20 re\nf\n")
        .expect("the clipped rectangle");
    let closed = clip + stream[clip..].find("Q\n").expect("the clip's end");
    let outside = stream
        .find("0 0 30 30 re\nf\n")
        .expect("the rectangle after");
    assert!(stream[..clip].ends_with("q\n"), "{stream}");
    assert!(
        clip < inside && inside < closed && closed < outside,
        "{stream}"
    );
}

/// The state a `Q` restores is the one from before its `q`, so whatever was
/// set inside is set again after it. The writer keeps a copy of that state
/// to skip an operator that would change nothing, and if the copy outlived
/// the `Q`, a mark after a group would take the paint of the group's last
/// mark with no operator to say so.
#[test]
fn what_a_group_sets_is_set_again_after_its_group_ends() {
    let mark = r##"<rect x="10" width="5" height="5" fill="#ff0000" fill-opacity="0.5" stroke="#0000ff" stroke-width="2" stroke-dasharray="1 1"/><text x="1" y="15" font-size="7" letter-spacing="1">a</text>"##;
    for group in [r#"clip-path="url(#c)""#, r#"transform="translate(1 1)""#] {
        let svg = format!(
            r#"<svg width="20" height="20"><defs><clipPath id="c"><rect width="10" height="10"/></clipPath></defs><g {group}>{mark}</g>{mark}</svg>"#
        );
        let pdf = convert(&svg);
        assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
        let stream = content(&pdf);
        let restored = stream.find("\nQ\n").expect("the group's end");
        let (inside, after) = stream.split_at(restored);
        for operator in [
            "1 0 0 rg\n",
            "/G1 gs\n",
            "0 0 1 RG\n",
            "2 w\n",
            "[1 1] 0 d\n",
            "/F1 7 Tf\n",
            "1 Tc\n",
        ] {
            assert!(
                inside.contains(operator),
                "{operator:?} inside {group}: {stream}"
            );
            assert!(
                after.contains(operator),
                "{operator:?} after {group}: {stream}"
            );
        }
    }
}

#[test]
fn a_panel_is_clipped_to_its_own_viewport() {
    let figure = crate::plot("chr1:1-1000")
        .unwrap()
        .add_coverage(vec![30.0; 1000])
        .into_figure();
    let (width, height) = figure.dimensions();
    let sheet = crate::Panels::new().push(&figure, "A");
    let stream = content(&sheet.to_pdf());
    let viewport = format!(
        "q\n0 0 {} {} re\nW n\n",
        crate::svg::num(width),
        crate::svg::num(height)
    );
    assert!(stream.contains(&viewport), "{viewport} in {stream}");
}

#[test]
fn a_stroke_of_width_nought_draws_nothing() {
    let mut svg = SvgWriter::new();
    svg.rect_outline(1.0, 1.0, 10.0, 10.0, "#000000", 0.0);
    svg.line(0.0, 0.0, 10.0, 10.0, "#000000", 0.0);
    let stream = content(&convert(&svg.finish(20.0, 20.0, "none", "sans-serif")));
    assert!(
        !stream.contains("\nS\n") && !stream.contains("\nB\n"),
        "{stream}"
    );
    // A width that is drawn is stroked.
    let mut svg = SvgWriter::new();
    svg.rect_outline(1.0, 1.0, 10.0, 10.0, "#000000", 0.5);
    let stream = content(&convert(&svg.finish(20.0, 20.0, "none", "sans-serif")));
    assert!(stream.contains("0.5 w\n1 1 10 10 re\nS\n"), "{stream}");
}

#[test]
fn every_spelling_of_a_colour_paints_the_same() {
    for spelling in ["blue", "#00f", "#0000ff", "rgb(0,0,255)", "rgb(0%,0%,100%)"] {
        let mut svg = SvgWriter::new();
        svg.rect(0.0, 0.0, 5.0, 5.0, spelling);
        let pdf = convert(&svg.finish(10.0, 10.0, "none", "sans-serif"));
        assert!(content(&pdf).contains("0 0 1 rg\n"), "{spelling}");
        assert!(pdf.notes.is_empty(), "{spelling}: {:?}", pdf.notes);
    }
    // An alpha in the colour is an alpha on the page.
    let mut svg = SvgWriter::new();
    svg.rect(0.0, 0.0, 5.0, 5.0, "rgba(0, 0, 255, 0.5)");
    let pdf = convert(&svg.finish(10.0, 10.0, "none", "sans-serif"));
    assert!(content(&pdf).contains("/G1 gs\n"));
    assert!(text(&pdf).contains("<< /Type /ExtGState /ca 0.5 /CA 1 >>"));
}

#[test]
fn an_unreadable_colour_fills_black_and_strokes_nothing_as_a_browser_does() {
    let mut svg = SvgWriter::new();
    svg.rect(0.0, 0.0, 5.0, 5.0, "#00f");
    svg.rect(0.0, 0.0, 6.0, 6.0, "chartreuse-ish");
    svg.line(0.0, 0.0, 9.0, 9.0, "lab(50 0 0)", 1.0);
    let pdf = convert(&svg.finish(10.0, 10.0, "none", "sans-serif"));
    let stream = content(&pdf);
    // Black, the fill a shape inherits from the page.
    assert!(
        stream.contains("0 0 1 rg\n0 0 5 5 re\nf\n0 0 0 rg\n0 0 6 6 re\nf\n"),
        "{stream}"
    );
    assert!(!stream.contains("\nS\n"), "{stream}");
    assert_eq!(pdf.notes.len(), 2, "{:?}", pdf.notes);
    assert!(pdf.notes[0].contains("chartreuse-ish"));
    assert!(pdf.notes[1].contains("lab(50 0 0)"));
}

#[test]
fn a_reference_to_a_missing_gradient_paints_nothing() {
    let mut svg = SvgWriter::new();
    svg.path("M0 0 L10 0 L10 10 Z", "url(#nowhere)", 1.0);
    let pdf = convert(&svg.finish(10.0, 10.0, "none", "sans-serif"));
    assert!(!content(&pdf).contains("\nf\n"));
    assert_eq!(pdf.notes.len(), 1, "{:?}", pdf.notes);
    assert!(pdf.notes[0].contains("url(#nowhere)"), "{:?}", pdf.notes);
}

#[test]
fn the_title_and_description_become_document_properties_and_alt_text() {
    let mut svg = SvgWriter::new();
    svg.describe("rpoB", "A figure");
    svg.begin_titled("a tooltip, which is not the title");
    svg.end_group();
    let written = text(&convert(&svg.finish(10.0, 10.0, "none", "sans-serif")));
    assert!(
        written.contains("/Title <FEFF00720070006F0042>"),
        "{written}"
    );
    assert!(written.contains("/Subject <FEFF00410020006600690067007500720065>"));
    // What aria-labelledby reads out: the title, then the description.
    assert!(
        written.contains(&format!("/Alt {}", utf16("rpoB. A figure"))),
        "{written}"
    );
    assert!(written.contains("/MarkInfo << /Marked true >>"));
    assert!(written.contains("/Figure << /MCID 0 >> BDC"));
}

#[test]
fn what_the_converter_does_not_read_is_named_once() {
    let svg = r##"<svg width="10" height="10"><rect width="5" height="5" opacity="0.5"/><rect width="5" height="5" opacity="0.5"/><foreignObject><p>x</p></foreignObject><text x="1" y="5">a<tspan>b</tspan></text><g transform="wobble(1)"/></svg>"##;
    let pdf = convert(svg);
    assert_eq!(pdf.notes.len(), 4, "{:#?}", pdf.notes);
    assert!(pdf.notes[0].contains("opacity attribute of <rect>"));
    assert!(pdf.notes[1].contains("<foreignObject>"));
    assert!(pdf.notes[2].contains("<tspan>"));
    assert!(pdf.notes[3].contains("wobble(1)"));
    // The rest is drawn: both rectangles, and the text before the tspan.
    let stream = content(&pdf);
    assert_eq!(stream.matches("0 0 5 5 re\nf\n").count(), 2, "{stream}");
    assert!(stream.contains("<0044> Tj"), "{stream}");
}

#[test]
fn a_group_s_paint_is_inherited_by_what_is_inside_it() {
    let svg = r##"<svg width="10" height="10"><g fill="#ff0000" stroke="#0000ff" stroke-width="2"><rect width="5" height="5"/><g fill="none"><rect width="6" height="6"/></g></g></svg>"##;
    let stream = content(&convert(svg));
    assert!(
        stream.contains("1 0 0 rg\n0 0 1 RG\n2 w\n0 0 5 5 re\nB\n"),
        "{stream}"
    );
    assert!(stream.contains("0 0 6 6 re\nS\n"), "{stream}");
}

#[test]
fn a_length_in_any_absolute_unit_is_read_by_the_drawing_as_by_the_page() {
    // A page 100 by 50 millimetres around a view 100 by 50: each unit of the
    // view is a millimetre, 3.78 pixels, so the rectangle fills the page.
    let pdf = convert(
        r##"<svg width="100mm" height="50mm" viewBox="0 0 100 50"><rect width="100" height="50" fill="#ff0000"/></svg>"##,
    );
    assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
    assert!(text(&pdf).contains("/MediaBox [0 0 283.465 141.732]"));
    let stream = content(&pdf);
    assert!(
        stream.contains("0.75 0 0 -0.75 0 141.732 cm\n4 M\n3.779528 0 0 3.779528 0 0 cm\n"),
        "{stream}"
    );
    // A size in points, and a rectangle in millimetres, inside a page in
    // pixels: 30 points is 40 pixels, and 10 millimetres 37.795.
    let pdf = convert(
        r##"<svg width="200" height="50"><text x="0" y="40" font-size="30pt">Hi</text><rect x="100" y="0" width="10mm" height="10mm" fill="red"/></svg>"##,
    );
    assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
    let stream = content(&pdf);
    assert!(stream.contains("/F1 40 Tf\n<002B004C> Tj"), "{stream}");
    assert!(stream.contains("100 0 37.795 37.795 re\nf\n"), "{stream}");
}

#[test]
fn a_length_measured_against_a_box_or_a_font_is_named_rather_than_dropped() {
    let pdf = convert(
        r##"<svg width="100%" height="50" viewBox="0 0 100 50"><rect width="50%" height="10"/><text x="1em" y="5" font-size="2em" stroke="#000" stroke-dasharray="1ex">a</text></svg>"##,
    );
    let said = pdf.notes.join("\n");
    for value in [
        "\"100%\" of <svg>",
        "\"50%\" of <rect>",
        "\"2em\"",
        "\"1em\" of <text>",
    ] {
        assert!(said.contains(value), "{value} in {said}");
    }
    assert!(said.contains("stroke-dasharray \"1ex\""), "{said}");
    assert_eq!(pdf.notes.len(), 5, "{said}");
    // The page is the view's size, as the drawing is.
    assert!(text(&pdf).contains("/MediaBox [0 0 75 37.5]"));
}

/// Every character the measure has a width of its own for: printable ASCII,
/// and the characters beyond it that karyon writes.
fn measured_characters() -> Vec<char> {
    (0x20u8..0x7f)
        .map(char::from)
        .chain(crate::svg::BEYOND_ASCII.iter().map(|(c, _, _)| *c))
        .collect()
}

/// How wide the PDF draws `letter` in `face`, in thousandths of an em.
fn drawn(letter: char, face: text::Face) -> f64 {
    let mut missing = Vec::new();
    let set = text::runs(&letter.to_string(), face, &mut missing);
    assert!(missing.is_empty(), "{letter:?} is in no embedded face");
    assert!(
        set.iter().all(|run| run.face == face),
        "{letter:?} is not in {face:?}"
    );
    set.iter().map(|run| run.width).sum()
}

/// The layout was made with Adobe's widths for Helvetica and Helvetica-Bold,
/// and Courier's 600, which is what a PDF set its text in before it embedded
/// its faces, with Symbol's for the five characters WinAnsi lacks. Liberation
/// is drawn to the same widths, so a label sits where it sat. Every character
/// karyon measures is held to that here, from the faces' own `hmtx` rounded
/// to thousandths as Adobe's are, and the three that differ are named with
/// both widths: the middle dot, which Liberation draws as Arial does rather
/// than as Helvetica does, and omega and the arrow, which are Liberation's
/// own rather than Symbol's. The tests after this one hold the layout's
/// measure to be at least as wide as each, which it is, so no table changes.
#[test]
fn the_embedded_faces_draw_every_measured_character_at_the_layout_s_widths() {
    // Each that differs, with Adobe's regular and bold widths and then
    // Liberation Sans's and Liberation Sans Bold's.
    let differs = [
        ('\u{b7}', (278, 278), (333, 333)),
        ('\u{3c9}', (686, 686), (781, 845)),
        ('\u{2192}', (987, 987), (1000, 1000)),
    ];
    let mut checked = 0;
    for c in measured_characters() {
        let adobe = match metrics::winansi(c) {
            Some(byte) => {
                let at = usize::from(byte - 32);
                (metrics::HELVETICA[at], metrics::HELVETICA_BOLD[at])
            }
            // Symbol has one weight.
            None => {
                let (_, width) = metrics::SYMBOL
                    .iter()
                    .find(|(known, _)| *known == c)
                    .unwrap_or_else(|| panic!("{c:?} has no Adobe width"));
                (*width, *width)
            }
        };
        let width = |face: text::Face| {
            let font = face.font();
            font.width(
                font.glyph(c)
                    .unwrap_or_else(|| panic!("{c:?} is not in {face:?}")),
            )
        };
        let liberation = (width(text::Face::Sans), width(text::Face::SansBold));
        match differs.iter().find(|(known, _, _)| *known == c) {
            Some((_, was, is)) => assert_eq!((adobe, liberation), (*was, *is), "{c:?}"),
            None => assert_eq!(liberation, adobe, "{c:?}: Liberation against Adobe"),
        }
        assert_eq!(width(text::Face::Mono), 600, "{c:?} in Liberation Mono");
        checked += 1;
    }
    assert_eq!(checked, 95 + 10);
}

/// The PDF sets regular text in Liberation Sans, and the layout was made
/// with [`text_width`](crate::svg::text_width): a measure narrower than the
/// face for any character would let a label that fitted in the SVG run past
/// its room on the page.
#[test]
fn the_regular_measure_is_never_narrower_than_liberation_sans() {
    for letter in measured_characters() {
        let measured = crate::svg::text_width(&letter.to_string(), 1000.0);
        let drawn = drawn(letter, text::Face::Sans);
        assert!(
            measured >= drawn - 1e-9,
            "{letter:?} is measured {measured} and Liberation Sans draws it {drawn}"
        );
    }
}

/// Bold and semibold text is set in Liberation Sans Bold and measured with
/// [`text_width_strong`](crate::svg::text_width_strong), which must not be
/// the narrower of the two for any character, or a bold panel letter "B" runs
/// past the room made for it.
#[test]
fn the_strong_measure_is_never_narrower_than_liberation_sans_bold() {
    for letter in measured_characters() {
        let measured = crate::svg::text_width_strong(&letter.to_string(), 1000.0);
        let drawn = drawn(letter, text::Face::SansBold);
        assert!(
            measured >= drawn - 1e-9,
            "{letter:?} is measured {measured} and Liberation Sans Bold draws it {drawn}"
        );
    }
}

/// The monospaced stack is set in Liberation Mono, every glyph of which is
/// 600 thousandths of an em, and measured at the same 600 by
/// [`mono_width`](crate::svg::mono_width), omega and the relations included.
#[test]
fn the_monospaced_measure_is_liberation_mono() {
    let stack = "JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace";
    assert_eq!(
        text::Face::choose(stack, "normal"),
        (text::Face::Mono, None)
    );
    for letter in measured_characters() {
        let measured = crate::svg::mono_width(&letter.to_string(), 1000.0);
        let drawn = drawn(letter, text::Face::Mono);
        assert!(
            (measured - drawn).abs() < 1e-9 && (drawn - 600.0).abs() < 1e-9,
            "{letter:?} is measured {measured} and Liberation Mono draws it {drawn}"
        );
    }
}

/// One embedded face, read back out of the file.
struct Embedded {
    /// The name the page gives it, subset tag and all.
    name: String,
    /// What is between the brackets of its `/W`.
    widths: String,
    /// Its `/FontFile2`.
    file: Vec<u8>,
    /// Its `/ToUnicode`, from glyph to text.
    unicode: std::collections::BTreeMap<u16, String>,
}

/// Where `needle` next is in `haystack`, from `from` on.
fn find(haystack: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    haystack
        .get(from..)?
        .windows(needle.len())
        .position(|window| window == needle)
        .map(|at| from + at)
}

/// Object `number`'s dictionary, as text.
fn dictionary(pdf: &Pdf, number: usize) -> String {
    let start = find(&pdf.bytes, format!("\n{number} 0 obj\n").as_bytes(), 0)
        .unwrap_or_else(|| panic!("no object {number}"));
    let end = find(&pdf.bytes, b"\nendobj", start).expect("its end");
    let body = String::from_utf8_lossy(&pdf.bytes[start..end]);
    body[..body.find("\nstream\n").unwrap_or(body.len())].to_string()
}

/// The object `key` refers to in `dictionary`, as in `/FontFile2 12 0 R`.
fn reference(dictionary: &str, key: &str) -> usize {
    let after =
        &dictionary[dictionary.find(key).unwrap_or_else(|| panic!("no {key}")) + key.len()..];
    after
        .trim_start_matches([' ', '['])
        .split(' ')
        .next()
        .and_then(|number| number.parse().ok())
        .unwrap_or_else(|| panic!("no number after {key}"))
}

/// Object `number`'s stream, by the length its dictionary gives.
fn object_stream(pdf: &Pdf, number: usize) -> Vec<u8> {
    let length = reference(&dictionary(pdf, number), "/Length ");
    let start = find(&pdf.bytes, format!("\n{number} 0 obj\n").as_bytes(), 0).unwrap();
    let begin = find(&pdf.bytes, b"stream\n", start).unwrap() + b"stream\n".len();
    pdf.bytes[begin..begin + length].to_vec()
}

/// Every composite font in the file, in the order the page numbers them.
fn embedded(pdf: &Pdf) -> Vec<Embedded> {
    let marker = b"/Subtype /Type0 /BaseFont /";
    let mut faces = Vec::new();
    let mut from = 0;
    while let Some(at) = find(&pdf.bytes, marker, from) {
        from = at + marker.len();
        // The object's number is on the line before its dictionary.
        let opened = pdf.bytes[..at]
            .windows(7)
            .rposition(|window| window == b" 0 obj\n")
            .expect("the font's object");
        let line = pdf.bytes[..opened]
            .iter()
            .rposition(|byte| *byte == b'\n')
            .unwrap()
            + 1;
        let number: usize = String::from_utf8_lossy(&pdf.bytes[line..opened])
            .parse()
            .expect("the font's object number");
        let font = dictionary(pdf, number);
        let name = font[font.find("/BaseFont /").unwrap() + 11..]
            .split(' ')
            .next()
            .unwrap()
            .to_string();
        let descendant = dictionary(pdf, reference(&font, "/DescendantFonts"));
        let widths = descendant[descendant.find("/W [").expect("/W") + 4..]
            .split(" /CIDToGIDMap")
            .next()
            .unwrap()
            .trim_end_matches(']')
            .to_string();
        let descriptor = dictionary(pdf, reference(&descendant, "/FontDescriptor"));
        let file = object_stream(pdf, reference(&descriptor, "/FontFile2"));
        let map = object_stream(pdf, reference(&font, "/ToUnicode"));
        let map = String::from_utf8(map).expect("a CMap is text");
        let mut unicode = std::collections::BTreeMap::new();
        // Every line of a pair, less the one range of codes there is.
        let pairs = map
            .lines()
            .filter(|line| line.starts_with('<') && *line != "<0000> <FFFF>");
        for line in pairs {
            let (glyph, text) = line.split_once("> <").expect("a pair");
            let glyph = u16::from_str_radix(&glyph[1..], 16).expect("a glyph");
            let units: Vec<u16> = text
                .trim_end_matches('>')
                .as_bytes()
                .chunks(4)
                .map(|unit| u16::from_str_radix(std::str::from_utf8(unit).unwrap(), 16).unwrap())
                .collect();
            unicode.insert(glyph, String::from_utf16(&units).expect("UTF-16"));
        }
        faces.push(Embedded {
            name,
            widths,
            file,
            unicode,
        });
    }
    faces
}

/// The glyphs a `/W` array lists, each with its width.
fn widths(array: &str) -> Vec<(u16, u16)> {
    let mut listed = Vec::new();
    for group in array.split(']').filter(|group| !group.trim().is_empty()) {
        let (first, list) = group.split_once('[').expect("a first glyph and its widths");
        let first: u16 = first.trim().parse().unwrap();
        for (at, width) in list.split_whitespace().enumerate() {
            listed.push((first + at as u16, width.parse().unwrap()));
        }
    }
    listed
}

/// The text of each label, as a reader copies it: every glyph the content
/// stream draws, through the map of the face it is drawn in.
fn copied(pdf: &Pdf, faces: &[Embedded]) -> Vec<String> {
    let stream = content(pdf);
    let mut labels = Vec::new();
    let mut face = 0;
    for line in stream.lines() {
        if line == "BT" {
            labels.push(String::new());
        } else if let Some(rest) = line.strip_prefix("/F").filter(|_| line.ends_with(" Tf")) {
            face = rest.split(' ').next().unwrap().parse::<usize>().unwrap() - 1;
        } else if let Some(hex) = line
            .strip_suffix("> Tj")
            .and_then(|line| line.strip_prefix('<'))
        {
            for glyph in hex.as_bytes().chunks(4) {
                let glyph = u16::from_str_radix(std::str::from_utf8(glyph).unwrap(), 16).unwrap();
                let text = faces[face]
                    .unicode
                    .get(&glyph)
                    .expect("every glyph is mapped");
                labels.last_mut().expect("inside BT").push_str(text);
            }
        }
    }
    labels
}
