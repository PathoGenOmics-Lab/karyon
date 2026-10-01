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
        svg.text_styled(
            10.0,
            10.0,
            "761,000 r\u{b2}",
            "#111111",
            9.0,
            anchor,
            TextStyle {
                family: Some("JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace"),
                weight: Some(600),
                tracking: 0.03,
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
    // Helvetica's d, e, p and h are 556 thousandths of an em and t is 278,
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
    // A is 722 in Helvetica-Bold: 14.44 at 20 pixels, asked to fill 7.22.
    svg.glyph(10.0, 40.0, 7.22, 20.0, "A", "#109648");
    let stream = content(&convert(&svg.finish(50.0, 50.0, "none", "sans-serif")));
    assert!(stream.contains("50 Tz\n"), "{stream}");
    assert!(stream.contains("/F1 20 Tf\n(A) Tj"), "{stream}");
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
fn omega_and_less_or_equal_are_set_in_symbol_and_the_ellipsis_in_winansi() {
    let mut svg = SvgWriter::new();
    svg.text(
        10.0,
        20.0,
        "\u{3c9} \u{2264} 1\u{2026} 2\u{d7}",
        "#000",
        10.0,
        Anchor::Start,
    );
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    let stream = content(&pdf);
    assert!(stream.contains("/F1 10 Tf\n(w) Tj\n"), "{stream}");
    assert!(stream.contains("/F1 10 Tf\n(\\243) Tj\n"), "{stream}");
    assert!(stream.contains("( 1\\205 2\\327) Tj"), "{stream}");
    let written = text(&pdf);
    assert!(written.contains("/BaseFont /Symbol >>"), "{written}");
    assert!(written.contains("/BaseFont /Helvetica /Encoding /WinAnsiEncoding"));
    assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
}

#[test]
fn a_character_no_builtin_font_has_is_a_question_mark_and_a_note() {
    let mut svg = SvgWriter::new();
    svg.describe("\u{4e2d} sample", "");
    svg.text(10.0, 20.0, "\u{4e2d}", "#000", 10.0, Anchor::Start);
    svg.text(10.0, 30.0, "a\u{4e2d}", "#000", 10.0, Anchor::Start);
    let pdf = convert(&svg.finish(100.0, 50.0, "none", "sans-serif"));
    assert!(content(&pdf).contains("(?) Tj"));
    assert_eq!(pdf.notes.len(), 1, "{:?}", pdf.notes);
    assert!(pdf.notes[0].contains("U+4E2D"), "{:?}", pdf.notes);
    // The title keeps the character, since a PDF string can hold any.
    assert!(text(&pdf).contains("/Title <FEFF4E2D002000730061006D0070006C0065>"));
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
    assert!(stream.contains("(a) Tj"), "{stream}");
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
