//! The PDF a user gets, read back the way a reader reads one.
//!
//! A PDF that a forgiving reader opens can still be wrong: Ghostscript and
//! poppler both repair a cross-reference table whose offsets are off, quietly,
//! and a stricter reader, or a journal's preflight, refuses the same file. So
//! the structure is checked here by a reader that repairs nothing: every
//! offset the table gives lands on its object, every stream is as long as it
//! says, every `q` has its `Q`, and every font, alpha and image the drawing
//! names is one the page defines. The checker is tested too, on files broken
//! on purpose, since a check that passes everything says nothing.

use std::fs;
use std::path::{Path, PathBuf};

use karyon::{
    AxisTrack, CoverageTrack, Feature, FeatureTrack, Figure, Panels, Pdf, Region, Strand, Variant,
    VariantTrack,
};

fn demo_figure() -> Figure {
    let region = Region::parse("NC_000962.3:761001-763000").unwrap();
    let depth: Vec<f64> = (0..2000).map(|i| 40.0 + (i % 31) as f64).collect();
    Figure::new(region)
        .title("integration figure, \u{3c9} \u{2264} 1")
        .push(CoverageTrack::new(761_000, depth).label("depth"))
        .push(
            FeatureTrack::new(vec![Feature::new(761_100, 761_900)
                .name("geneA")
                .strand(Strand::Forward)])
            .label("genes"),
        )
        .push(
            VariantTrack::new(vec![Variant::new(761_250).value(0.9).category("missense")])
                .label("variants"),
        )
        .push(AxisTrack::new())
}

fn find(bytes: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    bytes
        .get(from..)?
        .windows(needle.len())
        .position(|window| window == needle)
        .map(|at| at + from)
}

fn ascii(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The tokens of a content stream: string literals whole, everything else
/// split at whitespace.
fn tokens(content: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut chars = content.chars().peekable();
    while let Some(c) = chars.next() {
        if c.is_whitespace() {
            continue;
        }
        if c == '(' {
            let (mut depth, mut string) = (1, String::from("("));
            while let Some(c) = chars.next() {
                string.push(c);
                match c {
                    '\\' => string.extend(chars.next()),
                    '(' => depth += 1,
                    ')' => {
                        depth -= 1;
                        if depth == 0 {
                            break;
                        }
                    }
                    _ => {}
                }
            }
            out.push(string);
            continue;
        }
        let mut token = String::from(c);
        while let Some(&next) = chars.peek() {
            if next.is_whitespace() || next == '(' {
                break;
            }
            token.push(next);
            chars.next();
        }
        out.push(token);
    }
    out
}

/// Everything wrong with the structure of `pdf`, or nothing.
fn problems(pdf: &[u8]) -> Vec<String> {
    let mut found = Vec::new();
    if !pdf.starts_with(b"%PDF-1.") {
        found.push("no %PDF header".to_string());
    }
    if !pdf.ends_with(b"%%EOF\n") {
        found.push("no %%EOF at the end".to_string());
    }
    // The table, from where startxref says it is.
    let Some(start) = find(pdf, b"startxref\n", 0) else {
        found.push("no startxref".to_string());
        return found;
    };
    let after = ascii(&pdf[start + 10..]);
    let Some(xref) = after
        .lines()
        .next()
        .and_then(|line| line.parse::<usize>().ok())
    else {
        found.push("startxref names no offset".to_string());
        return found;
    };
    if !pdf
        .get(xref..)
        .is_some_and(|rest| rest.starts_with(b"xref\n0 "))
    {
        found.push(format!("startxref {xref} is not the table"));
        return found;
    }
    let header_end = find(pdf, b"\n", xref + 5).unwrap();
    let count: usize = ascii(&pdf[xref + 7..header_end]).parse().unwrap_or(0);
    let entries = header_end + 1;
    let mut offsets = Vec::new();
    for number in 0..count {
        let Some(entry) = pdf.get(entries + 20 * number..entries + 20 * (number + 1)) else {
            found.push(format!("the table stops before entry {number}"));
            return found;
        };
        let entry = ascii(entry);
        if number == 0 {
            if entry != "0000000000 65535 f \n" {
                found.push(format!("entry 0 is {entry:?}"));
            }
            continue;
        }
        if entry.len() != 20 || !entry.ends_with(" 00000 n \n") {
            found.push(format!("entry {number} is {entry:?}"));
            continue;
        }
        let offset: usize = entry[..10].parse().unwrap_or(usize::MAX);
        let opens = format!("{number} 0 obj\n");
        if !pdf
            .get(offset..)
            .is_some_and(|rest| rest.starts_with(opens.as_bytes()))
        {
            found.push(format!("object {number} is not at {offset}"));
        }
        offsets.push(offset);
    }
    let trailer = ascii(&pdf[entries + 20 * count..]);
    if !trailer.starts_with(&format!("trailer\n<< /Size {count} ")) {
        found.push(format!("the trailer does not say /Size {count}"));
    }
    // Every object's body, by number, for what the streams refer to.
    let body = |number: usize| -> Option<(usize, usize)> {
        let offset = *offsets.get(number.checked_sub(1)?)?;
        let start = offset + format!("{number} 0 obj\n").len();
        let end = find(pdf, b"\nendobj\n", start)?;
        Some((start, end))
    };
    for number in 1..count {
        let Some((start, end)) = body(number) else {
            found.push(format!("object {number} does not end"));
            continue;
        };
        let Some(stream) = find(&pdf[..end], b">>\nstream\n", start) else {
            continue;
        };
        let dictionary = ascii(&pdf[start..stream]);
        let words: Vec<&str> = dictionary.split_whitespace().collect();
        let Some(at) = words.iter().position(|word| *word == "/Length") else {
            found.push(format!("stream {number} has no /Length"));
            continue;
        };
        let length = if words.get(at + 3) == Some(&"R") {
            let target: usize = words[at + 1].parse().unwrap_or(0);
            body(target).and_then(|(s, e)| ascii(&pdf[s..e]).trim().parse::<usize>().ok())
        } else {
            words[at + 1].parse::<usize>().ok()
        };
        let data = stream + ">>\nstream\n".len();
        let ends = length.map(|length| data + length);
        if ends
            .and_then(|ends| pdf.get(ends..))
            .map_or(true, |rest| !rest.starts_with(b"\nendstream"))
        {
            found.push(format!("stream {number} is not {length:?} bytes long"));
        }
    }
    // The drawing: balanced, and naming only what the page defines.
    let Some((start, end)) = body(4) else {
        found.push("no content stream".to_string());
        return found;
    };
    let content_start = find(pdf, b"stream\n", start).unwrap() + 7;
    let content_end = find(pdf, b"\nendstream", content_start).unwrap().min(end);
    let content = ascii(&pdf[content_start..content_end]);
    let resources = body(6).map_or_else(String::new, |(s, e)| ascii(&pdf[s..e]));
    let mut depth = [0i64; 3];
    let pairs = [("q", "Q"), ("BT", "ET"), ("BDC", "EMC")];
    let words = tokens(&content);
    for (at, word) in words.iter().enumerate() {
        for (kind, (open, close)) in pairs.iter().enumerate() {
            if word == open {
                depth[kind] += 1;
            }
            if word == close {
                depth[kind] -= 1;
                if depth[kind] < 0 {
                    found.push(format!("{close} with no {open} at token {at}"));
                }
            }
        }
        let named = match word.as_str() {
            "Tf" => at.checked_sub(2),
            "gs" | "Do" => at.checked_sub(1),
            _ => None,
        };
        if let Some(name) = named.map(|n| &words[n]) {
            if !resources.contains(&format!("{name} ")) {
                found.push(format!("{name} is used and not defined"));
            }
        }
    }
    for (kind, (open, close)) in pairs.iter().enumerate() {
        if depth[kind] != 0 {
            found.push(format!("{} {open} left without {close}", depth[kind]));
        }
    }
    found
}

#[test]
fn the_demo_figure_and_a_sheet_are_well_formed_pdfs() {
    let figure = demo_figure();
    let sheet = Panels::new()
        .push(&figure, "A")
        .push_captioned(&figure, "B", "the same again");
    for pdf in [figure.to_pdf(), sheet.to_pdf()] {
        assert!(pdf.notes.is_empty(), "{:?}", pdf.notes);
        let found = problems(&pdf.bytes);
        assert!(found.is_empty(), "{found:#?}");
    }
}

#[test]
fn the_checker_catches_a_file_broken_in_each_way_it_checks() {
    let pdf = demo_figure().to_pdf().bytes;
    assert!(problems(&pdf).is_empty());
    let at = |needle: &str| find(&pdf, needle.as_bytes(), 0).unwrap_or_else(|| panic!("{needle}"));
    let breaks: Vec<(&str, Vec<u8>)> = vec![
        // The table off by a byte for the catalog, as a writer that counted
        // a line end it did not write would be.
        ("object 1 is not at", {
            let mut broken = pdf.clone();
            let entry = at("0000000000 65535 f \n") + 20;
            let offset: usize = ascii(&pdf[entry..entry + 10]).parse().unwrap();
            broken.splice(
                entry..entry + 10,
                format!("{:010}", offset + 1).into_bytes(),
            );
            broken
        }),
        // A stream whose length is counted one short. The digit is swapped
        // in place, so every offset after it stays true.
        ("stream 4 is not", {
            let mut broken = pdf.clone();
            let length = at("5 0 obj\n") + "5 0 obj\n".len();
            let end = length + pdf[length..].iter().position(|b| *b == b'\n').unwrap();
            let last = &mut broken[end - 1];
            *last = if *last == b'0' { b'1' } else { *last - 1 };
            broken
        }),
        // A `q` with no `Q`.
        ("left without Q", {
            let mut broken = pdf.clone();
            let content = at("4 M\n");
            broken.splice(content..content + 4, b"q  \n".iter().copied());
            broken
        }),
        // A font the page does not define.
        ("/F9 is used and not defined", {
            let mut broken = pdf.clone();
            let font = at("/F1 ");
            broken[font + 2] = b'9';
            broken
        }),
        ("no %%EOF", pdf[..pdf.len() - 2].to_vec()),
    ];
    for (said, broken) in breaks {
        let found = problems(&broken);
        assert!(
            found.iter().any(|problem| problem.contains(said)),
            "{said}: {found:?}"
        );
    }
}

/// Every figure the project commits, in order.
///
/// `assets/` is left out of the published crate, so a copy from crates.io
/// checks the site's figures alone; a checkout has to have both.
fn committed_figures() -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut figures = Vec::new();
    for (folder, fewest) in [("assets", 40), ("docs/assets/start", 20)] {
        let Ok(entries) = fs::read_dir(root.join(folder)) else {
            assert_eq!(folder, "assets", "{folder} is part of the crate");
            continue;
        };
        let before = figures.len();
        for entry in entries {
            let path = entry.unwrap().path();
            if path.extension().is_some_and(|extension| extension == "svg") {
                figures.push(path);
            }
        }
        assert!(
            figures.len() - before >= fewest,
            "{folder} holds {} figures",
            figures.len() - before
        );
    }
    figures.sort();
    figures
}

/// Every figure the project commits, converted: nothing the converter does
/// not read, and a well-formed file each time.
#[test]
fn every_committed_figure_converts_with_nothing_left_out() {
    for path in committed_figures() {
        let svg = fs::read_to_string(&path).unwrap();
        let pdf = Pdf::from_svg(&svg).unwrap_or_else(|| panic!("{}", path.display()));
        assert!(pdf.notes.is_empty(), "{}: {:?}", path.display(), pdf.notes);
        let found = problems(&pdf.bytes);
        assert!(found.is_empty(), "{}: {found:#?}", path.display());
    }
}

/// Every character beyond ASCII a committed figure sets as text has a width
/// of its own in the measure, rather than the flat six tenths of an em any
/// other is given. karyon writes such characters into legends and into labels
/// it cuts short, and an ellipsis given six tenths, which Helvetica and Arial
/// draw a whole em wide, let a cut label run past its room. The two let
/// through are the examples' own words, not karyon's: São Paulo on a map and
/// a fan of 250° in a title.
#[test]
fn every_character_a_committed_figure_sets_has_a_width_of_its_own() {
    // None of the widths the measure holds is exactly the flat one, so a
    // character measured at it is one the measure does not know.
    let flat = karyon::svg::text_width("\u{4e2d}", 1000.0);
    let mut seen = 0;
    for path in committed_figures() {
        let svg = fs::read_to_string(&path).unwrap();
        for text in svg
            .split("</text>")
            .filter_map(|part| part.rsplit('>').next())
        {
            for c in text.chars().filter(|c| !c.is_ascii()) {
                if matches!(c, '\u{e3}' | '\u{b0}') {
                    continue;
                }
                seen += 1;
                let width = karyon::svg::text_width(&c.to_string(), 1000.0);
                assert!(
                    width != flat,
                    "{} sets {c:?} in {text:?}, which the measure has no width for",
                    path.display()
                );
            }
        }
    }
    // The site's figures alone set 22, omega in the selection legends among
    // them, so a reading that found none would not pass for a clean one.
    assert!(seen >= 20, "{seen}");
}
