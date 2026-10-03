//! PDF, read from the SVG the crate already writes.
//!
//! A figure is drawn once, into an [`SvgWriter`](crate::SvgWriter), and a PDF
//! is that SVG converted: [`Pdf::from_svg`] reads the document back and writes
//! each element as the PDF operators that draw the same thing. There is no
//! second drawing backend, on purpose. Every track draws through a public
//! `&mut SvgWriter`, from hundreds of call sites in the crate and any
//! number outside it, and the writer's own interface already speaks SVG: a
//! track hands [`path`](crate::SvgWriter::path) a `d` string, and gets a fade
//! back from [`fade_down`](crate::SvgWriter::fade_down) as a `url(#id)` to
//! fill with. A backend behind a trait would still have to read path data,
//! colours and references, which is most of what this module does, and it
//! would break every one of those call sites. A [`Panels`](crate::Panels)
//! sheet is assembled from its panels' SVG, too, and comes out here as one
//! page with each panel clipped to its own box, exactly as it is drawn.
//!
//! The price is that the two have to agree on what an SVG can say, and that is
//! held by a test rather than by care: every method of the writer is called
//! once and the result converted, and any element or attribute the converter
//! does not read is a note that fails it by name.
//!
//! # A pixel is three quarters of a point
//!
//! A page is the SVG's size at 0.75 points to the pixel, which is CSS's 96
//! pixels to the inch and the size Inkscape and `rsvg-convert` give the same
//! file, so a figure drawn 900 pixels wide is 675 points, 9.4 inches, either
//! way it reaches a page. Inside, every coordinate is written as the SVG
//! wrote it, in pixels, under one matrix that scales and flips the page, so
//! the two documents can be read side by side. A side over 14,400 points,
//! two hundred inches, is the most a PDF page can be, and a wider figure is
//! drawn smaller with a `/UserUnit` that brings it back to size.
//!
//! # What the base fonts cost
//!
//! Text is set in the fonts every reader has, Helvetica, Courier and Symbol,
//! rather than in the SVG's Inter and JetBrains Mono, and none is embedded. A
//! face is 300 to 400 kilobytes, and embedding one means reading TrueType and
//! cutting a subset of it, which is most of a PDF library; the base fonts cost
//! nothing and the text stays text, to be searched and read aloud. The layout
//! survives the change because [`text_width`](crate::svg::text_width) is never
//! narrower than Helvetica, nor the strong measure than Helvetica-Bold, for
//! printable ASCII and the characters beyond it karyon writes, and each label
//! is anchored with Adobe's own widths. Characters outside what those
//! fonts encode are drawn as question marks and named in [`Pdf::notes`].
//!
//! # What does not carry over
//!
//! The tooltips a group's `<title>` gives in a browser have nowhere to go on
//! a page and are left out. The document's own title and description are
//! kept, as the PDF's title and subject and as the alternative text of the one
//! tagged figure on the page, which is where a screen reader looks for them.
//! The file is not compressed. Most figures come out about the size of their
//! SVG, and one of many dots up to four times it, since a circle is four
//! curves where the SVG writes one element: the association scan the guide
//! shows, 1,748 circles, is 3.7 times its SVG.

mod color;
mod metrics;
mod paint;
mod path;
mod text;
mod xml;

#[cfg(test)]
mod tests;

use std::fs;
use std::io;
use std::path::Path;

use self::paint::{Definitions, Painter};
use self::xml::Token;
use crate::svg::num;

/// A PDF document, converted from an SVG.
///
/// ```
/// use karyon::{plot, Pdf};
///
/// let svg = plot("chr1:1-1000").unwrap().add_coverage(vec![30.0; 1000]).to_svg();
/// let pdf = Pdf::from_svg(&svg).expect("the writer always gives a size");
/// assert!(pdf.bytes.starts_with(b"%PDF-1.4"));
/// assert!(pdf.notes.is_empty());
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pdf {
    /// The document, from `%PDF` to `%%EOF`.
    pub bytes: Vec<u8>,
    /// Each way the PDF differs from the SVG it came from, once each, for a
    /// person to read: a character no base font has, a colour or a transform
    /// that could not be read, an element or an attribute that is not drawn.
    /// Empty for everything this crate draws, unless a label holds a
    /// character outside what the base fonts encode, as a sample name in
    /// Cyrillic or Chinese does.
    pub notes: Vec<String>,
}

impl Pdf {
    /// Converts an SVG document into a one-page PDF.
    ///
    /// The document is one [`SvgWriter`](crate::SvgWriter) wrote, or one
    /// written in the same terms; anything this does not read is named in
    /// [`notes`](Pdf::notes) and drawn as well as it can be, since a figure
    /// with a note is still a figure. A length is read in pixels or in any
    /// absolute unit, so a root `width="100mm"` is a page 100 millimetres
    /// wide with its `viewBox` scaled to fill it; a percentage or an `em` is
    /// named in the notes. `None` only when there is no root `<svg>` to take
    /// a page size from: neither a width and a height it can read nor a
    /// `viewBox`.
    ///
    /// The same SVG always gives the same bytes. Nothing in the file depends
    /// on when or where it was written: there is no creation date and no
    /// document id.
    pub fn from_svg(svg: &str) -> Option<Pdf> {
        let (width, height) = page_size(svg)?;
        let mut notes = Notes::default();
        let definitions = Definitions::read(svg, &mut notes);

        // A reader's smallest page is 3 points a side and its largest 14,400.
        // Past that a page is written at a fraction of its size with a
        // `/UserUnit` saying how much larger a unit is, which every reader
        // since Acrobat 7 honours; one that does not shows the page smaller
        // and in proportion, rather than not at all.
        let points = |pixels: f64| (pixels * 0.75).max(3.0);
        let (page_w, page_h) = (points(width), points(height));
        let unit = (page_w.max(page_h) / 14_400.0).ceil().max(1.0);
        let (page_w, page_h) = (page_w / unit, page_h / unit);

        let mut out: Vec<u8> = Vec::with_capacity(svg.len() + svg.len() / 4 + 4096);
        let version = if unit > 1.0 { "1.6" } else { "1.4" };
        out.extend_from_slice(format!("%PDF-{version}\n").as_bytes());
        // Four bytes past 127, so a tool that sniffs a file treats it as
        // binary and does not rewrite its line ends.
        out.extend_from_slice(b"%\xe2\xe3\xcf\xd3\n");
        let mut offsets: Vec<(usize, usize)> = Vec::new();

        // The content stream goes first and straight into the file, with its
        // length in an object of its own, so the drawing is never held twice.
        offsets.push((4, out.len()));
        out.extend_from_slice(b"4 0 obj\n<< /Length 5 0 R >>\nstream\n");
        let start = out.len();
        // One figure, tagged as one, so a screen reader finds its alt text.
        out.extend_from_slice(b"/Figure << /MCID 0 >> BDC\nq\n");
        let scale = 0.75 / unit;
        path::Matrix([scale, 0.0, 0.0, -scale, 0.0, page_h]).write(&mut out);
        out.extend_from_slice(b"4 M\n");
        let (resources, title, description) = {
            let mut painter = Painter::new(&mut out, &definitions, &mut notes);
            painter.run(svg);
            (painter.resources, painter.title, painter.description)
        };
        out.extend_from_slice(b"Q\nEMC");
        let length = out.len() - start;
        // The line end before `endstream` is the file's, not the stream's.
        out.extend_from_slice(b"\nendstream\nendobj\n");
        offsets.push((5, out.len()));
        out.extend_from_slice(format!("5 0 obj\n{length}\nendobj\n").as_bytes());

        // Then everything that names an object by number, now the numbers
        // are known.
        let first_font = 7;
        let first_alpha = first_font + resources.fonts.len();
        let first_fade = first_alpha + resources.alphas.len();
        let info = first_fade + 2 * resources.fades.len();
        let (tree, figure, parents) = (info + 1, info + 2, info + 3);
        let translucent = !resources.fades.is_empty()
            || resources
                .alphas
                .iter()
                .any(|(fill, stroke)| *fill < 1.0 || *stroke < 1.0);

        object(
            &mut out,
            &mut offsets,
            1,
            &format!(
                "<< /Type /Catalog /Pages 2 0 R /MarkInfo << /Marked true >> \
                 /StructTreeRoot {tree} 0 R /Lang (en) \
                 /ViewerPreferences << /DisplayDocTitle true >> >>"
            ),
        );
        object(
            &mut out,
            &mut offsets,
            2,
            "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        );
        let mut page = format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {} {}] /Resources 6 0 R \
             /Contents 4 0 R /StructParents 0",
            num(page_w),
            num(page_h)
        );
        if unit > 1.0 {
            page.push_str(&format!(" /UserUnit {}", num(unit)));
        }
        if translucent {
            // Blend in RGB whatever the reader's device is, so a translucent
            // mark comes out the colour it is in the SVG.
            page.push_str(" /Group << /S /Transparency /CS /DeviceRGB >>");
        }
        page.push_str(" >>");
        object(&mut out, &mut offsets, 3, &page);

        let mut dictionary = String::from("<<");
        let named = |prefix: &str, first: usize, count: usize, step: usize| -> String {
            (0..count)
                .map(|at| format!(" /{prefix}{} {} 0 R", at + 1, first + at * step))
                .collect()
        };
        if !resources.fonts.is_empty() {
            dictionary.push_str(" /Font <<");
            dictionary.push_str(&named("F", first_font, resources.fonts.len(), 1));
            dictionary.push_str(" >>");
        }
        if !resources.alphas.is_empty() {
            dictionary.push_str(" /ExtGState <<");
            dictionary.push_str(&named("G", first_alpha, resources.alphas.len(), 1));
            dictionary.push_str(" >>");
        }
        if !resources.fades.is_empty() {
            dictionary.push_str(" /XObject <<");
            dictionary.push_str(&named("I", first_fade, resources.fades.len(), 2));
            dictionary.push_str(" >>");
        }
        dictionary.push_str(" >>");
        object(&mut out, &mut offsets, 6, &dictionary);

        for (at, face) in resources.fonts.iter().enumerate() {
            // Symbol carries its own encoding; the Latin faces are told to
            // read WinAnsi, which is what their strings are written in.
            let encoding = if *face == text::Face::Symbol {
                ""
            } else {
                " /Encoding /WinAnsiEncoding"
            };
            object(
                &mut out,
                &mut offsets,
                first_font + at,
                &format!(
                    "<< /Type /Font /Subtype /Type1 /BaseFont /{}{encoding} >>",
                    face.base_font()
                ),
            );
        }
        for (at, (fill, stroke)) in resources.alphas.iter().enumerate() {
            object(
                &mut out,
                &mut offsets,
                first_alpha + at,
                &format!(
                    "<< /Type /ExtGState /ca {} /CA {} >>",
                    num(*fill),
                    num(*stroke)
                ),
            );
        }
        for (at, fade) in resources.fades.iter().enumerate() {
            let (image, mask) = (first_fade + 2 * at, first_fade + 2 * at + 1);
            let rows = fade.alpha.len();
            // Smoothed by the reader between rows, so the fade has no steps
            // in a reader that honours /Interpolate and none it could see in
            // one that does not.
            offsets.push((image, out.len()));
            out.extend_from_slice(
                format!(
                    "{image} 0 obj\n<< /Type /XObject /Subtype /Image /Width 1 /Height {rows} \
                     /ColorSpace /DeviceRGB /BitsPerComponent 8 /Interpolate true \
                     /SMask {mask} 0 R /Length {} >>\nstream\n",
                    fade.rgb.len()
                )
                .as_bytes(),
            );
            out.extend_from_slice(&fade.rgb);
            out.extend_from_slice(b"\nendstream\nendobj\n");
            offsets.push((mask, out.len()));
            out.extend_from_slice(
                format!(
                    "{mask} 0 obj\n<< /Type /XObject /Subtype /Image /Width 1 /Height {rows} \
                     /ColorSpace /DeviceGray /BitsPerComponent 8 /Interpolate true \
                     /Length {} >>\nstream\n",
                    fade.alpha.len()
                )
                .as_bytes(),
            );
            out.extend_from_slice(&fade.alpha);
            out.extend_from_slice(b"\nendstream\nendobj\n");
        }

        let title = title.filter(|text| !text.is_empty());
        let description = description.filter(|text| !text.is_empty());
        let mut properties = String::from("<<");
        if let Some(title) = &title {
            properties.push_str(&format!(" /Title {}", utf16(title)));
        }
        if let Some(description) = &description {
            properties.push_str(&format!(" /Subject {}", utf16(description)));
        }
        properties.push_str(&format!(" /Producer (karyon {}) >>", crate::VERSION));
        object(&mut out, &mut offsets, info, &properties);
        object(
            &mut out,
            &mut offsets,
            tree,
            &format!(
                "<< /Type /StructTreeRoot /K {figure} 0 R /ParentTree {parents} 0 R \
                 /ParentTreeNextKey 1 >>"
            ),
        );
        // What a screen reader says for the figure is what the SVG's
        // `aria-labelledby` points it at: the title, then the description.
        let alt = match (&title, &description) {
            (Some(title), Some(description)) => {
                let stop = if title.ends_with(['.', '!', '?']) {
                    ""
                } else {
                    "."
                };
                Some(format!("{title}{stop} {description}"))
            }
            (Some(one), None) | (None, Some(one)) => Some(one.clone()),
            (None, None) => None,
        };
        let alt = alt.map_or_else(String::new, |alt| format!(" /Alt {}", utf16(&alt)));
        object(
            &mut out,
            &mut offsets,
            figure,
            &format!("<< /Type /StructElem /S /Figure /P {tree} 0 R /Pg 3 0 R /K 0{alt} >>"),
        );
        object(
            &mut out,
            &mut offsets,
            parents,
            &format!("<< /Nums [0 [{figure} 0 R]] >>"),
        );

        // Cross-reference: each entry exactly twenty bytes, ten digits of
        // offset, five of generation, the kind, and a two-byte line end.
        offsets.sort_unstable();
        debug_assert!(offsets.iter().enumerate().all(|(at, (n, _))| *n == at + 1));
        let xref = out.len();
        let count = parents + 1;
        out.extend_from_slice(format!("xref\n0 {count}\n0000000000 65535 f \n").as_bytes());
        for (_, offset) in &offsets {
            out.extend_from_slice(format!("{offset:010} 00000 n \n").as_bytes());
        }
        out.extend_from_slice(
            format!(
                "trailer\n<< /Size {count} /Root 1 0 R /Info {info} 0 R >>\nstartxref\n{xref}\n%%EOF\n"
            )
            .as_bytes(),
        );
        Some(Pdf {
            bytes: out,
            notes: notes.0,
        })
    }

    /// Writes the document to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, &self.bytes)
    }
}

/// The PDF of an SVG this crate drew.
///
/// The one place the drawing types turn their SVG into a PDF, so the reason
/// it cannot fail is written once: [`SvgWriter::finish`](crate::SvgWriter::finish)
/// opens every document with a root `<svg>` carrying a width and a height,
/// written through [`num`], which writes a number for anything, nought for
/// what is not one. A test in this module draws every kind of document the
/// crate makes, empty ones included, and converts it.
pub(crate) fn drawn(svg: &str) -> Pdf {
    Pdf::from_svg(svg).expect("SvgWriter::finish always writes a sized root")
}

/// Whether `path` names a PDF: its extension is `pdf` in any case.
///
/// ```
/// use karyon::pdf::named_pdf;
///
/// assert!(named_pdf("out/rpoB.pdf".as_ref()));
/// assert!(named_pdf("RPOB.PDF".as_ref()));
/// assert!(!named_pdf("rpoB.svg".as_ref()));
/// assert!(!named_pdf("pdf".as_ref()));
/// ```
pub fn named_pdf(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("pdf"))
}

/// Writes object `number` with `body`, and where it starts.
fn object(out: &mut Vec<u8>, offsets: &mut Vec<(usize, usize)>, number: usize, body: &str) {
    offsets.push((number, out.len()));
    out.extend_from_slice(format!("{number} 0 obj\n{body}\nendobj\n").as_bytes());
}

/// What the converter has to say, each thing once, in the order it came up.
#[derive(Debug, Default)]
pub(crate) struct Notes(Vec<String>);

impl Notes {
    pub(crate) fn add(&mut self, note: String) {
        if !self.0.contains(&note) {
            self.0.push(note);
        }
    }
}

/// The size of the root `<svg>`, in pixels, read as the painter reads it.
///
/// A width of nought, or one that came out of the writer as nought because it
/// was not a number, is still a size: the page is made as small as a reader
/// allows, rather than the conversion failing over a figure with nothing in it.
fn page_size(svg: &str) -> Option<(f64, f64)> {
    let root = xml::tokens(svg).find_map(|token| match token {
        Token::Start {
            name, attributes, ..
        } => Some((name, attributes)),
        _ => None,
    })?;
    let (name, attributes) = root;
    if name != "svg" {
        return None;
    }
    // Whatever the root holds that is not read is named when the painter
    // reaches it, so it is named once.
    paint::viewport(attributes, &mut Notes::default())
}

/// Text as a PDF string a reader shows in any script: UTF-16, big-endian,
/// after a byte order mark, in hexadecimal.
fn utf16(text: &str) -> String {
    let mut hex = String::from("<FEFF");
    for unit in text.encode_utf16() {
        hex.push_str(&format!("{unit:04X}"));
    }
    hex.push('>');
    hex
}
