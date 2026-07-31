//! A small SVG writer: exactly the elements the tracks need, and no more.
//!
//! Output is plain SVG 1.1 with no scripts, no external references and no
//! embedded fonts, so it opens unchanged in a browser, in Inkscape and in
//! Illustrator, and it survives being pasted into a manuscript figure.

use std::fmt::Write as _;

/// Horizontal anchoring of a text label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    /// `x` is the left edge of the text.
    Start,
    /// `x` is the centre of the text.
    Middle,
    /// `x` is the right edge of the text.
    End,
}

impl Anchor {
    fn as_str(self) -> &'static str {
        match self {
            Anchor::Start => "start",
            Anchor::Middle => "middle",
            Anchor::End => "end",
        }
    }
}

/// How a run of text is painted, bundled so the text calls stay narrow.
struct Ink<'a> {
    fill: &'a str,
    size: f64,
    anchor: Anchor,
}

/// Accumulates SVG elements and closes them into a finished document.
///
/// Tracks receive this through [`DrawContext`](crate::track::DrawContext) and
/// never build SVG strings themselves, which keeps escaping and number
/// formatting in one place.
#[derive(Debug, Default)]
pub struct SvgWriter {
    body: String,
    defs: String,
    next_id: usize,
    open_groups: usize,
}

impl SvgWriter {
    /// An empty document.
    pub fn new() -> Self {
        SvgWriter::default()
    }

    /// A filled rectangle.
    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        self.rect_opacity(x, y, w, h, fill, 1.0);
    }

    /// A filled rectangle with a fill opacity between 0 and 1.
    pub fn rect_opacity(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str, opacity: f64) {
        if w <= 0.0 || h <= 0.0 || !finite(&[x, y, w, h]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<rect x="{}" y="{}" width="{}" height="{}" fill="{}""#,
            num(x),
            num(y),
            num(w),
            num(h),
            fill
        );
        self.push_opacity("fill-opacity", opacity);
        self.body.push_str("/>");
    }

    /// A straight line.
    pub fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
        if !finite(&[x1, y1, x2, y2]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}"/>"#,
            num(x1),
            num(y1),
            num(x2),
            num(y2),
            stroke,
            num(width)
        );
    }

    /// A filled circle.
    pub fn circle(&mut self, cx: f64, cy: f64, r: f64, fill: &str) {
        if r <= 0.0 || !finite(&[cx, cy, r]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<circle cx="{}" cy="{}" r="{}" fill="{}"/>"#,
            num(cx),
            num(cy),
            num(r),
            fill
        );
    }

    /// A filled polygon through `points`.
    pub fn polygon(&mut self, points: &[(f64, f64)], fill: &str) {
        if points.len() < 3 {
            return;
        }
        let Some(list) = point_list(points) else {
            return;
        };
        let _ = write!(self.body, r#"<polygon points="{list}" fill="{fill}"/>"#);
    }

    /// An open stroked polyline through `points`.
    pub fn polyline(&mut self, points: &[(f64, f64)], stroke: &str, width: f64) {
        if points.len() < 2 {
            return;
        }
        let Some(list) = point_list(points) else {
            return;
        };
        let _ = write!(
            self.body,
            r#"<polyline points="{}" fill="none" stroke="{}" stroke-width="{}" stroke-linejoin="round"/>"#,
            list,
            stroke,
            num(width)
        );
    }

    /// A filled path from a ready-made `d` attribute.
    pub fn path(&mut self, d: &str, fill: &str, opacity: f64) {
        if d.is_empty() {
            return;
        }
        let _ = write!(self.body, r#"<path d="{d}" fill="{fill}""#);
        self.push_opacity("fill-opacity", opacity);
        self.body.push_str("/>");
    }

    /// A text label. `y` is the text baseline, not its centre or its top.
    pub fn text(&mut self, x: f64, y: f64, content: &str, fill: &str, size: f64, anchor: Anchor) {
        self.write_text(x, y, content, Ink { fill, size, anchor }, false);
    }

    /// A bold text label, for titles and for letters that must read at speed.
    pub fn text_bold(
        &mut self,
        x: f64,
        y: f64,
        content: &str,
        fill: &str,
        size: f64,
        anchor: Anchor,
    ) {
        self.write_text(x, y, content, Ink { fill, size, anchor }, true);
    }

    fn write_text(&mut self, x: f64, y: f64, content: &str, ink: Ink<'_>, bold: bool) {
        if content.is_empty() || !finite(&[x, y]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<text x="{}" y="{}" fill="{}" font-size="{}" text-anchor="{}""#,
            num(x),
            num(y),
            ink.fill,
            num(ink.size),
            ink.anchor.as_str()
        );
        if bold {
            self.body.push_str(r#" font-weight="bold""#);
        }
        let _ = write!(self.body, ">{}</text>", escape(content));
    }

    /// Opens a group clipped to a rectangle. Pair it with [`SvgWriter::end_group`].
    pub fn begin_clip(&mut self, x: f64, y: f64, w: f64, h: f64) {
        let id = format!("karyon-clip-{}", self.next_id);
        self.next_id += 1;
        let _ = write!(
            self.defs,
            r#"<clipPath id="{}"><rect x="{}" y="{}" width="{}" height="{}"/></clipPath>"#,
            id,
            num(x),
            num(y),
            num(w.max(0.0)),
            num(h.max(0.0))
        );
        let _ = write!(self.body, r#"<g clip-path="url(#{id})">"#);
        self.open_groups += 1;
    }

    /// Closes the innermost open group.
    pub fn end_group(&mut self) {
        if self.open_groups == 0 {
            return;
        }
        self.open_groups -= 1;
        self.body.push_str("</g>");
    }

    /// Wraps everything written so far into a standalone SVG document.
    pub fn finish(
        mut self,
        width: f64,
        height: f64,
        background: &str,
        font_family: &str,
    ) -> String {
        while self.open_groups > 0 {
            self.end_group();
        }
        let mut out = String::with_capacity(self.body.len() + self.defs.len() + 512);
        let _ = write!(
            out,
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}" font-family="{}">"#,
            num(width),
            num(height),
            num(width),
            num(height),
            escape(font_family)
        );
        if !self.defs.is_empty() {
            let _ = write!(out, "<defs>{}</defs>", self.defs);
        }
        if background != "none" {
            let _ = write!(
                out,
                r#"<rect x="0" y="0" width="{}" height="{}" fill="{}"/>"#,
                num(width),
                num(height),
                background
            );
        }
        out.push_str(&self.body);
        out.push_str("</svg>");
        out
    }

    fn push_opacity(&mut self, attr: &str, opacity: f64) {
        if opacity < 1.0 && opacity.is_finite() {
            let _ = write!(self.body, r#" {}="{}""#, attr, num(opacity.max(0.0)));
        }
    }
}

/// Formats a coordinate compactly: at most three decimals, no trailing zeros.
///
/// SVG files of a genome-wide figure are mostly digits, so this is the single
/// biggest lever on output size.
pub fn num(v: f64) -> String {
    if !v.is_finite() {
        return "0".to_string();
    }
    let rounded = (v * 1000.0).round() / 1000.0;
    if rounded == rounded.trunc() && rounded.abs() < 1e15 {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded}")
    }
}

/// Escapes the five XML characters that would otherwise break the document.
///
/// Sequence names and feature names come from user files and routinely contain
/// `&` and `<`, so nothing reaches the output without passing through here.
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

/// Rough advance width of a string, used to decide whether a label fits.
///
/// Deliberately an estimate: measuring properly would mean parsing font
/// metrics, and every renderer would still disagree slightly.
pub fn text_width(text: &str, font_size: f64) -> f64 {
    text.chars().count() as f64 * font_size * 0.55
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}

fn point_list(points: &[(f64, f64)]) -> Option<String> {
    let mut list = String::with_capacity(points.len() * 12);
    for (x, y) in points {
        if !x.is_finite() || !y.is_finite() {
            return None;
        }
        if !list.is_empty() {
            list.push(' ');
        }
        let _ = write!(list, "{},{}", num(*x), num(*y));
    }
    Some(list)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_lose_their_trailing_zeros() {
        assert_eq!(num(10.0), "10");
        assert_eq!(num(10.5), "10.5");
        assert_eq!(num(1.0 / 3.0), "0.333");
        assert_eq!(num(-0.0001), "0");
    }

    #[test]
    fn non_finite_numbers_never_reach_the_output() {
        assert_eq!(num(f64::NAN), "0");
        assert_eq!(num(f64::INFINITY), "0");
    }

    #[test]
    fn escaping_covers_every_xml_metacharacter() {
        assert_eq!(escape(r#"a&b<c>d"e'f"#), "a&amp;b&lt;c&gt;d&quot;e&apos;f");
    }

    #[test]
    fn text_content_is_escaped() {
        let mut svg = SvgWriter::new();
        svg.text(0.0, 0.0, "gene<A&B>", "#000", 10.0, Anchor::Start);
        let out = svg.finish(10.0, 10.0, "#fff", "sans-serif");
        assert!(out.contains("gene&lt;A&amp;B&gt;"));
        assert!(!out.contains("gene<A"));
    }

    #[test]
    fn degenerate_shapes_are_skipped() {
        let mut svg = SvgWriter::new();
        svg.rect(0.0, 0.0, 0.0, 10.0, "#000");
        svg.rect(0.0, 0.0, f64::NAN, 10.0, "#000");
        svg.circle(0.0, 0.0, -1.0, "#000");
        svg.polygon(&[(0.0, 0.0), (1.0, 1.0)], "#000");
        svg.polyline(&[(0.0, f64::NAN), (1.0, 1.0)], "#000", 1.0);
        let out = svg.finish(10.0, 10.0, "none", "sans-serif");
        assert!(!out.contains("<rect"), "{out}");
        assert!(!out.contains("<circle"));
        assert!(!out.contains("<polygon"));
        assert!(!out.contains("<polyline"));
    }

    #[test]
    fn unclosed_groups_are_closed_on_finish() {
        let mut svg = SvgWriter::new();
        svg.begin_clip(0.0, 0.0, 10.0, 10.0);
        svg.begin_clip(1.0, 1.0, 5.0, 5.0);
        let out = svg.finish(10.0, 10.0, "#fff", "sans-serif");
        assert_eq!(out.matches("<g ").count(), 2);
        assert_eq!(out.matches("</g>").count(), 2);
        assert!(out.ends_with("</svg>"));
    }

    #[test]
    fn transparent_background_writes_no_page_rect() {
        let svg = SvgWriter::new();
        let out = svg.finish(10.0, 10.0, "none", "sans-serif");
        assert!(!out.contains("<rect"));
    }
}
