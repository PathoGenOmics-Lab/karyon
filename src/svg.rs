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

    /// A filled rectangle with rounded corners.
    ///
    /// The radius is clamped to half the smaller side, so a short bar becomes a
    /// lozenge rather than losing its geometry to an oversized corner.
    pub fn rect_rounded(&mut self, x: f64, y: f64, w: f64, h: f64, radius: f64, fill: &str) {
        if w <= 0.0 || h <= 0.0 || !finite(&[x, y, w, h, radius]) {
            return;
        }
        let radius = radius.min(w / 2.0).min(h / 2.0).max(0.0);
        if radius <= 0.05 {
            return self.rect(x, y, w, h, fill);
        }
        let _ = write!(
            self.body,
            r#"<rect x="{}" y="{}" width="{}" height="{}" rx="{}" fill="{}"/>"#,
            num(x),
            num(y),
            num(w),
            num(h),
            num(radius),
            fill
        );
    }

    /// A filled circle with a ring in the surface colour around it.
    ///
    /// The ring is what keeps overlapping markers legible: two dots on top of
    /// each other read as two dots rather than one blob, and a dot sitting on
    /// its own stem stays a dot. It is spacing, not decoration, which is why it
    /// wears the page colour rather than a darker shade of the mark.
    pub fn circle_ringed(&mut self, cx: f64, cy: f64, r: f64, fill: &str, ring: &str, width: f64) {
        if r <= 0.0 || !finite(&[cx, cy, r]) {
            return;
        }
        if width > 0.0 && ring != "none" {
            self.circle(cx, cy, r + width, ring);
        }
        self.circle(cx, cy, r, fill);
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

    /// An outlined rectangle with no fill.
    pub fn rect_outline(&mut self, x: f64, y: f64, w: f64, h: f64, stroke: &str, width: f64) {
        if w <= 0.0 || h <= 0.0 || !finite(&[x, y, w, h]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="{}" stroke-width="{}"/>"#,
            num(x),
            num(y),
            num(w),
            num(h),
            stroke,
            num(width)
        );
    }

    /// An outlined path with no fill, for a shape that contains other shapes.
    pub fn path_stroked(&mut self, d: &str, stroke: &str, width: f64) {
        if d.is_empty() {
            return;
        }
        let _ = write!(
            self.body,
            r#"<path d="{}" fill="none" stroke="{}" stroke-width="{}"/>"#,
            d,
            stroke,
            num(width)
        );
    }

    /// Opens a group clipped to an arbitrary path, for shapes a rectangle
    /// cannot describe. Pair it with [`SvgWriter::end_group`].
    pub fn begin_clip_path(&mut self, d: &str) {
        let id = format!("karyon-clip-{}", self.next_id);
        self.next_id += 1;
        let _ = write!(
            self.defs,
            r#"<clipPath id="{id}"><path d="{d}"/></clipPath>"#
        );
        let _ = write!(self.body, r#"<g clip-path="url(#{id})">"#);
        self.open_groups += 1;
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

    /// A single glyph stretched to fill an exact box.
    ///
    /// This is what a sequence logo is made of, and it is the one place where
    /// text has to obey a geometry rather than a font size. `width` is enforced
    /// with `textLength`, so the renderer does the horizontal fitting and the
    /// result does not depend on this crate guessing font metrics. Height comes
    /// from the caller choosing `font_size` against the cap height ratio of the
    /// theme, and `baseline` is the bottom of the box.
    ///
    /// Symbols with descenders (`g`, `y`, `p`) hang below the baseline by
    /// design; a logo of uppercase letters sits exactly in its box.
    pub fn glyph(
        &mut self,
        x: f64,
        baseline: f64,
        width: f64,
        font_size: f64,
        content: &str,
        fill: &str,
    ) {
        if content.is_empty() || width <= 0.0 || font_size <= 0.0 || !finite(&[x, baseline]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<text x="{}" y="{}" font-size="{}" fill="{}" textLength="{}" lengthAdjust="spacingAndGlyphs" font-weight="bold">{}</text>"#,
            num(x),
            num(baseline),
            num(font_size),
            fill,
            num(width),
            escape(content)
        );
    }

    /// A text label turned about its own anchor point.
    ///
    /// `degrees` runs clockwise, so -90 stands a label on end reading upwards,
    /// which is how a column gets a caption wider than the column.
    pub fn text_rotated(
        &mut self,
        at: (f64, f64),
        degrees: f64,
        content: &str,
        fill: &str,
        size: f64,
        anchor: Anchor,
    ) {
        let (x, y) = at;
        if content.is_empty() || !finite(&[x, y, degrees]) {
            return;
        }
        let _ = write!(
            self.body,
            r#"<text x="0" y="0" transform="translate({} {}) rotate({})" fill="{}" font-size="{}" text-anchor="{}">{}</text>"#,
            num(x),
            num(y),
            num(degrees),
            fill,
            num(size),
            anchor.as_str(),
            escape(content)
        );
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
    fn a_corner_radius_cannot_outgrow_its_rectangle() {
        let mut svg = SvgWriter::new();
        svg.rect_rounded(0.0, 0.0, 10.0, 4.0, 20.0, "#000");
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        // Half the shorter side, not the twenty that was asked for.
        assert!(out.contains(r#"rx="2""#), "{out}");
    }

    #[test]
    fn a_zero_radius_falls_back_to_a_plain_rectangle() {
        let mut svg = SvgWriter::new();
        svg.rect_rounded(0.0, 0.0, 10.0, 4.0, 0.0, "#000");
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        assert!(out.contains("<rect"));
        assert!(!out.contains("rx="));
    }

    #[test]
    fn a_ringed_marker_draws_the_ring_underneath() {
        let mut svg = SvgWriter::new();
        svg.circle_ringed(10.0, 10.0, 4.0, "#111111", "#ffffff", 2.0);
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        let ring = out.find("#ffffff").unwrap();
        let fill = out.find("#111111").unwrap();
        assert!(ring < fill, "the ring has to be painted first");
        assert!(out.contains(r#"r="6""#), "ring radius is r + width");
    }

    #[test]
    fn a_ring_can_be_turned_off() {
        let mut svg = SvgWriter::new();
        svg.circle_ringed(10.0, 10.0, 4.0, "#111111", "none", 2.0);
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        assert_eq!(out.matches("<circle").count(), 1);
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
