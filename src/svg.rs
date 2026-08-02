//! A small SVG writer: exactly the elements the tracks need, and no more.
//!
//! Output is plain SVG 1.1 with no scripts, no external references and no
//! embedded fonts, so it opens unchanged in a browser, in Inkscape and in
//! Illustrator, and it survives being pasted into a manuscript figure.
//!
//! Writing it by hand is what keeps the dependency list empty. A general SVG
//! library brings a document model, a namespace machinery and a release
//! schedule of its own, in exchange for elements no track here draws. What a
//! plotting backend is asked for is narrower than that and sharper in two
//! places: it has to measure its own text, and it has to be frugal with the
//! coordinates it writes, which are the bulk of a large figure.
//!
//! # A layout has to know how wide its text is
//!
//! Where a track label sits, whether a feature name fits inside its box and how
//! much gutter to reserve are all settled before a single element is written,
//! and there is no font engine here to ask. [`text_width`] answers from
//! Helvetica's own advance widths, Helvetica being the first font in
//! [`Theme::font_family`](crate::Theme::font_family) and metrically the same as
//! Arial, so under the default theme the answer is exact rather than a guess.
//! Nearly every track calls it, because a label that overruns the room reserved
//! for it is a label that gets clipped.
//!
//! # Nothing malformed leaves here
//!
//! Names come out of user files, so every string is written through [`escape`]
//! and every number through [`num`]: an `&` in a sequence name cannot break the
//! document, and a coordinate that came out `NaN` cannot reach the output as
//! one. A shape with a non-finite corner or a width of zero is dropped rather
//! than written, and [`SvgWriter::finish`] closes whatever groups a track left
//! open. None of that is a judgement call at the call site, which is the point:
//! there is one way out of the crate and it is this one.

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
    id_prefix: String,
    next_id: usize,
    open_groups: usize,
    /// What the document as a whole is, for `<title>`.
    name: String,
    /// What it shows, for `<desc>`.
    description: String,
}

impl SvgWriter {
    /// An empty document.
    pub fn new() -> Self {
        SvgWriter::default()
    }

    /// An empty document whose generated ids all begin with `prefix`.
    ///
    /// An id in SVG belongs to the whole document rather than to the element
    /// carrying it, and `url(#id)` resolves to the first match anywhere in it.
    /// Two documents nested into one sheet would therefore both claim
    /// `karyon-clip-0`, and the second one's clips would silently resolve to
    /// the first one's rectangle, cropping it to the wrong band.
    /// [`Panels`](crate::Panels) hands every figure its own prefix; a document
    /// standing on its own needs none.
    pub fn with_id_prefix(prefix: impl Into<String>) -> Self {
        SvgWriter {
            id_prefix: prefix.into(),
            ..SvgWriter::default()
        }
    }

    /// The next unused id, carrying whatever prefix this document was given.
    fn next_clip_id(&mut self) -> String {
        let id = format!("{}karyon-clip-{}", self.id_prefix, self.next_id);
        self.next_id += 1;
        id
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
        let id = self.next_clip_id();
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
        // A track that derives its label size from the theme can arrive here
        // with a size at or below zero, and `font-size` is a length SVG does
        // not allow to be negative. Dropping matches `glyph`, and the element
        // would have drawn nothing anyway.
        if content.is_empty() || ink.size <= 0.0 || !finite(&[x, y, ink.size]) {
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
        // Same reason as `write_text`: `font-size` is a length that cannot be
        // negative, and this path writes its own.
        if content.is_empty() || size <= 0.0 || !finite(&[x, y, degrees, size]) {
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
        let id = self.next_clip_id();
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

    /// Opens a group carrying a tooltip, closed by [`SvgWriter::end_group`].
    ///
    /// A `<title>` as the first child of a group is what a browser shows when a
    /// pointer rests anywhere inside it, and what a screen reader announces.
    /// That is SVG 1.1 and not an extension: it needs no script, nothing is
    /// fetched to render it, and Inkscape and Illustrator both keep it, so a
    /// figure that gains one has given up none of what it promised.
    ///
    /// A group rather than a `<title>` inside the shape, because most glyphs
    /// here are more than one shape. A gene is an arrow and its name, a variant
    /// is a stem and a head, and a tooltip on half of one is worse than none.
    ///
    /// An empty `text` opens a plain group, so a caller with nothing to say
    /// stays balanced against its own [`SvgWriter::end_group`] without paying
    /// for an empty element.
    ///
    /// ```
    /// use karyon::SvgWriter;
    ///
    /// let mut svg = SvgWriter::new();
    /// svg.begin_titled("rpoB, 759,807 to 763,325, forward");
    /// svg.rect(10.0, 10.0, 100.0, 12.0, "#0072b2");
    /// svg.end_group();
    ///
    /// let out = svg.finish(200.0, 40.0, "#ffffff", "Helvetica");
    /// assert!(out.contains("<title>rpoB, 759,807 to 763,325, forward</title>"));
    /// ```
    pub fn begin_titled(&mut self, text: &str) {
        if text.is_empty() {
            self.body.push_str("<g>");
        } else {
            let _ = write!(self.body, "<g><title>{}</title>", escape(text));
        }
        self.open_groups += 1;
    }

    /// The same group, but transparent to the pointer.
    ///
    /// A `<title>` resolves to the innermost group under the pointer, so a
    /// translucent shape drawn over data takes every hover inside its own
    /// footprint and the data underneath can never be reached. That is the
    /// wrong way round when the shape on top is decoration: the region marker
    /// on an ideogram says where the figure is looking, and the bands beneath
    /// it are what a reader is pointing at.
    ///
    /// `pointer-events="none"` hands the hover back to the shapes below while
    /// leaving the title where a screen reader still finds it, which is why
    /// this is a group attribute rather than a dropped title.
    ///
    /// ```
    /// use karyon::SvgWriter;
    ///
    /// let mut svg = SvgWriter::new();
    /// svg.begin_titled_inert("region shown, 1 to 1,000");
    /// svg.rect(10.0, 10.0, 100.0, 12.0, "#d7263d");
    /// svg.end_group();
    ///
    /// let out = svg.finish(200.0, 40.0, "#ffffff", "Helvetica");
    /// assert!(out.contains(r#"<g pointer-events="none">"#));
    /// ```
    pub fn begin_titled_inert(&mut self, text: &str) {
        self.body.push_str(r#"<g pointer-events="none">"#);
        if !text.is_empty() {
            let _ = write!(self.body, "<title>{}</title>", escape(text));
        }
        self.open_groups += 1;
    }

    /// Names the document as a whole, for its `<title>` and `<desc>`.
    ///
    /// These are the first two children of the root element, which is where a
    /// screen reader looks and what a browser shows as the tooltip of the
    /// figure itself. Passing an empty string for either leaves that element
    /// out rather than writing an empty one.
    pub fn describe(&mut self, name: &str, description: &str) {
        self.name = name.to_string();
        self.description = description.to_string();
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
        // `role="img"` and `aria-labelledby` are what make an assistive
        // technology read the title and the description instead of walking
        // several thousand rectangles that mean nothing one at a time.
        let title_id = format!("{}karyon-title", self.id_prefix);
        let desc_id = format!("{}karyon-desc", self.id_prefix);
        let mut labels = String::new();
        if !self.name.is_empty() {
            labels.push_str(&title_id);
        }
        if !self.description.is_empty() {
            if !labels.is_empty() {
                labels.push(' ');
            }
            labels.push_str(&desc_id);
        }
        let accessible = if labels.is_empty() {
            String::new()
        } else {
            format!(r#" role="img" aria-labelledby="{labels}""#)
        };
        let _ = write!(
            out,
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}" font-family="{}"{}>"#,
            num(width),
            num(height),
            num(width),
            num(height),
            escape(font_family),
            accessible
        );
        // First two children, which is where the specification puts them.
        if !self.name.is_empty() {
            let _ = write!(
                out,
                r#"<title id="{}">{}</title>"#,
                title_id,
                escape(&self.name)
            );
        }
        if !self.description.is_empty() {
            let _ = write!(
                out,
                r#"<desc id="{}">{}</desc>"#,
                desc_id,
                escape(&self.description)
            );
        }
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

    /// The definitions and the elements written so far, without the wrapper.
    ///
    /// [`Panels`](crate::Panels) has to interleave text of its own with
    /// figures it did not draw, so it takes the pieces and assembles the
    /// document itself rather than splicing into a finished one.
    pub(crate) fn into_parts(mut self) -> (String, String) {
        while self.open_groups > 0 {
            self.end_group();
        }
        (self.defs, self.body)
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
    // Scaling by a thousand overflows to infinity above `f64::MAX / 1000`, and
    // `format!` would then write the literal `inf` into the document. A finite
    // value that large has an ulp far wider than one, so it is already an
    // integer and rounding it to three decimals is the identity.
    let scaled = v * 1000.0;
    let rounded = if scaled.is_finite() {
        scaled.round() / 1000.0
    } else {
        v
    };
    if rounded == rounded.trunc() && rounded.abs() < 1e15 {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded}")
    }
}

/// Escapes the five XML metacharacters, and drops the characters XML has no
/// way to write at all.
///
/// Sequence names and feature names come from user files and routinely contain
/// `&` and `<`, so nothing reaches the output without passing through here.
///
/// A name that carries a stray control byte is a harder case than an `&`: XML
/// 1.0 §2.2 admits only tab, newline and carriage return below `U+0020`, and a
/// numeric character reference cannot encode the others either. Such a
/// character is dropped, because the alternative is a document no parser will
/// open over a byte that would have drawn nothing.
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            // Outside XML 1.0's `Char` production, so not writable at all.
            // Surrogates need no arm here: a Rust `str` cannot hold one.
            c if (c < '\u{20}' && c != '\t' && c != '\n' && c != '\r')
                || c == '\u{fffe}'
                || c == '\u{ffff}' => {}
            _ => out.push(c),
        }
    }
    out
}

/// Advance width of a string, used to reserve room and to decide whether a
/// label fits where it is going.
///
/// The widths are Helvetica's own, which is the first font in
/// [`Theme::font_family`](crate::Theme::font_family) and metrically the same as
/// Arial, so for the default theme this is exact rather than approximate. That
/// matters more than it sounds: one flat width per character under-reserves for
/// a run of capitals by about a fifth, which is precisely what a column of
/// sample accessions is, and a label that overruns the space reserved for it
/// gets clipped.
///
/// A character outside printable ASCII falls back to a wide default, so an
/// accented name reserves a little too much rather than too little. Another
/// font stack will disagree in the third significant figure.
pub fn text_width(text: &str, font_size: f64) -> f64 {
    let per_mille: f64 = text
        .chars()
        .map(|c| {
            let index = c as u32;
            if (32..127).contains(&index) {
                HELVETICA_WIDTHS[index as usize - 32] as f64
            } else {
                600.0
            }
        })
        .sum();
    per_mille / 1000.0 * font_size
}

/// Helvetica advance widths for printable ASCII, in thousandths of an em,
/// starting at the space character.
const HELVETICA_WIDTHS: [u16; 95] = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278,
    278, // ' ' to '/'
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584,
    556, // '0' to '?'
    1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722,
    778, // '@' to 'O'
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469,
    556, // 'P' to '_'
    333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556,
    556, // '`' to 'o'
    556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584, // 'p' to '~'
];

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
    fn a_finite_number_too_large_to_scale_is_still_written_as_a_number() {
        // `v * 1000.0` overflows above f64::MAX / 1000 = 1.7976931348623157e305,
        // and the literal `inf` is not an SVG number.
        for v in [f64::MAX, 1e308, -1e308, 1.8e305, -1.8e305] {
            let out = num(v);
            assert!(!out.contains("inf"), "num({v}) = {out}");
            assert!(out.parse::<f64>().is_ok(), "num({v}) = {out}");
        }
        assert_eq!(num(1e308).len(), 309);
        // One ulp below the threshold keeps the answer it always gave.
        assert_eq!(num(1.7976931348623156e305).len(), 306);
    }

    #[test]
    fn escaping_covers_every_xml_metacharacter() {
        assert_eq!(escape(r#"a&b<c>d"e'f"#), "a&amp;b&lt;c&gt;d&quot;e&apos;f");
    }

    #[test]
    fn a_control_character_in_a_name_is_dropped_rather_than_written() {
        // XML 1.0 has no way to write these, escaped or not, so a name that
        // carries one would otherwise produce a document no parser will open.
        assert_eq!(escape("a\u{1}b"), "ab");
        assert_eq!(escape("rpoB\u{0}\u{7}\u{b}\u{c}\u{1b}"), "rpoB");
        assert_eq!(escape("a\u{fffe}b\u{ffff}c"), "abc");
        // Tab, newline and carriage return are characters, and so is everything
        // from U+0020 up.
        assert_eq!(escape("a\tb\nc\rd"), "a\tb\nc\rd");
        assert_eq!(escape("gene\u{200f}\u{1f9ec}"), "gene\u{200f}\u{1f9ec}");
    }

    #[test]
    fn a_font_size_at_or_below_zero_draws_nothing_instead_of_a_negative_length() {
        // A track that derives its label size from the theme can reach -1 from
        // a theme font size of 0, and `font-size` is not allowed to be negative.
        let mut svg = SvgWriter::new();
        svg.text(10.0, 10.0, "depth", "#6b7280", -1.0, Anchor::End);
        svg.text_bold(10.0, 10.0, "title", "#6b7280", 0.0, Anchor::Start);
        svg.text_rotated(
            (10.0, 10.0),
            -90.0,
            "column",
            "#6b7280",
            -3.0,
            Anchor::Middle,
        );
        svg.text(10.0, 10.0, "nan", "#6b7280", f64::NAN, Anchor::Start);
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        assert!(!out.contains("<text"), "{out}");
        // A positive size still writes its element.
        let mut svg = SvgWriter::new();
        svg.text(10.0, 10.0, "depth", "#6b7280", 1.0, Anchor::End);
        let out = svg.finish(20.0, 20.0, "none", "sans-serif");
        assert!(out.contains(r#"font-size="1""#), "{out}");
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
    fn text_width_matches_what_a_renderer_actually_draws() {
        // Measured with getComputedTextLength in a browser, at font-size 9.
        assert!((text_width("ERR5001", 9.0) - 39.015).abs() < 0.01);
        assert!((text_width("SNP distance", 9.0) - 54.522).abs() < 0.01);
        // Capitals are far wider than lowercase, which is the whole reason for
        // the table: one flat factor clips a column of accessions.
        assert!(text_width("MMMM", 10.0) > text_width("iiii", 10.0) * 3.0);
    }

    #[test]
    fn an_unknown_character_reserves_too_much_rather_than_too_little() {
        assert!(text_width("ñ", 10.0) > text_width("n", 10.0));
        assert_eq!(text_width("", 10.0), 0.0);
    }

    #[test]
    fn an_id_prefix_keeps_two_documents_from_claiming_the_same_id() {
        let one = {
            let mut svg = SvgWriter::with_id_prefix("p0-");
            svg.begin_clip(0.0, 0.0, 10.0, 10.0);
            svg.finish(10.0, 10.0, "none", "sans-serif")
        };
        let two = {
            let mut svg = SvgWriter::with_id_prefix("p1-");
            svg.begin_clip(0.0, 0.0, 10.0, 10.0);
            svg.finish(10.0, 10.0, "none", "sans-serif")
        };
        assert!(one.contains(r#"id="p0-karyon-clip-0""#), "{one}");
        assert!(two.contains(r#"id="p1-karyon-clip-0""#), "{two}");
        // Nested into one document these would otherwise collide, and the
        // second one's clip would resolve to the first one's rectangle.
        assert!(!two.contains("p0-"));
    }

    #[test]
    fn transparent_background_writes_no_page_rect() {
        let svg = SvgWriter::new();
        let out = svg.finish(10.0, 10.0, "none", "sans-serif");
        assert!(!out.contains("<rect"));
    }
}
