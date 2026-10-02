//! The interpreter: reads the document twice and writes what it draws.
//!
//! The first reading collects what a `url(#id)` can name, clips and
//! gradients, because SVG lets a reference come before what it refers to.
//! The second walks the elements in order and writes each as PDF operators
//! straight into the document, so a figure of a million marks is never held
//! twice.
//!
//! An SVG element carries its whole paint with it and a PDF operator does
//! not: `rg` sets the fill for everything after it. So the walk keeps the
//! state the PDF is in and writes an operator only when an element needs a
//! different value, and a thousand bars of one colour set it once. The state
//! is a stack, because `q` and `Q` save and restore all of it, and every
//! group that clips or moves is drawn inside a `q`.

use std::collections::HashMap;

use super::color::{self, Color};
use super::path::{self, Matrix, Scanner, Shape};
use super::text::{self, Face};
use super::xml::{self, Attributes, Token};
use super::Notes;
use crate::svg::num;

/// How many rows a fade's image has: one per level of an eight-bit alpha,
/// so a reader that does not smooth the image still shows no steps.
const FADE_ROWS: usize = 256;

/// Whether this reader acts on `attribute` of `element`.
///
/// Everything else is passed over and named in a note, so a writer that
/// starts emitting something new is heard about rather than quietly dropped:
/// the lock test in `tests.rs` holds every method of
/// [`SvgWriter`](crate::SvgWriter) to an empty list of notes.
pub(crate) fn reads(element: &str, attribute: &str) -> bool {
    // What the SVG says to a screen reader or a pointer. It draws nothing,
    // and the document's title and description reach the PDF's own
    // properties by another route.
    if matches!(
        attribute,
        "id" | "role" | "pointer-events" | "xmlns" | "xmlns:xlink" | "version" | "focusable"
    ) || attribute.starts_with("aria-")
    {
        return true;
    }
    let inherited = matches!(
        attribute,
        "fill"
            | "fill-opacity"
            | "stroke"
            | "stroke-opacity"
            | "stroke-width"
            | "stroke-dasharray"
            | "stroke-dashoffset"
            | "stroke-linejoin"
            | "stroke-linecap"
            | "stroke-miterlimit"
            | "font-family"
            | "font-size"
            | "font-weight"
            | "text-anchor"
            | "letter-spacing"
    );
    let placed = inherited || matches!(attribute, "transform" | "clip-path");
    match element {
        "svg" => inherited || matches!(attribute, "width" | "height" | "viewBox" | "x" | "y"),
        "g" => placed,
        "rect" => placed || matches!(attribute, "x" | "y" | "width" | "height" | "rx" | "ry"),
        "circle" => placed || matches!(attribute, "cx" | "cy" | "r"),
        "ellipse" => placed || matches!(attribute, "cx" | "cy" | "rx" | "ry"),
        "line" => placed || matches!(attribute, "x1" | "y1" | "x2" | "y2"),
        "polyline" | "polygon" => placed || attribute == "points",
        "path" => placed || attribute == "d",
        "text" => placed || matches!(attribute, "x" | "y" | "textLength" | "lengthAdjust"),
        "linearGradient" => matches!(attribute, "x1" | "y1" | "x2" | "y2"),
        "stop" => matches!(attribute, "offset" | "stop-color" | "stop-opacity"),
        _ => false,
    }
}

/// Elements this reader draws, or reads for what they define.
pub(crate) fn known(element: &str) -> bool {
    matches!(
        element,
        "svg"
            | "g"
            | "defs"
            | "title"
            | "desc"
            | "metadata"
            | "clipPath"
            | "linearGradient"
            | "stop"
            | "rect"
            | "circle"
            | "ellipse"
            | "line"
            | "polyline"
            | "polygon"
            | "path"
            | "text"
    )
}

/// Names each attribute of `element` this reader does not act on.
fn check(notes: &mut Notes, element: &str, attributes: Attributes<'_>) {
    for (attribute, _) in attributes.iter() {
        if !reads(element, attribute) {
            notes.add(format!(
                "the {attribute} attribute of <{element}> is not read, so the element is \
                 drawn as if it had none"
            ));
        }
    }
}

/// A gradient, in the box of the shape it fills.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Gradient {
    from: (f64, f64),
    to: (f64, f64),
    stops: Vec<(f64, Color)>,
}

impl Gradient {
    /// The colour at `t` along the gradient, padded past either end.
    fn at(&self, t: f64) -> Color {
        let first = self.stops.first().map_or(Color::BLACK, |(_, c)| *c);
        let mut before = (0.0, first);
        for (offset, color) in &self.stops {
            if t < *offset {
                let (start, from) = before;
                if *offset <= start {
                    return *color;
                }
                let share = ((t - start) / (offset - start)).clamp(0.0, 1.0);
                let mix = |a: f64, b: f64| a + (b - a) * share;
                return Color {
                    rgb: [
                        mix(from.rgb[0], color.rgb[0]),
                        mix(from.rgb[1], color.rgb[1]),
                        mix(from.rgb[2], color.rgb[2]),
                    ],
                    alpha: mix(from.alpha, color.alpha),
                };
            }
            before = (*offset, *color);
        }
        before.1
    }
}

/// What a `url(#id)` can name, collected before anything is drawn.
#[derive(Debug, Default)]
pub(crate) struct Definitions<'a> {
    clips: HashMap<&'a str, Vec<Shape>>,
    gradients: HashMap<&'a str, Gradient>,
}

impl<'a> Definitions<'a> {
    /// The first reading: every clip and every linear gradient, by id.
    pub(crate) fn read(document: &'a str, notes: &mut Notes) -> Definitions<'a> {
        let mut definitions = Definitions::default();
        let mut clip: Option<(&'a str, Vec<Shape>)> = None;
        let mut gradient: Option<(&'a str, Gradient)> = None;
        for token in xml::tokens(document) {
            match token {
                Token::Start {
                    name,
                    attributes,
                    empty,
                } => match name {
                    "clipPath" => {
                        check(notes, name, attributes);
                        let id = attributes.get("id").unwrap_or_default();
                        if empty {
                            definitions.clips.insert(id, Vec::new());
                        } else {
                            clip = Some((id, Vec::new()));
                        }
                    }
                    "linearGradient" => {
                        check(notes, name, attributes);
                        let fraction = |attribute: &str, default: f64| {
                            attributes
                                .get(attribute)
                                .and_then(proportion)
                                .unwrap_or(default)
                        };
                        let found = Gradient {
                            from: (fraction("x1", 0.0), fraction("y1", 0.0)),
                            to: (fraction("x2", 1.0), fraction("y2", 0.0)),
                            stops: Vec::new(),
                        };
                        let id = attributes.get("id").unwrap_or_default();
                        if empty {
                            definitions.gradients.insert(id, found);
                        } else {
                            gradient = Some((id, found));
                        }
                    }
                    "stop" => {
                        if let Some((_, found)) = &mut gradient {
                            check(notes, name, attributes);
                            let next = stop(attributes, found.stops.last(), notes);
                            found.stops.push(next);
                        }
                    }
                    _ => {
                        if let Some((_, shapes)) = &mut clip {
                            check(notes, name, attributes);
                            if let Some(shape) = geometry(name, attributes, notes) {
                                let shape = match attributes.get("transform") {
                                    Some(list) => match path::transform(&xml::unescape(list)) {
                                        Some(matrix) => shape.transformed(&matrix),
                                        None => shape,
                                    },
                                    None => shape,
                                };
                                shapes.push(shape);
                            } else if !known(name) {
                                notes.add(format!(
                                    "the <{name}> element is not read, so the clip it is in \
                                     leaves it out"
                                ));
                            }
                        }
                    }
                },
                Token::End("clipPath") => {
                    if let Some((id, shapes)) = clip.take() {
                        definitions.clips.insert(id, shapes);
                    }
                }
                Token::End("linearGradient") => {
                    if let Some((id, found)) = gradient.take() {
                        definitions.gradients.insert(id, found);
                    }
                }
                _ => {}
            }
        }
        definitions
    }
}

/// One stop of a gradient, its offset kept from going back past the stop
/// before it, as SVG requires.
fn stop(
    attributes: Attributes<'_>,
    before: Option<&(f64, Color)>,
    notes: &mut Notes,
) -> (f64, Color) {
    let offset = attributes
        .get("offset")
        .and_then(proportion)
        .unwrap_or(0.0)
        .clamp(0.0, 1.0)
        .max(before.map_or(0.0, |(offset, _)| *offset));
    let mut color = match attributes.get("stop-color") {
        Some(value) => color::parse(value).unwrap_or_else(|| {
            notes.add(unread_paint(value));
            Color::BLACK
        }),
        None => Color::BLACK,
    };
    if let Some(opacity) = attributes.get("stop-opacity").and_then(opacity) {
        color.alpha *= opacity;
    }
    (offset, color)
}

/// A number, or a percentage of one.
fn proportion(value: &str) -> Option<f64> {
    let value = value.trim();
    let number = match value.strip_suffix('%') {
        Some(percent) => percent.trim().parse::<f64>().ok()? / 100.0,
        None => value.parse::<f64>().ok()?,
    };
    number.is_finite().then_some(number)
}

/// An opacity, clamped to 0 to 1 as CSS clamps one.
fn opacity(value: &str) -> Option<f64> {
    proportion(value).map(|v| v.clamp(0.0, 1.0))
}

/// Pixels in one of each absolute unit CSS has, at its 96 pixels to the
/// inch.
const UNITS: [(&str, f64); 6] = [
    ("px", 1.0),
    ("pt", 4.0 / 3.0),
    ("pc", 16.0),
    ("in", 96.0),
    ("cm", 96.0 / 2.54),
    ("mm", 96.0 / 25.4),
];

/// A length in pixels, the unit a karyon figure is drawn in, from a bare
/// number or a number in any absolute unit CSS has.
///
/// This is the one reader of a length, the page's size included, so a root
/// `width="100mm"` scales the drawing by as much as it sizes the page. A
/// percentage, or a size relative to the font, is `None`: each needs a box or
/// a font to be measured against, and whoever asked names it in a note.
pub(crate) fn length(value: &str) -> Option<f64> {
    let value = value.trim();
    let (number, factor) = UNITS
        .iter()
        .find_map(|(unit, factor)| value.strip_suffix(unit).map(|number| (number, *factor)))
        .unwrap_or((value, 1.0));
    let pixels = number.trim_end().parse::<f64>().ok()? * factor;
    pixels.is_finite().then_some(pixels)
}

/// `value`, written as `attribute` of `element`, as a length, or `None` with
/// a note when it is in a unit this does not read, so the element is drawn
/// as if it had no such attribute.
fn placed(notes: &mut Notes, element: &str, attribute: &str, value: &str) -> Option<f64> {
    let read = length(value);
    if read.is_none() {
        notes.add(format!(
            "the {attribute} {value:?} of <{element}> is not read, so the element is drawn as \
             if it had none"
        ));
    }
    read
}

/// `value`, written as `attribute`, as a length an element can inherit, or
/// `None` with a note when it is in a unit this does not read, so the
/// element keeps the one it inherits.
fn own(notes: &mut Notes, attribute: &str, value: &str) -> Option<f64> {
    let read = length(value);
    if read.is_none() {
        notes.add(format!(
            "the {attribute} {value:?} is not read, so the element takes the one it inherits"
        ));
    }
    read
}

/// The width and height an `<svg>` gives what is inside it, in pixels: its
/// own, or its `viewBox`'s when it lacks either.
///
/// The page is sized by this too, so the page and the drawing on it are read
/// alike. A negative side is nought, as SVG draws nothing for one.
pub(crate) fn viewport(attributes: Attributes<'_>, notes: &mut Notes) -> Option<(f64, f64)> {
    let mut side = |name: &str| {
        let value = attributes.get(name)?;
        placed(notes, "svg", name, value).map(|pixels| pixels.max(0.0))
    };
    let (width, height) = (side("width"), side("height"));
    match (width, height, attributes.get("viewBox").and_then(view_box)) {
        (Some(w), Some(h), _) => Some((w, h)),
        (_, _, Some([_, _, w, h])) => Some((w, h)),
        _ => None,
    }
}

/// The id a `url(#id)` names, and whatever is written after it.
fn reference(value: &str) -> Option<(&str, &str)> {
    let inside = value.trim().strip_prefix("url(")?;
    let close = inside.find(')')?;
    let id = inside[..close]
        .trim()
        .trim_matches(|c: char| c == '"' || c == '\'');
    Some((id.strip_prefix('#')?, inside[close + 1..].trim()))
}

fn unread_paint(value: &str) -> String {
    format!("the colour {value:?} is not read, so what it paints takes the colour it inherits")
}

/// The outline of a shape element, or `None` when SVG would draw nothing.
fn geometry(element: &str, attributes: Attributes<'_>, notes: &mut Notes) -> Option<Shape> {
    let mut get = |name: &str| {
        let value = attributes.get(name)?;
        placed(notes, element, name, value)
    };
    let shape = match element {
        "rect" => {
            let (w, h) = (get("width")?, get("height")?);
            if !(w > 0.0 && h > 0.0) {
                return None;
            }
            let (x, y) = (get("x").unwrap_or(0.0), get("y").unwrap_or(0.0));
            Shape::rect(x, y, w, h, get("rx"), get("ry"))
        }
        "circle" => {
            let r = get("r")?;
            if r <= 0.0 {
                return None;
            }
            Shape::ellipse(get("cx").unwrap_or(0.0), get("cy").unwrap_or(0.0), r, r)
        }
        "ellipse" => {
            let (rx, ry) = (get("rx")?, get("ry")?);
            if rx <= 0.0 || ry <= 0.0 {
                return None;
            }
            Shape::ellipse(get("cx").unwrap_or(0.0), get("cy").unwrap_or(0.0), rx, ry)
        }
        "line" => {
            let mut at = |name: &str| get(name).unwrap_or(0.0);
            Shape::Path(vec![
                path::Segment::Move(at("x1"), at("y1")),
                path::Segment::Line(at("x2"), at("y2")),
            ])
        }
        "polyline" => Shape::points(attributes.get("points").unwrap_or_default(), false),
        "polygon" => Shape::points(attributes.get("points").unwrap_or_default(), true),
        "path" => Shape::Path(path::parse(attributes.get("d").unwrap_or_default())),
        _ => return None,
    };
    (!shape.is_empty()).then_some(shape)
}

/// What fills or strokes a shape, as written.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Paint<'a> {
    None,
    Color(Color),
    /// A gradient by id, and what to paint if the id names none.
    Url(&'a str, Option<Color>),
}

/// How far apart letters are set.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Spacing {
    Em(f64),
    Px(f64),
}

/// The paint and type an element inherits, overridden by its own attributes.
#[derive(Debug, Clone, Copy)]
struct Style<'a> {
    fill: Paint<'a>,
    fill_opacity: f64,
    stroke: Paint<'a>,
    stroke_opacity: f64,
    stroke_width: f64,
    dash: Option<&'a str>,
    dash_offset: f64,
    join: u8,
    cap: u8,
    miter: f64,
    family: &'a str,
    size: f64,
    weight: &'a str,
    anchor: u8,
    spacing: Spacing,
}

impl<'a> Style<'a> {
    /// SVG's initial values.
    fn initial() -> Style<'a> {
        Style {
            fill: Paint::Color(Color::BLACK),
            fill_opacity: 1.0,
            stroke: Paint::None,
            stroke_opacity: 1.0,
            stroke_width: 1.0,
            dash: None,
            dash_offset: 0.0,
            join: 0,
            cap: 0,
            miter: 4.0,
            family: "sans-serif",
            size: 16.0,
            weight: "normal",
            anchor: 0,
            spacing: Spacing::Px(0.0),
        }
    }
}

/// An element the walk is inside of.
struct Frame<'a> {
    name: &'a str,
    /// Whether it opened a `q` that its end has to close.
    restore: bool,
    style: Style<'a>,
}

/// The PDF graphics state, as far as this writer sets it.
#[derive(Debug, Clone, PartialEq)]
struct State {
    fill: [f64; 3],
    stroke: [f64; 3],
    width: f64,
    join: u8,
    cap: u8,
    miter: f64,
    dash: (Vec<f64>, f64),
    alpha: (f64, f64),
    font: Option<(usize, f64)>,
    spacing: f64,
    stretch: f64,
    render: u8,
}

impl State {
    /// A page's initial state, after the `4 M` every content stream here
    /// opens with: SVG's miter limit is 4 and PDF's is 10.
    fn initial() -> State {
        State {
            fill: [0.0; 3],
            stroke: [0.0; 3],
            width: 1.0,
            join: 0,
            cap: 0,
            miter: 4.0,
            dash: (Vec::new(), 0.0),
            alpha: (1.0, 1.0),
            font: None,
            spacing: 0.0,
            stretch: 100.0,
            render: 0,
        }
    }
}

/// A fade as two images: the colour, and the alpha it is masked by.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Fade {
    key: (String, u64),
    pub(crate) rgb: Vec<u8>,
    pub(crate) alpha: Vec<u8>,
}

/// What the content stream names and the page has to define.
#[derive(Debug, Default)]
pub(crate) struct Resources {
    /// Fonts, `/F1` onwards, in order of first use.
    pub(crate) fonts: Vec<Face>,
    /// Constant alphas for fill and stroke, `/G1` onwards.
    pub(crate) alphas: Vec<(f64, f64)>,
    /// Fades, `/I1` onwards.
    pub(crate) fades: Vec<Fade>,
}

/// A `<text>` being read, drawn when it closes.
struct Label<'a> {
    attributes: Attributes<'a>,
    style: Style<'a>,
    content: String,
}

/// Where the second reading is.
pub(crate) struct Painter<'a, 'w> {
    out: &'w mut Vec<u8>,
    definitions: &'w Definitions<'a>,
    pub(crate) notes: &'w mut Notes,
    pub(crate) resources: Resources,
    pub(crate) title: Option<String>,
    pub(crate) description: Option<String>,
    frames: Vec<Frame<'a>>,
    states: Vec<State>,
    missing: Vec<char>,
}

impl<'a, 'w> Painter<'a, 'w> {
    pub(crate) fn new(
        out: &'w mut Vec<u8>,
        definitions: &'w Definitions<'a>,
        notes: &'w mut Notes,
    ) -> Painter<'a, 'w> {
        Painter {
            out,
            definitions,
            notes,
            resources: Resources::default(),
            title: None,
            description: None,
            frames: Vec::new(),
            states: vec![State::initial()],
            missing: Vec::new(),
        }
    }

    /// The second reading.
    pub(crate) fn run(&mut self, document: &'a str) {
        // How deep the walk is inside something it does not draw.
        let mut skipping = 0usize;
        // The root's own title or description, while it is being read.
        let mut capture: Option<(bool, String)> = None;
        let mut label: Option<Label<'a>> = None;
        for token in xml::tokens(document) {
            match token {
                Token::Start {
                    name,
                    attributes,
                    empty,
                } => {
                    if skipping > 0 {
                        skipping += usize::from(!empty);
                        continue;
                    }
                    if label.is_some() {
                        self.notes.add(format!(
                            "a <{name}> inside <text> is not read, so its letters are left out"
                        ));
                        skipping = usize::from(!empty);
                        continue;
                    }
                    if !known(name) {
                        self.notes.add(format!("the <{name}> element is not drawn"));
                        skipping = usize::from(!empty);
                        continue;
                    }
                    check(self.notes, name, attributes);
                    match name {
                        "svg" => self.open_svg(attributes, empty),
                        "g" => self.open_group(attributes, empty),
                        "title" | "desc" => {
                            // Only the root's own: a group's title is a
                            // tooltip, and a PDF has nowhere to hover.
                            let root = self.frames.len() == 1;
                            let title = name == "title";
                            let wanted = if title {
                                self.title.is_none()
                            } else {
                                self.description.is_none()
                            };
                            if root && wanted && !empty {
                                capture = Some((title, String::new()));
                            }
                            skipping = usize::from(!empty);
                        }
                        "text" => {
                            if !empty {
                                label = Some(Label {
                                    attributes,
                                    style: self.inherit(attributes),
                                    content: String::new(),
                                });
                            }
                        }
                        "defs" | "metadata" | "clipPath" | "linearGradient" | "stop" => {
                            skipping = usize::from(!empty);
                        }
                        _ => {
                            self.draw_shape(name, attributes);
                            skipping = usize::from(!empty);
                        }
                    }
                }
                Token::End(name) => {
                    if skipping > 0 {
                        skipping -= 1;
                        if skipping == 0 {
                            if let Some((title, text)) = capture.take() {
                                let text = text::collapse(&text);
                                if title {
                                    self.title = Some(text);
                                } else {
                                    self.description = Some(text);
                                }
                            }
                        }
                        continue;
                    }
                    if name == "text" {
                        if let Some(label) = label.take() {
                            self.draw_text(label);
                        }
                        continue;
                    }
                    if label.is_none() {
                        self.close(name);
                    }
                }
                Token::Text(data) | Token::Raw(data) => {
                    let data = match token {
                        Token::Text(_) => xml::unescape(data),
                        _ => data.into(),
                    };
                    if skipping == 1 {
                        if let Some((_, text)) = &mut capture {
                            text.push_str(&data);
                        }
                    } else if skipping == 0 {
                        if let Some(label) = &mut label {
                            label.content.push_str(&data);
                        }
                    }
                }
            }
        }
        // Whatever never closed is closed, as the writer closes what a track
        // left open.
        while let Some(frame) = self.frames.pop() {
            if frame.restore {
                self.restore();
            }
        }
        if !self.missing.is_empty() {
            let named: Vec<String> = self
                .missing
                .iter()
                .map(|c| format!("{c} (U+{:04X})", u32::from(*c)))
                .collect();
            self.notes.add(format!(
                "{} {} in no font every PDF reader has, and {} drawn as a question mark",
                named.join(", "),
                if named.len() == 1 { "is" } else { "are" },
                if named.len() == 1 { "is" } else { "each is" },
            ));
        }
    }

    fn open_svg(&mut self, attributes: Attributes<'a>, empty: bool) {
        let style = self.inherit(attributes);
        let view_box = attributes.get("viewBox").and_then(view_box);
        let size = viewport(attributes, self.notes);
        let restore = if self.frames.is_empty() {
            // The root: the page is its viewport, so there is nothing to clip.
            if let (Some((w, h)), Some(view)) = (size, view_box) {
                let fit = fit(view, 0.0, 0.0, w, h);
                if fit != Matrix::IDENTITY {
                    fit.write(self.out);
                }
            }
            false
        } else {
            // A document inside another, as a sheet holds its panels: SVG
            // clips it to its own viewport, which is what keeps a panel's
            // marks from spilling into the one beside it.
            let mut get = |name: &str| {
                let value = attributes.get(name)?;
                placed(self.notes, "svg", name, value)
            };
            let (x, y) = (get("x").unwrap_or(0.0), get("y").unwrap_or(0.0));
            self.save();
            if let Some((w, h)) = size {
                Shape::Rect { x, y, w, h }.write(self.out);
                self.out.extend_from_slice(b"W n\n");
                let fit = match view_box {
                    Some(view) => fit(view, x, y, w, h),
                    None => Matrix::translate(x, y),
                };
                if fit != Matrix::IDENTITY {
                    fit.write(self.out);
                }
            }
            true
        };
        self.frames.push(Frame {
            name: "svg",
            restore,
            style,
        });
        if empty {
            self.close("svg");
        }
    }

    fn open_group(&mut self, attributes: Attributes<'a>, empty: bool) {
        let style = self.inherit(attributes);
        let restore = self.place(attributes);
        self.frames.push(Frame {
            name: "g",
            restore,
            style,
        });
        if empty {
            self.close("g");
        }
    }

    /// Ends the innermost open element if `name` is the one it is.
    fn close(&mut self, name: &str) {
        if self.frames.last().is_some_and(|frame| frame.name == name) {
            if let Some(frame) = self.frames.pop() {
                if frame.restore {
                    self.restore();
                }
            }
        }
    }

    /// Opens a `q` for an element's transform and clip, when it has either,
    /// and says whether it did.
    fn place(&mut self, attributes: Attributes<'a>) -> bool {
        let transform = attributes.get("transform").and_then(|list| {
            let matrix = path::transform(&xml::unescape(list));
            if matrix.is_none() {
                self.notes.add(format!(
                    "the transform {list:?} is not read, so its element is drawn where it \
                     would be without one"
                ));
            }
            matrix
        });
        let definitions = self.definitions;
        let clip = attributes
            .get("clip-path")
            .filter(|value| value.trim() != "none")
            .and_then(|value| {
                let found = reference(value).and_then(|(id, _)| definitions.clips.get(id));
                if found.is_none() {
                    self.notes.add(format!(
                        "the clip-path {value:?} names no clip in the document, so its \
                         element is not clipped"
                    ));
                }
                found
            });
        if transform.is_none() && clip.is_none() {
            return false;
        }
        self.save();
        if let Some(matrix) = transform {
            matrix.write(self.out);
        }
        if let Some(shapes) = clip {
            if shapes.is_empty() {
                // A clip with nothing in it hides everything it clips.
                self.out.extend_from_slice(b"0 0 0 0 re\n");
            }
            for shape in shapes {
                shape.write(self.out);
            }
            self.out.extend_from_slice(b"W n\n");
        }
        true
    }

    /// The style `attributes` give an element inside the innermost frame.
    fn inherit(&mut self, attributes: Attributes<'a>) -> Style<'a> {
        let mut style = self
            .frames
            .last()
            .map_or_else(Style::initial, |frame| frame.style);
        for (name, value) in attributes.iter() {
            if value.trim() == "inherit" {
                continue;
            }
            match name {
                "fill" => {
                    if let Some(paint) = self.paint(value) {
                        style.fill = paint;
                    }
                }
                "stroke" => {
                    if let Some(paint) = self.paint(value) {
                        style.stroke = paint;
                    }
                }
                "fill-opacity" => style.fill_opacity = opacity(value).unwrap_or(style.fill_opacity),
                "stroke-opacity" => {
                    style.stroke_opacity = opacity(value).unwrap_or(style.stroke_opacity);
                }
                "stroke-width" => {
                    if let Some(width) = own(self.notes, name, value).filter(|w| *w >= 0.0) {
                        style.stroke_width = width;
                    }
                }
                "stroke-dasharray" => {
                    style.dash = Some(value).filter(|value| value.trim() != "none");
                }
                "stroke-dashoffset" => {
                    style.dash_offset = own(self.notes, name, value).unwrap_or(style.dash_offset)
                }
                "stroke-linejoin" => {
                    style.join = match value.trim() {
                        "round" => 1,
                        "bevel" => 2,
                        "miter" | "miter-clip" | "arcs" => 0,
                        _ => style.join,
                    };
                }
                "stroke-linecap" => {
                    style.cap = match value.trim() {
                        "butt" => 0,
                        "round" => 1,
                        "square" => 2,
                        _ => style.cap,
                    };
                }
                "stroke-miterlimit" => {
                    if let Some(limit) = length(value).filter(|limit| *limit >= 1.0) {
                        style.miter = limit;
                    }
                }
                "font-family" => style.family = value,
                "font-size" => {
                    if let Some(size) = own(self.notes, name, value).filter(|size| *size >= 0.0) {
                        style.size = size;
                    }
                }
                "font-weight" => style.weight = value,
                "text-anchor" => {
                    style.anchor = match value.trim() {
                        "start" => 0,
                        "middle" => 1,
                        "end" => 2,
                        _ => style.anchor,
                    };
                }
                "letter-spacing" => {
                    let value = value.trim();
                    style.spacing = if value == "normal" {
                        Spacing::Px(0.0)
                    } else if let Some(em) = value.strip_suffix("em").and_then(|v| v.parse().ok()) {
                        Spacing::Em(em)
                    } else {
                        own(self.notes, name, value).map_or(style.spacing, Spacing::Px)
                    };
                }
                _ => {}
            }
        }
        style
    }

    /// A paint as written, or `None`, with a note, for one this cannot read.
    fn paint(&mut self, value: &'a str) -> Option<Paint<'a>> {
        let trimmed = value.trim();
        if trimmed == "none" {
            return Some(Paint::None);
        }
        if trimmed.eq_ignore_ascii_case("currentcolor") {
            // The `color` property is not read here, and its initial value
            // is black, which is what a browser paints with it unset.
            return Some(Paint::Color(Color::BLACK));
        }
        if let Some((id, fallback)) = reference(trimmed) {
            let fallback = match fallback {
                "" | "none" => None,
                other => match color::parse(other) {
                    Some(color) => Some(color),
                    None => {
                        self.notes.add(unread_paint(other));
                        None
                    }
                },
            };
            return Some(Paint::Url(id, fallback));
        }
        match color::parse(trimmed) {
            Some(color) => Some(Paint::Color(color)),
            None => {
                self.notes.add(unread_paint(trimmed));
                None
            }
        }
    }

    /// A paint brought down to one colour and its opacity, or `None` for no
    /// paint. A gradient where one colour has to do is its first stop.
    fn flat(&mut self, paint: Paint<'a>, opacity: f64, gradients: bool) -> Option<Color> {
        let definitions = self.definitions;
        let color = match paint {
            Paint::None => return None,
            Paint::Color(color) => color,
            Paint::Url(id, fallback) => match definitions.gradients.get(id) {
                Some(gradient) => {
                    if !gradients {
                        self.notes.add(
                            "a gradient on a stroke or on text is drawn in its first colour"
                                .to_string(),
                        );
                    }
                    gradient.stops.first().map(|(_, color)| *color)?
                }
                None => {
                    if fallback.is_none() {
                        self.notes.add(format!(
                            "url(#{id}) names no gradient in the document, so it paints nothing"
                        ));
                    }
                    fallback?
                }
            },
        };
        let alpha = color.alpha * opacity;
        (alpha > 0.0).then_some(Color { alpha, ..color })
    }

    fn draw_shape(&mut self, element: &'a str, attributes: Attributes<'a>) {
        let style = self.inherit(attributes);
        let Some(shape) = geometry(element, attributes, self.notes) else {
            return;
        };
        let restore = self.place(attributes);
        // A line has no inside to fill.
        let mut fill = if element == "line" {
            Paint::None
        } else {
            style.fill
        };
        if let Paint::Url(id, _) = fill {
            let definitions = self.definitions;
            if let Some(gradient) = definitions.gradients.get(id) {
                self.fade(&shape, id, gradient, style.fill_opacity);
                fill = Paint::None;
            }
        }
        let fill = self.flat(fill, style.fill_opacity, true);
        let stroke = if style.stroke_width > 0.0 {
            self.flat(style.stroke, style.stroke_opacity, false)
        } else {
            // PDF draws a width of nought as the thinnest line the device
            // has, and SVG draws nothing.
            None
        };
        if fill.is_some() || stroke.is_some() {
            if let Some(fill) = fill {
                self.set_fill(fill);
            }
            if let Some(stroke) = stroke {
                self.set_stroke(stroke, &style);
            }
            shape.write(self.out);
            self.out.extend_from_slice(match (fill, stroke) {
                (Some(_), Some(_)) => b"B\n",
                (Some(_), None) => b"f\n",
                _ => b"S\n",
            });
        }
        if restore {
            self.restore();
        }
    }

    /// Paints `shape` with a gradient: an image of the colour, its alpha
    /// in a soft mask, stretched over the shape's box and clipped to it.
    ///
    /// An image rather than a shading because a fade here changes its
    /// opacity, and opacity along a shading is a soft-mask group. That was
    /// measured on a coverage figure and failed twice: Quartz, which Preview
    /// draws with, cut it off two thirds of the way down, and Inkscape's
    /// import dropped it. An image with an alpha mask drew right in
    /// Ghostscript, poppler, Quartz and Inkscape.
    fn fade(&mut self, shape: &Shape, id: &str, gradient: &Gradient, opacity: f64) {
        let Some([left, top, right, bottom]) = shape.bounds() else {
            return;
        };
        let (w, h) = (right - left, bottom - top);
        // A gradient over the box of a shape with no width or no height has
        // no box to run over, and SVG draws nothing.
        if !(w > 0.0 && h > 0.0) || gradient.stops.is_empty() {
            return;
        }
        let (gx, gy) = (
            gradient.to.0 - gradient.from.0,
            gradient.to.1 - gradient.from.1,
        );
        let length_squared = gx * gx + gy * gy;
        if length_squared <= 0.0 || length_squared.is_nan() {
            // A gradient that goes nowhere is its last colour throughout.
            let last = gradient.stops.last().map(|(_, color)| *color);
            if let Some(color) = last.filter(|color| color.alpha * opacity > 0.0) {
                self.set_fill(Color {
                    alpha: color.alpha * opacity,
                    ..color
                });
                shape.write(self.out);
                self.out.extend_from_slice(b"f\n");
            }
            return;
        }
        // Where the four corners of the box fall along the gradient and
        // across it, in the unit square SVG lays a gradient out in.
        let length = length_squared.sqrt();
        let (nx, ny) = (-gy / length, gx / length);
        let (mut t_min, mut t_max, mut s_min, mut s_max) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            let (du, dv) = (u - gradient.from.0, v - gradient.from.1);
            let t = (du * gx + dv * gy) / length_squared;
            let s = du * nx + dv * ny;
            t_min = t_min.min(t);
            t_max = t_max.max(t);
            s_min = s_min.min(s);
            s_max = s_max.max(s);
        }
        let image = self.fade_image(id, gradient, opacity, t_min, t_max);
        self.save();
        // The alpha is in the mask, so the constant alpha has to be one.
        self.set_alpha(Some(1.0), None);
        shape.write(self.out);
        self.out.extend_from_slice(b"W n\n");
        Matrix([w, 0.0, 0.0, h, left, top]).write(self.out);
        // Image space is a unit square whose top row is drawn at y = 1. The
        // rows run along the gradient from its start, and the single column
        // is stretched across it.
        let span = t_max - t_min;
        let across = s_max - s_min;
        Matrix([
            across * nx,
            across * ny,
            -span * gx,
            -span * gy,
            gradient.from.0 + t_max * gx + s_min * nx,
            gradient.from.1 + t_max * gy + s_min * ny,
        ])
        .write(self.out);
        self.out
            .extend_from_slice(format!("/I{} Do\n", image + 1).as_bytes());
        self.restore();
    }

    /// The index of the image pair for a gradient at an opacity, made the
    /// first time it is asked for.
    fn fade_image(
        &mut self,
        id: &str,
        gradient: &Gradient,
        opacity: f64,
        from: f64,
        to: f64,
    ) -> usize {
        let key = (id.to_string(), opacity.to_bits());
        if let Some(at) = self.resources.fades.iter().position(|fade| fade.key == key) {
            return at;
        }
        let mut rgb = Vec::with_capacity(FADE_ROWS * 3);
        let mut alpha = Vec::with_capacity(FADE_ROWS);
        let byte = |v: f64| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
        for row in 0..FADE_ROWS {
            let t = from + (to - from) * (row as f64 + 0.5) / FADE_ROWS as f64;
            let color = gradient.at(t.clamp(0.0, 1.0));
            rgb.extend(color.rgb.map(byte));
            alpha.push(byte(color.alpha * opacity));
        }
        self.resources.fades.push(Fade { key, rgb, alpha });
        self.resources.fades.len() - 1
    }

    fn draw_text(&mut self, label: Label<'a>) {
        let style = label.style;
        let content = text::collapse(&label.content);
        if content.is_empty() || style.size <= 0.0 || style.size.is_nan() {
            return;
        }
        let face = Face::choose(&xml::unescape(style.family), style.weight.trim());
        let runs = text::runs(&content, face, &mut self.missing);
        let count = runs.iter().map(|run| run.bytes.len()).sum::<usize>() as f64;
        if count == 0.0 {
            return;
        }
        let natural = runs.iter().map(|run| run.width).sum::<f64>() / 1000.0 * style.size;
        let mut spacing = match style.spacing {
            Spacing::Em(em) => em * style.size,
            Spacing::Px(px) => px,
        };
        // CSS adds the spacing after every letter, the last one included, and
        // so does PDF's `Tc`, so the two agree on where a spaced label ends.
        let mut advance = natural + spacing * count;
        let mut stretch = 100.0;
        if let Some(target) = label
            .attributes
            .get("textLength")
            .and_then(|value| placed(self.notes, "text", "textLength", value))
            .filter(|target| *target > 0.0)
        {
            if label.attributes.get("lengthAdjust") == Some("spacingAndGlyphs") {
                if advance > 0.0 {
                    stretch = 100.0 * target / advance;
                }
            } else {
                spacing += (target - advance) / count;
            }
            advance = target;
        }
        let shift = match style.anchor {
            1 => -advance / 2.0,
            2 => -advance,
            _ => 0.0,
        };
        // A position on a label is a list of lengths, one for each letter,
        // and only the first is read.
        let mut first = |name: &str| {
            let list = label.attributes.get(name).unwrap_or_default();
            let mut items = list
                .split(|c: char| c == ',' || c.is_whitespace())
                .filter(|item| !item.is_empty());
            let value = items
                .next()
                .and_then(|item| placed(self.notes, "text", name, item));
            (value.unwrap_or(0.0), items.next().is_some())
        };
        let ((x, more_x), (y, more_y)) = (first("x"), first("y"));
        if more_x || more_y {
            self.notes.add(
                "a position for each letter of a label is not read, so the label is set from \
                 its first"
                    .to_string(),
            );
        }
        let fill = self.flat(style.fill, style.fill_opacity, false);
        let stroke = if style.stroke_width > 0.0 {
            self.flat(style.stroke, style.stroke_opacity, false)
        } else {
            None
        };
        // A halo is the label stroked as well as filled, which PDF does in
        // one pass with render mode 2.
        let render = match (fill.is_some(), stroke.is_some()) {
            (true, false) => 0,
            (false, true) => 1,
            (true, true) => 2,
            (false, false) => return,
        };
        let restore = self.place(label.attributes);
        if let Some(fill) = fill {
            self.set_fill(fill);
        }
        if let Some(stroke) = stroke {
            self.set_stroke(stroke, &style);
        }
        self.set_text_state(spacing, stretch, render);
        self.out.extend_from_slice(b"BT\n");
        // The label's own space is flipped back upright, and nothing else:
        // whatever turns or moves it was written as `cm` above, so the matrix
        // here never scales and a halo's width means the same in every reader.
        self.out.extend_from_slice(b"1 0 0 -1 ");
        path::push_numbers(self.out, &[x + shift, y]);
        self.out.extend_from_slice(b"Tm\n");
        for run in &runs {
            let font = self.font(run.face);
            let current = self.state().font;
            if current != Some((font, style.size)) {
                self.out.extend_from_slice(
                    format!("/F{} {} Tf\n", font + 1, num(style.size)).as_bytes(),
                );
                self.state_mut().font = Some((font, style.size));
            }
            text::write_string(self.out, &run.bytes);
            self.out.extend_from_slice(b" Tj\n");
        }
        self.out.extend_from_slice(b"ET\n");
        if restore {
            self.restore();
        }
    }

    fn font(&mut self, face: Face) -> usize {
        match self.resources.fonts.iter().position(|known| *known == face) {
            Some(at) => at,
            None => {
                self.resources.fonts.push(face);
                self.resources.fonts.len() - 1
            }
        }
    }

    fn state(&self) -> &State {
        self.states.last().expect("the page state is never popped")
    }

    fn state_mut(&mut self) -> &mut State {
        self.states
            .last_mut()
            .expect("the page state is never popped")
    }

    fn save(&mut self) {
        self.out.extend_from_slice(b"q\n");
        let state = self.state().clone();
        self.states.push(state);
    }

    fn restore(&mut self) {
        if self.states.len() > 1 {
            self.out.extend_from_slice(b"Q\n");
            self.states.pop();
        }
    }

    fn set_fill(&mut self, color: Color) {
        self.set_alpha(Some(color.alpha), None);
        if self.state().fill != color.rgb {
            path::push_numbers(self.out, &color.rgb);
            self.out.extend_from_slice(b"rg\n");
            self.state_mut().fill = color.rgb;
        }
    }

    fn set_stroke(&mut self, color: Color, style: &Style<'_>) {
        self.set_alpha(None, Some(color.alpha));
        if self.state().stroke != color.rgb {
            path::push_numbers(self.out, &color.rgb);
            self.out.extend_from_slice(b"RG\n");
            self.state_mut().stroke = color.rgb;
        }
        if self.state().width != style.stroke_width {
            path::push_numbers(self.out, &[style.stroke_width]);
            self.out.extend_from_slice(b"w\n");
            self.state_mut().width = style.stroke_width;
        }
        if self.state().join != style.join {
            self.out
                .extend_from_slice(format!("{} j\n", style.join).as_bytes());
            self.state_mut().join = style.join;
        }
        if self.state().cap != style.cap {
            self.out
                .extend_from_slice(format!("{} J\n", style.cap).as_bytes());
            self.state_mut().cap = style.cap;
        }
        if self.state().miter != style.miter {
            path::push_numbers(self.out, &[style.miter]);
            self.out.extend_from_slice(b"M\n");
            self.state_mut().miter = style.miter;
        }
        let dash = (dashes(style.dash, self.notes), style.dash_offset);
        let dash = if dash.0.is_empty() {
            (Vec::new(), 0.0)
        } else {
            dash
        };
        if self.state().dash != dash {
            self.out.push(b'[');
            for (at, value) in dash.0.iter().enumerate() {
                if at > 0 {
                    self.out.push(b' ');
                }
                self.out.extend_from_slice(num(*value).as_bytes());
            }
            self.out.extend_from_slice(b"] ");
            path::push_numbers(self.out, &[dash.1]);
            self.out.extend_from_slice(b"d\n");
            self.state_mut().dash = dash;
        }
    }

    /// Sets the constant alphas, leaving the one given as `None` as it is.
    fn set_alpha(&mut self, fill: Option<f64>, stroke: Option<f64>) {
        let current = self.state().alpha;
        // Rounded as it will be written, so two alphas that write the same
        // number share one graphics state.
        let round = |v: f64| (v * 1000.0).round() / 1000.0;
        let wanted = (
            fill.map_or(current.0, round),
            stroke.map_or(current.1, round),
        );
        if wanted == current {
            return;
        }
        let at = match self
            .resources
            .alphas
            .iter()
            .position(|known| *known == wanted)
        {
            Some(at) => at,
            None => {
                self.resources.alphas.push(wanted);
                self.resources.alphas.len() - 1
            }
        };
        self.out
            .extend_from_slice(format!("/G{} gs\n", at + 1).as_bytes());
        self.state_mut().alpha = wanted;
    }

    fn set_text_state(&mut self, spacing: f64, stretch: f64, render: u8) {
        if self.state().spacing != spacing {
            path::push_numbers(self.out, &[spacing]);
            self.out.extend_from_slice(b"Tc\n");
            self.state_mut().spacing = spacing;
        }
        if self.state().stretch != stretch {
            path::push_numbers(self.out, &[stretch]);
            self.out.extend_from_slice(b"Tz\n");
            self.state_mut().stretch = stretch;
        }
        if self.state().render != render {
            self.out
                .extend_from_slice(format!("{render} Tr\n").as_bytes());
            self.state_mut().render = render;
        }
    }
}

/// A dash pattern, or nothing for a solid line: SVG draws a pattern with a
/// negative length, or one that adds up to nought, as solid, and so is one
/// in a unit this does not read, with a note.
fn dashes(list: Option<&str>, notes: &mut Notes) -> Vec<f64> {
    let Some(list) = list else {
        return Vec::new();
    };
    let mut values = Vec::new();
    for piece in list.split(|c: char| c == ',' || c.is_whitespace()) {
        if piece.is_empty() {
            continue;
        }
        match length(piece) {
            Some(value) if value >= 0.0 => values.push(value),
            Some(_) => return Vec::new(),
            None => {
                notes.add(format!(
                    "the stroke-dasharray {list:?} is not read, so the line is drawn solid"
                ));
                return Vec::new();
            }
        }
    }
    if values.iter().sum::<f64>() <= 0.0 {
        return Vec::new();
    }
    // An odd list is repeated to make it even.
    if values.len() % 2 == 1 {
        values.extend(values.clone());
    }
    values
}

/// `[x, y, width, height]` from a `viewBox`, when it is one SVG draws.
fn view_box(value: &str) -> Option<[f64; 4]> {
    let mut scanner = Scanner::new(value);
    let mut numbers = [0.0; 4];
    for number in &mut numbers {
        *number = scanner.number()?;
    }
    (numbers[2] > 0.0 && numbers[3] > 0.0).then_some(numbers)
}

/// The map that fits `view` into the viewport at (`x`, `y`) of `w` by `h`,
/// as `preserveAspectRatio`'s default does: as large as fits, centred.
fn fit(view: [f64; 4], x: f64, y: f64, w: f64, h: f64) -> Matrix {
    let [vx, vy, vw, vh] = view;
    let scale = (w / vw).min(h / vh);
    Matrix([
        scale,
        0.0,
        0.0,
        scale,
        x + (w - vw * scale) / 2.0 - vx * scale,
        y + (h - vh * scale) / 2.0 - vy * scale,
    ])
}
