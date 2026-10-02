//! Geometry: SVG path data, the basic shapes, and the transforms that place
//! them, all brought down to the four things a PDF path is made of.
//!
//! A PDF path has a move, a straight line, a cubic Bezier and a close, and
//! nothing else. So a quadratic is raised to the cubic that draws the same
//! curve, a circle is four quarter arcs, and an elliptical arc is cut into
//! pieces of at most a right angle, each one cubic, by the endpoint to centre
//! conversion in appendix F.6 of SVG 1.1. A quarter circle drawn as one cubic
//! is off the true circle by 0.027 per cent of its radius at worst, a
//! hundredth of a pixel on a ring of forty: no reader can see it, and fewer
//! pieces would be.
//!
//! The whole grammar of path data is read, not only the absolute `M L A Q C Z`
//! the tracks write, because [`SvgWriter::path`](crate::SvgWriter::path) takes
//! a `d` from whoever calls it and a track outside the crate writes what it
//! likes. Path data that goes wrong part way is drawn up to where it went
//! wrong, which is the rule SVG gives a browser.

use crate::svg::num;

/// One piece of a path, in absolute coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Segment {
    Move(f64, f64),
    Line(f64, f64),
    Cubic(f64, f64, f64, f64, f64, f64),
    Close,
}

/// What an element outlines.
///
/// A plain rectangle stays a rectangle, because PDF writes one in a single
/// operator and more than a third of the elements in a typical figure are
/// rectangles: the background, every band, every block of the sequence.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Shape {
    Rect { x: f64, y: f64, w: f64, h: f64 },
    Path(Vec<Segment>),
}

/// The constant that makes a cubic Bezier follow a quarter circle:
/// 4/3 (sqrt 2 - 1).
const KAPPA: f64 = 0.552_284_749_830_793_4;

impl Shape {
    /// A rectangle, with its corners rounded by `rx` across and `ry` down.
    ///
    /// SVG's rules: a radius left out takes the other one's value, and each is
    /// at most half its side.
    pub(crate) fn rect(x: f64, y: f64, w: f64, h: f64, rx: Option<f64>, ry: Option<f64>) -> Shape {
        let (rx, ry) = match (rx, ry) {
            (Some(rx), Some(ry)) => (rx, ry),
            (Some(r), None) | (None, Some(r)) => (r, r),
            (None, None) => (0.0, 0.0),
        };
        let rx = rx.max(0.0).min(w / 2.0);
        let ry = ry.max(0.0).min(h / 2.0);
        if rx <= 0.0 || ry <= 0.0 {
            return Shape::Rect { x, y, w, h };
        }
        let (kx, ky) = (KAPPA * rx, KAPPA * ry);
        let (right, bottom) = (x + w, y + h);
        Shape::Path(vec![
            Segment::Move(x + rx, y),
            Segment::Line(right - rx, y),
            Segment::Cubic(right - rx + kx, y, right, y + ry - ky, right, y + ry),
            Segment::Line(right, bottom - ry),
            Segment::Cubic(
                right,
                bottom - ry + ky,
                right - rx + kx,
                bottom,
                right - rx,
                bottom,
            ),
            Segment::Line(x + rx, bottom),
            Segment::Cubic(x + rx - kx, bottom, x, bottom - ry + ky, x, bottom - ry),
            Segment::Line(x, y + ry),
            Segment::Cubic(x, y + ry - ky, x + rx - kx, y, x + rx, y),
            Segment::Close,
        ])
    }

    /// An ellipse about (`cx`, `cy`), as four quarter arcs.
    pub(crate) fn ellipse(cx: f64, cy: f64, rx: f64, ry: f64) -> Shape {
        let (kx, ky) = (KAPPA * rx, KAPPA * ry);
        Shape::Path(vec![
            Segment::Move(cx + rx, cy),
            Segment::Cubic(cx + rx, cy + ky, cx + kx, cy + ry, cx, cy + ry),
            Segment::Cubic(cx - kx, cy + ry, cx - rx, cy + ky, cx - rx, cy),
            Segment::Cubic(cx - rx, cy - ky, cx - kx, cy - ry, cx, cy - ry),
            Segment::Cubic(cx + kx, cy - ry, cx + rx, cy - ky, cx + rx, cy),
            Segment::Close,
        ])
    }

    /// The straight lines through a `points` list, closed for a polygon.
    ///
    /// An odd number of coordinates drops the last one, as SVG does.
    pub(crate) fn points(list: &str, close: bool) -> Shape {
        let mut scanner = Scanner::new(list);
        let mut segments = Vec::new();
        while let (Some(x), Some(y)) = (scanner.number(), scanner.number()) {
            segments.push(if segments.is_empty() {
                Segment::Move(x, y)
            } else {
                Segment::Line(x, y)
            });
        }
        if close && !segments.is_empty() {
            segments.push(Segment::Close);
        }
        Shape::Path(segments)
    }

    /// Whether there is anything to draw.
    pub(crate) fn is_empty(&self) -> bool {
        match self {
            Shape::Rect { w, h, .. } => !(*w > 0.0 && *h > 0.0),
            Shape::Path(segments) => segments.is_empty(),
        }
    }

    /// The same outline under `matrix`. A Bezier's control points move with
    /// the curve under any affine map, so the result is exact.
    pub(crate) fn transformed(&self, matrix: &Matrix) -> Shape {
        let segments = match self {
            Shape::Rect { x, y, w, h } => vec![
                Segment::Move(*x, *y),
                Segment::Line(x + w, *y),
                Segment::Line(x + w, y + h),
                Segment::Line(*x, y + h),
                Segment::Close,
            ],
            Shape::Path(segments) => segments.clone(),
        };
        let at = |x: f64, y: f64| matrix.apply(x, y);
        Shape::Path(
            segments
                .into_iter()
                .map(|segment| match segment {
                    Segment::Move(x, y) => {
                        let (x, y) = at(x, y);
                        Segment::Move(x, y)
                    }
                    Segment::Line(x, y) => {
                        let (x, y) = at(x, y);
                        Segment::Line(x, y)
                    }
                    Segment::Cubic(x1, y1, x2, y2, x, y) => {
                        let (x1, y1) = at(x1, y1);
                        let (x2, y2) = at(x2, y2);
                        let (x, y) = at(x, y);
                        Segment::Cubic(x1, y1, x2, y2, x, y)
                    }
                    Segment::Close => Segment::Close,
                })
                .collect(),
        )
    }

    /// The tight box around the outline, as `[left, top, right, bottom]`.
    ///
    /// Tight rather than the box of the control points, because that is the
    /// box SVG runs a gradient over: a curve's control points can stand well
    /// outside what it draws, and a fade stretched over them would start
    /// part way down its colour.
    pub(crate) fn bounds(&self) -> Option<[f64; 4]> {
        match self {
            Shape::Rect { x, y, w, h } => Some([*x, *y, x + w, y + h]),
            Shape::Path(segments) => {
                let mut bounds: Option<[f64; 4]> = None;
                let mut add = |x: f64, y: f64| {
                    let b = bounds.get_or_insert([x, y, x, y]);
                    b[0] = b[0].min(x);
                    b[1] = b[1].min(y);
                    b[2] = b[2].max(x);
                    b[3] = b[3].max(y);
                };
                let mut current = (0.0, 0.0);
                let mut start = (0.0, 0.0);
                for segment in segments {
                    match *segment {
                        Segment::Move(x, y) => {
                            add(x, y);
                            current = (x, y);
                            start = (x, y);
                        }
                        Segment::Line(x, y) => {
                            add(x, y);
                            current = (x, y);
                        }
                        Segment::Cubic(x1, y1, x2, y2, x, y) => {
                            add(x, y);
                            for t in extremes(current.0, x1, x2, x)
                                .into_iter()
                                .chain(extremes(current.1, y1, y2, y))
                                .flatten()
                            {
                                add(
                                    cubic_at(current.0, x1, x2, x, t),
                                    cubic_at(current.1, y1, y2, y, t),
                                );
                            }
                            current = (x, y);
                        }
                        Segment::Close => current = start,
                    }
                }
                bounds
            }
        }
    }

    /// Writes the outline as PDF path operators, in the coordinates it has.
    pub(crate) fn write(&self, out: &mut Vec<u8>) {
        match self {
            Shape::Rect { x, y, w, h } => {
                push_numbers(out, &[*x, *y, *w, *h]);
                out.extend_from_slice(b"re\n");
            }
            Shape::Path(segments) => {
                for segment in segments {
                    match *segment {
                        Segment::Move(x, y) => {
                            push_numbers(out, &[x, y]);
                            out.extend_from_slice(b"m\n");
                        }
                        Segment::Line(x, y) => {
                            push_numbers(out, &[x, y]);
                            out.extend_from_slice(b"l\n");
                        }
                        Segment::Cubic(x1, y1, x2, y2, x, y) => {
                            push_numbers(out, &[x1, y1, x2, y2, x, y]);
                            out.extend_from_slice(b"c\n");
                        }
                        Segment::Close => out.extend_from_slice(b"h\n"),
                    }
                }
            }
        }
    }
}

/// Each number followed by a space, in the form the SVG wrote it in.
pub(crate) fn push_numbers(out: &mut Vec<u8>, numbers: &[f64]) {
    for value in numbers {
        out.extend_from_slice(num(*value).as_bytes());
        out.push(b' ');
    }
}

/// Where in `0..1` a cubic's derivative on one axis is zero.
fn extremes(p0: f64, p1: f64, p2: f64, p3: f64) -> [Option<f64>; 2] {
    // The derivative is a quadratic a t^2 + b t + c.
    let a = 3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3);
    let b = 6.0 * (p0 - 2.0 * p1 + p2);
    let c = 3.0 * (p1 - p0);
    let inside = |t: f64| Some(t).filter(|t| *t > 0.0 && *t < 1.0);
    if a.abs() < 1e-12 {
        if b.abs() < 1e-12 {
            return [None, None];
        }
        return [inside(-c / b), None];
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return [None, None];
    }
    let root = discriminant.sqrt();
    [
        inside((-b + root) / (2.0 * a)),
        inside((-b - root) / (2.0 * a)),
    ]
}

fn cubic_at(p0: f64, p1: f64, p2: f64, p3: f64, t: f64) -> f64 {
    let u = 1.0 - t;
    u * u * u * p0 + 3.0 * u * u * t * p1 + 3.0 * u * t * t * p2 + t * t * t * p3
}

/// Reads SVG path data into absolute segments, stopping at the first thing
/// it cannot read.
pub(crate) fn parse(d: &str) -> Vec<Segment> {
    let mut scanner = Scanner::new(d);
    let mut out: Vec<Segment> = Vec::new();
    let (mut x, mut y) = (0.0f64, 0.0f64);
    let (mut start_x, mut start_y) = (0.0f64, 0.0f64);
    // The control point a smooth curve reflects, kept only while the segment
    // before it was the kind it reflects.
    let mut last_cubic: Option<(f64, f64)> = None;
    let mut last_quadratic: Option<(f64, f64)> = None;
    let mut command: Option<u8> = None;
    loop {
        scanner.skip_separators();
        let Some(next) = scanner.peek() else {
            break;
        };
        if next.is_ascii_alphabetic() {
            scanner.advance();
            if out.is_empty() && !matches!(next, b'M' | b'm') {
                // Path data has to open with a move, or it draws nothing.
                break;
            }
            command = Some(next);
            if matches!(next, b'Z' | b'z') {
                out.push(Segment::Close);
                x = start_x;
                y = start_y;
                last_cubic = None;
                last_quadratic = None;
                continue;
            }
        } else if matches!(command, None | Some(b'Z' | b'z')) {
            break;
        }
        let Some(letter) = command else {
            break;
        };
        let relative = letter.is_ascii_lowercase();
        let (ox, oy) = if relative { (x, y) } else { (0.0, 0.0) };
        let read = match letter.to_ascii_uppercase() {
            b'M' => scanner.numbers::<2>().map(|[px, py]| {
                x = px + ox;
                y = py + oy;
                start_x = x;
                start_y = y;
                out.push(Segment::Move(x, y));
                // The pairs after the first are lines.
                command = Some(if relative { b'l' } else { b'L' });
                last_cubic = None;
                last_quadratic = None;
            }),
            b'L' => scanner.numbers::<2>().map(|[px, py]| {
                x = px + ox;
                y = py + oy;
                out.push(Segment::Line(x, y));
                last_cubic = None;
                last_quadratic = None;
            }),
            b'H' => scanner.numbers::<1>().map(|[px]| {
                x = px + ox;
                out.push(Segment::Line(x, y));
                last_cubic = None;
                last_quadratic = None;
            }),
            b'V' => scanner.numbers::<1>().map(|[py]| {
                y = py + oy;
                out.push(Segment::Line(x, y));
                last_cubic = None;
                last_quadratic = None;
            }),
            b'C' => scanner.numbers::<6>().map(|[x1, y1, x2, y2, px, py]| {
                let (x2, y2) = (x2 + ox, y2 + oy);
                out.push(Segment::Cubic(x1 + ox, y1 + oy, x2, y2, px + ox, py + oy));
                x = px + ox;
                y = py + oy;
                last_cubic = Some((x2, y2));
                last_quadratic = None;
            }),
            b'S' => scanner.numbers::<4>().map(|[x2, y2, px, py]| {
                let (x1, y1) = match last_cubic {
                    Some((cx, cy)) => (2.0 * x - cx, 2.0 * y - cy),
                    None => (x, y),
                };
                let (x2, y2) = (x2 + ox, y2 + oy);
                out.push(Segment::Cubic(x1, y1, x2, y2, px + ox, py + oy));
                x = px + ox;
                y = py + oy;
                last_cubic = Some((x2, y2));
                last_quadratic = None;
            }),
            b'Q' => scanner.numbers::<4>().map(|[qx, qy, px, py]| {
                let (qx, qy) = (qx + ox, qy + oy);
                let (ex, ey) = (px + ox, py + oy);
                out.push(raise(x, y, qx, qy, ex, ey));
                x = ex;
                y = ey;
                last_quadratic = Some((qx, qy));
                last_cubic = None;
            }),
            b'T' => scanner.numbers::<2>().map(|[px, py]| {
                let (qx, qy) = match last_quadratic {
                    Some((cx, cy)) => (2.0 * x - cx, 2.0 * y - cy),
                    None => (x, y),
                };
                let (ex, ey) = (px + ox, py + oy);
                out.push(raise(x, y, qx, qy, ex, ey));
                x = ex;
                y = ey;
                last_quadratic = Some((qx, qy));
                last_cubic = None;
            }),
            b'A' => scanner
                .arc()
                .map(|(rx, ry, rotation, large, sweep, px, py)| {
                    let (ex, ey) = (px + ox, py + oy);
                    arc(&mut out, (x, y), (rx, ry), rotation, large, sweep, (ex, ey));
                    x = ex;
                    y = ey;
                    last_cubic = None;
                    last_quadratic = None;
                }),
            _ => None,
        };
        if read.is_none() {
            break;
        }
    }
    out
}

/// The cubic that draws the quadratic from (`x0`, `y0`) through the control
/// point (`qx`, `qy`) to (`x`, `y`).
fn raise(x0: f64, y0: f64, qx: f64, qy: f64, x: f64, y: f64) -> Segment {
    Segment::Cubic(
        x0 + 2.0 / 3.0 * (qx - x0),
        y0 + 2.0 / 3.0 * (qy - y0),
        x + 2.0 / 3.0 * (qx - x),
        y + 2.0 / 3.0 * (qy - y),
        x,
        y,
    )
}

/// An elliptical arc as cubics of at most a right angle each, by the
/// endpoint to centre conversion of SVG 1.1, appendix F.6.
fn arc(
    out: &mut Vec<Segment>,
    from: (f64, f64),
    radii: (f64, f64),
    rotation: f64,
    large: bool,
    sweep: bool,
    to: (f64, f64),
) {
    let ((x1, y1), (x2, y2)) = (from, to);
    if x1 == x2 && y1 == y2 {
        // An arc that ends where it starts is left out altogether.
        return;
    }
    let (mut rx, mut ry) = (radii.0.abs(), radii.1.abs());
    if rx == 0.0 || ry == 0.0 || !rx.is_finite() || !ry.is_finite() {
        out.push(Segment::Line(x2, y2));
        return;
    }
    let phi = rotation.to_radians();
    let (cos, sin) = (phi.cos(), phi.sin());
    let (dx, dy) = ((x1 - x2) / 2.0, (y1 - y2) / 2.0);
    let (x1p, y1p) = (cos * dx + sin * dy, -sin * dx + cos * dy);
    // Radii too small to reach are scaled up until they just do.
    let lambda = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
    if lambda > 1.0 {
        let scale = lambda.sqrt();
        rx *= scale;
        ry *= scale;
    }
    let numerator = rx * rx * ry * ry - rx * rx * y1p * y1p - ry * ry * x1p * x1p;
    let denominator = rx * rx * y1p * y1p + ry * ry * x1p * x1p;
    let mut coefficient = (numerator / denominator).max(0.0).sqrt();
    if large == sweep {
        coefficient = -coefficient;
    }
    let (cxp, cyp) = (coefficient * rx * y1p / ry, -coefficient * ry * x1p / rx);
    let cx = cos * cxp - sin * cyp + (x1 + x2) / 2.0;
    let cy = sin * cxp + cos * cyp + (y1 + y2) / 2.0;
    let angle = |ux: f64, uy: f64, vx: f64, vy: f64| (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
    let (ux, uy) = ((x1p - cxp) / rx, (y1p - cyp) / ry);
    let (vx, vy) = ((-x1p - cxp) / rx, (-y1p - cyp) / ry);
    let theta = angle(1.0, 0.0, ux, uy);
    let mut delta = angle(ux, uy, vx, vy);
    if !sweep && delta > 0.0 {
        delta -= std::f64::consts::TAU;
    } else if sweep && delta < 0.0 {
        delta += std::f64::consts::TAU;
    }
    let pieces = (delta.abs() / std::f64::consts::FRAC_PI_2 - 1e-9)
        .ceil()
        .max(1.0);
    let step = delta / pieces;
    let k = 4.0 / 3.0 * (step / 4.0).tan();
    let point = |t: f64| {
        (
            cx + rx * t.cos() * cos - ry * t.sin() * sin,
            cy + rx * t.cos() * sin + ry * t.sin() * cos,
        )
    };
    let tangent = |t: f64| {
        (
            -rx * t.sin() * cos - ry * t.cos() * sin,
            -rx * t.sin() * sin + ry * t.cos() * cos,
        )
    };
    let count = pieces as usize;
    for piece in 0..count {
        let a = theta + step * piece as f64;
        let b = a + step;
        let (ax, ay) = point(a);
        let (bx, by) = if piece + 1 == count { to } else { point(b) };
        let (dax, day) = tangent(a);
        let (dbx, dby) = tangent(b);
        out.push(Segment::Cubic(
            ax + k * dax,
            ay + k * day,
            bx - k * dbx,
            by - k * dby,
            bx,
            by,
        ));
    }
}

/// Reads numbers and flags out of path data or a points list.
pub(crate) struct Scanner<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Scanner<'a> {
    pub(crate) fn new(text: &'a str) -> Scanner<'a> {
        Scanner {
            bytes: text.as_bytes(),
            at: 0,
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.at).copied()
    }

    fn advance(&mut self) {
        self.at += 1;
    }

    fn skip_separators(&mut self) {
        while let Some(byte) = self.peek() {
            if byte.is_ascii_whitespace() || byte == b',' {
                self.advance();
            } else {
                break;
            }
        }
    }

    /// The next number, or `None` with nothing consumed.
    ///
    /// SVG's number grammar, in which `1.5.5` is two numbers and `-1-2` is
    /// two more: a second point or a sign ends the number before it.
    pub(crate) fn number(&mut self) -> Option<f64> {
        self.skip_separators();
        let begin = self.at;
        let mut end = self.at;
        let bytes = self.bytes;
        if matches!(bytes.get(end), Some(b'+' | b'-')) {
            end += 1;
        }
        let digits_from = end;
        while bytes.get(end).is_some_and(u8::is_ascii_digit) {
            end += 1;
        }
        let mut digits = end - digits_from;
        if bytes.get(end) == Some(&b'.') {
            end += 1;
            let fraction_from = end;
            while bytes.get(end).is_some_and(u8::is_ascii_digit) {
                end += 1;
            }
            digits += end - fraction_from;
        }
        if digits == 0 {
            return None;
        }
        if matches!(bytes.get(end), Some(b'e' | b'E')) {
            let mut exponent = end + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            if bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                    exponent += 1;
                }
                end = exponent;
            }
        }
        let text = std::str::from_utf8(&bytes[begin..end]).ok()?;
        let value = text.parse::<f64>().ok().filter(|v| v.is_finite())?;
        self.at = end;
        Some(value)
    }

    fn numbers<const N: usize>(&mut self) -> Option<[f64; N]> {
        let mut values = [0.0; N];
        for value in &mut values {
            *value = self.number()?;
        }
        Some(values)
    }

    /// One arc flag, which is a single `0` or `1` and may run straight on
    /// into the next number: `a5 5 0 015 5` has flags 0 and 1 and ends at 5,5.
    fn flag(&mut self) -> Option<bool> {
        self.skip_separators();
        let flag = match self.peek()? {
            b'0' => false,
            b'1' => true,
            _ => return None,
        };
        self.advance();
        Some(flag)
    }

    #[allow(clippy::type_complexity)]
    fn arc(&mut self) -> Option<(f64, f64, f64, bool, bool, f64, f64)> {
        let [rx, ry, rotation] = self.numbers::<3>()?;
        let large = self.flag()?;
        let sweep = self.flag()?;
        let [x, y] = self.numbers::<2>()?;
        Some((rx, ry, rotation, large, sweep, x, y))
    }
}

/// An affine map, `[a b c d e f]` as both SVG's `matrix()` and PDF's `cm`
/// write it: x' = a x + c y + e, y' = b x + d y + f.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Matrix(pub(crate) [f64; 6]);

impl Matrix {
    pub(crate) const IDENTITY: Matrix = Matrix([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);

    pub(crate) fn translate(x: f64, y: f64) -> Matrix {
        Matrix([1.0, 0.0, 0.0, 1.0, x, y])
    }

    pub(crate) fn scale(x: f64, y: f64) -> Matrix {
        Matrix([x, 0.0, 0.0, y, 0.0, 0.0])
    }

    /// This map applied after `inner`: the SVG list `self inner`.
    pub(crate) fn then_inner(&self, inner: &Matrix) -> Matrix {
        let [a, b, c, d, e, f] = self.0;
        let [p, q, r, s, t, u] = inner.0;
        Matrix([
            a * p + c * q,
            b * p + d * q,
            a * r + c * s,
            b * r + d * s,
            a * t + c * u + e,
            b * t + d * u + f,
        ])
    }

    pub(crate) fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let [a, b, c, d, e, f] = self.0;
        (a * x + c * y + e, b * x + d * y + f)
    }

    /// Writes the map as a `cm` operator.
    ///
    /// The four numbers that turn and scale get six decimals rather than the
    /// three a coordinate gets. A coordinate is rounded once; a rotation is
    /// multiplied into everything inside it, and at three decimals the cosine
    /// of a turned label's angle was off by enough to move its far end a
    /// twentieth of a pixel.
    pub(crate) fn write(&self, out: &mut Vec<u8>) {
        let [a, b, c, d, e, f] = self.0;
        for value in [a, b, c, d] {
            out.extend_from_slice(fine(value).as_bytes());
            out.push(b' ');
        }
        push_numbers(out, &[e, f]);
        out.extend_from_slice(b"cm\n");
    }
}

/// A number at six decimals, written as [`num`] writes one at three.
pub(crate) fn fine(value: f64) -> String {
    let scaled = (value * 1e6).round();
    if !scaled.is_finite() || scaled == 0.0 {
        return "0".to_string();
    }
    if scaled.abs() >= 1e15 {
        return num(value);
    }
    let rounded = scaled / 1e6;
    if rounded == rounded.trunc() {
        return format!("{}", rounded as i64);
    }
    let text = format!("{rounded:.6}");
    text.trim_end_matches('0').to_string()
}

/// Reads an SVG transform list, or `None` when any of it cannot be read, in
/// which case SVG draws the element as if it had none.
pub(crate) fn transform(list: &str) -> Option<Matrix> {
    let mut matrix = Matrix::IDENTITY;
    let mut rest = list.trim();
    while !rest.is_empty() {
        let open = rest.find('(')?;
        let close = rest.find(')')?;
        if close < open {
            return None;
        }
        let name = rest[..open].trim().trim_start_matches(',').trim();
        let mut scanner = Scanner::new(&rest[open + 1..close]);
        let mut values = Vec::with_capacity(6);
        while let Some(value) = scanner.number() {
            values.push(value);
        }
        scanner.skip_separators();
        if scanner.peek().is_some() {
            return None;
        }
        let step = match (name, values.as_slice()) {
            ("translate", [x]) => Matrix::translate(*x, 0.0),
            ("translate", [x, y]) => Matrix::translate(*x, *y),
            ("scale", [s]) => Matrix::scale(*s, *s),
            ("scale", [x, y]) => Matrix::scale(*x, *y),
            ("rotate", [degrees]) => rotation(*degrees),
            ("rotate", [degrees, cx, cy]) => Matrix::translate(*cx, *cy)
                .then_inner(&rotation(*degrees))
                .then_inner(&Matrix::translate(-cx, -cy)),
            ("skewX", [degrees]) => Matrix([1.0, 0.0, degrees.to_radians().tan(), 1.0, 0.0, 0.0]),
            ("skewY", [degrees]) => Matrix([1.0, degrees.to_radians().tan(), 0.0, 1.0, 0.0, 0.0]),
            ("matrix", [a, b, c, d, e, f]) => Matrix([*a, *b, *c, *d, *e, *f]),
            _ => return None,
        };
        matrix = matrix.then_inner(&step);
        rest = rest[close + 1..]
            .trim_start()
            .trim_start_matches(',')
            .trim_start();
    }
    Some(matrix)
}

fn rotation(degrees: f64) -> Matrix {
    let radians = degrees.to_radians();
    let (sin, cos) = radians.sin_cos();
    Matrix([cos, sin, -sin, cos, 0.0, 0.0])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn relative_path_commands_draw_what_absolute_ones_do() {
        let absolute = parse("M10 10 L20 10 H30 V20 C30 25 25 30 20 30 Q15 30 10 20 Z");
        let relative = parse("m10 10 l10 0 h10 v10 c0 5 -5 10 -10 10 q-5 0 -10 -10 z");
        assert_eq!(absolute.len(), relative.len());
        for (a, r) in absolute.iter().zip(&relative) {
            match (a, r) {
                (
                    Segment::Cubic(a1, a2, a3, a4, a5, a6),
                    Segment::Cubic(r1, r2, r3, r4, r5, r6),
                ) => {
                    for (x, y) in [(a1, r1), (a2, r2), (a3, r3), (a4, r4), (a5, r5), (a6, r6)] {
                        assert!(close(*x, *y), "{a:?} against {r:?}");
                    }
                }
                _ => assert_eq!(a, r),
            }
        }
    }

    #[test]
    fn arc_flags_written_without_separators_are_read() {
        let joined = parse("M0 0 A5 5 0 015 5");
        let spaced = parse("M0 0 A5 5 0 0 1 5 5");
        assert_eq!(joined, spaced);
        let Some(Segment::Cubic(.., x, y)) = joined.last() else {
            panic!("{joined:?}");
        };
        assert_eq!((*x, *y), (5.0, 5.0));
    }

    #[test]
    fn a_quarter_arc_is_one_cubic_within_a_thousandth_of_the_circle() {
        let segments = parse("M10 0 A10 10 0 0 1 0 10");
        assert_eq!(segments.len(), 2, "{segments:?}");
        let Segment::Cubic(x1, y1, x2, y2, x, y) = segments[1] else {
            panic!("{segments:?}");
        };
        // The circle is centred on the origin; every point of the cubic
        // should sit on it.
        for step in 0..=20 {
            let t = step as f64 / 20.0;
            let px = cubic_at(10.0, x1, x2, x, t);
            let py = cubic_at(0.0, y1, y2, y, t);
            let off = ((px * px + py * py).sqrt() - 10.0).abs() / 10.0;
            assert!(off < 0.001, "t {t}: off by {off}");
        }
        // A half turn is two pieces, a whole one drawn as two halves four.
        assert_eq!(parse("M10 0 A10 10 0 0 1 -10 0").len(), 3);
    }

    #[test]
    fn a_malformed_path_draws_up_to_the_error() {
        assert_eq!(
            parse("M0 0 L10 0 L10 x L0 10"),
            [Segment::Move(0.0, 0.0), Segment::Line(10.0, 0.0)]
        );
        // Path data has to start with a move.
        assert!(parse("L10 10").is_empty());
        // Numbers run together as SVG reads them.
        assert_eq!(
            parse("M1.5.5-2-3"),
            [Segment::Move(1.5, 0.5), Segment::Line(-2.0, -3.0)]
        );
    }

    #[test]
    fn the_box_of_a_curve_is_the_curve_s_and_not_its_control_points() {
        // A bump whose control points stand at y = -20 reaches only y = -15.
        let bump = Shape::Path(parse("M0 0 C0 -20 10 -20 10 0"));
        let [left, top, right, bottom] = bump.bounds().unwrap();
        assert!(close(left, 0.0) && close(right, 10.0) && close(bottom, 0.0));
        assert!(close(top, -15.0), "{top}");
    }

    #[test]
    fn a_transform_list_composes_left_to_right_as_svg_reads_it() {
        let m = transform("translate(10 20) rotate(-90)").unwrap();
        let (x, y) = m.apply(5.0, 0.0);
        assert!(close(x, 10.0) && close(y, 15.0), "{x} {y}");
        let m = transform("translate(3,4) scale(2)").unwrap();
        assert_eq!(m.apply(1.0, 1.0), (5.0, 6.0));
        assert!(transform("translate(1 2) wobble(3)").is_none());
        assert!(transform("scale(1 2 3)").is_none());
        assert_eq!(transform(""), Some(Matrix::IDENTITY));
    }

    #[test]
    fn a_rounded_corner_is_clamped_to_half_its_side() {
        let Shape::Path(segments) = Shape::rect(0.0, 0.0, 10.0, 4.0, Some(20.0), None) else {
            panic!("not rounded");
        };
        assert_eq!(segments[0], Segment::Move(5.0, 0.0));
        assert_eq!(
            Shape::rect(0.0, 0.0, 10.0, 4.0, Some(0.0), None),
            Shape::Rect {
                x: 0.0,
                y: 0.0,
                w: 10.0,
                h: 4.0
            }
        );
    }

    #[test]
    fn six_decimals_for_what_turns_and_scales() {
        assert_eq!(fine(0.123_456_789), "0.123457");
        assert_eq!(fine(6.123e-17), "0");
        assert_eq!(fine(-1.0), "-1");
        assert_eq!(fine(0.5), "0.5");
    }
}
