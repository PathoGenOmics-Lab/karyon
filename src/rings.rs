//! A circular genome, drawn as concentric rings.
//!
//! A bacterial chromosome, a plasmid and an organelle genome have no left end
//! and no right end, and everything else in the crate stacks bands over a
//! horizontal [`Scale`](crate::Scale), which assumes both. A ring maps position
//! to an angle instead, which is a different coordinate system and therefore a
//! different container: [`Rings`] is to [`Ring`] what
//! [`Figure`](crate::Figure) is to [`Track`](crate::Track), and the two have
//! nothing in common but the [`Drawing`] trait, which is all it takes to put a
//! stack and a circle on one [`Panels`](crate::Panels) sheet.
//!
//! # The rings carry what the bands carry
//!
//! Annotation ([`FeatureRing`]), a quantity in windows ([`SignalRing`]), points
//! ([`MarkerRing`]) and a ruler ([`AxisRing`]), each of them saying how thick it
//! wants to be and then drawing between the two radii it is given. What the
//! circle adds is the middle: [`Rings::link`] draws a chord across it between
//! the two ends of an inversion, or a duplication and its source, and it is the
//! one element of the plot that is not a ring.
//!
//! # Which ring a thing goes on is a choice
//!
//! An angle is the same all the way across the plot but an arc is not, so the
//! outer rings have pixels to spare and the inner ones do not.
//! [`Polar::bp_per_px`] is how a ring finds out where it stands, and
//! [`Rings::push`] puts each new ring inside the last, so the first thing
//! pushed gets the most room to say something in. Push the ring with detail in
//! it first.
//!
//! # A circle is the one place a label cannot go beside its mark
//!
//! A band has a row to write a name on and a margin to the right of it. A ring
//! has neither: an arc a third of the way round has no horizontal beside it,
//! and [`FeatureRing::show_names`] is off by default precisely because a whole
//! annotation drawn with names is a wheel of unreadable text. So the arcs and
//! the chords carry a `<title>` instead, and the plot as a whole names itself
//! the way a [`Figure`](crate::Figure) does. An arc under a pixel wide gets
//! none: at four megabases on a ring of that radius a gene is a hairline held
//! open by [`FeatureRing::min_degrees`], and a tooltip on it would belong to
//! whichever of several overlapping slivers the pointer happened to catch.
//!
//! # Where zero is
//!
//! At twelve o'clock, running clockwise, which is the convention every circular
//! genome viewer uses. [`Rings::origin_gap`] opens a few degrees there by
//! default, so the plot shows the join rather than hiding it. Because the
//! sequence closes, a span may arrive with its end below its start: 900 to 100
//! on a thousand-base circle is the two hundred bases across the origin, not
//! the eight hundred the other way round.

use std::f64::consts::PI;
use std::fs;
use std::io;
use std::path::Path;

use crate::svg::{num, Anchor, SvgWriter};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::feature::{feature_title, span_label, strand_color};
use crate::track::window::Window;
use crate::track::{Feature, Strand};

/// The mapping from a position on a circular sequence to a point on the page.
///
/// Position zero is at twelve o'clock and coordinates run clockwise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polar {
    length: u64,
    cx: f64,
    cy: f64,
    sweep: f64,
    start: f64,
}

impl Polar {
    /// A mapping of a sequence `length` bases long about `(cx, cy)`.
    ///
    /// `gap` is how many degrees are left blank at twelve o'clock. A little of
    /// it makes the origin visible as a seam; none of it closes the circle.
    pub fn new(length: u64, cx: f64, cy: f64, gap: f64) -> Self {
        let gap = gap.clamp(0.0, 90.0).to_radians();
        Polar {
            length: length.max(1),
            cx,
            cy,
            sweep: 2.0 * PI - gap,
            start: gap / 2.0,
        }
    }

    /// Angle of a position, in radians clockwise from twelve o'clock.
    pub fn angle(&self, pos: u64) -> f64 {
        self.angle_at(pos as f64)
    }

    /// Angle of a fractional position.
    pub fn angle_at(&self, pos: f64) -> f64 {
        self.start + (pos / self.length as f64).clamp(0.0, 1.0) * self.sweep
    }

    /// The point at an angle and a radius.
    pub fn at(&self, angle: f64, radius: f64) -> (f64, f64) {
        (
            self.cx + radius * angle.sin(),
            self.cy - radius * angle.cos(),
        )
    }

    /// The point at a position and a radius.
    pub fn point(&self, pos: u64, radius: f64) -> (f64, f64) {
        self.at(self.angle(pos), radius)
    }

    /// Centre of the circle.
    pub fn center(&self) -> (f64, f64) {
        (self.cx, self.cy)
    }

    /// Length of the sequence.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// How many bases one pixel covers at `radius`.
    ///
    /// The answer depends on the radius, which is the thing a circular plot
    /// does that a linear one does not: the same feature is half as long on a
    /// ring of half the radius.
    pub fn bp_per_px(&self, radius: f64) -> f64 {
        let circumference = radius.abs() * self.sweep;
        if circumference <= 0.0 {
            f64::INFINITY
        } else {
            self.length as f64 / circumference
        }
    }

    /// A closed annular sector: the region between two radii over a span of
    /// sequence. This is the shape a feature or a window is drawn as.
    ///
    /// A span whose end is before its start wraps through the origin, which on
    /// a circular sequence is a real thing and not an error: a feature either
    /// side of coordinate zero runs the short way, through the top of the circle,
    /// and not the long way round the rest of the sequence. It comes back as
    /// one path of two subpaths, so it is still a single element.
    pub fn sector(&self, from: u64, to: u64, inner: f64, outer: f64) -> String {
        if to < from {
            return format!(
                "{}{}",
                self.wedge(from, self.length, inner, outer),
                self.wedge(0, to, inner, outer)
            );
        }
        self.wedge(from, to, inner, outer)
    }

    /// One annular sector, start before end.
    fn wedge(&self, from: u64, to: u64, inner: f64, outer: f64) -> String {
        let (a0, a1) = self.span(from, to);
        let (inner, outer) = (inner.min(outer).max(0.0), inner.max(outer));
        let (x0, y0) = self.at(a0, outer);
        let (x1, y1) = self.at(a1, outer);
        let (x2, y2) = self.at(a1, inner);
        let (x3, y3) = self.at(a0, inner);
        let large = usize::from(a1 - a0 > PI);
        format!(
            "M{} {}{}L{} {}{}Z",
            num(x0),
            num(y0),
            arc(outer, large, 1, x1, y1),
            num(x2),
            num(y2),
            arc(inner, large, 0, x3, y3),
        )
    }

    /// The outline of a whole circle at `radius`, as two half arcs.
    ///
    /// Two, because one arc cannot return to where it started: an SVG arc with
    /// the same start and end point draws nothing at all.
    pub fn circle(&self, radius: f64) -> String {
        let (top_x, top_y) = self.at(0.0, radius);
        let (bottom_x, bottom_y) = self.at(PI, radius);
        format!(
            "M{} {}{}{}",
            num(top_x),
            num(top_y),
            arc(radius, 0, 1, bottom_x, bottom_y),
            arc(radius, 0, 1, top_x, top_y),
        )
    }

    /// A ribbon across the middle joining two spans of sequence.
    ///
    /// Both ends are arcs at `radius` and the two sides are curves pulled
    /// towards the centre, which is what keeps a chord readable when a dozen of
    /// them cross: a straight one would be a chord of noise.
    pub fn ribbon(&self, from: (u64, u64), to: (u64, u64), radius: f64) -> String {
        let (a0, a1) = self.chord_span(from.0, from.1);
        let (b0, b1) = self.chord_span(to.0, to.1);
        let (ax0, ay0) = self.at(a0, radius);
        let (ax1, ay1) = self.at(a1, radius);
        let (bx0, by0) = self.at(b0, radius);
        let (bx1, by1) = self.at(b1, radius);
        format!(
            "M{} {}{}Q{} {} {} {}{}Q{} {} {} {}Z",
            num(ax0),
            num(ay0),
            arc(radius, usize::from(a1 - a0 > PI), 1, ax1, ay1),
            num(self.cx),
            num(self.cy),
            num(bx0),
            num(by0),
            arc(radius, usize::from(b1 - b0 > PI), 1, bx1, by1),
            num(self.cx),
            num(self.cy),
            num(ax0),
            num(ay0),
        )
    }

    /// The angles of one end of a chord, which unlike a sector is a single
    /// closed shape and so cannot be cut in two at the origin: a span that
    /// wraps is carried past twelve o'clock as an angle beyond the sweep.
    fn chord_span(&self, from: u64, to: u64) -> (f64, f64) {
        if to < from {
            let a0 = self.angle(from);
            return (a0, (self.angle(to) + 2.0 * PI).max(a0 + 1e-4));
        }
        self.span(from, to)
    }

    /// The angles of a span, always at least a hair wide so that a single base
    /// is still a visible mark rather than a line of zero length.
    fn span(&self, from: u64, to: u64) -> (f64, f64) {
        let a0 = self.angle(from);
        let a1 = self.angle(to.max(from)).max(a0 + 1e-4);
        (a0, a1.min(self.start + self.sweep))
    }
}

/// One SVG elliptical arc command.
fn arc(radius: f64, large: usize, sweep: usize, x: f64, y: f64) -> String {
    format!(
        "A{} {} 0 {} {} {} {}",
        num(radius),
        num(radius),
        large,
        sweep,
        num(x),
        num(y)
    )
}

/// Everything a ring needs in order to draw itself.
pub struct RingContext<'a> {
    /// Where to write the SVG elements.
    pub svg: &'a mut SvgWriter,
    /// The shared mapping from position to angle.
    pub polar: &'a Polar,
    /// Radius of the inside of this ring.
    pub inner: f64,
    /// Radius of the outside of this ring.
    pub outer: f64,
    /// Shared colours and fonts.
    pub theme: &'a Theme,
}

impl RingContext<'_> {
    /// Halfway between the two radii.
    pub fn middle(&self) -> f64 {
        (self.inner + self.outer) / 2.0
    }

    /// How thick the ring is.
    pub fn thickness(&self) -> f64 {
        self.outer - self.inner
    }
}

/// One concentric ring of a circular plot.
///
/// The parallel of [`Track`](crate::Track): say how thick you want to be, then
/// draw between the two radii you are given.
pub trait Ring {
    /// How much radius this ring wants, in pixels.
    fn thickness(&self) -> f64;

    /// Blank radius left inside this ring before the next one starts.
    fn gap(&self) -> f64 {
        5.0
    }

    /// Draws the ring between `ctx.inner` and `ctx.outer`.
    fn draw(&self, ctx: &mut RingContext<'_>);
}

/// A chord across the middle of the plot.
struct Chord {
    from: (u64, u64),
    to: (u64, u64),
    color: Option<String>,
    opacity: f64,
}

impl Chord {
    /// What a reader hovering one chord is told: the two spans it joins.
    ///
    /// Both of them, because a chord is the one mark here that is in two places
    /// at once and neither end means anything without the other. The spans are
    /// written the way the ruler writes a position, 1-based and inclusive.
    ///
    /// Each end is labelled, the way an alignment block labels its query and
    /// its target. Four numbers strung together with `to` used it once as a
    /// connector between the ends and twice as a connector inside them, so
    /// which pair belonged to which end had to be inferred from the shape of
    /// the sentence rather than read off it.
    fn title(&self) -> String {
        format!(
            "link, source {}, target {}",
            chord_end_label(self.from.0, self.from.1),
            chord_end_label(self.to.0, self.to.1)
        )
    }
}

/// One end of a chord in words, 1-based and inclusive like the ruler.
///
/// A span whose end is below its start runs through the origin, and
/// [`span_label`] reads it as a backwards one and pushes the end up to the
/// start, so 900 to 100 came out as "901 to 901". The two numbers are still
/// the right ones the right way round; what has to be said as well is that the
/// span gets from the first to the second through twelve o'clock.
fn chord_end_label(from: u64, to: u64) -> String {
    if to < from {
        return format!(
            "{} to {} across the origin",
            group_thousands(from + 1),
            group_thousands(to.max(1))
        );
    }
    span_label(from, to)
}

/// A circular sequence and the rings drawn around it.
///
/// ```
/// use karyon::{AxisRing, Feature, FeatureRing, Rings};
///
/// let svg = Rings::new(4_411_532)
///     .title("H37Rv")
///     .push(AxisRing::new())
///     .push(FeatureRing::new(vec![
///         Feature::new(759_807, 763_325).name("rpoB"),
///     ]))
///     .link((10_000, 40_000), (2_000_000, 2_030_000))
///     .to_svg();
///
/// assert!(svg.starts_with("<svg"));
/// assert!(svg.contains("H37Rv"));
/// ```
pub struct Rings {
    length: u64,
    rings: Vec<Box<dyn Ring>>,
    chords: Vec<Chord>,
    diameter: f64,
    margin: f64,
    gap_degrees: f64,
    theme: Theme,
    title: Option<String>,
    subtitle: Option<String>,
    description: Option<String>,
}

impl Rings {
    /// A plot of a circular sequence `length` bases long.
    pub fn new(length: u64) -> Self {
        Rings {
            length: length.max(1),
            rings: Vec::new(),
            chords: Vec::new(),
            diameter: 640.0,
            margin: 14.0,
            gap_degrees: 2.0,
            theme: Theme::light(),
            title: None,
            subtitle: None,
            description: None,
        }
    }

    /// Sets the diameter of the outermost ring in pixels.
    pub fn diameter(mut self, diameter: f64) -> Self {
        self.diameter = diameter.max(80.0);
        self
    }

    /// Sets the whitespace around the circle.
    pub fn margin(mut self, margin: f64) -> Self {
        self.margin = margin.max(0.0);
        self
    }

    /// Sets how many degrees are left blank at twelve o'clock.
    ///
    /// A couple of degrees makes the origin visible as a seam, which matters:
    /// a closed circle hides the fact that a coordinate system has to start
    /// somewhere and that the choice was arbitrary. Zero closes it.
    pub fn origin_gap(mut self, degrees: f64) -> Self {
        self.gap_degrees = degrees.clamp(0.0, 90.0);
        self
    }

    /// Replaces the theme.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme;
        self
    }

    /// Sets the name written in the middle of the circle.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Sets a second, quieter line under the title.
    pub fn subtitle(mut self, subtitle: impl Into<String>) -> Self {
        self.subtitle = Some(subtitle.into());
        self
    }

    /// Sets the `<desc>` of the rendered document: what the plot shows.
    ///
    /// This is the alt text, read in place of several thousand arcs by a screen
    /// reader and shown in place of the image when it does not load. Without it
    /// the document still names itself, but a name says which sequence this is
    /// and an alt text says what happens on it, and only the person drawing the
    /// plot knows that.
    ///
    /// ```
    /// use karyon::{AxisRing, Rings};
    ///
    /// let svg = Rings::new(4_411_532)
    ///     .description("GC skew turns over at the origin and again at the terminus.")
    ///     .push(AxisRing::new())
    ///     .to_svg();
    ///
    /// assert!(svg.contains("<desc"));
    /// assert!(svg.contains("GC skew turns over"));
    /// ```
    pub fn description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }

    /// Adds a ring inside the ones already there.
    pub fn push(mut self, ring: impl Ring + 'static) -> Self {
        self.rings.push(Box::new(ring));
        self
    }

    /// Adds a boxed ring, for building a plot at runtime.
    pub fn push_boxed(mut self, ring: Box<dyn Ring>) -> Self {
        self.rings.push(ring);
        self
    }

    /// Joins two spans of sequence with a ribbon across the middle.
    ///
    /// Both spans are `(start, end)`, 0-based and half-open. This is the one
    /// thing a circle can say that a stack of bands cannot: that two places far
    /// apart in coordinates belong together.
    pub fn link(self, from: (u64, u64), to: (u64, u64)) -> Self {
        self.link_colored(from, to, None, 0.35)
    }

    /// A ribbon in a colour and opacity of its own.
    pub fn link_colored(
        mut self,
        from: (u64, u64),
        to: (u64, u64),
        color: Option<String>,
        opacity: f64,
    ) -> Self {
        self.chords.push(Chord {
            from,
            to,
            color,
            opacity: opacity.clamp(0.0, 1.0),
        });
        self
    }

    /// Length of the sequence.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// How many rings the plot holds.
    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    /// Width and height of the rendered image, which is square.
    pub fn dimensions(&self) -> (f64, f64) {
        let side = self.diameter + self.margin * 2.0;
        (side, side)
    }

    /// Radius of the innermost edge of the last ring, where chords start.
    pub fn inner_radius(&self) -> f64 {
        let mut radius = self.diameter / 2.0;
        for ring in &self.rings {
            radius -= ring.thickness().max(0.0);
            radius -= ring.gap().max(0.0);
        }
        radius.max(4.0)
    }

    /// What the document calls itself: the text in the middle of the circle,
    /// or the length of the sequence.
    ///
    /// One of them is always there, so the `<title>` is never empty. A circle
    /// has no locus to fall back on the way a figure does: the plot is the whole
    /// sequence, so how long it is is the only thing left that identifies it.
    fn document_name(&self) -> String {
        match (&self.title, &self.subtitle) {
            (Some(title), Some(subtitle)) => format!("{title}, {subtitle}"),
            (Some(title), None) => title.clone(),
            (None, Some(subtitle)) => subtitle.clone(),
            (None, None) => format!(
                "A circular sequence of {} bases",
                group_thousands(self.length)
            ),
        }
    }

    /// The alt text: whatever [`Rings::description`] was given, or a statement
    /// of what the plot is made of.
    ///
    /// The fallback is built only from what the plot knows for certain, since
    /// a ring is not asked for a label the way a track is and there is nothing
    /// here to name the rings with. What they mean is what
    /// [`Rings::description`] exists for.
    fn document_description(&self) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        let rings = match self.rings.len() {
            0 => "no rings".to_string(),
            1 => "one ring".to_string(),
            n => format!("{n} rings"),
        };
        let chords = match self.chords.len() {
            0 => String::new(),
            1 => " and one chord across the middle".to_string(),
            n => format!(" and {n} chords across the middle"),
        };
        format!(
            "A karyon plot of a circular sequence {} bases long, with {}{}.",
            group_thousands(self.length),
            rings,
            chords
        )
    }

    /// Renders the plot to a standalone SVG document.
    pub fn to_svg(&self) -> String {
        self.to_svg_with_id_prefix("")
    }

    /// Renders the plot with every id it generates carrying `prefix`.
    ///
    /// Needed only when nesting it beside another drawing; see
    /// [`Figure::to_svg_with_id_prefix`](crate::Figure::to_svg_with_id_prefix).
    pub fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        let (width, height) = self.dimensions();
        let centre = width / 2.0;
        let polar = Polar::new(self.length, centre, centre, self.gap_degrees);
        let mut svg = SvgWriter::with_id_prefix(prefix);
        // A prefix means this plot is going inside another document, and a
        // nested document must not name itself: `<title>` resolves to the
        // innermost element under the pointer, so a title here would shadow
        // the one the sheet puts on the panel over the panel's whole area.
        // The same rule as
        // [`Figure::to_svg_with_id_prefix`](crate::Figure::to_svg_with_id_prefix).
        if prefix.is_empty() {
            svg.describe(&self.document_name(), &self.document_description());
        }

        // Chords first, so that the rings sit over them rather than being
        // washed out by a dozen translucent ribbons crossing the middle.
        let inner = self.inner_radius();
        for chord in &self.chords {
            let color = chord
                .color
                .clone()
                .unwrap_or_else(|| self.theme.accent.clone());
            // One ribbon per link, always wide enough to point at, since the
            // middle of the circle is the one part of the plot with nothing
            // else in it.
            svg.begin_titled(&chord.title());
            svg.path(
                &polar.ribbon(chord.from, chord.to, inner),
                &color,
                chord.opacity,
            );
            svg.end_group();
        }

        let mut outer = self.diameter / 2.0;
        for ring in &self.rings {
            let thickness = ring.thickness().max(0.0);
            let mut ctx = RingContext {
                svg: &mut svg,
                polar: &polar,
                inner: (outer - thickness).max(0.0),
                outer,
                theme: &self.theme,
            };
            ring.draw(&mut ctx);
            outer -= thickness + ring.gap().max(0.0);
        }

        if let Some(title) = &self.title {
            let baseline = if self.subtitle.is_some() {
                centre
            } else {
                centre + self.theme.title_font_size * 0.35
            };
            svg.text_bold(
                centre,
                baseline,
                title,
                &self.theme.foreground,
                self.theme.title_font_size,
                Anchor::Middle,
            );
        }
        if let Some(subtitle) = &self.subtitle {
            svg.text(
                centre,
                centre + self.theme.font_size + 4.0,
                subtitle,
                &self.theme.muted,
                self.theme.font_size,
                Anchor::Middle,
            );
        }

        svg.finish(
            width,
            height,
            &self.theme.background,
            &self.theme.font_family,
        )
    }

    /// Renders the plot and writes it to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save_svg(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, self.to_svg())
    }
}

/// A ruler of positions around the outside.
#[derive(Debug, Clone)]
pub struct AxisRing {
    thickness: f64,
    ticks: usize,
    show_labels: bool,
}

impl AxisRing {
    /// A ruler with ten ticks.
    pub fn new() -> Self {
        AxisRing {
            thickness: 22.0,
            ticks: 10,
            show_labels: true,
        }
    }

    /// Sets how much radius the ruler takes.
    ///
    /// The coordinates are written inside this band, so a ruler thinner than
    /// the text it would carry is drawn as ticks alone. See
    /// [`AxisRing::show_labels`].
    pub fn thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness.max(4.0);
        self
    }

    /// Sets roughly how many ticks to draw.
    pub fn ticks(mut self, ticks: usize) -> Self {
        self.ticks = ticks.max(2);
        self
    }

    /// Draws or hides the coordinate labels.
    ///
    /// Asking for them is not enough on a ruler too thin to hold them: the
    /// labels sit inside the tick marks, so on a thin band they would be
    /// printed over the ring below. A ruler with no room for legible text
    /// draws its ticks and says nothing.
    pub fn show_labels(mut self, show: bool) -> Self {
        self.show_labels = show;
        self
    }
}

impl Default for AxisRing {
    fn default() -> Self {
        AxisRing::new()
    }
}

impl Ring for AxisRing {
    fn thickness(&self) -> f64 {
        self.thickness
    }

    fn draw(&self, ctx: &mut RingContext<'_>) {
        let polar = ctx.polar;
        let step = nice_step(polar.length() as f64 / self.ticks as f64);
        // Labels go inside the tick marks rather than outside the circle, so
        // the plot stays inside the square it said it would occupy. The size
        // has to fit the band the ring was given as well: the ruler is drawn
        // from the outer edge inwards, so a full-size label on a thin ring
        // lands on the ring below, over whatever that one is drawing. Below
        // the size a label can be read at there is no label, which is the rule
        // the rest of the crate uses when a band is too thin to say something.
        let size = (ctx.theme.font_size - 1.0).min(ctx.thickness() - 12.0);
        let label_radius = ctx.outer - 7.0 - size * 0.7;
        let show_labels = self.show_labels && size >= 4.0;

        ctx.svg
            .path_stroked(&polar.circle(ctx.outer), &ctx.theme.rule, 1.0);

        let mut pos = 0u64;
        while pos < polar.length() {
            let angle = polar.angle(pos);
            let (x0, y0) = polar.at(angle, ctx.outer);
            let (x1, y1) = polar.at(angle, (ctx.outer - 5.0).max(ctx.inner));
            ctx.svg.line(x0, y0, x1, y1, &ctx.theme.rule, 1.0);

            if show_labels {
                let (tx, ty) = polar.at(angle, label_radius);
                ctx.svg.text(
                    tx,
                    ty + size * 0.35,
                    &megabases(pos),
                    &ctx.theme.muted,
                    size,
                    Anchor::Middle,
                );
            }
            pos = match pos.checked_add(step) {
                Some(next) => next,
                None => break,
            };
        }
    }
}

/// Annotation drawn as arcs, forward strand outside and reverse inside.
#[derive(Debug, Clone)]
pub struct FeatureRing {
    features: Vec<Feature>,
    thickness: f64,
    color: Option<String>,
    reverse_color: Option<String>,
    split_strands: bool,
    show_names: bool,
    min_degrees: f64,
}

impl FeatureRing {
    /// A ring of `features`.
    pub fn new(features: impl Into<Vec<Feature>>) -> Self {
        FeatureRing {
            features: features.into(),
            thickness: 16.0,
            color: None,
            reverse_color: None,
            split_strands: true,
            show_names: false,
            min_degrees: 0.12,
        }
    }

    /// Sets how much radius the ring takes.
    pub fn thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness.max(2.0);
        self
    }

    /// Sets the colours of the forward and reverse halves.
    pub fn colors(mut self, forward: impl Into<String>, reverse: impl Into<String>) -> Self {
        self.color = Some(forward.into());
        self.reverse_color = Some(reverse.into());
        self
    }

    /// Whether the two strands get half the ring each.
    ///
    /// On by default, which is how a circular sequence is normally drawn: the
    /// two strands rarely carry features at the same density, and one lane hides
    /// that.
    pub fn split_strands(mut self, split: bool) -> Self {
        self.split_strands = split;
        self
    }

    /// Draws or hides names beside the features.
    ///
    /// Only worth it for a handful of named loci. A whole annotation drawn with
    /// names is a wheel of unreadable text.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Sets the smallest angle a feature is drawn at, in degrees.
    ///
    /// A feature of a thousand bases on a four megabase sequence is a tenth of a
    /// degree. Without a floor the whole annotation would be invisible, and
    /// with one a feature's width says nothing about its length.
    pub fn min_degrees(mut self, degrees: f64) -> Self {
        self.min_degrees = degrees.max(0.0);
        self
    }

    /// The features.
    pub fn features(&self) -> &[Feature] {
        &self.features
    }
}

impl Ring for FeatureRing {
    fn thickness(&self) -> f64 {
        self.thickness
    }

    fn draw(&self, ctx: &mut RingContext<'_>) {
        let forward = self
            .color
            .clone()
            .unwrap_or_else(|| strand_color(Strand::Forward, ctx.theme).to_string());
        let reverse = self
            .reverse_color
            .clone()
            .unwrap_or_else(|| strand_color(Strand::Reverse, ctx.theme).to_string());
        let middle = ctx.middle();
        let floor =
            (self.min_degrees.to_radians() / (2.0 * PI) * ctx.polar.length() as f64).round() as u64;

        for feature in &self.features {
            let (lane_inner, lane_outer) = match (self.split_strands, feature.strand) {
                (true, Strand::Reverse) => (ctx.inner, middle - 0.5),
                (true, _) => (middle + 0.5, ctx.outer),
                (false, _) => (ctx.inner, ctx.outer),
            };
            let color = if feature.strand == Strand::Reverse {
                &reverse
            } else {
                &forward
            };
            let end = feature.end.max(feature.start.saturating_add(floor.max(1)));

            // One arc per feature and nothing is binned, so a feature is named
            // when it is something a pointer can land on. The measure is the
            // arc as drawn, floor included, and the width comes from the lane's
            // own radius: the same gene is half as wide on a ring of half the
            // radius, which is the whole difference between a circle and a
            // stack of bands.
            let lane = (lane_inner + lane_outer) / 2.0;
            let pixels = (end - feature.start) as f64 / ctx.polar.bp_per_px(lane);
            let title = if pixels >= 1.0 {
                feature_title(feature)
            } else {
                String::new()
            };
            let named = !title.is_empty();
            if named {
                ctx.svg.begin_titled(&title);
            }

            ctx.svg.path(
                &ctx.polar.sector(feature.start, end, lane_inner, lane_outer),
                color,
                1.0,
            );

            if self.show_names {
                if let Some(name) = &feature.name {
                    // The middle of the arc as drawn, floor included, which is
                    // also the only subtraction that cannot run backwards: a
                    // span that wraps through the origin has its end below its
                    // start, and `feature.end` here would take a u64 below
                    // zero.
                    let mid = feature.start + (end - feature.start) / 2;
                    let angle = ctx.polar.angle(mid);
                    let (x, y) = ctx.polar.at(angle, ctx.inner - 4.0);
                    // Anchored away from the centre, so a name on the left of
                    // the circle runs left and one on the right runs right.
                    let anchor = if angle.sin() < 0.0 {
                        Anchor::End
                    } else {
                        Anchor::Start
                    };
                    ctx.svg.text(
                        x,
                        y,
                        name,
                        &ctx.theme.muted,
                        ctx.theme.font_size - 1.0,
                        anchor,
                    );
                }
            }

            // The name goes inside the group with the arc, because a feature
            // drawn with one is two shapes and a tooltip on half of it is
            // worse than none.
            if named {
                ctx.svg.end_group();
            }
        }
    }
}

/// A quantity in windows, drawn as a ring either side of a baseline circle.
#[derive(Debug, Clone)]
pub struct SignalRing {
    windows: Vec<Window>,
    thickness: f64,
    baseline: f64,
    above_color: Option<String>,
    below_color: Option<String>,
    extent: Option<f64>,
    show_baseline: bool,
}

impl SignalRing {
    /// A ring over `windows`, with the baseline at zero.
    pub fn new(windows: impl Into<Vec<Window>>) -> Self {
        SignalRing {
            windows: windows.into(),
            thickness: 40.0,
            baseline: 0.0,
            above_color: None,
            below_color: None,
            extent: None,
            show_baseline: true,
        }
    }

    /// Sets how much radius the ring takes.
    pub fn thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness.max(4.0);
        self
    }

    /// Moves the circle the quantity is read against.
    pub fn baseline(mut self, baseline: f64) -> Self {
        self.baseline = baseline;
        self
    }

    /// Sets the colours for windows outside and inside the baseline.
    pub fn colors(mut self, above: impl Into<String>, below: impl Into<String>) -> Self {
        self.above_color = Some(above.into());
        self.below_color = Some(below.into());
        self
    }

    /// Pins how far the ring reaches either side of the baseline.
    pub fn extent(mut self, extent: f64) -> Self {
        self.extent = Some(extent.abs());
        self
    }

    /// Draws or hides the baseline circle.
    pub fn show_baseline(mut self, show: bool) -> Self {
        self.show_baseline = show;
        self
    }

    /// The windows.
    pub fn windows(&self) -> &[Window] {
        &self.windows
    }

    /// How far the ring reaches either side of the baseline.
    pub fn reach(&self) -> f64 {
        if let Some(extent) = self.extent {
            return extent.max(f64::MIN_POSITIVE);
        }
        let reach = self
            .windows
            .iter()
            .filter(|window| window.value.is_finite())
            .map(|window| (window.value - self.baseline).abs())
            .fold(0.0f64, f64::max);
        if reach > 0.0 {
            reach * 1.06
        } else {
            1.0
        }
    }
}

impl Ring for SignalRing {
    fn thickness(&self) -> f64 {
        self.thickness
    }

    fn draw(&self, ctx: &mut RingContext<'_>) {
        let middle = ctx.middle();
        let half = ctx.thickness() / 2.0;
        let reach = self.reach();

        if self.show_baseline {
            ctx.svg.path_stroked(
                &ctx.polar.circle(middle),
                &mix(ctx.theme.surface(), &ctx.theme.rule, 0.8),
                0.8,
            );
        }

        let above = self
            .above_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(0).to_string());
        let below = self
            .below_color
            .clone()
            .unwrap_or_else(|| ctx.theme.color(1).to_string());

        for window in &self.windows {
            if !window.value.is_finite() {
                continue;
            }
            let offset = ((window.value - self.baseline) / reach).clamp(-1.0, 1.0) * half;
            if offset.abs() < 1e-9 {
                continue;
            }
            let (inner, outer, color) = if offset > 0.0 {
                (middle, middle + offset, &above)
            } else {
                (middle + offset, middle, &below)
            };
            ctx.svg.path(
                &ctx.polar.sector(window.start, window.end, inner, outer),
                color,
                1.0,
            );
        }
    }
}

/// Points on the sequence, drawn as radial ticks.
#[derive(Debug, Clone)]
pub struct MarkerRing {
    positions: Vec<(u64, usize)>,
    thickness: f64,
    width: f64,
    colors: Vec<String>,
}

impl MarkerRing {
    /// A ring of positions, all in one colour.
    pub fn new(positions: impl IntoIterator<Item = u64>) -> Self {
        MarkerRing {
            positions: positions.into_iter().map(|pos| (pos, 0)).collect(),
            thickness: 10.0,
            width: 1.2,
            colors: Vec::new(),
        }
    }

    /// A ring of positions, each carrying an index into the palette.
    ///
    /// For anything with a class attached: a variant and its consequence, a
    /// resistance mutation and its drug, an insertion sequence and its family.
    pub fn categorised(positions: impl IntoIterator<Item = (u64, usize)>) -> Self {
        MarkerRing {
            positions: positions.into_iter().collect(),
            thickness: 10.0,
            width: 1.2,
            colors: Vec::new(),
        }
    }

    /// Sets how much radius the ring takes.
    pub fn thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness.max(2.0);
        self
    }

    /// Sets how wide one tick is drawn, in pixels.
    pub fn width(mut self, width: f64) -> Self {
        self.width = width.max(0.3);
        self
    }

    /// Overrides the palette the category indices point into.
    pub fn colors(mut self, colors: impl Into<Vec<String>>) -> Self {
        self.colors = colors.into();
        self
    }

    /// The positions and their categories.
    pub fn positions(&self) -> &[(u64, usize)] {
        &self.positions
    }
}

impl Ring for MarkerRing {
    fn thickness(&self) -> f64 {
        self.thickness
    }

    fn gap(&self) -> f64 {
        3.0
    }

    fn draw(&self, ctx: &mut RingContext<'_>) {
        for (pos, category) in &self.positions {
            let color = if self.colors.is_empty() {
                ctx.theme.color(*category).to_string()
            } else {
                self.colors[*category % self.colors.len()].clone()
            };
            let angle = ctx.polar.angle(*pos);
            let (x0, y0) = ctx.polar.at(angle, ctx.inner);
            let (x1, y1) = ctx.polar.at(angle, ctx.outer);
            ctx.svg.line(x0, y0, x1, y1, &color, self.width);
        }
    }
}

/// A position as megabases, or bases when the sequence is short.
fn megabases(pos: u64) -> String {
    if pos == 0 {
        return "0".to_string();
    }
    if pos >= 1_000_000 {
        let mb = pos as f64 / 1e6;
        return format!("{}{}", trim(mb), " Mb");
    }
    if pos >= 1_000 {
        return format!("{}{}", trim(pos as f64 / 1e3), " kb");
    }
    format!("{pos}")
}

fn trim(value: f64) -> String {
    let rounded = (value * 100.0).round() / 100.0;
    if rounded == rounded.trunc() {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded}")
    }
}

/// Rounds a raw tick interval to 1, 2 or 5 times a power of ten.
fn nice_step(raw: f64) -> u64 {
    if !raw.is_finite() || raw <= 1.0 {
        return 1;
    }
    let magnitude = 10f64.powf(raw.log10().floor());
    let normalised = raw / magnitude;
    let multiplier = if normalised <= 1.5 {
        1.0
    } else if normalised <= 3.0 {
        2.0
    } else if normalised <= 7.0 {
        5.0
    } else {
        10.0
    };
    ((multiplier * magnitude).round() as u64).max(1)
}

/// Something that renders to a standalone SVG and can go on a sheet.
///
/// Implemented by [`Figure`](crate::Figure) and by [`Rings`], which is what
/// lets a linear stack and a circular plot sit on one
/// [`Panels`](crate::Panels) sheet despite having nothing else in common.
pub trait Drawing {
    /// Width and height of the rendered image.
    fn dimensions(&self) -> (f64, f64);

    /// Renders it, with every generated id carrying `prefix`.
    fn to_svg_with_id_prefix(&self, prefix: &str) -> String;
}

impl Drawing for Rings {
    fn dimensions(&self) -> (f64, f64) {
        Rings::dimensions(self)
    }

    fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        Rings::to_svg_with_id_prefix(self, prefix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn polar() -> Polar {
        Polar::new(1_000, 100.0, 100.0, 0.0)
    }

    #[test]
    fn zero_is_at_twelve_oclock_and_coordinates_run_clockwise() {
        let p = polar();
        let (x, y) = p.point(0, 50.0);
        assert!((x - 100.0).abs() < 1e-9, "straight up from the centre");
        assert!((y - 50.0).abs() < 1e-9);

        // A quarter of the way round is three o'clock, to the right.
        let (x, y) = p.point(250, 50.0);
        assert!((x - 150.0).abs() < 1e-6, "clockwise, not anticlockwise");
        assert!((y - 100.0).abs() < 1e-6);
    }

    #[test]
    fn an_origin_gap_leaves_a_seam_at_the_top() {
        let closed = Polar::new(1_000, 0.0, 0.0, 0.0);
        let seamed = Polar::new(1_000, 0.0, 0.0, 20.0);
        assert!((closed.angle(0) - 0.0).abs() < 1e-12);
        assert!((seamed.angle(0) - 10f64.to_radians()).abs() < 1e-12);
        // The last base stops short of where the first one starts.
        let end = seamed.angle(1_000);
        assert!(end < 2.0 * PI - 10f64.to_radians() + 1e-12);
    }

    #[test]
    fn a_position_past_the_end_is_clamped_rather_than_wrapped() {
        // Wrapping would put a feature that overruns the origin at the top of
        // the circle, which is a different claim from the one the data made.
        let p = polar();
        assert_eq!(p.angle(5_000), p.angle(1_000));
    }

    #[test]
    fn resolution_falls_with_the_radius() {
        let p = polar();
        let outer = p.bp_per_px(100.0);
        let inner = p.bp_per_px(50.0);
        assert!((inner / outer - 2.0).abs() < 1e-9, "half the circumference");
    }

    #[test]
    fn a_sector_closes_and_a_circle_is_two_arcs() {
        let p = polar();
        let sector = p.sector(0, 100, 40.0, 60.0);
        assert!(sector.starts_with('M'));
        assert!(sector.ends_with('Z'));
        assert_eq!(sector.matches('A').count(), 2, "one arc per radius");

        let circle = p.circle(50.0);
        // Two, because an arc cannot return to the point it started from.
        assert_eq!(circle.matches('A').count(), 2);
    }

    #[test]
    fn a_span_of_one_base_is_still_a_visible_mark() {
        let p = Polar::new(4_411_532, 0.0, 0.0, 0.0);
        let sector = p.sector(100, 101, 40.0, 60.0);
        assert!(!sector.contains("NaN"));
        // Not a shape of zero width, which would draw nothing.
        assert!(p.span(100, 101).1 > p.span(100, 101).0);
    }

    #[test]
    fn a_span_across_the_origin_takes_the_short_way_round() {
        // On a circular sequence "from 900 to 100" is two hundred bases through
        // the top of the circle, not the eight hundred the other way. Read as a
        // backwards span it would come out as most of the chromosome.
        let p = polar();
        let wrapped = p.sector(900, 100, 40.0, 60.0);
        assert_eq!(wrapped.matches('M').count(), 2, "a piece either side of 0");
        assert_eq!(wrapped.matches('Z').count(), 2);
        assert_eq!(
            wrapped.matches("A60 60 0 1").count(),
            0,
            "neither piece is the long way round"
        );
        assert_eq!(p.sector(100, 300, 40.0, 60.0).matches('M').count(), 1);
    }

    #[test]
    fn a_ribbon_joins_two_spans_through_the_middle() {
        let p = polar();
        let ribbon = p.ribbon((0, 50), (500, 550), 40.0);
        assert_eq!(ribbon.matches('Q').count(), 2, "one curve per side");
        assert_eq!(ribbon.matches('A').count(), 2, "one arc per end");
        assert!(ribbon.ends_with('Z'));
        assert!(ribbon.contains("100 100"), "pulled towards the centre");
    }

    #[test]
    fn a_chord_end_across_the_origin_keeps_all_two_hundred_of_its_bases() {
        // 900 to 100 on a thousand-base circle is the two hundred bases
        // through twelve o'clock. A sector is cut in two at the origin before
        // it is measured, but a ribbon is one closed shape and cannot be, so
        // its ends went through `span`, which reads a wrapped pair as
        // backwards and pushes the end back up to the start: 0.02 bases of the
        // 200, a hairline where one end of the chord should be.
        let p = polar();
        let (a0, a1) = p.chord_span(900, 100);
        let covered = (a1 - a0) / (2.0 * PI) * 1_000.0;
        assert!((covered - 200.0).abs() < 1e-9, "{covered} bases, want 200");

        // Carried past twelve o'clock rather than wrapped, and the sine and
        // the cosine do not care, so the end still lands on position 100.
        let (x, y) = p.at(a1, 40.0);
        let (px, py) = p.point(100, 40.0);
        assert!((x - px).abs() < 1e-9 && (y - py).abs() < 1e-9);

        // Which leaves the large-arc flag right on its own: a fifth of the
        // circle is the short way round and four fifths is not.
        let short = p.ribbon((900, 100), (400, 600), 40.0);
        assert_eq!(short.matches("A40 40 0 1").count(), 0, "{short}");
        let long = p.ribbon((600, 400), (400, 600), 40.0);
        assert_eq!(long.matches("A40 40 0 1").count(), 1, "{long}");
    }

    #[test]
    fn a_chord_end_across_the_origin_says_that_it_crosses_it() {
        // The same wrapped span read as a backwards one came out of the shared
        // span label as "901 to 901", a hundred and ninety nine bases short of
        // what the ribbon covers.
        let svg = Rings::new(1_000).link((900, 100), (400, 600)).to_svg();
        assert!(
            svg.contains(
                "<title>link, source 901 to 100 across the origin, target 401 to 600</title>"
            ),
            "{svg}"
        );
    }

    #[test]
    fn rings_stack_inwards_and_leave_room_in_the_middle() {
        let plot = Rings::new(1_000)
            .diameter(400.0)
            .push(AxisRing::new().thickness(20.0))
            .push(FeatureRing::new(Vec::new()).thickness(10.0));
        // 200 outer, less 20 and its gap, less 10 and its gap.
        assert_eq!(plot.inner_radius(), 200.0 - 20.0 - 5.0 - 10.0 - 5.0);
        assert_eq!(plot.ring_count(), 2);
        assert_eq!(plot.dimensions(), (400.0 + 28.0, 400.0 + 28.0));
    }

    #[test]
    fn too_many_rings_still_leave_a_middle_to_draw_in() {
        let mut plot = Rings::new(1_000).diameter(120.0);
        for _ in 0..20 {
            plot = plot.push(FeatureRing::new(Vec::new()).thickness(20.0));
        }
        assert!(plot.inner_radius() > 0.0);
        assert!(!plot.to_svg().contains("NaN"));
    }

    #[test]
    fn an_empty_plot_is_a_valid_document() {
        let svg = Rings::new(4_411_532).to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(svg.ends_with("</svg>"));
    }

    #[test]
    fn the_ruler_labels_the_origin_and_walks_round() {
        let svg = Rings::new(4_411_532).push(AxisRing::new()).to_svg();
        assert!(svg.contains(">0</text>"));
        assert!(svg.contains("Mb</text>"), "{svg}");
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_two_strands_take_half_the_ring_each() {
        let features = vec![
            Feature::new(0, 100_000).strand(Strand::Forward),
            Feature::new(200_000, 300_000).strand(Strand::Reverse),
        ];
        let split = Rings::new(1_000_000)
            .push(FeatureRing::new(features.clone()))
            .to_svg();
        let merged = Rings::new(1_000_000)
            .push(FeatureRing::new(features).split_strands(false))
            .to_svg();
        // Split, the two lanes use different radii and therefore different
        // arc commands; merged they use the same ones.
        assert_ne!(split, merged);
        assert!(split.contains(Theme::light().color(0)));
        assert!(split.contains(Theme::light().color(1)));
    }

    #[test]
    fn a_gene_too_small_to_see_is_widened_to_the_floor() {
        // A thousand bases on a four megabase chromosome is a tenth of a
        // degree, which is nothing at all without a floor.
        let hair = Rings::new(4_411_532)
            .push(FeatureRing::new(vec![Feature::new(10_000, 11_000)]).min_degrees(0.0))
            .to_svg();
        let visible = Rings::new(4_411_532)
            .push(FeatureRing::new(vec![Feature::new(10_000, 11_000)]).min_degrees(1.0))
            .to_svg();
        assert!(visible.len() > hair.len());
    }

    #[test]
    fn a_named_feature_across_the_origin_is_drawn_rather_than_panicking() {
        // A gene either side of coordinate zero is ordinary on a circular
        // sequence. Its name was anchored at `start + (end - start) / 2`,
        // which for a wrapped span takes a u64 below zero: a panic where the
        // rule is that nothing panics on data.
        let wrapped = Feature {
            start: 4_400_000,
            end: 1_000,
            name: Some("dnaA".to_string()),
            strand: Strand::Forward,
            color: None,
        };
        let svg = Rings::new(4_411_532)
            .push(FeatureRing::new(vec![wrapped]).show_names(true))
            .to_svg();
        assert!(svg.contains(">dnaA</text>"), "{svg}");
        assert!(!svg.contains("NaN"));

        // And it lands under the arc as drawn, a little anticlockwise of
        // twelve o'clock on a plot whose centre is (334, 334), rather than at
        // whatever angle the wrapped arithmetic happened to give.
        let (_, x, y) = &labels(&svg)[0];
        assert!(
            (x - 324.178).abs() < 0.5 && (y - 34.161).abs() < 0.5,
            "{x} {y}"
        );
    }

    /// Every text element of a document as `(content, x, y)`.
    fn labels(svg: &str) -> Vec<(String, f64, f64)> {
        svg.match_indices("<text ")
            .map(|(i, _)| {
                let element = &svg[i..i + svg[i..].find("</text>").unwrap()];
                let content = element[element.find('>').unwrap() + 1..].to_string();
                let read = |name: &str| -> f64 {
                    let key = format!("{name}=\"");
                    let at = element.find(&key).unwrap() + key.len();
                    element[at..at + element[at..].find('"').unwrap()]
                        .parse()
                        .unwrap()
                };
                (content, read("x"), read("y"))
            })
            .collect()
    }

    #[test]
    fn a_ruler_keeps_its_coordinates_inside_the_band_it_was_given() {
        // The label radius was a fixed fourteen pixels in from the outer edge
        // whatever the thickness, so a four pixel ruler on a 660 pixel plot
        // printed all eight of its coordinates at radius 316, inside a band of
        // [326, 330] and on top of the gene ring at [301, 321] below it.
        let centre = (660.0 + 28.0) / 2.0;
        let genes: Vec<Feature> = (0..400)
            .map(|i| Feature::new(i * 10_000, i * 10_000 + 9_000))
            .collect();

        for (thickness, wanted) in [(4.0, 0), (10.0, 0), (20.0, 8), (22.0, 8)] {
            let outer = 330.0;
            let inner = outer - thickness;
            let svg = Rings::new(4_000_000)
                .diameter(660.0)
                .push(AxisRing::new().thickness(thickness))
                .push(FeatureRing::new(genes.clone()).thickness(20.0))
                .to_svg();

            let drawn = labels(&svg);
            assert_eq!(drawn.len(), wanted, "thickness {thickness}: {drawn:?}");
            for (content, x, y) in drawn {
                let radius = ((x - centre).powi(2) + (y - centre).powi(2)).sqrt();
                assert!(
                    radius >= inner && radius <= outer,
                    "{content} at {radius} is outside [{inner}, {outer}]"
                );
            }
        }
    }

    #[test]
    fn a_signal_ring_goes_outward_and_inward_from_its_baseline() {
        let windows = vec![
            Window::new(0, 100_000, 0.4),
            Window::new(100_000, 200_000, -0.3),
        ];
        let ring = SignalRing::new(windows).colors("#111111", "#222222");
        let svg = Rings::new(1_000_000).push(ring).to_svg();
        assert!(svg.contains("#111111"), "nothing outside the baseline");
        assert!(svg.contains("#222222"), "nothing inside it");
    }

    #[test]
    fn a_flat_signal_does_not_divide_by_zero() {
        let ring = SignalRing::new(vec![Window::new(0, 100, 0.0)]);
        assert!(ring.reach() > 0.0);
        assert!(!Rings::new(1_000).push(ring).to_svg().contains("NaN"));
    }

    #[test]
    fn markers_carry_their_categories_into_the_palette() {
        let theme = Theme::light();
        let svg = Rings::new(1_000_000)
            .push(MarkerRing::categorised(vec![(1_000, 0), (2_000, 1)]))
            .to_svg();
        assert!(svg.contains(theme.color(0)));
        assert!(svg.contains(theme.color(1)));
    }

    #[test]
    fn a_link_is_drawn_under_the_rings_rather_than_over_them() {
        let svg = Rings::new(1_000_000)
            .push(AxisRing::new())
            .link((0, 10_000), (500_000, 510_000))
            .to_svg();
        let ribbon = svg.find("<path").unwrap();
        let ruler = svg.find("stroke-width").unwrap();
        assert!(ribbon < ruler, "a dozen ribbons would wash the rings out");
    }

    #[test]
    fn the_title_sits_in_the_middle_of_the_circle() {
        let svg = Rings::new(4_411_532)
            .title("H37Rv")
            .subtitle("4.41 Mb")
            .to_svg();
        assert!(svg.contains(">H37Rv</text>"));
        assert!(svg.contains(">4.41 Mb</text>"));
    }

    #[test]
    fn the_document_names_itself_the_way_a_figure_does() {
        let svg = Rings::new(4_411_532)
            .title("H37Rv")
            .subtitle("4.41 Mb")
            .push(AxisRing::new())
            .to_svg();
        assert!(
            svg.contains(r#"<title id="karyon-title">H37Rv, 4.41 Mb</title>"#),
            "{svg}"
        );
        assert!(svg.contains(r#"role="img""#));
    }

    #[test]
    fn a_plot_with_no_title_falls_back_to_how_long_the_sequence_is() {
        // A circle has no locus to print: the plot is the whole sequence.
        let svg = Rings::new(4_411_532).to_svg();
        assert!(
            svg.contains(
                r#"<title id="karyon-title">A circular sequence of 4,411,532 bases</title>"#
            ),
            "{svg}"
        );
    }

    #[test]
    fn the_description_says_what_the_plot_is_made_of() {
        let svg = Rings::new(1_000_000)
            .push(AxisRing::new())
            .push(FeatureRing::new(Vec::new()))
            .link((0, 100), (500, 600))
            .to_svg();
        assert!(
            svg.contains(
                "A karyon plot of a circular sequence 1,000,000 bases long, with 2 rings and one chord across the middle."
            ),
            "{svg}"
        );
    }

    #[test]
    fn a_description_of_the_authors_own_replaces_the_inventory() {
        let svg = Rings::new(1_000)
            .description("Skew turns over at the terminus.")
            .to_svg();
        assert!(svg.contains("Skew turns over at the terminus."));
        assert!(!svg.contains("A karyon plot of a circular sequence"));
    }

    #[test]
    fn a_gene_wide_enough_to_point_at_carries_its_name_and_its_span() {
        let svg = Rings::new(1_000_000)
            .diameter(600.0)
            .push(FeatureRing::new(vec![Feature::new(100_000, 200_000)
                .name("rpoB")
                .strand(Strand::Forward)]))
            .to_svg();
        // 1-based and inclusive, the coordinates the ruler prints.
        assert!(
            svg.contains("<title>rpoB, 100,001 to 200,000, forward</title>"),
            "{svg}"
        );
    }

    #[test]
    fn a_nameless_gene_is_still_told_where_it_is() {
        let svg = Rings::new(1_000_000)
            .diameter(600.0)
            .push(FeatureRing::new(vec![
                Feature::new(400_000, 500_000).strand(Strand::Reverse)
            ]))
            .to_svg();
        assert!(
            svg.contains("<title>feature, 400,001 to 500,000, reverse</title>"),
            "{svg}"
        );
    }

    #[test]
    fn a_gene_thinner_than_a_pixel_is_not_named() {
        // A thousand bases on a four megabase chromosome is drawn as a hairline
        // held open by `min_degrees`. Naming it would name the floor rather
        // than the gene, and there is no pointing at it to read the name.
        let plot = || {
            Rings::new(4_411_532)
                .diameter(600.0)
                .push(FeatureRing::new(vec![
                    Feature::new(10_000, 11_000).name("x")
                ]))
        };
        assert!(!plot().to_svg().contains("<title>"), "{}", plot().to_svg());

        // The same gene on a sequence short enough for it to be an arc.
        let wide = Rings::new(20_000)
            .diameter(600.0)
            .push(FeatureRing::new(vec![
                Feature::new(10_000, 11_000).name("x")
            ]))
            .to_svg();
        assert!(
            wide.contains("<title>x, 10,001 to 11,000</title>"),
            "{wide}"
        );
    }

    #[test]
    fn a_chord_says_both_of_the_places_it_joins() {
        let svg = Rings::new(1_000_000)
            .link((100_000, 200_000), (600_000, 700_000))
            .to_svg();
        assert!(
            svg.contains(
                "<title>link, source 100,001 to 200,000, target 600,001 to 700,000</title>"
            ),
            "{svg}"
        );
    }

    #[test]
    fn every_group_a_tooltip_opens_is_closed_again() {
        let svg = Rings::new(1_000_000)
            .diameter(600.0)
            .push(AxisRing::new())
            .push(
                FeatureRing::new(vec![
                    Feature::new(0, 100_000).name("a").strand(Strand::Forward),
                    Feature::new(200_000, 300_000)
                        .name("b")
                        .strand(Strand::Reverse),
                    // Too thin to name, so it opens no group at all.
                    Feature::new(400_000, 400_100).name("c"),
                ])
                .show_names(true),
            )
            .push(MarkerRing::new([1_000, 2_000]))
            .link((0, 1_000), (500_000, 501_000))
            .to_svg();
        let open = svg.matches("<g ").count() + svg.matches("<g>").count();
        assert_eq!(open, svg.matches("</g>").count(), "{svg}");
        assert_eq!(svg.matches("<title>").count(), 3, "two arcs and one chord");
    }

    #[test]
    fn positions_are_written_in_the_unit_that_suits_them() {
        assert_eq!(megabases(0), "0");
        assert_eq!(megabases(500), "500");
        assert_eq!(megabases(2_500), "2.5 kb");
        assert_eq!(megabases(2_000_000), "2 Mb");
        assert_eq!(megabases(4_411_532), "4.41 Mb");
    }

    #[test]
    fn tick_steps_are_round_numbers() {
        assert_eq!(nice_step(441_153.2), 500_000);
        assert_eq!(nice_step(0.5), 1);
        assert_eq!(nice_step(f64::NAN), 1);
    }

    #[test]
    fn a_circle_and_a_stack_can_share_one_sheet() {
        use crate::figure::Figure;
        use crate::panels::Panels;
        use crate::region::Region;
        use crate::track::{AxisTrack, CoverageTrack};

        // The whole reason for the Drawing trait: the two have nothing in
        // common except that both render to an SVG of a known size.
        let stack = Figure::new(Region::new("chr1", 0, 1_000).unwrap())
            .push(CoverageTrack::new(0, vec![3.0; 1_000]))
            .push(AxisTrack::new());
        let circle = Rings::new(4_411_532).push(AxisRing::new()).title("H37Rv");

        let sheet = Panels::new().push(&stack, "A").push(&circle, "B").to_svg();
        assert_eq!(sheet.matches("<svg ").count(), 3);
        assert!(sheet.contains(">H37Rv</text>"));

        // The prefix is for ids the writer generates, which is what `url(#id)`
        // resolves against. Checking for the prefix on its own would pass on
        // any id at all, so it is checked only where there is one to prefix.
        if circle.to_svg().contains("karyon-clip-") {
            assert!(sheet.contains("p1-karyon-clip-"), "{sheet}");
        }

        // A nested document does not name itself: two `<title>` elements
        // inside one panel would leave the outer one unreachable everywhere
        // the inner one covers.
        assert_eq!(sheet.matches("<title id=").count(), 1, "{sheet}");
    }
}
