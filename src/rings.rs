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
//! # From the command line
//!
//! `karyon NC_000962.3 --circular genes.gff3 calls.vcf.gz depth.bedgraph` draws
//! the place, one whole sequence, as one of these plots, a ring a track in the
//! order written and a key under it naming each ring, through
//! [`build_circle`](crate::cli::stack::build_circle).
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

use crate::pdf::Pdf;
use crate::region::Region;
use crate::style::{Density, LinePattern, RenderProfile};
use crate::svg::{fit_text, num, text_rounded, text_width, Anchor, SvgWriter};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::coverage::{Aggregate, CoverageTrack};
use crate::track::feature::{feature_title, span_label, strand_color};
use crate::track::legend::Legend;
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
    ///
    /// An end narrower than [`CHORD_END`] is drawn that wide, centred on the
    /// span it stands for. The hair [`Polar::span`] holds a single base open
    /// with is 1e-4 of a radian, and a breakend join drawn at that width had
    /// ends 0.03 px across at a radius of 293: a ribbon nobody could see or
    /// point at. Only the drawing widens, so the chord's tooltip still says
    /// the base the file named.
    fn chord_span(&self, from: u64, to: u64) -> (f64, f64) {
        let (a0, a1) = if to < from {
            let a0 = self.angle(from);
            (a0, (self.angle(to) + 2.0 * PI).max(a0 + 1e-4))
        } else {
            self.span(from, to)
        };
        let floor = CHORD_END.to_radians();
        if a1 - a0 >= floor {
            return (a0, a1);
        }
        // Kept inside the sweep, so an end at the origin widens away from
        // the seam rather than into it.
        let middle = (a0 + a1) / 2.0;
        let lo = (middle - floor / 2.0)
            .max(self.start)
            .min(self.start + self.sweep - floor);
        (lo, lo + floor)
    }

    /// The angles of a span, always at least a hair wide so that a single base
    /// is still a visible mark rather than a line of zero length.
    fn span(&self, from: u64, to: u64) -> (f64, f64) {
        let a0 = self.angle(from);
        let a1 = self.angle(to.max(from)).max(a0 + 1e-4);
        (a0, a1.min(self.start + self.sweep))
    }
}

/// The narrowest a chord's end is drawn, in degrees.
///
/// About a pixel and a half at the radius chords leave from on a figure the
/// default size, which is wide enough to see and to put a pointer on, and
/// narrow enough that two breakends a few kilobases apart on a chromosome stay
/// two ends. A chord of 50 kb on a sequence of 4.4 Mb is four degrees, so
/// anything a reader drew as a span keeps its own width.
const CHORD_END: f64 = 0.3;

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
    /// Scale applied to ring thickness, gaps and fixed pixel measurements.
    pub visual_scale: f64,
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

    /// Scales one fixed pixel measurement for the active profile and density.
    pub fn px(&self, value: f64) -> f64 {
        value * self.visual_scale
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

    /// What the ring is called, where it is called anything.
    ///
    /// A band says its name in the gutter beside it, and a ring has nowhere
    /// to write one: a word at twelve o'clock runs over the arcs left of the
    /// seam, and anywhere else it runs over the ring's own data. So a named
    /// ring is one tooltip, its name, wherever a pointer lands on it and no
    /// mark of its own answers first, and [`Rings::key`] is where it is
    /// written out. Two rings of depth were otherwise two rings of arcs that
    /// nothing on the plot told apart.
    fn label(&self) -> Option<&str> {
        None
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

/// One line of the key under a circle, laid out.
struct KeyLine {
    /// Its top, below the foot of the square the circle is drawn in.
    top: f64,
    /// How tall it is: one row, or as many as its legend wraps into.
    height: f64,
    /// The ring's name, shortened to the column of names.
    name: String,
    /// The height of one row, which the name is centred on.
    row: f64,
    /// Where its legend starts, and how wide it may run.
    legend_x: f64,
    legend_width: f64,
}

/// One end of a chord in words, 1-based and inclusive like the ruler.
///
/// A span whose end is below its start runs through the origin, and
/// [`span_label`] reads it as a backwards one and pushes the end up to the
/// start, so 900 to 100 came out as "901 to 901". The two numbers are still
/// the right ones the right way round; what has to be said as well is that the
/// span gets from the first to the second through twelve o'clock.
fn chord_end_label(from: u64, to: u64) -> String {
    // One base, as a breakend is, is one position: `40,001 to 40,001` said
    // it twice.
    if to >= from && to - from <= 1 {
        return group_thousands(from + 1);
    }
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
    visual_scale: f64,
    density: Density,
    key: Vec<(String, Legend)>,
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
            visual_scale: 1.0,
            density: Density::Balanced,
            key: Vec::new(),
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

    /// Applies a named palette, type scale and ring density together.
    pub fn profile(mut self, profile: RenderProfile) -> Self {
        self.theme = if profile.is_dark() {
            Theme::dark()
        } else {
            Theme::light()
        };
        self.visual_scale = profile.visual_scale();
        self.density = profile.density();
        self
    }

    /// Scales typography, marks, margins and ring geometry together.
    pub fn visual_scale(mut self, factor: f64) -> Self {
        self.visual_scale = if factor.is_finite() {
            factor.max(0.25)
        } else {
            1.0
        };
        self
    }

    /// Sets how tightly the concentric data bands are packed.
    pub fn density(mut self, density: Density) -> Self {
        self.density = density;
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

    /// Adds a line to the key under the circle: a ring's name, then what its
    /// colours mean.
    ///
    /// The lines are written in the order they are added, one under the
    /// other, so added in the order the rings were pushed they read outside
    /// in, as the circle does. The colours repeat from ring to ring, since
    /// each ring takes the palette from its start, and what tells two rings
    /// apart is where they sit: a key that merged every colour into one list
    /// would say that one blue meant a forward strand, a missense call and a
    /// depth above the median all at once. A line's legend may be empty, and
    /// then the line is the name alone.
    ///
    /// The image grows downwards to hold the key, and a plot with none is
    /// the square it always was.
    ///
    /// ```
    /// use karyon::{FeatureRing, Feature, Legend, Rings};
    ///
    /// let plain = Rings::new(10_000);
    /// let keyed = Rings::new(10_000)
    ///     .push(FeatureRing::new(vec![Feature::new(100, 900)]).label("genes"))
    ///     .key("genes", Legend::new().key("forward strand", "#0072b2"));
    /// assert!(keyed.dimensions().1 > plain.dimensions().1);
    /// assert!(keyed.to_svg().contains("forward strand"));
    /// ```
    pub fn key(mut self, name: impl Into<String>, legend: Legend) -> Self {
        self.key.push((name.into(), legend));
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

    /// Width and height of the rendered image: a square, and the key under
    /// it where [`Rings::key`] was given one.
    pub fn dimensions(&self) -> (f64, f64) {
        let side = self.diameter + self.margin * self.visual_scale * 2.0;
        let theme = self.theme.clone().scaled(self.visual_scale);
        let key = self.key_lines(side, &theme);
        match key.last() {
            Some(last) => (
                side,
                side + last.top + last.height + self.margin * self.visual_scale,
            ),
            None => (side, side),
        }
    }

    /// Where each line of the key goes under the circle: its top, measured
    /// from the foot of the square, and how tall it is.
    ///
    /// Worked out in one place for [`Rings::dimensions`] and for the drawing,
    /// since a legend wraps at the width it is given and the image has to be
    /// as tall as the wrapped key, not as the key on one line.
    fn key_lines(&self, side: f64, theme: &Theme) -> Vec<KeyLine> {
        if self.key.is_empty() {
            return Vec::new();
        }
        let margin = self.margin * self.visual_scale;
        let font = theme.font_size;
        // The column of names is as wide as the widest of them, up to two
        // fifths of the plot, so a long file name is shortened rather than
        // pushing every legend off the right of the square.
        let widest = self
            .key
            .iter()
            .map(|(name, _)| text_width(name, font))
            .fold(0.0f64, f64::max);
        let column = widest.min((side - 2.0 * margin) * 0.4);
        let gap = theme.tokens.label_gap.max(8.0);
        let legend_x = margin + column + gap;
        let legend_width = (side - margin - legend_x).max(1.0);
        // The height of one row of a legend, which `Legend::height` gives a
        // legend with items in it: a name with nothing after it takes as
        // much, so the lines are evenly spaced whatever they hold.
        let row = Legend::new()
            .key("", "")
            .height(f64::INFINITY, theme)
            .max(font + 4.0);
        let mut top = 0.0;
        let mut lines = Vec::with_capacity(self.key.len());
        for (name, legend) in &self.key {
            let height = legend.height(legend_width, theme).max(row);
            lines.push(KeyLine {
                top,
                height,
                name: fit_text(name, column, font),
                row,
                legend_x,
                legend_width,
            });
            top += height;
        }
        lines
    }

    /// Radius of the innermost edge of the last ring, where chords start.
    pub fn inner_radius(&self) -> f64 {
        let scale = self.visual_scale * self.density.scale();
        let mut radius = self.diameter / 2.0;
        for ring in &self.rings {
            radius -= ring.thickness().max(0.0) * scale;
            radius -= ring.gap().max(0.0) * scale;
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
    /// The fallback is built only from what the plot knows for certain: how
    /// many rings, and what the ones given a [`Ring::label`] are called,
    /// outside in. What they mean is what [`Rings::description`] exists for.
    fn document_description(&self) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        let mut rings = match self.rings.len() {
            0 => "no rings".to_string(),
            1 => "one ring".to_string(),
            n => format!("{n} rings"),
        };
        // Said only where a ring has a name, so a plot of unnamed rings says
        // what it always said.
        let names: Vec<&str> = self
            .rings
            .iter()
            .filter_map(|ring| ring.label())
            .filter(|label| !label.is_empty())
            .collect();
        if let Some((last, rest)) = names.split_last() {
            let named = if rest.is_empty() {
                (*last).to_string()
            } else {
                format!("{} and {last}", rest.join(", "))
            };
            // A comma closes the list where a chord follows it, so the
            // last ring's name is not read as joined to the chord.
            let close = if self.chords.is_empty() { "" } else { "," };
            rings = format!("{rings}, outside in: {named}{close}");
        }
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
        let theme = self.theme.clone().scaled(self.visual_scale);
        let content_scale = self.visual_scale * self.density.scale();
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
            let color = chord.color.clone().unwrap_or_else(|| theme.accent.clone());
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
            let thickness = ring.thickness().max(0.0) * content_scale;
            // A named ring is one group under its name, so a pointer between
            // its marks, on a baseline or in a gap of the annotation, still
            // says which ring it is on. A mark with a tooltip of its own is a
            // group inside it and answers first.
            let label = ring.label().filter(|label| !label.is_empty());
            if let Some(label) = label {
                svg.begin_titled(label);
            }
            let mut ctx = RingContext {
                svg: &mut svg,
                polar: &polar,
                inner: (outer - thickness).max(0.0),
                outer,
                theme: &theme,
                visual_scale: content_scale,
            };
            ring.draw(&mut ctx);
            if label.is_some() {
                svg.end_group();
            }
            outer -= thickness + ring.gap().max(0.0) * content_scale;
        }

        // The key, under the square the circle is drawn in: each ring's name
        // in the colour of the text, and what its colours mean beside it.
        let margin = self.margin * self.visual_scale;
        for (line, (_, legend)) in self.key_lines(width, &theme).iter().zip(&self.key) {
            let top = width + line.top;
            if !line.name.is_empty() {
                svg.text(
                    margin,
                    top + line.row / 2.0 + theme.font_size * 0.35,
                    &line.name,
                    &theme.foreground,
                    theme.font_size,
                    Anchor::Start,
                );
            }
            legend.draw(&mut svg, line.legend_x, top, line.legend_width, &theme);
        }

        if let Some(title) = &self.title {
            let baseline = if self.subtitle.is_some() {
                centre
            } else {
                centre + theme.title_font_size * 0.35
            };
            let title = fit_text(title, self.inner_radius() * 1.65, theme.title_font_size);
            // A chord crosses the middle wherever it pleases, and the name of
            // the molecule is written in the middle, so the name carries a
            // halo of the page: the ribbon goes behind the letters instead of
            // through them.
            svg.text_haloed(
                centre,
                baseline,
                &title,
                &theme.foreground,
                theme.surface(),
                theme.title_font_size,
                Anchor::Middle,
                true,
            );
        }
        if let Some(subtitle) = &self.subtitle {
            let subtitle = fit_text(subtitle, self.inner_radius() * 1.65, theme.font_size);
            svg.text_haloed(
                centre,
                centre + theme.font_size + theme.tokens.row_gap,
                &subtitle,
                &theme.muted,
                theme.surface(),
                theme.font_size,
                Anchor::Middle,
                false,
            );
        }

        svg.finish(width, height, &theme.background, &theme.font_family)
    }

    /// Renders the plot and writes it to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save_svg(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, self.to_svg())
    }

    /// Renders the plot as a one-page PDF, converted from
    /// [`Rings::to_svg`]; see [`Pdf`] for what carries over.
    pub fn to_pdf(&self) -> Pdf {
        crate::pdf::drawn(&self.to_svg())
    }

    /// Renders the plot as PDF and writes it to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save_pdf(&self, path: impl AsRef<Path>) -> io::Result<()> {
        self.to_pdf().save(path)
    }
}

/// A ruler of positions around the outside.
#[derive(Debug, Clone)]
pub struct AxisRing {
    thickness: f64,
    ticks: usize,
    show_labels: bool,
    label: Option<String>,
}

impl AxisRing {
    /// A ruler with ten ticks.
    pub fn new() -> Self {
        AxisRing {
            thickness: 22.0,
            ticks: 10,
            show_labels: true,
            label: None,
        }
    }

    /// Names the ring, for its tooltip; see [`Ring::label`].
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
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

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
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
        let size = (ctx.theme.font_size - 1.0).min(ctx.thickness() - ctx.px(12.0));
        let label_radius = ctx.outer - ctx.px(7.0) - size * 0.7;
        let show_labels = self.show_labels && size >= 4.0;

        ctx.svg.path_stroked(
            &polar.circle(ctx.outer),
            &ctx.theme.rule,
            ctx.theme.tokens.stroke,
        );

        let mut pos = 0u64;
        while pos < polar.length() {
            let angle = polar.angle(pos);
            let (x0, y0) = polar.at(angle, ctx.outer);
            let (x1, y1) = polar.at(
                angle,
                (ctx.outer - ctx.theme.tokens.tick_length).max(ctx.inner),
            );
            ctx.svg
                .line(x0, y0, x1, y1, &ctx.theme.rule, ctx.theme.tokens.stroke);

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
    label: Option<String>,
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
            label: None,
        }
    }

    /// Names the ring, for its tooltip; see [`Ring::label`].
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets how much radius the ring takes.
    pub fn thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness.max(2.0);
        self
    }

    /// Sets the colours of the forward and reverse halves, for the features
    /// with no [`Feature::color`] of their own.
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

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
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
        let mut left_labels: Vec<(f64, f64)> = Vec::new();
        let mut right_labels: Vec<(f64, f64)> = Vec::new();
        let mut centre_labels: Vec<(f64, f64)> = Vec::new();

        for feature in &self.features {
            let (lane_inner, lane_outer) = match (self.split_strands, feature.strand) {
                (true, Strand::Reverse) => (ctx.inner, middle - 0.5),
                (true, _) => (middle + 0.5, ctx.outer),
                (false, _) => (ctx.inner, ctx.outer),
            };
            // A feature's own colour first, then the ring's, then the
            // strand's: the order a band of features takes them in, so one
            // annotation is painted alike drawn either way.
            let color = match &feature.color {
                Some(own) => own,
                None if feature.strand == Strand::Reverse => &reverse,
                None => &forward,
            };
            let end = if feature.end < feature.start {
                feature.end
            } else {
                feature.end.max(feature.start.saturating_add(floor.max(1)))
            };

            // One arc per feature and nothing is binned, so a feature is named
            // when it is something a pointer can land on. The measure is the
            // arc as drawn, floor included, and the width comes from the lane's
            // own radius: the same gene is half as wide on a ring of half the
            // radius, which is the whole difference between a circle and a
            // stack of bands.
            let lane = (lane_inner + lane_outer) / 2.0;
            let span = if end < feature.start {
                ctx.polar.length() - feature.start.min(ctx.polar.length()) + end
            } else {
                end - feature.start
            };
            let pixels = span as f64 / ctx.polar.bp_per_px(lane);
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
                    let mid = feature.start.saturating_add(span / 2) % ctx.polar.length();
                    let angle = ctx.polar.angle(mid);
                    let (x, y) = ctx.polar.at(angle, ctx.inner - ctx.px(4.0));
                    let (cx, _) = ctx.polar.center();
                    let font = ctx.theme.font_size - 1.0;
                    let bounds = (y - font * 0.8, y + font * 0.35);
                    // Beside the centre line a name runs towards it and has
                    // only the gap to the line to run in. Just past six
                    // o'clock that is a few pixels, and the name came out as
                    // a lone ellipsis, which names nothing: one with no room
                    // for a letter beside the line is centred on its arc.
                    let beside = (cx - x).abs() - ctx.theme.tokens.label_gap;
                    let lettered = !matches!(
                        fit_text(name, beside.max(0.0), font).as_str(),
                        "" | "\u{2026}"
                    );
                    let centred = (cx - x).abs() < ctx.theme.tokens.label_gap || !lettered;
                    let room = if centred {
                        ctx.inner * 1.5
                    } else {
                        (cx - x).abs() - ctx.theme.tokens.label_gap
                    };
                    let visible = fit_text(name, room.max(0.0), font);
                    let (anchor, occupied) = if centred {
                        (Anchor::Middle, &mut centre_labels)
                    } else if x < cx {
                        (Anchor::Start, &mut left_labels)
                    } else {
                        (Anchor::End, &mut right_labels)
                    };
                    let collides = occupied
                        .iter()
                        .any(|(top, bottom)| bounds.0 < *bottom && bounds.1 > *top);
                    if !visible.is_empty() && !collides {
                        ctx.svg.text(x, y, &visible, &ctx.theme.muted, font, anchor);
                        occupied.push(bounds);
                    }
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
    label: Option<String>,
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
            label: None,
        }
    }

    /// A ring of a per-base signal over a sequence `length` bases long,
    /// reduced to `bins` equal arcs, each the `aggregate` of the bases under
    /// it.
    ///
    /// What a depth file states, taken as [`CoverageTrack::from_spans`]
    /// takes it: half-open `(start, end, value)` spans, and a base no span
    /// covers at nought. A ring has no pixel columns to reduce a signal into
    /// the way a band does, and `samtools depth` over a chromosome is four
    /// million lines, which drawn as they stand were four million sectors.
    /// Held as runs and cut into arcs, it costs what its changes of value
    /// cost, and the ring draws at most `bins` sectors; neighbouring arcs of
    /// one value are drawn as one.
    ///
    /// ```
    /// use karyon::{Aggregate, SignalRing};
    ///
    /// // A thousand bases of depth 30 with one base of 300 among them.
    /// let spans = [(0, 1_000, 30.0), (500, 501, 300.0)];
    /// let max = SignalRing::from_spans(1_000, spans, 10, Aggregate::Max);
    /// let mean = SignalRing::from_spans(1_000, spans, 10, Aggregate::Mean);
    /// assert_eq!(max.windows().iter().map(|w| w.value).fold(0.0, f64::max), 300.0);
    /// assert_eq!(mean.windows().iter().map(|w| w.value).fold(0.0, f64::max), 32.7);
    /// ```
    pub fn from_spans(
        length: u64,
        spans: impl IntoIterator<Item = (u64, u64, f64)>,
        bins: usize,
        aggregate: Aggregate,
    ) -> Self {
        let whole = Region::new("ring", 0, length.max(1)).expect("a sequence of one base or more");
        let track = CoverageTrack::from_spans(&whole, spans).aggregate(aggregate);
        SignalRing::new(track.binned(bins))
    }

    /// Moves the baseline circle to the median of the windows, each counted
    /// for as many bases as it spans.
    ///
    /// The circle's counterpart of reading a depth against its usual level: a
    /// loss dips inside the line and a gain stands outside it, where against
    /// nought every window of a sequenced genome stands outside and a loss is
    /// only a shorter one. Windows with no value are left out, and a ring of
    /// none keeps its baseline.
    pub fn baseline_at_median(mut self) -> Self {
        let mut held: Vec<(f64, u64)> = self
            .windows
            .iter()
            .filter(|window| window.value.is_finite())
            .map(|window| (window.value, window.end.saturating_sub(window.start).max(1)))
            .collect();
        held.sort_by(|a, b| a.0.total_cmp(&b.0));
        let total: u64 = held.iter().map(|(_, bases)| bases).sum();
        let mut counted = 0u64;
        for (value, bases) in held {
            counted += bases;
            if counted.saturating_mul(2) >= total {
                self.baseline = value;
                break;
            }
        }
        self
    }

    /// Names the ring, for its tooltip; see [`Ring::label`].
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Where the baseline circle sits, in the units of the windows.
    pub fn baseline_value(&self) -> f64 {
        self.baseline
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

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut RingContext<'_>) {
        let middle = ctx.middle();
        let half = ctx.thickness() / 2.0;
        let reach = self.reach();

        if self.show_baseline {
            ctx.svg.path_stroked(
                &ctx.polar.circle(middle),
                &mix(ctx.theme.surface(), &ctx.theme.rule, 0.8),
                ctx.theme.tokens.hairline,
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
    label: Option<String>,
}

impl MarkerRing {
    /// A ring of positions, all in one colour.
    pub fn new(positions: impl IntoIterator<Item = u64>) -> Self {
        MarkerRing {
            positions: positions.into_iter().map(|pos| (pos, 0)).collect(),
            thickness: 10.0,
            width: 1.2,
            colors: Vec::new(),
            label: None,
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
            label: None,
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

    /// Names the ring, for its tooltip; see [`Ring::label`].
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
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

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
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
            let pattern = match category % 3 {
                1 => LinePattern::Dashed,
                2 => LinePattern::Dotted,
                _ => LinePattern::Solid,
            };
            ctx.svg.line_pattern(
                x0,
                y0,
                x1,
                y1,
                &color,
                self.width * ctx.visual_scale,
                pattern,
            );
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
        return format!("{}{}", text_rounded(mb, 2), " Mb");
    }
    if pos >= 1_000 {
        return format!("{}{}", text_rounded(pos as f64 / 1e3, 2), " kb");
    }
    format!("{pos}")
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
/// [`Panels`](crate::Panels) sheet despite having nothing else in common. The
/// maps implement it too, and so does [`Panels`](crate::Panels) itself, so a
/// sheet can be a panel of another sheet or one of several drawings inlined
/// into the same page.
pub trait Drawing {
    /// Width and height of the rendered image.
    fn dimensions(&self) -> (f64, f64);

    /// Renders it, with every generated id carrying `prefix`.
    ///
    /// The prefix goes in front of each id and each reference to one as it is
    /// given, neither escaped nor checked: a [`Panels`](crate::Panels) sheet
    /// renders a panel before it knows the prefix the panel's ids will need,
    /// and finds them again afterwards by the one it handed over.
    fn to_svg_with_id_prefix(&self, prefix: &str) -> String;

    /// Horizontal origin of the data area, when drawings can be aligned on a
    /// panel sheet. Circular and other free-form drawings return no anchor.
    fn content_anchor(&self) -> Option<f64> {
        None
    }

    /// The region a coordinate axis runs along underneath it, when one does.
    ///
    /// A figure whose tracks lay their marks out along its region answers
    /// with that region, and it is the window a viewer can pan and zoom.
    /// Everything else answers `None`: a circle maps position to an angle, a
    /// map and a sheet have no window of their own, and a stack of
    /// phylogenies holds one only because every figure is given one.
    fn region(&self) -> Option<&Region> {
        None
    }
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
            ..Feature::new(0, 1)
        };
        let svg = Rings::new(4_411_532)
            .push(FeatureRing::new(vec![wrapped]).show_names(true))
            .to_svg();
        assert!(svg.contains(">dnaA</text>"), "{svg}");
        assert!(!svg.contains("NaN"));

        // And it lands under the midpoint of the wrapped arc as drawn rather
        // than at whichever end happened not to underflow.
        let (_, x, y) = &labels(&svg)[0];
        let polar = Polar::new(4_411_532, 334.0, 334.0, 2.0);
        let span = 4_411_532 - 4_400_000 + 1_000;
        let midpoint = (4_400_000 + span / 2) % 4_411_532;
        let (expected_x, expected_y) = polar.point(midpoint, 300.0);
        assert!(
            (x - expected_x).abs() < 0.5 && (y - expected_y).abs() < 0.5,
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
    fn a_named_profile_scales_ring_geometry_and_typography_together() {
        let base = Rings::new(1_000).push(AxisRing::new());
        let large = Rings::new(1_000)
            .profile(RenderProfile::Presentation)
            .push(AxisRing::new());
        assert!(large.dimensions().0 > base.dimensions().0);
        assert!(large.inner_radius() < base.inner_radius());
        assert!(large.to_svg().contains(r#"font-size=""#));
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

    /// The prefix is only as good as its reach. An id it missed is an id two
    /// drawings in one page can both claim, and a reference it missed points
    /// at the other drawing's element, so every id each kind of drawing writes
    /// is checked, and everything that points at one: `url(#...)` for a clip,
    /// and `aria-labelledby` for the title and description.
    #[test]
    fn every_drawing_puts_its_prefix_on_every_id_it_writes() {
        use crate::figure::Figure;
        use crate::map::{GeoLocation, Map, PhyloMap};
        use crate::panels::Panels;
        use crate::region::Region;
        use crate::track::{AxisTrack, CoverageTrack};
        use crate::tree::Tree;

        let places = || {
            [
                GeoLocation::new("Peru", -9.19, -75.0152),
                GeoLocation::new("Spain", 40.4637, -3.7492),
            ]
        };
        let stack = Figure::new(Region::new("chr1", 0, 1_000).unwrap())
            .push(CoverageTrack::new(0, vec![3.0; 1_000]).label("depth"))
            .push(AxisTrack::new());
        let circle = Rings::new(4_411_532)
            .push(AxisRing::new())
            .push(SignalRing::new(vec![Window::new(0, 2_000_000, 0.4)]))
            .title("H37Rv");
        let map = Map::new().extend(places());
        let tree = Tree::parse_annotated_newick(
            "((A[&country=Peru]:1,B[&country=Peru]:1):1,C[&country=Spain]:2);",
        )
        .unwrap();
        let phylo = PhyloMap::new(tree)
            .location_by("country")
            .coordinates(places());
        let sheet = Panels::new().push(&stack, "A").push(&map, "B");

        let drawings: [(&str, &dyn Drawing); 5] = [
            ("a figure", &stack),
            ("a circle", &circle),
            ("a map", &map),
            ("a phylogeny on a map", &phylo),
            ("a sheet", &sheet),
        ];
        for (what, drawing) in drawings {
            let svg = drawing.to_svg_with_id_prefix("x-");
            let after = |key: &str, end: char| -> Vec<String> {
                svg.match_indices(key)
                    .map(|(at, found)| {
                        let rest = &svg[at + found.len()..];
                        rest[..rest.find(end).unwrap()].to_string()
                    })
                    .collect()
            };
            let written = after(r#" id=""#, '"');
            let pointed: Vec<String> = after("url(#", ')')
                .into_iter()
                .chain(
                    after(r#"aria-labelledby=""#, '"')
                        .iter()
                        .flat_map(|list| list.split(' ').map(str::to_string)),
                )
                .collect();
            // A circle may have nothing to clip, and nested it names nothing,
            // so it can pass with no ids at all. Everything else here clips.
            assert!(
                !written.is_empty() || what == "a circle",
                "{what} wrote no ids, so this checked nothing"
            );
            for id in written.iter().chain(&pointed) {
                assert!(
                    id.starts_with("x-"),
                    "{what} wrote {id:?} without the prefix"
                );
            }
            for id in &pointed {
                assert!(
                    written.contains(id),
                    "{what} points at {id:?}, which is not there"
                );
            }
            // The same drawing under two prefixes shares nothing, which is the
            // whole of what the prefix is for.
            let other = drawing.to_svg_with_id_prefix("y-");
            for id in &written {
                assert!(!other.contains(&format!(r#" id="{id}""#)), "{what}: {id}");
            }
        }
    }

    /// The angle between the two ends of the outer arc of the first end of
    /// the ribbon in `svg`, in degrees.
    fn first_chord_end(svg: &str) -> f64 {
        let path = svg
            .split("<title>link")
            .nth(1)
            .and_then(|rest| rest.split(" d=\"").nth(1))
            .and_then(|rest| rest.split('"').next())
            .expect("a ribbon");
        // `M x0 y0 A r r 0 large sweep x1 y1 Q ...`
        let numbers: Vec<f64> = path
            .split(|c: char| c.is_ascii_alphabetic() || c == ' ')
            .filter(|word| !word.is_empty())
            .take(9)
            .map(|word| word.parse().unwrap())
            .collect();
        let (x0, y0, x1, y1) = (numbers[0], numbers[1], numbers[7], numbers[8]);
        let radius = numbers[2];
        let chord = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
        (2.0 * (chord / (2.0 * radius)).asin()).to_degrees()
    }

    #[test]
    fn a_chord_between_two_points_is_wide_enough_to_see() {
        // Two breakends, one base each, on a megabase: 0.03 px across at the
        // radius chords leave from, before the floor.
        let svg = Rings::new(1_000_000)
            .link((100_000, 100_000), (600_000, 600_000))
            .to_svg();
        let wide = first_chord_end(&svg);
        assert!((wide - CHORD_END).abs() < 0.01, "{wide} degrees");
        // The tooltip still says the base the join is at, once.
        assert!(
            svg.contains("<title>link, source 100,001, target 600,001</title>"),
            "{svg}"
        );
        // And an end wider than the floor keeps its own width.
        let svg = Rings::new(1_000_000)
            .origin_gap(0.0)
            .link((100_000, 110_000), (600_000, 610_000))
            .to_svg();
        let own = first_chord_end(&svg);
        assert!((own - 3.6).abs() < 0.01, "{own} degrees");
        // A point at the origin widens away from the seam, not into it.
        let at_origin = Rings::new(1_000_000)
            .link((0, 0), (600_000, 600_000))
            .to_svg();
        let polar = Polar::new(1_000_000, 0.0, 0.0, 2.0);
        let (a0, a1) = polar.chord_span(0, 0);
        assert!(a0 >= 1f64.to_radians() - 1e-12, "{a0}");
        assert!((a1 - a0 - CHORD_END.to_radians()).abs() < 1e-12);
        assert!(at_origin.contains("source 1, target 600,001"));
    }

    #[test]
    fn a_labelled_ring_is_one_titled_group_and_an_unlabelled_one_renders_as_before() {
        let ring = || {
            SignalRing::new(vec![
                Window::new(0, 500, 2.0),
                Window::new(500, 1_000, -1.0),
            ])
        };
        let plain = Rings::new(1_000).push(ring()).to_svg();
        let named = Rings::new(1_000).push(ring().label("depth")).to_svg();
        assert!(!plain.contains("<title>depth</title>"));
        assert!(named.contains("<g><title>depth</title><path"), "{named}");
        // The same marks either way, inside the group.
        let marks = |svg: &str| svg.matches("<path").count();
        assert_eq!(marks(&plain), marks(&named));
        // An empty name is no name.
        assert_eq!(Rings::new(1_000).push(ring().label("")).to_svg(), plain);
        // And every ring takes one.
        let all = Rings::new(1_000)
            .push(AxisRing::new().label("ruler"))
            .push(FeatureRing::new(vec![Feature::new(10, 400)]).label("genes"))
            .push(MarkerRing::new([100, 200]).label("calls"))
            .to_svg();
        for name in ["ruler", "genes", "calls"] {
            assert!(all.contains(&format!("<g><title>{name}</title>")), "{name}");
        }
        assert!(all.contains("outside in: ruler, genes and calls."), "{all}");
    }

    #[test]
    fn signal_ring_from_spans_aggregates_each_bin_with_max_mean_and_min() {
        // Ten bases of 10 with one of 100 at base 3, and nothing past 20.
        let spans = [(0, 20, 10.0), (3, 4, 100.0)];
        let values = |aggregate| -> Vec<(u64, u64, f64)> {
            SignalRing::from_spans(40, spans, 4, aggregate)
                .windows()
                .iter()
                .map(|window| (window.start, window.end, window.value))
                .collect()
        };
        // Neighbouring arcs of one value are one window.
        assert_eq!(
            values(Aggregate::Max),
            [(0, 10, 100.0), (10, 20, 10.0), (20, 40, 0.0)]
        );
        assert_eq!(
            values(Aggregate::Mean),
            [(0, 10, 19.0), (10, 20, 10.0), (20, 40, 0.0)]
        );
        assert_eq!(values(Aggregate::Min), [(0, 20, 10.0), (20, 40, 0.0)]);
        // Never more arcs than bases, and never more than asked for.
        assert_eq!(
            SignalRing::from_spans(3, [(0, 3, 1.0)], 1_000, Aggregate::Max)
                .windows()
                .len(),
            1
        );
        let jagged: Vec<(u64, u64, f64)> = (0..10_000)
            .map(|at| (at, at + 1, (at % 5) as f64))
            .collect();
        assert_eq!(
            SignalRing::from_spans(10_000, jagged, 100, Aggregate::Max).windows(),
            [Window::new(0, 10_000, 4.0)]
        );
    }

    #[test]
    fn the_baseline_at_the_median_counts_each_window_for_its_bases() {
        // Most of the sequence at 30, a short stretch at 0 and two at 60:
        // by window the middle one of five is 30 either way, and by base too;
        // with the short windows the many, by window alone it would be 60.
        let ring = SignalRing::new(vec![
            Window::new(0, 900, 30.0),
            Window::new(900, 910, 60.0),
            Window::new(910, 920, 60.0),
            Window::new(920, 930, 60.0),
            Window::new(930, 1_000, f64::NAN),
        ])
        .baseline_at_median();
        assert_eq!(ring.baseline_value(), 30.0);
        // No window with a value keeps the baseline where it was.
        let none = SignalRing::new(vec![Window::new(0, 10, f64::NAN)])
            .baseline(2.0)
            .baseline_at_median();
        assert_eq!(none.baseline_value(), 2.0);
    }

    #[test]
    fn a_feature_ring_paints_a_feature_its_own_colour() {
        let svg = Rings::new(1_000)
            .push(
                FeatureRing::new(vec![
                    Feature::new(10, 400).name("own").color("#123456"),
                    Feature::new(500, 900).name("plain"),
                ])
                .colors("#aaaaaa", "#bbbbbb"),
            )
            .to_svg();
        let fill = |name: &str| -> String {
            let at = svg.find(&format!("<title>{name}, ")).unwrap();
            svg[at..]
                .split("fill=\"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap()
                .to_string()
        };
        assert_eq!(fill("own"), "#123456");
        assert_eq!(fill("plain"), "#aaaaaa");
    }

    #[test]
    fn a_name_with_no_room_beside_the_centre_line_is_centred_on_its_arc() {
        // Just past six o'clock, a few pixels from the centre line: beside
        // it the name had room for an ellipsis and nothing else.
        let svg = Rings::new(1_000_000)
            .push(
                FeatureRing::new(vec![Feature::new(503_000, 506_000).name("katG")])
                    .show_names(true),
            )
            .to_svg();
        assert!(!svg.contains(">\u{2026}</text>"), "{svg}");
        assert!(svg.contains("text-anchor=\"middle\">katG</text>"), "{svg}");
    }

    #[test]
    fn the_key_under_a_circle_is_a_line_a_ring_in_the_order_given() {
        let plain = Rings::new(1_000).push(AxisRing::new());
        let (side, height) = plain.dimensions();
        assert_eq!(side, height, "a circle with no key is a square");
        let keyed = Rings::new(1_000)
            .push(AxisRing::new())
            .key(
                "outer",
                Legend::new().key("one", "#111111").key("two", "#222222"),
            )
            .key("inner", Legend::new());
        let (width, tall) = keyed.dimensions();
        assert_eq!(width, side);
        assert!(tall > side + 30.0, "{tall}");
        let svg = keyed.to_svg();
        assert!(svg.contains(&format!("height=\"{}\"", num(tall))), "{svg}");
        let y = |text: &str| -> f64 {
            let end = svg.find(&format!(">{text}</text>")).unwrap();
            let start = svg[..end].rfind("<text").unwrap();
            svg[start..end]
                .split(" y=\"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        assert!(y("outer") > side && y("outer") < y("inner"));
        assert!(y("one") > side && y("one") < y("inner"));
        assert!(
            (y("outer") - y("one")).abs() < 1.0,
            "a name sits on its legend's row"
        );
        // A legend too wide for one row wraps, and the image is as tall as
        // the wrapped key.
        let long: Legend = (0..30).fold(Legend::new(), |legend, at| {
            legend.key(format!("consequence number {at}"), "#333333")
        });
        let wrapped = Rings::new(1_000).key("calls", long);
        assert!(wrapped.dimensions().1 > keyed.dimensions().1 + 40.0);
    }
}
