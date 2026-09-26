//! The figure: a stack of tracks over one shared coordinate axis.
//!
//! Layout is the whole of this module: one pass that turns a region, a width
//! and a list of tracks into rectangles, and it is kept away from the tracks on
//! purpose. A track is handed a band and the shared [`Scale`], and is never
//! told which band it is, how many others there are, or what any of them asked
//! for, so a track type has nothing to negotiate with and no way to disturb its
//! neighbours.
//!
//! # Why one track's value axis moves every other track
//!
//! Two strips come out of the width before anything is drawn, and each is
//! settled by asking every track and taking one answer for all of them. The
//! label gutter is there when any track has a [`label`](Track::label) to put in
//! it. The value axis is the widest [`y_axis_width`](Track::y_axis_width) any
//! one track asks for, and that width is then taken from all of them: a depth
//! profile with room for tick labels would otherwise begin a tick label's width
//! to the right of the ruler beneath it, and two tracks whose x axes disagree
//! are worse than either alone, because nothing in the picture says they do.
//!
//! # Heights are computed, and computed last
//!
//! The horizontal decisions come first because the heights depend on them. The
//! two strips are reserved, what is left of the width becomes the [`Scale`],
//! and only then is each track asked how tall it wants to be, since some answer
//! with the scale in hand: a pileup packs the reads that are in view, so it is
//! a different height at a different zoom. Nothing after that may change the
//! width, which is why the width is a setting on the figure that no track can
//! influence, while the height is not a setting at all: it is whatever the
//! tracks came to.
//!
//! # Nothing a track draws can leave its band
//!
//! [`Scale`] does not clamp, so a feature beginning before the window has a
//! negative x and a read running past the end has one off the right edge. Each
//! track is drawn inside a clip over its band and the axis strip it asked for
//! itself, which is what makes that overhang free: a track draws the whole of a
//! partly visible thing and the clip decides how much shows. A track that asked
//! for no axis is clipped to its band alone, so it cannot paint left of the
//! plot origin however wide a neighbour's axis is. The label is the exception,
//! drawn by the figure outside the clip and to the left of the widest axis
//! strip in the figure, so names line up whether or not a track drew an axis,
//! and cut down to the gutter when there is not room for the whole of it.

use std::fs;
use std::io;
use std::path::Path;

use crate::region::Region;
use crate::scale::Scale;
use crate::style::{Density, RenderProfile};
use crate::svg::{
    fit_text_by, mono_width, text_width, text_width_strong, Anchor, SvgWriter, TextStyle,
};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Rect, Track};

const DEFAULT_LABEL_WIDTH: f64 = 84.0;
const MIN_AUTO_LABEL_WIDTH: f64 = 48.0;
const MAX_AUTO_LABEL_WIDTH: f64 = 160.0;

/// Whitespace around the plotting area, in pixels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Margin {
    /// Space above the title.
    pub top: f64,
    /// Space to the right of the plotting area.
    pub right: f64,
    /// Space below the last track.
    pub bottom: f64,
    /// Space to the left of the track labels.
    pub left: f64,
}

impl Default for Margin {
    fn default() -> Self {
        Margin {
            top: 14.0,
            right: 18.0,
            bottom: 14.0,
            left: 16.0,
        }
    }
}

/// A stack of tracks sharing one horizontal coordinate system.
///
/// The figure owns the layout. Tracks say how tall they want to be, and the
/// figure decides where each band sits, reserves a gutter for labels when any
/// track asks for one, clips every track to its band, and works out the total
/// image height. Nothing about a track depends on its position in the stack, so
/// reordering them is a matter of reordering the [`Figure::push`] calls.
///
/// ```
/// use karyon::{AxisTrack, CoverageTrack, Figure, Region};
///
/// let region = Region::parse("NC_000962.3:761000-761500").unwrap();
/// let depth: Vec<f64> = (0..500).map(|i| 40.0 - (i as f64 / 25.0)).collect();
///
/// let svg = Figure::new(region)
///     .title("rpoB promoter")
///     .push(CoverageTrack::new(760_999, depth).label("depth"))
///     .push(AxisTrack::new())
///     .to_svg();
///
/// assert!(svg.starts_with("<svg"));
/// assert!(svg.ends_with("</svg>"));
/// ```
pub struct Figure {
    region: Region,
    width: f64,
    title: Option<String>,
    theme: Theme,
    tracks: Vec<Box<dyn Track>>,
    margin: Margin,
    label_width: Option<f64>,
    track_gap: f64,
    visual_scale: f64,
    density: Density,
    show_region_label: bool,
    description: Option<String>,
}

impl Figure {
    /// An empty figure over `region`, 900 pixels wide with the light theme.
    pub fn new(region: Region) -> Self {
        Figure {
            region,
            width: 900.0,
            title: None,
            theme: Theme::light(),
            tracks: Vec::new(),
            margin: Margin::default(),
            label_width: None,
            track_gap: 12.0,
            visual_scale: 1.0,
            density: Density::Balanced,
            show_region_label: true,
            description: None,
        }
    }

    /// Sets the image width in pixels.
    ///
    /// Widths that would leave no plotting area are raised to the smallest one
    /// that does, so a figure is always renderable. That smallest one depends
    /// on the margins and the label gutter, so it is worked out when the
    /// figure is laid out, and a [`Figure::margin`] or [`Figure::label_width`]
    /// counts the same written before this call as after it.
    pub fn width(mut self, width: f64) -> Self {
        self.width = width;
        self
    }

    /// Sets the title drawn above the tracks.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Sets the `<desc>` of the rendered document: what the figure shows.
    ///
    /// This is the alt text. It is read by a screen reader in place of the
    /// several thousand rectangles a figure is made of, and it is what a reader
    /// gets when the image does not load. Without it the document still carries
    /// a `<title>`, but a title says which locus this is and an alt text says
    /// what happens in it, and only the person drawing the figure knows that.
    ///
    /// ```
    /// use karyon::{AxisTrack, Figure, Region};
    ///
    /// let svg = Figure::new(Region::parse("chr7:1-1000").unwrap())
    ///     .description("Read depth falls to zero across the deleted exon.")
    ///     .push(AxisTrack::new())
    ///     .to_svg();
    ///
    /// assert!(svg.contains("<desc"));
    /// assert!(svg.contains("Read depth falls to zero"));
    /// ```
    pub fn description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }

    /// Replaces the theme.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme;
        self
    }

    /// Applies a named output profile as one coherent starting point.
    ///
    /// Profiles select palette, type scale and track density together. Any
    /// later [`Figure::theme`], [`Figure::visual_scale`] or [`Figure::density`]
    /// call may still override the relevant part.
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

    /// Scales the visual chrome without changing the canvas width or genomic
    /// coordinates.
    ///
    /// Fonts, margins, the label gutter, gaps and rounded corners move together,
    /// so a slide-sized figure can be made more assertive without retuning each
    /// value independently. Data marks retain their meaning and coordinates.
    /// The one width that moves is the floor: a width too narrow for the
    /// scaled margins and gutter is raised to what they need, as it is at the
    /// plain scale.
    pub fn visual_scale(mut self, factor: f64) -> Self {
        self.visual_scale = if factor.is_finite() {
            factor.max(0.25)
        } else {
            1.0
        };
        self
    }

    /// Sets the packing of repeated rows and track-internal marks.
    pub fn density(mut self, density: Density) -> Self {
        self.density = density;
        self
    }

    /// Replaces the margins.
    ///
    /// A side that is negative or not a number is taken as zero, since the
    /// margins are added into the total height and a total height that is not a
    /// number is written out as a document zero pixels tall.
    pub fn margin(mut self, margin: Margin) -> Self {
        let side = |value: f64| {
            if value.is_finite() {
                value.max(0.0)
            } else {
                0.0
            }
        };
        self.margin = Margin {
            top: side(margin.top),
            right: side(margin.right),
            bottom: side(margin.bottom),
            left: side(margin.left),
        };
        self
    }

    /// Sets the width of the left gutter holding track labels.
    ///
    /// The gutter is only reserved when at least one track returns a label, so
    /// a figure of unlabelled tracks uses the full width whatever this says.
    /// Without this override the gutter is measured from the widest label and
    /// capped so a single long sample name cannot consume the figure.
    ///
    /// A width that is negative or not a number is taken as zero, since the
    /// gutter is part of the floor the image width is held to, and an image
    /// width that is not a number is written out as a document zero pixels
    /// wide.
    pub fn label_width(mut self, width: f64) -> Self {
        self.label_width = Some(if width.is_finite() {
            width.max(0.0)
        } else {
            0.0
        });
        self
    }

    /// Sets the vertical gap between tracks.
    pub fn track_gap(mut self, gap: f64) -> Self {
        self.track_gap = gap.max(0.0);
        self
    }

    /// Shows or hides the locus string in the top right corner.
    ///
    /// A figure where no track [shows the window](Track::shows_region) never
    /// shows it, whatever this says: a locus above a phylogeny names a window
    /// the tree is not drawn in, which every figure is given, and reads as a
    /// claim about where the tree is.
    pub fn show_region_label(mut self, show: bool) -> Self {
        self.show_region_label = show;
        self
    }

    /// Whether the locus is drawn: asked for, and something in the figure
    /// shows the window it names.
    fn shows_region_label(&self) -> bool {
        self.show_region_label && self.names_region()
    }

    /// Whether the window means anything to a reader of this figure: empty,
    /// or holding a track that [shows where it is](Track::shows_region).
    fn names_region(&self) -> bool {
        self.tracks.is_empty() || self.tracks.iter().any(|track| track.shows_region())
    }

    /// Appends a track below the ones already added.
    pub fn push(mut self, track: impl Track + 'static) -> Self {
        self.tracks.push(Box::new(track));
        self
    }

    /// Appends a boxed track, for building a stack at runtime.
    pub fn push_boxed(mut self, track: Box<dyn Track>) -> Self {
        self.tracks.push(track);
        self
    }

    /// Puts a ruler under the last track measured against the coordinates,
    /// or at the bottom of a figure that has none.
    ///
    /// A ruler numbers what is above it. At the bottom of a stack that ends in
    /// a tree, or in a panel of sites laid out by their own index, it sat a
    /// track away from the coverage it numbered, under something it does not
    /// measure.
    pub(crate) fn push_ruler(mut self, ruler: impl Track + 'static) -> Self {
        let at = self
            .tracks
            .iter()
            .rposition(|track| track.on_coordinates())
            .map_or(self.tracks.len(), |last| last + 1);
        self.tracks.insert(at, Box::new(ruler));
        self
    }

    /// Whether a ruler would be measuring anything.
    ///
    /// True unless everything in the figure says otherwise. A stack of
    /// phylogenies says otherwise: a ruler under one measures a window that
    /// exists because every figure is given one, and not because the tree is
    /// anywhere in it.
    ///
    /// An empty figure is a window with nothing in it yet, and a window is
    /// worth showing, so it keeps its ruler.
    pub fn measures_coordinates(&self) -> bool {
        self.tracks.is_empty() || self.tracks.iter().any(|track| track.on_coordinates())
    }

    /// The region on display.
    pub fn region(&self) -> &Region {
        &self.region
    }

    /// How many tracks the figure holds.
    pub fn track_count(&self) -> usize {
        self.tracks.len()
    }

    /// Width and height of the rendered image in pixels.
    ///
    /// The height is computed, not configured: it follows from the tracks, so
    /// asking for it means laying the figure out.
    pub fn dimensions(&self) -> (f64, f64) {
        let layout = self.layout();
        (layout.width, layout.total_height)
    }

    /// What the document calls itself: the visible title, and the locus where
    /// anything is measured against it.
    ///
    /// Never empty, so a reader hovering the image is never told nothing. A
    /// figure with no title and nothing on coordinates, a phylogeny drawn on
    /// its own, is named by what its tracks are called; it was named by the
    /// window every figure is given, which for a tree was a placeholder the
    /// command line made up.
    fn document_name(&self) -> String {
        let region = self.names_region().then(|| self.region.to_string());
        match (&self.title, region) {
            (Some(title), Some(region)) => format!("{title}, {region}"),
            (Some(title), None) => title.clone(),
            (None, Some(region)) => region,
            (None, None) => {
                let labels: Vec<&str> = self.tracks.iter().filter_map(|t| t.label()).collect();
                if labels.is_empty() {
                    "A karyon figure".to_string()
                } else {
                    labels.join(", ")
                }
            }
        }
    }

    /// The alt text: whatever [`Figure::description`] was given, or a
    /// statement of what the figure is made of.
    ///
    /// The fallback is composed only from things the figure knows for certain,
    /// the region and the tracks in the order they are drawn, each by its
    /// label or, where it has none, by what it is, so it can say what is here
    /// without claiming anything about what it shows. That is the part
    /// [`Figure::description`] exists for.
    ///
    /// It counted every track and named only the labelled ones, so a figure
    /// of four tracks was three names long, the ruler said by nothing.
    fn document_description(&self) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        let named: Vec<&str> = self
            .tracks
            .iter()
            .map(|track| track.label().unwrap_or_else(|| track.noun()))
            .collect();
        // The window only where something shows it.
        let over = if self.names_region() {
            format!(" over {}", self.region)
        } else {
            String::new()
        };
        match named.as_slice() {
            [] => format!("A karyon figure{over}, with no tracks."),
            [one] => format!("A karyon figure{over}, with one track: {one}."),
            [first @ .., last] => format!(
                "A karyon figure{over}, with {} tracks, drawn top to bottom: {} and {last}.",
                named.len(),
                first.join(", ")
            ),
        }
    }

    /// Renders the figure to a standalone SVG document.
    pub fn to_svg(&self) -> String {
        self.to_svg_with_id_prefix("")
    }

    /// Renders the figure with every id it generates carrying `prefix`.
    ///
    /// Only needed when the result is going to be nested inside another SVG
    /// alongside a second figure. Ids are document-wide in SVG, so two figures
    /// in one document would otherwise both claim `karyon-clip-0` and the
    /// second one's clips would resolve to the first one's rectangles, cropping
    /// its tracks to the wrong bands. [`Panels`](crate::Panels) does this for
    /// you; do it yourself if you assemble a sheet by hand.
    pub fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        let theme = self.theme.clone().scaled(self.visual_scale);
        let layout = self.layout_with_theme(&theme);
        let mut svg = SvgWriter::with_id_prefix(prefix);
        // A prefix means this document is going inside another one, and a
        // nested document must not name itself. `<title>` resolves to the
        // innermost element under the pointer, so a title here would shadow
        // the one the sheet puts on the panel over the panel's whole area,
        // and `role="img"` inside another `role="img"` hides its contents
        // from a screen reader rather than describing them. The sheet is the
        // root and names itself; see [`Panels`](crate::Panels).
        if prefix.is_empty() {
            svg.describe(&self.document_name(), &self.document_description());
        }

        let locus = self.region.to_string();
        let (pill_width, pill_height) = locus_pill(&locus, &theme);
        if let Some(title) = &self.title {
            let room = if self.shows_region_label() {
                layout.width
                    - layout.margin_right
                    - pill_width
                    - theme.tokens.label_gap
                    - layout.margin_left
            } else {
                layout.width - layout.margin_right - layout.margin_left
            };
            // A title too long for its line is set a little smaller before it
            // is cut short, as a category too long for its strip is: a point
            // less type costs far less than the end of the sentence.
            let room = room.max(0.0);
            let largest = theme.title_font_size;
            let smallest = largest * 0.85;
            let mut size = largest;
            while size > smallest && text_width_strong(title, size) > room {
                size = (size - 0.5).max(smallest);
            }
            let visible_title = fit_text_by(title, room, |text| text_width_strong(text, size));
            svg.text_bold(
                layout.margin_left,
                layout.header_baseline,
                &visible_title,
                &theme.foreground,
                size,
                Anchor::Start,
            );
        }
        if self.shows_region_label() {
            // The locus in a pill, in the monospaced stack: it is a pair of
            // coordinates to be read digit by digit against the ruler, not a
            // phrase, and the pill keeps it from reading as the end of the
            // title beside it.
            let right = layout.width - layout.margin_right;
            let middle = layout.header_baseline - theme.font_size * 0.35;
            svg.rect_rounded(
                right - pill_width,
                middle - pill_height / 2.0,
                pill_width,
                pill_height,
                pill_height / 2.0,
                &mix(&theme.foreground, theme.surface(), 0.9),
            );
            svg.text_styled(
                right - pill_width / 2.0,
                layout.header_baseline,
                &locus,
                &theme.muted,
                theme.font_size,
                Anchor::Middle,
                TextStyle {
                    family: Some(&theme.mono_family),
                    ..TextStyle::default()
                },
            );
        }

        let mut y = layout.margin_top + layout.header_height;
        for (track, height) in self.tracks.iter().zip(&layout.track_heights) {
            let band = Rect {
                x: layout.plot_x,
                y,
                w: layout.plot_width,
                h: *height,
            };

            // The strip is what this track asked for, laid against the plot
            // area, and not the widest strip in the figure: a track that asked
            // for no axis gets none, and is clipped to its band alone rather
            // than to a neighbour's room.
            let axis_width = track.y_axis_width(&theme).max(0.0);
            let axis = Rect {
                x: band.x - axis_width,
                y,
                w: axis_width,
                h: *height,
            };

            // What the value axis measures goes under the name, a line of
            // its own in the ticks' ink, so the name and the title read as one
            // block centred on the band.
            let title = track.axis_title().filter(|title| !title.is_empty());
            let title_size = theme.font_size - 1.0;
            let lift = match (track.label(), title) {
                (Some(_), Some(_)) => (title_size + 3.0 * self.visual_scale) / 2.0,
                _ => 0.0,
            };
            if let Some(title) = title {
                let right = band.x - layout.axis_width - 10.0 * self.visual_scale;
                let visible = fit_text_by(title, right - layout.margin_left, |text| {
                    text_width(text, title_size)
                });
                let baseline = if track.label().is_some() {
                    band.mid_y() + theme.label_font_size * 0.35 + lift
                } else {
                    band.mid_y() + title_size * 0.35
                };
                svg.text(
                    right,
                    baseline,
                    &visible,
                    &theme.muted,
                    title_size,
                    Anchor::End,
                );
            }

            if let Some(label) = track.label() {
                // Labels sit to the left of the widest value axis, so a track
                // with an axis and one without still line their names up, and
                // they are cut down to the gutter, since a name wider than the
                // room reserved for it would start off the left edge of the
                // image and lose its first characters.
                //
                // They are set semibold and a little spaced, in the quieter
                // ink: the weight is what sets a track's name apart from the
                // tick labels beside it, so the ink can step back and leave
                // the data the loudest thing in the figure. The name keeps the
                // case it was given, since `katG` and `KATG` are not the same
                // gene.
                let right = band.x - layout.axis_width - 10.0 * self.visual_scale;
                let size = theme.label_font_size;
                let visible = fit_text_by(label, right - layout.margin_left, |text| {
                    label_width(text, size)
                });
                svg.text_styled(
                    right,
                    band.mid_y() + size * 0.35 - lift,
                    &visible,
                    &theme.muted,
                    size,
                    Anchor::End,
                    TextStyle {
                        weight: Some(600),
                        tracking: LABEL_TRACKING,
                        ..TextStyle::default()
                    },
                );
            }

            svg.begin_clip(axis.x, band.y, axis.w + band.w, band.h);
            let mut ctx = DrawContext {
                svg: &mut svg,
                scale: &layout.scale,
                theme: &theme,
                band,
                axis,
                region: &self.region,
                visual_scale: self.visual_scale * self.density.scale(),
            };
            track.draw(&mut ctx);
            svg.end_group();

            // No rule between one track and the next. The gap separates them
            // and each name sits level with the middle of its own band; a
            // hairline across the full width was one more line in a figure
            // made of lines, and the one that measured nothing.
            y += height + layout.track_gap;
        }

        svg.finish(
            layout.width,
            layout.total_height,
            &theme.background,
            &theme.font_family,
        )
    }

    /// Renders the figure and writes it to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save_svg(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, self.to_svg())
    }

    /// How many pixels one base is drawn across, once the gutter and the
    /// value axes are taken out of the width: what decides whether a base is
    /// a letter, a block of colour or nothing.
    pub fn px_per_bp(&self) -> f64 {
        self.layout().scale.px_per_bp()
    }

    /// The keys this figure's tracks need at the zoom it is drawn at, each
    /// once, for a [`LegendTrack`](crate::LegendTrack) to explain: the colours
    /// of bases drawn as blocks, for one. Empty where every mark names itself.
    pub fn key(&self) -> crate::track::legend::Legend {
        let px_per_bp = self.px_per_bp();
        self.tracks
            .iter()
            .filter_map(|track| track.key(&self.region, px_per_bp, &self.theme))
            .fold(crate::track::legend::Legend::new(), |key, more| {
                key.and(&more)
            })
    }

    fn layout(&self) -> Layout {
        let theme = self.theme.clone().scaled(self.visual_scale);
        self.layout_with_theme(&theme)
    }

    fn layout_with_theme(&self, theme: &Theme) -> Layout {
        // The floor is settled here and not in `width`, which only ever saw
        // the settings written before it: a gutter or a margin set after the
        // width never reached the floor, so the same settings written in two
        // orders drew two figures. The default is held to it as well, since
        // 900 is a width like any other and was set before everything else.
        // It is scaled as the layout below scales what it adds up: added up
        // unscaled, a figure at twice the visual scale was held to the width
        // of one at the plain scale, and its plotting area started past the
        // right edge of the image.
        let spacing = self.visual_scale;
        let floor = (self.margin.left
            + self.margin.right
            + self.label_width.unwrap_or(DEFAULT_LABEL_WIDTH)
            + 50.0)
            * spacing;
        let width = if self.width.is_finite() {
            self.width.max(floor)
        } else {
            floor
        };
        let margin_top = self.margin.top * spacing;
        let margin_right = self.margin.right * spacing;
        let margin_bottom = self.margin.bottom * spacing;
        let margin_left = self.margin.left * spacing;
        let track_gap = self.track_gap * spacing;
        let has_header = self.title.is_some() || self.shows_region_label();
        let header_height = if has_header {
            theme.title_font_size + 12.0 * spacing
        } else {
            0.0
        };
        let header_baseline = margin_top + theme.title_font_size;

        let named = self
            .tracks
            .iter()
            .any(|t| t.label().is_some() || t.axis_title().is_some());
        let gutter = if named {
            self.label_width.map_or_else(
                || self.automatic_label_width(theme),
                |width| width * spacing,
            )
        } else {
            0.0
        };
        // The widest axis any track asks for, reserved for all of them, so that
        // every plotting area still starts at the same x.
        let axis_width = self
            .tracks
            .iter()
            .map(|t| t.y_axis_width(theme).max(0.0))
            .fold(0.0f64, f64::max);
        let plot_x = margin_left + gutter + axis_width;
        let plot_width = (width - plot_x - margin_right).max(1.0);
        let scale = Scale::new(&self.region, plot_x, plot_width);

        // A height that is not a number is not a height. It would reach
        // `total_height`, which is written out as `height="0"`, and on a sheet
        // it would stack every later panel back at the top.
        let content_scale = self.visual_scale * self.density.scale();
        let track_heights: Vec<f64> = self
            .tracks
            .iter()
            .map(|t| {
                let height = t.height(&scale) * content_scale;
                if height.is_finite() {
                    height.max(1.0)
                } else {
                    1.0
                }
            })
            .collect();
        let content_height: f64 = track_heights.iter().sum::<f64>()
            + track_gap * (self.tracks.len().saturating_sub(1)) as f64;
        // Checking each height on its own is not enough: two of them can each
        // be a number and still add to one that is not, and then the total
        // goes out as `height="0"` while `dimensions` keeps saying infinity.
        // The ceiling is where an f64 stops holding consecutive integers, so a
        // figure taller than this could not state its own height exactly even
        // if something were willing to draw it.
        const TALLEST: f64 = (1u64 << 53) as f64;

        Layout {
            width,
            scale,
            plot_x,
            axis_width,
            plot_width,
            header_height,
            header_baseline,
            margin_top,
            margin_right,
            margin_left,
            track_gap,
            track_heights,
            total_height: (margin_top + header_height + content_height + margin_bottom)
                .min(TALLEST),
        }
    }

    /// Room for the widest label plus the quiet gap between labels and axes.
    fn automatic_label_width(&self, theme: &Theme) -> f64 {
        let names = self
            .tracks
            .iter()
            .filter_map(|track| track.label())
            .map(|label| label_width(label, theme.label_font_size));
        let titles = self
            .tracks
            .iter()
            .filter_map(|track| track.axis_title())
            .map(|title| text_width(title, theme.font_size - 1.0));
        let widest = names.chain(titles).fold(0.0f64, f64::max);
        (widest + 14.0 * self.visual_scale).clamp(
            MIN_AUTO_LABEL_WIDTH * self.visual_scale,
            MAX_AUTO_LABEL_WIDTH * self.visual_scale,
        )
    }
}

/// Extra space after every letter of a track's name, in ems.
const LABEL_TRACKING: f64 = 0.03;

/// How wide a track's name is as it is set: semibold, and spaced out.
fn label_width(label: &str, size: f64) -> f64 {
    text_width_strong(label, size) + LABEL_TRACKING * size * label.chars().count() as f64
}

/// The width and height of the pill the locus sits in beside the title.
fn locus_pill(locus: &str, theme: &Theme) -> (f64, f64) {
    let height = theme.font_size + 9.0 * theme.font_size / 11.5;
    (mono_width(locus, theme.font_size) + height, height)
}

impl crate::rings::Drawing for Figure {
    fn dimensions(&self) -> (f64, f64) {
        Figure::dimensions(self)
    }

    fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        Figure::to_svg_with_id_prefix(self, prefix)
    }

    fn content_anchor(&self) -> Option<f64> {
        Some(self.layout().plot_x)
    }

    fn region(&self) -> Option<&Region> {
        // The test that decides whether a plot gets a ruler.
        // A window nothing is measured against is not one worth moving.
        self.measures_coordinates().then_some(&self.region)
    }
}

struct Layout {
    /// The width the figure is drawn at: the one it was given, held to the
    /// floor.
    width: f64,
    scale: Scale,
    plot_x: f64,
    axis_width: f64,
    plot_width: f64,
    header_height: f64,
    header_baseline: f64,
    margin_top: f64,
    margin_right: f64,
    margin_left: f64,
    track_gap: f64,
    track_heights: Vec<f64>,
    total_height: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::track::{AxisTrack, CoverageTrack, Feature, FeatureTrack};

    /// A track that answers with whatever height it was built with.
    struct Tall(f64);

    impl Track for Tall {
        fn height(&self, _scale: &Scale) -> f64 {
            self.0
        }
        fn draw(&self, ctx: &mut DrawContext<'_>) {
            let _ = ctx.scale;
        }
    }

    #[test]
    fn two_heights_that_are_numbers_cannot_add_to_one_that_is_not() {
        // Every height is checked for being a number on its own, and that is
        // one track at a time. Two of them can each pass and still overflow
        // the sum, and then the total went out as `height="0"` while
        // `dimensions` reported infinity: the file and the API disagreed about
        // the same figure.
        let figure = Figure::new(region())
            .show_region_label(false)
            .push(Tall(1e300))
            .push(Tall(f64::MAX));
        let (width, height) = figure.dimensions();
        assert!(
            width.is_finite() && height.is_finite() && height > 0.0,
            "dimensions {width} x {height}"
        );
        let svg = figure.to_svg();
        let stated: f64 = svg
            .split("height=\"")
            .nth(1)
            .and_then(|rest| rest.split('"').next())
            .and_then(|value| value.parse().ok())
            .expect("the document states a height");
        assert!(
            (stated - height).abs() < 0.01,
            "the document says {stated} and dimensions says {height}"
        );
    }

    fn region() -> Region {
        Region::parse("chr1:1-1000").unwrap()
    }

    #[test]
    fn an_empty_figure_still_renders_a_valid_document() {
        let svg = Figure::new(region()).to_svg();
        assert!(svg.starts_with("<svg"));
        assert!(svg.ends_with("</svg>"));
        assert!(svg.contains("chr1:1-1000"));
    }

    #[test]
    fn height_grows_with_each_track() {
        let one = Figure::new(region()).push(AxisTrack::new().height(20.0));
        let two = Figure::new(region())
            .push(AxisTrack::new().height(20.0))
            .push(AxisTrack::new().height(20.0));
        let (_, h1) = one.dimensions();
        let (_, h2) = two.dimensions();
        assert_eq!(h2 - h1, 20.0 + 12.0);
    }

    #[test]
    fn named_profiles_move_type_marks_and_density_as_one_system() {
        let manuscript = Figure::new(region())
            .title("profile")
            .push(FeatureTrack::new(vec![Feature::new(0, 10)]));
        let presentation = Figure::new(region())
            .title("profile")
            .profile(RenderProfile::Presentation)
            .push(FeatureTrack::new(vec![Feature::new(0, 10)]));
        assert!(presentation.dimensions().1 > manuscript.dimensions().1);
        let svg = presentation.to_svg();
        assert!(
            svg.contains(r#"font-size="22.95""#),
            "a 17 px title at 1.35"
        );
    }

    #[test]
    fn a_long_header_keeps_its_exact_accessible_name_when_the_visible_line_is_fitted() {
        let title = "A deliberately long title that would otherwise run through the locus label";
        let svg = Figure::new(region()).width(260.0).title(title).to_svg();
        assert!(svg.contains("\u{2026}</text>"), "{svg}");
        assert!(svg.contains(&format!(
            "<title id=\"karyon-title\">{title}, chr1:1-1000</title>"
        )));
    }

    #[test]
    fn the_label_gutter_is_only_reserved_when_a_track_wants_it() {
        let bare = Figure::new(region()).push(AxisTrack::new());
        let labelled = Figure::new(region()).push(AxisTrack::new().label("pos"));
        assert_eq!(bare.layout().plot_x, 16.0);
        assert_eq!(labelled.layout().plot_x, 16.0 + MIN_AUTO_LABEL_WIDTH);
        assert!(labelled.layout().plot_width < bare.layout().plot_width);
    }

    #[test]
    fn a_value_axis_is_only_reserved_when_a_track_asks_for_one() {
        let bare = Figure::new(region()).push(AxisTrack::new());
        let quantitative = Figure::new(region()).push(CoverageTrack::new(0, vec![1.0; 1000]));
        assert_eq!(bare.layout().axis_width, 0.0);
        assert!(quantitative.layout().axis_width > 0.0);
        assert!(quantitative.layout().plot_x > bare.layout().plot_x);
    }

    #[test]
    fn one_track_asking_for_an_axis_moves_every_plotting_area() {
        // The whole point of reserving the widest request for all of them: two
        // tracks in one figure must still share an x axis.
        let figure = Figure::new(region())
            .push(CoverageTrack::new(0, vec![1.0; 1000]).label("depth"))
            .push(AxisTrack::new());
        let layout = figure.layout();
        assert!(layout.axis_width > 0.0);
        // Both bands start at plot_x, because there is only one plot_x.
        let gutter = figure.automatic_label_width(&Theme::light());
        assert_eq!(layout.plot_x, 16.0 + gutter + layout.axis_width);
    }

    #[test]
    fn turning_the_axis_off_gives_the_room_back() {
        let with = Figure::new(region()).push(CoverageTrack::new(0, vec![1.0; 1000]));
        let without =
            Figure::new(region()).push(CoverageTrack::new(0, vec![1.0; 1000]).show_max(false));
        assert!(without.layout().plot_width > with.layout().plot_width);
        assert_eq!(without.layout().axis_width, 0.0);
    }

    #[test]
    fn every_track_is_clipped_to_its_own_band() {
        let svg = Figure::new(region())
            .push(AxisTrack::new())
            .push(CoverageTrack::new(0, vec![1.0; 1000]))
            .to_svg();
        assert_eq!(svg.matches("<clipPath").count(), 2);
        assert_eq!(svg.matches("clip-path=").count(), 2);
        assert_eq!(svg.matches("</g>").count(), 2);
    }

    #[test]
    fn tracks_are_stacked_top_to_bottom_in_push_order() {
        let figure = Figure::new(region())
            .push(AxisTrack::new().height(20.0))
            .push(AxisTrack::new().height(40.0));
        let layout = figure.layout();
        assert_eq!(layout.track_heights, vec![20.0, 40.0]);
    }

    #[test]
    fn an_absurd_width_is_raised_to_something_renderable() {
        let figure = Figure::new(region()).width(1.0).push(AxisTrack::new());
        let (width, _) = figure.dimensions();
        assert!(width >= 50.0);
        assert!(figure.layout().plot_width >= 1.0);
        assert!(figure.to_svg().starts_with("<svg"));
    }

    #[test]
    fn the_settings_a_width_is_floored_against_count_before_it_or_after_it() {
        // The smallest width that leaves a plotting area depends on the
        // margins and the label gutter, and `width` used to work it out on the
        // spot, from whatever those two were when it was called. A gutter or a
        // margin written after it never reached the floor, so the same
        // settings written in two orders drew two different figures.
        let draw = |figure: Figure| {
            let figure = figure
                .push(CoverageTrack::new(0, vec![10.0; 1000]).label("depth"))
                .push(AxisTrack::new());
            (figure.dimensions(), figure.to_svg())
        };
        let wide = Margin {
            left: 120.0,
            right: 120.0,
            ..Margin::default()
        };
        let orders = [
            (
                Figure::new(region()).width(200.0).label_width(150.0),
                Figure::new(region()).label_width(150.0).width(200.0),
            ),
            (
                Figure::new(region()).width(200.0).margin(wide),
                Figure::new(region()).margin(wide).width(200.0),
            ),
        ];
        for (width_first, width_last) in orders {
            let (first, last) = (draw(width_first), draw(width_last));
            assert_eq!(first.0, last.0, "the width written first, then last");
            assert_eq!(first.1, last.1);
        }
    }

    #[test]
    fn the_default_width_is_held_to_the_same_floor_as_one_that_was_asked_for() {
        // 900 is a width like any other, and it is set when the figure is
        // made, so before every other setting. A gutter wider than it used to
        // leave the plotting area off the right of the image, unless `width`
        // was called afterwards, even with 900, to have it raised.
        let gutter = |figure: Figure| {
            figure
                .label_width(1_000.0)
                .push(AxisTrack::new().label("position"))
        };
        let unset = gutter(Figure::new(region()));
        let set = gutter(Figure::new(region())).width(900.0);
        assert_eq!(unset.dimensions(), set.dimensions());
        assert_eq!(unset.to_svg(), set.to_svg());
        let (width, _) = unset.dimensions();
        let layout = unset.layout();
        assert!(
            layout.plot_x + layout.plot_width <= width,
            "the plotting area runs from {} to {} on an image {width} wide",
            layout.plot_x,
            layout.plot_x + layout.plot_width
        );
    }

    #[test]
    fn the_floor_grows_with_the_visual_scale_the_margins_and_gutter_grow_with() {
        // The layout multiplies the margins and the gutter by the visual
        // scale, and the floor added them up unmultiplied. At twice the scale
        // a figure held to its floor was 234 pixels wide and started its
        // plotting area at 332, so the track and its name were both drawn
        // off the right of the image.
        let floored = |scale: f64| {
            Figure::new(region())
                .visual_scale(scale)
                .label_width(150.0)
                .width(1.0)
                .push(CoverageTrack::new(0, vec![10.0; 1000]).label("depth"))
        };
        let unscaled = floored(1.0).layout().plot_width;
        for scale in [1.0, 2.0, 3.5] {
            let figure = floored(scale);
            let (width, _) = figure.dimensions();
            let layout = figure.layout();
            assert!(
                layout.plot_x + layout.plot_width + layout.margin_right <= width,
                "at {scale} the plotting area runs from {} to {} on an image {width} wide",
                layout.plot_x,
                layout.plot_x + layout.plot_width
            );
            // The area grows with everything else, rather than being what is
            // left of an image that did not.
            assert!(
                layout.plot_width >= unscaled * scale * 0.99,
                "at {scale} the plotting area is {} wide, against {unscaled} unscaled",
                layout.plot_width
            );
        }
        // Nothing moves at the scale every committed figure is drawn at.
        assert_eq!(floored(1.0).dimensions().0, 234.0);
    }

    #[test]
    fn a_gutter_that_is_not_finite_does_not_make_a_document_zero_pixels_wide() {
        // The gutter is part of the floor a width is held to, so an infinite
        // one made the floor infinite, and that went out as `width="0"`. It
        // did so already when `width` came after it, and with the floor
        // settled at layout it would have in every order.
        for figure in [
            Figure::new(region()).label_width(f64::INFINITY),
            Figure::new(region())
                .label_width(f64::INFINITY)
                .width(900.0),
            Figure::new(region())
                .width(900.0)
                .label_width(f64::INFINITY),
        ] {
            let figure = figure.push(AxisTrack::new().label("position"));
            let (width, _) = figure.dimensions();
            assert_eq!(width, 900.0);
            let layout = figure.layout();
            assert!(
                layout.plot_x + layout.plot_width <= width,
                "the plotting area runs from {} to {} on an image {width} wide",
                layout.plot_x,
                layout.plot_x + layout.plot_width
            );
        }
    }

    #[test]
    fn a_non_finite_width_does_not_leak_into_the_output() {
        let svg = Figure::new(region())
            .width(f64::NAN)
            .push(AxisTrack::new())
            .to_svg();
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_region_label_can_be_turned_off() {
        // What is turned off is the drawn label. The document still calls
        // itself by its locus, since a figure with no accessible name is
        // worse than one with a name nobody can see.
        let svg = Figure::new(region()).show_region_label(false).to_svg();
        assert!(!drawn_text(&svg).contains("chr1:1-1000"));
        assert!(svg.contains("<title id=\"karyon-title\">chr1:1-1000</title>"));
    }

    #[test]
    fn a_figure_nothing_is_measured_in_does_not_print_or_speak_its_window() {
        // A phylogeny is not drawn in the window every figure is given, so a
        // locus above it, or in its accessible name, is a claim about where
        // the tree is. The command line had to make one up for a tree, and
        // printed it: "x:1-1".
        let tree = crate::Tree::parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let bare = Figure::new(region())
            .push(crate::TreeTrack::new(tree.clone()))
            .to_svg();
        assert!(!bare.contains("chr1:1-1000"), "{bare}");
        assert!(bare.contains("<title id=\"karyon-title\">A karyon figure</title>"));
        assert!(bare.contains("A karyon figure, with one track: a phylogeny."));

        let titled = Figure::new(region())
            .title("Outbreak")
            .push(crate::TreeTrack::new(tree.clone()).label("phylogeny"))
            .to_svg();
        assert!(!titled.contains("chr1:1-1000"), "{titled}");
        assert!(titled.contains("<title id=\"karyon-title\">Outbreak</title>"));

        let labelled = Figure::new(region())
            .push(crate::TreeTrack::new(tree).label("phylogeny"))
            .to_svg();
        assert!(labelled.contains("<title id=\"karyon-title\">phylogeny</title>"));

        // An ideogram is not on the coordinates either, and it marks the
        // window on its chromosome, so the locus names something on the page.
        let banded = Figure::new(region())
            .push(crate::IdeogramTrack::new(5_000, Vec::new()))
            .to_svg();
        assert!(drawn_text(&banded).contains("chr1:1-1000"), "{banded}");

        // One track on coordinates and the window is back, pill and all.
        let mixed = Figure::new(region())
            .push(crate::TreeTrack::new(
                crate::Tree::parse_newick("(A:1,B:1);").unwrap(),
            ))
            .push(AxisTrack::new())
            .to_svg();
        assert!(drawn_text(&mixed).contains("chr1:1-1000"));
    }

    /// Everything the figure actually draws as text, with the title and the
    /// description left out.
    fn drawn_text(svg: &str) -> String {
        svg.split("<text")
            .skip(1)
            .filter_map(|piece| piece.split_once('>'))
            .filter_map(|(_, rest)| rest.split_once("</text>"))
            .map(|(content, _)| content)
            .collect::<Vec<_>>()
            .join(" ")
    }

    #[test]
    fn a_label_wider_than_the_gutter_is_shortened_rather_than_run_off_the_canvas() {
        // An explicitly narrow gutter keeps long names from silently running
        // off the canvas: they retain their recognisable start and an ellipsis.
        let name = "NC_000962.3 read depth";
        let svg = Figure::new(region())
            .label_width(DEFAULT_LABEL_WIDTH)
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 1000]).label(name))
            .to_svg();

        let drawn = drawn_text(&svg);
        assert!(!drawn.contains(name), "the whole name was drawn: {drawn}");
        let shown = drawn
            .split_whitespace()
            .find(|word| word.ends_with('\u{2026}'))
            .unwrap_or_else(|| panic!("no shortened name in {drawn}"));
        assert!(
            name.starts_with(shown.trim_end_matches('\u{2026}')),
            "{shown}"
        );
        assert!(shown.len() > "NC_0".len(), "kept too little of it: {shown}");
        // Right aligned at x = 90, so the ink starts inside the left margin,
        // measured as it is set: semibold, and spaced.
        assert!(90.0 - label_width(shown, Theme::light().label_font_size) >= 16.0);
    }

    #[test]
    fn a_label_that_fits_the_gutter_is_left_exactly_as_it_was() {
        let name = "enrich / deplete";
        let svg = Figure::new(region())
            .label_width(130.0)
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 1000]).label(name))
            .to_svg();
        assert!(drawn_text(&svg).contains(name), "{svg}");
    }

    #[test]
    fn a_gutter_with_no_room_draws_no_label_rather_than_one_off_the_canvas() {
        // With the gutter turned off the label has no drawable room.
        let svg = Figure::new(region())
            .label_width(0.0)
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 1000]).label("dp"))
            .to_svg();
        assert!(!drawn_text(&svg).contains("dp"), "{svg}");
    }

    #[test]
    fn a_track_that_asked_for_no_axis_is_clipped_to_its_band_alone() {
        // The feature track asks for no axis, and a gene overhanging the left
        // of the window has to stop at the plot origin rather than inside the
        // quantitative track's y-axis strip.
        let figure = Figure::new(Region::new("chr1", 1000, 2000).unwrap())
            .show_region_label(false)
            .push(CoverageTrack::new(1000, vec![40.0; 1000]))
            .push(FeatureTrack::new(vec![Feature::new(0, 1500)]));
        let layout = figure.layout();
        assert!(layout.axis_width > 0.0, "{}", layout.axis_width);
        assert_eq!(layout.plot_x, 16.0 + layout.axis_width);

        let svg = figure.to_svg();
        let clips: Vec<f64> = svg
            .match_indices("<clipPath")
            .map(|(at, _)| {
                let rest = &svg[at..];
                let x = rest.find(r#"x=""#).unwrap() + 3;
                rest[x..].split('"').next().unwrap().parse().unwrap()
            })
            .collect();
        // The document writes three decimals, so that is as close as the
        // comparison can be.
        assert_eq!(clips.len(), 2);
        assert!(
            (clips[0] - 16.0).abs() < 1e-3,
            "the strip this track asked for: {clips:?}"
        );
        assert!(
            (clips[1] - layout.plot_x).abs() < 1e-3,
            "and no strip at all for the track that asked for none: {clips:?}"
        );
    }

    #[test]
    fn a_label_still_lines_up_with_its_neighbours_when_the_strips_differ() {
        // The strips are now the track's own, but the names are not: they hang
        // off the widest strip in the figure so that a track with an axis and
        // one without still read down a single edge.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![10.0; 1000]).label("depth"))
            .push(AxisTrack::new().label("position"))
            .to_svg();
        let anchors: Vec<&str> = svg
            .match_indices(r#"<text x=""#)
            .map(|(at, prefix)| {
                let rest = &svg[at + prefix.len()..];
                &rest[..rest.find('"').unwrap()]
            })
            .collect();
        assert!(anchors.windows(2).any(|pair| pair[0] == pair[1]), "{svg}");
    }

    #[test]
    fn a_non_finite_margin_is_taken_as_no_margin_rather_than_no_document() {
        // Margin is added into the total height, and a total height that is
        // not a number is written out as a document zero pixels tall.
        let figure = Figure::new(region())
            .margin(Margin {
                top: f64::NAN,
                right: 16.0,
                bottom: -4.0,
                left: 12.0,
            })
            .push(AxisTrack::new());
        assert_eq!(figure.margin.top, 0.0);
        assert_eq!(figure.margin.bottom, 0.0);
        let (_, height) = figure.dimensions();
        assert!(height.is_finite() && height > 0.0, "{height}");
        assert!(!figure.to_svg().contains(r#"height="0""#));
    }

    #[test]
    fn a_non_finite_track_height_is_floored_rather_than_flattening_the_document() {
        let broken = Figure::new(region())
            .show_region_label(false)
            .push(AxisTrack::new().height(f64::INFINITY));
        let (_, height) = broken.dimensions();
        assert_eq!(height, 29.0, "14 of margin, one pixel of band, 14 more");
        let svg = broken.to_svg();
        assert!(!svg.contains(r#"height="0""#), "{svg}");
        assert!(!svg.contains("inf"), "{svg}");
    }

    #[test]
    fn no_rule_runs_between_one_track_and_the_next() {
        // The gap separates the bands. A line across the whole width between
        // them was the one line in the figure that measured nothing.
        let theme = Theme::light();
        let two = Figure::new(region())
            .show_region_label(false)
            .push(FeatureTrack::new(vec![Feature::new(0, 10)]).label("genes"))
            .push(FeatureTrack::new(vec![Feature::new(0, 10)]).label("more"))
            .to_svg();
        let full_width = format!(r#"<line x1="{}" "#, crate::svg::num(Margin::default().left));
        assert!(!two.contains(&full_width), "{two}");
        assert!(!two.contains(&theme.rule), "{two}");
    }

    #[test]
    fn a_track_is_named_semibold_in_the_quiet_ink_and_keeps_its_case() {
        // The weight sets a name apart from the tick labels beside it, so the
        // ink can step back; and a name is never re-cased, since `katG` and
        // `KATG` are different genes.
        let theme = Theme::light();
        let svg = Figure::new(region())
            .push(FeatureTrack::new(vec![Feature::new(0, 10)]).label("katG"))
            .to_svg();
        let named = format!(
            r#"fill="{}" font-size="{}" text-anchor="end" font-weight="600" letter-spacing="0.03em">katG<"#,
            theme.muted,
            crate::svg::num(theme.label_font_size)
        );
        assert!(svg.contains(&named), "{svg}");
    }

    /// A track that says what its axis measures has it written under its
    /// name, the two lines centred on the band together, and a title alone
    /// still has a gutter to go in.
    #[test]
    fn an_axis_title_goes_under_the_track_s_name() {
        struct Titled(Option<&'static str>);
        impl Track for Titled {
            fn height(&self, _scale: &Scale) -> f64 {
                60.0
            }
            fn label(&self) -> Option<&str> {
                self.0
            }
            fn axis_title(&self) -> Option<&str> {
                Some("-log10 p")
            }
            fn draw(&self, _ctx: &mut DrawContext<'_>) {}
        }
        let y_of = |svg: &str, text: &str| -> f64 {
            let at = svg.find(&format!(">{text}<")).expect("the text");
            let tag = &svg[svg[..at].rfind("<text").unwrap()..at];
            let y = &tag[tag.find(" y=\"").unwrap() + 4..];
            y[..y.find('"').unwrap()].parse().unwrap()
        };
        let both = Figure::new(region()).push(Titled(Some("scan"))).to_svg();
        let (name, title) = (y_of(&both, "scan"), y_of(&both, "-log10 p"));
        assert!(
            title > name,
            "the title at {title} is not under the name at {name}"
        );
        // Alone, the title takes the name's place, and the gutter is kept for it.
        let alone = Figure::new(region()).push(Titled(None));
        assert!(alone.layout().plot_x > Margin::default().left);
        assert!(alone.to_svg().contains(">-log10 p<"));
        // A title wider than any name widens the gutter to fit it.
        struct Wide;
        impl Track for Wide {
            fn height(&self, _scale: &Scale) -> f64 {
                40.0
            }
            fn label(&self) -> Option<&str> {
                Some("a")
            }
            fn axis_title(&self) -> Option<&str> {
                Some("a title a good deal wider than the name")
            }
            fn draw(&self, _ctx: &mut DrawContext<'_>) {}
        }
        let wide = Figure::new(region()).push(Wide);
        let named = Figure::new(region()).push(Titled(Some("a")));
        assert!(wide.layout().plot_x > named.layout().plot_x);
    }

    /// The keys are the ones the tracks need at the zoom the figure draws at,
    /// which only the figure knows once it has laid out its gutter.
    #[test]
    fn a_figure_asks_its_tracks_for_the_keys_its_zoom_needs() {
        use crate::track::SequenceTrack;
        let bases = b"ACGT".repeat(1_000);
        let blocks = Figure::new(Region::new("chr1", 0, 400).unwrap())
            .push(SequenceTrack::new(0, bases.clone()).label("reference"));
        assert_eq!(blocks.key().len(), 4);
        let letters = Figure::new(Region::new("chr1", 0, 40).unwrap())
            .push(SequenceTrack::new(0, bases).label("reference"));
        assert!(letters.key().is_empty());
    }

    /// The description counted every track and named only the labelled
    /// ones, so four tracks were three names long; each is named now, by its
    /// label or by what it is.
    #[test]
    fn the_description_names_every_track_it_counts() {
        use crate::track::legend::{Legend, LegendTrack};
        use crate::track::CoverageTrack;
        let svg = Figure::new(region())
            .push(CoverageTrack::new(0, vec![1.0; 1000]).label("depth"))
            .push(AxisTrack::new())
            .push(LegendTrack::new(Legend::new().key("A", "#111111")))
            .to_svg();
        assert!(
            svg.contains(
                "with 3 tracks, drawn top to bottom: depth, a ruler and a key to the colours."
            ),
            "{svg}"
        );
        // A track that says nothing of itself is still counted and named.
        struct Quiet;
        impl Track for Quiet {
            fn height(&self, _scale: &Scale) -> f64 {
                10.0
            }
            fn draw(&self, _ctx: &mut DrawContext<'_>) {}
        }
        let svg = Figure::new(region()).push(Quiet).push(Quiet).to_svg();
        assert!(
            svg.contains("with 2 tracks, drawn top to bottom: a track and a track."),
            "{svg}"
        );
    }

    #[test]
    fn the_locus_sits_in_a_pill_in_the_monospaced_stack() {
        let theme = Theme::light();
        let svg = Figure::new(region()).title("a locus").to_svg();
        let locus = format!(r#"font-family="{}">chr1:1-1000</text>"#, theme.mono_family);
        assert!(svg.contains(&locus), "{svg}");
        let pill = mix(&theme.foreground, theme.surface(), 0.9);
        assert!(svg.contains(&format!(r#"fill="{pill}""#)), "{svg}");
    }

    #[test]
    fn the_title_is_drawn_when_given() {
        let svg = Figure::new(region()).title("H37Rv rpoB").to_svg();
        assert!(svg.contains("H37Rv rpoB"));
        assert!(svg.contains("font-weight=\"bold\""));
    }

    #[test]
    fn the_dark_theme_paints_a_dark_page() {
        let svg = Figure::new(region()).theme(Theme::dark()).to_svg();
        assert!(svg.contains(&Theme::dark().background));
    }

    #[test]
    fn automatic_label_width_tracks_content_and_has_a_cap() {
        let short = Figure::new(region()).push(AxisTrack::new().label("pos"));
        let medium = Figure::new(region()).push(AxisTrack::new().label("chromosome position"));
        let huge = Figure::new(region()).push(AxisTrack::new().label("x".repeat(200)));

        assert!(medium.layout().plot_x > short.layout().plot_x);
        assert_eq!(
            huge.layout().plot_x,
            Margin::default().left + MAX_AUTO_LABEL_WIDTH
        );
    }

    #[test]
    fn visual_scale_moves_type_and_spacing_together() {
        let normal = Figure::new(region())
            .title("Scaled")
            .push(AxisTrack::new().label("position"));
        let scaled = Figure::new(region())
            .title("Scaled")
            .visual_scale(1.5)
            .push(AxisTrack::new().label("position"));

        assert!(scaled.layout().plot_x > normal.layout().plot_x);
        assert!(scaled.dimensions().1 > normal.dimensions().1);
        assert!(
            scaled.to_svg().contains(r#"font-size="25.5""#),
            "a 17 px title at 1.5"
        );
    }

    #[test]
    fn a_figure_offers_its_region_only_where_something_is_measured_against_it() {
        use crate::rings::Drawing;
        use crate::track::TreeTrack;
        use crate::tree::Tree;

        let depth = Figure::new(region()).push(CoverageTrack::new(0, vec![5.0; 1_000]));
        assert_eq!(Drawing::region(&depth), Some(&region()));

        // A window a phylogeny is handed because every figure is handed one is
        // not a window a viewer can move along, whatever it says.
        let tree = Tree::parse_newick("((a:1,b:1):1,c:2);").unwrap();
        let alone = Figure::new(Region::parse("tree:1-1").unwrap()).push(TreeTrack::new(tree));
        assert_eq!(Drawing::region(&alone), None);

        // And an empty figure is a window with nothing in it yet, which is
        // still a window.
        assert_eq!(Drawing::region(&Figure::new(region())), Some(&region()));
    }
}
