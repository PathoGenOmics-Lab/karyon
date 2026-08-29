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
use crate::svg::{fit_text, text_width, Anchor, SvgWriter};
use crate::theme::Theme;
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
    /// that does, so a figure is always renderable.
    pub fn width(mut self, width: f64) -> Self {
        let floor = self.margin.left
            + self.margin.right
            + self.label_width.unwrap_or(DEFAULT_LABEL_WIDTH)
            + 50.0;
        self.width = if width.is_finite() {
            width.max(floor)
        } else {
            floor
        };
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
    pub fn label_width(mut self, width: f64) -> Self {
        self.label_width = Some(width.max(0.0));
        self
    }

    /// Sets the vertical gap between tracks.
    pub fn track_gap(mut self, gap: f64) -> Self {
        self.track_gap = gap.max(0.0);
        self
    }

    /// Shows or hides the locus string in the top right corner.
    pub fn show_region_label(mut self, show: bool) -> Self {
        self.show_region_label = show;
        self
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

    /// Whether a ruler along the bottom would be measuring anything.
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
        (self.width, layout.total_height)
    }

    /// What the document calls itself: the visible title, or the locus.
    ///
    /// A figure always has one of the two, so the `<title>` is never empty and
    /// a reader hovering the image is never told nothing.
    fn document_name(&self) -> String {
        match &self.title {
            Some(title) => format!("{}, {}", title, self.region),
            None => self.region.to_string(),
        }
    }

    /// The alt text: whatever [`Figure::description`] was given, or a
    /// statement of what the figure is made of.
    ///
    /// The fallback is composed only from things the figure knows for certain,
    /// the region and the labels of the tracks in the order they are drawn, so
    /// it can say what is here without claiming anything about what it shows.
    /// That is the part [`Figure::description`] exists for.
    fn document_description(&self) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        let labels: Vec<&str> = self.tracks.iter().filter_map(|t| t.label()).collect();
        let count = self.tracks.len();
        let stack = match count {
            0 => "no tracks".to_string(),
            1 => "one track".to_string(),
            n => format!("{n} tracks"),
        };
        if labels.is_empty() {
            format!("A karyon figure over {}, with {}.", self.region, stack)
        } else {
            format!(
                "A karyon figure over {}, with {}, drawn top to bottom: {}.",
                self.region,
                stack,
                labels.join(", ")
            )
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

        if let Some(title) = &self.title {
            let room = if self.show_region_label {
                let locus_width = text_width(&self.region.to_string(), theme.font_size);
                self.width
                    - layout.margin_right
                    - locus_width
                    - theme.tokens.label_gap
                    - layout.margin_left
            } else {
                self.width - layout.margin_right - layout.margin_left
            };
            let visible_title = fit_text(title, room.max(0.0), theme.title_font_size);
            svg.text_bold(
                layout.margin_left,
                layout.header_baseline,
                &visible_title,
                &theme.foreground,
                theme.title_font_size,
                Anchor::Start,
            );
        }
        if self.show_region_label {
            svg.text(
                self.width - layout.margin_right,
                layout.header_baseline,
                &self.region.to_string(),
                &theme.muted,
                theme.font_size,
                Anchor::End,
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

            if let Some(label) = track.label() {
                // Labels sit to the left of the widest value axis, so a track
                // with an axis and one without still line their names up, and
                // they are cut down to the gutter, since a name wider than the
                // room reserved for it would start off the left edge of the
                // image and lose its first characters.
                let right = band.x - layout.axis_width - 10.0 * self.visual_scale;
                let visible = fit_text(label, right - layout.margin_left, theme.label_font_size);
                svg.text(
                    right,
                    band.mid_y() + theme.label_font_size * 0.35,
                    &visible,
                    &theme.muted,
                    theme.label_font_size,
                    Anchor::End,
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

            y += height + layout.track_gap;
        }

        svg.finish(
            self.width,
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

    fn layout(&self) -> Layout {
        let theme = self.theme.clone().scaled(self.visual_scale);
        self.layout_with_theme(&theme)
    }

    fn layout_with_theme(&self, theme: &Theme) -> Layout {
        let spacing = self.visual_scale;
        let margin_top = self.margin.top * spacing;
        let margin_right = self.margin.right * spacing;
        let margin_bottom = self.margin.bottom * spacing;
        let margin_left = self.margin.left * spacing;
        let track_gap = self.track_gap * spacing;
        let has_header = self.title.is_some() || self.show_region_label;
        let header_height = if has_header {
            theme.title_font_size + 12.0 * spacing
        } else {
            0.0
        };
        let header_baseline = margin_top + theme.title_font_size;

        let gutter = if self.tracks.iter().any(|t| t.label().is_some()) {
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
        let plot_width = (self.width - plot_x - margin_right).max(1.0);
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
        let widest = self
            .tracks
            .iter()
            .filter_map(|track| track.label())
            .map(|label| text_width(label, theme.label_font_size))
            .fold(0.0f64, f64::max);
        (widest + 14.0 * self.visual_scale).clamp(
            MIN_AUTO_LABEL_WIDTH * self.visual_scale,
            MAX_AUTO_LABEL_WIDTH * self.visual_scale,
        )
    }
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
}

struct Layout {
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
        assert!(svg.contains(r#"font-size="24.3""#), "{svg}");
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
        assert!(drawn.contains("NC_000962\u{2026}"), "{drawn}");
        // Right aligned at x = 90, so the ink starts inside the left margin.
        assert!(90.0 - text_width("NC_000962\u{2026}", 12.0) >= 16.0);
    }

    #[test]
    fn a_label_that_fits_the_gutter_is_left_exactly_as_it_was() {
        let name = "enrich / deplete";
        let svg = Figure::new(region())
            .label_width(100.0)
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
        assert_eq!(clips.len(), 2);
        assert!(
            (clips[0] - 16.0).abs() < 1e-9,
            "the strip this track asked for: {clips:?}"
        );
        assert!(
            (clips[1] - layout.plot_x).abs() < 1e-9,
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
        assert!(scaled.to_svg().contains(r#"font-size="27""#));
    }
}
