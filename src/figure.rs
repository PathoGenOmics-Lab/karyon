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
//! track is drawn inside a clip over its band and its axis strip, which is what
//! makes that overhang free: a track draws the whole of a partly visible thing
//! and the clip decides how much shows. The label is the exception, drawn by
//! the figure outside the clip and to the left of the axis strip, so names line
//! up whether or not a track drew an axis.

use std::fs;
use std::io;
use std::path::Path;

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{Anchor, SvgWriter};
use crate::theme::Theme;
use crate::track::{DrawContext, Rect, Track};

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
            top: 10.0,
            right: 16.0,
            bottom: 10.0,
            left: 12.0,
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
    label_width: f64,
    track_gap: f64,
    show_region_label: bool,
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
            label_width: 84.0,
            track_gap: 10.0,
            show_region_label: true,
        }
    }

    /// Sets the image width in pixels.
    ///
    /// Widths that would leave no plotting area are raised to the smallest one
    /// that does, so a figure is always renderable.
    pub fn width(mut self, width: f64) -> Self {
        let floor = self.margin.left + self.margin.right + self.label_width + 50.0;
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

    /// Replaces the theme.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme;
        self
    }

    /// Replaces the margins.
    pub fn margin(mut self, margin: Margin) -> Self {
        self.margin = margin;
        self
    }

    /// Sets the width of the left gutter holding track labels.
    ///
    /// The gutter is only reserved when at least one track returns a label, so
    /// a figure of unlabelled tracks uses the full width whatever this says.
    pub fn label_width(mut self, width: f64) -> Self {
        self.label_width = width.max(0.0);
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
        let layout = self.layout();
        let mut svg = SvgWriter::with_id_prefix(prefix);

        if let Some(title) = &self.title {
            svg.text_bold(
                self.margin.left,
                layout.header_baseline,
                title,
                &self.theme.foreground,
                self.theme.title_font_size,
                Anchor::Start,
            );
        }
        if self.show_region_label {
            svg.text(
                self.width - self.margin.right,
                layout.header_baseline,
                &self.region.to_string(),
                &self.theme.muted,
                self.theme.font_size,
                Anchor::End,
            );
        }

        let mut y = self.margin.top + layout.header_height;
        for (track, height) in self.tracks.iter().zip(&layout.track_heights) {
            let band = Rect {
                x: layout.plot_x,
                y,
                w: layout.plot_width,
                h: *height,
            };

            let axis = Rect {
                x: band.x - layout.axis_width,
                y,
                w: layout.axis_width,
                h: *height,
            };

            if let Some(label) = track.label() {
                // Labels sit to the left of the value axes, so a track with an
                // axis and one without still line their names up.
                svg.text(
                    axis.x - 8.0,
                    band.mid_y() + self.theme.label_font_size * 0.35,
                    label,
                    &self.theme.muted,
                    self.theme.label_font_size,
                    Anchor::End,
                );
            }

            svg.begin_clip(axis.x, band.y, axis.w + band.w, band.h);
            let mut ctx = DrawContext {
                svg: &mut svg,
                scale: &layout.scale,
                theme: &self.theme,
                band,
                axis,
                region: &self.region,
            };
            track.draw(&mut ctx);
            svg.end_group();

            y += height + self.track_gap;
        }

        svg.finish(
            self.width,
            layout.total_height,
            &self.theme.background,
            &self.theme.font_family,
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
        let has_header = self.title.is_some() || self.show_region_label;
        let header_height = if has_header {
            self.theme.title_font_size + 10.0
        } else {
            0.0
        };
        let header_baseline = self.margin.top + self.theme.title_font_size;

        let gutter = if self.tracks.iter().any(|t| t.label().is_some()) {
            self.label_width
        } else {
            0.0
        };
        // The widest axis any track asks for, reserved for all of them, so that
        // every plotting area still starts at the same x.
        let axis_width = self
            .tracks
            .iter()
            .map(|t| t.y_axis_width(&self.theme).max(0.0))
            .fold(0.0f64, f64::max);
        let plot_x = self.margin.left + gutter + axis_width;
        let plot_width = (self.width - plot_x - self.margin.right).max(1.0);
        let scale = Scale::new(&self.region, plot_x, plot_width);

        let track_heights: Vec<f64> = self
            .tracks
            .iter()
            .map(|t| t.height(&scale).max(1.0))
            .collect();
        let content_height: f64 = track_heights.iter().sum::<f64>()
            + self.track_gap * (self.tracks.len().saturating_sub(1)) as f64;

        Layout {
            scale,
            plot_x,
            axis_width,
            plot_width,
            header_height,
            header_baseline,
            track_heights,
            total_height: self.margin.top + header_height + content_height + self.margin.bottom,
        }
    }
}

impl crate::rings::Drawing for Figure {
    fn dimensions(&self) -> (f64, f64) {
        Figure::dimensions(self)
    }

    fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        Figure::to_svg_with_id_prefix(self, prefix)
    }
}

struct Layout {
    scale: Scale,
    plot_x: f64,
    axis_width: f64,
    plot_width: f64,
    header_height: f64,
    header_baseline: f64,
    track_heights: Vec<f64>,
    total_height: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::track::{AxisTrack, CoverageTrack};

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
        assert_eq!(h2 - h1, 20.0 + 10.0);
    }

    #[test]
    fn the_label_gutter_is_only_reserved_when_a_track_wants_it() {
        let bare = Figure::new(region()).push(AxisTrack::new());
        let labelled = Figure::new(region()).push(AxisTrack::new().label("pos"));
        assert_eq!(bare.layout().plot_x, 12.0);
        assert_eq!(labelled.layout().plot_x, 12.0 + 84.0);
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
        assert_eq!(layout.plot_x, 12.0 + figure.label_width + layout.axis_width);
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
        let svg = Figure::new(region()).show_region_label(false).to_svg();
        assert!(!svg.contains("chr1:1-1000"));
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
}
