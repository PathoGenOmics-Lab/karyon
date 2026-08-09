//! Geographic sample maps rendered with the same visual contracts as the
//! genomic and phylogenetic drawings.
//!
//! A map is not a genomic track. Latitude and longitude cannot share a scale
//! with bases, and pretending that they can would make the alignment between
//! neighbouring tracks meaningless. [`Map`] therefore implements
//! [`Drawing`] directly, like [`Rings`](crate::Rings), and can sit beside other
//! complete drawings on a [`Panels`](crate::Panels) sheet.
//!
//! # Coordinates are data
//!
//! [`GeoLocation`] keeps the coordinates it was given. A finite latitude
//! outside -90 to 90 degrees or a longitude outside -180 to 180 degrees is not
//! clamped into a plausible place: the location is counted as invalid and the
//! rendered map says how many inputs it could not draw. An orthographic map
//! likewise reports locations on the hidden hemisphere.
//!
//! # The base map is local and deterministic
//!
//! The simplified land polygons are derived from Natural Earth 1:110m land
//! data, which Natural Earth places in the public domain. They are compiled
//! into the crate; rendering performs no network request and writes no raster
//! image into the SVG.

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::fs;
use std::io;
use std::path::Path;
use std::sync::OnceLock;

use crate::rings::Drawing;
use crate::style::{RenderProfile, Symbol};
use crate::svg::{fit_text, num, Anchor, SvgWriter};
use crate::theme::{mix, Theme};
use crate::track::TreeShape;
use crate::tree::{AnnotationValue, Placement, TimeDirection, Tree};

const WORLD_110M: &str = include_str!("data/world_110m.txt");

mod draw;
mod phylo;
mod project;

#[cfg(test)]
mod tests;

pub use self::phylo::{PhyloConnector, PhyloMap};

use self::draw::*;
use self::project::*;

/// One point on Earth in decimal degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoPosition {
    latitude: f64,
    longitude: f64,
}

impl GeoPosition {
    /// Stores a latitude and longitude in decimal degrees.
    pub fn new(latitude: f64, longitude: f64) -> Self {
        GeoPosition {
            latitude,
            longitude,
        }
    }

    /// Latitude in decimal degrees.
    pub fn latitude(self) -> f64 {
        self.latitude
    }

    /// Longitude in decimal degrees.
    pub fn longitude(self) -> f64 {
        self.longitude
    }

    /// Whether both coordinates are finite and inside their geographic bounds.
    pub fn is_valid(self) -> bool {
        self.latitude.is_finite()
            && self.longitude.is_finite()
            && (-90.0..=90.0).contains(&self.latitude)
            && (-180.0..=180.0).contains(&self.longitude)
    }

    fn radians(self) -> (f64, f64) {
        (self.latitude.to_radians(), self.longitude.to_radians())
    }
}

/// Geographic projection used by [`Map`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum GeoProjection {
    /// Whole world on a two-to-one plate carrée rectangle.
    #[default]
    Equirectangular,
    /// Whole world with conformal local angles and polar latitudes clipped.
    Mercator,
    /// One visible hemisphere projected onto a circle.
    Orthographic {
        /// Latitude at the centre of the globe, in decimal degrees.
        center_latitude: f64,
        /// Longitude at the centre of the globe, in decimal degrees.
        center_longitude: f64,
    },
}

impl GeoProjection {
    /// A globe centred on `latitude`, `longitude`.
    pub fn orthographic(latitude: f64, longitude: f64) -> Self {
        GeoProjection::Orthographic {
            center_latitude: if latitude.is_finite() {
                latitude.clamp(-90.0, 90.0)
            } else {
                0.0
            },
            center_longitude: if longitude.is_finite() {
                longitude.clamp(-180.0, 180.0)
            } else {
                0.0
            },
        }
    }

    fn is_round(self) -> bool {
        matches!(self, GeoProjection::Orthographic { .. })
    }
}

/// One named sample location on a [`Map`].
#[derive(Debug, Clone, PartialEq)]
pub struct GeoLocation {
    name: String,
    position: GeoPosition,
    category: Option<String>,
    value: Option<f64>,
    count: u64,
}

impl GeoLocation {
    /// Creates a named location at `latitude`, `longitude` decimal degrees.
    pub fn new(name: impl Into<String>, latitude: f64, longitude: f64) -> Self {
        GeoLocation {
            name: name.into(),
            position: GeoPosition::new(latitude, longitude),
            category: None,
            value: None,
            count: 1,
        }
    }

    /// Assigns a categorical group used for colour and marker shape.
    pub fn category(mut self, category: impl Into<String>) -> Self {
        self.category = Some(category.into());
        self
    }

    /// Attaches one finite quantitative value to the location.
    pub fn value(mut self, value: f64) -> Self {
        self.value = value.is_finite().then_some(value);
        self
    }

    /// Sets how many observations the location represents.
    pub fn count(mut self, count: u64) -> Self {
        self.count = count.max(1);
        self
    }

    /// Location name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Geographic position.
    pub fn position(&self) -> GeoPosition {
        self.position
    }

    /// Optional categorical group.
    pub fn category_name(&self) -> Option<&str> {
        self.category.as_deref()
    }

    /// Optional quantitative value.
    pub fn numeric_value(&self) -> Option<f64> {
        self.value
    }

    /// Number of observations represented by this mark.
    pub fn observation_count(&self) -> u64 {
        self.count
    }

    fn title(&self) -> String {
        let mut parts = vec![
            self.name.clone(),
            format!("latitude {}", data_number(self.position.latitude)),
            format!("longitude {}", data_number(self.position.longitude)),
        ];
        if self.count > 1 {
            parts.push(format!("{} observations", self.count));
        }
        if let Some(category) = &self.category {
            parts.push(format!("category {category}"));
        }
        if let Some(value) = self.value {
            parts.push(format!("value {}", data_number(value)));
        }
        parts.join("; ")
    }
}

/// A connection between two named [`GeoLocation`] entries.
#[derive(Debug, Clone, PartialEq)]
pub struct GeoFlow {
    from: String,
    to: String,
    category: Option<String>,
    weight: f64,
    directed: bool,
}

impl GeoFlow {
    /// Connects the uniquely named locations `from` and `to`.
    pub fn new(from: impl Into<String>, to: impl Into<String>) -> Self {
        GeoFlow {
            from: from.into(),
            to: to.into(),
            category: None,
            weight: 1.0,
            directed: true,
        }
    }

    /// Assigns a categorical group used for colour.
    pub fn category(mut self, category: impl Into<String>) -> Self {
        self.category = Some(category.into());
        self
    }

    /// Sets a positive visual weight for the connection.
    pub fn weight(mut self, weight: f64) -> Self {
        self.weight = if weight.is_finite() {
            weight.max(0.0)
        } else {
            1.0
        };
        self
    }

    /// Draws the connection without a directional arrowhead.
    pub fn undirected(mut self) -> Self {
        self.directed = false;
        self
    }

    /// Source location name.
    pub fn source(&self) -> &str {
        &self.from
    }

    /// Target location name.
    pub fn target(&self) -> &str {
        &self.to
    }

    /// Whether the connection carries a directional arrowhead.
    pub fn is_directed(&self) -> bool {
        self.directed
    }

    /// Visual weight.
    pub fn flow_weight(&self) -> f64 {
        self.weight
    }

    fn title(&self) -> String {
        let noun = if self.directed { "flow" } else { "link" };
        let mut title = format!(
            "{noun}, source {}, target {}, weight {}",
            self.from,
            self.to,
            data_number(self.weight)
        );
        if let Some(category) = &self.category {
            title.push_str(&format!(", category {category}"));
        }
        title
    }
}

/// A standalone geographic drawing.
///
/// ```
/// use karyon::{GeoFlow, GeoLocation, GeoProjection, Map};
///
/// let svg = Map::new()
///     .title("Surveillance sites")
///     .projection(GeoProjection::orthographic(18.0, 0.0))
///     .push(GeoLocation::new("Lima", -12.046, -77.043).category("South America"))
///     .push(GeoLocation::new("Madrid", 40.417, -3.704).category("Europe"))
///     .push_flow(GeoFlow::new("Lima", "Madrid").weight(4.0))
///     .to_svg();
///
/// assert!(svg.contains("latitude -12.046"));
/// assert!(svg.contains("source Lima, target Madrid"));
/// ```
pub struct Map {
    locations: Vec<GeoLocation>,
    flows: Vec<GeoFlow>,
    projection: GeoProjection,
    width: f64,
    height: f64,
    margin: f64,
    theme: Theme,
    title: Option<String>,
    subtitle: Option<String>,
    description: Option<String>,
    show_graticule: bool,
    show_labels: bool,
    show_legend: bool,
    visual_scale: f64,
}

impl Map {
    /// An empty world map.
    pub fn new() -> Self {
        Map {
            locations: Vec::new(),
            flows: Vec::new(),
            projection: GeoProjection::Equirectangular,
            width: 760.0,
            height: 480.0,
            margin: 18.0,
            theme: Theme::light(),
            title: None,
            subtitle: None,
            description: None,
            show_graticule: true,
            show_labels: false,
            show_legend: true,
            visual_scale: 1.0,
        }
    }

    /// Changes the geographic projection.
    pub fn projection(mut self, projection: GeoProjection) -> Self {
        self.projection = projection;
        self
    }

    /// Sets the output width in pixels.
    pub fn width(mut self, width: f64) -> Self {
        self.width = finite_at_least(width, 180.0, 760.0);
        self
    }

    /// Sets the output height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = finite_at_least(height, 160.0, 480.0);
        self
    }

    /// Sets whitespace around the drawing.
    pub fn margin(mut self, margin: f64) -> Self {
        self.margin = finite_at_least(margin, 0.0, 18.0);
        self
    }

    /// Replaces the visual theme.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme;
        self
    }

    /// Applies a named output profile.
    pub fn profile(mut self, profile: RenderProfile) -> Self {
        self.theme = if profile.is_dark() {
            Theme::dark()
        } else {
            Theme::light()
        };
        self.visual_scale = profile.visual_scale();
        self
    }

    /// Scales typography and marks without changing the map extent.
    pub fn visual_scale(mut self, scale: f64) -> Self {
        self.visual_scale = finite_at_least(scale, 0.25, 1.0);
        self
    }

    /// Sets the visible title and document name.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Sets a quieter line below the title.
    pub fn subtitle(mut self, subtitle: impl Into<String>) -> Self {
        self.subtitle = Some(subtitle.into());
        self
    }

    /// Sets the accessible description of the map.
    pub fn description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }

    /// Draws or hides the latitude-longitude graticule.
    pub fn show_graticule(mut self, show: bool) -> Self {
        self.show_graticule = show;
        self
    }

    /// Draws or hides location names beside their marks.
    pub fn show_labels(mut self, show: bool) -> Self {
        self.show_labels = show;
        self
    }

    /// Draws or hides the categorical legend.
    pub fn show_legend(mut self, show: bool) -> Self {
        self.show_legend = show;
        self
    }

    /// Adds one location in input order.
    pub fn push(mut self, location: GeoLocation) -> Self {
        self.locations.push(location);
        self
    }

    /// Adds several locations in input order.
    pub fn extend(mut self, locations: impl IntoIterator<Item = GeoLocation>) -> Self {
        self.locations.extend(locations);
        self
    }

    /// Adds one geographic connection.
    pub fn push_flow(mut self, flow: GeoFlow) -> Self {
        self.flows.push(flow);
        self
    }

    /// Locations in input order, including invalid entries.
    pub fn locations(&self) -> &[GeoLocation] {
        &self.locations
    }

    /// Connections in input order.
    pub fn flows(&self) -> &[GeoFlow] {
        &self.flows
    }

    /// Number of locations whose coordinates are not drawable.
    pub fn invalid_location_count(&self) -> usize {
        self.locations
            .iter()
            .filter(|location| !location.position.is_valid())
            .count()
    }

    /// Number of flows whose endpoint name is missing or ambiguous.
    pub fn unresolved_flow_count(&self) -> usize {
        let names = location_names(&self.locations);
        self.flows
            .iter()
            .filter(|flow| {
                !unique_endpoint(&names, &flow.from) || !unique_endpoint(&names, &flow.to)
            })
            .count()
    }

    /// Output width and height.
    pub fn dimensions(&self) -> (f64, f64) {
        (self.width, self.height)
    }

    /// Renders the map to a standalone SVG document.
    pub fn to_svg(&self) -> String {
        self.to_svg_with_id_prefix("")
    }

    /// Renders the map with generated SVG ids prefixed for nesting.
    pub fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        let theme = self.theme.clone().scaled(self.visual_scale);
        let mut svg = SvgWriter::with_id_prefix(prefix);
        svg.describe(&self.document_name(), &self.document_description());
        let layout = self.layout(&theme);

        if let Some(title) = &self.title {
            svg.text_bold(
                self.margin,
                self.margin + theme.title_font_size,
                title,
                &theme.foreground,
                theme.title_font_size,
                Anchor::Start,
            );
        }
        if let Some(subtitle) = &self.subtitle {
            let y = self.margin
                + if self.title.is_some() {
                    theme.title_font_size + theme.font_size + 5.0
                } else {
                    theme.font_size
                };
            svg.text(
                self.margin,
                y,
                subtitle,
                &theme.muted,
                theme.font_size,
                Anchor::Start,
            );
        }

        let categories = category_indices(&self.locations, &self.flows);
        draw_base_map(
            &mut svg,
            layout.map,
            self.projection,
            &theme,
            self.show_graticule,
        );
        let state = draw_flows(
            &mut svg,
            layout.map,
            self.projection,
            &theme,
            &categories,
            &self.locations,
            &self.flows,
        );
        let hidden = draw_locations(
            &mut svg,
            layout.map,
            self.projection,
            &theme,
            &categories,
            &self.locations,
            self.show_labels,
        );
        if self.show_legend {
            draw_legend(&mut svg, layout, &theme, &categories);
        }
        draw_warning(
            &mut svg,
            layout.warning_y,
            &theme,
            self.invalid_location_count(),
            hidden,
            state.unresolved,
            state.hidden,
        );
        svg.finish(
            self.width,
            self.height,
            &theme.background,
            &theme.font_family,
        )
    }

    /// Saves the SVG document to `path`.
    pub fn save_svg(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, self.to_svg())
    }

    fn document_name(&self) -> String {
        match (&self.title, &self.subtitle) {
            (Some(title), Some(subtitle)) => format!("{title}, {subtitle}"),
            (Some(title), None) => title.clone(),
            (None, Some(subtitle)) => subtitle.clone(),
            (None, None) => "Geographic sample map".to_string(),
        }
    }

    fn document_description(&self) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        format!(
            "A karyon geographic map with {} locations and {} flows; {} invalid locations and {} unresolved flows.",
            self.locations.len(),
            self.flows.len(),
            self.invalid_location_count(),
            self.unresolved_flow_count()
        )
    }

    fn layout(&self, theme: &Theme) -> MapLayout {
        let title_room = match (self.title.is_some(), self.subtitle.is_some()) {
            (true, true) => theme.title_font_size + theme.font_size + 12.0,
            (true, false) => theme.title_font_size + 8.0,
            (false, true) => theme.font_size + 8.0,
            (false, false) => 0.0,
        };
        let warning_room = if self.invalid_location_count() > 0
            || self.unresolved_flow_count() > 0
            || self.projection.is_round()
        {
            theme.font_size + 8.0
        } else {
            0.0
        };
        let categories = category_indices(&self.locations, &self.flows);
        let legend_width = if self.show_legend && !categories.is_empty() {
            154.0 * self.visual_scale
        } else {
            0.0
        };
        let available = MapRect {
            x: self.margin,
            y: self.margin + title_room,
            w: (self.width - self.margin * 2.0 - legend_width).max(20.0),
            h: (self.height - self.margin * 2.0 - title_room - warning_room).max(20.0),
        };
        let map = if self.projection.is_round() {
            let side = available.w.min(available.h);
            MapRect {
                x: available.x + (available.w - side) / 2.0,
                y: available.y + (available.h - side) / 2.0,
                w: side,
                h: side,
            }
        } else {
            let ratio = match self.projection {
                GeoProjection::Equirectangular => 2.0,
                GeoProjection::Mercator => 1.0,
                GeoProjection::Orthographic { .. } => unreachable!(),
            };
            let (w, h) = if available.w / available.h > ratio {
                (available.h * ratio, available.h)
            } else {
                (available.w, available.w / ratio)
            };
            MapRect {
                x: available.x + (available.w - w) / 2.0,
                y: available.y + (available.h - h) / 2.0,
                w,
                h,
            }
        };
        MapLayout {
            map,
            legend_x: self.width - self.margin - legend_width + 14.0,
            legend_y: available.y + 4.0,
            legend_width: (legend_width - 18.0).max(0.0),
            warning_y: self.height - self.margin,
        }
    }
}

impl Default for Map {
    fn default() -> Self {
        Self::new()
    }
}

impl Drawing for Map {
    fn dimensions(&self) -> (f64, f64) {
        Map::dimensions(self)
    }

    fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        Map::to_svg_with_id_prefix(self, prefix)
    }
}
