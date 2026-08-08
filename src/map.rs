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

const WORLD_110M: &str = include_str!("map/data/world_110m.txt");

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

/// How terminal taxa are connected to their geographic location in a
/// [`PhyloMap`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PhyloConnector {
    /// One connector per location, with its exact terminal count in a tooltip.
    #[default]
    Aggregated,
    /// One connector per mapped terminal taxon.
    Individual,
    /// Draw the tree and map marks without joining them.
    None,
}

/// A circular phylogeny surrounding a geographic map.
///
/// The tree supplies one annotation per terminal taxon, such as `country` or
/// `location`, and the coordinate table supplies the matching named places.
/// Karyon does not infer ancestral locations, migrations or transmission.
/// Branches only inherit an ancestor annotation that already exists in the
/// source tree.
///
/// ```
/// use karyon::{GeoLocation, GeoProjection, PhyloMap, Tree};
///
/// let tree = Tree::parse_annotated_newick(
///     "(sample_A[&country=Peru]:1,sample_B[&country=Spain]:1);",
/// )?;
/// let svg = PhyloMap::new(tree)
///     .location_by("country")
///     .projection(GeoProjection::orthographic(15.0, -10.0))
///     .coordinates([
///         GeoLocation::new("Peru", -9.19, -75.015),
///         GeoLocation::new("Spain", 40.464, -3.749),
///     ])
///     .to_svg();
///
/// assert!(svg.contains("Peru; 1 mapped tip"));
/// # Ok::<(), karyon::Error>(())
/// ```
pub struct PhyloMap {
    tree: Tree,
    locations: Vec<GeoLocation>,
    location_key: String,
    projection: GeoProjection,
    connector: PhyloConnector,
    shape: TreeShape,
    time: Option<PhyloTime>,
    diameter: f64,
    margin: f64,
    start_degrees: f64,
    sweep_degrees: f64,
    theme: Theme,
    visual_scale: f64,
    title: Option<String>,
    subtitle: Option<String>,
    description: Option<String>,
    show_graticule: bool,
    show_legend: bool,
    show_nodes: bool,
    show_tip_labels: bool,
    show_time_guides: bool,
}

#[derive(Debug, Clone)]
struct PhyloTime {
    key: String,
    direction: TimeDirection,
    unit: Option<String>,
}

impl PhyloMap {
    /// Starts a circular phylogeographic composition from `tree`.
    pub fn new(tree: Tree) -> Self {
        PhyloMap {
            tree,
            locations: Vec::new(),
            location_key: "location".to_string(),
            projection: GeoProjection::orthographic(15.0, -5.0),
            connector: PhyloConnector::Aggregated,
            shape: TreeShape::Phylogram,
            time: None,
            diameter: 760.0,
            margin: 18.0,
            start_degrees: -95.0,
            sweep_degrees: 350.0,
            theme: Theme::light(),
            visual_scale: 1.0,
            title: None,
            subtitle: None,
            description: None,
            show_graticule: true,
            show_legend: true,
            show_nodes: true,
            show_tip_labels: false,
            show_time_guides: true,
        }
    }

    /// Selects the terminal annotation that names a coordinate-table entry.
    pub fn location_by(mut self, key: impl Into<String>) -> Self {
        self.location_key = key.into();
        self
    }

    /// Adds the named coordinate table used by terminal annotations.
    pub fn coordinates(mut self, locations: impl IntoIterator<Item = GeoLocation>) -> Self {
        self.locations.extend(locations);
        self
    }

    /// Adds one named coordinate-table entry.
    pub fn coordinate(mut self, location: GeoLocation) -> Self {
        self.locations.push(location);
        self
    }

    /// Changes the central map projection.
    pub fn projection(mut self, projection: GeoProjection) -> Self {
        self.projection = projection;
        self
    }

    /// Chooses individual, aggregated or absent tree-to-map connectors.
    pub fn connector(mut self, connector: PhyloConnector) -> Self {
        self.connector = connector;
        self
    }

    /// Chooses a branch-length phylogram or branch-count cladogram.
    pub fn shape(mut self, shape: TreeShape) -> Self {
        self.shape = shape;
        self
    }

    /// Places nodes by the finite numeric annotation `key` when every tip has it.
    pub fn time(mut self, key: impl Into<String>) -> Self {
        self.time = Some(PhyloTime {
            key: key.into(),
            direction: TimeDirection::Increasing,
            unit: None,
        });
        self
    }

    /// Changes the direction of the configured time values.
    pub fn time_direction(mut self, direction: TimeDirection) -> Self {
        if let Some(time) = &mut self.time {
            time.direction = direction;
        }
        self
    }

    /// Adds a unit to concentric time-guide labels.
    pub fn time_unit(mut self, unit: impl Into<String>) -> Self {
        if let Some(time) = &mut self.time {
            time.unit = Some(unit.into());
        }
        self
    }

    /// Sets the diameter of the tree-bearing circle in pixels.
    pub fn diameter(mut self, diameter: f64) -> Self {
        self.diameter = finite_at_least(diameter, 320.0, 760.0);
        self
    }

    /// Sets whitespace around the drawing.
    pub fn margin(mut self, margin: f64) -> Self {
        self.margin = finite_at_least(margin, 0.0, 18.0);
        self
    }

    /// Rotates the first terminal taxon in clockwise degrees.
    pub fn radial_start(mut self, degrees: f64) -> Self {
        if degrees.is_finite() {
            self.start_degrees = degrees;
        }
        self
    }

    /// Sets the clockwise angular span from 120 to 360 degrees.
    pub fn radial_sweep(mut self, degrees: f64) -> Self {
        self.sweep_degrees = if degrees.is_finite() {
            degrees.clamp(120.0, 360.0)
        } else {
            350.0
        };
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

    /// Scales typography and marks without changing the tree diameter.
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

    /// Sets the accessible description of the composition.
    pub fn description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }

    /// Draws or hides the central graticule.
    pub fn show_graticule(mut self, show: bool) -> Self {
        self.show_graticule = show;
        self
    }

    /// Draws or hides the location legend.
    pub fn show_legend(mut self, show: bool) -> Self {
        self.show_legend = show;
        self
    }

    /// Draws or hides internal and terminal node marks.
    pub fn show_nodes(mut self, show: bool) -> Self {
        self.show_nodes = show;
        self
    }

    /// Draws or hides terminal names beside the inner circumference.
    pub fn show_tip_labels(mut self, show: bool) -> Self {
        self.show_tip_labels = show;
        self
    }

    /// Draws or hides concentric guides when a time layout succeeds.
    pub fn show_time_guides(mut self, show: bool) -> Self {
        self.show_time_guides = show;
        self
    }

    /// The source tree, unchanged by rendering.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Coordinate-table entries in input order.
    pub fn locations(&self) -> &[GeoLocation] {
        &self.locations
    }

    /// Number of tips with no unique, valid coordinate-table match.
    pub fn unlocated_tip_count(&self) -> usize {
        let names = location_names(&self.locations);
        self.tree
            .leaves()
            .into_iter()
            .filter(|tip| {
                resolve_tip_location(
                    &self.tree,
                    *tip,
                    &self.location_key,
                    &self.locations,
                    &names,
                )
                .is_none()
            })
            .count()
    }

    /// Width and height of the rendered document.
    pub fn dimensions(&self) -> (f64, f64) {
        let theme = self.theme.clone().scaled(self.visual_scale);
        let scene = PhyloScene::new(&self.tree, self.shape, self.time.as_ref());
        let categories = phylo_categories(&self.tree, &scene, &self.location_key);
        let layout = self.layout(&theme, !categories.is_empty());
        (layout.width, layout.height)
    }

    /// Renders the complete composition to a standalone SVG document.
    pub fn to_svg(&self) -> String {
        self.to_svg_with_id_prefix("")
    }

    /// Renders with generated SVG ids prefixed for nesting in [`Panels`](crate::Panels).
    pub fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        let theme = self.theme.clone().scaled(self.visual_scale);
        let scene = PhyloScene::new(&self.tree, self.shape, self.time.as_ref());
        let categories = phylo_categories(&self.tree, &scene, &self.location_key);
        let layout = self.layout(&theme, !categories.is_empty());
        let names = location_names(&self.locations);
        let mapped = mapped_tip_groups(
            &self.tree,
            &scene,
            &self.location_key,
            &self.locations,
            &names,
            self.projection,
            layout.map,
        );
        let unlocated = self.unlocated_tip_count();
        let hidden = mapped
            .iter()
            .filter(|group| group.projected.is_none())
            .count();
        let mut svg = SvgWriter::with_id_prefix(prefix);
        svg.describe(
            &self.document_name(),
            &self.document_description(&scene, unlocated, hidden),
        );
        draw_phylo_title(self, &mut svg, &theme);
        draw_base_map(
            &mut svg,
            layout.map,
            self.projection,
            &theme,
            self.show_graticule,
        );
        if self.show_time_guides && scene.temporal {
            draw_phylo_time_guides(self, &mut svg, &theme, &scene, layout);
        }
        draw_phylo_connectors(self, &mut svg, &theme, &scene, layout, &categories, &mapped);
        draw_phylo_tree(self, &mut svg, &theme, &scene, layout, &categories);
        draw_phylo_locations(&mut svg, &theme, &categories, &mapped);
        if self.show_legend {
            draw_phylo_legend(&mut svg, &theme, layout, &categories, &mapped);
        }
        draw_phylo_warning(
            &mut svg,
            &theme,
            layout,
            unlocated,
            hidden,
            scene.time_failed,
        );
        svg.finish(
            layout.width,
            layout.height,
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
            (None, None) => "Circular phylogeographic map".to_string(),
        }
    }

    fn document_description(&self, scene: &PhyloScene, unlocated: usize, hidden: usize) -> String {
        if let Some(description) = &self.description {
            return description.clone();
        }
        format!(
            "A karyon circular phylogeny of {} tips around a geographic map, connected by {}; {} tips lack a drawable location and {} mapped locations are outside the projection.{}",
            scene.leaves.len(),
            match self.connector {
                PhyloConnector::Aggregated => "one line per location",
                PhyloConnector::Individual => "one line per mapped tip",
                PhyloConnector::None => "no lines",
            },
            unlocated,
            hidden,
            if scene.time_failed {
                " The requested time layout was unavailable."
            } else {
                ""
            }
        )
    }

    fn layout(&self, theme: &Theme, has_categories: bool) -> PhyloLayout {
        let title_room = match (self.title.is_some(), self.subtitle.is_some()) {
            (true, true) => theme.title_font_size + theme.font_size + 12.0,
            (true, false) => theme.title_font_size + 8.0,
            (false, true) => theme.font_size + 8.0,
            (false, false) => 0.0,
        };
        let warning_room = theme.font_size + 8.0;
        let legend_width = if self.show_legend && has_categories {
            172.0 * self.visual_scale
        } else {
            0.0
        };
        let width = self.diameter + self.margin * 2.0 + legend_width;
        let height = self.diameter + self.margin * 2.0 + title_room + warning_room;
        let cx = self.margin + self.diameter / 2.0;
        let cy = self.margin + title_room + self.diameter / 2.0;
        let map_radius = self.diameter * 0.235;
        let tree_inner = map_radius + self.diameter * 0.075;
        let tree_outer = self.diameter / 2.0 - theme.font_size * 1.7;
        PhyloLayout {
            width,
            height,
            cx,
            cy,
            map: MapRect {
                x: cx - map_radius,
                y: cy - map_radius,
                w: map_radius * 2.0,
                h: map_radius * 2.0,
            },
            tree_inner,
            tree_outer: tree_outer.max(tree_inner + 10.0),
            start: self.start_degrees.to_radians(),
            sweep: self.sweep_degrees.to_radians(),
            legend_x: self.margin + self.diameter + 14.0,
            legend_y: self.margin + title_room + 12.0,
            legend_width: (legend_width - 22.0).max(0.0),
            warning_y: height - self.margin,
        }
    }
}

impl Drawing for PhyloMap {
    fn dimensions(&self) -> (f64, f64) {
        PhyloMap::dimensions(self)
    }

    fn to_svg_with_id_prefix(&self, prefix: &str) -> String {
        PhyloMap::to_svg_with_id_prefix(self, prefix)
    }
}

#[derive(Debug, Clone, Copy)]
struct PhyloLayout {
    width: f64,
    height: f64,
    cx: f64,
    cy: f64,
    map: MapRect,
    tree_inner: f64,
    tree_outer: f64,
    start: f64,
    sweep: f64,
    legend_x: f64,
    legend_y: f64,
    legend_width: f64,
    warning_y: f64,
}

struct PhyloScene {
    placements: Vec<Placement>,
    leaves: Vec<usize>,
    minimum: f64,
    maximum: f64,
    temporal: bool,
    time_failed: bool,
    direction: TimeDirection,
}

impl PhyloScene {
    fn new(tree: &Tree, shape: TreeShape, time: Option<&PhyloTime>) -> Self {
        let timed = time.and_then(|time| tree.time_layout(&time.key, time.direction));
        let temporal = timed.is_some();
        let time_failed = time.is_some() && timed.is_none();
        let placements = timed.unwrap_or_else(|| tree.layout(shape == TreeShape::Cladogram));
        let (minimum, maximum) = placements
            .iter()
            .map(|placement| placement.depth)
            .filter(|value| value.is_finite())
            .fold((f64::MAX, f64::MIN), |(minimum, maximum), value| {
                (minimum.min(value), maximum.max(value))
            });
        let (minimum, maximum) = if minimum.is_finite() && maximum.is_finite() {
            (minimum, maximum)
        } else {
            (0.0, 1.0)
        };
        PhyloScene {
            placements,
            leaves: tree.leaves(),
            minimum,
            maximum,
            temporal,
            time_failed,
            direction: time.map_or(TimeDirection::Increasing, |time| time.direction),
        }
    }

    fn fraction(&self, value: f64) -> f64 {
        let span = self.maximum - self.minimum;
        if span <= 0.0 {
            return 0.0;
        }
        let fraction = match self.direction {
            TimeDirection::Increasing => (value - self.minimum) / span,
            TimeDirection::Decreasing if self.temporal => (self.maximum - value) / span,
            TimeDirection::Decreasing => (value - self.minimum) / span,
        };
        fraction.clamp(0.0, 1.0)
    }

    fn radius(&self, layout: PhyloLayout, value: f64) -> f64 {
        layout.tree_outer - self.fraction(value) * (layout.tree_outer - layout.tree_inner)
    }

    fn angle(&self, layout: PhyloLayout, row: f64) -> f64 {
        if self.leaves.len() <= 1 {
            return layout.start + layout.sweep / 2.0;
        }
        let denominator = if layout.sweep >= TAU - 1e-6 {
            self.leaves.len() as f64
        } else {
            (self.leaves.len() - 1) as f64
        };
        layout.start + layout.sweep * row / denominator
    }
}

struct MappedTipGroup<'a> {
    location: &'a GeoLocation,
    tips: Vec<usize>,
    projected: Option<(f64, f64)>,
}

fn resolve_tip_location<'a>(
    tree: &Tree,
    tip: usize,
    key: &str,
    locations: &'a [GeoLocation],
    names: &BTreeMap<&str, Vec<usize>>,
) -> Option<&'a GeoLocation> {
    let value = inherited_tree_annotation(tree, tip, key)?.to_string();
    let indices = names.get(value.as_str())?;
    if indices.len() != 1 {
        return None;
    }
    let location = &locations[indices[0]];
    location.position.is_valid().then_some(location)
}

fn mapped_tip_groups<'a>(
    tree: &Tree,
    scene: &PhyloScene,
    key: &str,
    locations: &'a [GeoLocation],
    names: &BTreeMap<&str, Vec<usize>>,
    projection: GeoProjection,
    area: MapRect,
) -> Vec<MappedTipGroup<'a>> {
    let mut groups: Vec<MappedTipGroup<'a>> = Vec::new();
    let mut by_name: BTreeMap<&str, usize> = BTreeMap::new();
    for tip in &scene.leaves {
        let Some(location) = resolve_tip_location(tree, *tip, key, locations, names) else {
            continue;
        };
        let index = if let Some(index) = by_name.get(location.name.as_str()) {
            *index
        } else {
            let index = groups.len();
            by_name.insert(location.name.as_str(), index);
            groups.push(MappedTipGroup {
                location,
                tips: Vec::new(),
                projected: project(location.position, projection, area),
            });
            index
        };
        groups[index].tips.push(*tip);
    }
    groups
}

fn phylo_categories(tree: &Tree, scene: &PhyloScene, key: &str) -> BTreeMap<String, usize> {
    let mut categories = BTreeMap::new();
    for placement in &scene.placements {
        let Some(value) = inherited_tree_annotation(tree, placement.node, key) else {
            continue;
        };
        let value = value.to_string();
        let next = categories.len();
        categories.entry(value).or_insert(next);
    }
    categories
}

fn inherited_tree_annotation<'a>(
    tree: &'a Tree,
    node: usize,
    key: &str,
) -> Option<&'a AnnotationValue> {
    tree.annotation(node, key).or_else(|| {
        tree.ancestors(node)
            .into_iter()
            .find_map(|ancestor| tree.annotation(ancestor, key))
    })
}

fn phylo_point(layout: PhyloLayout, radius: f64, angle: f64) -> (f64, f64) {
    (
        layout.cx + radius * angle.cos(),
        layout.cy + radius * angle.sin(),
    )
}

fn draw_phylo_title(map: &PhyloMap, svg: &mut SvgWriter, theme: &Theme) {
    if let Some(title) = &map.title {
        svg.text_bold(
            map.margin,
            map.margin + theme.title_font_size,
            title,
            &theme.foreground,
            theme.title_font_size,
            Anchor::Start,
        );
    }
    if let Some(subtitle) = &map.subtitle {
        let y = map.margin
            + if map.title.is_some() {
                theme.title_font_size + theme.font_size + 5.0
            } else {
                theme.font_size
            };
        svg.text(
            map.margin,
            y,
            subtitle,
            &theme.muted,
            theme.font_size,
            Anchor::Start,
        );
    }
}

fn draw_phylo_time_guides(
    map: &PhyloMap,
    svg: &mut SvgWriter,
    theme: &Theme,
    scene: &PhyloScene,
    layout: PhyloLayout,
) {
    let Some(time) = &map.time else {
        return;
    };
    for index in 0..=2 {
        let fraction = index as f64 / 2.0;
        let value = scene.minimum + fraction * (scene.maximum - scene.minimum);
        let radius = scene.radius(layout, value);
        let path = radial_arc(layout, radius, layout.start, layout.start + layout.sweep);
        svg.path_stroked(&path, &theme.rule, theme.tokens.hairline * 0.75);
        let label = match &time.unit {
            Some(unit) => format!("{} {unit}", num(value)),
            None => num(value),
        };
        let angle = layout.start - 0.018;
        let (x, y) = phylo_point(layout, radius, angle);
        svg.text_rotated(
            (x, y),
            upright_tangent(angle),
            &label,
            &theme.muted,
            theme.font_size - 2.0,
            Anchor::End,
        );
    }
}

fn draw_phylo_connectors(
    map: &PhyloMap,
    svg: &mut SvgWriter,
    theme: &Theme,
    scene: &PhyloScene,
    layout: PhyloLayout,
    categories: &BTreeMap<String, usize>,
    groups: &[MappedTipGroup<'_>],
) {
    if map.connector == PhyloConnector::None {
        return;
    }
    for group in groups {
        let Some(map_point) = group.projected else {
            continue;
        };
        let (base, _) = phylo_category_style(Some(&group.location.name), categories, theme);
        let color = mix(base, theme.surface(), 0.35);
        match map.connector {
            PhyloConnector::Aggregated => {
                let angle = circular_mean(
                    group
                        .tips
                        .iter()
                        .map(|tip| scene.angle(layout, scene.placements[*tip].row)),
                );
                let tree_point = phylo_point(layout, layout.tree_inner - 2.0, angle);
                let path = connector_path(map_point, tree_point, layout);
                let title = format!(
                    "{}; {}",
                    group.location.name,
                    plural(group.tips.len(), "mapped tip", "mapped tips")
                );
                svg.begin_titled(&title);
                svg.path_stroked(
                    &path,
                    &color,
                    theme.tokens.stroke
                        * (0.8 + (group.tips.len() as f64).sqrt().clamp(1.0, 5.0) * 0.35),
                );
                svg.end_group();
            }
            PhyloConnector::Individual => {
                for tip in &group.tips {
                    let angle = scene.angle(layout, scene.placements[*tip].row);
                    let tree_point = phylo_point(layout, layout.tree_inner - 2.0, angle);
                    let path = connector_path(map_point, tree_point, layout);
                    let sample = map.tree.nodes()[*tip]
                        .name
                        .as_deref()
                        .unwrap_or("unnamed tip");
                    svg.begin_titled(&format!("{sample}; location {}", group.location.name));
                    svg.path_stroked(&path, &color, theme.tokens.hairline);
                    svg.end_group();
                }
            }
            PhyloConnector::None => {}
        }
    }
}

fn connector_path(start: (f64, f64), end: (f64, f64), layout: PhyloLayout) -> String {
    let midpoint = ((start.0 + end.0) / 2.0, (start.1 + end.1) / 2.0);
    let from_centre = (midpoint.0 - layout.cx, midpoint.1 - layout.cy);
    let length = from_centre.0.hypot(from_centre.1).max(1.0);
    let lift = (end.0 - start.0).hypot(end.1 - start.1) * 0.08;
    let control = (
        midpoint.0 + from_centre.0 / length * lift,
        midpoint.1 + from_centre.1 / length * lift,
    );
    format!(
        "M {} {} Q {} {} {} {}",
        num(start.0),
        num(start.1),
        num(control.0),
        num(control.1),
        num(end.0),
        num(end.1)
    )
}

fn draw_phylo_tree(
    map: &PhyloMap,
    svg: &mut SvgWriter,
    theme: &Theme,
    scene: &PhyloScene,
    layout: PhyloLayout,
    categories: &BTreeMap<String, usize>,
) {
    for placement in &scene.placements {
        let node = &map.tree.nodes()[placement.node];
        let angle = scene.angle(layout, placement.row);
        let radius = scene.radius(layout, placement.depth);
        let value = inherited_tree_annotation(&map.tree, placement.node, &map.location_key)
            .map(ToString::to_string);
        let (color, symbol) = phylo_category_style(value.as_deref(), categories, theme);
        if let Some(parent) = node.parent {
            let parent_radius = scene.radius(layout, scene.placements[parent].depth);
            let from = phylo_point(layout, parent_radius, angle);
            let to = phylo_point(layout, radius, angle);
            let title = phylo_branch_title(map, placement.node);
            if let Some(title) = &title {
                svg.begin_titled(title);
            }
            svg.line(from.0, from.1, to.0, to.1, color, theme.tokens.stroke);
            if title.is_some() {
                svg.end_group();
            }
        }
        if node.children.len() >= 2 {
            let first = scene.angle(layout, scene.placements[node.children[0]].row);
            let last = scene.angle(layout, scene.placements[*node.children.last().unwrap()].row);
            svg.path_stroked(
                &radial_arc(layout, radius, first, last),
                color,
                theme.tokens.stroke,
            );
        }
        if node.is_leaf() {
            let boundary = phylo_point(layout, layout.tree_inner, angle);
            let tip = phylo_point(layout, radius, angle);
            if (radius - layout.tree_inner).abs() > 0.3 {
                svg.line(
                    tip.0,
                    tip.1,
                    boundary.0,
                    boundary.1,
                    color,
                    theme.tokens.hairline,
                );
            }
            if map.show_nodes {
                svg.begin_titled(
                    &phylo_branch_title(map, placement.node).unwrap_or_else(|| "tip".into()),
                );
                svg.symbol_ringed(
                    boundary.0,
                    boundary.1,
                    theme.tokens.marker_radius * 0.85,
                    symbol,
                    color,
                    theme.surface(),
                    theme.tokens.hairline,
                );
                svg.end_group();
            }
            if map.show_tip_labels {
                draw_phylo_tip_label(map, svg, theme, layout, angle, placement.node);
            }
        } else if map.show_nodes {
            let point = phylo_point(layout, radius, angle);
            svg.symbol_ringed(
                point.0,
                point.1,
                theme.tokens.marker_radius * 0.48,
                Symbol::Circle,
                color,
                theme.surface(),
                theme.tokens.hairline,
            );
        }
    }
}

fn phylo_branch_title(map: &PhyloMap, node: usize) -> Option<String> {
    let clade = &map.tree.nodes()[node];
    let mut parts = Vec::new();
    if clade.is_leaf() {
        if let Some(name) = &clade.name {
            if !name.is_empty() {
                parts.push(name.clone());
            }
        }
    }
    if let Some(support) = clade.support.filter(|support| support.is_finite()) {
        parts.push(format!("clade support {}", num(support)));
    }
    if let Some(value) = inherited_tree_annotation(&map.tree, node, &map.location_key) {
        parts.push(format!("{} {value}", map.location_key));
    }
    (!parts.is_empty()).then(|| parts.join("; "))
}

fn draw_phylo_tip_label(
    map: &PhyloMap,
    svg: &mut SvgWriter,
    theme: &Theme,
    layout: PhyloLayout,
    angle: f64,
    node: usize,
) {
    let label = map.tree.nodes()[node].name.as_deref().unwrap_or("");
    if label.is_empty() {
        return;
    }
    let radius = layout.tree_inner - theme.tokens.label_gap - 3.0;
    let (x, y) = phylo_point(layout, radius, angle);
    let degrees = angle.to_degrees().rem_euclid(360.0);
    let right = angle.cos() >= 0.0;
    let (rotation, anchor) = if right {
        (degrees, Anchor::End)
    } else {
        (degrees + 180.0, Anchor::Start)
    };
    svg.text_rotated(
        (x, y + theme.font_size * 0.3),
        rotation,
        label,
        &theme.muted,
        theme.font_size - 1.0,
        anchor,
    );
}

fn draw_phylo_locations(
    svg: &mut SvgWriter,
    theme: &Theme,
    categories: &BTreeMap<String, usize>,
    groups: &[MappedTipGroup<'_>],
) {
    for group in groups {
        let Some((x, y)) = group.projected else {
            continue;
        };
        let (color, symbol) = phylo_category_style(Some(&group.location.name), categories, theme);
        let radius =
            theme.tokens.marker_radius * 1.2 * (group.tips.len() as f64).sqrt().clamp(1.0, 3.6);
        let title = format!(
            "{}; latitude {}; longitude {}; {}",
            group.location.name,
            data_number(group.location.position.latitude),
            data_number(group.location.position.longitude),
            plural(group.tips.len(), "mapped tip", "mapped tips")
        );
        svg.begin_titled(&title);
        svg.symbol_ringed(
            x,
            y,
            radius,
            symbol,
            color,
            theme.surface(),
            theme.tokens.hairline * 1.6,
        );
        svg.end_group();
    }
}

fn draw_phylo_legend(
    svg: &mut SvgWriter,
    theme: &Theme,
    layout: PhyloLayout,
    categories: &BTreeMap<String, usize>,
    groups: &[MappedTipGroup<'_>],
) {
    if categories.is_empty() || layout.legend_width <= 0.0 {
        return;
    }
    svg.text_bold(
        layout.legend_x,
        layout.legend_y + theme.font_size,
        "Locations",
        &theme.foreground,
        theme.font_size,
        Anchor::Start,
    );
    let mut ordered = categories.iter().collect::<Vec<_>>();
    ordered.sort_by_key(|(_, index)| **index);
    for (row, (category, index)) in ordered.into_iter().enumerate() {
        let y = layout.legend_y + theme.font_size + 12.0 + row as f64 * (theme.font_size + 8.0);
        let count = groups
            .iter()
            .find(|group| group.location.name == *category)
            .map_or(0, |group| group.tips.len());
        svg.symbol_ringed(
            layout.legend_x + theme.tokens.marker_radius,
            y - theme.font_size * 0.28,
            theme.tokens.marker_radius,
            theme.symbol(*index),
            theme.color(*index),
            theme.surface(),
            theme.tokens.hairline,
        );
        let full = format!("{category}  n={count}");
        let label = fit_text(
            &full,
            (layout.legend_width - theme.tokens.marker_radius * 3.0).max(1.0),
            theme.font_size - 1.0,
        );
        if label != full {
            svg.begin_titled(&full);
        }
        svg.text(
            layout.legend_x + theme.tokens.marker_radius * 2.5,
            y,
            &label,
            &theme.muted,
            theme.font_size - 1.0,
            Anchor::Start,
        );
        if label != full {
            svg.end_group();
        }
    }
}

fn draw_phylo_warning(
    svg: &mut SvgWriter,
    theme: &Theme,
    layout: PhyloLayout,
    unlocated: usize,
    hidden: usize,
    time_failed: bool,
) {
    let mut parts = Vec::new();
    if unlocated > 0 {
        parts.push(plural(unlocated, "unlocated tip", "unlocated tips"));
    }
    if hidden > 0 {
        parts.push(plural(
            hidden,
            "mapped location outside projection",
            "mapped locations outside projection",
        ));
    }
    if time_failed {
        parts.push("requested time layout unavailable".to_string());
    }
    if parts.is_empty() {
        return;
    }
    svg.text(
        18.0,
        layout.warning_y,
        &parts.join("; "),
        &theme.muted,
        theme.font_size - 1.0,
        Anchor::Start,
    );
}

fn radial_arc(layout: PhyloLayout, radius: f64, start: f64, end: f64) -> String {
    let delta = (end - start).abs();
    let first = phylo_point(layout, radius, start);
    if delta >= TAU - 1e-6 {
        let middle = phylo_point(layout, radius, start + PI);
        return format!(
            "M {} {} A {} {} 0 0 1 {} {} A {} {} 0 0 1 {} {}",
            num(first.0),
            num(first.1),
            num(radius),
            num(radius),
            num(middle.0),
            num(middle.1),
            num(radius),
            num(radius),
            num(first.0),
            num(first.1)
        );
    }
    let last = phylo_point(layout, radius, end);
    format!(
        "M {} {} A {} {} 0 {} 1 {} {}",
        num(first.0),
        num(first.1),
        num(radius),
        num(radius),
        usize::from(delta > PI),
        num(last.0),
        num(last.1)
    )
}

fn upright_tangent(angle: f64) -> f64 {
    let mut degrees = (angle.to_degrees() + 90.0).rem_euclid(360.0);
    if degrees > 90.0 && degrees < 270.0 {
        degrees += 180.0;
    }
    degrees
}

fn circular_mean(angles: impl IntoIterator<Item = f64>) -> f64 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut count = 0usize;
    for angle in angles {
        x += angle.cos();
        y += angle.sin();
        count += 1;
    }
    if count == 0 || (x.abs() < 1e-12 && y.abs() < 1e-12) {
        0.0
    } else {
        y.atan2(x)
    }
}

#[derive(Debug, Clone, Copy)]
struct MapRect {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
}

#[derive(Debug, Clone, Copy)]
struct MapLayout {
    map: MapRect,
    legend_x: f64,
    legend_y: f64,
    legend_width: f64,
    warning_y: f64,
}

#[derive(Debug)]
struct LandRing {
    hole: bool,
    points: Vec<GeoPosition>,
}

#[derive(Debug, Default)]
struct FlowState {
    unresolved: usize,
    hidden: usize,
}

fn finite_at_least(value: f64, minimum: f64, fallback: f64) -> f64 {
    if value.is_finite() {
        value.max(minimum)
    } else {
        fallback
    }
}

fn data_number(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "not finite".to_string()
    }
}

fn location_names(locations: &[GeoLocation]) -> BTreeMap<&str, Vec<usize>> {
    let mut names: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (index, location) in locations.iter().enumerate() {
        names.entry(&location.name).or_default().push(index);
    }
    names
}

fn unique_endpoint(names: &BTreeMap<&str, Vec<usize>>, name: &str) -> bool {
    names.get(name).is_some_and(|indices| indices.len() == 1)
}

fn category_indices(locations: &[GeoLocation], flows: &[GeoFlow]) -> BTreeMap<String, usize> {
    let mut categories = BTreeMap::new();
    for category in locations
        .iter()
        .filter_map(|location| location.category.as_ref())
        .chain(flows.iter().filter_map(|flow| flow.category.as_ref()))
    {
        let next = categories.len();
        categories.entry(category.clone()).or_insert(next);
    }
    categories
}

fn category_style<'a>(
    category: Option<&str>,
    categories: &BTreeMap<String, usize>,
    theme: &'a Theme,
) -> (&'a str, Symbol) {
    let index = category
        .and_then(|category| categories.get(category))
        .copied()
        .unwrap_or(0);
    (theme.color(index), theme.symbol(index))
}

fn phylo_category_style<'a>(
    category: Option<&str>,
    categories: &BTreeMap<String, usize>,
    theme: &'a Theme,
) -> (&'a str, Symbol) {
    match category.and_then(|category| categories.get(category)) {
        Some(index) => (theme.color(*index), theme.symbol(*index)),
        None => (&theme.muted, Symbol::Circle),
    }
}

fn world_land() -> &'static [LandRing] {
    static LAND: OnceLock<Vec<LandRing>> = OnceLock::new();
    LAND.get_or_init(|| {
        WORLD_110M
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .filter_map(|line| {
                let (kind, points) = line.split_once(' ')?;
                let points = points
                    .split_whitespace()
                    .filter_map(|point| {
                        let (longitude, latitude) = point.split_once(',')?;
                        Some(GeoPosition::new(
                            latitude.parse().ok()?,
                            longitude.parse().ok()?,
                        ))
                    })
                    .collect::<Vec<_>>();
                (points.len() >= 3).then_some(LandRing {
                    hole: kind == "H",
                    points,
                })
            })
            .collect()
    })
}

fn draw_base_map(
    svg: &mut SvgWriter,
    area: MapRect,
    projection: GeoProjection,
    theme: &Theme,
    show_graticule: bool,
) {
    let water = mix(theme.surface(), &theme.accent, 0.07);
    let land = mix(theme.surface(), &theme.muted, 0.30);
    let coast = mix(&theme.rule, &theme.muted, 0.30);
    match projection {
        GeoProjection::Orthographic { .. } => {
            let circle = circle_path(area);
            svg.path(&circle, &water, 1.0);
            svg.begin_clip_path(&circle);
            for ring in world_land() {
                let fill = if ring.hole { &water } else { &land };
                for path in orthographic_ring_paths(ring, projection, area) {
                    svg.path(&path, fill, 1.0);
                    svg.path_stroked(&path, &coast, theme.tokens.hairline * 0.65);
                }
            }
            if show_graticule {
                draw_graticule(svg, area, projection, theme);
            }
            svg.end_group();
            svg.path_stroked(&circle, &theme.rule, theme.tokens.hairline);
        }
        GeoProjection::Equirectangular | GeoProjection::Mercator => {
            svg.rect(area.x, area.y, area.w, area.h, &water);
            svg.begin_clip(area.x, area.y, area.w, area.h);
            for ring in world_land() {
                let fill = if ring.hole { &water } else { &land };
                for path in rectangular_ring_paths(ring, projection, area) {
                    svg.path(&path, fill, 1.0);
                    svg.path_stroked(&path, &coast, theme.tokens.hairline * 0.65);
                }
            }
            if show_graticule {
                draw_graticule(svg, area, projection, theme);
            }
            svg.end_group();
            svg.rect_outline(
                area.x,
                area.y,
                area.w,
                area.h,
                &theme.rule,
                theme.tokens.hairline,
            );
        }
    }
}

fn rectangular_ring_paths(
    ring: &LandRing,
    projection: GeoProjection,
    area: MapRect,
) -> Vec<String> {
    let mut unwrapped = Vec::with_capacity(ring.points.len());
    let mut previous = ring.points[0].longitude;
    unwrapped.push(previous);
    for point in ring.points.iter().skip(1) {
        let mut longitude = point.longitude;
        while longitude - previous > 180.0 {
            longitude -= 360.0;
        }
        while longitude - previous < -180.0 {
            longitude += 360.0;
        }
        unwrapped.push(longitude);
        previous = longitude;
    }
    [-area.w, 0.0, area.w]
        .into_iter()
        .filter_map(|shift| {
            let points = ring
                .points
                .iter()
                .zip(&unwrapped)
                .filter_map(|(point, longitude)| {
                    project_rectangular(
                        GeoPosition::new(point.latitude, *longitude),
                        projection,
                        area,
                        false,
                    )
                    .map(|(x, y)| (x + shift, y))
                })
                .collect::<Vec<_>>();
            path_from_points(&points, true)
        })
        .collect()
}

fn orthographic_ring_paths(
    ring: &LandRing,
    projection: GeoProjection,
    area: MapRect,
) -> Vec<String> {
    let visible = ring
        .points
        .iter()
        .filter(|point| orthographic_unit(**point, projection).is_some())
        .count();
    if visible == ring.points.len() {
        let points = ring
            .points
            .iter()
            .filter_map(|point| project(*point, projection, area))
            .collect::<Vec<_>>();
        return path_from_points(&points, true).into_iter().collect();
    }
    if visible == 0 {
        return Vec::new();
    }
    let start = ring
        .points
        .iter()
        .position(|point| orthographic_unit(*point, projection).is_none())
        .unwrap_or(0);
    let mut runs: Vec<Vec<GeoPosition>> = Vec::new();
    let mut run = Vec::new();
    for offset in 0..ring.points.len() {
        let a = ring.points[(start + offset) % ring.points.len()];
        let b = ring.points[(start + offset + 1) % ring.points.len()];
        let a_visible = orthographic_visibility(a, projection) >= 0.0;
        let b_visible = orthographic_visibility(b, projection) >= 0.0;
        match (a_visible, b_visible) {
            (false, true) => {
                if let Some(intersection) = horizon_intersection(a, b, projection) {
                    run.push(intersection);
                }
                run.push(b);
            }
            (true, true) => {
                if run.is_empty() {
                    run.push(a);
                }
                run.push(b);
            }
            (true, false) => {
                if run.is_empty() {
                    run.push(a);
                }
                if let Some(intersection) = horizon_intersection(a, b, projection) {
                    run.push(intersection);
                }
                if run.len() >= 3 {
                    runs.push(std::mem::take(&mut run));
                } else {
                    run.clear();
                }
            }
            (false, false) => {}
        }
    }
    if run.len() >= 3 {
        runs.push(run);
    }
    runs.into_iter()
        .filter_map(|run| {
            let mut points = run
                .iter()
                .filter_map(|point| project(*point, projection, area))
                .collect::<Vec<_>>();
            close_along_horizon(&mut points, area);
            path_from_points(&points, true)
        })
        .collect()
}

fn close_along_horizon(points: &mut Vec<(f64, f64)>, area: MapRect) {
    let (Some(first), Some(last)) = (points.first().copied(), points.last().copied()) else {
        return;
    };
    let cx = area.x + area.w / 2.0;
    let cy = area.y + area.h / 2.0;
    let radius = area.w.min(area.h) / 2.0;
    let start = (last.1 - cy).atan2(last.0 - cx);
    let end = (first.1 - cy).atan2(first.0 - cx);
    let mut delta = (end - start).rem_euclid(TAU);
    if delta > PI {
        delta -= TAU;
    }
    let steps = ((delta.abs() * 12.0).ceil() as usize).clamp(2, 28);
    for index in 1..steps {
        let angle = start + delta * index as f64 / steps as f64;
        points.push((cx + radius * angle.cos(), cy + radius * angle.sin()));
    }
}

fn draw_graticule(svg: &mut SvgWriter, area: MapRect, projection: GeoProjection, theme: &Theme) {
    let color = mix(&theme.rule, theme.surface(), 0.18);
    for latitude in [-60.0, -30.0, 0.0, 30.0, 60.0] {
        let points = (-180..=180)
            .step_by(3)
            .map(|longitude| GeoPosition::new(latitude, longitude as f64));
        draw_projected_line(
            svg,
            points,
            area,
            projection,
            &color,
            theme.tokens.hairline * 0.55,
        );
    }
    for longitude in (-150..=180).step_by(30) {
        let points = (-90..=90)
            .step_by(2)
            .map(|latitude| GeoPosition::new(latitude as f64, longitude as f64));
        draw_projected_line(
            svg,
            points,
            area,
            projection,
            &color,
            theme.tokens.hairline * 0.55,
        );
    }
}

fn draw_projected_line(
    svg: &mut SvgWriter,
    points: impl IntoIterator<Item = GeoPosition>,
    area: MapRect,
    projection: GeoProjection,
    color: &str,
    width: f64,
) {
    let mut run = Vec::new();
    for point in points {
        if let Some(projected) = project(point, projection, area) {
            run.push(projected);
        } else {
            if run.len() >= 2 {
                svg.polyline(&run, color, width);
            }
            run.clear();
        }
    }
    if run.len() >= 2 {
        svg.polyline(&run, color, width);
    }
}

fn draw_locations(
    svg: &mut SvgWriter,
    area: MapRect,
    projection: GeoProjection,
    theme: &Theme,
    categories: &BTreeMap<String, usize>,
    locations: &[GeoLocation],
    show_labels: bool,
) -> usize {
    let mut hidden = 0;
    for location in locations {
        if !location.position.is_valid() {
            continue;
        }
        let Some((x, y)) = project(location.position, projection, area) else {
            hidden += 1;
            continue;
        };
        let (color, symbol) = category_style(location.category.as_deref(), categories, theme);
        let radius =
            theme.tokens.marker_radius * 1.15 * (location.count as f64).sqrt().clamp(1.0, 3.2);
        svg.begin_titled(&location.title());
        svg.symbol_ringed(
            x,
            y,
            radius,
            symbol,
            color,
            theme.surface(),
            theme.tokens.hairline * 1.5,
        );
        svg.end_group();
        if show_labels {
            svg.text(
                x + radius + theme.tokens.label_gap,
                y + theme.font_size * 0.32,
                &location.name,
                &theme.foreground,
                theme.font_size - 1.0,
                Anchor::Start,
            );
        }
    }
    hidden
}

fn draw_flows(
    svg: &mut SvgWriter,
    area: MapRect,
    projection: GeoProjection,
    theme: &Theme,
    categories: &BTreeMap<String, usize>,
    locations: &[GeoLocation],
    flows: &[GeoFlow],
) -> FlowState {
    let names = location_names(locations);
    let mut state = FlowState::default();
    for flow in flows {
        let (Some(from), Some(to)) = (names.get(flow.from.as_str()), names.get(flow.to.as_str()))
        else {
            state.unresolved += 1;
            continue;
        };
        if from.len() != 1 || to.len() != 1 {
            state.unresolved += 1;
            continue;
        }
        let from = &locations[from[0]];
        let to = &locations[to[0]];
        if !from.position.is_valid() || !to.position.is_valid() {
            state.unresolved += 1;
            continue;
        }
        let (Some(start), Some(end)) = (
            project(from.position, projection, area),
            project(to.position, projection, area),
        ) else {
            state.hidden += 1;
            continue;
        };
        let (base, _) = category_style(flow.category.as_deref(), categories, theme);
        let color = mix(base, theme.surface(), 0.26);
        let dx = end.0 - start.0;
        let dy = end.1 - start.1;
        let distance = dx.hypot(dy).max(1.0);
        let bend = distance.min(area.w.min(area.h) * 0.45) * 0.17;
        let control = (
            (start.0 + end.0) / 2.0 - dy / distance * bend,
            (start.1 + end.1) / 2.0 + dx / distance * bend,
        );
        let path = format!(
            "M {} {} Q {} {} {} {}",
            num(start.0),
            num(start.1),
            num(control.0),
            num(control.1),
            num(end.0),
            num(end.1)
        );
        let width = theme.tokens.stroke * (0.72 + flow.weight.sqrt().clamp(0.0, 4.0) * 0.34);
        svg.begin_titled(&flow.title());
        svg.path_stroked(&path, &color, width);
        if flow.directed && distance > theme.tokens.arrow_size * 2.0 {
            draw_flow_arrow(svg, start, control, end, base, theme.tokens.arrow_size);
        }
        svg.end_group();
    }
    state
}

fn draw_flow_arrow(
    svg: &mut SvgWriter,
    start: (f64, f64),
    control: (f64, f64),
    end: (f64, f64),
    color: &str,
    size: f64,
) {
    let t = 0.82;
    let point = quadratic(start, control, end, t);
    let tangent = (
        2.0 * (1.0 - t) * (control.0 - start.0) + 2.0 * t * (end.0 - control.0),
        2.0 * (1.0 - t) * (control.1 - start.1) + 2.0 * t * (end.1 - control.1),
    );
    let length = tangent.0.hypot(tangent.1).max(1.0);
    let ux = tangent.0 / length;
    let uy = tangent.1 / length;
    let px = -uy;
    let py = ux;
    svg.polygon(
        &[
            (point.0 + ux * size, point.1 + uy * size),
            (
                point.0 - ux * size * 0.65 + px * size * 0.55,
                point.1 - uy * size * 0.65 + py * size * 0.55,
            ),
            (
                point.0 - ux * size * 0.65 - px * size * 0.55,
                point.1 - uy * size * 0.65 - py * size * 0.55,
            ),
        ],
        color,
    );
}

fn quadratic(start: (f64, f64), control: (f64, f64), end: (f64, f64), t: f64) -> (f64, f64) {
    let one = 1.0 - t;
    (
        one * one * start.0 + 2.0 * one * t * control.0 + t * t * end.0,
        one * one * start.1 + 2.0 * one * t * control.1 + t * t * end.1,
    )
}

fn draw_legend(
    svg: &mut SvgWriter,
    layout: MapLayout,
    theme: &Theme,
    categories: &BTreeMap<String, usize>,
) {
    if categories.is_empty() || layout.legend_width <= 0.0 {
        return;
    }
    svg.text_bold(
        layout.legend_x,
        layout.legend_y + theme.font_size,
        "Categories",
        &theme.foreground,
        theme.font_size,
        Anchor::Start,
    );
    let mut ordered = categories.iter().collect::<Vec<_>>();
    ordered.sort_by_key(|(_, index)| **index);
    for (row, (category, index)) in ordered.into_iter().enumerate() {
        let y = layout.legend_y + theme.font_size + 12.0 + row as f64 * (theme.font_size + 7.0);
        let color = theme.color(*index);
        svg.symbol_ringed(
            layout.legend_x + theme.tokens.marker_radius,
            y - theme.font_size * 0.28,
            theme.tokens.marker_radius,
            theme.symbol(*index),
            color,
            theme.surface(),
            theme.tokens.hairline,
        );
        let label = fit_text(
            category,
            (layout.legend_width - theme.tokens.marker_radius * 3.0).max(1.0),
            theme.font_size - 1.0,
        );
        if label != *category {
            svg.begin_titled(category);
        }
        svg.text(
            layout.legend_x + theme.tokens.marker_radius * 2.5,
            y,
            &label,
            &theme.muted,
            theme.font_size - 1.0,
            Anchor::Start,
        );
        if label != *category {
            svg.end_group();
        }
    }
}

fn draw_warning(
    svg: &mut SvgWriter,
    y: f64,
    theme: &Theme,
    invalid_locations: usize,
    hidden_locations: usize,
    unresolved_flows: usize,
    hidden_flows: usize,
) {
    let mut parts = Vec::new();
    if invalid_locations > 0 {
        parts.push(plural(
            invalid_locations,
            "invalid location",
            "invalid locations",
        ));
    }
    if hidden_locations > 0 {
        parts.push(plural(
            hidden_locations,
            "location outside projection",
            "locations outside projection",
        ));
    }
    if unresolved_flows > 0 {
        parts.push(plural(
            unresolved_flows,
            "unresolved flow",
            "unresolved flows",
        ));
    }
    if hidden_flows > 0 {
        parts.push(plural(
            hidden_flows,
            "flow outside projection",
            "flows outside projection",
        ));
    }
    if parts.is_empty() {
        return;
    }
    svg.text(
        18.0,
        y,
        &parts.join("; "),
        &theme.muted,
        theme.font_size - 1.0,
        Anchor::Start,
    );
}

fn plural(count: usize, one: &str, many: &str) -> String {
    format!("{count} {}", if count == 1 { one } else { many })
}

fn project(position: GeoPosition, projection: GeoProjection, area: MapRect) -> Option<(f64, f64)> {
    if !position.is_valid() {
        return None;
    }
    match projection {
        GeoProjection::Equirectangular | GeoProjection::Mercator => {
            project_rectangular(position, projection, area, true)
        }
        GeoProjection::Orthographic { .. } => {
            let (x, y) = orthographic_unit(position, projection)?;
            let radius = area.w.min(area.h) / 2.0;
            Some((
                area.x + area.w / 2.0 + x * radius,
                area.y + area.h / 2.0 - y * radius,
            ))
        }
    }
}

fn project_rectangular(
    position: GeoPosition,
    projection: GeoProjection,
    area: MapRect,
    require_bounds: bool,
) -> Option<(f64, f64)> {
    if !position.latitude.is_finite() || !position.longitude.is_finite() {
        return None;
    }
    if require_bounds && !position.is_valid() {
        return None;
    }
    let x = area.x + (position.longitude + 180.0) / 360.0 * area.w;
    let y_fraction = match projection {
        GeoProjection::Equirectangular => (90.0 - position.latitude) / 180.0,
        GeoProjection::Mercator => {
            let latitude = position.latitude.clamp(-85.051_128_78, 85.051_128_78);
            let radians = latitude.to_radians();
            (1.0 - (radians.tan() + 1.0 / radians.cos()).ln() / PI) / 2.0
        }
        GeoProjection::Orthographic { .. } => return None,
    };
    Some((x, area.y + y_fraction * area.h))
}

fn orthographic_visibility(position: GeoPosition, projection: GeoProjection) -> f64 {
    let GeoProjection::Orthographic {
        center_latitude,
        center_longitude,
    } = projection
    else {
        return 1.0;
    };
    let (latitude, longitude) = position.radians();
    let centre_latitude = center_latitude.to_radians();
    let delta = longitude - center_longitude.to_radians();
    centre_latitude.sin() * latitude.sin() + centre_latitude.cos() * latitude.cos() * delta.cos()
}

fn orthographic_unit(position: GeoPosition, projection: GeoProjection) -> Option<(f64, f64)> {
    let GeoProjection::Orthographic {
        center_latitude,
        center_longitude,
    } = projection
    else {
        return None;
    };
    if orthographic_visibility(position, projection) < -1e-9 {
        return None;
    }
    let (latitude, longitude) = position.radians();
    let centre_latitude = center_latitude.to_radians();
    let delta = longitude - center_longitude.to_radians();
    Some((
        latitude.cos() * delta.sin(),
        centre_latitude.cos() * latitude.sin()
            - centre_latitude.sin() * latitude.cos() * delta.cos(),
    ))
}

fn horizon_intersection(
    a: GeoPosition,
    b: GeoPosition,
    projection: GeoProjection,
) -> Option<GeoPosition> {
    let da = orthographic_visibility(a, projection);
    let db = orthographic_visibility(b, projection);
    let denominator = da - db;
    if denominator.abs() < 1e-12 {
        return None;
    }
    let t = (da / denominator).clamp(0.0, 1.0);
    let a = sphere(a);
    let b = sphere(b);
    let mut point = (
        a.0 + (b.0 - a.0) * t,
        a.1 + (b.1 - a.1) * t,
        a.2 + (b.2 - a.2) * t,
    );
    let length = (point.0 * point.0 + point.1 * point.1 + point.2 * point.2).sqrt();
    if length <= 1e-12 {
        return None;
    }
    point.0 /= length;
    point.1 /= length;
    point.2 /= length;
    Some(GeoPosition::new(
        point.2.asin().to_degrees(),
        point.1.atan2(point.0).to_degrees(),
    ))
}

fn sphere(position: GeoPosition) -> (f64, f64, f64) {
    let (latitude, longitude) = position.radians();
    (
        latitude.cos() * longitude.cos(),
        latitude.cos() * longitude.sin(),
        latitude.sin(),
    )
}

fn path_from_points(points: &[(f64, f64)], close: bool) -> Option<String> {
    let first = points.first()?;
    if points.len() < if close { 3 } else { 2 } {
        return None;
    }
    let mut path = format!("M {} {}", num(first.0), num(first.1));
    for point in points.iter().skip(1) {
        path.push_str(&format!(" L {} {}", num(point.0), num(point.1)));
    }
    if close {
        path.push_str(" Z");
    }
    Some(path)
}

fn circle_path(area: MapRect) -> String {
    let cx = area.x + area.w / 2.0;
    let cy = area.y + area.h / 2.0;
    let radius = area.w.min(area.h) / 2.0;
    format!(
        "M {} {} A {} {} 0 1 1 {} {} A {} {} 0 1 1 {} {} Z",
        num(cx - radius),
        num(cy),
        num(radius),
        num(radius),
        num(cx + radius),
        num(cy),
        num(radius),
        num(radius),
        num(cx - radius),
        num(cy)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn phylo_tree() -> Tree {
        Tree::parse_annotated_newick(concat!(
            "((A[&date=2023,country=Peru]:1,B[&date=2024,country=Peru]:1)",
            "[&country=Peru]:1,(C[&date=2023.5,country=Spain]:1,",
            "D[&date=2024.2,country=Kenya]:1):1)[&date=2022];"
        ))
        .unwrap()
    }

    fn phylo_coordinates() -> [GeoLocation; 3] {
        [
            GeoLocation::new("Peru", -9.19, -75.0152),
            GeoLocation::new("Spain", 40.4637, -3.7492),
            GeoLocation::new("Kenya", -0.0236, 37.9062),
        ]
    }

    #[test]
    fn equirectangular_places_the_origin_at_the_centre() {
        let area = MapRect {
            x: 10.0,
            y: 20.0,
            w: 360.0,
            h: 180.0,
        };
        assert_eq!(
            project(
                GeoPosition::new(0.0, 0.0),
                GeoProjection::Equirectangular,
                area
            ),
            Some((190.0, 110.0))
        );
    }

    #[test]
    fn orthographic_projection_hides_the_back_hemisphere() {
        let area = MapRect {
            x: 0.0,
            y: 0.0,
            w: 200.0,
            h: 200.0,
        };
        let projection = GeoProjection::orthographic(0.0, 0.0);
        assert_eq!(
            project(GeoPosition::new(0.0, 0.0), projection, area),
            Some((100.0, 100.0))
        );
        assert!(project(GeoPosition::new(0.0, 179.0), projection, area).is_none());
    }

    #[test]
    fn invalid_coordinates_are_counted_and_named_on_the_page() {
        let map = Map::new()
            .push(GeoLocation::new("valid", 0.0, 0.0))
            .push(GeoLocation::new("invalid", 120.0, 0.0));
        assert_eq!(map.invalid_location_count(), 1);
        let svg = map.to_svg();
        assert!(svg.contains(">1 invalid location</text>"), "{svg}");
        assert!(svg.contains("1 invalid locations and 0 unresolved flows"));
    }

    #[test]
    fn location_tooltips_keep_exact_coordinates_and_values() {
        let svg = Map::new()
            .push(
                GeoLocation::new("Lima", -12.046, -77.043)
                    .category("Peru")
                    .value(42.5)
                    .count(7),
            )
            .to_svg();
        assert!(svg.contains(
            "<title>Lima; latitude -12.046; longitude -77.043; 7 observations; category Peru; value 42.5</title>"
        ));
    }

    #[test]
    fn categories_use_shape_as_well_as_colour() {
        let svg = Map::new()
            .push(GeoLocation::new("A", 0.0, 0.0).category("one"))
            .push(GeoLocation::new("B", 10.0, 10.0).category("two"))
            .to_svg();
        assert!(svg.contains("<circle"), "{svg}");
        assert!(svg.contains("<rect"), "{svg}");
        assert!(svg.contains(">one</text>"), "{svg}");
        assert!(svg.contains(">two</text>"), "{svg}");
    }

    #[test]
    fn flows_are_named_and_direction_is_visible() {
        let svg = Map::new()
            .push(GeoLocation::new("Lima", -12.046, -77.043))
            .push(GeoLocation::new("Madrid", 40.417, -3.704))
            .push_flow(GeoFlow::new("Lima", "Madrid").weight(4.0))
            .to_svg();
        assert!(svg.contains("<title>flow, source Lima, target Madrid, weight 4</title>"));
        assert!(svg.contains(" Q "), "{svg}");
        assert!(svg.contains("<polygon"), "arrowhead: {svg}");
    }

    #[test]
    fn a_missing_or_ambiguous_flow_endpoint_is_not_guessed() {
        let missing = Map::new()
            .push(GeoLocation::new("A", 0.0, 0.0))
            .push_flow(GeoFlow::new("A", "B"));
        assert_eq!(missing.unresolved_flow_count(), 1);
        let ambiguous = Map::new()
            .push(GeoLocation::new("A", 0.0, 0.0))
            .push(GeoLocation::new("A", 1.0, 1.0))
            .push(GeoLocation::new("B", 2.0, 2.0))
            .push_flow(GeoFlow::new("A", "B"));
        assert_eq!(ambiguous.unresolved_flow_count(), 1);
    }

    #[test]
    fn base_land_is_compiled_into_every_projection_without_non_finite_output() {
        for projection in [
            GeoProjection::Equirectangular,
            GeoProjection::Mercator,
            GeoProjection::orthographic(20.0, 0.0),
        ] {
            let svg = Map::new().projection(projection).to_svg();
            assert!(svg.matches("<path").count() > 100, "{projection:?}");
            assert!(!svg.contains("NaN"), "{projection:?}");
            assert!(!svg.contains("inf"), "{projection:?}");
        }
    }

    #[test]
    fn rendering_is_deterministic_and_accessible() {
        let map = Map::new()
            .title("Sites")
            .description("Three surveillance sites.")
            .push(GeoLocation::new("A", 0.0, 0.0));
        assert_eq!(map.to_svg(), map.to_svg());
        let svg = map.to_svg();
        assert!(svg.contains("role=\"img\""));
        assert!(svg.contains("<title id=\"karyon-title\">Sites</title>"));
        assert!(svg.contains("<desc id=\"karyon-desc\">Three surveillance sites.</desc>"));
    }

    #[test]
    fn phylo_map_preserves_the_tree_and_aggregates_exact_location_counts() {
        let tree = phylo_tree();
        let leaves = tree.leaf_names();
        let map = PhyloMap::new(tree)
            .location_by("country")
            .coordinates(phylo_coordinates());
        assert_eq!(map.tree().leaf_names(), leaves);
        let svg = map.to_svg();
        assert!(svg.contains("<title>Peru; 2 mapped tips</title>"), "{svg}");
        assert!(svg.contains("<title>Spain; 1 mapped tip</title>"), "{svg}");
        assert!(svg.contains("<title>Kenya; 1 mapped tip</title>"), "{svg}");
        assert!(
            svg.contains("<title>Peru; latitude -9.19; longitude -75.0152; 2 mapped tips</title>")
        );
    }

    #[test]
    fn individual_phylo_connectors_name_every_sample() {
        let svg = PhyloMap::new(phylo_tree())
            .location_by("country")
            .coordinates(phylo_coordinates())
            .connector(PhyloConnector::Individual)
            .to_svg();
        for title in [
            "A; location Peru",
            "B; location Peru",
            "C; location Spain",
            "D; location Kenya",
        ] {
            assert!(svg.contains(&format!("<title>{title}</title>")), "{title}");
        }
    }

    #[test]
    fn unlocated_phylo_tips_are_counted_instead_of_dropped() {
        let map = PhyloMap::new(phylo_tree())
            .location_by("country")
            .coordinates(phylo_coordinates().into_iter().take(2));
        assert_eq!(map.unlocated_tip_count(), 1);
        let svg = map.to_svg();
        assert!(svg.contains(">1 unlocated tip</text>"), "{svg}");
        assert!(svg.contains("<title>D; country Kenya</title>"), "{svg}");
    }

    #[test]
    fn ambiguous_coordinate_names_are_not_guessed() {
        let map = PhyloMap::new(phylo_tree())
            .location_by("country")
            .coordinates([
                GeoLocation::new("Peru", -9.19, -75.0152),
                GeoLocation::new("Peru", -12.046, -77.043),
                GeoLocation::new("Spain", 40.4637, -3.7492),
                GeoLocation::new("Kenya", -0.0236, 37.9062),
            ]);
        assert_eq!(map.unlocated_tip_count(), 2);
        assert!(map.to_svg().contains(">2 unlocated tips</text>"));
    }

    #[test]
    fn hidden_phylo_locations_are_reported() {
        let tree = Tree::parse_annotated_newick("(A[&place=Far]:1);").unwrap();
        let svg = PhyloMap::new(tree)
            .location_by("place")
            .projection(GeoProjection::orthographic(0.0, 0.0))
            .coordinate(GeoLocation::new("Far", 0.0, 179.0))
            .to_svg();
        assert!(
            svg.contains(">1 mapped location outside projection</text>"),
            "{svg}"
        );
    }

    #[test]
    fn phylo_time_guides_keep_exact_values_and_units() {
        let tree = Tree::parse_annotated_newick(
            "(A[&date=2023,place=Near]:1,B[&date=2025,place=Near]:1)[&date=2022];",
        )
        .unwrap();
        let svg = PhyloMap::new(tree)
            .location_by("place")
            .coordinate(GeoLocation::new("Near", 0.0, 0.0))
            .time("date")
            .time_unit("year")
            .to_svg();
        for value in ["2022 year", "2023.5 year", "2025 year"] {
            assert!(svg.contains(&format!(">{value}</text>")), "{value}: {svg}");
        }
    }

    #[test]
    fn an_incomplete_phylo_time_layout_is_explicit() {
        let tree =
            Tree::parse_annotated_newick("(A[&date=2023,place=Near]:1,B[&place=Near]:1);").unwrap();
        let svg = PhyloMap::new(tree)
            .location_by("place")
            .coordinate(GeoLocation::new("Near", 0.0, 0.0))
            .time("date")
            .to_svg();
        assert!(
            svg.contains(">requested time layout unavailable</text>"),
            "{svg}"
        );
    }

    #[test]
    fn phylo_map_rendering_is_deterministic_accessible_and_finite() {
        let map = PhyloMap::new(phylo_tree())
            .title("Circular surveillance")
            .description("Four synthetic samples at three supplied locations.")
            .location_by("country")
            .coordinates(phylo_coordinates());
        assert_eq!(map.to_svg(), map.to_svg());
        let svg = map.to_svg();
        assert!(svg.contains("role=\"img\""));
        assert!(svg.contains("<title id=\"karyon-title\">Circular surveillance</title>"));
        assert!(svg.contains(
            "<desc id=\"karyon-desc\">Four synthetic samples at three supplied locations.</desc>"
        ));
        assert!(!svg.contains("NaN"));
        assert!(!svg.contains("inf"));
    }

    #[test]
    fn profiles_scale_map_themes_once() {
        let map = Map::new().profile(RenderProfile::Presentation);
        let phylo = PhyloMap::new(phylo_tree()).profile(RenderProfile::Presentation);
        assert_eq!(map.theme.title_font_size, Theme::light().title_font_size);
        assert_eq!(phylo.theme.title_font_size, Theme::light().title_font_size);
        assert_eq!(map.visual_scale, 1.35);
        assert_eq!(phylo.visual_scale, 1.35);
    }
}
