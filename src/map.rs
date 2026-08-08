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
            format!("latitude {}", num(self.position.latitude)),
            format!("longitude {}", num(self.position.longitude)),
        ];
        if self.count > 1 {
            parts.push(format!("{} observations", self.count));
        }
        if let Some(category) = &self.category {
            parts.push(format!("category {category}"));
        }
        if let Some(value) = self.value {
            parts.push(format!("value {}", num(value)));
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
            num(self.weight)
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
        self.theme = Theme::for_profile(profile);
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
}
