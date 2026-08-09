//! A phylogeny drawn around a globe, with its tips joined to where they came
//! from.
//!
//! [`PhyloMap`] is not a [`Map`](super::Map) with a tree on it. The tree is
//! the frame and the map sits inside it, so the composition owns its own
//! layout, its own legend and its own accounting of what it could not place.
//! What it shares with a plain map is the projection and the marks, which is
//! why both live one module away.

use super::*;

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
    pub(super) theme: Theme,
    pub(super) visual_scale: f64,
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
pub(super) struct PhyloTime {
    pub(super) key: String,
    pub(super) direction: TimeDirection,
    pub(super) unit: Option<String>,
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
        self.diameter = finite_within(diameter, 320.0, f64::INFINITY, 760.0);
        self
    }

    /// Sets whitespace around the drawing.
    pub fn margin(mut self, margin: f64) -> Self {
        self.margin = finite_within(margin, 0.0, f64::INFINITY, 18.0);
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
        self.visual_scale = finite_within(scale, 0.25, f64::INFINITY, 1.0);
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

    pub(super) fn document_name(&self) -> String {
        match (&self.title, &self.subtitle) {
            (Some(title), Some(subtitle)) => format!("{title}, {subtitle}"),
            (Some(title), None) => title.clone(),
            (None, Some(subtitle)) => subtitle.clone(),
            (None, None) => "Circular phylogeographic map".to_string(),
        }
    }

    pub(super) fn document_description(
        &self,
        scene: &PhyloScene,
        unlocated: usize,
        hidden: usize,
    ) -> String {
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

    pub(super) fn layout(&self, theme: &Theme, has_categories: bool) -> PhyloLayout {
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
pub(super) struct PhyloLayout {
    pub(super) width: f64,
    pub(super) height: f64,
    pub(super) cx: f64,
    pub(super) cy: f64,
    pub(super) map: MapRect,
    pub(super) tree_inner: f64,
    pub(super) tree_outer: f64,
    pub(super) start: f64,
    pub(super) sweep: f64,
    pub(super) legend_x: f64,
    pub(super) legend_y: f64,
    pub(super) legend_width: f64,
    pub(super) warning_y: f64,
}

pub(super) struct PhyloScene {
    pub(super) placements: Vec<Placement>,
    pub(super) leaves: Vec<usize>,
    pub(super) minimum: f64,
    pub(super) maximum: f64,
    pub(super) temporal: bool,
    pub(super) time_failed: bool,
    pub(super) direction: TimeDirection,
}

impl PhyloScene {
    pub(super) fn new(tree: &Tree, shape: TreeShape, time: Option<&PhyloTime>) -> Self {
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

    pub(super) fn fraction(&self, value: f64) -> f64 {
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

    pub(super) fn radius(&self, layout: PhyloLayout, value: f64) -> f64 {
        layout.tree_outer - self.fraction(value) * (layout.tree_outer - layout.tree_inner)
    }

    pub(super) fn angle(&self, layout: PhyloLayout, row: f64) -> f64 {
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

pub(super) struct MappedTipGroup<'a> {
    pub(super) location: &'a GeoLocation,
    pub(super) tips: Vec<usize>,
    pub(super) projected: Option<(f64, f64)>,
}

pub(super) fn resolve_tip_location<'a>(
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

pub(super) fn mapped_tip_groups<'a>(
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

pub(super) fn phylo_categories(
    tree: &Tree,
    scene: &PhyloScene,
    key: &str,
) -> BTreeMap<String, usize> {
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

pub(super) fn inherited_tree_annotation<'a>(
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

pub(super) fn phylo_point(layout: PhyloLayout, radius: f64, angle: f64) -> (f64, f64) {
    (
        layout.cx + radius * angle.cos(),
        layout.cy + radius * angle.sin(),
    )
}

pub(super) fn draw_phylo_title(map: &PhyloMap, svg: &mut SvgWriter, theme: &Theme) {
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

pub(super) fn draw_phylo_time_guides(
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

pub(super) fn draw_phylo_connectors(
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

pub(super) fn connector_path(start: (f64, f64), end: (f64, f64), layout: PhyloLayout) -> String {
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

pub(super) fn draw_phylo_tree(
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

pub(super) fn phylo_branch_title(map: &PhyloMap, node: usize) -> Option<String> {
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

pub(super) fn draw_phylo_tip_label(
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

pub(super) fn draw_phylo_locations(
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

pub(super) fn draw_phylo_legend(
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

pub(super) fn draw_phylo_warning(
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

pub(super) fn radial_arc(layout: PhyloLayout, radius: f64, start: f64, end: f64) -> String {
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

pub(super) fn upright_tangent(angle: f64) -> f64 {
    let mut degrees = (angle.to_degrees() + 90.0).rem_euclid(360.0);
    if degrees > 90.0 && degrees < 270.0 {
        degrees += 180.0;
    }
    degrees
}

pub(super) fn circular_mean(angles: impl IntoIterator<Item = f64>) -> f64 {
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
