//! The marks a map is made of, shared by both compositions.
//!
//! Land, graticule, locations, flows, legend and the notice that says what
//! could not be drawn. Nothing here decides where a point goes, which is
//! [`project`](super::project)'s job; everything here decides what it looks
//! like once it is there, so a change to a symbol cannot move a sample and a
//! change to the projection cannot restyle one.

use super::*;
use crate::svg::text_exact;

#[derive(Debug, Clone, Copy)]
pub(super) struct MapRect {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) w: f64,
    pub(super) h: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct MapLayout {
    pub(super) map: MapRect,
    pub(super) legend_x: f64,
    pub(super) legend_y: f64,
    pub(super) legend_width: f64,
    pub(super) warning_y: f64,
}

#[derive(Debug)]
pub(super) struct LandRing {
    pub(super) hole: bool,
    pub(super) points: Vec<GeoPosition>,
}

#[derive(Debug, Default)]
pub(super) struct FlowState {
    pub(super) unresolved: usize,
    pub(super) hidden: usize,
}

/// A measured value as map tooltip text.
///
/// The formatting is [`text_exact`](crate::svg::text_exact); the wording is
/// this module's. A map tooltip reads as a sentence, and "value not finite"
/// belongs in one where "value NaN" does not. Sharing the arithmetic and
/// keeping the phrasing local is the whole of what this is for.
pub(super) fn data_number(value: f64) -> String {
    if !value.is_finite() {
        return "not finite".to_string();
    }
    text_exact(value)
}

pub(super) fn location_names(locations: &[GeoLocation]) -> BTreeMap<&str, Vec<usize>> {
    let mut names: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (index, location) in locations.iter().enumerate() {
        names.entry(&location.name).or_default().push(index);
    }
    names
}

pub(super) fn unique_endpoint(names: &BTreeMap<&str, Vec<usize>>, name: &str) -> bool {
    names.get(name).is_some_and(|indices| indices.len() == 1)
}

pub(super) fn category_indices(
    locations: &[GeoLocation],
    flows: &[GeoFlow],
) -> BTreeMap<String, usize> {
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

pub(super) fn category_style<'a>(
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

pub(super) fn phylo_category_style<'a>(
    category: Option<&str>,
    categories: &BTreeMap<String, usize>,
    theme: &'a Theme,
) -> (&'a str, Symbol) {
    match category.and_then(|category| categories.get(category)) {
        Some(index) => (theme.color(*index), theme.symbol(*index)),
        None => (&theme.muted, Symbol::Circle),
    }
}

pub(super) fn world_land() -> &'static [LandRing] {
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

pub(super) fn draw_base_map(
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
                    if !touches(&path, area, theme.tokens.hairline) {
                        continue;
                    }
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
                    if !touches(&path, area, theme.tokens.hairline) {
                        continue;
                    }
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

/// Whether a finished path can put any ink inside `area`.
///
/// Two thirds of the coastline of a world map lies outside the window a figure
/// draws, and every ring of it was written out in full and then thrown away by
/// the clip: measured on the committed map figure, 996 of the 1,536 paths
/// inside the clipped groups had a bounding box that missed the clip entirely.
/// They cost bytes in the document and drew nothing.
///
/// The test is the path's bounding box against the area, widened by the stroke
/// so a coastline that only grazes the edge keeps its hairline. A bounding box
/// is enough because it can only ever be too generous: a path it rejects has no
/// point inside the area, so nothing is dropped that would have drawn.
pub(super) fn touches(path: &str, area: MapRect, margin: f64) -> bool {
    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut number = String::new();
    let mut coords: usize = 0;
    let mut pending: Option<f64> = None;
    let mut flush = |number: &mut String, coords: &mut usize, pending: &mut Option<f64>| {
        if number.is_empty() {
            return;
        }
        if let Ok(value) = number.parse::<f64>() {
            if *coords % 2 == 0 {
                *pending = Some(value);
            } else if let Some(x) = pending.take() {
                min_x = min_x.min(x);
                max_x = max_x.max(x);
                min_y = min_y.min(value);
                max_y = max_y.max(value);
            }
            *coords += 1;
        }
        number.clear();
    };
    for c in path.chars() {
        if c.is_ascii_digit() || c == '.' || (c == '-' && number.is_empty()) {
            number.push(c);
        } else {
            flush(&mut number, &mut coords, &mut pending);
            if c.is_ascii_alphabetic() {
                coords = 0;
                pending = None;
            }
        }
    }
    flush(&mut number, &mut coords, &mut pending);
    if !min_x.is_finite() || !min_y.is_finite() {
        return false;
    }
    max_x >= area.x - margin
        && min_x <= area.x + area.w + margin
        && max_y >= area.y - margin
        && min_y <= area.y + area.h + margin
}

pub(super) fn rectangular_ring_paths(
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

pub(super) fn orthographic_ring_paths(
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

pub(super) fn close_along_horizon(points: &mut Vec<(f64, f64)>, area: MapRect) {
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

pub(super) fn draw_graticule(
    svg: &mut SvgWriter,
    area: MapRect,
    projection: GeoProjection,
    theme: &Theme,
) {
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

pub(super) fn draw_projected_line(
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

pub(super) fn draw_locations(
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

pub(super) fn draw_flows(
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

pub(super) fn draw_flow_arrow(
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

pub(super) fn quadratic(
    start: (f64, f64),
    control: (f64, f64),
    end: (f64, f64),
    t: f64,
) -> (f64, f64) {
    let one = 1.0 - t;
    (
        one * one * start.0 + 2.0 * one * t * control.0 + t * t * end.0,
        one * one * start.1 + 2.0 * one * t * control.1 + t * t * end.1,
    )
}

pub(super) fn draw_legend(
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

pub(super) fn draw_warning(
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

pub(super) fn plural(count: usize, one: &str, many: &str) -> String {
    format!("{count} {}", if count == 1 { one } else { many })
}
