//! The circular projection, and the annotation rings that only it has.
//!
//! A radial tree is the rectangular one in polar coordinates: the same scene,
//! the same rows, with depth becoming radius and row becoming angle. What is
//! genuinely different is here, and the largest part of it is the trait rings,
//! which have no rectangular equivalent because a ring is a column bent round
//! until its ends meet.

use super::*;

#[derive(Debug, Clone, Copy)]
pub(super) struct RadialGeometry {
    pub(super) cx: f64,
    pub(super) cy: f64,
    pub(super) tree_inner: f64,
    pub(super) tree_outer: f64,
    pub(super) ring_outer: f64,
    pub(super) label_radius: f64,
    pub(super) start: f64,
    pub(super) sweep: f64,
    pub(super) terminals: usize,
    pub(super) direction: RadialDirection,
}

impl RadialGeometry {
    pub(super) fn new(track: &TreeTrack, theme: &Theme, scene: &TreeScene, area: Rect) -> Self {
        let label_extent = if track.show_tips {
            scene
                .terminals
                .iter()
                .map(|node| {
                    text_width(
                        &terminal_label(&track.tree, *node, &track.collapsed),
                        theme.font_size - 1.0,
                    )
                })
                .fold(0.0f64, f64::max)
                + 6.0
        } else {
            4.0
        };
        let gap = theme.tokens.legend_gap.clamp(1.0, 4.0);
        let ring_room = if track.trait_columns.is_empty() {
            0.0
        } else {
            track
                .trait_columns
                .iter()
                .map(|column| column.ring_width)
                .sum::<f64>()
                + gap * track.trait_columns.len() as f64
        };
        let half = (area.w.min(area.h) / 2.0 - 4.0).max(2.0);
        let (tree_outer, tree_inner, ring_outer, label_radius) = match track.radial.direction {
            RadialDirection::Outward => {
                let ring_outer = (half - label_extent).max(4.0);
                let tree_outer = (ring_outer - ring_room).max(2.0);
                let tree_inner = tree_outer * track.radial.inner_radius;
                (tree_outer, tree_inner, ring_outer, ring_outer + 4.0)
            }
            RadialDirection::Inward => {
                let ring_outer = half;
                let tree_outer = (ring_outer - ring_room).max(2.0);
                let requested = tree_outer * track.radial.inner_radius;
                let tree_inner = if track.show_tips {
                    requested.max(label_extent + 4.0).min(tree_outer * 0.9)
                } else {
                    requested
                };
                (
                    tree_outer,
                    tree_inner,
                    ring_outer,
                    (tree_inner - 4.0).max(0.0),
                )
            }
        };
        RadialGeometry {
            cx: area.x + area.w / 2.0,
            cy: area.y + area.h / 2.0,
            tree_inner,
            tree_outer,
            ring_outer,
            label_radius,
            start: track.radial.start_degrees.to_radians(),
            sweep: track.radial.sweep_degrees.to_radians(),
            terminals: scene.terminals.len().max(1),
            direction: track.radial.direction,
        }
    }

    pub(super) fn angle(&self, row: f64) -> f64 {
        if self.terminals == 1 {
            return if self.full_circle() {
                self.start
            } else {
                self.start + self.sweep / 2.0
            };
        }
        let denominator = if self.full_circle() {
            self.terminals as f64
        } else {
            (self.terminals - 1) as f64
        };
        self.start + self.sweep * row / denominator
    }

    pub(super) fn angular_step(&self) -> f64 {
        if self.terminals <= 1 {
            self.sweep
        } else if self.full_circle() {
            self.sweep / self.terminals as f64
        } else {
            self.sweep / (self.terminals - 1) as f64
        }
    }

    pub(super) fn full_circle(&self) -> bool {
        self.sweep >= std::f64::consts::TAU - 1e-6
    }

    pub(super) fn point(&self, radius: f64, angle: f64) -> (f64, f64) {
        (
            self.cx + angle.cos() * radius,
            self.cy + angle.sin() * radius,
        )
    }

    pub(super) fn radius(&self, scene: &TreeScene, value: f64) -> f64 {
        let fraction = scene.fraction(value);
        match self.direction {
            RadialDirection::Outward => {
                self.tree_inner + fraction * (self.tree_outer - self.tree_inner)
            }
            RadialDirection::Inward => {
                self.tree_outer - fraction * (self.tree_outer - self.tree_inner)
            }
        }
    }

    pub(super) fn terminal_boundary(&self) -> f64 {
        match self.direction {
            RadialDirection::Outward => self.tree_outer,
            RadialDirection::Inward => self.tree_inner,
        }
    }
}

pub(super) fn draw_radial_track(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    let color = track
        .color
        .clone()
        .unwrap_or_else(|| ctx.theme.foreground.clone());
    let scene = TreeScene::new(
        &track.tree,
        track.shape,
        track.time.as_ref(),
        &track.collapsed,
    );
    let header_room = track.annotation_header_room();
    let area = Rect {
        x: ctx.band.x,
        y: ctx.band.y + header_room,
        w: ctx.band.w,
        h: (ctx.band.h - header_room).max(1.0),
    };
    let geometry = RadialGeometry::new(track, ctx.theme, &scene, area);
    let colors = branch_colors(
        &track.tree,
        &scene,
        track.color_by.as_deref(),
        ctx.theme,
        &color,
    );
    let styles = branch_styles(
        &track.tree,
        &colors,
        track.dnds.as_ref(),
        ctx.theme,
        track.line_width,
    );

    draw_radial_clade_highlights(track, ctx, &scene, &geometry);
    if !track.homoplasy_layers.is_empty() {
        let points: Vec<(usize, (f64, f64))> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| {
                let parent = track.tree.nodes()[placement.node].parent?;
                let parent_placement = scene.placements[parent]?;
                let angle = geometry.angle(placement.row);
                let middle_depth = (parent_placement.depth + placement.depth) / 2.0;
                Some((
                    placement.node,
                    geometry.point(geometry.radius(&scene, middle_depth), angle),
                ))
            })
            .collect();
        draw_homoplasy_links(
            ctx,
            &track.tree,
            &track.homoplasy_layers,
            &points,
            LinkGeometry::Centred {
                centre: (geometry.cx, geometry.cy),
            },
        );
    }
    if let Some(time) = track.time.as_ref().filter(|time| time.show_axis) {
        draw_radial_time_axis(ctx, &scene, &geometry, time);
    }
    draw_radial_padding(track, ctx, &scene, &geometry);
    draw_radial_branches(track, ctx, &scene, &geometry, &styles, &colors);
    draw_radial_node_glyphs(track, ctx, &scene, &geometry);
    if track.show_root {
        if let Some(root) = scene.placements[track.tree.root()] {
            let (x, y) = geometry.point(
                geometry.radius(&scene, root.depth),
                geometry.angle(root.row),
            );
            draw_root_marker(ctx, x, y);
        }
    }
    draw_radial_collapsed(track, ctx, &scene, &geometry, &styles);
    draw_trait_rings(track, ctx, &scene, &geometry);
    draw_radial_labels(track, ctx, &scene, &geometry);
    draw_trait_ring_headings(track, ctx);
    if let Some(bar) = track.branch_scale() {
        draw_radial_scale_bar(ctx, &scene, &geometry, area, bar);
    }
    draw_annotation_legend(track, ctx);
}

pub(super) fn draw_radial_padding(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if !track.show_tips && track.trait_columns.is_empty() {
        return;
    }
    let boundary = geometry.terminal_boundary();
    for (row, node) in scene.terminals.iter().enumerate() {
        if !track.tree.nodes()[*node].is_leaf() {
            continue;
        }
        let placement = scene.placements[*node].unwrap();
        let radius = geometry.radius(scene, placement.depth);
        if (radius - boundary).abs() <= 0.5 {
            continue;
        }
        let angle = geometry.angle(row as f64);
        let (x0, y0) = geometry.point(radius, angle);
        let (x1, y1) = geometry.point(boundary, angle);
        ctx.svg
            .line(x0, y0, x1, y1, &ctx.theme.rule, ctx.theme.tokens.hairline);
    }
}

pub(super) fn draw_radial_branches(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    styles: &BranchStyles<'_>,
    colors: &PerNode<String>,
) {
    for placement in scene.placements.iter().flatten() {
        let node = &track.tree.nodes()[placement.node];
        let Some(parent) = node.parent else {
            continue;
        };
        let Some(parent_placement) = scene.placements[parent] else {
            continue;
        };
        let angle = geometry.angle(placement.row);
        let (x0, y0) = geometry.point(geometry.radius(scene, parent_placement.depth), angle);
        let (x1, y1) = geometry.point(geometry.radius(scene, placement.depth), angle);
        let title = branch_title(
            &track.tree,
            placement.node,
            track.color_by.as_deref(),
            track.dnds.as_ref(),
            track
                .branch_labels
                .as_ref()
                .map(|labels| labels.key.as_str()),
            !track.show_tips,
            false,
        );
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        let style = styles.get(placement.node);
        ctx.svg
            .line_pattern(x0, y0, x1, y1, &style.color, style.width, style.pattern);
        if title.is_some() {
            ctx.svg.end_group();
        }
        if let Some(labels) = &track.branch_labels {
            draw_branch_annotation(ctx, &track.tree, placement.node, labels, (x0, y0), (x1, y1));
        }
        draw_branch_rate_mixtures(
            ctx,
            &track.tree,
            placement.node,
            &track.rate_mixtures,
            (x0, y0),
            (x1, y1),
        );
        draw_branch_event_layers(
            ctx,
            &track.tree,
            placement.node,
            &track.branch_event_layers,
            (x0, y0),
            (x1, y1),
        );
        draw_branch_intervals(
            ctx,
            &track.tree,
            placement.node,
            &track.branch_interval_layers,
            (x0, y0),
            (x1, y1),
        );
        draw_ancestral_transitions(
            ctx,
            &track.tree,
            placement.node,
            &track.ancestral_state_layers,
            (x0, y0),
            (x1, y1),
        );
    }

    for placement in scene.placements.iter().flatten() {
        let node = &track.tree.nodes()[placement.node];
        if node.is_leaf() || scene.terminals.contains(&placement.node) {
            continue;
        }
        let angles: Vec<f64> = node
            .children
            .iter()
            .filter_map(|child| scene.placements[*child])
            .map(|child| geometry.angle(child.row))
            .collect();
        let (Some(start), Some(end)) = (angles.first(), angles.last()) else {
            continue;
        };
        let radius = geometry.radius(scene, placement.depth);
        if radius > 0.5 && (end - start).abs() > 1e-9 {
            let title = branch_title(
                &track.tree,
                placement.node,
                track.color_by.as_deref(),
                None,
                None,
                false,
                true,
            );
            if let Some(title) = &title {
                ctx.svg.begin_titled(title);
            }
            let connector = connector_style(
                track.dnds.as_ref(),
                ctx.theme,
                colors.get(placement.node),
                track.line_width,
            );
            ctx.svg.path_stroked_pattern(
                &radial_arc_path(geometry, radius, *start, *end),
                &connector.color,
                connector.width,
                connector.pattern,
            );
            if title.is_some() {
                ctx.svg.end_group();
            }
        }
        let angle = geometry.angle(placement.row);
        let (x, y) = geometry.point(radius, angle);
        if let Some(support) = node.support.filter(|value| {
            track.support_style != SupportStyle::None
                && support_fraction(*value).is_some_and(|value| value >= track.support_threshold)
        }) {
            draw_support(
                ctx,
                x,
                y,
                support,
                &styles.get(placement.node).color,
                track.support_style,
            );
        } else if track.show_nodes {
            ctx.svg.circle_ringed(
                x,
                y,
                ctx.theme.tokens.marker_radius * 0.65,
                &styles.get(placement.node).color,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
        }
    }
}

pub(super) fn draw_radial_labels(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if !track.show_tips {
        return;
    }
    let size = ctx.theme.font_size - 1.0;
    for (row, node) in scene.terminals.iter().enumerate() {
        let angle = geometry.angle(row as f64);
        let (x, y) = geometry.point(geometry.label_radius, angle);
        let degrees = angle.to_degrees().rem_euclid(360.0);
        let right = angle.cos() >= 0.0;
        let (rotation, anchor) = match (geometry.direction, right) {
            (RadialDirection::Outward, true) => (degrees, crate::svg::Anchor::Start),
            (RadialDirection::Outward, false) => (degrees + 180.0, crate::svg::Anchor::End),
            (RadialDirection::Inward, true) => (degrees, crate::svg::Anchor::End),
            (RadialDirection::Inward, false) => (degrees + 180.0, crate::svg::Anchor::Start),
        };
        ctx.svg.text_rotated(
            (x, y + size * 0.32),
            rotation,
            &terminal_label(&track.tree, *node, &track.collapsed),
            &ctx.theme.muted,
            size,
            anchor,
        );
    }
}

pub(super) fn radial_arc_path(
    geometry: &RadialGeometry,
    radius: f64,
    start: f64,
    end: f64,
) -> String {
    let delta = (end - start).abs();
    let (x0, y0) = geometry.point(radius, start);
    if delta >= std::f64::consts::TAU - 1e-6 {
        let middle = start + std::f64::consts::PI;
        let (xm, ym) = geometry.point(radius, middle);
        return format!(
            "M {} {} A {} {} 0 0 1 {} {} A {} {} 0 0 1 {} {}",
            num(x0),
            num(y0),
            num(radius),
            num(radius),
            num(xm),
            num(ym),
            num(radius),
            num(radius),
            num(x0),
            num(y0)
        );
    }
    let (x1, y1) = geometry.point(radius, end);
    format!(
        "M {} {} A {} {} 0 {} 1 {} {}",
        num(x0),
        num(y0),
        num(radius),
        num(radius),
        usize::from(delta > std::f64::consts::PI),
        num(x1),
        num(y1)
    )
}

pub(super) fn draw_radial_time_axis(
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    time: &TimeAxis,
) {
    if !scene.temporal {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    for index in 0..=2 {
        let fraction = index as f64 / 2.0;
        let value = scene.minimum + fraction * (scene.maximum - scene.minimum);
        let radius = geometry.radius(scene, value);
        if radius > 0.5 {
            ctx.svg.path_stroked(
                &radial_arc_path(
                    geometry,
                    radius,
                    geometry.start,
                    geometry.start + geometry.sweep,
                ),
                &ctx.theme.rule,
                ctx.theme.tokens.hairline,
            );
        }
        let label = match &time.unit {
            Some(unit) => format!("{} {unit}", text_rounded(value, 3)),
            None => text_rounded(value, 3),
        };
        let (x, y) = geometry.point((radius - 4.0).max(0.0), geometry.start);
        let rotation = upright_tangent(geometry.start);
        ctx.svg.text_rotated(
            (x, y - 2.0),
            rotation,
            &label,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );
    }
}

pub(super) fn draw_radial_collapsed(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    styles: &BranchStyles<'_>,
) {
    for (row, node) in scene.terminals.iter().enumerate() {
        if track.tree.nodes()[*node].is_leaf() {
            continue;
        }
        let placement = scene.placements[*node].unwrap();
        let start_radius = geometry.radius(scene, placement.depth);
        let descendant_radii: Vec<f64> = track
            .tree
            .descendants(*node)
            .into_iter()
            .map(|descendant| geometry.radius(scene, scene.source_placements[descendant].depth))
            .collect();
        let far_radius = match geometry.direction {
            RadialDirection::Outward => descendant_radii
                .into_iter()
                .fold(start_radius + 2.0, f64::max),
            RadialDirection::Inward => descendant_radii
                .into_iter()
                .fold((start_radius - 2.0).max(0.0), f64::min),
        };
        let angle = geometry.angle(row as f64);
        let half = (geometry.angular_step() * 0.34).min(std::f64::consts::PI * 0.24);
        let (tip_x, tip_y) = geometry.point(start_radius, angle);
        let (left_x, left_y) = geometry.point(far_radius, angle - half);
        let (right_x, right_y) = geometry.point(far_radius, angle + half);
        let d = format!(
            "M {} {} L {} {} A {} {} 0 0 1 {} {} Z",
            num(tip_x),
            num(tip_y),
            num(left_x),
            num(left_y),
            num(far_radius),
            num(far_radius),
            num(right_x),
            num(right_y)
        );
        let title = collapsed_title(&track.tree, *node);
        ctx.svg.begin_titled(&title);
        ctx.svg.path(&d, &styles.get(*node).color, 0.28);
        ctx.svg.end_group();
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn draw_trait_sector(
    ctx: &mut DrawContext<'_>,
    column: &TraitColumn,
    value: Option<&AnnotationValue>,
    domain: &TraitDomain,
    title: &str,
    centre: (f64, f64),
    inner: f64,
    outer: f64,
    start: f64,
    end: f64,
    angle: f64,
) {
    let path = annular_sector_path(centre.0, centre.1, inner, outer, start, end);
    let fill = domain.color(column, value, ctx.theme);
    let middle_radius = (inner + outer) / 2.0;
    let (x, y) = (
        centre.0 + angle.cos() * middle_radius,
        centre.1 + angle.sin() * middle_radius,
    );
    let thickness = (outer - inner).max(0.0);
    let arc_room = middle_radius * (end - start).abs();
    let marker_radius = (thickness.min(arc_room) * 0.28).clamp(1.4, 5.5);
    ctx.svg.begin_titled(title);
    match column.style {
        TraitStyle::Strip => {
            if let Some(fill) = &fill {
                ctx.svg.path(&path, fill, 1.0);
            } else {
                ctx.svg
                    .path_stroked(&path, &ctx.theme.rule, ctx.theme.tokens.hairline);
            }
        }
        TraitStyle::Bar => {
            ctx.svg
                .path_stroked(&path, &ctx.theme.rule, ctx.theme.tokens.hairline);
            if let Some(fraction) = domain.fraction(value) {
                let bar_outer = inner + thickness * fraction;
                if bar_outer > inner + 0.1 {
                    let bar = annular_sector_path(centre.0, centre.1, inner, bar_outer, start, end);
                    ctx.svg
                        .path(&bar, fill.as_deref().unwrap_or(&ctx.theme.accent), 0.92);
                }
            }
        }
        TraitStyle::Binary => match binary_state(value) {
            Some(true) => ctx.svg.circle_ringed(
                x,
                y,
                marker_radius,
                &ctx.theme.accent,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            ),
            Some(false) => ctx.svg.circle_ringed(
                x,
                y,
                (marker_radius * 0.42).max(1.0),
                &ctx.theme.rule,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            ),
            None => ctx
                .svg
                .path_stroked(&path, &ctx.theme.rule, ctx.theme.tokens.hairline),
        },
        TraitStyle::Symbol => {
            if let Some(index) = domain.category(value) {
                ctx.svg.symbol_ringed(
                    x,
                    y,
                    marker_radius,
                    ctx.theme.symbol(index),
                    fill.as_deref().unwrap_or(&ctx.theme.accent),
                    &ctx.theme.background,
                    ctx.theme.tokens.hairline,
                );
            } else {
                ctx.svg
                    .path_stroked(&path, &ctx.theme.rule, ctx.theme.tokens.hairline);
            }
        }
    }
    if column.show_values && matches!(column.style, TraitStyle::Strip | TraitStyle::Bar) {
        let text = value
            .map(ToString::to_string)
            .unwrap_or_else(|| crate::tree::ABSENT.to_string());
        let size = (ctx.theme.font_size - 3.0).max(6.0);
        if thickness >= size + 1.0 && arc_room >= text_width(&text, size) + 4.0 {
            let ink = fill
                .as_deref()
                .map(contrast_ink)
                .unwrap_or(ctx.theme.muted.as_str());
            ctx.svg.text_rotated(
                (x, y + size * 0.3),
                upright_tangent(angle),
                &text,
                ink,
                size,
                crate::svg::Anchor::Middle,
            );
        }
    }
    ctx.svg.end_group();
}

pub(super) fn draw_trait_rings(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    if track.trait_columns.is_empty() {
        return;
    }
    let gap = ctx.theme.tokens.legend_gap.clamp(1.0, 4.0);
    let mut inner = geometry.tree_outer + gap;
    for column in &track.trait_columns {
        let outer = (inner + column.ring_width).min(geometry.ring_outer);
        let values: Vec<Option<&AnnotationValue>> = scene
            .terminals
            .iter()
            .map(|node| inherited_annotation(&track.tree, *node, &column.key))
            .collect();
        let domain = TraitDomain::new(scene.placements.iter().flatten().filter_map(|placement| {
            inherited_annotation(&track.tree, placement.node, &column.key)
        }));
        for (row, node) in scene.terminals.iter().enumerate() {
            let angle = geometry.angle(row as f64);
            let gap_angle = if outer > 0.0 { 0.8 / outer } else { 0.0 };
            let half = (geometry.angular_step() / 2.0 - gap_angle)
                .max(geometry.angular_step() * 0.12)
                .min(std::f64::consts::PI * 0.45);
            let start = if geometry.full_circle() {
                angle - half
            } else {
                (angle - half).max(geometry.start)
            };
            let end = if geometry.full_circle() {
                angle + half
            } else {
                (angle + half).min(geometry.start + geometry.sweep)
            };
            let value = values[row];
            let displayed = value.map(ToString::to_string);
            let name = terminal_label(&track.tree, *node, &track.collapsed);
            let title = match &displayed {
                Some(value) => format!("{name}; {} {value}", column.key),
                None => format!("{name}; {} missing", column.key),
            };
            draw_trait_sector(
                ctx,
                column,
                value,
                &domain,
                &title,
                (geometry.cx, geometry.cy),
                inner,
                outer,
                start,
                end,
                angle,
            );
        }
        inner = outer + gap;
    }
}

pub(super) fn annular_sector_path(
    cx: f64,
    cy: f64,
    inner: f64,
    outer: f64,
    start: f64,
    end: f64,
) -> String {
    let point = |radius: f64, angle: f64| (cx + angle.cos() * radius, cy + angle.sin() * radius);
    let (x0, y0) = point(outer, start);
    let (x1, y1) = point(outer, end);
    let (x2, y2) = point(inner, end);
    let (x3, y3) = point(inner, start);
    let large = usize::from((end - start).abs() > std::f64::consts::PI);
    format!(
        "M {} {} A {} {} 0 {} 1 {} {} L {} {} A {} {} 0 {} 0 {} {} Z",
        num(x0),
        num(y0),
        num(outer),
        num(outer),
        large,
        num(x1),
        num(y1),
        num(x2),
        num(y2),
        num(inner),
        num(inner),
        large,
        num(x3),
        num(y3)
    )
}

pub(super) fn draw_trait_ring_headings(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    if track.trait_columns.is_empty() {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let slot = ctx.band.w / track.trait_columns.len() as f64;
    for (index, column) in track.trait_columns.iter().enumerate() {
        let x = ctx.band.x + index as f64 * slot;
        let visible = fit_text(&column.label, (slot - 18.0).max(1.0), size);
        if visible != column.label {
            ctx.svg.begin_titled(&column.label);
        }
        match column.style {
            TraitStyle::Strip => {
                ctx.svg
                    .rect_rounded(x + 3.0, ctx.band.y + 3.0, 8.0, 8.0, 1.5, &ctx.theme.accent)
            }
            TraitStyle::Bar => {
                ctx.svg.rect_outline(
                    x + 2.0,
                    ctx.band.y + 3.0,
                    10.0,
                    8.0,
                    &ctx.theme.rule,
                    ctx.theme.tokens.hairline,
                );
                ctx.svg
                    .rect(x + 2.0, ctx.band.y + 6.0, 7.0, 5.0, &ctx.theme.accent);
            }
            TraitStyle::Binary => ctx.svg.circle_ringed(
                x + 7.0,
                ctx.band.y + 7.0,
                3.6,
                &ctx.theme.accent,
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            ),
            TraitStyle::Symbol => ctx.svg.symbol_ringed(
                x + 7.0,
                ctx.band.y + 7.0,
                3.8,
                ctx.theme.symbol(index),
                ctx.theme.color(index),
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            ),
        }
        ctx.svg.text(
            x + 16.0,
            ctx.band.y + size + 2.0,
            &visible,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Start,
        );
        if visible != column.label {
            ctx.svg.end_group();
        }
    }
}

pub(super) fn upright_tangent(angle: f64) -> f64 {
    let mut degrees = (angle.to_degrees() + 90.0).rem_euclid(360.0);
    if degrees > 90.0 && degrees < 270.0 {
        degrees += 180.0;
    }
    degrees
}
