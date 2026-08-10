//! What is drawn on top of a tree rather than being part of one.
//!
//! Clade highlights and node glyphs are the two decorations that have to exist
//! three times over, once per projection, because a band across rows, a sector
//! of a circle and a hull around a cloud of points are the same statement in
//! three coordinate systems. Keeping them together is what makes it visible
//! when one of the three drifts from the other two.

use super::*;

pub(super) fn node_in_clade(tree: &Tree, node: usize, clade: usize) -> bool {
    node == clade || tree.ancestors(node).contains(&clade)
}

pub(super) fn highlight_title(tree: &Tree, highlight: &CladeHighlight) -> String {
    let name = highlight
        .label
        .as_deref()
        .or_else(|| tree.nodes()[highlight.node].name.as_deref())
        .unwrap_or("clade");
    format!("{name}; {} tips", tree.clade_size(highlight.node))
}

pub(super) fn highlight_color(highlight: &CladeHighlight, index: usize, theme: &Theme) -> String {
    highlight
        .color
        .clone()
        .unwrap_or_else(|| theme.color(index).to_string())
}

pub(super) fn draw_highlight_label(
    ctx: &mut DrawContext<'_>,
    highlight: &CladeHighlight,
    at: (f64, f64),
    available: f64,
) {
    let Some(label) = &highlight.label else {
        return;
    };
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let visible = fit_text(label, available.max(0.0), size);
    ctx.svg.text_bold(
        at.0,
        at.1,
        &visible,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
}

pub(super) fn draw_rectangular_clade_highlights(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    area: Rect,
) {
    for (index, highlight) in track.clade_highlights.iter().enumerate() {
        let Some(placement) = scene
            .placements
            .get(highlight.node)
            .and_then(|placement| *placement)
        else {
            continue;
        };
        let rows: Vec<usize> = scene
            .terminals
            .iter()
            .enumerate()
            .filter_map(|(row, node)| {
                node_in_clade(&track.tree, *node, highlight.node).then_some(row)
            })
            .collect();
        let (Some(first), Some(last)) = (rows.first(), rows.last()) else {
            continue;
        };
        let x = (scene.x(area, placement.depth) - 3.0).max(area.x);
        let y = area.y + *first as f64 * track.row_height + 1.0;
        let bottom = area.y + (*last + 1) as f64 * track.row_height - 1.0;
        let color = highlight_color(highlight, index, ctx.theme);
        ctx.svg
            .begin_titled(&highlight_title(&track.tree, highlight));
        ctx.svg.rect_rounded_opacity(
            x,
            y,
            (area.right() - x).max(1.0),
            (bottom - y).max(1.0),
            ctx.theme.corner_radius * 1.6,
            &color,
            highlight.opacity,
        );
        ctx.svg.rect_rounded(
            x + 1.5,
            y + 3.0,
            2.5,
            (bottom - y - 6.0).max(1.0),
            1.25,
            &color,
        );
        draw_highlight_label(
            ctx,
            highlight,
            (x + 8.0, y + (ctx.theme.font_size - 2.0).max(6.0) + 2.0),
            area.right() - x - 12.0,
        );
        ctx.svg.end_group();
    }
}

pub(super) fn draw_radial_clade_highlights(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    for (index, highlight) in track.clade_highlights.iter().enumerate() {
        let Some(placement) = scene
            .placements
            .get(highlight.node)
            .and_then(|placement| *placement)
        else {
            continue;
        };
        let rows: Vec<usize> = scene
            .terminals
            .iter()
            .enumerate()
            .filter_map(|(row, node)| {
                node_in_clade(&track.tree, *node, highlight.node).then_some(row)
            })
            .collect();
        let (Some(first), Some(last)) = (rows.first(), rows.last()) else {
            continue;
        };
        let half = geometry.angular_step() * 0.48;
        let start = if geometry.full_circle() && rows.len() == scene.terminals.len() {
            geometry.start
        } else {
            (geometry.angle(*first as f64) - half).max(geometry.start)
        };
        let mut end = if geometry.full_circle() && rows.len() == scene.terminals.len() {
            geometry.start + geometry.sweep - 1e-5
        } else {
            (geometry.angle(*last as f64) + half).min(geometry.start + geometry.sweep)
        };
        if end <= start {
            end = start + 1e-5;
        }
        let root_radius = geometry.radius(scene, placement.depth);
        let terminal = geometry.terminal_boundary();
        let inner = root_radius.min(terminal).max(0.5);
        let outer = root_radius.max(terminal).max(inner + 0.5);
        let path = annular_sector_path(
            geometry.cx,
            geometry.cy,
            (inner - 3.0).max(0.5),
            outer + 3.0,
            start,
            end,
        );
        let color = highlight_color(highlight, index, ctx.theme);
        let outline = mix(
            ctx.theme.surface(),
            &color,
            (highlight.opacity * 3.0).min(0.5),
        );
        ctx.svg
            .begin_titled(&highlight_title(&track.tree, highlight));
        ctx.svg.path(&path, &color, highlight.opacity);
        ctx.svg
            .path_stroked(&path, &outline, ctx.theme.tokens.hairline.max(0.8));
        let angle = (start + end) / 2.0;
        let radius = (inner + outer) / 2.0;
        let (x, y) = geometry.point(radius, angle);
        if let Some(label) = &highlight.label {
            let size = (ctx.theme.font_size - 2.0).max(6.0);
            ctx.svg.text_rotated(
                (x, y),
                upright_tangent(angle),
                label,
                &ctx.theme.muted,
                size,
                crate::svg::Anchor::Middle,
            );
        }
        ctx.svg.end_group();
    }
}

pub(super) fn draw_unrooted_clade_highlights(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &UnrootedScene,
    geometry: &UnrootedGeometry,
) {
    for (index, highlight) in track.clade_highlights.iter().enumerate() {
        if !scene.visible.contains(&highlight.node) {
            continue;
        }
        let points: Vec<(f64, f64)> = scene
            .visible
            .iter()
            .filter(|node| node_in_clade(&track.tree, **node, highlight.node))
            .filter_map(|node| scene.positions[*node].map(|point| geometry.node(point)))
            .collect();
        if points.is_empty() {
            continue;
        }
        let mut hull = convex_hull(points);
        let centre = hull
            .iter()
            .fold((0.0, 0.0), |sum, point| (sum.0 + point.0, sum.1 + point.1));
        let centre = (centre.0 / hull.len() as f64, centre.1 / hull.len() as f64);
        for point in &mut hull {
            let dx = point.0 - centre.0;
            let dy = point.1 - centre.1;
            let distance = dx.hypot(dy).max(1.0);
            point.0 += dx / distance * 7.0;
            point.1 += dy / distance * 7.0;
        }
        let color = highlight_color(highlight, index, ctx.theme);
        let fill = mix(ctx.theme.surface(), &color, highlight.opacity);
        let outline = mix(
            ctx.theme.surface(),
            &color,
            (highlight.opacity * 3.0).min(0.5),
        );
        ctx.svg
            .begin_titled(&highlight_title(&track.tree, highlight));
        if hull.len() >= 3 {
            let mut path = format!("M {} {}", num(hull[0].0), num(hull[0].1));
            for point in hull.iter().skip(1) {
                path.push_str(&format!(" L {} {}", num(point.0), num(point.1)));
            }
            path.push_str(" Z");
            ctx.svg.path(&path, &color, highlight.opacity);
            ctx.svg
                .path_stroked(&path, &outline, ctx.theme.tokens.hairline.max(0.8));
        } else {
            ctx.svg.circle_ringed(
                centre.0,
                centre.1,
                8.0,
                &fill,
                &outline,
                ctx.theme.tokens.hairline.max(0.8),
            );
        }
        let top = hull
            .iter()
            .min_by(|left, right| left.1.total_cmp(&right.1))
            .copied()
            .unwrap_or(centre);
        draw_highlight_label(ctx, highlight, (top.0 + 3.0, top.1 + 10.0), 90.0);
        ctx.svg.end_group();
    }
}

pub(super) fn convex_hull(mut points: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    points.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });
    points.dedup_by(|left, right| left.0 == right.0 && left.1 == right.1);
    if points.len() <= 2 {
        return points;
    }
    let cross = |origin: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - origin.0) * (b.1 - origin.1) - (a.1 - origin.1) * (b.0 - origin.0)
    };
    let mut lower = Vec::new();
    for point in &points {
        while lower.len() >= 2
            && cross(lower[lower.len() - 2], lower[lower.len() - 1], *point) <= 0.0
        {
            lower.pop();
        }
        lower.push(*point);
    }
    let mut upper = Vec::new();
    for point in points.iter().rev() {
        while upper.len() >= 2
            && cross(upper[upper.len() - 2], upper[upper.len() - 1], *point) <= 0.0
        {
            upper.pop();
        }
        upper.push(*point);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

pub(super) fn glyph_matches(tree: &Tree, node: usize, target: NodeGlyphTarget) -> bool {
    match target {
        NodeGlyphTarget::All => true,
        NodeGlyphTarget::Internal => !tree.nodes()[node].is_leaf(),
        NodeGlyphTarget::Leaves => tree.nodes()[node].is_leaf(),
    }
}

pub(super) fn glyph_values(tree: &Tree, node: usize, glyph: &NodeGlyph) -> Option<Vec<f64>> {
    glyph
        .keys
        .iter()
        .map(|key| {
            tree.annotation(node, key)
                .and_then(AnnotationValue::as_number)
                .filter(|value| value.is_finite() && *value >= 0.0)
        })
        .collect()
}

pub(super) fn node_title(tree: &Tree, node: usize) -> String {
    tree.nodes()[node].name.clone().unwrap_or_else(|| {
        if node == tree.root() {
            "root".to_string()
        } else {
            format!("node {node}")
        }
    })
}

pub(super) fn glyph_title(tree: &Tree, node: usize, glyph: &NodeGlyph, values: &[f64]) -> String {
    let values = glyph
        .keys
        .iter()
        .zip(values)
        .map(|(key, value)| format!("{key} {}", num(*value)))
        .collect::<Vec<_>>()
        .join("; ");
    format!("{}; {values}", node_title(tree, node))
}

pub(super) fn bubble_max<'a>(
    tree: &Tree,
    nodes: impl IntoIterator<Item = &'a usize>,
    glyph: &NodeGlyph,
) -> f64 {
    nodes
        .into_iter()
        .filter(|node| glyph_matches(tree, **node, glyph.target))
        .filter_map(|node| glyph_values(tree, *node, glyph))
        .filter_map(|values| values.first().copied())
        .fold(0.0f64, f64::max)
}

pub(super) fn draw_node_glyph(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    glyph: &NodeGlyph,
    glyph_index: usize,
    at: (f64, f64),
    maximum: f64,
) {
    if !glyph_matches(tree, node, glyph.target) {
        return;
    }
    let Some(values) = glyph_values(tree, node, glyph) else {
        return;
    };
    let total: f64 = values.iter().sum();
    if glyph.style != NodeGlyphStyle::Bubble && total <= 0.0 {
        return;
    }
    ctx.svg
        .begin_titled(&glyph_title(tree, node, glyph, &values));
    match glyph.style {
        NodeGlyphStyle::Bubble => {
            let value = values[0];
            let radius = if value <= 0.0 || maximum <= 0.0 {
                glyph.minimum_size * 0.45
            } else {
                glyph.minimum_size + (value / maximum).sqrt() * (glyph.size - glyph.minimum_size)
            };
            ctx.svg.circle_ringed(
                at.0,
                at.1,
                radius,
                ctx.theme.color(glyph_index),
                ctx.theme.surface(),
                ctx.theme.tokens.stroke.max(1.1),
            );
        }
        NodeGlyphStyle::Pie | NodeGlyphStyle::Donut => {
            ctx.svg
                .circle(at.0, at.1, glyph.size + 1.2, ctx.theme.surface());
            let positive = values.iter().filter(|value| **value > 0.0).count();
            if positive == 1 {
                let index = values.iter().position(|value| *value > 0.0).unwrap_or(0);
                ctx.svg.circle_ringed(
                    at.0,
                    at.1,
                    glyph.size,
                    ctx.theme.color(index),
                    ctx.theme.surface(),
                    0.7,
                );
            } else {
                let mut start = -std::f64::consts::FRAC_PI_2;
                for (index, value) in values.iter().enumerate() {
                    if *value <= 0.0 {
                        continue;
                    }
                    let end = start + std::f64::consts::TAU * *value / total;
                    let path = pie_slice_path(at.0, at.1, glyph.size, start, end);
                    ctx.svg.path(&path, ctx.theme.color(index), 1.0);
                    ctx.svg.path_stroked(&path, ctx.theme.surface(), 0.7);
                    start = end;
                }
            }
            if glyph.style == NodeGlyphStyle::Donut {
                ctx.svg
                    .circle(at.0, at.1, glyph.size * 0.48, ctx.theme.surface());
            }
        }
        NodeGlyphStyle::StackedBar => {
            let width = glyph.size * 3.0;
            let height = (glyph.size * 0.78).max(3.0);
            let left = at.0 - width / 2.0;
            let top = at.1 - height / 2.0;
            ctx.svg.rect_rounded(
                left - 1.0,
                top - 1.0,
                width + 2.0,
                height + 2.0,
                2.0,
                &ctx.theme.background,
            );
            let mut cursor = left;
            for (index, value) in values.iter().enumerate() {
                let segment = width * *value / total;
                if segment > 0.0 {
                    ctx.svg
                        .rect(cursor, top, segment, height, ctx.theme.color(index));
                }
                cursor += segment;
                if cursor < left + width - 0.5 {
                    ctx.svg.line(
                        cursor,
                        top + 0.5,
                        cursor,
                        top + height - 0.5,
                        ctx.theme.surface(),
                        0.8,
                    );
                }
            }
        }
    }
    ctx.svg.end_group();
}

pub(super) fn pie_slice_path(cx: f64, cy: f64, radius: f64, start: f64, end: f64) -> String {
    let (x0, y0) = (cx + start.cos() * radius, cy + start.sin() * radius);
    let (x1, y1) = (cx + end.cos() * radius, cy + end.sin() * radius);
    format!(
        "M {} {} L {} {} A {} {} 0 {} 1 {} {} Z",
        num(cx),
        num(cy),
        num(x0),
        num(y0),
        num(radius),
        num(radius),
        usize::from((end - start).abs() > std::f64::consts::PI),
        num(x1),
        num(y1)
    )
}

pub(super) fn draw_rectangular_node_glyphs(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    area: Rect,
) {
    let visible: Vec<usize> = scene
        .placements
        .iter()
        .flatten()
        .map(|placement| placement.node)
        .collect();
    for (glyph_index, glyph) in track.node_glyphs.iter().enumerate() {
        let maximum = bubble_max(&track.tree, &visible, glyph);
        for node in &visible {
            let placement = scene.placements[*node].unwrap();
            draw_node_glyph(
                ctx,
                &track.tree,
                *node,
                glyph,
                glyph_index,
                (
                    scene.x(area, placement.depth),
                    area.y + track.row_height / 2.0 + placement.row * track.row_height,
                ),
                maximum,
            );
        }
    }
}

pub(super) fn draw_radial_node_glyphs(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
) {
    let visible: Vec<usize> = scene
        .placements
        .iter()
        .flatten()
        .map(|placement| placement.node)
        .collect();
    for (glyph_index, glyph) in track.node_glyphs.iter().enumerate() {
        let maximum = bubble_max(&track.tree, &visible, glyph);
        for node in &visible {
            let placement = scene.placements[*node].unwrap();
            let at = geometry.point(
                geometry.radius(scene, placement.depth),
                geometry.angle(placement.row),
            );
            draw_node_glyph(ctx, &track.tree, *node, glyph, glyph_index, at, maximum);
        }
    }
}

pub(super) fn draw_unrooted_node_glyphs(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &UnrootedScene,
    geometry: &UnrootedGeometry,
) {
    for (glyph_index, glyph) in track.node_glyphs.iter().enumerate() {
        let maximum = bubble_max(&track.tree, &scene.visible, glyph);
        for node in &scene.visible {
            let Some(point) = scene.positions[*node] else {
                continue;
            };
            draw_node_glyph(
                ctx,
                &track.tree,
                *node,
                glyph,
                glyph_index,
                geometry.node(point),
                maximum,
            );
        }
    }
}

pub(super) fn draw_annotation_legend(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    if track.node_glyphs.is_empty()
        && track.dnds.is_none()
        && track.rate_mixtures.is_empty()
        && track.homoplasy_layers.is_empty()
    {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let mut x = ctx.band.x + 2.0;
    let top = ctx.band.y + 1.0;
    let height = size + 7.0;
    let y = top + height / 2.0 + size * 0.34;
    let chip = mix(ctx.theme.surface(), &ctx.theme.rule, 0.32);
    if let Some(dnds) = &track.dnds {
        x = draw_dnds_legend(ctx, dnds, x, top, height, size, &chip);
    }
    for mixture in &track.rate_mixtures {
        if x >= ctx.band.right() - 10.0 {
            break;
        }
        x = draw_rate_mixture_legend(ctx, mixture, x, top, height, size, &chip);
    }
    for layer in &track.branch_event_layers {
        if x >= ctx.band.right() - 10.0 {
            break;
        }
        x = draw_branch_event_legend(ctx, layer, x, top, height, size, &chip);
    }
    for layer in &track.branch_interval_layers {
        if x >= ctx.band.right() - 10.0 {
            break;
        }
        x = draw_branch_interval_legend(ctx, layer, x, top, height, size, &chip);
    }
    for layer in &track.homoplasy_layers {
        if x >= ctx.band.right() - 10.0 {
            break;
        }
        x = draw_homoplasy_legend(ctx, layer, x, top, height, size, &chip);
    }
    for (glyph_index, glyph) in track.node_glyphs.iter().enumerate() {
        if x >= ctx.band.right() - 10.0 {
            break;
        }
        match glyph.style {
            NodeGlyphStyle::Bubble => {
                let available = (ctx.band.right() - x).min(110.0);
                let label = fit_text(&glyph.label, (available - 20.0).max(0.0), size);
                let width = (text_width(&label, size) + 24.0).min(available);
                ctx.svg
                    .rect_rounded(x, top, width, height, height / 2.0, &chip);
                ctx.svg.circle_ringed(
                    x + 9.0,
                    top + height / 2.0,
                    3.0,
                    ctx.theme.color(glyph_index),
                    ctx.theme.surface(),
                    0.7,
                );
                ctx.svg.text_bold(
                    x + 16.0,
                    y,
                    &label,
                    &ctx.theme.muted,
                    size,
                    crate::svg::Anchor::Start,
                );
                x += width + 6.0;
            }
            _ => {
                let label = fit_text(&glyph.label, 90.0, size);
                let keys: Vec<String> = glyph
                    .keys
                    .iter()
                    .map(|key| fit_text(key, 58.0, size))
                    .collect();
                let natural = 16.0
                    + text_width(&label, size)
                    + keys
                        .iter()
                        .map(|key| 16.0 + text_width(key, size))
                        .sum::<f64>();
                let width = natural.min(ctx.band.right() - x);
                ctx.svg
                    .rect_rounded(x, top, width, height, height / 2.0, &chip);
                let mut cursor = x + 8.0;
                ctx.svg.text_bold(
                    cursor,
                    y,
                    &label,
                    &ctx.theme.muted,
                    size,
                    crate::svg::Anchor::Start,
                );
                cursor += text_width(&label, size) + 8.0;
                for (key_index, key) in keys.iter().enumerate() {
                    if cursor + 12.0 >= x + width {
                        break;
                    }
                    ctx.svg.circle(
                        cursor + 3.0,
                        top + height / 2.0,
                        3.0,
                        ctx.theme.color(key_index),
                    );
                    cursor += 9.0;
                    let key = fit_text(key, (x + width - cursor - 5.0).max(0.0), size);
                    ctx.svg.text(
                        cursor,
                        y,
                        &key,
                        &ctx.theme.muted,
                        size,
                        crate::svg::Anchor::Start,
                    );
                    cursor += text_width(&key, size) + 7.0;
                }
                x += width + 6.0;
            }
        }
    }
}

fn draw_branch_event_legend(
    ctx: &mut DrawContext<'_>,
    layer: &BranchEventLayer,
    x: f64,
    top: f64,
    height: f64,
    size: f64,
    chip: &str,
) -> f64 {
    let label = fit_text(&layer.label, 88.0, size);
    let width = (text_width(&label, size) + 51.0).min((ctx.band.right() - x).max(0.0));
    if width <= 14.0 {
        return x;
    }
    let y = top + height / 2.0;
    ctx.svg.begin_titled(&format!(
        "{}; ordered direct branch events use colour and shape",
        layer.label
    ));
    ctx.svg
        .rect_rounded(x, top, width, height, height / 2.0, chip);
    for (index, symbol) in [
        crate::style::Symbol::Diamond,
        crate::style::Symbol::Circle,
        crate::style::Symbol::Square,
    ]
    .into_iter()
    .enumerate()
    {
        ctx.svg.symbol(
            x + 9.0 + index as f64 * 9.0,
            y,
            2.8,
            symbol,
            ctx.theme.color(index),
        );
    }
    ctx.svg.text_bold(
        x + 38.0,
        y + size * 0.34,
        &label,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
    ctx.svg.end_group();
    x + width + 6.0
}

fn draw_branch_interval_legend(
    ctx: &mut DrawContext<'_>,
    layer: &BranchIntervalLayer,
    x: f64,
    top: f64,
    height: f64,
    size: f64,
    chip: &str,
) -> f64 {
    let label = fit_text(&layer.label, 88.0, size);
    let width = (text_width(&label, size) + 49.0).min((ctx.band.right() - x).max(0.0));
    if width <= 14.0 {
        return x;
    }
    let y = top + height / 2.0;
    ctx.svg.begin_titled(&format!(
        "{}; point estimate with lower and upper bounds",
        layer.label
    ));
    ctx.svg
        .rect_rounded(x, top, width, height, height / 2.0, chip);
    ctx.svg
        .line(x + 8.0, y, x + 29.0, y, ctx.theme.color(0), 2.4);
    for at in [x + 8.0, x + 29.0] {
        ctx.svg
            .line(at, y - 2.5, at, y + 2.5, ctx.theme.color(0), 1.0);
    }
    ctx.svg
        .circle_ringed(x + 22.0, y, 2.4, ctx.theme.color(0), chip, 0.7);
    ctx.svg.text_bold(
        x + 37.0,
        y + size * 0.34,
        &label,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
    ctx.svg.end_group();
    x + width + 6.0
}

fn draw_rate_mixture_legend(
    ctx: &mut DrawContext<'_>,
    mixture: &BranchRateMixture,
    x: f64,
    top: f64,
    height: f64,
    size: f64,
    chip: &str,
) -> f64 {
    let label = fit_text(&mixture.label, 92.0, size);
    let width = (text_width(&label, size) + 43.0).min((ctx.band.right() - x).max(0.0));
    if width <= 14.0 {
        return x;
    }
    ctx.svg.begin_titled(&format!(
        "{}; segment width is fitted class weight and colour is omega",
        mixture.label
    ));
    ctx.svg
        .rect_rounded(x, top, width, height, height / 2.0, chip);
    let y = top + height / 2.0;
    let left = x + 8.0;
    let segment_widths = [5.0, 7.0, 10.0];
    let values = [0.25, 1.0, mixture.saturation];
    let mut cursor = left;
    ctx.svg
        .line(left, y, left + 22.0, y, ctx.theme.surface(), 6.2);
    for (segment_width, value) in segment_widths.into_iter().zip(values) {
        ctx.svg.line(
            cursor,
            y,
            cursor + segment_width,
            y,
            &omega_color(
                ctx.theme,
                value,
                mixture.neutral_lower,
                mixture.neutral_upper,
                mixture.saturation,
            ),
            4.3,
        );
        cursor += segment_width;
    }
    ctx.svg.text_bold(
        x + 35.0,
        y + size * 0.34,
        &label,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
    ctx.svg.end_group();
    x + width + 6.0
}

fn draw_homoplasy_legend(
    ctx: &mut DrawContext<'_>,
    layer: &HomoplasyLayer,
    x: f64,
    top: f64,
    height: f64,
    size: f64,
    chip: &str,
) -> f64 {
    let label = fit_text(&layer.label, 92.0, size);
    let width = (text_width(&label, size) + 42.0).min((ctx.band.right() - x).max(0.0));
    if width <= 14.0 {
        return x;
    }
    let y = top + height / 2.0;
    let color = mix(ctx.theme.surface(), ctx.theme.color(2), 0.72);
    ctx.svg.begin_titled(&format!(
        "{}; dashed curves connect recurrent direct branch events",
        layer.label
    ));
    ctx.svg
        .rect_rounded(x, top, width, height, height / 2.0, chip);
    ctx.svg.line_pattern(
        x + 8.0,
        y,
        x + 28.0,
        y,
        &color,
        layer.width,
        LinePattern::Dashed,
    );
    ctx.svg.circle(x + 8.0, y, 2.1, &color);
    ctx.svg.circle(x + 28.0, y, 2.1, &color);
    ctx.svg.text_bold(
        x + 34.0,
        y + size * 0.34,
        &label,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
    ctx.svg.end_group();
    x + width + 6.0
}

fn draw_dnds_legend(
    ctx: &mut DrawContext<'_>,
    dnds: &DnDsLayer,
    x: f64,
    top: f64,
    height: f64,
    size: f64,
    chip: &str,
) -> f64 {
    let label = fit_text(&dnds.label, 86.0, size);
    let labels = ["purifying", "near neutral", "diversifying"];
    let values = [1.0 / dnds.saturation, 1.0, dnds.saturation];
    let significance = dnds
        .significance
        .as_ref()
        .map(|test| format!("{} ≤ {}", test.key, text_rounded(test.maximum, 3)));
    let natural = 16.0
        + text_width(&label, size)
        + labels
            .iter()
            .map(|label| 15.0 + text_width(label, size))
            .sum::<f64>()
        + 20.0
        + text_width("missing", size)
        + significance
            .as_ref()
            .map_or(0.0, |label| 48.0 + text_width(label, size));
    let available = (ctx.band.right() - x).max(0.0);
    let width = natural.min(available);
    if width <= 12.0 {
        return x;
    }
    let y = top + height / 2.0 + size * 0.34;
    ctx.svg.begin_titled(&format!(
        "{}; cool branches ω < {}; neutral {}–{}; warm branches ω > {}",
        dnds.label,
        text_rounded(dnds.neutral_lower, 3),
        text_rounded(dnds.neutral_lower, 3),
        text_rounded(dnds.neutral_upper, 3),
        text_rounded(dnds.neutral_upper, 3)
    ));
    ctx.svg
        .rect_rounded(x, top, width, height, height / 2.0, chip);
    let mut cursor = x + 8.0;
    ctx.svg.text_bold(
        cursor,
        y,
        &label,
        &ctx.theme.muted,
        size,
        crate::svg::Anchor::Start,
    );
    cursor += text_width(&label, size) + 9.0;
    for (label, value) in labels.into_iter().zip(values) {
        if cursor + 15.0 >= x + width {
            break;
        }
        ctx.svg.line(
            cursor,
            top + height / 2.0,
            cursor + 9.0,
            top + height / 2.0,
            &dnds_color(dnds, ctx.theme, value),
            3.0,
        );
        cursor += 13.0;
        let visible = fit_text(label, (x + width - cursor - 4.0).max(0.0), size);
        ctx.svg.text(
            cursor,
            y,
            &visible,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Start,
        );
        cursor += text_width(&visible, size) + 7.0;
    }
    if cursor + 25.0 < x + width {
        ctx.svg.line_pattern(
            cursor,
            top + height / 2.0,
            cursor + 9.0,
            top + height / 2.0,
            &ctx.theme.rule,
            1.2,
            LinePattern::Dotted,
        );
        cursor += 13.0;
        let visible = fit_text("missing", (x + width - cursor - 4.0).max(0.0), size);
        ctx.svg.text(
            cursor,
            y,
            &visible,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Start,
        );
        cursor += text_width(&visible, size) + 7.0;
    }
    if let Some(significance) = significance {
        if cursor + 20.0 < x + width {
            ctx.svg.line(
                cursor,
                top + height / 2.0,
                cursor + 9.0,
                top + height / 2.0,
                &ctx.theme.muted,
                2.4,
            );
            cursor += 13.0;
            let visible = fit_text(&significance, (x + width - cursor - 4.0).max(0.0), size);
            ctx.svg.text(
                cursor,
                y,
                &visible,
                &ctx.theme.muted,
                size,
                crate::svg::Anchor::Start,
            );
        }
    }
    ctx.svg.end_group();
    x + width + 6.0
}
