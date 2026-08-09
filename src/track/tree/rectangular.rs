//! The default projection: depth across, rows down.
//!
//! This is the one a phylogeny is usually drawn in, and the one the other two
//! are defined against. The trait columns live here rather than in
//! [`decorate`](super::decorate) because a column beside a rectangular tree is
//! a strip of cells on the same rows, which is not a decoration but the
//! rectangular answer to what a ring is in polar coordinates.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn draw_tree_scene(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    scene: &TreeScene,
    area: Rect,
    row_pitch: f64,
    default_color: &str,
    width: f64,
    color_by: Option<&str>,
    show_nodes: bool,
    support_style: SupportStyle,
    support_threshold: f64,
    branch_labels: Option<&BranchLabels>,
    name_leaves: bool,
) {
    let colors = branch_colors(tree, scene, color_by, ctx.theme, default_color);
    let y_of = |row: f64| area.y + row_pitch / 2.0 + row * row_pitch;

    for placement in scene.placements.iter().flatten() {
        let node = &tree.nodes()[placement.node];
        let Some(parent) = node.parent else {
            continue;
        };
        let Some(parent_placement) = scene.placements[parent] else {
            continue;
        };
        let x0 = scene.x(area, parent_placement.depth);
        let x1 = scene.x(area, placement.depth);
        let y = y_of(placement.row);
        let title = branch_title(
            tree,
            placement.node,
            color_by,
            branch_labels.map(|labels| labels.key.as_str()),
            name_leaves,
            false,
        );
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        ctx.svg.line(x0, y, x1, y, &colors[placement.node], width);
        if title.is_some() {
            ctx.svg.end_group();
        }
        if let Some(labels) = branch_labels {
            draw_branch_annotation(ctx, tree, placement.node, labels, (x0, y), (x1, y));
        }
    }

    for placement in scene.placements.iter().flatten() {
        let node = &tree.nodes()[placement.node];
        if node.is_leaf() || node.children.is_empty() || scene.terminals.contains(&placement.node) {
            continue;
        }
        let rows: Vec<f64> = node
            .children
            .iter()
            .filter_map(|child| scene.placements[*child].map(|value| value.row))
            .collect();
        if rows.is_empty() {
            continue;
        }
        let (top, bottom) = rows.iter().fold((f64::MAX, f64::MIN), |(lo, hi), row| {
            (lo.min(*row), hi.max(*row))
        });
        let x = scene.x(area, placement.depth);
        let title = branch_title(tree, placement.node, color_by, None, false, true);
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        ctx.svg.line(
            x,
            y_of(top),
            x,
            y_of(bottom),
            &colors[placement.node],
            width,
        );
        if title.is_some() {
            ctx.svg.end_group();
        }
        if let Some(support) = node.support.filter(|value| {
            support_style != SupportStyle::None
                && support_fraction(*value).is_some_and(|value| value >= support_threshold)
        }) {
            draw_support(
                ctx,
                x,
                y_of(placement.row),
                support,
                &colors[placement.node],
                support_style,
            );
        } else if show_nodes {
            ctx.svg.circle_ringed(
                x,
                y_of(placement.row),
                ctx.theme.tokens.marker_radius * 0.65,
                &colors[placement.node],
                &ctx.theme.background,
                ctx.theme.tokens.hairline,
            );
        }
    }

    for node in &scene.terminals {
        if !tree.nodes()[*node].is_leaf() {
            let placement = scene.placements[*node].unwrap();
            let start = scene.x(area, placement.depth);
            let far = tree
                .descendants(*node)
                .into_iter()
                .map(|descendant| scene.x(area, scene.source_placements[descendant].depth))
                .fold(start, f64::max);
            let y = y_of(placement.row);
            let half = row_pitch * 0.34;
            let d = format!(
                "M {} {} L {} {} L {} {} Z",
                num(start),
                num(y),
                num(far.max(start + 2.0)),
                num(y - half),
                num(far.max(start + 2.0)),
                num(y + half)
            );
            let title = format!(
                "{} ({} tips)",
                tree.nodes()[*node].name.as_deref().unwrap_or("clade"),
                tree.clade_size(*node)
            );
            ctx.svg.begin_titled(&title);
            ctx.svg.path(&d, &colors[*node], 0.28);
            ctx.svg.end_group();
        }
    }
}

pub(super) fn branch_colors(
    tree: &Tree,
    scene: &TreeScene,
    key: Option<&str>,
    theme: &Theme,
    default_color: &str,
) -> Vec<String> {
    let mut colors = vec![default_color.to_string(); tree.nodes().len()];
    let Some(key) = key else {
        return colors;
    };
    let values: Vec<Option<&AnnotationValue>> = (0..tree.nodes().len())
        .map(|node| inherited_annotation(tree, node, key))
        .collect();
    let visible: Vec<usize> = scene
        .placements
        .iter()
        .flatten()
        .map(|placement| placement.node)
        .collect();
    let numeric: Vec<f64> = visible
        .iter()
        .filter_map(|node| values[*node].and_then(AnnotationValue::as_number))
        .collect();
    let all_numeric = numeric.len()
        == visible
            .iter()
            .filter(|node| values[**node].is_some())
            .count()
        && !numeric.is_empty();
    if all_numeric {
        let minimum = numeric.iter().copied().fold(f64::MAX, f64::min);
        let maximum = numeric.iter().copied().fold(f64::MIN, f64::max);
        for node in scene
            .placements
            .iter()
            .flatten()
            .map(|placement| placement.node)
        {
            if let Some(value) = values[node].and_then(AnnotationValue::as_number) {
                let fraction = if maximum <= minimum {
                    1.0
                } else {
                    (value - minimum) / (maximum - minimum)
                };
                colors[node] = mix(&theme.muted, &theme.accent, fraction);
            }
        }
    } else {
        let mut categories = BTreeMap::new();
        for node in scene
            .placements
            .iter()
            .flatten()
            .map(|placement| placement.node)
        {
            let Some(value) = values[node] else {
                continue;
            };
            let value = value.to_string();
            let next = categories.len();
            let index = *categories.entry(value).or_insert(next);
            colors[node] = theme.color(index).to_string();
        }
    }
    colors
}

pub(super) fn inherited_annotation<'a>(
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

pub(super) fn branch_title(
    tree: &Tree,
    node: usize,
    color_by: Option<&str>,
    branch_label: Option<&str>,
    name_leaf: bool,
    include_support: bool,
) -> Option<String> {
    let clade = &tree.nodes()[node];
    let mut parts = Vec::new();
    if name_leaf && clade.is_leaf() {
        if let Some(name) = &clade.name {
            if !name.is_empty() {
                parts.push(name.clone());
            }
        }
    }
    if include_support {
        if let Some(support) = clade.support.filter(|value| value.is_finite()) {
            parts.push(format!("clade support {}", text_rounded(support, 3)));
        }
    }
    if let Some(key) = color_by {
        if let Some(value) = inherited_annotation(tree, node, key) {
            parts.push(format!("{key} {value}"));
        }
    }
    if let Some(key) = branch_label.filter(|key| Some(*key) != color_by) {
        if let Some(value) = tree.annotation(node, key) {
            parts.push(format!("{key} {value}"));
        }
    }
    (!parts.is_empty()).then(|| parts.join("; "))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn draw_trait_columns(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    scene: &TreeScene,
    collapsed: &BTreeSet<usize>,
    area: Rect,
    tip_width: f64,
    columns: &[TraitColumn],
    row_pitch: f64,
) {
    if columns.is_empty() {
        return;
    }
    let size = (ctx.theme.font_size - 2.0).max(6.0);
    let mut x = area.right() + tip_width + ctx.theme.tokens.label_gap;

    for column in columns {
        let values: Vec<Option<&AnnotationValue>> = scene
            .terminals
            .iter()
            .map(|node| inherited_annotation(tree, *node, &column.key))
            .collect();
        let domain =
            TraitDomain::new(
                scene.placements.iter().flatten().filter_map(|placement| {
                    inherited_annotation(tree, placement.node, &column.key)
                }),
            );

        let heading = fit_text(&column.label, column.width, size);
        ctx.svg.text(
            x + column.width / 2.0,
            area.y - 5.0,
            &heading,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );

        for (row, node) in scene.terminals.iter().enumerate() {
            let y = area.y + row as f64 * row_pitch + 1.0;
            let height = (row_pitch - 2.0).max(1.0);
            let name = terminal_label(tree, *node, collapsed);
            let value = values[row];
            let fill = domain.color(column, value, ctx.theme);
            let displayed = value.map(ToString::to_string);
            let title = match &displayed {
                Some(value) => format!("{name}; {} {value}", column.key),
                None => format!("{name}; {} missing", column.key),
            };
            ctx.svg.begin_titled(&title);
            match column.style {
                TraitStyle::Strip => {
                    if let Some(fill) = &fill {
                        ctx.svg.rect_rounded(
                            x,
                            y,
                            column.width,
                            height,
                            ctx.theme.corner_radius.min(2.0),
                            fill,
                        );
                    } else {
                        ctx.svg.rect_outline(
                            x,
                            y,
                            column.width,
                            height,
                            &ctx.theme.rule,
                            ctx.theme.tokens.hairline,
                        );
                    }
                }
                TraitStyle::Bar => {
                    ctx.svg.rect_outline(
                        x,
                        y,
                        column.width,
                        height,
                        &ctx.theme.rule,
                        ctx.theme.tokens.hairline,
                    );
                    if let Some(fraction) = domain.fraction(value) {
                        ctx.svg.rect_rounded(
                            x,
                            y,
                            column.width * fraction,
                            height,
                            ctx.theme.corner_radius.min(2.0),
                            fill.as_deref().unwrap_or(&ctx.theme.accent),
                        );
                    }
                }
                TraitStyle::Binary => match binary_state(value) {
                    Some(true) => ctx.svg.circle_ringed(
                        x + column.width / 2.0,
                        y + height / 2.0,
                        (height * 0.28).clamp(1.4, 5.0),
                        &ctx.theme.accent,
                        &ctx.theme.background,
                        ctx.theme.tokens.hairline,
                    ),
                    Some(false) => ctx.svg.circle_ringed(
                        x + column.width / 2.0,
                        y + height / 2.0,
                        (height * 0.12).clamp(0.8, 2.0),
                        &ctx.theme.rule,
                        &ctx.theme.background,
                        ctx.theme.tokens.hairline,
                    ),
                    None => ctx.svg.rect_outline(
                        x,
                        y,
                        column.width,
                        height,
                        &ctx.theme.rule,
                        ctx.theme.tokens.hairline,
                    ),
                },
                TraitStyle::Symbol => {
                    if let Some(index) = domain.category(value) {
                        ctx.svg.symbol_ringed(
                            x + column.width / 2.0,
                            y + height / 2.0,
                            (height * 0.28).clamp(1.4, 5.0),
                            ctx.theme.symbol(index),
                            fill.as_deref().unwrap_or(&ctx.theme.accent),
                            &ctx.theme.background,
                            ctx.theme.tokens.hairline,
                        );
                    } else {
                        ctx.svg.rect_outline(
                            x,
                            y,
                            column.width,
                            height,
                            &ctx.theme.rule,
                            ctx.theme.tokens.hairline,
                        );
                    }
                }
            }
            if column.show_values && matches!(column.style, TraitStyle::Strip | TraitStyle::Bar) {
                let text = displayed.as_deref().unwrap_or("—");
                let visible = fit_text(text, column.width - 4.0, size);
                let ink = fill
                    .as_deref()
                    .filter(|_| column.style == TraitStyle::Strip)
                    .map(contrast_ink)
                    .unwrap_or(ctx.theme.muted.as_str());
                ctx.svg.text(
                    x + column.width / 2.0,
                    y + height / 2.0 + size * 0.35,
                    &visible,
                    ink,
                    size,
                    crate::svg::Anchor::Middle,
                );
            }
            ctx.svg.end_group();
        }
        x += column.width + ctx.theme.tokens.legend_gap;
    }
}
