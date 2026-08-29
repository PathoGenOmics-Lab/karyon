//! The default projection: depth across, rows down.
//!
//! This is the one a phylogeny is usually drawn in, and the one the other two
//! are defined against. The trait columns live here rather than in
//! [`decorate`](super::decorate) because a column beside a rectangular tree is
//! a strip of cells on the same rows, which is not a decoration but the
//! rectangular answer to what a ring is in polar coordinates.

use super::*;

#[derive(Debug, Clone)]
pub(super) struct BranchStyle {
    pub(super) color: String,
    pub(super) width: f64,
    pub(super) pattern: LinePattern,
}

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
    dnds: Option<&DnDsLayer>,
    show_nodes: bool,
    support_style: SupportStyle,
    support_threshold: f64,
    branch_labels: Option<&BranchLabels>,
    rate_mixtures: &[BranchRateMixture],
    homoplasy_layers: &[HomoplasyLayer],
    branch_event_layers: &[BranchEventLayer],
    branch_interval_layers: &[BranchIntervalLayer],
    ancestral_state_layers: &[AncestralStateLayer],
    branch_geometry: BranchGeometry,
    name_leaves: bool,
) {
    let colors = branch_colors(tree, scene, color_by, ctx.theme, default_color);
    let styles = branch_styles(tree, &colors, dnds, ctx.theme, width);
    let y_of = |row: f64| area.y + row_pitch / 2.0 + row * row_pitch;

    if !homoplasy_layers.is_empty() {
        let points: Vec<(usize, (f64, f64))> = scene
            .placements
            .iter()
            .flatten()
            .filter_map(|placement| {
                let parent = tree.nodes()[placement.node].parent?;
                let parent_placement = scene.placements[parent]?;
                let x0 = scene.x(area, parent_placement.depth);
                let x1 = scene.x(area, placement.depth);
                let (start, end) = rectangular_branch_endpoints(
                    branch_geometry,
                    x0,
                    y_of(parent_placement.row),
                    x1,
                    y_of(placement.row),
                );
                Some((
                    placement.node,
                    ((start.0 + end.0) / 2.0, (start.1 + end.1) / 2.0),
                ))
            })
            .collect();
        draw_homoplasy_links(
            ctx,
            tree,
            homoplasy_layers,
            &points,
            LinkGeometry::Rectangular {
                right: area.right(),
            },
        );
    }

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
        let (start, end) = rectangular_branch_endpoints(
            branch_geometry,
            x0,
            y_of(parent_placement.row),
            x1,
            y_of(placement.row),
        );
        let title = branch_title(
            tree,
            placement.node,
            color_by,
            dnds,
            branch_labels.map(|labels| labels.key.as_str()),
            name_leaves,
            false,
        );
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        let style = &styles[placement.node];
        match branch_geometry {
            BranchGeometry::Curved => ctx.svg.path_stroked_pattern(
                &rectangular_curve_path(start, end),
                &style.color,
                style.width,
                style.pattern,
            ),
            _ => ctx.svg.line_pattern(
                start.0,
                start.1,
                end.0,
                end.1,
                &style.color,
                style.width,
                style.pattern,
            ),
        }
        if title.is_some() {
            ctx.svg.end_group();
        }
        if let Some(labels) = branch_labels {
            draw_branch_annotation(ctx, tree, placement.node, labels, start, end);
        }
        draw_branch_rate_mixtures(ctx, tree, placement.node, rate_mixtures, start, end);
        draw_branch_event_layers(ctx, tree, placement.node, branch_event_layers, start, end);
        draw_branch_intervals(
            ctx,
            tree,
            placement.node,
            branch_interval_layers,
            start,
            end,
        );
        draw_ancestral_transitions(
            ctx,
            tree,
            placement.node,
            ancestral_state_layers,
            start,
            end,
        );
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
        if branch_geometry == BranchGeometry::Orthogonal {
            let title = branch_title(tree, placement.node, color_by, None, None, false, true);
            if let Some(title) = &title {
                ctx.svg.begin_titled(title);
            }
            let connector = connector_style(dnds, ctx.theme, &colors[placement.node], width);
            ctx.svg.line_pattern(
                x,
                y_of(top),
                x,
                y_of(bottom),
                &connector.color,
                connector.width,
                connector.pattern,
            );
            if title.is_some() {
                ctx.svg.end_group();
            }
        }
        if let Some(support) = node.support.filter(|value| {
            support_style != SupportStyle::None
                && support_fraction(*value).is_some_and(|value| value >= support_threshold)
        }) {
            let title = format!("clade support {}", text_rounded(support, 3));
            ctx.svg.begin_titled(&title);
            draw_support(
                ctx,
                x,
                y_of(placement.row),
                support,
                &styles[placement.node].color,
                support_style,
            );
            ctx.svg.end_group();
        } else if show_nodes {
            ctx.svg.circle_ringed(
                x,
                y_of(placement.row),
                ctx.theme.tokens.marker_radius * 0.65,
                &styles[placement.node].color,
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
                "{} ({})",
                tree.nodes()[*node].name.as_deref().unwrap_or("clade"),
                tip_count(tree.clade_size(*node))
            );
            ctx.svg.begin_titled(&title);
            ctx.svg.path(&d, &styles[*node].color, 0.28);
            ctx.svg.end_group();
        }
    }
}

fn rectangular_branch_endpoints(
    geometry: BranchGeometry,
    parent_x: f64,
    parent_y: f64,
    child_x: f64,
    child_y: f64,
) -> ((f64, f64), (f64, f64)) {
    match geometry {
        BranchGeometry::Orthogonal => ((parent_x, child_y), (child_x, child_y)),
        BranchGeometry::Diagonal | BranchGeometry::Curved => {
            ((parent_x, parent_y), (child_x, child_y))
        }
    }
}

fn rectangular_curve_path(start: (f64, f64), end: (f64, f64)) -> String {
    let middle_x = (start.0 + end.0) / 2.0;
    format!(
        "M {} {} C {} {} {} {} {} {}",
        num(start.0),
        num(start.1),
        num(middle_x),
        num(start.1),
        num(middle_x),
        num(end.1),
        num(end.0),
        num(end.1)
    )
}

pub(super) fn branch_styles(
    tree: &Tree,
    colors: &[String],
    dnds: Option<&DnDsLayer>,
    theme: &Theme,
    width: f64,
) -> Vec<BranchStyle> {
    (0..tree.nodes().len())
        .map(|node| match dnds {
            Some(dnds) => dnds_branch_style(tree, node, dnds, theme, width),
            None => BranchStyle {
                color: colors
                    .get(node)
                    .cloned()
                    .unwrap_or_else(|| theme.foreground.clone()),
                width,
                pattern: LinePattern::Solid,
            },
        })
        .collect()
}

pub(super) fn connector_style(
    dnds: Option<&DnDsLayer>,
    theme: &Theme,
    default_color: &str,
    width: f64,
) -> BranchStyle {
    BranchStyle {
        color: if dnds.is_some() {
            mix(theme.surface(), &theme.muted, 0.58)
        } else {
            default_color.to_string()
        },
        width,
        pattern: LinePattern::Solid,
    }
}

pub(super) fn dnds_branch_style(
    tree: &Tree,
    node: usize,
    dnds: &DnDsLayer,
    theme: &Theme,
    width: f64,
) -> BranchStyle {
    let Some(value) = tree
        .annotation(node, &dnds.key)
        .and_then(AnnotationValue::as_number)
        .filter(|value| value.is_finite() && *value >= 0.0)
    else {
        return BranchStyle {
            color: mix(theme.surface(), &theme.rule, 0.88),
            width: (width * 0.9).max(theme.tokens.hairline),
            pattern: LinePattern::Dotted,
        };
    };
    let significant = dnds.significance.as_ref().is_some_and(|test| {
        tree.annotation(node, &test.key)
            .and_then(AnnotationValue::as_number)
            .is_some_and(|score| score.is_finite() && score >= 0.0 && score <= test.maximum)
    });
    BranchStyle {
        color: dnds_color(dnds, theme, value),
        width: if significant {
            (width * 1.85).max(width + 0.9)
        } else {
            width
        },
        pattern: LinePattern::Solid,
    }
}

pub(super) fn dnds_color(dnds: &DnDsLayer, theme: &Theme, value: f64) -> String {
    omega_color(
        theme,
        value,
        dnds.neutral_lower,
        dnds.neutral_upper,
        dnds.saturation,
    )
}

pub(super) fn omega_color(
    theme: &Theme,
    value: f64,
    neutral_lower: f64,
    neutral_upper: f64,
    saturation: f64,
) -> String {
    let neutral = mix(theme.surface(), &theme.muted, 0.58);
    if value >= neutral_lower && value <= neutral_upper {
        return neutral;
    }
    let logarithmic_span = saturation.log2().max(f64::EPSILON);
    if value < neutral_lower {
        let strength = if value <= 0.0 {
            1.0
        } else {
            (-value.log2() / logarithmic_span).clamp(0.0, 1.0)
        };
        mix(&neutral, theme.color(0), 0.25 + strength * 0.75)
    } else {
        let strength = (value.log2() / logarithmic_span).clamp(0.0, 1.0);
        mix(&neutral, theme.color(1), 0.25 + strength * 0.75)
    }
}

pub(super) fn dnds_regime(dnds: &DnDsLayer, value: f64) -> &'static str {
    if value < dnds.neutral_lower {
        "purifying"
    } else if value > dnds.neutral_upper {
        "diversifying"
    } else {
        "approximately neutral"
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
    dnds: Option<&DnDsLayer>,
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
    if let Some(dnds) = dnds {
        match tree
            .annotation(node, &dnds.key)
            .and_then(AnnotationValue::as_number)
            .filter(|value| value.is_finite() && *value >= 0.0)
        {
            Some(value) => parts.push(format!(
                "dN/dS ω {} ({})",
                text_rounded(value, 3),
                dnds_regime(dnds, value)
            )),
            None => parts.push("dN/dS missing".to_string()),
        }
        if let Some(test) = &dnds.significance {
            if let Some(value) = tree
                .annotation(node, &test.key)
                .and_then(AnnotationValue::as_number)
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                let relation = if value <= test.maximum { "≤" } else { ">" };
                parts.push(format!(
                    "{} {} ({relation} {})",
                    test.key,
                    text_rounded(value, 3),
                    text_rounded(test.maximum, 3)
                ));
            }
        }
    }
    if let Some(key) = branch_label
        .filter(|key| Some(*key) != color_by && dnds.map_or(true, |dnds| *key != dnds.key))
    {
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
    let mut x = area.right() + tip_width + ctx.theme.tokens.label_gap;
    let names: Vec<String> = scene
        .terminals
        .iter()
        .map(|node| terminal_label(tree, *node, collapsed))
        .collect();

    for column in columns {
        // The domain is every placement and not only the terminals, because a
        // value inherited from an internal node is a value this column has and
        // a level nobody drew still has to keep its colour.
        let domain =
            TraitDomain::new(
                scene.placements.iter().flatten().filter_map(|placement| {
                    inherited_annotation(tree, placement.node, &column.key)
                }),
            );
        let rows: Vec<TraitRow<'_>> = names
            .iter()
            .zip(&scene.terminals)
            .enumerate()
            .map(|(row, (name, node))| TraitRow {
                name,
                top: area.y + row as f64 * row_pitch + 1.0,
                height: (row_pitch - 2.0).max(1.0),
                value: inherited_annotation(tree, *node, &column.key),
            })
            .collect();
        draw_column(ctx, column, &domain, x, Some(area.y - 5.0), &rows);
        x += column.width + ctx.theme.tokens.legend_gap;
    }
}
