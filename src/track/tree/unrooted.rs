//! The equal-angle projection, which has no root and therefore no rows.
//!
//! Every other projection reads [`Placement`](crate::tree::Placement), which
//! is a depth and a row measured from a root. An unrooted tree has neither, so
//! this module computes its own positions by walking out from a starting node
//! and dividing the remaining angle among the leaves below each branch.

use super::*;

#[derive(Debug, Clone)]
pub(super) struct UnrootedScene {
    pub(super) positions: Vec<Option<(f64, f64)>>,
    pub(super) parents: Vec<Option<usize>>,
    pub(super) angles: Vec<Option<f64>>,
    pub(super) terminals: Vec<usize>,
    pub(super) visible: Vec<usize>,
    pub(super) radius: f64,
}

impl UnrootedScene {
    pub(super) fn new(
        tree: &Tree,
        shape: TreeShape,
        collapsed: &BTreeSet<usize>,
        start_degrees: f64,
    ) -> Self {
        let visibility = visible_nodes(tree, collapsed);
        let visible: Vec<usize> = visibility
            .iter()
            .enumerate()
            .filter_map(|(node, visible)| visible.then_some(node))
            .collect();
        let terminal_set: BTreeSet<usize> =
            visible_terminals(tree, collapsed).into_iter().collect();
        let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); tree.nodes().len()];
        for (node, clade) in tree.nodes().iter().enumerate() {
            let Some(parent) = clade.parent else {
                continue;
            };
            if !visibility[node] || !visibility[parent] {
                continue;
            }
            let length = if shape == TreeShape::Cladogram {
                1.0
            } else {
                clade
                    .branch_length
                    .filter(|value| value.is_finite())
                    .unwrap_or(1.0)
                    .max(0.0)
            };
            adjacency[parent].push((node, length));
            adjacency[node].push((parent, length));
        }

        let candidates: Vec<usize> = visible
            .iter()
            .copied()
            .filter(|node| adjacency[*node].len() > 1)
            .collect();
        let centre = candidates
            .into_iter()
            .min_by_key(|candidate| {
                let largest = adjacency[*candidate]
                    .iter()
                    .map(|(next, _)| {
                        component_terminal_count(*next, *candidate, &adjacency, &terminal_set)
                    })
                    .max()
                    .unwrap_or(0);
                (
                    largest,
                    std::cmp::Reverse(adjacency[*candidate].len()),
                    *candidate,
                )
            })
            .or_else(|| visible.first().copied())
            .unwrap_or(tree.root());

        let mut positions = vec![None; tree.nodes().len()];
        let mut parents = vec![None; tree.nodes().len()];
        let mut angles = vec![None; tree.nodes().len()];
        positions[centre] = Some((0.0, 0.0));
        let start = start_degrees.to_radians();
        if terminal_set.contains(&centre) {
            angles[centre] = Some(start);
        }

        #[derive(Debug, Clone, Copy)]
        struct Task {
            node: usize,
            parent: Option<usize>,
            start: f64,
            end: f64,
        }

        let mut stack = vec![Task {
            node: centre,
            parent: None,
            start,
            end: start + std::f64::consts::TAU,
        }];
        while let Some(task) = stack.pop() {
            let children: Vec<(usize, f64, usize)> = adjacency[task.node]
                .iter()
                .filter(|(next, _)| Some(*next) != task.parent)
                .map(|(next, length)| {
                    (
                        *next,
                        *length,
                        component_terminal_count(*next, task.node, &adjacency, &terminal_set)
                            .max(1),
                    )
                })
                .collect();
            let total: usize = children.iter().map(|(_, _, count)| *count).sum();
            if total == 0 {
                continue;
            }
            let parent_position = positions[task.node].unwrap_or((0.0, 0.0));
            let mut cursor = task.start;
            let mut pending = Vec::with_capacity(children.len());
            for (child, length, count) in children {
                let span = (task.end - task.start) * count as f64 / total as f64;
                let end = cursor + span;
                let angle = cursor + span / 2.0;
                positions[child] = Some((
                    parent_position.0 + angle.cos() * length,
                    parent_position.1 + angle.sin() * length,
                ));
                parents[child] = Some(task.node);
                angles[child] = Some(angle);
                pending.push(Task {
                    node: child,
                    parent: Some(task.node),
                    start: cursor,
                    end,
                });
                cursor = end;
            }
            stack.extend(pending.into_iter().rev());
        }

        let mut terminals: Vec<usize> = terminal_set.into_iter().collect();
        terminals.sort_by(|left, right| {
            angles[*left]
                .unwrap_or(start)
                .total_cmp(&angles[*right].unwrap_or(start))
        });
        let radius = positions
            .iter()
            .flatten()
            .map(|(x, y)| x.hypot(*y))
            .fold(0.0f64, f64::max)
            .max(1e-9);
        UnrootedScene {
            positions,
            parents,
            angles,
            terminals,
            visible,
            radius,
        }
    }
}

pub(super) fn component_terminal_count(
    start: usize,
    blocked: usize,
    adjacency: &[Vec<(usize, f64)>],
    terminals: &BTreeSet<usize>,
) -> usize {
    let mut count = 0usize;
    let mut stack = vec![(start, blocked)];
    while let Some((node, parent)) = stack.pop() {
        count += usize::from(terminals.contains(&node));
        for (next, _) in &adjacency[node] {
            if *next != parent {
                stack.push((*next, node));
            }
        }
    }
    count
}

#[derive(Debug, Clone, Copy)]
pub(super) struct UnrootedGeometry {
    pub(super) cx: f64,
    pub(super) cy: f64,
    pub(super) scale: f64,
    pub(super) branch_radius: f64,
    pub(super) ring_outer: f64,
    pub(super) label_radius: f64,
}

impl UnrootedGeometry {
    pub(super) fn new(track: &TreeTrack, theme: &Theme, scene: &UnrootedScene, area: Rect) -> Self {
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
        let ring_outer = (half - label_extent).max(4.0);
        let branch_radius = (ring_outer - ring_room).max(2.0);
        UnrootedGeometry {
            cx: area.x + area.w / 2.0,
            cy: area.y + area.h / 2.0,
            scale: branch_radius * 0.88 / scene.radius,
            branch_radius,
            ring_outer,
            label_radius: ring_outer + 4.0,
        }
    }

    pub(super) fn node(&self, point: (f64, f64)) -> (f64, f64) {
        (
            self.cx + point.0 * self.scale,
            self.cy + point.1 * self.scale,
        )
    }

    pub(super) fn point(&self, radius: f64, angle: f64) -> (f64, f64) {
        (
            self.cx + angle.cos() * radius,
            self.cy + angle.sin() * radius,
        )
    }
}

pub(super) fn draw_unrooted_track(track: &TreeTrack, ctx: &mut DrawContext<'_>) {
    let color = track
        .color
        .clone()
        .unwrap_or_else(|| ctx.theme.foreground.clone());
    let scene = UnrootedScene::new(
        &track.tree,
        track.shape,
        &track.collapsed,
        track.radial.start_degrees,
    );
    let header_room = track.annotation_header_room();
    let area = Rect {
        x: ctx.band.x,
        y: ctx.band.y + header_room,
        w: ctx.band.w,
        h: (ctx.band.h - header_room).max(1.0),
    };
    let geometry = UnrootedGeometry::new(track, ctx.theme, &scene, area);
    let colors = unrooted_branch_colors(
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

    draw_unrooted_clade_highlights(track, ctx, &scene, &geometry);

    if !track.homoplasy_layers.is_empty() {
        let points: Vec<(usize, (f64, f64))> = scene
            .visible
            .iter()
            .filter_map(|node| {
                let parent = scene.parents[*node]?;
                let (from, to) = (scene.positions[parent]?, scene.positions[*node]?);
                let owner = if track.tree.nodes()[*node].parent == Some(parent) {
                    *node
                } else {
                    parent
                };
                let (x0, y0) = geometry.node(from);
                let (x1, y1) = geometry.node(to);
                Some((owner, ((x0 + x1) / 2.0, (y0 + y1) / 2.0)))
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

    // Terminal leaders align labels and annotation rings without pretending
    // unequal branch lengths all end at the same evolutionary distance.
    for node in &scene.terminals {
        let Some(raw) = scene.positions[*node] else {
            continue;
        };
        let angle = scene.angles[*node].unwrap_or(track.radial.start_degrees.to_radians());
        let (x0, y0) = geometry.node(raw);
        let (x1, y1) = geometry.point(geometry.branch_radius, angle);
        ctx.svg
            .line(x0, y0, x1, y1, &ctx.theme.rule, ctx.theme.tokens.hairline);
    }

    for node in &scene.visible {
        let Some(parent) = scene.parents[*node] else {
            continue;
        };
        let (Some(from), Some(to)) = (scene.positions[parent], scene.positions[*node]) else {
            continue;
        };
        let owner = if track.tree.nodes()[*node].parent == Some(parent) {
            *node
        } else {
            parent
        };
        let title = branch_title(
            &track.tree,
            owner,
            track.color_by.as_deref(),
            track.dnds.as_ref(),
            track
                .branch_labels
                .as_ref()
                .map(|labels| labels.key.as_str()),
            !track.show_tips && scene.terminals.contains(node),
            true,
        );
        if let Some(title) = &title {
            ctx.svg.begin_titled(title);
        }
        let (x0, y0) = geometry.node(from);
        let (x1, y1) = geometry.node(to);
        let style = &styles[owner];
        ctx.svg
            .line_pattern(x0, y0, x1, y1, &style.color, style.width, style.pattern);
        if title.is_some() {
            ctx.svg.end_group();
        }
        if let Some(labels) = &track.branch_labels {
            draw_branch_annotation(ctx, &track.tree, owner, labels, (x0, y0), (x1, y1));
        }
        draw_branch_rate_mixtures(
            ctx,
            &track.tree,
            owner,
            &track.rate_mixtures,
            (x0, y0),
            (x1, y1),
        );
        draw_branch_event_layers(
            ctx,
            &track.tree,
            owner,
            &track.branch_event_layers,
            (x0, y0),
            (x1, y1),
        );
        draw_branch_intervals(
            ctx,
            &track.tree,
            owner,
            &track.branch_interval_layers,
            (x0, y0),
            (x1, y1),
        );
        draw_ancestral_transitions(
            ctx,
            &track.tree,
            owner,
            &track.ancestral_state_layers,
            (x0, y0),
            (x1, y1),
        );
    }

    if track.show_nodes || track.support_style != SupportStyle::None {
        for node in &scene.visible {
            if scene.terminals.contains(node) {
                continue;
            }
            let Some(raw) = scene.positions[*node] else {
                continue;
            };
            let (x, y) = geometry.node(raw);
            if let Some(support) = track.tree.nodes()[*node].support.filter(|value| {
                track.support_style != SupportStyle::None
                    && support_fraction(*value)
                        .is_some_and(|value| value >= track.support_threshold)
            }) {
                draw_support(
                    ctx,
                    x,
                    y,
                    support,
                    &styles[*node].color,
                    track.support_style,
                );
            } else if track.show_nodes {
                ctx.svg.circle_ringed(
                    x,
                    y,
                    ctx.theme.tokens.marker_radius * 0.65,
                    &styles[*node].color,
                    &ctx.theme.background,
                    ctx.theme.tokens.hairline,
                );
            }
        }
    }

    draw_unrooted_node_glyphs(track, ctx, &scene, &geometry);

    draw_unrooted_trait_rings(track, ctx, &scene, &geometry);

    if track.show_tips {
        let size = ctx.theme.font_size - 1.0;
        for node in &scene.terminals {
            let angle = scene.angles[*node].unwrap_or(track.radial.start_degrees.to_radians());
            let (x, y) = geometry.point(geometry.label_radius, angle);
            let right = angle.cos() >= 0.0;
            ctx.svg.text_rotated(
                (x, y + size * 0.32),
                if right {
                    angle.to_degrees()
                } else {
                    angle.to_degrees() + 180.0
                },
                &terminal_label(&track.tree, *node, &track.collapsed),
                &ctx.theme.muted,
                size,
                if right {
                    crate::svg::Anchor::Start
                } else {
                    crate::svg::Anchor::End
                },
            );
        }
    }
    draw_trait_ring_headings(track, ctx);
    if let Some(bar) = track.branch_scale() {
        draw_unrooted_scale_bar(ctx, &scene, &geometry, area, bar);
    }
    draw_annotation_legend(track, ctx);
}

pub(super) fn unrooted_branch_colors(
    tree: &Tree,
    scene: &UnrootedScene,
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
    let numeric: Vec<f64> = scene
        .visible
        .iter()
        .filter_map(|node| values[*node].and_then(AnnotationValue::as_number))
        .collect();
    let present = scene
        .visible
        .iter()
        .filter(|node| values[**node].is_some())
        .count();
    if !numeric.is_empty() && numeric.len() == present {
        let minimum = numeric.iter().copied().fold(f64::MAX, f64::min);
        let maximum = numeric.iter().copied().fold(f64::MIN, f64::max);
        for node in &scene.visible {
            if let Some(value) = values[*node].and_then(AnnotationValue::as_number) {
                let fraction = if maximum <= minimum {
                    1.0
                } else {
                    (value - minimum) / (maximum - minimum)
                };
                colors[*node] = mix(&theme.muted, &theme.accent, fraction);
            }
        }
    } else {
        let mut categories = BTreeMap::new();
        for node in &scene.visible {
            let Some(value) = values[*node] else {
                continue;
            };
            let value = value.to_string();
            let next = categories.len();
            let index = *categories.entry(value).or_insert(next);
            colors[*node] = theme.color(index).to_string();
        }
    }
    colors
}

pub(super) fn draw_unrooted_trait_rings(
    track: &TreeTrack,
    ctx: &mut DrawContext<'_>,
    scene: &UnrootedScene,
    geometry: &UnrootedGeometry,
) {
    if track.trait_columns.is_empty() || scene.terminals.is_empty() {
        return;
    }
    let gap = ctx.theme.tokens.legend_gap.clamp(1.0, 4.0);
    let mut inner = geometry.branch_radius + gap;
    let step = std::f64::consts::TAU / scene.terminals.len() as f64;
    for column in &track.trait_columns {
        let outer = (inner + column.ring_width).min(geometry.ring_outer);
        let values: Vec<Option<&AnnotationValue>> = scene
            .terminals
            .iter()
            .map(|node| inherited_annotation(&track.tree, *node, &column.key))
            .collect();
        let domain = TraitDomain::new(
            scene
                .visible
                .iter()
                .filter_map(|node| inherited_annotation(&track.tree, *node, &column.key)),
        );
        for (row, node) in scene.terminals.iter().enumerate() {
            let angle = scene.angles[*node]
                .unwrap_or(track.radial.start_degrees.to_radians() + row as f64 * step);
            let gap_angle = if outer > 0.0 { 0.8 / outer } else { 0.0 };
            let half = (step / 2.0 - gap_angle).max(step * 0.12);
            let value = values[row];
            let name = terminal_label(&track.tree, *node, &track.collapsed);
            let title = match value {
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
                angle - half,
                angle + half,
                angle,
            );
        }
        inner = outer + gap;
    }
}
