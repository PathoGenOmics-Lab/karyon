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
    /// The box the drawing actually occupies, as `(min x, min y, max x, max
    /// y)`. It is not centred on the origin: the origin is whichever node the
    /// walk started from, and a tree hangs off it however its branches fall.
    pub(super) bounds: (f64, f64, f64, f64),
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

        // Every edge asks how many tips lie beyond it, and the answers all
        // come from this one pass.
        let components = ComponentTerminals::new(&adjacency, &terminal_set);

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
                    .map(|(next, _)| components.beyond(*next, *candidate))
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
                .map(|(next, length)| (*next, *length, components.beyond(*next, task.node).max(1)))
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
        let bounds = positions.iter().flatten().fold(
            (f64::MAX, f64::MAX, f64::MIN, f64::MIN),
            |(x0, y0, x1, y1), (x, y)| (x0.min(*x), y0.min(*y), x1.max(*x), y1.max(*y)),
        );
        let bounds = if bounds.0 > bounds.2 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            bounds
        };
        UnrootedScene {
            positions,
            parents,
            angles,
            terminals,
            visible,
            bounds,
        }
    }
}

/// How many terminals lie on each side of each edge, worked out once.
///
/// An unrooted layout asks, for every edge it walks, how many tips are out
/// that way, and it asked by walking out that way and counting: order n work
/// at every one of n nodes. Measured on a twenty thousand tip tree capped to
/// rows, it took 79 ms at 500 rows, 94 at 1000, 162 at 2000 and 596 at 4000,
/// which is the shape of a square. The other two projections stayed flat.
///
/// The adjacency is a tree, so the question has a closed answer. Root it
/// anywhere, count the terminals under every node once from the bottom up, and
/// then the tips beyond an edge are either that subtree or everything else.
pub(super) struct ComponentTerminals {
    /// Terminals under each node, in the arbitrary rooting used here.
    under: Vec<usize>,
    /// That rooting's parent for each node.
    parent: Vec<usize>,
    /// Every node's share of the whole, which is its own component's total
    /// and not the tree's: hiding a clade can leave the visible edges in
    /// several disconnected pieces, and a count taken across pieces would be
    /// wrong in a way nothing else would catch.
    whole: Vec<usize>,
}

impl ComponentTerminals {
    pub(super) fn new(adjacency: &[Vec<(usize, f64)>], terminals: &BTreeSet<usize>) -> Self {
        let n = adjacency.len();
        let mut under = vec![0usize; n];
        let mut parent = vec![usize::MAX; n];
        let mut whole = vec![0usize; n];
        let mut seen = vec![false; n];
        for root in 0..n {
            if seen[root] {
                continue;
            }
            // One depth-first walk fixes a rooting and an order for this
            // piece, then the counts come back up that order without walking
            // anything a second time.
            let mut order = Vec::new();
            let mut stack = vec![root];
            seen[root] = true;
            while let Some(node) = stack.pop() {
                order.push(node);
                for (next, _) in &adjacency[node] {
                    if !seen[*next] {
                        seen[*next] = true;
                        parent[*next] = node;
                        stack.push(*next);
                    }
                }
            }
            for node in order.iter().rev() {
                under[*node] += usize::from(terminals.contains(node));
                if parent[*node] != usize::MAX {
                    let up = parent[*node];
                    under[up] += under[*node];
                }
            }
            let total = under[root];
            for node in &order {
                whole[*node] = total;
            }
        }
        ComponentTerminals {
            under,
            parent,
            whole,
        }
    }

    /// The terminals reachable from `start` without going back through
    /// `blocked`.
    pub(super) fn beyond(&self, start: usize, blocked: usize) -> usize {
        if self.parent.get(start).copied() == Some(blocked) {
            self.under.get(start).copied().unwrap_or(0)
        } else {
            // `blocked` is below `start`, so what lies beyond is this whole
            // piece less that side of it.
            self.whole
                .get(start)
                .copied()
                .unwrap_or(0)
                .saturating_sub(self.under.get(blocked).copied().unwrap_or(0))
        }
    }
}

impl UnrootedScene {
    /// Half the longer side of the drawing, in tree units. It is what a scale
    /// bar sizes itself against, and it replaced the distance from the walk's
    /// starting node, which said more about where the walk began than about
    /// how big the tree is.
    pub(super) fn span(&self) -> f64 {
        let (x0, y0, x1, y1) = self.bounds;
        ((x1 - x0).max(y1 - y0) / 2.0).max(1e-9)
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct UnrootedGeometry {
    /// Where the layout's own origin lands on the page.
    pub(super) cx: f64,
    pub(super) cy: f64,
    pub(super) scale: f64,
    pub(super) branch_radius: f64,
    pub(super) ring_outer: f64,
    pub(super) label_radius: f64,
    /// The centre of the page area, which the drawing is centred on and the
    /// rings are drawn around.
    pub(super) mx: f64,
    pub(super) my: f64,
    /// Whether anything needs the tips gathered onto a shared circle.
    pub(super) ringed: bool,
}

impl UnrootedGeometry {
    pub(super) fn new(track: &TreeTrack, theme: &Theme, scene: &UnrootedScene, area: Rect) -> Self {
        let size = theme.font_size - 1.0;
        let name = |node: &usize| terminal_label(&track.tree, *node, &track.collapsed);
        let label_extent = if track.show_tips {
            scene
                .terminals
                .iter()
                .map(|node| text_width(&name(node), size))
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
        // Two ways to label an unrooted tree, and the tips decide which. A
        // name belongs at the end of its own branch, and on a small tree that
        // is how it reads. Past some number of tips the names start landing on
        // each other, and then gathering them onto a circle with a short
        // leader each is the only way to read any of them. Annotation rings
        // force the circle whatever the count, since they need every tip at
        // one radius to line up against.
        let mut ringed = !track.trait_columns.is_empty();
        let half = (area.w.min(area.h) / 2.0 - 4.0).max(2.0);
        let ring_outer = (half - label_extent).max(4.0);
        let branch_radius = (ring_outer - ring_room).max(2.0);
        let (mx, my) = (area.x + area.w / 2.0, area.y + area.h / 2.0);

        let (x0, y0, x1, y1) = scene.bounds;
        let (span_x, span_y) = ((x1 - x0).max(1e-9), (y1 - y0).max(1e-9));
        let ring = |ringed: bool| {
            // The tips have to land on the circle, so the fit is the box half
            // span against the radius, and the box centre goes to the middle
            // of the area rather than the origin, which is only whichever node
            // the walk happened to start from.
            let scale = branch_radius / (span_x / 2.0).max(span_y / 2.0);
            UnrootedGeometry {
                cx: mx - (x0 + x1) / 2.0 * scale,
                cy: my - (y0 + y1) / 2.0 * scale,
                scale,
                branch_radius,
                ring_outer,
                label_radius: ring_outer + 4.0,
                mx,
                my,
                ringed,
            }
        };
        if ringed {
            return ring(true);
        }

        // With the names at the tips, what has to fit is the drawing plus a
        // label hanging off each tip along its own branch. A label is a fixed
        // number of pixels whatever the scale, so the two cannot be fitted in
        // one step: pick a scale, measure the box the labels make, shrink, and
        // repeat until it stops moving.
        let room_x = (area.w - 8.0).max(8.0);
        let room_y = (area.h - 8.0).max(8.0);
        let labels: Vec<(f64, f64, f64, f64)> = if track.show_tips {
            scene
                .terminals
                .iter()
                .filter_map(|node| {
                    let (x, y) = scene.positions[*node]?;
                    let angle = scene.angles[*node].unwrap_or(0.0);
                    let reach = text_width(&name(node), size) + 6.0;
                    Some((x, y, angle.cos() * reach, angle.sin() * reach))
                })
                .collect()
        } else {
            Vec::new()
        };
        let extent = |scale: f64| {
            let mut box_ = (x0 * scale, y0 * scale, x1 * scale, y1 * scale);
            for (x, y, dx, dy) in &labels {
                let (px, py) = (x * scale + dx, y * scale + dy);
                box_ = (
                    box_.0.min(px),
                    box_.1.min(py),
                    box_.2.max(px),
                    box_.3.max(py),
                );
            }
            box_
        };
        let mut scale = (room_x / span_x).min(room_y / span_y);
        for _ in 0..12 {
            let (a, b, c, d) = extent(scale);
            let factor = (room_x / (c - a).max(1e-9)).min(room_y / (d - b).max(1e-9));
            scale *= factor;
            if (factor - 1.0).abs() < 1e-4 {
                break;
            }
        }
        // Does any name land on its neighbour? Labels radiate, so the pairs
        // that can touch are the ones next to each other going round, and each
        // label is checked against the next one along by sampling both down
        // their length. Anything closer than a line of text is a collision,
        // and one collision is enough: a figure with two names written over
        // each other is a figure with two names nobody can read.
        let (a, b, c, d) = extent(scale);
        /// One name to check against its neighbours: the angle it sits at
        /// going round, where it starts and where it ends.
        struct Reading {
            angle: f64,
            from: (f64, f64),
            to: (f64, f64),
        }
        let mut order: Vec<Reading> = labels
            .iter()
            .map(|(x, y, dx, dy)| {
                let (px, py) = (x * scale, y * scale);
                Reading {
                    angle: (py - (b + d) / 2.0).atan2(px - (a + c) / 2.0),
                    from: (px, py),
                    to: (px + dx, py + dy),
                }
            })
            .collect();
        order.sort_by(|left, right| left.angle.total_cmp(&right.angle));
        let sample = |from: (f64, f64), to: (f64, f64)| {
            [0.0, 0.25, 0.5, 0.75, 1.0]
                .map(|t| (from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t))
        };
        for pair in order.windows(2) {
            let (left, right) = (&pair[0], &pair[1]);
            let near = sample(left.from, left.to)
                .iter()
                .flat_map(|p| {
                    sample(right.from, right.to)
                        .into_iter()
                        .map(move |q| (p.0 - q.0).hypot(p.1 - q.1))
                })
                .fold(f64::MAX, f64::min);
            // Eight tenths of the body size, measured rather than picked.
            // Rendering forty five trees of eight to forty tips with the names
            // at the tips and intersecting the text boxes, every figure that
            // really had two names touching came in under 0.73 bodies and
            // every clean one at or above it, so this clears the dirty side
            // with a margin and still keeps twenty four of the forty five on
            // the readable rung, where a full body kept only nineteen.
            if near < size * 0.8 {
                ringed = true;
                break;
            }
        }
        if ringed {
            return ring(true);
        }

        UnrootedGeometry {
            // Centred on the box the names make, not the one the branches
            // make, so a long name on one side does not push the tree off the
            // other.
            cx: mx - (a + c) / 2.0,
            cy: my - (b + d) / 2.0,
            scale,
            branch_radius,
            ring_outer,
            label_radius: ring_outer + 4.0,
            mx,
            my,
            ringed,
        }
    }

    /// The direction a tip's name reads in, and for a ringed drawing the
    /// direction its leader runs.
    ///
    /// This is the angle the layout handed the tip, not the angle its drawn
    /// position happens to sit at. The equal-angle walk gives every tip its
    /// own slice of the circle, so those angles are spread evenly and the
    /// names come out evenly spaced on the ring, while the drawn positions
    /// bunch wherever the branches bunch. Using the drawn ones shortened every
    /// leader and made seventeen of seventeen test figures write names over
    /// each other.
    pub(super) fn outward(&self, scene: &UnrootedScene, node: usize) -> f64 {
        scene.angles[node].unwrap_or(0.0)
    }

    /// Where a tip's name starts.
    pub(super) fn tip_anchor(&self, scene: &UnrootedScene, node: usize, angle: f64) -> (f64, f64) {
        if self.ringed {
            return self.point(self.label_radius, angle);
        }
        let (x, y) = scene.positions[node].unwrap_or((0.0, 0.0));
        let (px, py) = self.node((x, y));
        (px + angle.cos() * 3.0, py + angle.sin() * 3.0)
    }

    pub(super) fn node(&self, point: (f64, f64)) -> (f64, f64) {
        (
            self.cx + point.0 * self.scale,
            self.cy + point.1 * self.scale,
        )
    }

    /// A point on a circle around the middle of the area, which is where the
    /// rings live. The drawing's own origin is somewhere else entirely.
    pub(super) fn point(&self, radius: f64, angle: f64) -> (f64, f64) {
        (
            self.mx + angle.cos() * radius,
            self.my + angle.sin() * radius,
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

    // Terminal leaders align the annotation rings without pretending unequal
    // branch lengths all end at the same evolutionary distance. With no rings
    // there is nothing to line up and no leader is drawn: the name sits at the
    // end of its own branch instead.
    if geometry.ringed {
        for node in &scene.terminals {
            let Some(raw) = scene.positions[*node] else {
                continue;
            };
            let angle = geometry.outward(&scene, *node);
            let (x0, y0) = geometry.node(raw);
            let (x1, y1) = geometry.point(geometry.branch_radius, angle);
            ctx.svg
                .line(x0, y0, x1, y1, &ctx.theme.rule, ctx.theme.tokens.hairline);
        }
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
        let style = styles.get(owner);
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
                    &styles.get(*node).color,
                    track.support_style,
                );
            } else if track.show_nodes {
                ctx.svg.circle_ringed(
                    x,
                    y,
                    ctx.theme.tokens.marker_radius * 0.65,
                    &styles.get(*node).color,
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
            let angle = geometry.outward(&scene, *node);
            let (x, y) = geometry.tip_anchor(&scene, *node, angle);
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
) -> PerNode<String> {
    let mut colors = PerNode::shared(default_color.to_string());
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
                colors.set(*node, mix(&theme.muted, &theme.accent, fraction));
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
            colors.set(*node, theme.color(index).to_string());
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
            .map(|node| row_annotation(&track.tree, *node, &column.key, &track.collapsed))
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
