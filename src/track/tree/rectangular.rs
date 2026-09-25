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
    color_levels: Dealt<'_>,
    dnds: Option<&DnDsLayer>,
    show_nodes: bool,
    support_style: SupportStyle,
    support: SupportReading,
    branch_labels: Option<&BranchLabels>,
    rate_mixtures: &[BranchRateMixture],
    homoplasy_layers: &[HomoplasyLayer],
    branch_event_layers: &[BranchEventLayer],
    branch_interval_layers: &[BranchIntervalLayer],
    ancestral_state_layers: &[AncestralStateLayer],
    branch_geometry: BranchGeometry,
    name_leaves: bool,
) {
    let colors = branch_colors(
        tree,
        scene,
        color_by,
        color_levels,
        ctx.theme,
        default_color,
    );
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
        let style = styles.get(placement.node);
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
            let connector = connector_style(dnds, ctx.theme, colors.get(placement.node), width);
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
        if let Some(value) = node
            .support
            .filter(|value| support_style != SupportStyle::None && support.shown(*value))
        {
            let title = format!("clade support {}", text_rounded(value, 3));
            ctx.svg.begin_titled(&title);
            draw_support(
                ctx,
                x,
                y_of(placement.row),
                value,
                support,
                &styles.get(placement.node).color,
                support_style,
            );
            ctx.svg.end_group();
        } else if show_nodes {
            ctx.svg.circle_ringed(
                x,
                y_of(placement.row),
                ctx.theme.tokens.marker_radius * 0.65,
                &styles.get(placement.node).color,
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
            let title = collapsed_title(tree, *node);
            ctx.svg.begin_titled(&title);
            ctx.svg.path(&d, &styles.get(*node).color, 0.28);
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

/// One value for every node, where almost every node has the same one.
///
/// A million tip tree is two million nodes, and giving each its own `String`
/// is two million allocations to say the same colour two million times, which
/// was 130 of the 190 ms such a tree spent drawing sixty rows. The nodes that
/// differ are the ones the figure was asked to colour, and a figure that
/// colours every node of a million tip tree has nothing to say either way.
#[derive(Debug, Clone)]
pub(super) struct PerNode<T> {
    common: T,
    apart: BTreeMap<usize, T>,
}

impl<T> PerNode<T> {
    pub(super) fn shared(common: T) -> Self {
        PerNode {
            common,
            apart: BTreeMap::new(),
        }
    }

    pub(super) fn set(&mut self, node: usize, value: T) {
        self.apart.insert(node, value);
    }

    pub(super) fn get(&self, node: usize) -> &T {
        self.apart.get(&node).unwrap_or(&self.common)
    }
}

/// The styles branches are drawn with, worked out as each branch is reached.
///
/// Held rather than built because building it meant one `BranchStyle` per node
/// whether or not the node is drawn, and a folded tree draws a few hundred of
/// its two million.
pub(super) struct BranchStyles<'a> {
    tree: &'a Tree,
    colors: &'a PerNode<String>,
    dnds: Option<&'a DnDsLayer>,
    theme: &'a Theme,
    width: f64,
}

impl BranchStyles<'_> {
    pub(super) fn get(&self, node: usize) -> BranchStyle {
        match self.dnds {
            Some(dnds) => dnds_branch_style(self.tree, node, dnds, self.theme, self.width),
            None => BranchStyle {
                color: self.colors.get(node).clone(),
                width: self.width,
                pattern: LinePattern::Solid,
            },
        }
    }
}

pub(super) fn branch_styles<'a>(
    tree: &'a Tree,
    colors: &'a PerNode<String>,
    dnds: Option<&'a DnDsLayer>,
    theme: &'a Theme,
    width: f64,
) -> BranchStyles<'a> {
    BranchStyles {
        tree,
        colors,
        dnds,
        theme,
        width,
    }
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
    levels: Dealt<'_>,
    theme: &Theme,
    default_color: &str,
) -> PerNode<String> {
    let mut colors = PerNode::shared(default_color.to_string());
    let Some(key) = key else {
        return colors;
    };
    let values = branch_values(tree, key);
    let domain = TraitDomain::ordered(
        levels.first,
        levels.levels,
        values.iter().flatten().copied(),
    );
    let continuous = is_continuous(&values);
    for node in scene
        .placements
        .iter()
        .flatten()
        .map(|placement| placement.node)
    {
        let color = if continuous {
            domain
                .fraction(values[node])
                .map(|fraction| mix(&theme.muted, &theme.accent, fraction))
        } else {
            domain
                .category(values[node])
                .map(|index| theme.color(index).to_string())
        };
        if let Some(color) = color {
            colors.set(node, color);
        }
    }
    colors
}

/// Whether `--color-by` reads an annotation as a ramp rather than as levels:
/// when every value it has anywhere in the tree is a number.
///
/// Over the whole tree, as the domain is, so folding a clade of the only
/// words cannot turn the rest of the tree into a ramp.
pub(super) fn is_continuous(values: &[Option<&AnnotationValue>]) -> bool {
    let mut stated = values.iter().flatten().peekable();
    stated.peek().is_some() && stated.all(|value| value.as_number().is_some())
}

/// The levels and the range one annotation covers, over the whole tree.
///
/// Every colour the track paints for `key` is read off this one count: the
/// branches `--color-by` colours, the strips beside the tips in all three
/// projections, [`TreeTrack::strips`] and [`TreeTrack::legend`]. It used to be
/// worked out from the nodes on screen, once per picture, so the key a caller
/// built from the sample sheet named the levels in the order the sheet sorts
/// its rows while the tree numbered them in the order it met them, and a
/// figure of two countries printed each one's colour beside the other's name.
/// Counting from what is on screen also meant that folding a clade could
/// repaint the rest of the tree. Counting over every node, in the tree's own
/// order, gives one answer per tree and key whatever is folded or scrolled
/// away, and the values folded rows show are their tips' values, so none is
/// left without a colour.
///
/// `levels` is the order a column carries from its sample sheet, dealt the
/// palette first, so a level is the colour here that it is in every other strip
/// the sheet is drawn in. Empty, the tree deals the palette in its own order.
pub(super) fn tree_domain(tree: &Tree, key: &str, levels: Dealt<'_>) -> TraitDomain {
    TraitDomain::ordered(
        levels.first,
        levels.levels,
        branch_values(tree, key).into_iter().flatten(),
    )
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

/// The value every node takes its colour from, worked out in one pass.
///
/// A node's own annotation first, then one inherited from an ancestor, which is
/// how a tree out of BEAST carries a state down a lineage. Failing both, what
/// its descendants agree on, which is how a sheet keyed by sample name reaches
/// the branches above the tips: without it a tree coloured by lineage had its
/// two hundred terminal branches in six colours and all three hundred and
/// ninety seven internal ones left black, which says the clades are of unknown
/// lineage when every tip in them says otherwise.
///
/// One post-order pass and not a walk per node. Asking each node separately
/// what its descendants agree on is a traversal per node, which is the shape
/// this crate has spent a good deal of effort removing.
pub(super) fn branch_values<'a>(tree: &'a Tree, key: &str) -> Vec<Option<&'a AnnotationValue>> {
    let nodes = tree.nodes();
    let mut agreed: Vec<Option<&AnnotationValue>> = vec![None; nodes.len()];
    let mut settled = vec![false; nodes.len()];
    for node in postorder_nodes(tree) {
        if nodes[node].is_leaf() {
            agreed[node] = tree.annotation(node, key);
            settled[node] = agreed[node].is_some();
            continue;
        }
        let mut value: Option<&AnnotationValue> = None;
        let mut uniform = !nodes[node].children.is_empty();
        for child in &nodes[node].children {
            if !settled[*child] {
                uniform = false;
                break;
            }
            match value {
                None => value = agreed[*child],
                Some(had) if Some(had) == agreed[*child] => {}
                Some(_) => {
                    uniform = false;
                    break;
                }
            }
        }
        agreed[node] = if uniform { value } else { None };
        settled[node] = uniform && value.is_some();
    }
    (0..nodes.len())
        .map(|node| inherited_annotation(tree, node, key).or(agreed[node]))
        .collect()
}

/// What a folded clade can honestly show in a metadata strip.
///
/// A folded row is an internal node, and an internal node carries no sample's
/// metadata, so reading it the ordinary way gives nothing: a tree drawn to
/// sixty rows had every cell of every strip empty, which is what a figure looks
/// like when it read the wrong file. What it can say is what its tips agree on.
/// A clade whose samples are all L4 is an L4 clade and the strip says so; a
/// clade holding two lineages is not either of them and the strip stays empty
/// rather than picking one.
///
/// One tip with nothing recorded is enough to withhold it. A clade cannot be
/// called uniform on the strength of the members that happen to have been
/// typed, and a strip that quietly ignored the untyped ones would be claiming
/// more than the sheet says.
pub(super) fn agreed_annotation<'a>(
    tree: &'a Tree,
    node: usize,
    key: &str,
) -> Option<&'a AnnotationValue> {
    let mut agreed: Option<&AnnotationValue> = None;
    let mut stack = vec![node];
    while let Some(at) = stack.pop() {
        let clade = &tree.nodes()[at];
        if clade.is_leaf() {
            let value = inherited_annotation(tree, at, key)?;
            match agreed {
                None => agreed = Some(value),
                Some(had) if had == value => {}
                Some(_) => return None,
            }
            continue;
        }
        for child in &clade.children {
            stack.push(*child);
        }
    }
    agreed
}

/// The value a drawn row stands for, folded or not.
pub(super) fn row_annotation<'a>(
    tree: &'a Tree,
    node: usize,
    key: &str,
    collapsed: &BTreeSet<usize>,
) -> Option<&'a AnnotationValue> {
    if collapsed.contains(&node) {
        agreed_annotation(tree, node, key)
    } else {
        inherited_annotation(tree, node, key)
    }
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
        // Over the whole tree and not the rows on screen, for the reasons
        // `tree_domain` gives. A folded row shows what its tips agree on, and
        // that is one of their values, so it has a colour here: counting only
        // the nodes the walk placed once left forty rows of lineage as forty
        // empty outlines.
        let domain = tree_domain(tree, &column.key, column.dealt());
        let rows: Vec<TraitRow<'_>> = names
            .iter()
            .zip(&scene.terminals)
            .enumerate()
            .map(|(row, (name, node))| TraitRow {
                name,
                top: area.y + row as f64 * row_pitch + 1.0,
                height: (row_pitch - 2.0).max(1.0),
                value: row_annotation(tree, *node, &column.key, collapsed),
            })
            .collect();
        draw_column(ctx, column, &domain, x, Some(area.y - 5.0), &rows);
        x += column.width + ctx.theme.tokens.legend_gap;
    }
}
