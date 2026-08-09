//! The layout every projection reads, and the marks every projection draws.
//!
//! A [`TreeScene`] is the tree after collapsing, ladderising and rerooting have
//! been applied and before anything is drawn: which nodes are visible, where
//! each one sits, and what its terminal is called. The projections differ in
//! how they turn that into coordinates, not in what they turn.

use super::*;

pub(super) struct TreeScene {
    pub(super) placements: Vec<Option<Placement>>,
    pub(super) source_placements: Vec<Placement>,
    pub(super) terminals: Vec<usize>,
    pub(super) minimum: f64,
    pub(super) maximum: f64,
    pub(super) temporal: bool,
    pub(super) direction: TimeDirection,
}

impl TreeScene {
    pub(super) fn new(
        tree: &Tree,
        shape: TreeShape,
        time: Option<&TimeAxis>,
        collapsed: &BTreeSet<usize>,
    ) -> Self {
        let timed = time.and_then(|axis| tree.time_layout(&axis.key, axis.direction));
        let temporal = timed.is_some();
        let source_placements = timed.unwrap_or_else(|| tree.layout(shape == TreeShape::Cladogram));
        let direction = time.map_or(TimeDirection::Increasing, |axis| axis.direction);
        let terminals = visible_terminals(tree, collapsed);
        let mut rows = vec![None; tree.nodes().len()];
        for (row, node) in terminals.iter().enumerate() {
            rows[*node] = Some(row as f64);
        }
        let visible = visible_nodes(tree, collapsed);
        for node in postorder_nodes(tree) {
            if !visible[node] || rows[node].is_some() {
                continue;
            }
            let children: Vec<f64> = tree.nodes()[node]
                .children
                .iter()
                .filter(|child| visible[**child])
                .filter_map(|child| rows[*child])
                .collect();
            if !children.is_empty() {
                rows[node] = Some(children.iter().sum::<f64>() / children.len() as f64);
            }
        }
        let placements: Vec<Option<Placement>> = source_placements
            .iter()
            .map(|placement| {
                visible[placement.node].then_some(Placement {
                    row: rows[placement.node].unwrap_or(0.0),
                    ..*placement
                })
            })
            .collect();
        let (minimum, maximum) = source_placements
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
        TreeScene {
            placements,
            source_placements,
            terminals,
            minimum,
            maximum,
            temporal,
            direction,
        }
    }

    pub(super) fn x(&self, area: Rect, value: f64) -> f64 {
        area.x + self.fraction(value) * area.w
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
}

pub(super) fn visible_nodes(tree: &Tree, collapsed: &BTreeSet<usize>) -> Vec<bool> {
    let mut visible = vec![false; tree.nodes().len()];
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        visible[node] = true;
        if collapsed.contains(&node) {
            continue;
        }
        for child in tree.nodes()[node].children.iter().rev() {
            stack.push(*child);
        }
    }
    visible
}

pub(super) fn visible_terminals(tree: &Tree, collapsed: &BTreeSet<usize>) -> Vec<usize> {
    let mut terminals = Vec::new();
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        let clade = &tree.nodes()[node];
        if clade.is_leaf() || collapsed.contains(&node) {
            terminals.push(node);
            continue;
        }
        for child in clade.children.iter().rev() {
            stack.push(*child);
        }
    }
    terminals
}

pub(super) fn postorder_nodes(tree: &Tree) -> Vec<usize> {
    let mut order = Vec::with_capacity(tree.nodes().len());
    let mut stack = vec![tree.root()];
    while let Some(node) = stack.pop() {
        order.push(node);
        stack.extend(tree.nodes()[node].children.iter().copied());
    }
    order.reverse();
    order
}

pub(super) fn terminal_label(tree: &Tree, node: usize, collapsed: &BTreeSet<usize>) -> String {
    let name = tree.nodes()[node].name.as_deref().unwrap_or("clade");
    if collapsed.contains(&node) {
        format!("{} ({} tips)", name, tree.clade_size(node))
    } else {
        name.to_string()
    }
}

pub(super) fn support_fraction(value: f64) -> Option<f64> {
    if !value.is_finite() {
        return None;
    }
    let fraction = if value > 1.0 { value / 100.0 } else { value };
    Some(fraction.clamp(0.0, 1.0))
}

pub(super) fn draw_root_marker(ctx: &mut DrawContext<'_>, x: f64, y: f64) {
    ctx.svg.begin_titled("selected root");
    ctx.svg.symbol_ringed(
        x,
        y,
        ctx.theme.tokens.marker_radius * 1.15,
        crate::style::Symbol::Diamond,
        &ctx.theme.accent,
        &ctx.theme.background,
        ctx.theme.tokens.hairline.max(1.0),
    );
    ctx.svg.end_group();
}

pub(super) fn draw_support(
    ctx: &mut DrawContext<'_>,
    x: f64,
    y: f64,
    support: f64,
    color: &str,
    style: SupportStyle,
) {
    let fraction = support_fraction(support).unwrap_or(0.0);
    let radius = ctx.theme.tokens.marker_radius * (0.45 + fraction * 0.55);
    if style.symbols() {
        ctx.svg.circle_ringed(
            x,
            y,
            radius,
            color,
            &ctx.theme.background,
            ctx.theme.tokens.hairline,
        );
    }
    if style.labels() {
        let label = text_rounded(support, 3);
        let size = (ctx.theme.font_size - 3.0).max(6.0);
        let offset = if style.symbols() { radius + 2.5 } else { 3.0 };
        let width = text_width(&label, size) + 4.0;
        ctx.svg.rect_rounded(
            x + offset - 2.0,
            y - size - 1.0,
            width,
            size + 3.0,
            2.0,
            &ctx.theme.background,
        );
        ctx.svg.text(
            x + offset,
            y - 1.5,
            &label,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Start,
        );
    }
}

pub(super) fn draw_branch_annotation(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    labels: &BranchLabels,
    start: (f64, f64),
    end: (f64, f64),
) {
    let Some(value) = tree.annotation(node, &labels.key) else {
        return;
    };
    let exact = value.to_string();
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    let visible = fit_text(&exact, (length - 6.0).max(0.0), labels.size);
    if visible.is_empty() {
        return;
    }
    let angle = dy.atan2(dx);
    let mut rotation = angle.to_degrees().rem_euclid(360.0);
    if rotation > 90.0 && rotation < 270.0 {
        rotation += 180.0;
    }
    let offset = labels.size * 0.55 + 1.0;
    let x = (start.0 + end.0) / 2.0 + angle.sin() * offset;
    let y = (start.1 + end.1) / 2.0 - angle.cos() * offset;
    if visible != exact {
        ctx.svg.begin_titled(&format!("{} {exact}", labels.key));
    }
    ctx.svg.text_rotated(
        (x, y),
        rotation,
        &visible,
        &ctx.theme.muted,
        labels.size,
        crate::svg::Anchor::Middle,
    );
    if visible != exact {
        ctx.svg.end_group();
    }
}
