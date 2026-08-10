//! Branch-local rate mixtures and cross-branch recurrent-event links.

use super::*;

type Point = (f64, f64);
type BranchPoint = (usize, Point);

#[derive(Debug, Clone, Copy)]
pub(super) enum LinkGeometry {
    Rectangular { right: f64 },
    Centred { centre: (f64, f64) },
}

pub(super) fn draw_branch_rate_mixtures(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    mixtures: &[BranchRateMixture],
    start: (f64, f64),
    end: (f64, f64),
) {
    if mixtures.is_empty() {
        return;
    }
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let branch_length = dx.hypot(dy);
    if branch_length < 4.0 || !branch_length.is_finite() {
        return;
    }
    let direction = (dx / branch_length, dy / branch_length);
    let normal = (-direction.1, direction.0);
    let middle = ((start.0 + end.0) / 2.0, (start.1 + end.1) / 2.0);

    for (layer_index, mixture) in mixtures.iter().enumerate() {
        let Some(classes) = mixture_values(tree, node, mixture) else {
            continue;
        };
        let visual_width = mixture.width.min((branch_length * 0.78).max(4.0));
        let offset = layer_index as f64 * (mixture.thickness + 2.0);
        let centre = (middle.0 + normal.0 * offset, middle.1 + normal.1 * offset);
        let from = (
            centre.0 - direction.0 * visual_width / 2.0,
            centre.1 - direction.1 * visual_width / 2.0,
        );
        let to = (
            centre.0 + direction.0 * visual_width / 2.0,
            centre.1 + direction.1 * visual_width / 2.0,
        );
        let title = mixture_title(mixture, &classes);
        ctx.svg.begin_titled(&title);
        ctx.svg.line(
            from.0,
            from.1,
            to.0,
            to.1,
            ctx.theme.surface(),
            mixture.thickness + 2.2,
        );
        let mut cursor = 0.0;
        for class in &classes {
            let next = (cursor + visual_width * class.normalised_weight).min(visual_width);
            if next > cursor {
                let segment_start = (from.0 + direction.0 * cursor, from.1 + direction.1 * cursor);
                let segment_end = (from.0 + direction.0 * next, from.1 + direction.1 * next);
                ctx.svg.line(
                    segment_start.0,
                    segment_start.1,
                    segment_end.0,
                    segment_end.1,
                    &omega_color(
                        ctx.theme,
                        class.rate,
                        mixture.neutral_lower,
                        mixture.neutral_upper,
                        mixture.saturation,
                    ),
                    mixture.thickness,
                );
            }
            cursor = next;
        }
        ctx.svg.circle(
            from.0,
            from.1,
            mixture.thickness * 0.50,
            &omega_color(
                ctx.theme,
                classes[0].rate,
                mixture.neutral_lower,
                mixture.neutral_upper,
                mixture.saturation,
            ),
        );
        if let Some(last) = classes.last() {
            ctx.svg.circle(
                to.0,
                to.1,
                mixture.thickness * 0.50,
                &omega_color(
                    ctx.theme,
                    last.rate,
                    mixture.neutral_lower,
                    mixture.neutral_upper,
                    mixture.saturation,
                ),
            );
        }
        ctx.svg.end_group();
    }
}

#[derive(Debug, Clone, Copy)]
struct MixtureClass {
    rate: f64,
    source_weight: f64,
    normalised_weight: f64,
}

fn mixture_values(
    tree: &Tree,
    node: usize,
    mixture: &BranchRateMixture,
) -> Option<Vec<MixtureClass>> {
    if mixture.rate_keys.len() != mixture.weight_keys.len() {
        return None;
    }
    let raw: Vec<(f64, f64)> = mixture
        .rate_keys
        .iter()
        .zip(&mixture.weight_keys)
        .filter_map(|(rate_key, weight_key)| {
            let rate = tree.annotation(node, rate_key)?.as_number()?;
            let weight = tree.annotation(node, weight_key)?.as_number()?;
            (rate.is_finite() && rate >= 0.0 && weight.is_finite() && weight > 0.0)
                .then_some((rate, weight))
        })
        .collect();
    let total: f64 = raw.iter().map(|(_, weight)| *weight).sum();
    if raw.is_empty() || !total.is_finite() || total <= 0.0 {
        return None;
    }
    Some(
        raw.into_iter()
            .map(|(rate, source_weight)| MixtureClass {
                rate,
                source_weight,
                normalised_weight: source_weight / total,
            })
            .collect(),
    )
}

fn mixture_title(mixture: &BranchRateMixture, classes: &[MixtureClass]) -> String {
    let details = classes
        .iter()
        .enumerate()
        .map(|(index, class)| {
            format!(
                "class {} omega {} weight {}",
                index + 1,
                text_rounded(class.rate, 4),
                text_rounded(class.source_weight, 4)
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("{} | {details}", mixture.label)
}

pub(super) fn draw_homoplasy_links(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    layers: &[HomoplasyLayer],
    points: &[BranchPoint],
    geometry: LinkGeometry,
) {
    for (layer_index, layer) in layers.iter().enumerate() {
        let mut groups: BTreeMap<String, Vec<BranchPoint>> = BTreeMap::new();
        for (node, point) in points {
            let Some(value) = tree.annotation(*node, &layer.key) else {
                continue;
            };
            let label = value.to_string();
            if label.is_empty() {
                continue;
            }
            groups.entry(label).or_default().push((*node, *point));
        }

        let mut emitted = 0usize;
        for (event_index, (event, mut occurrences)) in groups.into_iter().enumerate() {
            if occurrences.len() < layer.minimum_occurrences {
                continue;
            }
            match geometry {
                LinkGeometry::Rectangular { .. } => {
                    occurrences.sort_by(|left, right| left.1 .1.total_cmp(&right.1 .1));
                }
                LinkGeometry::Centred { centre } => occurrences.sort_by(|left, right| {
                    let left_angle = (left.1 .1 - centre.1).atan2(left.1 .0 - centre.0);
                    let right_angle = (right.1 .1 - centre.1).atan2(right.1 .0 - centre.0);
                    left_angle.total_cmp(&right_angle)
                }),
            }
            let count = occurrences.len();
            let color = mix(
                ctx.theme.surface(),
                ctx.theme.color(layer_index + event_index),
                0.72,
            );
            let title = format!(
                "recurrent event {} = {}; {} branches",
                layer.key, event, count
            );
            for pair in occurrences.windows(2) {
                if emitted >= layer.maximum_connections {
                    break;
                }
                let from = pair[0].1;
                let to = pair[1].1;
                let path = link_path(from, to, geometry, emitted);
                ctx.svg.begin_titled(&title);
                ctx.svg
                    .path_stroked_pattern(&path, &color, layer.width, LinePattern::Dashed);
                ctx.svg
                    .circle_ringed(from.0, from.1, 2.2, &color, ctx.theme.surface(), 0.8);
                ctx.svg
                    .circle_ringed(to.0, to.1, 2.2, &color, ctx.theme.surface(), 0.8);
                ctx.svg.end_group();
                emitted += 1;
            }
            if emitted >= layer.maximum_connections {
                break;
            }
        }
    }
}

fn link_path(from: (f64, f64), to: (f64, f64), geometry: LinkGeometry, index: usize) -> String {
    match geometry {
        LinkGeometry::Rectangular { right } => {
            let bend = (from.0.max(to.0) + 12.0 + index as f64 * 2.0).min(right - 1.0);
            format!(
                "M {} {} C {} {} {} {} {} {}",
                num(from.0),
                num(from.1),
                num(bend),
                num(from.1),
                num(bend),
                num(to.1),
                num(to.0),
                num(to.1)
            )
        }
        LinkGeometry::Centred { centre } => {
            let pull = 0.72;
            let c1 = (
                from.0 + (centre.0 - from.0) * pull,
                from.1 + (centre.1 - from.1) * pull,
            );
            let c2 = (
                to.0 + (centre.0 - to.0) * pull,
                to.1 + (centre.1 - to.1) * pull,
            );
            format!(
                "M {} {} C {} {} {} {} {} {}",
                num(from.0),
                num(from.1),
                num(c1.0),
                num(c1.1),
                num(c2.0),
                num(c2.1),
                num(to.0),
                num(to.1)
            )
        }
    }
}
