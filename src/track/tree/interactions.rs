//! Branch-local rate mixtures and cross-branch recurrent-event links.

use super::*;
use crate::svg::text_exact;

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

/// Names every class the capsule drew, with the two numbers it was drawn from.
///
/// Both numbers are written exactly, and the two lines below are one decision
/// rather than two, because the capsule keeps neither of them. It normalises the
/// weights for its geometry, so weights of 0.00001 and 0.00001 paint the picture
/// that 0.5 and 0.5 paint and this string is the only place a mixture that does
/// not sum to one can be seen at all; it paints the rate as a colour on a band
/// around neutral, so this string is likewise the only place a rate of 0.00001
/// differs from a rate of nought. Rounded to four decimals, both of those
/// printed `0`: a class described by the one number `mixture_values` refuses to
/// draw, under a capsule split precisely in half. `text_exact` is the crate's
/// answer to a value whose precision belongs to the caller, and the type's own
/// documentation already promised the supplied weights survive into the tooltip.
fn mixture_title(mixture: &BranchRateMixture, classes: &[MixtureClass]) -> String {
    let details = classes
        .iter()
        .enumerate()
        .map(|(index, class)| {
            format!(
                "class {} omega {} weight {}",
                index + 1,
                text_exact(class.rate),
                text_exact(class.source_weight)
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("{} | {details}", mixture.label)
}

pub(super) fn draw_branch_event_layers(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    layers: &[BranchEventLayer],
    start: Point,
    end: Point,
) {
    let Some((direction, normal, branch_length)) = branch_vectors(start, end) else {
        return;
    };
    if branch_length < 6.0 {
        return;
    }
    for (layer_index, layer) in layers.iter().enumerate() {
        let Some(annotation) = tree.annotation(node, &layer.key) else {
            continue;
        };
        let events = event_values(annotation);
        if events.is_empty() {
            continue;
        }
        let visible = events.len().min(layer.maximum_events);
        let hidden = events.len() - visible;
        let offset = -(layer.size + 4.0 + layer_index as f64 * (layer.size * 2.0 + 2.0));
        for (index, event) in events.iter().take(visible).enumerate() {
            let fraction = (index + 1) as f64 / (visible + 1) as f64;
            let at = (
                start.0 + direction.0 * branch_length * fraction + normal.0 * offset,
                start.1 + direction.1 * branch_length * fraction + normal.1 * offset,
            );
            let category = event_category(layer, event);
            let symbol = match category % 4 {
                0 => crate::style::Symbol::Diamond,
                1 => crate::style::Symbol::Circle,
                2 => crate::style::Symbol::Square,
                _ => crate::style::Symbol::Triangle,
            };
            let mut title = format!("{} | {} = {}", layer.label, layer.key, event);
            if hidden > 0 && index + 1 == visible {
                title.push_str(&format!(" | {hidden} additional events not drawn"));
            }
            ctx.svg.begin_titled(&title);
            ctx.svg.symbol_ringed(
                at.0,
                at.1,
                layer.size,
                symbol,
                ctx.theme.color(category),
                ctx.theme.surface(),
                0.9,
            );
            ctx.svg.end_group();
        }
    }
}

fn event_values(value: &AnnotationValue) -> Vec<String> {
    match value {
        AnnotationValue::List(values) => values
            .iter()
            .map(ToString::to_string)
            .filter(|value| !value.is_empty())
            .collect(),
        _ => {
            let value = value.to_string();
            (!value.is_empty()).then_some(value).into_iter().collect()
        }
    }
}

fn event_category(layer: &BranchEventLayer, event: &str) -> usize {
    // An explicit hash keeps event identity stable when taxa or other events
    // are appended, and avoids rescanning the whole tree for every mark.
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in layer.key.bytes().chain([0]).chain(event.bytes()) {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash as usize
}

pub(super) fn draw_branch_intervals(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    layers: &[BranchIntervalLayer],
    start: Point,
    end: Point,
) {
    let Some((direction, normal, branch_length)) = branch_vectors(start, end) else {
        return;
    };
    for (index, layer) in layers.iter().enumerate() {
        let Some(estimate) = direct_number(tree, node, &layer.estimate_key) else {
            continue;
        };
        let Some(lower) = direct_number(tree, node, &layer.lower_key) else {
            continue;
        };
        let Some(upper) = direct_number(tree, node, &layer.upper_key) else {
            continue;
        };
        if upper < lower {
            continue;
        }
        let width = layer.width.min(branch_length * 0.82);
        if width < 7.0 {
            continue;
        }
        let offset = 7.0 + index as f64 * 7.0;
        let centre = (
            (start.0 + end.0) / 2.0 + normal.0 * offset,
            (start.1 + end.1) / 2.0 + normal.1 * offset,
        );
        let origin = (
            centre.0 - direction.0 * width / 2.0,
            centre.1 - direction.1 * width / 2.0,
        );
        let point = |value: f64| {
            let fraction = interval_fraction(value, layer.minimum, layer.maximum);
            (
                origin.0 + direction.0 * width * fraction,
                origin.1 + direction.1 * width * fraction,
            )
        };
        let lower_at = point(lower);
        let upper_at = point(upper);
        let estimate_at = point(estimate);
        let highlighted = layer
            .threshold
            .is_some_and(|threshold| estimate >= threshold);
        let color = if highlighted {
            ctx.theme.color(1)
        } else {
            ctx.theme.color(0)
        };
        // Exactly, not to five decimals. The tooltip is the last surviving copy
        // of the number the caller handed over: the drawing keeps only where it
        // fell on the axis. An estimate of 2.2e-308 rounded to five decimals is
        // the string "0", and then the mark sitting correctly at the top of its
        // own domain reads as a nought drawn at the ceiling. The same mistake
        // was found in the selection tooltips, and this is the same answer.
        let title = format!(
            "{} | estimate {} | interval {} to {}{}",
            layer.label,
            text_exact(estimate),
            text_exact(lower),
            text_exact(upper),
            layer
                .threshold
                .map_or_else(String::new, |threshold| format!(
                    " | threshold {} ({})",
                    text_exact(threshold),
                    if highlighted {
                        "reached"
                    } else {
                        "not reached"
                    }
                ))
        );
        ctx.svg.begin_titled(&title);
        ctx.svg.line(
            origin.0,
            origin.1,
            origin.0 + direction.0 * width,
            origin.1 + direction.1 * width,
            ctx.theme.surface(),
            5.2,
        );
        ctx.svg.line(
            origin.0,
            origin.1,
            origin.0 + direction.0 * width,
            origin.1 + direction.1 * width,
            &ctx.theme.rule,
            1.0,
        );
        ctx.svg
            .line(lower_at.0, lower_at.1, upper_at.0, upper_at.1, color, 2.6);
        for at in [lower_at, upper_at] {
            ctx.svg.line(
                at.0 - normal.0 * 2.4,
                at.1 - normal.1 * 2.4,
                at.0 + normal.0 * 2.4,
                at.1 + normal.1 * 2.4,
                color,
                1.0,
            );
        }
        ctx.svg.circle_ringed(
            estimate_at.0,
            estimate_at.1,
            2.5,
            color,
            ctx.theme.surface(),
            0.8,
        );
        ctx.svg.end_group();
    }
}

/// How far along the layer's domain `value` sits, as a fraction of the way.
///
/// The halving is the whole function, and it has a name of its own because the
/// width of a domain is the same arithmetic `mixture_values` and `maximum_state`
/// already guard a few lines apart: a total that can leave the range a fraction
/// can be taken from. `BranchIntervalLayer::range` takes any finite pair, so
/// `range(f64::MIN, f64::MAX)` is what somebody writes when they mean "any
/// number at all", and `maximum - minimum` is then infinity. Divided by that, an
/// ordinary estimate came back nought and was drawn on the left end of its own
/// axis while the tooltip beside it still read the number it was drawn from, and
/// an estimate near the end of the domain came back infinity over infinity,
/// which is NaN, which the clamp keeps and the writer then refuses: a caption
/// over the blank paper where the mark should be. Halving both ends first cannot
/// overflow for any finite pair, and it runs only when the plain span was not a
/// number an axis can divide by, so no mark that was already placed correctly
/// moves by so much as an ulp.
///
/// The numerator needs no such care. It overflows only when the value is
/// further from the domain than the domain is wide, which is a value off the
/// end of the axis, and the clamp puts it on the end it ran off.
fn interval_fraction(value: f64, minimum: f64, maximum: f64) -> f64 {
    let span = maximum - minimum;
    let fraction = if span.is_finite() {
        (value - minimum) / span
    } else {
        (value * 0.5 - minimum * 0.5) / (maximum * 0.5 - minimum * 0.5)
    };
    fraction.clamp(0.0, 1.0)
}

pub(super) fn draw_ancestral_transitions(
    ctx: &mut DrawContext<'_>,
    tree: &Tree,
    node: usize,
    layers: &[AncestralStateLayer],
    start: Point,
    end: Point,
) {
    let Some(parent) = tree.nodes()[node].parent else {
        return;
    };
    let Some((direction, normal, _)) = branch_vectors(start, end) else {
        return;
    };
    for (layer_index, layer) in layers.iter().enumerate() {
        if !layer.show_transitions {
            continue;
        }
        let Some((parent_state, parent_probability)) = maximum_state(tree, parent, layer) else {
            continue;
        };
        let Some((child_state, child_probability)) = maximum_state(tree, node, layer) else {
            continue;
        };
        if parent_state == child_state
            || parent_probability < layer.confidence
            || child_probability < layer.confidence
        {
            continue;
        }
        let offset = -(14.0 + layer_index as f64 * 7.0);
        let middle = (
            (start.0 + end.0) / 2.0 + normal.0 * offset,
            (start.1 + end.1) / 2.0 + normal.1 * offset,
        );
        let from = (middle.0 - direction.0 * 3.0, middle.1 - direction.1 * 3.0);
        let to = (middle.0 + direction.0 * 3.0, middle.1 + direction.1 * 3.0);
        let title = format!(
            "{} transition {} ({}) to {} ({})",
            layer.label,
            layer.keys[parent_state],
            text_rounded(parent_probability, 4),
            layer.keys[child_state],
            text_rounded(child_probability, 4)
        );
        ctx.svg.begin_titled(&title);
        ctx.svg
            .line(from.0, from.1, to.0, to.1, ctx.theme.surface(), 5.2);
        ctx.svg
            .line(from.0, from.1, to.0, to.1, &ctx.theme.muted, 1.0);
        ctx.svg.circle_ringed(
            from.0,
            from.1,
            2.2,
            ctx.theme.color(parent_state),
            ctx.theme.surface(),
            0.7,
        );
        ctx.svg.symbol_ringed(
            to.0,
            to.1,
            2.8,
            crate::style::Symbol::Diamond,
            ctx.theme.color(child_state),
            ctx.theme.surface(),
            0.7,
        );
        ctx.svg.end_group();
    }
}

fn maximum_state(tree: &Tree, node: usize, layer: &AncestralStateLayer) -> Option<(usize, f64)> {
    let values: Vec<f64> = layer
        .keys
        .iter()
        .map(|key| direct_number(tree, node, key).filter(|value| *value >= 0.0))
        .collect::<Option<Vec<_>>>()?;
    let total: f64 = values.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return None;
    }
    values
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1))
        .map(|(index, value)| (index, *value / total))
}

fn direct_number(tree: &Tree, node: usize, key: &str) -> Option<f64> {
    tree.annotation(node, key)
        .and_then(AnnotationValue::as_number)
        .filter(|value| value.is_finite())
}

fn branch_vectors(start: Point, end: Point) -> Option<(Point, Point, f64)> {
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if !length.is_finite() || length <= f64::EPSILON {
        return None;
    }
    let direction = (dx / length, dy / length);
    Some((direction, (-direction.1, direction.0), length))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("tree", 0, 1_000).unwrap()
    }

    /// The value of `name` in the fragment, which the drawing wrote itself.
    fn attribute(fragment: &str, name: &str) -> f64 {
        let at = fragment
            .find(name)
            .unwrap_or_else(|| panic!("no {name} in {fragment}"));
        let rest = &fragment[at + name.len()..];
        let end = rest.find('"').unwrap_or_else(|| panic!("{fragment}"));
        rest[..end].parse().unwrap_or_else(|_| panic!("{fragment}"))
    }

    /// The one titled group whose tooltip starts with `opening`.
    fn group<'a>(svg: &'a str, opening: &str) -> &'a str {
        let at = svg
            .find(opening)
            .unwrap_or_else(|| panic!("no group opening {opening:?} in {svg}"));
        let rest = &svg[at..];
        &rest[..rest.find("</g>").unwrap_or_else(|| panic!("{rest}"))]
    }

    fn interval(source: &str, minimum: f64, maximum: f64) -> String {
        Figure::new(region())
            .width(700.0)
            .show_region_label(false)
            .push(
                TreeTrack::new(Tree::parse_annotated_newick(source).unwrap()).branch_interval(
                    BranchIntervalLayer::new("cf", "cf_lo", "cf_hi").range(minimum, maximum),
                ),
            )
            .to_svg()
    }

    #[test]
    fn an_interval_is_placed_on_a_domain_wider_than_f64_can_subtract() {
        // `range(f64::MIN, f64::MAX)` is what a caller writes for "any number at
        // all", and its width is the one span the subtraction cannot return. The
        // estimate used to be divided by that infinity: an ordinary value came
        // back nought and was drawn on the left end of its own axis under a
        // tooltip reading 0.5, and a value near the end of the domain came back
        // NaN and was not drawn at all, leaving the tooltip over blank paper.
        let source = "((A[&cf=0.5,cf_lo=0.25,cf_hi=0.75]:1,B:1):1,C:2);";
        let wide = interval(source, f64::MIN, f64::MAX);
        let marks = group(&wide, "<title>cf | estimate 0.5");
        let axis = attribute(marks, "x1=\"");
        let far_end = attribute(marks, "x2=\"");
        let estimate = attribute(marks, "<circle cx=\"");
        assert!(
            (estimate - (axis + far_end) / 2.0).abs() < 0.001,
            "half way along a domain symmetric about nought is half way along \
             the axis, not {estimate} on an axis of {axis} to {far_end}: {marks}"
        );

        // The same numbers on a domain the subtraction can hold, to show that
        // the repair is in the arithmetic and not in the picture: an estimate
        // half way up nought to one is drawn where an estimate half way up the
        // whole of f64 is.
        let readable = interval(source, 0.0, 1.0);
        let ordinary = group(&readable, "<title>cf | estimate 0.5");
        assert_eq!(
            attribute(ordinary, "<circle cx=\""),
            estimate,
            "{ordinary}\n{marks}"
        );

        // And the case that drew nothing: infinity over infinity is NaN, the
        // clamp of a NaN is a NaN, and the writer refuses a shape it cannot
        // place. The tooltip promised an estimate, so there has to be one.
        let enormous = interval(
            "((A[&cf=1e308,cf_lo=0,cf_hi=1e308]:1,B:1):1,C:2);",
            f64::MIN,
            f64::MAX,
        );
        let marks = group(&enormous, "<title>cf | estimate 1e308");
        let estimate = attribute(marks, "<circle cx=\"");
        let axis = attribute(marks, "x1=\"");
        let far_end = attribute(marks, "x2=\"");
        assert!(
            estimate > (axis + far_end) / 2.0 && estimate <= far_end,
            "an estimate past the middle of its domain belongs past the middle \
             of its axis, not at {estimate} on an axis of {axis} to {far_end}: {marks}"
        );
        assert!(!enormous.contains("NaN"), "{enormous}");
    }

    #[test]
    fn a_class_too_small_to_round_keeps_the_numbers_it_was_drawn_from() {
        // The capsule normalises the weights for its geometry, so 0.00001 and
        // 0.00001 paint exactly the picture 0.5 and 0.5 paint, and the tooltip
        // is the only place the difference survives. Rounded to four decimals it
        // did not survive either: both classes printed `weight 0`, which is a
        // weight the drawing had already refused to draw, beside a capsule split
        // in half. The rate on the line above rounded away the same way.
        let source = "((A[&o1=0.00001,w1=0.00001,o2=6,w2=0.00001]:1,B:1):1,C:2);";
        let svg = Figure::new(region())
            .width(700.0)
            .show_region_label(false)
            .push(
                TreeTrack::new(Tree::parse_annotated_newick(source).unwrap())
                    .branch_rate_mixture(BranchRateMixture::new(["o1", "o2"], ["w1", "w2"])),
            )
            .to_svg();
        assert!(
            svg.contains("class 1 omega 1e-5 weight 1e-5; class 2 omega 6 weight 1e-5"),
            "{svg}"
        );
    }
}
