//! Scale bars and time axes: the two ways a branch length is given a unit.
//!
//! Both exist because a drawn branch is a distance in pixels and the reader
//! needs it in substitutions or in years. Both have to work in all three
//! projections, and in the circular and unrooted ones a bar is a straight line
//! through a space that has no straight lines, so it is drawn where the
//! distortion is least rather than wherever there is room.

use super::*;

pub(super) fn nice_scale_length(span: f64) -> Option<f64> {
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    let target = span * 0.2;
    let magnitude = 10.0f64.powf(target.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|step| step * magnitude)
        .rfind(|value| *value <= target * (1.0 + 1e-12))
        .or(Some(magnitude))
}

pub(super) fn scale_bar_value(bar: &ScaleBar, span: f64) -> Option<f64> {
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    Some(bar.length.unwrap_or(nice_scale_length(span)?).min(span))
}

pub(super) fn draw_scale_bar(
    ctx: &mut DrawContext<'_>,
    start: (f64, f64),
    pixels: f64,
    value: f64,
    bar: &ScaleBar,
) {
    if !pixels.is_finite() || pixels <= 0.0 {
        return;
    }
    let label = match &bar.unit {
        Some(unit) => format!("{} {unit}", text_rounded(value, 3)),
        None => text_rounded(value, 3),
    };
    let title = format!("branch length scale {label}");
    let tick = ctx.theme.tokens.tick_length.max(3.0);
    let end = start.0 + pixels;
    ctx.svg.begin_titled(&title);
    ctx.svg.line(
        start.0,
        start.1,
        end,
        start.1,
        &ctx.theme.foreground,
        ctx.theme.tokens.hairline.max(1.0),
    );
    for x in [start.0, end] {
        ctx.svg.line(
            x,
            start.1 - tick / 2.0,
            x,
            start.1 + tick / 2.0,
            &ctx.theme.foreground,
            ctx.theme.tokens.hairline.max(1.0),
        );
    }
    ctx.svg.text(
        start.0,
        start.1 + tick / 2.0 + ctx.theme.font_size,
        &label,
        &ctx.theme.muted,
        (ctx.theme.font_size - 1.0).max(6.0),
        crate::svg::Anchor::Start,
    );
    ctx.svg.end_group();
}

pub(super) fn draw_rectangular_scale_bar(
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    area: Rect,
    bar: &ScaleBar,
) {
    let span = scene.maximum - scene.minimum;
    let Some(value) = scale_bar_value(bar, span) else {
        return;
    };
    draw_scale_bar(
        ctx,
        (area.x, area.bottom() + 4.0),
        area.w * value / span,
        value,
        bar,
    );
}

pub(super) fn draw_radial_scale_bar(
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    geometry: &RadialGeometry,
    area: Rect,
    bar: &ScaleBar,
) {
    let span = scene.maximum - scene.minimum;
    let Some(value) = scale_bar_value(bar, span) else {
        return;
    };
    let radial_pixels = (geometry.tree_outer - geometry.tree_inner).abs();
    draw_scale_bar(
        ctx,
        (
            area.x + 8.0,
            area.bottom() - ctx.theme.font_size - ctx.theme.tokens.tick_length - 9.0,
        ),
        radial_pixels * value / span,
        value,
        bar,
    );
}

pub(super) fn draw_unrooted_scale_bar(
    ctx: &mut DrawContext<'_>,
    scene: &UnrootedScene,
    geometry: &UnrootedGeometry,
    area: Rect,
    bar: &ScaleBar,
) {
    let Some(value) = scale_bar_value(bar, scene.radius) else {
        return;
    };
    draw_scale_bar(
        ctx,
        (
            area.x + 8.0,
            area.bottom() - ctx.theme.font_size - ctx.theme.tokens.tick_length - 9.0,
        ),
        geometry.scale * value,
        value,
        bar,
    );
}

pub(super) fn draw_time_axis(
    ctx: &mut DrawContext<'_>,
    scene: &TreeScene,
    area: Rect,
    time: &TimeAxis,
) {
    if !scene.temporal {
        return;
    }
    let y = area.bottom() + 2.0;
    ctx.svg.line(
        area.x,
        y,
        area.right(),
        y,
        &ctx.theme.foreground,
        ctx.theme.tokens.hairline,
    );
    let size = ctx.theme.font_size - 1.0;
    for index in 0..=2 {
        let fraction = index as f64 / 2.0;
        let value = scene.minimum + fraction * (scene.maximum - scene.minimum);
        let x = scene.x(area, value);
        ctx.svg.line(
            x,
            y,
            x,
            y + ctx.theme.tokens.tick_length,
            &ctx.theme.foreground,
            ctx.theme.tokens.hairline,
        );
        let label = match &time.unit {
            Some(unit) => format!("{} {unit}", text_rounded(value, 3)),
            None => text_rounded(value, 3),
        };
        ctx.svg.text(
            x,
            y + ctx.theme.tokens.tick_length + size,
            &label,
            &ctx.theme.muted,
            size,
            match index {
                0 => crate::svg::Anchor::Start,
                2 => crate::svg::Anchor::End,
                _ => crate::svg::Anchor::Middle,
            },
        );
    }
}
