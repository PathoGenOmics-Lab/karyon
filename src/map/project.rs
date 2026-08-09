//! Turning a latitude and a longitude into a point, or refusing to.
//!
//! Every function here can answer that it will not: a coordinate off the Earth
//! and a coordinate on the hidden side of a globe are both `None` rather than
//! a plausible-looking point somewhere. That is the whole reason the module is
//! separate. A projection that silently returns a number puts a sample
//! somewhere nobody sampled, and the map looks right.

use super::*;

pub(super) fn project(
    position: GeoPosition,
    projection: GeoProjection,
    area: MapRect,
) -> Option<(f64, f64)> {
    if !position.is_valid() {
        return None;
    }
    match projection {
        GeoProjection::Equirectangular | GeoProjection::Mercator => {
            project_rectangular(position, projection, area, true)
        }
        GeoProjection::Orthographic { .. } => {
            let (x, y) = orthographic_unit(position, projection)?;
            let radius = area.w.min(area.h) / 2.0;
            Some((
                area.x + area.w / 2.0 + x * radius,
                area.y + area.h / 2.0 - y * radius,
            ))
        }
    }
}

pub(super) fn project_rectangular(
    position: GeoPosition,
    projection: GeoProjection,
    area: MapRect,
    require_bounds: bool,
) -> Option<(f64, f64)> {
    if !position.latitude.is_finite() || !position.longitude.is_finite() {
        return None;
    }
    if require_bounds && !position.is_valid() {
        return None;
    }
    let x = area.x + (position.longitude + 180.0) / 360.0 * area.w;
    let y_fraction = match projection {
        GeoProjection::Equirectangular => (90.0 - position.latitude) / 180.0,
        GeoProjection::Mercator => {
            let latitude = position.latitude.clamp(-85.051_128_78, 85.051_128_78);
            let radians = latitude.to_radians();
            (1.0 - (radians.tan() + 1.0 / radians.cos()).ln() / PI) / 2.0
        }
        GeoProjection::Orthographic { .. } => return None,
    };
    Some((x, area.y + y_fraction * area.h))
}

pub(super) fn orthographic_visibility(position: GeoPosition, projection: GeoProjection) -> f64 {
    let GeoProjection::Orthographic {
        center_latitude,
        center_longitude,
    } = projection
    else {
        return 1.0;
    };
    let (latitude, longitude) = position.radians();
    let centre_latitude = center_latitude.to_radians();
    let delta = longitude - center_longitude.to_radians();
    centre_latitude.sin() * latitude.sin() + centre_latitude.cos() * latitude.cos() * delta.cos()
}

pub(super) fn orthographic_unit(
    position: GeoPosition,
    projection: GeoProjection,
) -> Option<(f64, f64)> {
    let GeoProjection::Orthographic {
        center_latitude,
        center_longitude,
    } = projection
    else {
        return None;
    };
    if orthographic_visibility(position, projection) < -1e-9 {
        return None;
    }
    let (latitude, longitude) = position.radians();
    let centre_latitude = center_latitude.to_radians();
    let delta = longitude - center_longitude.to_radians();
    Some((
        latitude.cos() * delta.sin(),
        centre_latitude.cos() * latitude.sin()
            - centre_latitude.sin() * latitude.cos() * delta.cos(),
    ))
}

pub(super) fn horizon_intersection(
    a: GeoPosition,
    b: GeoPosition,
    projection: GeoProjection,
) -> Option<GeoPosition> {
    let da = orthographic_visibility(a, projection);
    let db = orthographic_visibility(b, projection);
    let denominator = da - db;
    if denominator.abs() < 1e-12 {
        return None;
    }
    let t = (da / denominator).clamp(0.0, 1.0);
    let a = sphere(a);
    let b = sphere(b);
    let mut point = (
        a.0 + (b.0 - a.0) * t,
        a.1 + (b.1 - a.1) * t,
        a.2 + (b.2 - a.2) * t,
    );
    let length = (point.0 * point.0 + point.1 * point.1 + point.2 * point.2).sqrt();
    if length <= 1e-12 {
        return None;
    }
    point.0 /= length;
    point.1 /= length;
    point.2 /= length;
    Some(GeoPosition::new(
        point.2.asin().to_degrees(),
        point.1.atan2(point.0).to_degrees(),
    ))
}

pub(super) fn sphere(position: GeoPosition) -> (f64, f64, f64) {
    let (latitude, longitude) = position.radians();
    (
        latitude.cos() * longitude.cos(),
        latitude.cos() * longitude.sin(),
        latitude.sin(),
    )
}

pub(super) fn path_from_points(points: &[(f64, f64)], close: bool) -> Option<String> {
    let first = points.first()?;
    if points.len() < if close { 3 } else { 2 } {
        return None;
    }
    let mut path = format!("M {} {}", num(first.0), num(first.1));
    for point in points.iter().skip(1) {
        path.push_str(&format!(" L {} {}", num(point.0), num(point.1)));
    }
    if close {
        path.push_str(" Z");
    }
    Some(path)
}

pub(super) fn circle_path(area: MapRect) -> String {
    let cx = area.x + area.w / 2.0;
    let cy = area.y + area.h / 2.0;
    let radius = area.w.min(area.h) / 2.0;
    format!(
        "M {} {} A {} {} 0 1 1 {} {} A {} {} 0 1 1 {} {} Z",
        num(cx - radius),
        num(cy),
        num(radius),
        num(radius),
        num(cx + radius),
        num(cy),
        num(radius),
        num(radius),
        num(cx - radius),
        num(cy)
    )
}
