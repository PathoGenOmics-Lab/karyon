//! Renders the deterministic geographic gallery used in the documentation.
//!
//! ```text
//! cargo run --example maps -- assets
//! ```

use std::env;
use std::path::PathBuf;

use karyon::{GeoFlow, GeoLocation, GeoProjection, Map, Panels};

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let global = Map::new()
        .title("Global surveillance sites")
        .subtitle("Synthetic counts; exact coordinates")
        .width(650.0)
        .height(430.0)
        .show_labels(true)
        .extend(sites());

    let mercator = Map::new()
        .title("The same observations in Mercator")
        .subtitle("Projection changes; locations do not")
        .projection(GeoProjection::Mercator)
        .width(650.0)
        .height(430.0)
        .extend(sites());

    let flows = Map::new()
        .title("Directed introductions")
        .subtitle("User-provided links, not inferred transmission")
        .projection(GeoProjection::orthographic(15.0, -5.0))
        .width(650.0)
        .height(540.0)
        .show_labels(true)
        .extend(
            sites()
                .into_iter()
                .filter(|site| !matches!(site.name(), "Tokyo" | "Sydney" | "Bangkok")),
        )
        .push_flow(
            GeoFlow::new("Lima", "Madrid")
                .category("Introduction")
                .weight(4.0),
        )
        .push_flow(
            GeoFlow::new("Cape Town", "Madrid")
                .category("Introduction")
                .weight(2.0),
        )
        .push_flow(
            GeoFlow::new("Madrid", "Reykjavik")
                .category("Onward movement")
                .weight(3.0),
        )
        .push_flow(
            GeoFlow::new("Madrid", "Nairobi")
                .category("Onward movement")
                .weight(1.0),
        );

    let gallery = Panels::new()
        .title("Geographic genomics (synthetic)")
        .columns(2)
        .gap(24.0)
        .push_captioned(
            &global,
            "A",
            "A full-world occurrence map with categorical shape and colour",
        )
        .push_captioned(
            &mercator,
            "B",
            "One dataset under a second explicit projection",
        )
        .push_captioned(
            &flows,
            "C",
            "A hemisphere view with weighted, directed geographic links",
        );

    gallery.save_svg(out.join("example-maps.svg"))?;
    let (width, height) = gallery.dimensions();
    println!("example-maps.svg {width:.0} x {height:.0}, 3 map views");
    Ok(())
}

fn sites() -> Vec<GeoLocation> {
    vec![
        GeoLocation::new("Lima", -12.0464, -77.0428)
            .category("South America")
            .value(0.18)
            .count(12),
        GeoLocation::new("São Paulo", -23.5505, -46.6333)
            .category("South America")
            .value(0.31)
            .count(7),
        GeoLocation::new("New York", 40.7128, -74.0060)
            .category("North America")
            .value(0.24)
            .count(9),
        GeoLocation::new("Madrid", 40.4168, -3.7038)
            .category("Europe")
            .value(0.42)
            .count(16),
        GeoLocation::new("Reykjavik", 64.1466, -21.9426)
            .category("Europe")
            .value(0.12)
            .count(4),
        GeoLocation::new("Nairobi", -1.2921, 36.8219)
            .category("Africa")
            .value(0.36)
            .count(10),
        GeoLocation::new("Cape Town", -33.9249, 18.4241)
            .category("Africa")
            .value(0.28)
            .count(6),
        GeoLocation::new("Bangkok", 13.7563, 100.5018)
            .category("Asia-Pacific")
            .value(0.33)
            .count(8),
        GeoLocation::new("Tokyo", 35.6762, 139.6503)
            .category("Asia-Pacific")
            .value(0.39)
            .count(11),
        GeoLocation::new("Sydney", -33.8688, 151.2093)
            .category("Asia-Pacific")
            .value(0.15)
            .count(5),
    ]
}
