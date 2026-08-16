# Geographic genomics

!!! info "Library only"
    Everything on this page is the Rust library. There is no map flag on the
    command line.

Karyon can draw supplied sample locations, geographic links and a circular
phylogeny around a map. Geography has its own coordinate system, so `Map` and
`PhyloMap` are complete `Drawing` implementations rather than genomic tracks.
They can share a `Panels` sheet with figures and rings, but they never pretend
that longitude aligns with a base position.

![Three deterministic geographic views of the same synthetic observations: an equirectangular occurrence map, a Mercator occurrence map and an orthographic globe with weighted directed links](../assets/figures/example-maps.svg)

Generate the geographic gallery with:

```bash
cargo run --example maps -- assets
```

## Draw locations and explicit links

```rust
use karyon::{GeoFlow, GeoLocation, GeoProjection, Map};

let map = Map::new()
    .title("Surveillance sites")
    .subtitle("Counts and coordinates supplied by the analysis")
    .projection(GeoProjection::orthographic(15.0, -5.0))
    .show_labels(true)
    .push(
        GeoLocation::new("Lima", -12.0464, -77.0428)
            .category("South America")
            .count(12),
    )
    .push(
        GeoLocation::new("Madrid", 40.4168, -3.7038)
            .category("Europe")
            .count(16),
    )
    .push_flow(GeoFlow::new("Lima", "Madrid").weight(4.0));

map.save_svg("sites.svg")?;
```

A location category is encoded by both colour and shape. `count` controls mark
area and remains exact in the tooltip. `value` attaches an additional numeric
value without silently turning it into a second visual scale. A flow weight
controls stroke width, and its arrowhead makes direction visible. Call
`undirected` when the input only supports a link.

Flow endpoints match exact location names. A missing name or a duplicate name
is unresolved and counted in the rendered warning; Karyon never chooses one.

### Choose the projection explicitly

| Projection | Extent | Behaviour |
|:--|:--|:--|
| `Equirectangular` | Whole world, 2:1 | Direct longitude and latitude spacing; useful when area and angle preservation are not the claim. |
| `Mercator` | Whole world, square | Preserves local angles; polar latitude is clipped rather than sent to infinity. |
| `orthographic(lat, lon)` | One hemisphere | A globe centred on the supplied decimal degrees; locations on the back are counted as outside the projection. |

`GeoPosition::is_valid` accepts only finite latitude from -90 to 90 and
longitude from -180 to 180. Invalid values are not clamped into a plausible
place. Coordinates and attached values retain their input precision in SVG
titles; only the path geometry is rounded to keep the document compact.

## Put a phylogeny around the map

![Two circular views of one synthetic outbreak phylogeny: an inward calendar tree with one connector per location and a partial cladogram with one connector per sample](../assets/figures/example-phylo-map.svg)

`PhyloMap` resolves one tree annotation against a named coordinate table. The
tree stays a tree: the map neither changes its topology nor estimates a
geographic history.

```rust
use karyon::{GeoLocation, GeoProjection, PhyloMap, Tree};

let tree = Tree::parse_annotated_newick(
    "((A[&date=2023.1,country=Peru]:0.2,\
        B[&date=2023.4,country=Peru]:0.3):0.4,\
       C[&date=2023.8,country=Spain]:0.7);",
)?;

PhyloMap::new(tree)
    .title("Circular phylogeography")
    .location_by("country")
    .coordinates([
        GeoLocation::new("Peru", -9.19, -75.0152),
        GeoLocation::new("Spain", 40.4637, -3.7492),
    ])
    .projection(GeoProjection::orthographic(15.0, -18.0))
    .time("date")
    .time_unit("year")
    .save_svg("phylogeography.svg")?;
```

The default is one aggregated connector per mapped location. Its width follows
the square root of the exact terminal count, and its title states that count.
`PhyloConnector::Individual` draws one named connector per terminal sample;
`PhyloConnector::None` keeps only the tree and map marks. These are display
choices over the same mapping table, not different analyses.

| Builder | Effect |
|:--|:--|
| `location_by(key)` | Matches the exact text of a terminal or nearest ancestral annotation to a unique coordinate name. |
| `coordinates(items)` | Supplies the named coordinate table; no geocoding or lookup is performed. |
| `time(key)` | Uses a complete finite terminal annotation as the radial scale. |
| `shape(Phylogram | Cladogram)` | Uses branch length or branch count when no time layout is selected. |
| `connector(Aggregated | Individual | None)` | Chooses the visible tree-to-place links. |
| `radial_start(degrees)` | Rotates the first terminal clockwise. |
| `radial_sweep(degrees)` | Uses a complete circle or a 120-to-360-degree fan. |
| `show_tip_labels(bool)` | Places terminal names between the tree and central map. |
| `show_time_guides(bool)` | Shows or hides concentric temporal guides. |

An annotation without a unique valid coordinate remains an unlocated tip and
is counted on the page. A valid point on the hidden orthographic hemisphere is
reported separately. The tree tip is still drawn in either case, so missing
geography cannot silently delete a sample. If a requested time layout is
incomplete, the drawing falls back to the selected phylogram or cladogram and
says that the time layout was unavailable.

Branches use only annotations already present on a node or inherited from an
annotated ancestor. Karyon does not reconstruct ancestral states, infer
migration events, geocode labels or claim transmission. Use an upstream method
for those analyses and pass its explicit result as annotations or `GeoFlow`
records.

## Compose maps with other drawings

`Map` and `PhyloMap` implement the same `Drawing` interface as `Rings` and can
be placed directly in `Panels`:

```rust
use karyon::Panels;

Panels::new()
    .columns(2)
    .push_captioned(&map, "A", "Occurrence and supplied links")
    .push_captioned(&phylo_map, "B", "Tree and reported sampling location")
    .save_svg("geographic-genomics.svg")?;
```

The land geometry is compiled into the crate from Natural Earth 1:110m land
data. [Natural Earth makes its vector data public domain](https://www.naturalearthdata.com/about/terms-of-use/).
Rendering performs no network request, depends on no tile server and is
deterministic across repeated runs.
