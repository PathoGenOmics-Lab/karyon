# Annotated phylogenetics

Karyon can keep metadata inside a phylogeny, manipulate its topology and draw
the result on evolutionary distance or calendar time. Rectangular and radial
projections use the same topology and values. The same metadata can colour
branches and form categorical or continuous columns, or rings, aligned to the
terminal taxa.

![A synthetic dated outbreak phylogeny with branches coloured by country, aligned country and sequencing-depth columns, and a second view with two named clades collapsed](../assets/figures/example-phylogenetics.svg)

The figure is deterministic and its outbreak is synthetic. Generate it with:

```bash
cargo run --example phylogenetics -- assets
```

## Read annotations instead of flattening them

There are three entry points, separated deliberately:

| Method | Input | Metadata behaviour |
|:--|:--|:--|
| `Tree::parse_newick` | Newick | Compatibility parser; comments are ignored. |
| `Tree::parse_annotated_newick` | Newick, BEAST or NHX | Preserves typed node annotations and `[&R]` / `[&U]`. |
| `Tree::parse_nexus` | Nexus trees block | Reads the first tree, applies its `translate` table and preserves annotations. |

BEAST numbers, text, booleans and brace-delimited lists become
`AnnotationValue::Number`, `Text`, `Boolean` and `List`. Access them without
re-parsing strings:

```rust
use karyon::{AnnotationValue, Tree};

let tree = Tree::parse_annotated_newick(
    "[&R] (sample_A[&date=2024.25,country=Peru,selected=true]:0.2,\
            sample_B[&date=2024.50,country=Spain]:0.3);",
)?;
let sample = tree.node_named("sample_A").unwrap();

assert_eq!(
    tree.annotation(sample, "date").and_then(AnnotationValue::as_number),
    Some(2024.25),
);
assert_eq!(tree.rooted(), Some(true));
```

`annotations_mut` and `tree_annotations_mut` add or replace metadata after
parsing. Node indices are stable across rotation, ladderising and rerooting;
operations that extract or delete nodes return a compact tree with new indices.

## Draw time, branches and sample traits together

```rust
use karyon::{Figure, Region, TraitColumn, TreeTrack};

let track = TreeTrack::new(tree)
    .time("date")
    .time_unit("year")
    .color_by("country")
    .show_nodes(true)
    .trait_column(
        TraitColumn::categorical("country")
            .label("Country")
            .width(62.0),
    )
    .trait_column(
        TraitColumn::continuous("coverage")
            .label("Depth")
            .width(46.0),
    );

Figure::new(Region::new("phylogeny", 0, 1)?)
    .show_region_label(false)
    .push(track)
    .save_svg("outbreak.svg")?;
```

`color_by` uses a continuous ramp when every visible value is numeric and the
categorical palette otherwise. A branch without its own value inherits the
nearest annotated ancestor. Its exact value remains in the SVG title, so the
colour is not the only way to retrieve it.

Trait columns use the same categorical mapping as branches in that view.
Continuous cells span the visible numeric range. Missing or non-numeric values
in a continuous column are outlined and marked with an em dash rather than
invented as zero. `show_values(false)` keeps only the colour blocks when space
is tight.

## Change the projection, not the tree

![Four views of the same synthetic outbreak phylogeny: an outward circular time tree with two trait rings, a partial fan with a collapsed clade, an inward time tree and a circular cladogram](../assets/figures/example-phylo-layouts.svg)

`circular` wraps the terminal order around a complete circumference. Branch
lengths, dates, annotations, labels and node order are not recomputed. Calendar
ticks become concentric guides, trait columns become annular rings and a
visually collapsed clade becomes a wedge.

```rust
use karyon::{RadialDirection, TraitColumn, TreeTrack};

let outward = TreeTrack::new(tree.clone())
    .time("date")
    .color_by("country")
    .trait_column(
        TraitColumn::categorical("country")
            .label("Country")
            .ring_width(12.0),
    )
    .circular()
    .radial_start(-90.0)
    .radial_size(520.0);

let inward_fan = TreeTrack::new(tree)
    .time("date")
    .fan(250.0)
    .radial_start(-215.0)
    .radial_direction(RadialDirection::Inward)
    .inner_radius(0.32);
```

| Builder | Visual effect |
|:--|:--|
| `circular()` | Complete 360-degree projection. |
| `fan(degrees)` | Partial clockwise sweep from 10 to 359 degrees. |
| `radial_start(degrees)` | Rotates the first terminal; `-90` starts at twelve o'clock. |
| `radial_sweep(degrees)` | Sets a 10-to-360-degree sweep directly. |
| `radial_direction(Outward | Inward)` | Places terminals at the outside or towards the centre. |
| `inner_radius(fraction)` | Reserves a central gap from 0 to 0.85 of the tree radius. |
| `radial_size(pixels)` | Sets radial drawing height independently of the tip count. |
| `projection(TreeProjection::Circular)` | Selects circular coordinates without changing the other radial settings. |
| `TraitColumn::ring_width(pixels)` | Sets a trait ring's thickness from 2 to 24 pixels. |

A full circle is useful when topology and metadata are the subject. A fan leaves
a quiet sector for labels or an external annotation. An inward projection keeps
early branches around the circumference and later nodes near a controlled
central gap. For a dense tree, turn off exact ring text with
`show_values(false)`, omit tip labels with `show_tips(false)`, or hide the
concentric time guides with `show_time_axis(false)`. Tooltips still retain exact
sample and annotation values.

Radial layouts occupy a standalone square and therefore do not share rows with
`SnpTrack`, `MatrixTrack` or `CladeTrack`. Use the rectangular projection when
leaf-to-row alignment is the analytical claim.

### Requirements for a time tree

Every tip must carry a finite numeric value for the key passed to `time`.
Annotated internal values are used exactly. Missing internal values are
inferred from the children and branch lengths: lengths are subtracted for
ordinary calendar dates and added for heights before present.

```rust
use karyon::TimeDirection;

let track = TreeTrack::new(tree)
    .time("height")
    .time_direction(TimeDirection::Decreasing)
    .time_unit("years BP");
```

`Tree::time_layout` returns `None` when a tip value is missing or non-finite.
`TreeTrack` then falls back to its ordinary phylogram and omits the temporal
axis; validate with `time_layout` first when incomplete dates must be an error
in an analysis pipeline.

## Work with clades and topology

All operations are iterative, including deep trees.

| Operation | Effect |
|:--|:--|
| `ancestors`, `descendants`, `clade_size` | Query the rooted topology. |
| `mrca` | Find the most recent common ancestor of a non-empty node set. |
| `rotate` | Reverse one split without changing its clades. |
| `ladderize` | Order every split by descendant tip count. |
| `reroot` | Reorient around an internal node while preserving tips and edge lengths. |
| `subtree` | Copy one clade into a compact standalone tree. |
| `collapse` | Replace descendants in the data with one terminal node. |
| `TreeTrack::collapse` | Draw a clade as a triangle without modifying the source tree. |

Visual collapse uses the selected node's own or inherited annotations. It does
not calculate a mean, majority category or any other summary for descendants;
if the collapsed node has no sequencing depth, the aligned depth cell is
explicitly missing.

```rust
let outbreak = tree.node_named("PER_outbreak").unwrap();
let track = TreeTrack::new(tree).collapse(outbreak);

assert_eq!(track.tree().clade_size(outbreak), 4);
```

## Scope

The renderer is for rooted rectangular, circular and fan trees. It does not
infer trees, fit clocks, reconstruct ancestral states or claim epidemiological
transmission. [`PhyloMap`](maps.md#put-a-phylogeny-around-the-map) can place
terminal annotations at explicitly supplied coordinates, but it does not infer
those locations or the movement between them. Those analyses belong upstream;
Karyon preserves their topology, lengths, support and annotations and makes the
chosen encodings explicit. Nexus support is intentionally limited to the
portable first-tree and translation-table subset.
