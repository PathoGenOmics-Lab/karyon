# Annotated phylogenetics

Karyon can keep metadata inside a phylogeny, manipulate its topology and draw
the result on evolutionary distance or calendar time. The same metadata can
colour branches and form categorical or continuous columns aligned to the
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

The renderer is for rooted rectangular trees. It does not infer trees, fit
clocks, reconstruct ancestral states or claim epidemiological transmission.
Those analyses belong upstream; Karyon preserves their topology, lengths,
support and annotations and makes the chosen encodings explicit. Nexus support
is intentionally limited to the portable first-tree and translation-table
subset.
