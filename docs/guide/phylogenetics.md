# Annotated phylogenetics

Karyon can keep metadata inside a phylogeny, manipulate its topology and draw
the result on evolutionary distance or calendar time. Rectangular, radial and
unrooted projections use the same topology and values. The same metadata can
colour branches and form colour strips, heatmaps, bars, binary marks or shaped
categories aligned to the terminal taxa.

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

## Choose a tree geometry for the reading task

![Eight synthetic evolutionary views beginning with orthogonal, diagonal, curved, circular and unrooted tree geometries](../assets/figures/example-evolutionary-surveillance.svg)

Geometry is a reading aid, not a transformation of the owned tree. The three
rectangular branch shapes keep the same node coordinates and terminal order;
the radial and unrooted projections change the coordinate system while
retaining topology, branch values and annotation ownership.

```rust
use karyon::{BranchGeometry, TreeTrack};

let aligned = TreeTrack::new(tree.clone())
    .branch_geometry(BranchGeometry::Orthogonal);
let topology_forward = TreeTrack::new(tree.clone())
    .branch_geometry(BranchGeometry::Diagonal);
let presentation = TreeTrack::new(tree)
    .branch_geometry(BranchGeometry::Curved);
```

| Geometry | Best reading | Constraint |
|:--|:--|:--|
| orthogonal rectangular | aligned tip rows, events and dense metadata columns | parent risers can dominate very unbalanced trees |
| diagonal rectangular | topology and branch-length direction | weaker visual alignment between a node and its descendants |
| curved rectangular | annotated internal nodes and presentation figures | use restrained node glyph sizes to avoid crossings |
| circular | many terminal taxa plus metadata rings | root and terminal order remain meaningful |
| fan | radial context with a quiet sector for annotation | partial sweep gives taxa unequal screen directions, not unequal evolutionary weight |
| unrooted | split structure without privileging the source root | cannot show rooted time direction |

`branch_geometry` is intentionally ignored outside rectangular coordinates.
Rerooting, ladderising and rotating are separate topology operations; choosing
a path shape never performs one of them implicitly.

## Draw topology without privileging the Newick root

![An unrooted tree with a layered metadata halo beside a circular cladogram carrying the same four annotation datasets](../assets/figures/example-phylo-annotations.svg)

`unrooted()` treats every branch as an undirected edge, chooses a centre that
minimises the largest terminal-bearing component and assigns equal angular
space to terminal taxa. The root stored in the source tree is not moved or
deleted. A phylogram keeps branch lengths; a cladogram gives every edge one
unit.

```rust
use karyon::{TraitColumn, TreeShape, TreeTrack};

let view = TreeTrack::new(tree)
    .shape(TreeShape::Phylogram)
    .unrooted()
    .unrooted_start(-104.0)
    .unrooted_size(560.0)
    .color_by("country")
    .show_nodes(true)
    .trait_column(TraitColumn::categorical("country").label("Country"))
    .trait_column(TraitColumn::bar("coverage").label("Depth"))
    .trait_column(TraitColumn::binary("resistant").label("AMR"))
    .trait_column(TraitColumn::symbol("host").label("Host"));
```

Terminal leaders connect unequal branch endpoints to one common annotation
halo. They are guides, not extra evolutionary distance. Branch and ring
colours share one categorical domain, so a country cannot change colour
between the topology and its metadata.

`unrooted_size` fixes the standalone height and `unrooted_start` rotates the
equal-angle sectors. Time axes are rooted quantities and are therefore not
drawn in this projection. Use rectangular or circular coordinates when root
age, direction or calendar time is part of the claim.

## Put support, events and distance on the branches

![The same synthetic phylogram drawn in rectangular, circular and unrooted coordinates with scaled support markers, exact support labels, mutation labels and branch-length scale bars](../assets/figures/example-phylo-evidence.svg)

Support, an event and branch length answer different questions, so Karyon gives
each one an independent channel. Support uses node markers and optional text;
an event follows the branch that owns it; evolutionary distance gets a scale
bar rather than being inferred from panel width.

```rust
use karyon::{SupportStyle, TreeTrack};

let view = TreeTrack::new(tree)
    .support_style(SupportStyle::SymbolsAndLabels)
    .support_threshold(0.70)
    .branch_labels("mutation")
    .branch_label_size(7.0)
    .scale_bar()
    .scale_bar_length(0.1)
    .scale_bar_unit("substitutions/site");
```

`SupportStyle::Symbols`, `Labels`, `SymbolsAndLabels` and `None` control only
the visible encoding. Exact support remains in branch tooltips. Thresholds can
use either the 0–1 convention (`0.70`) or the percentage convention (`70.0`);
labels retain the value as supplied rather than silently converting it.

`branch_labels` reads only the annotation attached to the incoming branch. It
does not inherit ancestral values as `color_by` does, because a mutation, gain
or loss must not be repeated on every descendant. Labels rotate with circular
and unrooted edges. When an edge is too short, visible text is ellipsised and
the complete key and value remain in the SVG tooltip.

`scale_bar()` chooses a 1–2–5 length near one fifth of the visible branch span.
`scale_bar_length` requests a value explicitly and clamps it to that span;
`scale_bar_unit` prints its unit exactly. Scale bars are omitted from
cladograms and explicitly time-scaled trees, where a branch-length ruler would
make the wrong claim.

## Render ancestral states, events and branch uncertainty

Ancestral reconstruction usually contains at least three different objects:
a probability distribution at a node, an inferred event on an edge and an
uncertainty measure for an estimate. Combining them into one branch colour
loses both ownership and uncertainty. These layers keep them separate in every
tree projection.

```rust
use karyon::{
    AncestralStateLayer, BranchEventLayer, BranchGeometry,
    BranchIntervalLayer, TreeTrack,
};

let reconstruction = TreeTrack::new(tree)
    .branch_geometry(BranchGeometry::Curved)
    .ancestral_states(
        AncestralStateLayer::new(["state_human", "state_animal", "state_water"])
            .label("ancestral host posterior")
            .confidence(0.72),
    )
    .branch_event_layer(
        BranchEventLayer::new("mutations")
            .label("ancestral mutations")
            .maximum_events(6),
    )
    .branch_interval(
        BranchIntervalLayer::new("gcf", "gcf_low", "gcf_high")
            .label("gene concordance")
            .range(0.0, 1.0)
            .threshold(0.70),
    );
```

The ancestral layer normalises non-negative supplied probabilities only for
donut geometry. Exact values remain in its tooltip. A transition cue appears
only when the maximum-posterior state changes and both endpoint maxima reach
the confidence floor. It is a visualisation of the supplied reconstruction,
not a newly inferred transition.

Text, numbers and booleans under a `BranchEventLayer` key become one direct
event; an annotated-Newick list becomes ordered marks, capped per branch by
`maximum_events`. `BranchIntervalLayer` draws the supplied estimate and valid
lower/upper bounds on a fixed compact scale. Reversed or non-finite intervals
are omitted rather than repaired. None of these layers inherits an ancestor's
value onto descendants.

## Show branch-wise dN/dS without moving the neutral point

![Four views of one synthetic codon-model tree: a phylogram, a circular tree with genomic metadata rings, an unrooted tree and a cladogram, all with branch-wise dN/dS](../assets/figures/example-phylo-dnds.svg)

`dnds` is a dedicated branch encoding for the nonsynonymous-to-synonymous rate
ratio, ω. It differs from a generic continuous `color_by` in two important
ways: the diverging scale is fixed around the biologically meaningful value
ω = 1, and an estimate belongs only to its incoming branch. It is never
inherited by descendants.

```rust
use karyon::TreeTrack;

let view = TreeTrack::new(tree)
    .dnds("omega")
    .dnds_label("Branch dN/dS (ω)")
    .dnds_neutral_band(0.9, 1.1)
    .dnds_saturation(4.0)
    .dnds_significance("q", 0.05)
    .branch_labels("amino_acid_change");
```

The cool side represents ω below the neutral band, grey represents values near
one and the warm side represents ω above it. Colour strength follows
`abs(log2(ω))` and saturates symmetrically: with `dnds_saturation(4.0)`, ω ≤
0.25 and ω ≥ 4 use the strongest colours. Zero is retained as the strongest
purifying value; negative, non-finite and missing estimates are drawn as quiet
dotted edges rather than converted to zero.

`dnds_significance(key, maximum)` adds an independent evidence channel. A
branch whose direct numeric `p`, `q` or other test value is at most the chosen
threshold becomes thicker, while its colour continues to describe effect
size. Exact ω and test values, the selected regime and the threshold comparison
remain in the SVG tooltip in rectangular, circular and unrooted projections.

The renderer visualises estimates fitted upstream; it does not calculate dN,
dS, likelihood-ratio tests or multiple-testing corrections. In particular,
ω > 1 alone is shown as a diversifying regime, not presented as proof of
positive selection. Generate the synthetic gallery above with:

```bash
cargo run --example phylo_dnds -- assets
```

## Build a branch-to-codon selection atlas

![A four-panel synthetic molecular-selection atlas with weighted branch rate classes, recurrent-event connections, a circular dN/dS tree, a frequentist site scan and a posterior site scan](../assets/figures/example-selection-atlas.svg)

A single mean ω is sometimes the result, and sometimes the summary that hides
the result. Branch-site models can fit several rate classes to one branch,
while site models report evidence and effect at coding positions. Karyon keeps
those quantities in separate visual channels and lets them meet only through
their shared biological interpretation.

### Preserve fitted rate classes on branches

`BranchRateMixture` pairs any number of direct rate annotations with their
weight annotations. Segment length is the fitted class weight; segment colour
is the class ω on the same neutral-centred logarithmic scale as `dnds`. Weights
are normalised only to fill the capsule. Their original values remain exact in
the SVG tooltip.

```rust
use karyon::{BranchRateMixture, HomoplasyLayer, TreeTrack};

let rates = BranchRateMixture::new(
    ["omega_1", "omega_2", "omega_3"],
    ["weight_1", "weight_2", "weight_3"],
)
.label("aBSREL ω classes")
.neutral_band(0.9, 1.1)
.saturation(6.0);

let view = TreeTrack::new(tree)
    .branch_rate_mixture(rates)
    .homoplasy_layer(
        HomoplasyLayer::new("amino_acid_change")
            .label("recurrent amino-acid change"),
    );
```

The rate and weight keys are paired in iterator order. A class with a missing,
negative or non-finite rate, or with a non-positive weight, is omitted. A
branch with no valid class is left untouched rather than receiving a zero-rate
capsule. Values are read from the node that owns the incoming edge and are
never inherited.

`HomoplasyLayer` groups equal direct branch annotations. Events appearing on
at least two branches are joined with dashed curves: contained arcs in the
rectangular tree and centre-seeking chords in circular and unrooted trees.
`minimum_occurrences` raises the recurrence threshold and
`maximum_connections` prevents a common event from turning a dense tree into
an all-to-all web. The renderer calls these *recurrent events*, not proven
homoplasies: convergence, reversal and ancestral-state uncertainty have to be
settled by the upstream analysis.

### Separate site evidence from rate direction

`SelectionTrack` uses genomic x coordinates, so it can sit under protein
domains, codons, variants or an axis. Its upper tier draws p-values as
`-log10(p)` or posterior probability on `0..1`; its lower tier draws signed
`log2(ω)` around the explicit neutral baseline. A diamond means the chosen
evidence threshold was crossed. Colour still means purifying, near-neutral or
diversifying rate, so a strongly supported purifying site is not painted as a
positive-selection hit.

```rust
use karyon::{SelectionEvidence, SelectionSite, SelectionTrack};

let sites = vec![
    SelectionSite::new(44)
        .rates(0.18, 1.52)
        .p_value(0.0014)
        .episodic_rates(0.05, 3.8, 0.18)
        .label("surface loop"),
    SelectionSite::new(103)
        .rates(0.50, 0.07)
        .p_value(0.008),
];

let scan = SelectionTrack::new(sites)
    .evidence(SelectionEvidence::PValue)
    .p_threshold(0.05)
    .neutral_band(0.85, 1.15)
    .saturation(8.0)
    .label("FEL / MEME");
```

Switch to `SelectionEvidence::Posterior` and set
`posterior_threshold` for FUBAR-like posterior scans. A site can carry both
forms of evidence; the track-level mode decides which one is drawn, preventing
p-values and posterior probabilities from sharing a false common axis.
`episodic_rates(beta_minus, beta_plus, positive_weight)` adds a compact
two-class capsule above the point while preserving every supplied value in its
tooltip.

Positions are 0-based coordinates, consistent with every genomic track in the
crate. `dS = 0, dN > 0` remains an infinite ratio in the tooltip and saturates
safely in geometry; a missing rate pair remains missing. The renderer does not
run FEL, MEME, FUBAR, aBSREL, ancestral reconstruction or multiple-testing
correction. It renders their supplied results without silently converting one
statistical quantity into another.

Generate the complete synthetic atlas with:

```bash
cargo run --example selection_atlas -- assets
```

## Layer annotation rings like iTOL datasets

`TraitColumn` uses the same dataset in rectangular columns, circular rings and
the halo around an unrooted tree. The mark changes with the projection; the
annotation key and its exact SVG tooltip do not.

| Builder | Rectangular mark | Circular or unrooted mark | Accepted value |
|:--|:--|:--|:--|
| `categorical(key)` | colour strip | annular colour strip | any typed value |
| `continuous(key)` | heatmap cell | annular heatmap sector | finite number |
| `bar(key)` | horizontal bar | outward radial bar | finite number |
| `binary(key)` | presence marker | ring marker | boolean or finite number; zero is absent |
| `symbol(key)` | coloured shape | coloured ring shape | any typed value |

Missing values remain outlined. Text is never guessed as binary, numeric
values are scaled only across the visible dataset, and symbols repeat category
identity with shape as well as colour. `ring_width` controls each annular
dataset independently; `show_values(false)` removes in-cell text without
removing tooltips.

## Attach data graphics to nodes and clades

![Four synthetic phylogenetic data views: a rectangular tree with abundance bubbles, stacked host bars and a highlighted clade; a circular tree with ancestral-state donuts and a clade sector; a tree-aligned nucleotide alignment; and tree-aligned protein domain architectures](../assets/figures/example-phylo-faces.svg)

`NodeGlyph` turns numeric node annotations into small plots without flattening
them into labels. Bubble area follows one value; pie, donut and stacked-bar
segments follow several keys in the order supplied. Composition geometry is
normalised locally, while the tooltip retains every original value.

```rust
use karyon::{CladeHighlight, NodeGlyph, NodeGlyphTarget, TreeTrack};

let outbreak = tree.node_named("outbreak").unwrap();
let track = TreeTrack::new(tree)
    .node_glyph(
        NodeGlyph::bubble("isolates")
            .label("Isolate count")
            .target(NodeGlyphTarget::Internal),
    )
    .node_glyph(
        NodeGlyph::donut(["human", "animal", "environment"])
            .label("Host probability")
            .target(NodeGlyphTarget::Internal),
    )
    .clade_highlight(
        CladeHighlight::new(outbreak)
            .label("Transmission cluster")
            .opacity(0.12),
    );
```

| Constructor | Data requirement | Encoding |
|:--|:--|:--|
| `NodeGlyph::bubble(key)` | one finite, non-negative number | circle area |
| `NodeGlyph::pie(keys)` | one finite, non-negative number per key | filled sectors |
| `NodeGlyph::donut(keys)` | one finite, non-negative number per key | annular sectors |
| `NodeGlyph::stacked_bar(keys)` | one finite, non-negative number per key | compact horizontal composition |

`NodeGlyphTarget::All`, `Internal` and `Leaves` prevent a dataset from being
repeated where it has no biological meaning. A missing key suppresses that
node's glyph rather than treating absence as zero. `CladeHighlight` becomes a
descendant band in rectangular coordinates, an annular sector in radial
coordinates and a topology-following field in an unrooted view. Its tooltip
always reports the exact descendant-tip count.

## Align sequences and domain architectures to descent

`MsaTrack::tree` and `DomainTrack::tree` match rows to leaves by exact name,
sort them by descent and draw the tree in the same gutter. A row not named by
the tree remains at the bottom instead of disappearing.

```rust
use karyon::{DomainArchitecture, DomainFeature, DomainTrack, MsaTrack};

let alignment = MsaTrack::new(sequences)
    .tree(tree.clone())
    .tree_width(110.0);

let architectures = vec![
    DomainArchitecture::new("sample_A", 300)
        .feature(DomainFeature::new(20, 110).label("sensor"))
        .feature(DomainFeature::new(170, 260).label("kinase")),
];
let domains = DomainTrack::new(architectures)
    .tree(tree)
    .tree_width(110.0);
```

Domain and motif boundaries remain 0-based and half-open. Colours are stable by
feature label, explicit colours override the palette, and full names and
boundaries remain in tooltips when visible text must be shortened. The renderer
does not infer domains or ancestral states; it displays intervals and numeric
probabilities supplied by an upstream analysis.

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

## Choose the root explicitly

![The same synthetic phylogeny with its source root, a checked monophyletic outgroup root and a weighted midpoint root; a diamond identifies each selected root](../assets/figures/example-phylo-reroot.svg)

Rerooting changes orientation, not the undirected tree. Karyon preserves every
tip-to-tip distance, keeps support on the same split and appends a new root only
when the chosen position lies inside an edge. A diamond marks the selected root
in rectangular and circular projections.

```rust
use karyon::TreeTrack;

let by_clade = TreeTrack::new(tree.clone()).reroot_named("lineage_4");
let by_outgroup = TreeTrack::new(tree.clone())
    .reroot_outgroup(["outgroup_A", "outgroup_B"]);
let by_midpoint = TreeTrack::new(tree).reroot_midpoint();
```

| Builder | Validation and result |
|:--|:--|
| `reroot(node)` | Accepts an internal node index; a sampled tip or invalid index leaves the tree unchanged. |
| `reroot_named(name)` | Finds one exact internal label and uses that node. |
| `reroot_outgroup(names)` | Requires existing, distinct leaf names that are exactly one monophyletic clade; inserts a root halfway along its incoming edge. |
| `reroot_midpoint()` | Requires every edge to have a finite, non-negative length; bisects the longest weighted tip-to-tip path. |
| `show_root(false)` | Hides the diamond without undoing the reroot. |

The builder API deliberately leaves an invalid request unchanged so it remains
composable. Pipelines that must treat failure as an error should call
`Tree::reroot`, `Tree::reroot_outgroup` or `Tree::reroot_midpoint` first and
inspect their `bool` or `Option<usize>` result before constructing the track.
An unrooted projection never draws the diamond because its geometry explicitly
discards the source root.

## Work with clades and topology

All operations are iterative, including deep trees.

| Operation | Effect |
|:--|:--|
| `ancestors`, `descendants`, `clade_size` | Query the rooted topology. |
| `mrca` | Find the most recent common ancestor of a non-empty node set. |
| `rotate` | Reverse one split without changing its clades. |
| `ladderize` | Order every split by descendant tip count. |
| `reroot` | Reorient around an internal node while preserving tips, edge lengths and split support. |
| `reroot_outgroup` | Validate a monophyletic leaf set and insert a root on its incoming edge. |
| `reroot_midpoint` | Bisect the weighted diameter when every branch length is valid. |
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

## Make tree-to-tree disagreement traceable

![Core and accessory genome trees with duplicated terminal labels, ward-coloured ties, dashed crossings and a before-to-after untangling summary](../assets/figures/example-tanglegram.svg)

A tanglegram compares terminal order, not branch coordinates. `untangle`
rotates free clades on both trees and accepts a rotation only when the crossing
count strictly falls. It preserves every clade, annotation and branch length;
the result is deterministic but is not presented as a global optimum.

```rust
use karyon::{TangleLabels, TangleTieStyle, TanglegramTrack};

let comparison = TanglegramTrack::new(core, accessory)
    .names("core genome", "accessory genome")
    .labels(TangleLabels::Both)
    .tie_style(TangleTieStyle::Curved)
    .color_by("ward")
    .untangle();

assert!(comparison.crossings() <= comparison.initial_crossings());
```

The central summary reports initial and final crossings, linked taxa and tips
present in only one tree. Crossing ties use a dashed pattern, so colour remains
available for metadata. When the two trees give a matched taxon different
values for the selected annotation, the endpoint marks retain both colours and
the tooltip states the exact mismatch.

`TangleLabels::Left`, `Right`, `Both` and `None` control repeated terminal
names. `TangleTieStyle::Straight` is compact; `Curved` is easiest to trace; and
`Ribbon` remains visible after reduction for print. `tie_widths`, `tree_width`,
`label_width` and `row_height` control density without changing the comparison.

## Align phylodynamics and surveillance over time

The bottom of the evolutionary-surveillance atlas deliberately places an
inferred process above observed composition. Their x pivots agree; their y
quantities, uncertainty and evidential status do not.

```rust
use karyon::{
    AxisTrack, Figure, PhylodynamicPoint, PhylodynamicScale,
    PhylodynamicTrack, Region, SurveillanceObservation, SurveillanceTrack,
};

let skyline = PhylodynamicTrack::new(vec![
    PhylodynamicPoint::new(2020, 120.0).interval(70.0, 210.0),
    PhylodynamicPoint::new(2021, 430.0).interval(250.0, 760.0),
])
.scale(PhylodynamicScale::Log10)
.unit("Ne");

let observed = SurveillanceTrack::new(vec![
    SurveillanceObservation::new(2020, "L1", 38, 100),
    SurveillanceObservation::new(2020, "L2", 62, 100),
    SurveillanceObservation::new(2021, "L1", 73, 120),
    SurveillanceObservation::new(2021, "L2", 47, 120),
])
.minimum_total(20)
.frequency_alert(0.50)
.growth_alert(0.15);

Figure::new(Region::new("year", 2020, 2022)?)
    .push(skyline)
    .push(observed)
    .push(AxisTrack::new())
    .save_svg("evolution-through-time.svg")?;
```

`PhylodynamicTrack` retains point intervals as a ribbon and supports linear or
base-ten logarithmic y geometry. Non-positive estimates are absent in log mode
instead of being nudged above zero. `SurveillanceTrack` retains each count and
denominator, can switch between frequencies and raw counts, and draws stacked
composition or independent lines. A sampling floor omits underpowered
observations explicitly; alert markers report whether frequency, stepwise
growth or both crossed the chosen rule. Missing lineage/time pairs break a
line and make a stacked pivot explicitly incomplete; zero must be supplied as
an observed zero. Duplicate pairs are marked as ambiguous rather than summed.

Neither track performs inference, smoothing, interpolation, forecasting or
outbreak detection. Their purpose is to align already computed results while
keeping their provenance inspectable in the SVG.

Generate the integrated sheet with:

```bash
cargo run --example evolutionary_surveillance -- assets
```

## Scope

The renderer is for orthogonal, diagonal and curved rectangular trees,
circular and fan layouts, and equal-angle unrooted trees. It does not infer
trees, fit clocks or population models, reconstruct ancestral states, estimate
selection, smooth surveillance observations or claim epidemiological
transmission. [`PhyloMap`](maps.md#put-a-phylogeny-around-the-map) can place
terminal annotations at explicitly supplied coordinates, but it does not infer
those locations or the movement between them. Those analyses belong upstream;
Karyon preserves their topology, lengths, support, intervals and annotations
and makes the chosen encodings explicit. Nexus support is intentionally
limited to the portable first-tree and translation-table subset.
