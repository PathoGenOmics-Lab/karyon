# Figure

[`plot()`](plot.md) is the short way to write a figure down. `Figure` is what it
builds, and it is the layer to reach for when a track comes from an alternative
constructor, when a track has to be read back before it is drawn, or when the
stack is assembled somewhere other than a single expression. This page is that
layer, how it lays a stack out, and the three things around it: a `Panels`
sheet, the circular `Rings`, and the `Genome` that puts several sequences on one
axis.

```rust
use karyon::{AxisTrack, CoverageTrack, Figure, Region};

Figure::new(Region::parse("chr7:5,220,001-5,222,000")?)
    .title("A locus")
    .push(CoverageTrack::new(5_220_000, depth).label("depth"))
    .push(AxisTrack::new())
    .save_svg("locus.svg")?;
```

`Region::parse` reads the 1-based inclusive locus string samtools and IGV use,
so `chr7:5,220,001-5,222,000` becomes the 0-based half-open span
`[5220000, 5222000)`. That is why the coverage array starts at `5_220_000` and
not at `5_220_001`. Every other number in the API is 0-based; see
[Coordinates](../how-it-works/coordinates.md).

`Plot::into_figure` hands a `Figure` back from the short form, which is how a
plot ends up on a [`Panels`](#panels-a-sheet-of-figures) sheet.

## Building a stack

`push` appends a track below the ones already there, so **the order of the
`push` calls is the order of the stack**. Nothing about a track depends on where
it sits, so reordering a figure is reordering its calls.

`push_boxed` takes a `Box<dyn Track>` instead, which is what a stack built in a
loop or behind a condition needs: every arm of the loop produces a different
concrete track type, and a `Vec<Box<dyn Track>>` is the one thing they have in
common.

```rust
use karyon::{AxisTrack, CoverageTrack, Figure, Region, Track};

let mut tracks: Vec<Box<dyn Track>> = Vec::new();
for (name, depth) in samples {
    tracks.push(Box::new(CoverageTrack::new(0, depth).label(name)));
}
tracks.push(Box::new(AxisTrack::new()));

let mut figure = Figure::new(Region::new("chr2", 0, 4_000)?);
for track in tracks {
    figure = figure.push_boxed(track);
}
```

## How the layout works

A track knows two things: how tall it wants to be, and how to draw itself. It
never knows where it sits. The figure works out the rest.

```text
+------------------------------------------------------------------+
|  A locus                                   chr7:5220001-5222000  |  header
|                                                                  |
|       depth |  40 |    ____        ___                           |  band 0
|             |   0 |___/    \______/   \____                      |
|                                                                  |  track_gap
|       genes |     |   >>>>>>>>>>       <<<<<<<                   |  band 1
|                                                                  |  track_gap
|             |     |  |      |      |      |                      |  band 2
+------------------------------------------------------------------+
   \__________/\____/\____________________________________________/
   label gutter axis                 plotting area
```

The box is the image; the space between it and the content is `margin`.

The pieces, in the order the figure decides them:

**The header.** Present when there is a title, a region label, or both, and
`title_font_size + 12` pixels tall. The title is drawn bold at the left margin,
the locus string in the muted colour at the right. `show_region_label(false)`
turns the locus off; a figure with neither has no header at all.

**The label gutter.** Reserved only when at least one track returns a label from
[`Track::label`](../how-it-works/extending.md). A figure of unlabelled tracks
uses the full width whatever `label_width` says, so the setting costs nothing
until it is used. Labels are drawn right-aligned against the left edge of the
value axis strip, which is what makes a track with an axis and a track without
one line their names up. By default the gutter follows the widest label, between
48 and 160 pixels; `label_width` is the explicit override for fixed layouts.

**The value axis strip.** Each track is asked, through `Track::y_axis_width`,
how much room it wants for its own ticks. The figure then gives **every track
the widest of those requests**, rather than sizing each strip to its own track,
because two tracks in one figure must still start their plotting areas at the
same x. A coverage track asks for enough to print its ceiling; a pileup or a
sequence has no value to put a number on and asks for nothing, which is the
default.

**The bands.** Each track gets a rectangle as wide as the plotting area, as tall
as it asked for, stacked top to bottom with `track_gap` between them. The band
and its axis strip are wrapped in a clip, so a track that miscalculates cannot
paint into its neighbour.

**The height.** It is computed rather than configured. It is
`margin.top + header + sum of the track heights + track_gap * (n - 1) + margin.bottom`,
which is why `dimensions()` has to lay the figure out to answer.

The heights come from `Track::height(&scale)`, and the scale is passed in
because some tracks change height with zoom: a feature track needs more rows
when features overlap on screen, so the same track is shorter zoomed in than
zoomed out.

!!! note "The one track that ignores the shared axis"
    `IdeogramTrack` draws the whole chromosome across the plotting area and
    marks the region on it. Every other track maps its data through the same
    [`Scale`](../how-it-works/scale.md), which is what keeps the x axes aligned.

## Size and spacing

Every setting is a consuming builder method, so they chain in any order.

| Method | Default | What it does |
|:-------|:--------|:-------------|
| `width(f64)` | `900.0` | Image width in pixels. |
| `title(impl Into<String>)` | none | Bold line above the tracks. |
| `theme(Theme)` | `Theme::light()` | Colours, fonts and corner radius. See [Theming](theming.md). |
| `profile(RenderProfile)` | `Manuscript` geometry | Selects a named palette, scale and density together. |
| `visual_scale(f64)` | `1.0` | Scales type, marks, bands, margins, gaps and corners without changing genomic coordinates. |
| `density(Density)` | `Balanced` | Packs repeated data rows compactly or spaciously. |
| `margin(Margin)` | `14, 18, 14, 16` | Whitespace around the plotting area. |
| `label_width(f64)` | automatic | Explicit width of the left gutter holding track labels. |
| `track_gap(f64)` | `12.0` | Vertical gap between bands. |
| `show_region_label(bool)` | `true` | The locus string in the top right corner. |

```rust
use karyon::{Figure, Margin, Region, Theme};

let figure = Figure::new(Region::new("contig_01", 0, 50_000)?)
    .width(1_200.0)
    .theme(Theme::dark())
    .visual_scale(1.15)
    .margin(Margin { top: 14.0, right: 20.0, bottom: 14.0, left: 16.0 })
    .label_width(120.0)
    .track_gap(6.0)
    .show_region_label(false);
```

`Margin` has four public `f64` fields, `top`, `right`, `bottom` and `left`, and
a `Default` of `14, 18, 14, 16`. `top` is the space above the title, `left` the
space to the left of the track labels.

A width that would leave no plotting area is raised to the smallest one that
does, `margin.left + margin.right + label_width + 50`, and a width that is not
finite is replaced by the same floor. A figure is always renderable, and `NaN`
never reaches the output.

## Reading a figure back

```rust
let (width, height) = figure.dimensions();
println!("{width} x {height}, {} tracks", figure.track_count());
println!("{}", figure.region());
```

`dimensions()` is the one that does work: the height follows from the tracks, so
asking for it means laying the figure out. `region()` prints back in the 1-based
inclusive form, the same string the region label shows.

## Rendering

`to_svg()` returns a standalone SVG 1.1 document as a `String`. `save_svg(path)`
renders and writes it, returning whatever `fs::write` returns. The `save_svg`
methods on `Figure`, `Panels` and `Rings` are the only I/O the library does.

`to_svg_with_id_prefix(prefix)` is the same document with every id it generates
carrying `prefix`. Ids in SVG belong to the whole document rather than to the
element carrying them, and `url(#id)` resolves to the first match anywhere in
it. Two figures nested into one document would both claim `karyon-clip-0`, and
the second one's clips would silently resolve to the first one's rectangles,
cropping its tracks to bands belonging to a different figure. `Panels` hands
every figure a prefix of its own; do it yourself only if you assemble a sheet by
hand.

```rust
let left = figure.to_svg_with_id_prefix("left-");
```

## Panels: a sheet of figures

A paper figure is usually several plots with letters on them. `Panels` stacks
finished drawings into one SVG without any of them having to know about the
others.

![Every kind of plot karyon draws, on one sheet of twenty-two panels in three columns: a genomic stack, a read pileup, sequence logos, association statistics with a genotype matrix, a dotplot and synteny ribbons, a multiple sequence alignment, variable sites with a phylogeny, a tree, windowed statistics read against a baseline, a circular chromosome, raw nanopore signal, one locus compared across three genomes, Dam methylation across the E. coli origin of replication, an association scan across a whole draft assembly, structural variants as arcs between their breakpoints, the six reading frames, two trees face to face, a human imprinting control region read one molecule at a time, a coding sequence ruled in codons, one molecule aligned in three pieces, SARS-CoV-2 lineage deletions painted onto a phylogeny, and transcription units from start site to terminator](../assets/figures/gallery.svg)

```rust
use karyon::Panels;

let sheet = Panels::new()
    .title("Figure 1")
    .columns(2)
    .gap(20.0)
    .push_captioned(&stack, "A", "Read depth and annotation across the locus")
    .push_captioned(&circle, "B", "The same sample as a closed sequence");

sheet.save_svg("figure-1.svg")?;
```

There are three ways to add one: `push` for a figure with a letter,
`push_captioned` for a letter and a caption underneath, and `push_bare` for
neither. All three render the figure **now** and store the result, so later
changes to it do not reach the sheet.

| Method | Default | What it does |
|:-------|:--------|:-------------|
| `gap(f64)` | `18.0` | Vertical gap between panels. |
| `columns(usize)` | `1` | How many columns to lay the panels out in. |
| `row_major()` | off | Fills a comparison grid left-to-right before starting the next row. |
| `column_gap(f64)` | `26.0` | Horizontal gap between columns. |
| `margin(f64)` | `14.0` | Whitespace around the whole sheet. |
| `theme(Theme)` | `Theme::light()` | Used for the sheet's own title and letters, not for the panels, each of which keeps its own. |
| `profile(RenderProfile)` | manuscript scale | Selects a named sheet palette and scale. |
| `visual_scale(f64)` | `1.0` | Scales nested drawings and all sheet chrome together. |
| `align_plot_areas(bool)` | `true` | Aligns the data origins of linear figures despite different gutters. |
| `title(impl Into<String>)` | none | A line across the top of the sheet. |

`len()`, `is_empty()` and `dimensions()` read the sheet back; `to_svg()` and
`save_svg(path)` render it.

Panels fill each column top to bottom before starting the next, so the letters
still read in order: a reader following A, B, C down a column and back up to the
top of the next is doing what a two-column paper makes them do anyway. Where the
columns break is chosen by weighing every possible cut rather than by giving
each column an equal count of unequal panels: the tallest column is what the
sheet is as tall as, so it is minimised first, and the sum of squares breaks the
tie so the remaining columns come out level rather than ragged.

For a dashboard-like comparison where A and B belong on the first row and C
and D on the second, add `row_major()`. It keeps that left-to-right reading
order and gives every row the height of its tallest panel. Use
`align_plot_areas(false)` when the grid deliberately mixes unrelated coordinate
systems, such as a rectangular tree, a radial tree and an alignment.

The letters get a strip of their own on the left rather than being drawn over
the panels. A letter on top of a panel is a letter on top of data, and on top of
an opaque page colour it is a letter that is not there at all.

Each panel goes into the sheet as a nested `<svg>` inside a translated group.
Nothing is reparsed and nothing is rewritten, so a panel on a sheet is the same
picture as the panel on its own. The one thing a panel does not keep to itself
is its ids, which is why every one is rendered with a prefix.

## Rings: a circular sequence

A plasmid, an organelle genome, a viral genome and most bacterial chromosomes
have no ends. Drawing one as a line puts an edge where the biology has none,
straight through whatever happens to sit at coordinate zero. `Rings` maps
position to an angle instead, which is a different coordinate system and
therefore a different container: `Rings` is to `Ring` what `Figure` is to
`Track`.

![A 4.41 Mb chromosome drawn as concentric rings, its name and length in the middle: a position ruler around the outside, a dense ring of genes with the forward strand outside the reverse, three named loci, point mutations as tick marks coloured by category, then GC content and GC skew as rings rising above and falling below their baselines, and two translucent ribbons across the middle joining the ends of a rearrangement](../assets/figures/example-circular.svg)

```rust
use karyon::{AxisRing, FeatureRing, MarkerRing, Rings, SignalRing};

Rings::new(154_478)
    .title("chloroplast")
    .subtitle("154.5 kb")
    .push(AxisRing::new())
    .push(FeatureRing::new(genes).thickness(18.0))
    .push(MarkerRing::categorised(sites).thickness(12.0))
    .push(SignalRing::new(composition).thickness(40.0))
    .link((10_000, 34_000), (110_000, 134_000))
    .save_svg("plastome.svg")?;
```

Position zero is at twelve o'clock and coordinates run clockwise, which is the
convention every circular genome viewer uses. Rings are added outside in: the
first `push` takes the outermost radius, and each one after it sits inside the
last.

Four ring types ship, and they are the circular counterparts of what goes in a
band: `AxisRing` is a ruler of positions, `FeatureRing` draws annotation as arcs
with the forward strand outside and the reverse inside, `SignalRing` draws a
quantity in windows either side of a baseline circle, and `MarkerRing` draws
points as radial ticks, optionally categorised so each category gets its own
colour.

What a circle adds over a stack of bands is the middle. `link` draws a
translucent ribbon between two spans that belong together, the two ends of an
inversion, a duplication and its source, a pair of sequences a rearrangement
joined. Both spans are `(start, end)`, 0-based and half-open. `link_colored`
takes a colour and an opacity of its own. Chords are drawn before the rings, so
a dozen ribbons crossing the middle do not wash the data out.

| Method | Default | What it does |
|:-------|:--------|:-------------|
| `diameter(f64)` | `640.0` | Diameter of the outermost ring, in pixels. |
| `margin(f64)` | `14.0` | Whitespace around the circle. |
| `origin_gap(f64)` | `2.0` | Degrees left blank at twelve o'clock. |
| `theme(Theme)` | `Theme::light()` | Colours and fonts. |
| `profile(RenderProfile)` | manuscript scale | Selects palette, type scale and ring density together. |
| `visual_scale(f64)` | `1.0` | Scales typography, marks, margins and ring geometry. |
| `density(Density)` | `Balanced` | Packs concentric bands compactly or spaciously. |
| `title(impl Into<String>)` | none | Name written in the middle. |
| `subtitle(impl Into<String>)` | none | A second, quieter line under it. |

The default two-degree origin gap is there on purpose. A closed circle hides the
fact that a coordinate system has to start somewhere and that the choice was
arbitrary; a seam says so. `origin_gap(0.0)` closes it.

`length()`, `ring_count()`, `dimensions()` and `inner_radius()` read the plot
back, the last being the radius where chords start. `push_boxed`, `to_svg`,
`to_svg_with_id_prefix` and `save_svg` work exactly as they do on `Figure`.

## Genome: several sequences end to end

A figure is one region on one sequence, which is right for a locus and wrong for
an assembly. Two hundred contigs, a couple of dozen chromosomes, or a
chromosome and its plasmid do not fit in one region. `Genome` lays the sequences
end to end and hands back the single region that covers them all, and from then
on every track in the crate works across all of them at once.

![An association scan across a whole draft assembly: twelve contigs laid end to end on one axis, association points banded by contig with one peak crossing the genome-wide line, a depth profile that falls away at every contig end, and a track naming the contigs underneath](../assets/figures/example-genomewide.svg)

```rust
use karyon::{Figure, Genome, GenomeTrack};

let genome = Genome::new([
    ("chrI", 230_218u64),
    ("chrII", 813_184),
    ("chrIII", 316_620),
]);

assert_eq!(genome.total(), 1_360_022);

// A position on one sequence becomes a position on the shared axis.
let at = genome.at("chrII", 1_000).unwrap();
assert_eq!(at, 231_218);
assert_eq!(genome.locate(at), Some(("chrII", 1_000)));

Figure::new(genome.region())
    .push(GenomeTrack::new(genome.clone()).label("chromosomes"))
    .save_svg("genome.svg")?;
```

The order given is the drawing order and is not sorted: a reference is
conventionally ordered by chromosome number and an assembly by descending contig
length, and neither is something to guess at.

| Method | What it gives back |
|:-------|:-------------------|
| `total()` | Length of the whole axis, gaps included. |
| `region()` | The region covering every sequence, named `genome`. |
| `sequences()` | The `Chromosome` list, each a `name` and a `length`. |
| `len()`, `is_empty()` | How many sequences there are. |
| `offset(name)` | Where a sequence starts on the shared axis, or `None`. |
| `at(name, pos)` | A per-sequence position as a shared-axis position. |
| `locate(pos)` | The inverse: the sequence and offset a position falls on. |
| `boundaries()` | Where each sequence starts, for `ManhattanTrack::bands`. |
| `spans()` | Each sequence as its name and both ends of its span. |
| `map(items)` | Maps `(name, pos, value)` triples onto the axis. |

`at` returns `None` when the genome has no such sequence, which is what a name
mismatch between a variant file and a reference looks like. `map` does the same
for a whole list and reports how many it dropped, so a caller can say how many
went missing rather than wondering why the figure came out thin.

```rust
let (mapped, dropped) = genome.map(hits);
if dropped > 0 {
    eprintln!("{dropped} points named a sequence the genome does not have");
}
```

`Genome::checked` refuses a repeated name with `Error::DuplicateSequence`, since
`at` would otherwise place everything silently on the first of the two.
`gap(bases)` puts blank axis between one sequence and the next, which makes the
joins obvious without any drawing, at the cost of coordinates belonging to no
sequence: `locate` returns `None` inside a gap, which is the honest answer.

!!! warning "A concatenated axis is not a ruler"
    A distance measured across a boundary is not a distance. Two points a pixel
    apart may be the last base of one contig and the first of the next, which
    are not neighbours in any sense that matters. `GenomeTrack` draws where the
    joins are, and it labels each sequence rather than printing global
    coordinates, because a ruler of the concatenated coordinate system would be
    a ruler of a system nothing else uses.

## Drawing: what a sheet asks for

`Figure` and `Rings` have nothing in common. One stacks bands over a horizontal
scale, the other wraps rings around an angle. What lets them share a sheet is
that both can state a size and render themselves at it, and that is the whole of
the `Drawing` trait:

```rust
pub trait Drawing {
    fn dimensions(&self) -> (f64, f64);
    fn to_svg_with_id_prefix(&self, prefix: &str) -> String;
}
```

`Panels::push` and its two siblings take `&impl Drawing`, so anything that
implements those two methods goes on a sheet, including a type of your own:

```rust
use karyon::{Drawing, Panels};

fn sheet(a: &impl Drawing, b: &impl Drawing) -> Panels {
    Panels::new().push(a, "A").push(b, "B")
}
```

Panel J of the gallery above is a `Rings` plot and every other panel is a
`Figure`. The sheet does not know or care which is which.

## Next

- [Plot API](plot.md), for the short form that builds one of these.
- [Theming](theming.md), for the colours and fonts a figure hands every track.
- [Writing a track](../how-it-works/extending.md), for what a figure asks of a
  track and what it does on the track's behalf.
