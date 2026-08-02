# Plot API

`plot()` is the short way to write a figure down: one call per track, in the
order they stack, and the figure comes out the other end. What it remembers
between calls, what it fills in without being asked, and where the short form
stops are what this page covers.

A figure is a stack of tracks over one shared coordinate axis. Built with
[`Figure`](figure.md) directly, a good deal of what gets typed is plumbing: a
`Region` to unwrap, a `use` line naming every track type, a `push` around each
track, and the window start repeated on every track that takes one even though
the figure already holds it. `plot()` holds the region and names the track
types, so none of that is written twice.

![A coverage profile with a dropout, the reference sequence, two gene models and variants coloured by consequence, all over one coordinate axis](../assets/figures/example.svg)

```rust
use karyon::{plot, Aggregate, Feature, Strand, Variant};

plot("NC_000962.3:761000-762999")?
    .title("rpoB locus, resistance determining region")
    .add_coverage(depth)
    .label("depth")
    .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
    .add_sequence(bases)
    .label("reference")
    .add_features(vec![
        Feature::new(759_806, 763_325).name("rpoB").strand(Strand::Forward),
        Feature::new(761_081, 761_162).name("RRDR").strand(Strand::Forward),
    ])
    .label("annotation")
    .add_variants(vec![
        Variant::new(761_108).value(0.98).category("missense"),
        Variant::new(761_154).value(1.00).category("missense"),
        Variant::new(761_155).value(0.21).category("synonymous"),
    ])
    .label("variants")
    .save("example.svg")?;
```

Three things in that chain are worth naming now. The coordinate ruler along the
bottom was never asked for: a plot fills it in. `add_coverage` takes the start
of its array from the region, so the window start is written once rather than on
every track that needs it. And `label` and `adjust` talk about the track the
line above them added, not about the figure.

The figure above comes from `cargo run --example locus`, which generates the
depth and the bases, colours the RRDR box, and carries two more variants.

## Starting a plot

`plot(locus)` parses a locus string and hands back a `Plot` over that region.

```rust
use karyon::plot;

let figure = plot("chr7:140,753,000-140,754,000")?;
```

The string is the 1-based inclusive form samtools and IGV use. Thousands
separators are ignored, and the sequence name may itself contain colons, since
the split happens at the last one. Anything that does not parse is an
`Error::InvalidLocus` rather than a panic: no colon, no dash after it,
coordinates that are not numbers, a start below 1, or an end before the start.
The error converts into `io::Error`, so the region and the file it renders to
share one `?`.

!!! warning "Coordinates"
    The locus string is 1-based and inclusive. Everything else in the API is
    **0-based and half-open**, the BED convention, so a VCF `POS` or a GFF
    `start` is `pos - 1` on the way in. `plot("chr1:101-200")` gives a region
    whose `start()` is `100` and whose `end()` is `200`. Tick labels are printed
    back in the 1-based form. See [Coordinates](../how-it-works/coordinates.md).

`Plot::over` is the same thing from a `Region` already in hand, which is the
form to use when the coordinates are computed rather than typed:

```rust
use karyon::{Plot, Region};

let region = Region::new("NC_000962.3", 761_000, 763_000)?;
let figure = Plot::over(region);
```

`region()` reads the region back at any point in the chain.

## The track a plot remembers

A plot holds the track the last `add_` built rather than handing it straight to
the figure, so the call after it can still reach it. `label` names it for the
left gutter and `adjust` hands it to a closure for everything else:

```rust
use karyon::{plot, Aggregate, CoverageStyle};

plot("chr7:140,753,000-140,754,000")?
    .add_coverage(depth)
    .label("depth")
    .adjust(|track| {
        track
            .aggregate(Aggregate::Min)
            .style(CoverageStyle::Area)
            .height(70.0)
    })
```

The track is held in the plot's type parameter: `add_coverage` returns a
`Plot<CoverageTrack>`, `add_features` a `Plot<FeatureTrack>`, and so on down the
chain. That is what lets the closure be handed the concrete type, with every
builder method the track has in reach and no downcasting anywhere. A method that
is not on that track is a compile error rather than something ignored at render
time:

```rust
plot("chr7:140,753,000-140,754,000")?
    .add_coverage(depth)
    .add_features(genes)
    .adjust(|track| track.aggregate(Aggregate::Min))
```

```text
error[E0599]: no method named `aggregate` found for struct `FeatureTrack` in the current scope
  |
  |     .adjust(|track| track.aggregate(Aggregate::Min))
  |                           ^^^^^^^^^ method not found in `FeatureTrack`
```

A fresh plot is a `Plot<Empty>`, and `Empty` is not a track, so `label` and
`adjust` are not there to be called before anything has been added:

```text
error[E0599]: the method `label` exists for struct `Plot`, but its trait bounds were not satisfied
  |
  |     plot("chr1:1-1000")?.label("depth")
  |                          ^^^^^ method cannot be called on `Plot` due to unsatisfied trait bounds
  |
  = note: the following trait bounds were not satisfied:
          `karyon::plot::Empty: karyon::plot::Named`
```

### What is not checked

**Which track a call reaches is a matter of where the call sits, and that is not
checked.** An `adjust` written one `add_` too late configures the next track
instead, and compiles:

```rust
plot("chr7:140,753,000-140,754,000")?
    .add_coverage(depth)
    .label("depth")
    .add_sequence(bases)
    .adjust(|track| track.height(70.0))   // the sequence track, not the depth
```

`height` is a builder on many of the track types, so the closure type-checks
against either one and the figure renders with the wrong band resized. `label`
is worse: every track has it, which is exactly what makes `Plot::label` work
without being told the type, so a misplaced `label` always compiles and always
names the wrong row. Nothing in the type system can catch either. Read a chain
by its `add_` calls and treat what follows each one as belonging to it.

### What puts the track away

The pending track reaches the figure when anything else needs the plot to move
on. The next `add_` puts it away and holds the new one instead. `add_track`,
`add_boxed`, `done` and `save` hand back a plot with nothing pending, so `label`
and `adjust` are again not callable, which is the same compile error as on a
fresh plot. `to_svg` and `into_figure` consume the plot altogether.

!!! note "Passing a plot around"
    `Plot` on its own means `Plot<Empty>` and will not take a plot mid-chain,
    since the type changes with every `add_`. A function that wants one is
    generic over `T: Slot`, or takes the `Figure` that `into_figure` hands back.
    `Slot` is sealed: a track type from outside the crate goes in through
    `add_track`, not by implementing it.

## Every `add_` method

One per track type the crate ships. Each returns a plot holding that track.

| Method | Builds | Takes |
|:-------|:-------|:------|
| `add_axis()` | `AxisTrack` | nothing; the coordinate ruler |
| `add_bisulfite(sites, molecules)` | `BisulfiteTrack` | `Vec<u64>`, `Vec<Molecule>` |
| `add_clades(tree, blocks)` | `CladeTrack` | `Tree`, `Vec<CladeBlock>` |
| `add_codons(start, end, strand)` | `CodonTrack` | `u64`, `u64`, `Strand` |
| `add_coverage(values)` | `CoverageTrack` | `Vec<f64>`, one per base |
| `add_dotplot(blocks)` | `DotplotTrack` | `Vec<AlignmentBlock>` |
| `add_features(features)` | `FeatureTrack` | `Vec<Feature>` |
| `add_genome(genome)` | `GenomeTrack` | `Genome` |
| `add_ideogram(length, bands)` | `IdeogramTrack` | `u64`, `Vec<Band>` |
| `add_legend(legend)` | `LegendTrack` | `Legend` |
| `add_loci(loci)` | `LocusTrack` | `Vec<Locus>` |
| `add_logo(columns)` | `LogoTrack` | `Vec<LogoColumn>` |
| `add_manhattan(points)` | `ManhattanTrack` | `Vec<Association>` |
| `add_matrix(sites, rows)` | `MatrixTrack` | `Vec<u64>`, `Vec<MatrixRow>` |
| `add_methylation(sites)` | `MethylationTrack` | `Vec<MethylSite>` |
| `add_msa(sequences)` | `MsaTrack` | `Vec<MsaSequence>` |
| `add_orfs(seq)` | `OrfTrack` | `Vec<u8>` |
| `add_pileup(reads)` | `PileupTrack` | `Vec<Read>` |
| `add_sequence(seq)` | `SequenceTrack` | `Vec<u8>` |
| `add_snps(names, sites)` | `SnpTrack` | anything iterable of names, `Vec<SnpSite>` |
| `add_split_reads(reads)` | `SplitReadTrack` | `Vec<SplitRead>` |
| `add_squiggle(signal)` | `SquiggleTrack` | `Vec<f64>`, one per sample |
| `add_structural(variants)` | `StructuralTrack` | `Vec<StructuralVariant>` |
| `add_synteny(blocks)` | `SyntenyTrack` | `Vec<AlignmentBlock>` |
| `add_tanglegram(left, right)` | `TanglegramTrack` | `Tree`, `Tree` |
| `add_transcription_units(units)` | `TranscriptionUnitTrack` | `Vec<TranscriptionUnit>` |
| `add_tree(tree)` | `TreeTrack` | `Tree` |
| `add_variants(variants)` | `VariantTrack` | `Vec<Variant>` |
| `add_windows(windows)` | `WindowTrack` | `Vec<Window>` |

Five of them have an `_at` form as well: `add_coverage_at`, `add_logo_at`,
`add_orfs_at`, `add_sequence_at` and `add_squiggle_at`. Those are the tracks
that lay an array along the axis, and the extra argument is where the first
element sits. See [what a plot fills in](#what-a-plot-fills-in) below.

Every argument that says `Vec<T>` is really `impl Into<Vec<T>>`, so an array or
a slice's `to_vec()` goes in unchanged. What each track draws and what its
builder methods do is on the [Tracks](../tracks.md) page.

## What a plot fills in

Two things happen without being asked. Both can be undone, and neither is a
default that has to be lived with.

### The axis at the bottom

A figure with no coordinates along it is rarely what anyone meant, so an
`AxisTrack` is appended at the bottom when the plot is turned into a figure.

`add_axis()` puts the ruler somewhere else instead, and the automatic one is
then not appended, so there is still exactly one. This is how a ruler goes
between two tracks, or above them:

```rust
plot("chr1:1-1000")?
    .add_axis()
    .add_coverage(depth)
    .label("depth")
```

`remove_axis()` leaves it out altogether, which is what a tanglegram or a bare
tree wants: those have leaf rows, not coordinates, and a ruler under them would
be measuring nothing.

![Two trees of the same eight isolates face to face, core genome on the left and accessory genome on the right, tips joined by lines across the middle and the crossings coloured](../assets/figures/example-tanglegram.svg)

```rust
plot("taxa:1-8")?
    .title("K. pneumoniae: core and accessory genome trees over one collection")
    .remove_region_label()
    .remove_axis()
    .add_track(track)
```

`remove_axis` suppresses the automatic ruler and only that. An axis put
somewhere with `add_axis`, or handed over with `add_track`, stays where it was
put. It can be called from either side of the tracks, since it is a decision
about the finished figure rather than a position in the stack.

!!! warning "`add_track` cannot tell a ruler from anything else"
    `add_track(AxisTrack::new())` adds a ruler **in addition** to the automatic
    one, and the figure gets two. `add_axis()` is the call that knows about the
    one at the bottom.

### The region start for tracks that carry an array

A coverage profile, a reference sequence, a set of reading frames, a logo and a
raw signal are all arrays laid along the axis, and each needs to know where its
first element sits. Almost always that is the left edge of the region, which the
plot already holds:

```rust
plot("NC_000962.3:761001-762000")?
    .add_coverage(depth)      // depth[0] is at 0-based 761,000
```

When the array starts somewhere else, the `_at` form takes the start
explicitly. A profile covering the second half of the window:

```rust
plot("NC_000962.3:761001-762000")?
    .add_coverage_at(761_500, depth)
```

The two forms build the same track: `add_coverage(values)` is
`add_coverage_at(region.start(), values)`, and the tests pin an implicit start
against the same start written out as identical renderings. For a squiggle the
start is a sample index rather than a base, so `add_squiggle_at` takes a `usize`
and `add_squiggle` casts the region start to one.

### What the axis counts

An alignment is indexed by column and a raw signal by sample, so a plot holding
one of those wants a region in that unit: `plot("alignment:1-320")`. The region
is a coordinate system, not a claim about a genome. What the ruler counts is
whatever the data is indexed by.

![A conservation logo above a multiple sequence alignment of twelve isolates, with only the disagreements painted](../assets/figures/example-msa.svg)

The locus string printed at the top right is the one place that still reads as
genomic, so `remove_region_label()` leaves it out when it would mislead:

```rust
plot(&format!("alignment:1-{columns}"))?
    .title("An alignment, and what disagrees in it")
    .remove_region_label()
```

## Figure settings

These are the figure's own settings rather than any track's. They return the
plot unchanged in every other respect, including its pending track, so they can
sit anywhere in the chain.

| Method | Default | Sets |
|:-------|:--------|:-----|
| `width(f64)` | `900.0` | image width in pixels |
| `title(impl Into<String>)` | none | the title drawn above the tracks |
| `theme(Theme)` | `Theme::light()` | colours and fonts |
| `margin(Margin)` | `10 / 16 / 10 / 12` | whitespace top, right, bottom, left |
| `label_width(f64)` | `84.0` | width of the left gutter holding track labels |
| `track_gap(f64)` | `10.0` | vertical gap between tracks |
| `remove_region_label()` | shown | drops the locus string at the top right |

The height is not in the table because it is not configurable: it follows from
the tracks, each of which reports how tall it wants to be for the scale in
force. A width too small to leave a plotting area is raised to the smallest one
that does, so a figure always renders. The label gutter is only reserved when at
least one track carries a label, so `label_width` on a figure of unlabelled
tracks changes nothing.

`Theme::dark()` is the other theme that ships, and it is a selected set of
colours rather than an inversion of the light one. See
[Theming](theming.md).

## Rendering

`to_svg()` returns a standalone SVG document as a `String`. `save(path)` writes
it and hands the plot back, which is what lets one stack be rendered twice:

```rust
use karyon::{plot, Theme};

plot("plasmid:1-2000")?
    .add_coverage(depth)
    .save("example.svg")?
    .theme(Theme::dark())
    .save("example-dark.svg")?;
```

A save closes the stack. It puts the pending track away and fills in the axis,
so saving twice does not draw two rulers, and a track added after a save sits
**below** the axis rather than above it.

`into_figure()` hands over the `Figure`, with the pending track pushed and the
axis in place. That is the way out to anything that wants a figure:
`Figure::dimensions` to find out how tall the thing came out,
`Figure::to_svg_with_id_prefix` to nest it in another document by hand, or
`Panels` to put it on a sheet beside other figures, including a circular `Rings`
plot.

```rust
use karyon::{plot, Panels};

let top = plot("chr1:1-1000")?.add_coverage(depth).into_figure();
let bottom = plot("chr1:1-1000")?.add_manhattan(points).into_figure();

Panels::new().push(&top, "A").push(&bottom, "B").save_svg("sheet.svg")?;
```

## Where it stops

A plot is a stack of tracks over one axis and nothing more. Three situations ask
for more than the `add_` methods can express, and all three are answered inside
the plot rather than by abandoning it.

### A track from an alternative constructor

Several track types have constructors that compute the track from raw input:
`WindowTrack::gc_skew`, `WindowTrack::ratios`, `SnpTrack::from_alignment`,
`LogoTrack::from_sequences`, `IdeogramTrack::bare`. There is no `add_` for each
of those, because the list would be as long as the constructors are. `add_track`
takes a finished track and pushes it straight through:

![Two windowed statistics over the same window: pN/pS coloured by which side of the baseline each window falls on, and GC skew as a line that turns over a third of the way along](../assets/figures/example-selection.svg)

```rust
use karyon::{Plot, Region, Window, WindowStyle, WindowTrack};

Plot::over(Region::new("NC_000962.3", start, start + 40_000)?)
    .title("Selection and strand composition, read against their baselines")
    .add_track(WindowTrack::ratios(windows).label("pN/pS").height(70.0))
    .add_track(
        WindowTrack::gc_skew(start, &bases, 1_000)
            .style(WindowStyle::Line)
            .label("GC skew"),
    )
```

Name the track before handing it over. `add_track` leaves nothing pending, so
`Plot::label` and `Plot::adjust` do not reach it afterwards. The track's own
`label` builder is the one to use, and it is the same method `Plot::label` ends
up calling.

`add_boxed` is the same thing for a `Box<dyn Track>`, which is what a stack
assembled at runtime from a list of choices ends up holding. Both take any
implementation of the [`Track`](../how-it-works/extending.md) trait, including
one from outside the crate.

### A track read back before it is drawn

Some tracks answer questions about their own data, and the answer is often what
the rest of the figure is positioned against. `CodonTrack::span_of` gives the
reference span of a codon, `SnpTrack::sites` lists the variable sites the panel
kept, `TanglegramTrack::crossings` counts the crossing ties in the drawing.
Reading any of those needs the track as a variable, so it is built first and
added with `add_track`:

![A coding sequence drawn as numbered codons with their translated residues, two variant lollipops sitting over the codons they change, and a base ruler underneath](../assets/figures/example-codons.svg)

```rust
use karyon::{CodonTrack, Plot, Region, Strand, Variant};

let ruler = CodonTrack::new(759_806, 763_325, Strand::Forward)
    .sequence(view_from, bases)
    .label("rpoB");
let s450 = ruler.span_of(450).expect("rpoB has a codon 450");

Plot::over(Region::new("NC_000962.3", view_from, view_to)?)
    .add_variants(vec![Variant::new(s450.0 + 1).category("S450L")])
    .label("variants")
    .add_track(ruler)
```

The variant lands on the codon the ruler says it does, because both come from
the same object rather than from the same arithmetic done twice.

### A stack built in a loop

The plot's type changes with every `add_`, and every arm of a loop or of an
`if` has to have one type. `done()` puts the pending track away and gives back a
`Plot<Empty>`, which is a type that does not change:

```rust
let mut figure = plot("chr2:1-4,000")?;
for (name, depth) in samples {
    figure = figure.add_coverage(depth).label(name).done();
}
figure.save("samples.svg")?;
```

`done` renders nothing and decides nothing. It is only the point at which the
plot stops naming the track it is holding.

### The layer underneath

For anything past those three, [`Figure`](figure.md) takes tracks that are
already built and is what a plot produces anyway. Neither layer can draw
anything the other cannot: the same stack written both ways renders the same
document, byte for byte, and a test in `src/plot.rs` asserts it.

=== "With `plot`"

    ```rust
    plot("NC_000962.3:761001-762000")?
        .title("a locus")
        .add_coverage(depth)
        .label("depth")
        .add_features(genes)
        .label("genes")
        .to_svg()
    ```

=== "With `Figure`"

    ```rust
    Figure::new(Region::parse("NC_000962.3:761001-762000")?)
        .title("a locus")
        .push(CoverageTrack::new(761_000, depth).label("depth"))
        .push(FeatureTrack::new(genes).label("genes"))
        .push(AxisTrack::new())
        .to_svg()
    ```

The difference is only how much has to be written down.

## Next

- [Figure](figure.md), for the layer underneath and the containers a finished
  figure goes into.
- [Tracks](../tracks.md), for what each `add_` builds and what it refuses to do.
- [Command line](cli.md), for the same grammar with spaces instead of dots.
