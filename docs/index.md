---
hide:
  - navigation
  - toc
---

<div class="hero" markdown>

![The karyon logo: a marker pen drawing a DNA double helix, half of it still blue and half already redrawn in orange](assets/karyon_logo.svg){ .hero-logo }

# karyon

Genomic track plots for Rust. A stack of tracks over one shared coordinate
axis, rendered to standalone SVG.

</div>

![A coverage profile with a dropout, the reference sequence, two gene models and variants coloured by consequence, all over one coordinate axis](assets/figures/example.svg)

A general plotting library knows about points and lines. It does not know that a
position is a base, that a gene has a strand, or that a figure stops being worth
anything the moment its tracks stop lining up. `karyon` is the small amount of
code that does. Thirty-three track types, one shared coordinate axis, plain SVG
1.1, and no dependencies at all.

<div class="grid cards" markdown>

-   :material-download: **[Installation](getting-started/installation.md)**

    Point Cargo at the repository for the library, `cargo install` for the
    command. Nothing else to install alongside it.

-   :material-rocket-launch: **[Quickstart](getting-started/quickstart.md)**

    A first figure in a dozen lines, and the same figure without writing any
    Rust.

-   :material-view-grid-outline: **[Plot catalogue](plots/index.md)**

    Eight visual routes through all thirty-three tracks plus circular and
    geographic drawings. Choose by biological question or data shape.

-   :material-code-braces: **[Plot API](guide/plot.md)**

    One call per track, in the order they stack. What the plot remembers
    between calls, what it fills in, and where the short form stops.

-   :material-console: **[Command line](guide/cli.md)**

    The same grammar with spaces instead of dots. Flag order is stack order,
    and any track file may be `-`, so the pipeline is the parser.

-   :material-file-document-outline: **[File formats](guide/formats.md)**

    bedGraph, BED, GFF3, VCF, FASTA, Newick, cytoBand and SAM text, each read
    with the coordinate convention its own specification defines.

-   :material-family-tree: **[Annotated phylogenetics](guide/phylogenetics.md)**

    BEAST, NHX and Nexus metadata; dated trees, topology operations, branch
    colours, collapsed clades and aligned sample traits.

-   :material-map: **[Geographic genomics](guide/maps.md)**

    Occurrence maps, explicit geographic links and circular phylogenies around
    an offline, deterministic world map.

-   :material-ruler: **[Coordinates](how-it-works/coordinates.md)**

    0-based and half-open everywhere, with two deliberate exceptions, both of
    them places where a person reads the number.

-   :material-magnify-scan: **[Scale](how-it-works/scale.md)**

    What a track does when one pixel covers two thousand bases, and why a
    genome-wide figure is still a small file.

-   :material-book-open-variant: **[Recipes](recipes.md)**

    Short worked figures for the questions that come up repeatedly, each one a
    program that runs.

</div>

## The figure above, in full

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
        Feature::new(759_806, 763_325)
            .name("rpoB")
            .strand(Strand::Forward),
        Feature::new(761_081, 761_162)
            .name("RRDR")
            .strand(Strand::Forward)
            .color("#d55e00"),
    ])
    .label("annotation")
    .add_variants(vec![
        Variant::new(761_108).value(0.98).category("missense"),
        Variant::new(761_138).value(0.55).category("missense"),
        Variant::new(761_154).value(1.00).category("missense"),
        Variant::new(761_155).value(0.21).category("synonymous"),
        Variant::new(761_051).value(0.12).category("synonymous"),
    ])
    .label("variants")
    .save("example.svg")?;
```

The thing to notice is that nothing in there positions anything: the box over
the resistance determining region and the lollipops standing inside it line up
without being told to, and the dropout in the depth track falls exactly where
the depth values put it, because every track maps its data through the one scale
the region defines. The ruler along the bottom was added without being asked
for. This is `examples/locus.rs`, which renders the figure above.

!!! note "Two coordinate systems, on purpose"
    The locus string is the 1-based inclusive form samtools and IGV use, and so
    are the tick labels, because those are the two numbers a person reads.
    Everything else, `Feature::new` and `Variant::new` among them, is 0-based
    and half-open like BED, so a VCF `POS` goes in as `POS - 1`.
    [Coordinates](how-it-works/coordinates.md) has the whole of it.
