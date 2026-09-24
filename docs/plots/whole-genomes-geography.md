# Whole genomes and maps

Plots for context beyond one locus: where a window sits in its chromosome, a
signal across a whole assembly, a circular genome from end to end, and where
samples were taken.
{ .k-lead }

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| Where in the chromosome is the window I am showing? | [Ideogram](../tracks/whole-genome.md#ideogramtrack) | `IdeogramTrack`, `--ideogram` |
| How does a signal run across every contig of an assembly? | [Whole assembly](../tracks/whole-genome.md#genometrack) | `GenomeTrack` over a `Genome`, Rust only |
| What does a circular chromosome or plasmid look like all the way round? | [Circular genome](../guide/plot.md) | `Rings`, Rust only |
| Where were the samples taken, and what moved between places? | [Map](../guide/maps.md) | `Map`, Rust only |

A circular genome and a map are drawings rather than tracks: each has a
coordinate system of its own, and a `Panels` sheet puts them beside a figure.
For a tree around a map, see [phylogeny and clades](phylogeny-clades.md).

## Plots

<div class="k-plots" markdown>

-   [![A banded chromosome 7 with a marker at the sixty kilobases shown beneath it as a depth profile, two genes and variant calls](../assets/figures/example-ideogram.svg){ width="900" height="274" loading="lazy" }](../tracks/whole-genome.md#ideogramtrack)

    **[Ideogram](../tracks/whole-genome.md#ideogramtrack)**
    The whole chromosome across the plot with a marker on the region in view, banded from a UCSC cytoBand table.

-   [![The M. tuberculosis H37Rv chromosome as a bare outline with the rpoB window marked on it, above the rpoB gene](../assets/figures/example-ideogram-bacterial.svg){ width="760" height="140" loading="lazy" .k-wide }](../tracks/whole-genome.md#ideogramtrack)

    **[Bare ideogram](../tracks/whole-genome.md#ideogramtrack)**
    An outline for the sequences with no bands to draw, such as bacteria, plasmids, organelles, viruses and draft assemblies, which still shows where the window is.

-   [![An association scan and a depth profile across every contig of a draft assembly, with the contigs beneath as alternating named blocks and a count of the names that did not fit](../assets/figures/example-genomewide.svg){ width="940" height="288" loading="lazy" }](../tracks/whole-genome.md#genometrack)

    **[Whole assembly](../tracks/whole-genome.md#genometrack)**
    Every sequence of an assembly laid end to end as named blocks, so any track can run across all of them on one axis.

-   [![The 4.41 Mb H37Rv chromosome as concentric rings: a ruler, genes on both strands, three named resistance genes, resistance mutations coloured by drug, GC content and GC skew, with two chords across the middle joining the ends of rearrangements](../assets/figures/example-circular.svg){ width="688" height="688" loading="lazy" }](../guide/plot.md)

    **[Circular genome](../guide/plot.md)**
    Position as an angle, with rings for features, markers, signals and a ruler, and chords across the middle for rearrangements.

-   [![Synthetic sampling sites on an equirectangular world map and again in Mercator, sized by count and coded by region in colour and shape, and an orthographic hemisphere with directed introductions drawn as arrows](../assets/figures/example-maps.svg){ width="1398" height="972" loading="lazy" }](../guide/maps.md)

    **[Map](../guide/maps.md)**
    Locations and flows you supply under an equirectangular, Mercator or orthographic projection; a coordinate off the Earth is counted, never moved.

</div>

## Related

<div class="grid cards" markdown>

-   **[Maps](../guide/maps.md)**

    Projections, flows, locations that cannot be drawn, and a tree around a
    map.

-   **[The Rust API](../guide/plot.md)**

    `Genome`, `Rings` and `Panels`, and how a figure is put together.

-   **[Phylogeny and clades](phylogeny-clades.md)**

    A circular tree framing a map, with its tips linked to where they were
    sampled.

-   **[Comparisons and alignments](comparisons-alignments.md)**

    One genome against another, rather than one genome in context.

</div>
