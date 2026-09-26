# Variation and association

Plots for calls and the patterns they make across samples: point and
structural variants, copy number, variable sites, genotype matrices and
heatmaps, association scans, linkage and other pairs, and site-wise
selection.
{ .k-lead }

## How to choose

Choose by what carries the result: one position, two breakpoints, a fitted
segment, a variable column, a cell for each sample and site, or one test per
position.

| Your question | Plot | Build it with |
|:--|:--|:--|
| Where are the calls, and how strong is each one? | [Point variants](../tracks/variation.md#varianttrack) | `VariantTrack`, `--variants` |
| Which two positions does a rearrangement join? | [Structural variants](../tracks/variation.md#structuraltrack) | `StructuralTrack`, `--structural` |
| How many copies are there, and has one allele been lost? | [Copy number](../tracks/variation.md#copynumbertrack) | `CopyNumberTrack`, `--copy-number` with `--ploidy` |
| Which sites tell closely related samples apart? | [Variable sites](../tracks/variation.md#snptrack) | `SnpTrack`, `--snps` |
| Which samples carry what, site by site? | [Genotype matrix](../tracks/variation.md#matrixtrack) | `MatrixTrack`, `--matrix` |
| Which samples lost or gained a stretch, window by window? | [Heatmap of samples](../tracks/variation.md#matrixtrack) | `MatrixTrack::windows`, `--heatmap` |
| Where does a scan cross its significance line? | [Association scan](../tracks/variation.md#manhattantrack) | `ManhattanTrack`, `--manhattan` |
| Which markers of a peak travel with its strongest? | [Scan coloured by linkage](../tracks/variation.md#manhattantrack) | `ManhattanTrack::linkage`, `--manhattan` with `--ld` |
| Which variants are inherited together, or which places are in contact? | [Pairs of positions](../tracks/variation.md#pairtrack) | `PairTrack`, `--pairs` |
| Which codons are under selection, and in which direction? | [Site-wise selection](../tracks/variation.md#selectiontrack) | `SelectionTrack`, `--selection` |

## Plots

<div class="k-plots" markdown>

-   [![Variant lollipops coloured by consequence, below a depth profile with a dropout, a reference band and the rpoB gene model](../assets/figures/example.svg){ width="900" height="305" loading="lazy" }](../tracks/variation.md#varianttrack)

    **[Point variants](../tracks/variation.md#varianttrack)**
    Point events as lollipops whose height is a value, or as plain ticks when there are too many for heads, with categories coloured in the order they first appear.

-   [![Five structural calls as arcs between their breakpoints, a deletion, a duplication, an inversion, an insertion and a translocation leaving the view, above a depth profile that drops under the deletion and steps up under the duplication](../assets/figures/example-structural.svg){ width="880" height="288" loading="lazy" }](../tracks/variation.md#structuraltrack)

    **[Structural variants](../tracks/variation.md#structuraltrack)**
    Each call as an arc between its two breakpoints, heavier with more supporting reads; a coverage track beneath shows whether the depth agrees.

-   [![Copy number gains and losses across a cohort along one arm of chromosome 8, and beneath them one tumour's segments on a ladder of whole copies, with lost heterozygosity marked along the foot](../assets/figures/example-copy-number.svg){ width="900" height="341" loading="lazy" }](../tracks/variation.md#copynumbertrack)

    **[Copy number](../tracks/variation.md#copynumbertrack)**
    Fitted segments drawn at their level on a ladder of whole copies, loss of heterozygosity in a lane of its own, and balanced wherever you say it is.

-   [![A phylogeny of twelve isolates beside lineage, resistance and year strips and a panel of thirty-four variable sites, each column labelled with its position and each row ending in its count of differences](../assets/figures/example-snps.svg){ width="900" height="388" loading="lazy" }](../tracks/variation.md#snptrack)

    **[Variable sites](../tracks/variation.md#snptrack)**
    Only the columns that vary, evenly spaced and each labelled with its position; a tree beside the rows lines a clade's shared changes up into a block.

-   [![An association scan over rpoB whose peak crosses the threshold line, the gene beneath it, and a genotype matrix showing which isolates carry the associated alleles](../assets/figures/example-association.svg){ width="940" height="348" loading="lazy" }](../tracks/variation.md#matrixtrack)

    **[Genotype matrix](../tracks/variation.md#matrixtrack)**
    One row per sample and one cell per site at its real coordinate, where a sample without the allele, a sample never typed and a stretch with no site all look different.

-   [![Presence and absence of twenty-six genes across nine Klebsiella isolates, the rows sorted by the phylogeny beside them so two accessory islands come out as solid blocks](../assets/figures/example-pangenome.svg){ width="940" height="264" loading="lazy" .k-wide }](../tracks/variation.md#matrixtrack)

    **[Presence and absence by descent](../tracks/variation.md#matrixtrack)**
    The same matrix sorted by a tree drawn beside it, which turns a speckle into blocks a clade carries.

-   [![An association scan and a depth profile across every contig of a draft assembly, with the contigs beneath as alternating named blocks](../assets/figures/example-genomewide.svg){ width="940" height="289" loading="lazy" }](../tracks/variation.md#manhattantrack)

    **[Association scan](../tracks/variation.md#manhattantrack)**
    One point per test, a threshold you set and the points above it ringed; laid over a `Genome`, the scan runs across a whole assembly.

-   [![A peak of an association scan, each marker coloured from grey to blue by its linkage with the strongest, which is a diamond labelled with its name, with a recombination rate with two hotspots laid behind it](../assets/start/locus.svg){ .k-light width="720" height="214" loading="lazy" }![A peak of an association scan coloured by linkage](../assets/start/locus-dark.svg){ .k-dark width="720" height="214" loading="lazy" }](../tracks/variation.md#manhattantrack)

    **[Scan coloured by linkage](../tracks/variation.md#manhattantrack)**
    A peak as LocusZoom draws one: every marker coloured by its r² with the lead, so a second signal beside the first stands out grey.

-   [![Forty samples in the order of a tree, each a row of cells along a chromosome, with one clade missing a stretch and three samples carrying another twice](../assets/start/heatmap.svg){ .k-light width="720" height="602" loading="lazy" }![A heatmap of forty samples in windows](../assets/start/heatmap-dark.svg){ .k-dark width="720" height="602" loading="lazy" }](../tracks/variation.md#matrixtrack)

    **[Heatmap of samples](../tracks/variation.md#matrixtrack)**
    A value per sample per window, each sample read against its own usual value, in the order of a tree.

-   [![A triangle under a gene in which each pair of variants is a cell coloured by its linkage, with three dark blocks of variants inherited together](../assets/start/pairs.svg){ .k-light width="720" height="386" loading="lazy" }![A triangle of linkage under a gene](../assets/start/pairs-dark.svg){ .k-dark width="720" height="386" loading="lazy" }](../tracks/variation.md#pairtrack)

    **[Pairs of positions](../tracks/variation.md#pairtrack)**
    Linkage or contacts as a triangle hung under the axis, and a few pairs far apart as arcs.

-   [![A molecular selection atlas: rate classes and recurrent changes on a rectangular tree, mean branch omega on a circular tree, and two site-wise scans over protein domains with evidence above signed omega effects](../assets/figures/example-selection-atlas.svg){ width="1508" height="1075" loading="lazy" }](../tracks/variation.md#selectiontrack)

    **[Site-wise selection](../tracks/variation.md#selectiontrack)**
    Evidence, as a p-value or a posterior, in one tier and the signed log2(ω) effect in another, so a significant purifying site still reads as purifying.

</div>

## Related

<div class="grid cards" markdown>

-   **[Reads and molecules](reads-molecules.md)**

    The reads and molecules behind a call.

-   **[Phylogeny and clades](phylogeny-clades.md)**

    Calls read on a tree, and spans painted onto the clades that carry them.

-   **[Annotation and coordinates](annotation-coordinates.md)**

    The genes and codons a call lands in.

-   **[Whole genomes and maps](whole-genomes-geography.md)**

    A whole assembly laid end to end under a genome-wide scan.

</div>
