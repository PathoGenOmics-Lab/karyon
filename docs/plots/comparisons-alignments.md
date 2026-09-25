# Comparisons and alignments

Plots that compare two or more sequences without losing orientation:
alignments column by column, domain architectures, dotplots, synteny ribbons
and gene neighbourhoods.
{ .k-lead }

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| Where do aligned sequences disagree, column by column? | [Multiple alignment](../tracks/comparison.md#msatrack) | `MsaTrack`, `--msa` |
| Which domains does each protein carry, and in what order? | [Domain architectures](../tracks/comparison.md#domaintrack) | `DomainTrack`, `--domains` |
| What shape is the rearrangement between two sequences? | [Dotplot](../tracks/comparison.md#dotplottrack) | `DotplotTrack`, `--dotplot` |
| Where did each block of one sequence go in the other, in a band that stacks with the rest? | [Synteny ribbons](../tracks/comparison.md#syntenytrack) | `SyntenyTrack`, `--synteny` |
| What does one genome's locus have that another's lacks? | [Gene neighbourhoods](../tracks/comparison.md#locustrack) | `LocusTrack`, `--loci` with `--links` |

The dotplot and the ribbons take the same alignment blocks, so the same PAF
file can be drawn either way. For two trees over the same taxa, use a
[tanglegram](phylogeny-clades.md).

## Plots

<div class="k-plots" markdown>

-   [![A conservation logo above an alignment of twelve isolates over 120 columns, each row a quiet bar with only the bases that disagree with the first isolate painted](../assets/figures/example-msa.svg){ width="940" height="326" loading="lazy" }](../tracks/comparison.md#msatrack)

    **[Multiple alignment](../tracks/comparison.md#msatrack)**
    Rows over alignment columns with only the disagreements painted, against the consensus or a row you name.

-   [![Twenty-two columns of a KatG protein alignment across five strains, every residue shown and coloured by physicochemical class, with one gap](../assets/figures/example-msa-protein.svg){ width="720" height="174" loading="lazy" .k-wide }](../tracks/comparison.md#msatrack)

    **[Protein alignment](../tracks/comparison.md#msatrack)**
    Every residue shown and coloured by one of six physicochemical classes.

-   [![Node bubbles and host bars on a rectangular tree, ancestral-host donuts on a radial tree, and an alignment and four protein architectures of sensor, repeat and kinase domains with rows sorted by the tree beside them](../assets/figures/example-phylo-faces.svg){ width="1386" height="660" loading="lazy" }](../tracks/comparison.md#domaintrack)

    **[Domain architectures](../tracks/comparison.md#domaintrack)**
    Domains, motifs or repeats as labelled intervals, one protein per row and one colour per label; a tree beside the rows makes a gain or a loss read by clade.

-   [![Two M. tuberculosis chromosomes compared as a dotplot, where an inversion runs as an anti-diagonal and a translocated block sits off the main diagonal, above the same alignment drawn as ribbons](../assets/figures/example-synteny.svg){ width="900" height="437" loading="lazy" }](../tracks/comparison.md#dotplottrack)

    **[Dotplot](../tracks/comparison.md#dotplottrack)**
    Alignment blocks as diagonals on two axes: an inversion is an anti-diagonal, and a translocation is a diagonal that has slid sideways.

-   [![An inversion close up: two sequence bars joined by ribbons that cross where the alignment reverses, above a feature marking the inverted segment](../assets/figures/example-synteny-inversion.svg){ width="760" height="251" loading="lazy" }](../tracks/comparison.md#syntenytrack)

    **[Synteny ribbons](../tracks/comparison.md#syntenytrack)**
    The same blocks as ribbons between two bars, compact enough to stack under coverage and genes; an inversion is a twist.

-   [![The ESX-1 locus in H37Rv, M. bovis and BCG, genes drawn as arrows and joined by identity ribbons, with the block deleted in BCG left without ribbons](../assets/figures/example-cluster.svg){ width="880" height="258" loading="lazy" }](../tracks/comparison.md#locustrack)

    **[Gene neighbourhoods](../tracks/comparison.md#locustrack)**
    One row per genome, genes as arrows and identity ribbons to the row below; genes with no match are outlined, because the absence is usually the finding.

</div>

## Related

<div class="grid cards" markdown>

-   **[Variation and association](variation-association.md)**

    Only the variable sites, when most of an alignment agrees.

-   **[Phylogeny and clades](phylogeny-clades.md)**

    Alignments and domains sorted by a tree, and two trees face to face.

-   **[Whole genomes and maps](whole-genomes-geography.md)**

    Context across a whole chromosome or assembly.

-   **[File formats](../guide/formats.md)**

    How PAF, aligned FASTA and InterProScan tables are read.

</div>
