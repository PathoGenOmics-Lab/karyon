# Reads and molecules

Plots that keep the evidence one molecule at a time, before reads are
collapsed into depth, a call or an average per site.
{ .k-lead }

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| Does a call hold up in the reads behind it? | [Read pileup](../tracks/reads-molecules.md#pileuptrack) | `PileupTrack`, `--pileup` |
| Did one molecule align to several places, and in what order? | [Split reads](../tracks/reads-molecules.md#splitreadtrack) | `SplitReadTrack`, `--split-reads` |
| Are the methylated sites on the same molecules, or scattered across them? | [Methylation by molecule](../tracks/reads-molecules.md#bisulfitetrack) | `BisulfiteTrack`, `--bisulfite` |
| How many reads crossed each intron? | [Splice junctions](../tracks/reads-molecules.md#junctiontrack) | `JunctionTrack`, `--junctions` |
| What did the current look like before basecalling? | [Nanopore signal](../tracks/reads-molecules.md#squiggletrack) | `SquiggleTrack`, Rust only |

## Plots

<div class="k-plots" markdown>

-   [![Reads packed into rows under a depth profile, a candidate SNV and the reference, coloured by strand, with mismatches painted against the reference, a deletion, an insertion, faded low-quality reads and a count of the reads not shown](../assets/figures/example-pileup.svg){ width="920" height="473" loading="lazy" }](../tracks/reads-molecules.md#pileuptrack)

    **[Read pileup](../tracks/reads-molecules.md#pileuptrack)**
    Reads placed by their real CIGAR and packed into rows, with mismatches found against the reference and the reads past the row limit counted on the band.

-   [![Eight molecules each aligned in three segments, both ends at a new insertion site and the middle on the reference IS6110 copy on the reverse strand, with the connectors drawn under the rows, above a depth profile that doubles over that copy](../assets/figures/example-split.svg){ width="880" height="330" loading="lazy" }](../tracks/reads-molecules.md#splitreadtrack)

    **[Split reads](../tracks/reads-molecules.md#splitreadtrack)**
    One row per molecule and one bar per alignment, joined in the order the molecule ran, with any step back down the reference drawn under the row.

-   [![Sixteen molecules across the H19/IGF2 imprinting control region, one row each, filled circles for methylated CpGs and open circles for unmethylated ones: some reads are methylated almost throughout and others almost nowhere](../assets/figures/example-bisulfite.svg){ width="880" height="303" loading="lazy" }](../tracks/reads-molecules.md#bisulfitetrack)

    **[Methylation by molecule](../tracks/reads-molecules.md#bisulfitetrack)**
    One row per read, a filled circle where a site is methylated, an open one where it is not and nothing where the read did not reach, so two alleles show as stripes.

-   [![Splice junction arcs labelled with read counts above an RNA depth profile, with one junction held back and counted, and per-base model attribution beneath](../assets/figures/example-regulation.svg){ width="811" height="788" loading="lazy" }](../tracks/reads-molecules.md#junctiontrack)

    **[Splice junctions](../tracks/reads-molecules.md#junctiontrack)**
    One arc per intron, weighted and labelled by the reads that crossed it; arcs sit in lanes, so their height means nothing.

-   [![Raw nanopore current for one read against sample number, stepping between levels, with the base the basecaller assigned to each stretch written above it](../assets/figures/example-squiggle.svg){ width="880" height="185" loading="lazy" .k-wide }](../tracks/reads-molecules.md#squiggletrack)

    **[Nanopore signal](../tracks/reads-molecules.md#squiggletrack)**
    Raw current against sample number, an envelope of the extremes when zoomed out and the samples themselves when zoomed in; a move table adds the called bases.

</div>

## Related

<div class="grid cards" markdown>

-   **[Signal and sequence](signal-sequence.md)**

    After the reads have become depth or a methylated fraction per site.

-   **[Variation and association](variation-association.md)**

    After the evidence has become a call.

-   **[File formats](../guide/formats.md)**

    SAM and its `SA` tag, STAR junction tables and Bismark extractor files,
    and how BAM comes in through `samtools`.

-   **[Recipes](../recipes.md)**

    Complete programs that stack several tracks into one figure.

</div>
