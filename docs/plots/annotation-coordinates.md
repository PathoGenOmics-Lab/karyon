# Annotation and coordinates

Plots that name what lies along a sequence and show how to read positions on
it: gene models, transcription units, reading frames, rulers in bases or
codons, and a key to the colours.
{ .k-lead }

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| Where are the genes, exons, repeats or primers? | [Features](../tracks/annotation.md#featuretrack) | `FeatureTrack`, `--features` |
| Which genes share one transcript, and where does it start and stop? | [Transcription units](../tracks/annotation.md#transcriptionunittrack) | `TranscriptionUnitTrack`, Rust only |
| Could an unannotated stretch be coding, and on which strand? | [Six reading frames](../tracks/annotation.md#orftrack) | `OrfTrack`, `--orfs` |
| Which position is this, as a genome browser would number it? | [Coordinate ruler](../tracks/scales-keys.md#axistrack) | `AxisTrack`, added at the bottom by `plot()` and by the command line; `--axis` places it |
| Which residue is this, as in S450L or V600E? | [Codon ruler](../tracks/scales-keys.md#codontrack) | `CodonTrack`, Rust only |
| What do the colours and shapes stand for? | [Legend](../tracks/scales-keys.md#legendtrack) | `LegendTrack`, Rust only |

## Plots

<div class="k-plots" markdown>

-   [![A depth profile, a reference band, the rpoB gene with its resistance-determining region as a second, stranded feature on a row of its own, and variant calls](../assets/figures/example.svg){ width="900" height="304" loading="lazy" }](../tracks/annotation.md#featuretrack)

    **[Features](../tracks/annotation.md#featuretrack)**
    Intervals from BED or GFF3 packed into as few rows as the zoom allows, with strand shown by an arrowhead and by the strand colours the other tracks use.

-   [![Three transcription units over the ESX-1 genes of M. tuberculosis, each a bent arrow at its start site with a hairpin or a bar at its end, one of them leaderless, above the gene models](../assets/figures/example-transcripts.svg){ width="880" height="216" loading="lazy" .k-wide }](../tracks/annotation.md#transcriptionunittrack)

    **[Transcription units](../tracks/annotation.md#transcriptionunittrack)**
    One RNA from its start site through its leader to its terminator, so the genes in a feature track beneath it read as transcribed together.

-   [![Six reading frames across three and a half kilobases, stop codons as ticks and open stretches as bars, forward frames above the line and reverse frames below](../assets/figures/example-frames.svg){ width="880" height="170" loading="lazy" .k-wide }](../tracks/annotation.md#orftrack)

    **[Six reading frames](../tracks/annotation.md#orftrack)**
    Stops and open stretches in all six frames, worked out from the sequence alone; requiring a start codon is a separate switch, off by default.

-   [![Two figures over the rpoB locus, each ending in a coordinate ruler labelled in kilobases: depth, annotation and allele frequencies above, a centred strand-composition score below](../assets/figures/example-visual-system.svg){ width="832" height="561" loading="lazy" }](../tracks/scales-keys.md#axistrack)

    **[Coordinate ruler](../tracks/scales-keys.md#axistrack)**
    Round, 1-based tick labels in one unit for the whole ruler: the numbers you would type into a genome browser.

-   [![rpoB codons 439 to 465 numbered and translated, with variant lollipops over codons 445 and 450 and a base ruler underneath](../assets/figures/example-codons.svg){ width="880" height="169" loading="lazy" .k-wide }](../tracks/scales-keys.md#codontrack)

    **[Codon ruler](../tracks/scales-keys.md#codontrack)**
    Codons numbered and translated along a coding sequence, counting from the far end on the reverse strand, so a change can be named by its residue.

-   [![One gene cluster in three genomes joined by identity ribbons, with a legend band beneath it: a filled swatch for the deleted block, an identity ramp and an outline for genes with no match](../assets/figures/example-cluster.svg){ width="880" height="257" loading="lazy" }](../tracks/scales-keys.md#legendtrack)

    **[Legend](../tracks/scales-keys.md#legendtrack)**
    A key as a band of its own, placed where you push it and wrapped onto another row rather than dropping an entry.

</div>

## Related

<div class="grid cards" markdown>

-   **[Signal and sequence](signal-sequence.md)**

    The measurements an annotation sits beside.

-   **[Variation and association](variation-association.md)**

    Calls that land inside genes and codons.

-   **[Whole genomes and maps](whole-genomes-geography.md)**

    Where the locus sits in its chromosome or its assembly.

-   **[Coordinates](../how-it-works/coordinates.md)**

    Why positions are 0-based inside karyon and 1-based wherever a reader
    sees them.

</div>
