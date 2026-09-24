# Signal and sequence

Plots for a value or a symbol at every position: read depth, a statistic in
windows, methylation by strand, the reference bases, sequence logos and
per-base model attribution.
{ .k-lead }

## How to choose

Decide first what zero means for your values: the floor they rise from, as for
depth; a line they fall either side of, as for a statistic in windows; or
nothing at all, as for bases and letters.

| Your question | Plot | Build it with |
|:--|:--|:--|
| How much is there at each base, counting up from zero? | [Coverage](../tracks/signal-sequence.md#coveragetrack) | `CoverageTrack`, `--coverage` |
| Which side of its baseline does each window fall on? | [Windowed statistic](../tracks/signal-sequence.md#windowtrack) | `WindowTrack`, `--windows` |
| How methylated is each site, strand by strand? | [Methylation by strand](../tracks/signal-sequence.md#methylationtrack) | `MethylationTrack`, `--methylation` |
| What are the reference bases here? | [Reference sequence](../tracks/signal-sequence.md#sequencetrack) | `SequenceTrack`, `--sequence` |
| What does a motif look like across aligned sequences? | [Sequence logo](../tracks/signal-sequence.md#logotrack) | `LogoTrack`, `--logo` |
| Which bases did a model rely on, and in which direction? | [Per-base attribution](../tracks/signal-sequence.md#dynseqtrack) | `DynseqTrack`, `--dynseq` with `--with-sequence` |

## Plots

<div class="k-plots" markdown>

-   [![A read depth profile over the rpoB locus with a dropout, above a reference band, the rpoB gene model and variant calls](../assets/figures/example.svg){ width="900" height="304" loading="lazy" }](../tracks/signal-sequence.md#coveragetrack)

    **[Coverage](../tracks/signal-sequence.md#coveragetrack)**
    One value per base drawn up from a floor of zero, each pixel column showing its maximum, or its minimum when you are hunting dropouts.

-   [![pN/pS on a log2 scale and GC skew in windows along forty kilobases, each drawn either side of its own baseline, with the two sides of each line in different colours](../assets/figures/example-selection.svg){ width="880" height="232" loading="lazy" .k-wide }](../tracks/signal-sequence.md#windowtrack)

    **[Windowed statistic](../tracks/signal-sequence.md#windowtrack)**
    A statistic per window drawn either side of the line it is read against, such as pN/pS, GC skew or Tajima's D.

-   [![Dam methylation at GATC sites around the E. coli origin of replication, forward-strand calls above the line and reverse-strand calls below, with the reverse strand close to zero inside oriC](../assets/figures/example-methylation.svg){ width="880" height="196" loading="lazy" .k-wide }](../tracks/signal-sequence.md#methylationtrack)

    **[Methylation by strand](../tracks/signal-sequence.md#methylationtrack)**
    The methylated fraction at each site with each strand in its own lane, calls from too few reads dropped and counted, and the rest faded by depth.

-   [![Sixty bases of the rpoB locus at base resolution: a depth profile, the reference drawn as coloured letters, and variant calls](../assets/figures/example-zoom.svg){ width="900" height="221" loading="lazy" .k-wide }](../tracks/signal-sequence.md#sequencetrack)

    **[Reference sequence](../tracks/signal-sequence.md#sequencetrack)**
    The bases as letters when there is room, coloured blocks when there is not, and a prompt to zoom in once they are too thin to draw.

-   [![One eight-column DNA motif drawn three ways: as probabilities, as information in bits, and as enrichment above the line with depletion below it](../assets/figures/example-logo.svg){ width="900" height="378" loading="lazy" }](../tracks/signal-sequence.md#logotrack)

    **[Sequence logo](../tracks/signal-sequence.md#logotrack)**
    A motif from aligned sequences or a weight matrix, scored as probability, information or enrichment against a background, with depletion hanging below the line.

-   [![Four columns scored five ways against a background, as log odds, KL divergence, difference, ratio and odds ratio, each putting the emphasis on a different column](../assets/figures/example-logo-scores.svg){ width="760" height="542" loading="lazy" }](../tracks/signal-sequence.md#logotrack)

    **[Logo scores compared](../tracks/signal-sequence.md#logotrack)**
    Seven scores to choose from, and the choice decides which column reads loudest, not only how the logo looks.

-   [![The same motif proportions from five, fifty and five hundred sequences, each drawn raw and shrunk: the raw logos look alike, while the shrunk ones grow with the number of sequences](../assets/figures/example-logo-stability.svg){ width="700" height="484" loading="lazy" }](../tracks/signal-sequence.md#logotrack)

    **[Stabilised logo](../tracks/signal-sequence.md#logotrack)**
    Each column shrunk towards the background, more when it rests on fewer sequences, so five sequences do not look as certain as five hundred.

-   [![A sequence logo whose symbols are three-letter amino acid codes](../assets/figures/example-logo-protein.svg){ width="640" height="204" loading="lazy" }](../tracks/signal-sequence.md#logotrack)

    **[Logo of any alphabet](../tracks/signal-sequence.md#logotrack)**
    Symbols are strings, so three-letter residue codes and k-mers stack the way bases do.

-   [![Junction arcs over an RNA depth profile with per-base model attribution beneath, and a close-up of a promoter motif where each base is a letter sized by its score, above or below the line](../assets/figures/example-regulation.svg){ width="810" height="786" loading="lazy" }](../tracks/signal-sequence.md#dynseqtrack)

    **[Per-base attribution](../tracks/signal-sequence.md#dynseqtrack)**
    One signed score per base drawn as the base itself: tall where the model relied on it, hanging below the line where it pulled the prediction down.

</div>

## Related

<div class="grid cards" markdown>

-   **[Reads and molecules](reads-molecules.md)**

    Methylation one read at a time, raw nanopore current, and the reads behind
    a depth profile.

-   **[Variation and association](variation-association.md)**

    When a position carries a call rather than a value.

-   **[Annotation and coordinates](annotation-coordinates.md)**

    Genes, reading frames and rulers to put under a signal.

-   **[File formats](../guide/formats.md)**

    What bedGraph, `samtools depth` and bedMethyl are read as, and which
    convention each one counts positions in.

</div>
