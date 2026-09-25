# Gallery

Find the plot for the question you are asking of your data: start from what the
data describes, or from the file you already have.
{ .k-lead }

Every figure here is drawn by karyon in your browser, in this page's light or
dark. On the pages they lead to, drag a figure along its genome, zoom it with
the buttons under it or a pinch, and hover a mark for what it holds.

## Browse by what you want to show

<div class="k-plots" markdown>

-   [![A read depth profile with a dropout, a reference band, the rpoB gene with its resistance-determining region, and variant lollipops, all on one axis](../assets/figures/example.svg){ width="900" height="304" loading="lazy" }](signal-sequence.md)

    **[Signal and sequence](signal-sequence.md)**
    A value or a symbol at every position: depth, windowed statistics, methylation by strand, the reference bases, logos and per-base attribution.

-   [![Three transcription units over the ESX-1 genes of M. tuberculosis, each a bent arrow at its start site with a hairpin or a bar at its end, above the gene models they carry](../assets/figures/example-transcripts.svg){ width="880" height="216" loading="lazy" .k-wide }](annotation-coordinates.md)

    **[Annotation and coordinates](annotation-coordinates.md)**
    Genes, transcription units and reading frames, with rulers in bases or codons and a key to the colours.

-   [![An association scan over rpoB whose peak crosses the threshold line, the gene beneath it, and a genotype matrix of the isolates at the associated sites](../assets/figures/example-association.svg){ width="940" height="347" loading="lazy" }](variation-association.md)

    **[Variation and association](variation-association.md)**
    From one call to a cohort: point and structural variants, copy number, variable sites, genotype matrices, association scans and site-wise selection.

-   [![Reads packed into rows under a depth profile and a candidate SNV, coloured by strand, with mismatches painted against the reference, a deletion, an insertion, faded low-quality reads and a count of the reads not shown](../assets/figures/example-pileup.svg){ width="920" height="472" loading="lazy" }](reads-molecules.md)

    **[Reads and molecules](reads-molecules.md)**
    The evidence one molecule at a time: pileups, split reads, per-read methylation, splice junctions and raw nanopore current.

-   [![Two M. tuberculosis chromosomes compared as a dotplot, where an inversion runs as an anti-diagonal and a translocated block sits off the main diagonal, above the same alignment drawn as ribbons](../assets/figures/example-synteny.svg){ width="900" height="436" loading="lazy" }](comparisons-alignments.md)

    **[Comparisons and alignments](comparisons-alignments.md)**
    Two or more sequences compared: alignments, domain architectures, dotplots, synteny ribbons and gene neighbourhoods.

-   [![Four radial views of one outbreak tree: a circular time tree with country and depth rings, a fan with a collapsed clade, time radiating inwards, and a circular cladogram](../assets/figures/example-phylo-layouts.svg){ width="1398" height="1226" loading="lazy" }](phylogeny-clades.md)

    **[Phylogeny and clades](phylogeny-clades.md)**
    Trees in every layout with their metadata, two trees face to face, spans painted onto clades, and a tree around a map.

-   [![An eight-panel synthetic atlas of tree geometries, ancestral states, selection on branches and at sites, a tanglegram, a population trajectory and lineage frequencies through time](../assets/figures/example-evolutionary-surveillance.svg){ width="1406" height="2047" loading="lazy" }](evolution-surveillance.md)

    **[Evolution and surveillance](evolution-surveillance.md)**
    Results fitted upstream: ancestral states, selection on branches and at sites, population trajectories and lineage counts through time.

-   [![The 4.41 Mb H37Rv chromosome as concentric rings of genes, GC content and GC skew, with resistance mutations marked and chords across the middle joining the ends of two rearrangements](../assets/figures/example-circular.svg){ width="688" height="688" loading="lazy" }](whole-genomes-geography.md)

    **[Whole genomes and maps](whole-genomes-geography.md)**
    Context beyond one locus: an ideogram, a whole assembly, a circular genome and a world map.

</div>

## Choose from the shape of your data

Find the file you have: the plot links to its reference and names the command
line flag that draws it, and the route opens its page in this gallery.

| Your data | Start with | Route |
|:--|:--|:--|
| Aligned reads: SAM text, as `samtools view` writes it | [Read pileup](../tracks/reads-molecules.md#pileuptrack), `--pileup` | [Reads and molecules](reads-molecules.md) |
| Reads aligned in pieces: SAM with an `SA` tag | [Split reads](../tracks/reads-molecules.md#splitreadtrack), `--split-reads` | [Reads and molecules](reads-molecules.md) |
| Splice junction counts: STAR's `SJ.out.tab` | [Splice junctions](../tracks/reads-molecules.md#junctiontrack), `--junctions` | [Reads and molecules](reads-molecules.md) |
| Methylation calls read by read: a Bismark methylation extractor file | [Methylation by molecule](../tracks/reads-molecules.md#bisulfitetrack), `--bisulfite` | [Reads and molecules](reads-molecules.md) |
| Raw nanopore current and a basecaller move table | [Nanopore signal](../tracks/reads-molecules.md#squiggletrack), Rust only | [Reads and molecules](reads-molecules.md) |
| A value per base: bedGraph, `samtools depth` or a column of values | [Coverage](../tracks/signal-sequence.md#coveragetrack), `--coverage` | [Signal and sequence](signal-sequence.md) |
| A statistic in windows: bedGraph | [Windowed statistic](../tracks/signal-sequence.md#windowtrack), `--windows` | [Signal and sequence](signal-sequence.md) |
| Modified bases per strand: bedMethyl from `modkit pileup` | [Methylation by strand](../tracks/signal-sequence.md#methylationtrack), `--methylation` | [Signal and sequence](signal-sequence.md) |
| Per-base model attribution: bedGraph, plus the reference as FASTA | [Per-base attribution](../tracks/signal-sequence.md#dynseqtrack), `--dynseq` with `--with-sequence` | [Signal and sequence](signal-sequence.md) |
| Aligned sequences of a motif: aligned FASTA | [Sequence logo](../tracks/signal-sequence.md#logotrack), `--logo` | [Signal and sequence](signal-sequence.md) |
| A reference sequence: FASTA | [Reference sequence](../tracks/signal-sequence.md#sequencetrack), `--sequence`, or [six reading frames](../tracks/annotation.md#orftrack), `--orfs` | [Signal and sequence](signal-sequence.md), [Annotation and coordinates](annotation-coordinates.md) |
| Genes and other intervals: BED or GFF3 | [Features](../tracks/annotation.md#featuretrack), `--features` | [Annotation and coordinates](annotation-coordinates.md) |
| Point calls: VCF | [Point variants](../tracks/variation.md#varianttrack), `--variants` | [Variation and association](variation-association.md) |
| Structural calls: VCF with symbolic alleles or `SVTYPE` | [Structural variants](../tracks/variation.md#structuraltrack), `--structural` | [Variation and association](variation-association.md) |
| Copy number segments: CNVkit `.cns`, ASCAT or `.seg` | [Copy number](../tracks/variation.md#copynumbertrack), `--copy-number` with `--ploidy` | [Variation and association](variation-association.md) |
| Closely related genomes: aligned FASTA | [Variable sites](../tracks/variation.md#snptrack), `--snps`, or a [multiple alignment](../tracks/comparison.md#msatrack), `--msa` | [Variation and association](variation-association.md), [Comparisons and alignments](comparisons-alignments.md) |
| A value per sample per site: a table | [Genotype matrix](../tracks/variation.md#matrixtrack), `--matrix` | [Variation and association](variation-association.md) |
| Association statistics: a table of position and value | [Association scan](../tracks/variation.md#manhattantrack), `--manhattan` | [Variation and association](variation-association.md) |
| Site-wise results from a codon model | [Site-wise selection](../tracks/variation.md#selectiontrack), Rust only | [Variation and association](variation-association.md), [Evolution and surveillance](evolution-surveillance.md) |
| A pairwise alignment: PAF from `minimap2` | [Dotplot](../tracks/comparison.md#dotplottrack), `--dotplot`, or [synteny ribbons](../tracks/comparison.md#syntenytrack), `--synteny` | [Comparisons and alignments](comparisons-alignments.md) |
| Protein domains: an InterProScan table | [Domain architectures](../tracks/comparison.md#domaintrack), `--domains` | [Comparisons and alignments](comparisons-alignments.md) |
| Genes from several genomes, and the hits between them | [Gene neighbourhoods](../tracks/comparison.md#locustrack), `--loci` with `--links` | [Comparisons and alignments](comparisons-alignments.md) |
| A phylogeny: Newick, with or without BEAST and NHX annotations | [Annotated tree](../tracks/phylogeny.md#treetrack), `--tree` | [Phylogeny and clades](phylogeny-clades.md) |
| Two trees over the same taxa | [Tanglegram](../tracks/phylogeny.md#tanglegramtrack), `--tanglegram` with `--against` | [Phylogeny and clades](phylogeny-clades.md) |
| Spans carried by named taxa: Gubbins GFF3, and the tree | [Clade blocks](../tracks/phylogeny.md#cladetrack), `--clades` with `--with-tree` | [Phylogeny and clades](phylogeny-clades.md) |
| Ancestral states, branch events or branch ω, as tree annotations | [Ancestral states and branch events](../tracks/phylogeny.md#treetrack), Rust only | [Evolution and surveillance](evolution-surveillance.md) |
| Estimates through time, such as effective population size or R | [Phylodynamic trajectory](../tracks/evolution-surveillance.md#phylodynamictrack), Rust only | [Evolution and surveillance](evolution-surveillance.md) |
| Lineage counts and totals through time | [Lineage surveillance](../tracks/evolution-surveillance.md#surveillancetrack), Rust only | [Evolution and surveillance](evolution-surveillance.md) |
| Cytogenetic bands: a UCSC cytoBand table | [Ideogram](../tracks/whole-genome.md#ideogramtrack), `--ideogram` | [Whole genomes and maps](whole-genomes-geography.md) |
| An assembly of many contigs or chromosomes | [Whole assembly](../tracks/whole-genome.md#genometrack), Rust only | [Whole genomes and maps](whole-genomes-geography.md) |
| A circular chromosome or plasmid | [Circular genome](../guide/plot.md), Rust only | [Whole genomes and maps](whole-genomes-geography.md) |
| A latitude and longitude per sample | [Map](../guide/maps.md), Rust only | [Whole genomes and maps](whole-genomes-geography.md) |

!!! tip "BAM, CRAM and BCF"
    Binary formats come in through a pipe: `samtools view`, `samtools depth`
    and `bcftools view` already write the text these readers take, and a track
    file can be `-` for standard input. See [file formats](../guide/formats.md).

!!! note "Tracks and drawings"
    Most plots here are tracks: bands that stack over one shared coordinate
    axis, so any of them can sit above or below any other. The circular
    genome, the map and the tree around a map are drawings instead (`Rings`,
    `Map` and `PhyloMap`), each with a coordinate system of its own, and a
    `Panels` sheet puts drawings and figures side by side.

## Every track in one place

<div class="grid cards" markdown>

-   **[Track catalogue](../tracks/index.md)**

    All 36 track types, family by family, with the options of each and the
    files it reads.

-   **[Recipes](../recipes.md)**

    Complete programs that stack several of these plots into one figure.

-   **[Command line](../guide/cli.md)**

    Every flag, and how the order of the flags becomes the order of the stack.

</div>

Every figure in this gallery is drawn by a program in the repository's
`examples/` directory, so the code behind each one is there to read and copy.
