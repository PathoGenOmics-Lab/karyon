---
title: Track catalogue
description: Every track type karyon draws, by family, with the Rust call that builds it and the command line flag that reads it.
---

# Track catalogue

Every track type karyon draws, grouped by what it shows, with the Rust call that builds it and the command line flag that reads it. Each name links to its full entry: options, defaults and pitfalls.
{ .k-lead }

All 36 work the same way. A track says how tall it wants to be and draws inside the band the figure gives it, on the shared horizontal scale, so any of them stacks with any other. Each has an `add_` method on [`plot()`](../guide/plot.md), and each is also a type you can build yourself and hand over with `add_track`, or push onto a `Figure`, which is the way in for an alternative constructor or a track you want to query before it is drawn:

```rust
use karyon::{plot, read, CoverageTrack, Region};

let region = Region::parse("NC_000962.3:761,001-763,000")?;
let text = std::fs::read_to_string("depth.bedgraph")?;
let depth = CoverageTrack::from_spans(&region, read::signal::spans(&text, &region, None)?);

// genes: Vec<Feature>
plot("NC_000962.3:761,001-763,000")?
    .add_track(depth.label("depth"))
    .add_features(genes)
    .label("genes")
    .save("rpoB.svg")?;
```

On the [command line](../guide/cli.md), 28 of the 36 have a flag. Each flag starts a track, the options after it describe that track, and every track takes `--label`.

## Signal and sequence

[Signal and sequence tracks](signal-sequence.md): values along the sequence, and the sequence itself.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [CoverageTrack](signal-sequence.md#coveragetrack) | A value per base: depth, GC content, mappability | `.add_coverage(values)` | `--coverage` |
| [WindowTrack](signal-sequence.md#windowtrack) | A statistic in windows, either side of a baseline | `.add_windows(windows)` | `--windows` |
| [MethylationTrack](signal-sequence.md#methylationtrack) | Methylation per site, one lane per strand | `.add_methylation(sites)` | `--methylation` |
| [SequenceTrack](signal-sequence.md#sequencetrack) | The reference bases | `.add_sequence(seq)` | `--sequence` |
| [LogoTrack](signal-sequence.md#logotrack) | A sequence logo, scored seven ways | `.add_logo(columns)` | `--logo` |
| [DynseqTrack](signal-sequence.md#dynseqtrack) | Model attribution as bases sized by their score | `.add_dynseq(start, seq, scores)` | `--dynseq` with `--with-sequence` |

## Annotation

[Annotation tracks](annotation.md): what is annotated on the sequence.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [FeatureTrack](annotation.md#featuretrack) | Genes and other intervals, packed into rows | `.add_features(features)` | `--features` |
| [TranscriptionUnitTrack](annotation.md#transcriptionunittrack) | Transcripts from start site to terminator | `.add_transcription_units(units)` | none |
| [OrfTrack](annotation.md#orftrack) | Stops and open reading frames in six frames | `.add_orfs(seq)` | `--orfs` |

## Variation

[Variation tracks](variation.md): how samples differ from the reference and from each other.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [VariantTrack](variation.md#varianttrack) | Point calls as lollipops or ticks | `.add_variants(variants)` | `--variants` |
| [StructuralTrack](variation.md#structuraltrack) | Structural calls as arcs between breakpoints | `.add_structural(variants)` | `--structural` |
| [CopyNumberTrack](variation.md#copynumbertrack) | Segmented copy number and lost heterozygosity | `.add_copy_number(segments, ploidy)` | `--copy-number` with `--ploidy` |
| [SnpTrack](variation.md#snptrack) | The variable sites of an alignment, one row per sample | `.add_snps(names, sites)` | `--snps` |
| [MatrixTrack](variation.md#matrixtrack) | A value per sample per site | `.add_matrix(sites, rows)` | `--matrix` |
| [ManhattanTrack](variation.md#manhattantrack) | Association statistics against a threshold | `.add_manhattan(points)` | `--manhattan` |
| [SelectionTrack](variation.md#selectiontrack) | Selection evidence above, the ω effect below | `.add_selection(sites)` | none |

## Reads and molecules

[Read and molecule tracks](reads-molecules.md): the evidence behind a call.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [PileupTrack](reads-molecules.md#pileuptrack) | Aligned reads, with mismatches painted | `.add_pileup(reads)` | `--pileup` |
| [SplitReadTrack](reads-molecules.md#splitreadtrack) | Molecules that aligned in pieces | `.add_split_reads(reads)` | `--split-reads` |
| [BisulfiteTrack](reads-molecules.md#bisulfitetrack) | Methylation one molecule at a time | `.add_bisulfite(sites, molecules)` | `--bisulfite` |
| [JunctionTrack](reads-molecules.md#junctiontrack) | Splice junctions as arcs with read counts | `.add_junctions(junctions)` | `--junctions` |
| [SquiggleTrack](reads-molecules.md#squiggletrack) | Raw nanopore current for one read | `.add_squiggle(signal)` | none |

## Comparison

[Comparison tracks](comparison.md): one sequence against others.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [MsaTrack](comparison.md#msatrack) | A multiple alignment, differences painted | `.add_msa(sequences)` | `--msa` |
| [DomainTrack](comparison.md#domaintrack) | Domain architectures, one protein per row | `.add_domains(rows)` | `--domains` |
| [DotplotTrack](comparison.md#dotplottrack) | Two sequences on two axes | `.add_dotplot(blocks)` | `--dotplot` |
| [SyntenyTrack](comparison.md#syntenytrack) | Alignment ribbons between two bars | `.add_synteny(blocks)` | `--synteny` |
| [LocusTrack](comparison.md#locustrack) | Gene neighbourhoods joined by homology | `.add_loci(loci)` | `--loci` with `--links` |

## Phylogeny

[Phylogeny tracks](phylogeny.md): trees, and what is painted on them.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [TreeTrack](phylogeny.md#treetrack) | A phylogeny in three projections | `.add_tree(tree)` | `--tree` |
| [TanglegramTrack](phylogeny.md#tanglegramtrack) | Two trees face to face, shared tips joined | `.add_tanglegram(left, right)` | `--tanglegram` with `--against` |
| [CladeTrack](phylogeny.md#cladetrack) | Genomic spans painted onto the clades carrying them | `.add_clades(tree, blocks)` | `--clades` with `--with-tree` |

## Evolution and surveillance

[Evolution and surveillance tracks](evolution-surveillance.md): time on the shared axis.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [PhylodynamicTrack](evolution-surveillance.md#phylodynamictrack) | A trajectory through time, with its interval | `.add_phylodynamics(points)` | none |
| [SurveillanceTrack](evolution-surveillance.md#surveillancetrack) | Lineage counts or frequencies through time | `.add_surveillance(observations)` | none |

## Whole genome

[Whole genome tracks](whole-genome.md): where a region sits, and figures across an assembly.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [IdeogramTrack](whole-genome.md#ideogramtrack) | A whole chromosome, with the region marked | `.add_ideogram(length, bands)` | `--ideogram` |
| [GenomeTrack](whole-genome.md#genometrack) | The sequences of an assembly, end to end | `.add_genome(genome)` | none |

## Scales and keys

[Scale and key tracks](scales-keys.md): what a figure is read by.

| Track | Draws | Rust | Command line |
|:--|:--|:--|:--|
| [AxisTrack](scales-keys.md#axistrack) | The coordinate ruler | `.add_axis()` | `--axis` |
| [CodonTrack](scales-keys.md#codontrack) | A ruler in codons, with the residues | `.add_codons(start, end, strand)` | none |
| [LegendTrack](scales-keys.md#legendtrack) | A key to the colours, as a band | `.add_legend(legend)` | none |

!!! note "Coordinates"
    Positions in Rust are 0-based and half-open, the BED convention, so a GFF3 interval `759806..763325` is `Feature::new(759_805, 763_325)` and a VCF `POS` is `POS - 1`; the readers convert for you. The exceptions are what a reader sees: locus strings and ruler labels are 1-based and inclusive, as samtools and IGV write them. Not every axis is genomic either. An alignment counts columns, a variable-site panel counts sites, a squiggle counts samples, a domain track counts residues, the time tracks count time points, and a tree measures branch length across. See [Coordinates](../how-it-works/coordinates.md).

!!! note "Circular sequences"
    A plasmid, an organelle genome or a circular chromosome is drawn with `Rings`, which is not a track: it maps position to an angle rather than to a horizontal scale, puts annotation, signal, markers and a ruler on concentric rings, and draws chords across the middle to join the two ends of a rearrangement. A `Rings` plot and a `Figure` can share one `Panels` sheet. See [The Rust API](../guide/plot.md).

## Metadata columns

Seven tracks draw a row per named thing, and each can carry columns of metadata beside its rows: [MatrixTrack](variation.md#matrixtrack), [MsaTrack](comparison.md#msatrack), [SnpTrack](variation.md#snptrack), [CladeTrack](phylogeny.md#cladetrack), [DomainTrack](comparison.md#domaintrack) and [LocusTrack](comparison.md#locustrack) through `.traits(...)`, and [TreeTrack](phylogeny.md#treetrack), which draws the same columns from its own annotations through `trait_column`. The track answers which ones; the columns say what they were.

```rust
use karyon::read;
use karyon::track::traits::Traits;
use karyon::{plot, MatrixRow, MatrixTrack};

let sheet = read::sheet::sheet(
    "sample\tlineage\thost\tdepth\n\
     S1\tL4\thuman\t72.5\n\
     S2\tL2\tbovine\t61\n\
     S3\tL4\t\t48.2\n",
)?;
let traits = Traits::new(sheet.rows).spread(sheet.columns);

let rows = vec![
    MatrixRow::new("S1", vec![1.0, 0.0]),
    MatrixRow::new("S2", vec![0.0, 1.0]),
    MatrixRow::new("S3", vec![1.0, 1.0]), // no host: its cell is an empty outline
];
plot("chr1:1-1,000")?
    .add_track(MatrixTrack::new(vec![120, 340], rows).traits(traits))
    .save("matrix.svg")?;
```

- **The join is by name**, so the strips follow whatever order the rows are in, including the order a phylogeny put them in. A row the sheet says nothing about gets an empty outline, the one mark in a strip that cannot be mistaken for a level.
- **`Traits::spread` picks the mark.** A column whose every stated value is a number gets a ramp; anything else gets the categorical palette, and a column with more levels than the palette's six gets `TraitStyle::Symbol`, which carries the level in a shape as well as a hue. `TraitColumn::categorical`, `continuous`, `bar`, `binary` and `symbol` build a column by hand, for `Traits::column`.
- **Levels are numbered as they are first met**, never sorted, so a figure redrawn from the same file colours the same way, and one more sample does not repaint the others.
- **One vocabulary.** The strips beside a tree and beside a matrix are drawn by the same code, so the same lineage is the same colour in both, which is most of the reason to put them in one figure.
- **Nothing adds a key for you.** `Traits::legend(&theme)` builds one naming every level and both ends of every ramp; where it goes, and whether the figure needs it, is yours to decide.

From the command line this is `--traits FILE`, and `--columns A,B` to choose and order the columns, after `--matrix`, `--msa`, `--snps`, `--clades`, `--domains`, `--loci` or `--tree`. The sample sheet format is in [File formats](../guide/formats.md).

## From the command line

The Command line column above maps each flag to its track. A few flags need company:

- **A second file.** `--dynseq` needs `--with-sequence`, `--tanglegram` needs `--against`, `--clades` needs `--with-tree` and `--loci` needs `--links`, and each is refused without it. `--pileup` takes `--with-sequence` optionally, to find mismatches.
- **A choice inside the file.** `--methylation` takes `--modification`, `--bisulfite` takes `--context` and `--domains` takes `--analysis`, for a file that holds several datasets; the command refuses to pick one for you.
- **A number the file does not hold.** `--copy-number` needs `--ploidy`, since where balanced sits is not in the file.
- **Standard input.** Any track file may be `-`, for one track per command, which is how BAM, CRAM and BCF get in: `samtools` and `bcftools` already write the text these readers take.

Eight tracks are library only. `TranscriptionUnitTrack`, `SelectionTrack`, `PhylodynamicTrack`, `SurveillanceTrack`, `CodonTrack` and `GenomeTrack` would need a table with no single standard behind it; `SquiggleTrack` reads raw current, which comes in binary formats; and `LegendTrack` is built from what the other tracks drew rather than from a file. The whole grammar is in [Command line](../guide/cli.md).

## Where next

<div class="grid cards" markdown>

-   **[The Rust API](../guide/plot.md)**

    The `plot()` builder and `Figure`, which every track here goes through.

-   **[Command line](../guide/cli.md)**

    The grammar the flags above follow, and every option they take.

-   **[File formats](../guide/formats.md)**

    What each reader accepts, column by column.

-   **[Writing a track](../how-it-works/extending.md)**

    What a new track type has to do, and why every track lives on the shared axis.

</div>

<script>
  // Entries used to live on this page; send an old deep link on to the
  // family page that holds its entry now, keeping the anchor.
  (function () {
    var family = {
      coveragetrack: "signal-sequence",
      windowtrack: "signal-sequence",
      methylationtrack: "signal-sequence",
      sequencetrack: "signal-sequence",
      logotrack: "signal-sequence",
      dynseqtrack: "signal-sequence",
      featuretrack: "annotation",
      transcriptionunittrack: "annotation",
      orftrack: "annotation",
      varianttrack: "variation",
      structuraltrack: "variation",
      copynumbertrack: "variation",
      snptrack: "variation",
      matrixtrack: "variation",
      manhattantrack: "variation",
      selectiontrack: "variation",
      pileuptrack: "reads-molecules",
      splitreadtrack: "reads-molecules",
      bisulfitetrack: "reads-molecules",
      junctiontrack: "reads-molecules",
      squiggletrack: "reads-molecules",
      msatrack: "comparison",
      domaintrack: "comparison",
      dotplottrack: "comparison",
      syntenytrack: "comparison",
      locustrack: "comparison",
      treetrack: "phylogeny",
      tanglegramtrack: "phylogeny",
      cladetrack: "phylogeny",
      phylodynamictrack: "evolution-surveillance",
      surveillancetrack: "evolution-surveillance",
      ideogramtrack: "whole-genome",
      genometrack: "whole-genome",
      axistrack: "scales-keys",
      codontrack: "scales-keys",
      legendtrack: "scales-keys"
    };
    var anchor = window.location.hash.slice(1);
    if (Object.prototype.hasOwnProperty.call(family, anchor)) {
      window.location.replace(family[anchor] + "/#" + anchor);
    }
  })();
</script>
