---
title: Variation and association plots
description: Point variants, structural events, variable sites, selection scans, genotype matrices and association statistics.
---

<div class="plot-hero plot-hero--variation" markdown>

<span class="plot-eyebrow">Plot catalogue · Variation and association</span>

# From one call to a cohort pattern

Choose by the unit that carries the result: one position, two breakpoints, a
variable alignment column, a sample-by-site cell or a statistic tested across
the genome.

<div class="plot-stats"><span><strong>7</strong> track types</span><span><strong>5</strong> data shapes</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#selectiontrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-selection-atlas.svg" alt="Site-wise molecular selection evidence and signed omega effects aligned to protein domains" loading="lazy" width="1510" height="1057"></span>
    <span class="plot-card__body"><small>Tested coding position</small><strong>SelectionTrack</strong><span>Separate p-value or posterior evidence from the direction and magnitude of the synonymous-to-nonsynonymous rate ratio.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#varianttrack">
    <span class="plot-card__media"><img src="../../assets/figures/example.svg" alt="Variant calls aligned to coverage and genes" loading="lazy" width="900" height="306"></span>
    <span class="plot-card__body"><small>Point event</small><strong>VariantTrack</strong><span>SNPs, indels or other named positions shown as lollipops or dense ticks with value and category.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#structuraltrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-structural.svg" alt="Structural variants drawn as arcs between breakpoints" loading="lazy" width="880" height="289"></span>
    <span class="plot-card__body"><small>Two breakpoints</small><strong>StructuralTrack</strong><span>Deletions, inversions, duplications or translocations as span-aware arcs weighted by support.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#copynumbertrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-copy-number.svg" alt="A cohort copy number landscape over one arm, and one tumour's segmentation under it" loading="lazy" width="900" height="342"></span>
    <span class="plot-card__body"><small>Segments a caller fitted</small><strong>CopyNumberTrack</strong><span>Total and minor allele copies on a ladder, with lost heterozygosity in a lane of its own.</span><b>Open reference <span aria-hidden="true">&rarr;</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#snptrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-snps.svg" alt="Variable alignment sites ordered by a phylogeny, with lineage, resistance and year beside them" loading="lazy" width="900" height="387"></span>
    <span class="plot-card__body"><small>Variable columns only</small><strong>SnpTrack</strong><span>Invariant columns removed so shared substitutions become visible blocks across ordered samples.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#matrixtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-association.svg" alt="A genotype matrix aligned below association statistics" loading="lazy" width="940" height="349"></span>
    <span class="plot-card__body"><small>Sample × site</small><strong>MatrixTrack</strong><span>Genotype or presence cells where zero, missing and every categorical state remain distinct.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#manhattantrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-genomewide.svg" alt="Association statistics across every contig of an assembly" loading="lazy" width="940" height="290"></span>
    <span class="plot-card__body"><small>Genome-wide statistic</small><strong>ManhattanTrack</strong><span>Association values, explicit threshold and highlighted hits on a shared genomic or assembly axis.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose the unit of evidence

| Evidence | Start with | Pair it with |
|:--|:--|:--|
| one named position | `VariantTrack` | `FeatureTrack` or `CodonTrack` |
| two linked breakpoints | `StructuralTrack` | `SplitReadTrack` for molecule-level support |
| differences in an alignment | `SnpTrack` | `TreeTrack` to order the rows |
| a state for each sample and site | `MatrixTrack` | `ManhattanTrack` above it |
| one statistic per tested position | `ManhattanTrack` | `GenomeTrack` for multi-contig context |
| molecular-selection evidence and rate effect | `SelectionTrack` | `FeatureTrack` or `CodonTrack` |

## Related routes

- [Reads and molecules](reads-molecules.md) when the evidence needs to stay visible below the call.
- [Phylogeny and clades](phylogeny-clades.md) when variants are interpreted on a tree.
- [Annotation and coordinates](annotation-coordinates.md) to name the genes, codons and positions involved.
