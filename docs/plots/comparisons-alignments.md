---
title: Comparison and alignment plots
description: Multiple alignments, pairwise dotplots, synteny ribbons and homologous loci.
---

<div class="plot-hero plot-hero--comparison" markdown>

<span class="plot-eyebrow">Plot catalogue · Comparisons and alignments</span>

# Compare sequences without losing orientation

Choose by whether the comparison is column-wise, pairwise, block-wise or
anchored on gene neighbourhood. Reversals remain reversals rather than being
normalised into forward-looking matches.

<div class="plot-stats"><span><strong>5</strong> track types</span><span><strong>2+</strong> sequences</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#msatrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-msa.svg" alt="Multiple nucleotide sequence alignment" loading="lazy"></span>
    <span class="plot-card__body"><small>Column-wise comparison</small><strong>MsaTrack</strong><span>Aligned nucleotides or residues coloured against a reference, consensus or residue class.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#domaintrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-faces.svg" alt="Protein domain architectures aligned to the leaves of a phylogeny" loading="lazy"></span>
    <span class="plot-card__body"><small>Interval architecture</small><strong>DomainTrack</strong><span>Domains, motifs, exons or repeats as labelled intervals, optionally sorted and justified by an adjacent tree.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#dotplottrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-synteny.svg" alt="Pairwise sequence dotplot with forward and reverse blocks" loading="lazy"></span>
    <span class="plot-card__body"><small>Two coordinate axes</small><strong>DotplotTrack</strong><span>Alignment blocks as diagonals or anti-diagonals, preserving exact position on both sequences.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#syntenytrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-synteny-inversion.svg" alt="Synteny ribbons crossing at an inversion" loading="lazy"></span>
    <span class="plot-card__body"><small>Block correspondence</small><strong>SyntenyTrack</strong><span>Two sequence bars joined by ribbons that cross where order or orientation changes.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#locustrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-pangenome.svg" alt="One homologous locus compared across several genomes" loading="lazy"></span>
    <span class="plot-card__body"><small>Gene neighbourhood</small><strong>LocusTrack</strong><span>Homologous loci as gene arrows and identity ribbons, with absent homologues kept explicit.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose the comparison unit

| The shared unit is… | Start with |
|:--|:--|
| an aligned column across many sequences | `MsaTrack` |
| domains or motifs along many sequences | `DomainTrack` |
| a coordinate pair between two sequences | `DotplotTrack` |
| an aligned block whose order matters | `SyntenyTrack` |
| a gene and its neighbourhood across genomes | `LocusTrack` |

For two trees over the same taxa, use the [TanglegramTrack route](phylogeny-clades.md#tree-to-tree-comparison).

## Related routes

- [Variation and association](variation-association.md) to keep only variable alignment columns.
- [Phylogeny and clades](phylogeny-clades.md) when evolutionary relationships organise the comparison.
- [Whole genomes and geography](whole-genomes-geography.md) for assembly-wide context.
