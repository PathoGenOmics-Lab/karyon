---
title: Phylogeny and clade plots
description: Annotated trees, tanglegrams, clade-linked intervals and circular phylogeography.
---

<div class="plot-hero plot-hero--phylo" markdown>

<span class="plot-eyebrow">Plot catalogue · Phylogeny and clades</span>

# Put topology beside the evidence

These views preserve terminal order, branch length, support and typed metadata.
They render an existing analysis; they do not infer a tree, clock, ancestral
state, migration event or transmission chain.

<div class="plot-stats"><span><strong>3</strong> tracks</span><span><strong>1</strong> standalone drawing</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#treetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-annotations.svg" alt="Unrooted and circular trees with layered metadata rings" loading="lazy"></span>
    <span class="plot-card__body"><small>One annotated topology</small><strong>TreeTrack</strong><span>Rectangular, circular, fan or unrooted trees with branch metadata and iTOL-style strips, bars, binary marks and symbols.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" id="tree-to-tree-comparison" href="../../tracks/#tanglegramtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-tanglegram.svg" alt="Two trees joined through their matching terminal taxa" loading="lazy"></span>
    <span class="plot-card__body"><small>Two topologies</small><strong>TanglegramTrack</strong><span>Matching taxa joined between trees, with crossing minimisation that changes order but never topology.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#cladetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-clades.svg" alt="Genomic intervals painted across the clades that carry them" loading="lazy"></span>
    <span class="plot-card__body"><small>Tree × genomic interval</small><strong>CladeTrack</strong><span>Coordinate blocks whose height is a clade, with non-carrier rows cut out rather than hidden.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/maps/#put-a-phylogeny-around-the-map">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-map.svg" alt="A circular time tree and cladogram around a map" loading="lazy"></span>
    <span class="plot-card__body"><small>Tree × supplied geography</small><strong>PhyloMap</strong><span>A circular tree around a central projection with aggregated, individual or absent sample-to-place links.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose what topology is aligned to

| Analytical structure | Start with |
|:--|:--|
| one tree and terminal metadata | `TreeTrack` |
| two trees over matching terminal names | `TanglegramTrack` |
| genomic intervals carried by named taxa | `CladeTrack` |
| terminal annotations matched to supplied coordinates | `PhyloMap` |

## Related guides

- [Annotated phylogenetics](../guide/phylogenetics.md) covers parsing, dates, radial geometry and topology operations.
- [Geographic genomics](../guide/maps.md) covers coordinate matching, projections and missing-location rules.
- [Variation and association](variation-association.md) covers variable sites and genotype patterns organised by a tree.
