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

<div class="plot-stats"><span><strong>3</strong> tree tracks</span><span><strong>4</strong> selection views</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../evolution-surveillance/">
    <span class="plot-card__media"><img src="../../assets/figures/example-evolutionary-surveillance.svg" alt="Tree geometries, ancestral reconstruction, molecular evolution, comparison, phylodynamics and lineage surveillance" loading="lazy"></span>
    <span class="plot-card__body"><small>Topology × genome × time</small><strong>Evolutionary surveillance atlas</strong><span>Move from tree geometry through ancestral evidence and molecular selection to phylodynamic estimates and observed lineage change.</span><b>Open catalogue <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/phylogenetics/#build-a-branch-to-codon-selection-atlas">
    <span class="plot-card__media"><img src="../../assets/figures/example-selection-atlas.svg" alt="A molecular-selection atlas with branch rate mixtures, recurrent-event links and site-wise frequentist and posterior scans" loading="lazy"></span>
    <span class="plot-card__body"><small>Tree × codon model</small><strong>Molecular selection atlas</strong><span>Keep weighted branch rate classes, recurrence, statistical evidence and signed site effects in separate visual channels.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/phylogenetics/#attach-data-graphics-to-nodes-and-clades">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-faces.svg" alt="Node bubbles, ancestral composition donuts, clade fields, an aligned MSA and aligned domain architectures" loading="lazy"></span>
    <span class="plot-card__body"><small>Tree × node and row data</small><strong>Phylogenetic data faces</strong><span>Attach abundance, composition and clade context to nodes, then align sequences or domains by descent.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#treetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-evidence.svg" alt="Rectangular, circular and unrooted trees with support, branch events and evolutionary scale bars" loading="lazy"></span>
    <span class="plot-card__body"><small>One annotated topology</small><strong>TreeTrack</strong><span>Root by node, outgroup or midpoint; then draw rectangular, radial or unrooted with support, events, scales and iTOL-style metadata.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
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
| ancestral probabilities, direct events or branch uncertainty | `TreeTrack` with `AncestralStateLayer`, `BranchEventLayer` or `BranchIntervalLayer` |
| abundance or composition on internal nodes | `TreeTrack` with `NodeGlyph` |
| aligned residues or domains ordered by descent | `MsaTrack::tree` or `DomainTrack::tree` |
| two trees over matching terminal names | `TanglegramTrack` |
| genomic intervals carried by named taxa | `CladeTrack` |
| terminal annotations matched to supplied coordinates | `PhyloMap` |
| inferred population trajectory through time | `PhylodynamicTrack` |
| observed lineage counts and denominators through time | `SurveillanceTrack` |

## Related guides

- [Annotated phylogenetics](../guide/phylogenetics.md) covers parsing, dates, radial geometry and topology operations.
- [Geographic genomics](../guide/maps.md) covers coordinate matching, projections and missing-location rules.
- [Variation and association](variation-association.md) covers variable sites and genotype patterns organised by a tree.
