---
title: Evolution and surveillance plots
description: Tree geometry, ancestral reconstruction, molecular selection, phylodynamics and lineage surveillance.
---

<div class="plot-hero plot-hero--variation" markdown>

<span class="plot-eyebrow">Plot catalogue · Evolution and surveillance</span>

# Keep estimates, uncertainty and observations distinct

These views join topology, molecular evolution and time without flattening
their different evidence into one colour. Karyon renders supplied results: it
does not infer ancestral states, fit codon or coalescent models, smooth lineage
frequencies or assert transmission.

<div class="plot-stats"><span><strong>6</strong> tree geometries</span><span><strong>2</strong> temporal tracks</span><span><strong>1</strong> genomic selection track</span></div>

</div>

![Eight synthetic views covering orthogonal, diagonal, curved, circular and unrooted trees, a tanglegram, genomic selection layers, a phylodynamic trajectory and lineage surveillance](../assets/figures/example-evolutionary-surveillance.svg){ width="1410" height="2057" loading="lazy" }

## Choose the biological question

<div class="plot-card-grid">
  <a class="plot-card" href="../../guide/phylogenetics/#choose-a-tree-geometry-for-the-reading-task">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-layouts.svg" alt="Four views of one outbreak tree: an outward circular time tree with two trait rings, a partial fan with a collapsed clade, an inward time tree and a circular cladogram" loading="lazy" width="1402" height="1232"></span>
    <span class="plot-card__body"><small>Topology × geometry</small><strong>Tree shape</strong><span>Choose row alignment, compact radial context or topology-balanced unrooted space without changing the owned tree.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/phylogenetics/#render-ancestral-states-events-and-branch-uncertainty">
    <span class="plot-card__media"><img src="../../assets/figures/example-phylo-evidence.svg" alt="Rectangular, circular and unrooted phylograms carrying support, branch events and evolutionary distance scales" loading="lazy" width="1742" height="632"></span>
    <span class="plot-card__body"><small>Reconstruction × uncertainty</small><strong>Ancestral evidence</strong><span>Place posterior state composition on internal nodes and direct mutation or interval evidence on the branch that owns it.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/phylogenetics/#build-a-branch-to-codon-selection-atlas">
    <span class="plot-card__media"><img src="../../assets/figures/example-selection-atlas.svg" alt="Branch rate mixtures, recurrence links and genomic site-wise selection evidence" loading="lazy" width="1510" height="1057"></span>
    <span class="plot-card__body"><small>Branch × codon × genome</small><strong>Molecular evolution</strong><span>Separate fitted rate classes, effect direction, statistical evidence, recurrence and observed genomic variants.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#phylodynamictrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-evolutionary-surveillance.svg" alt="An eight-panel evolutionary atlas whose sixth panel is a log-scale effective population trajectory with its uncertainty ribbon" loading="lazy" width="1410" height="2057"></span>
    <span class="plot-card__body"><small>Time × inferred population process</small><strong>PhylodynamicTrack</strong><span>Draw skyline, effective population, reproductive-number or growth estimates with their interval still visible.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#surveillancetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-evolutionary-surveillance.svg" alt="An eight-panel evolutionary atlas whose last panel is stacked lineage frequencies with sampling-aware alert markers" loading="lazy" width="1410" height="2057"></span>
    <span class="plot-card__body"><small>Time × observed lineage counts</small><strong>SurveillanceTrack</strong><span>Compare lineage composition or trajectories while keeping counts, denominators and alert reasons inspectable.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#tanglegramtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-tanglegram.svg" alt="Core and accessory trees joined through matching taxa" loading="lazy" width="760" height="236"></span>
    <span class="plot-card__body"><small>Topology × topology</small><strong>Evolutionary comparison</strong><span>Expose gene-tree, species-tree, core or accessory disagreement as traceable crossings.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Data contract

| Supplied result | Karyon encoding | Deliberate boundary |
|:--|:--|:--|
| ancestral state probabilities | internal-node donuts and confident transition cues | no ancestral reconstruction |
| direct branch events | ordered branch symbols and optional recurrence links | no event inference or convergence claim |
| point estimate with bounds | branch whisker or temporal uncertainty ribbon | no interval estimation |
| branch or site ω results | neutral-centred colour, rate mixtures and separate evidence marks | no codon-model fitting |
| effective size, R or growth through time | linear or log trajectory | no coalescent or clock fitting |
| lineage counts and denominators | stacked composition or comparable lines | no smoothing or extrapolation |

Generate the complete synthetic sheet with:

```bash
cargo run --example evolutionary_surveillance -- assets
```

The visual contracts follow the structure of results produced by established
upstream workflows: [HyPhy methods](https://hyphy.org/methods/) for molecular
selection, [Augur ancestral](https://docs.nextstrain.org/projects/augur/en/latest/usage/cli/ancestral.html)
and [traits](https://docs.nextstrain.org/projects/augur/en/latest/usage/cli/traits.html)
for node and branch reconstruction, [Augur frequencies](https://docs.nextstrain.org/projects/augur/en/latest/usage/cli/frequencies.html)
for temporal observations, and [TreeTime coalescent
models](https://treetime.readthedocs.io/en/latest/merger_models.html) for
phylodynamic estimates. Karyon accepts typed values derived from any upstream
tool; it does not require or invoke these programs.

Continue with the [annotated phylogenetics guide](../guide/phylogenetics.md)
for builders and failure semantics, or open the [track API
reference](../tracks.md) for every option.
