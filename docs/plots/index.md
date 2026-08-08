---
title: Plot catalogue
description: Choose a Karyon plot by biological question, data shape or scale.
---

<div class="plot-hero plot-hero--catalog" markdown>

<span class="plot-eyebrow">Visual catalogue</span>

# Find the plot that matches the question

Start from the biological structure in the data, not from a Rust type name.
Each route below leads to a small family of plots, a visual example and the
exact reference entry for every component.

<div class="plot-stats" aria-label="Catalogue summary">
  <span><strong>30</strong> genomic tracks</span>
  <span><strong>3</strong> standalone drawings</span>
  <span><strong>7</strong> biological routes</span>
</div>

</div>

## Browse by biological task

<div class="plot-category-grid">
  <a class="plot-category-card plot-category-card--signal" href="signal-sequence/">
    <img src="../assets/figures/example.svg" alt="Coverage, sequence and annotation over one locus" loading="lazy">
    <span class="plot-category-card__body">
      <small>5 track types</small>
      <strong>Signal and sequence</strong>
      <span>Depth, signed windows, bases, methylation and motifs.</span>
      <b>Explore signal plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--annotation" href="annotation-coordinates/">
    <img src="../assets/figures/example-transcripts.svg" alt="Transcription units and genes over genomic coordinates" loading="lazy">
    <span class="plot-category-card__body">
      <small>6 track types</small>
      <strong>Annotation and coordinates</strong>
      <span>Genes, transcripts, reading frames, rulers and keys.</span>
      <b>Explore annotation plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--variation" href="variation-association/">
    <img src="../assets/figures/example-association.svg" alt="Association statistics aligned to genes and a genotype matrix" loading="lazy">
    <span class="plot-category-card__body">
      <small>5 track types</small>
      <strong>Variation and association</strong>
      <span>Variants, structural events, variable sites and genotypes.</span>
      <b>Explore variation plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--reads" href="reads-molecules/">
    <img src="../assets/figures/example-pileup.svg" alt="Read pileup with strand, mismatches and gaps" loading="lazy">
    <span class="plot-category-card__body">
      <small>4 track types</small>
      <strong>Reads and molecules</strong>
      <span>Pileups, split alignments, methylation patterns and raw signal.</span>
      <b>Explore molecule plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--comparison" href="comparisons-alignments/">
    <img src="../assets/figures/example-synteny.svg" alt="A dotplot and synteny ribbons comparing two sequences" loading="lazy">
    <span class="plot-category-card__body">
      <small>4 track types</small>
      <strong>Comparisons and alignments</strong>
      <span>MSAs, dotplots, synteny and homologous loci.</span>
      <b>Explore comparison plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--phylo" href="phylogeny-clades/">
    <img src="../assets/figures/example-phylo-map.svg" alt="Circular phylogenies around an orthographic map" loading="lazy">
    <span class="plot-category-card__body">
      <small>3 tracks + 1 drawing</small>
      <strong>Phylogeny and clades</strong>
      <span>Trees, tanglegrams, clade intervals and phylogeography.</span>
      <b>Explore phylogeny plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
  <a class="plot-category-card plot-category-card--world" href="whole-genomes-geography/">
    <img src="../assets/figures/example-maps.svg" alt="World maps under three geographic projections" loading="lazy">
    <span class="plot-category-card__body">
      <small>2 tracks + 2 drawings</small>
      <strong>Whole genomes and geography</strong>
      <span>Assemblies, ideograms, circular genomes and world maps.</span>
      <b>Explore context plots <span aria-hidden="true">→</span></b>
    </span>
  </a>
</div>

## Choose from the shape of the data

<div class="plot-route-grid">
  <a href="reads-molecules/"><strong>I have SAM alignments or molecules</strong><span>Start with pileups, split reads or per-molecule calls.</span></a>
  <a href="variation-association/"><strong>I have VCF calls, sites or genotypes</strong><span>Start with variants, variable sites, matrices or association statistics.</span></a>
  <a href="phylogeny-clades/"><strong>I have Newick, BEAST, NHX or Nexus</strong><span>Start with a tree, aligned clade blocks or circular phylogeography.</span></a>
  <a href="comparisons-alignments/"><strong>I have two or more sequences</strong><span>Start with an MSA, dotplot, synteny view or homologous locus.</span></a>
  <a href="signal-sequence/"><strong>I have a value for every base or window</strong><span>Start with coverage, a signed window track or per-site methylation.</span></a>
  <a href="whole-genomes-geography/"><strong>I need global context</strong><span>Start with an assembly track, ideogram, circular genome or map.</span></a>
</div>

## Tracks, drawings and sheets

<div class="plot-concept-strip">
  <div><span class="plot-chip">Track</span><strong>Shares a figure band</strong><p>Usually maps genomic position through the same horizontal scale as its neighbours.</p></div>
  <div><span class="plot-chip">Drawing</span><strong>Own coordinate system</strong><p>`Rings`, `Map` and `PhyloMap` render a complete circular or geographic document.</p></div>
  <div><span class="plot-chip">Panels</span><strong>Composes finished views</strong><p>Places figures and drawings on one labelled, aligned manuscript sheet.</p></div>
</div>

Need the exhaustive API behaviour rather than a visual route? Open the
[alphabetical track reference](../tracks.md), or jump to the [worked
recipes](../recipes.md) when the desired output is already clear.
