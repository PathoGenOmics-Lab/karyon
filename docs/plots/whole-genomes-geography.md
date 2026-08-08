---
title: Whole-genome and geographic plots
description: Ideograms, concatenated assemblies, circular chromosomes and geographic maps.
---

<div class="plot-hero plot-hero--world" markdown>

<span class="plot-eyebrow">Plot catalogue · Whole genomes and geography</span>

# Preserve context beyond one locus

These views answer where a window sits in a chromosome, where a signal sits in
an assembly, how features wrap around a circular replicon, or where supplied
samples sit on Earth.

<div class="plot-stats"><span><strong>2</strong> tracks</span><span><strong>2</strong> standalone drawings</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#ideogramtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-ideogram.svg" alt="Chromosome ideogram with the visible locus marked" loading="lazy"></span>
    <span class="plot-card__body"><small>Window in chromosome</small><strong>IdeogramTrack</strong><span>A banded or bare whole chromosome with an explicit marker for the region shown below.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#genometrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-genomewide.svg" alt="Signals and association statistics across a multi-contig assembly" loading="lazy"></span>
    <span class="plot-card__body"><small>Multi-sequence assembly</small><strong>GenomeTrack</strong><span>Contigs placed end to end as named blocks so every neighbouring track can use one global scale.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/figure/#rings-a-circular-sequence">
    <span class="plot-card__media"><img src="../../assets/figures/example-circular.svg" alt="Circular bacterial chromosome with concentric genomic rings" loading="lazy"></span>
    <span class="plot-card__body"><small>Circular sequence</small><strong>Rings</strong><span>Features, markers, signal and ruler rings around one sequence, with structural chords through the centre.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../guide/maps/#draw-locations-and-explicit-links">
    <span class="plot-card__media"><img src="../../assets/figures/example-maps.svg" alt="Locations and flows under three geographic projections" loading="lazy"></span>
    <span class="plot-card__body"><small>Supplied geography</small><strong>Map</strong><span>Locations, categories, counts and explicit flows under equirectangular, Mercator or orthographic projection.</span><b>Open guide <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose the coordinate system

| Context | Coordinate system | Start with |
|:--|:--|:--|
| one visible locus inside a chromosome | whole-chromosome overview beside genomic position | `IdeogramTrack` |
| many contigs on one figure | concatenated assembly position | `GenomeTrack` |
| one circular replicon | clockwise angle around sequence length | `Rings` |
| named latitude and longitude | explicit geographic projection | `Map` |

## Related guides

- [Figure composition](../guide/figure.md) explains `Genome`, `Rings` and `Panels`.
- [Geographic genomics](../guide/maps.md) covers projections, flows and circular phylogeography.
- [Comparison and alignments](comparisons-alignments.md) covers genome-to-genome correspondence rather than global context.
