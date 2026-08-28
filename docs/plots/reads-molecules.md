---
title: Reads and molecule plots
description: Aligned reads, split molecules, bisulfite patterns and raw nanopore signal.
---

<div class="plot-hero plot-hero--reads" markdown>

<span class="plot-eyebrow">Plot catalogue · Reads and molecules</span>

# Keep the evidence at molecule resolution

Use these tracks before reads have been collapsed into depth, a variant call or
a site average. Each row or trace remains attributable to the molecule that
produced it.

<div class="plot-stats"><span><strong>4</strong> track types</span><span><strong>1</strong> molecule per row</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#pileuptrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-pileup.svg" alt="Read pileup with strand, mismatches, insertions and deletions" loading="lazy"></span>
    <span class="plot-card__body"><small>Contiguous alignment</small><strong>PileupTrack</strong><span>Real CIGAR operations packed into rows, with strand, mapping quality and reference-aware mismatches.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#splitreadtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-split.svg" alt="One molecule aligned in several genomic segments" loading="lazy"></span>
    <span class="plot-card__body"><small>Segmented alignment</small><strong>SplitReadTrack</strong><span>One row per molecule, one bar per segment and connectors that preserve visit order and orientation.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#bisulfitetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-bisulfite.svg" alt="Single-molecule methylation calls at cytosines" loading="lazy"></span>
    <span class="plot-card__body"><small>Modification pattern</small><strong>BisulfiteTrack</strong><span>Filled and open calls per covered site, leaving truly uncovered positions empty rather than unmodified.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#junctiontrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-regulation.svg" alt="Junction arcs over a depth profile, with per-base attribution under them" loading="lazy"></span>
    <span class="plot-card__body"><small>Introns reads stepped over</small><strong>JunctionTrack</strong><span>Sashimi arcs weighted and labelled by the reads that crossed each junction.</span><b>Open reference <span aria-hidden="true">&rarr;</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#squiggletrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-squiggle.svg" alt="Raw nanopore current resolving from an envelope into a trace" loading="lazy"></span>
    <span class="plot-card__body"><small>Raw current</small><strong>SquiggleTrack</strong><span>Nanopore signal as a min–max envelope at overview scale and the original trace when resolution allows.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose by what was measured

| Input still contains… | Start with |
|:--|:--|
| one reference position and CIGAR per read | `PileupTrack` |
| several alignments belonging to one molecule | `SplitReadTrack` |
| modified and unmodified calls along each molecule | `BisulfiteTrack` |
| current samples and a basecaller move table | `SquiggleTrack` |

## Related routes

- [Signal and sequence](signal-sequence.md) after reads have become coverage or per-site fractions.
- [Variation and association](variation-association.md) after molecule evidence has become a call.
- [Recipes](../recipes.md) for piping SAM text from `samtools` and attaching a reference in Rust.
