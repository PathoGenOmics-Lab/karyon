---
title: Annotation and coordinate plots
description: Features, transcription units, reading frames, coordinate rulers, codon rulers and legends.
---

<div class="plot-hero plot-hero--annotation" markdown>

<span class="plot-eyebrow">Plot catalogue · Annotation and coordinates</span>

# Name the structures and the scale

These tracks turn coordinates into biological objects a reader can point at:
genes, transcripts, open reading frames, codons, base positions and the keys
needed to decode the rest of a figure.

<div class="plot-stats"><span><strong>6</strong> track types</span><span><strong>2</strong> coordinate rulers</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#featuretrack">
    <span class="plot-card__media"><img src="../../assets/figures/example.svg" alt="Gene and interval annotations under a genomic signal" loading="lazy"></span>
    <span class="plot-card__body"><small>Intervals</small><strong>FeatureTrack</strong><span>Genes, exons, repeats or primers packed into non-overlapping rows with strand-aware shapes.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#transcriptionunittrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-transcripts.svg" alt="Transcription units with starts, leaders and terminators" loading="lazy"></span>
    <span class="plot-card__body"><small>RNA span</small><strong>TranscriptionUnitTrack</strong><span>Start site, 5′ leader and terminator shown as one transcriptional claim rather than separate features.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#orftrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-frames.svg" alt="Six reading frames with stops and open stretches" loading="lazy"></span>
    <span class="plot-card__body"><small>Six reading frames</small><strong>OrfTrack</strong><span>Forward and reverse open stretches partitioned by stop codons without pretending they are annotations.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#axistrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-visual-system.svg" alt="Aligned figures with shared quantitative and coordinate axes" loading="lazy"></span>
    <span class="plot-card__body"><small>Base coordinates</small><strong>AxisTrack</strong><span>A 1-based human-readable ruler with round steps and one consistent bp, kb or Mb unit.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#codontrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-codons.svg" alt="Numbered codons and translated residues under variant calls" loading="lazy"></span>
    <span class="plot-card__body"><small>Protein coordinates</small><strong>CodonTrack</strong><span>A strand-aware codon ruler that lets a genomic figure be referenced as S450L or V600E.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#legendtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-visual-system.svg" alt="A manuscript sheet using shared categorical symbols and legends" loading="lazy"></span>
    <span class="plot-card__body"><small>Visual key</small><strong>LegendTrack</strong><span>A stackable legend band that wraps instead of covering data or dropping entries.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose the biological object

| You need to show… | Start with |
|:--|:--|
| an interval that exists independently | `FeatureTrack` |
| one RNA molecule from start to termination | `TranscriptionUnitTrack` |
| possible coding stretches in unannotated sequence | `OrfTrack` |
| genomic positions a reader can type elsewhere | `AxisTrack` |
| residue numbering and codon boundaries | `CodonTrack` |
| the meaning of colours, symbols or line styles | `LegendTrack` |

## Related routes

- [Signal and sequence](signal-sequence.md) for the measurements annotation sits beside.
- [Variation and association](variation-association.md) for calls that land inside genes or codons.
- [Whole genomes and geography](whole-genomes-geography.md) when the context is an assembly or complete chromosome.
