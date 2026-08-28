---
title: Signal and sequence plots
description: Coverage, signed windows, sequence, per-site methylation and sequence logos.
---

<div class="plot-hero plot-hero--signal" markdown>

<span class="plot-eyebrow">Plot catalogue · Signal and sequence</span>

# Values and symbols along a sequence

Use these tracks when each genomic position, window or aligned column carries a
measurement or symbol. The important choice is whether zero is a floor, a
baseline to cross, or not part of the scale at all.

<div class="plot-stats"><span><strong>5</strong> track types</span><span><strong>3</strong> scale contracts</span></div>

</div>

<div class="plot-card-grid">
  <a class="plot-card" href="../../tracks/#coveragetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example.svg" alt="Coverage profile over a locus" loading="lazy"></span>
    <span class="plot-card__body"><small>Dense or sparse signal</small><strong>CoverageTrack</strong><span>Per-base depth, GC content or mappability reduced honestly when several bases share a pixel.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#windowtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-selection.svg" alt="Signed statistics in genomic windows" loading="lazy"></span>
    <span class="plot-card__body"><small>Signed window statistic</small><strong>WindowTrack</strong><span>pN/pS, GC skew or Tajima's D drawn on both sides of the baseline they are interpreted against.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#methylationtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-methylation.svg" alt="Per-site methylation on forward and reverse strands" loading="lazy"></span>
    <span class="plot-card__body"><small>Per-site fraction</small><strong>MethylationTrack</strong><span>Forward and reverse methylation calls kept separate, with depth filtering and hemimethylation queries.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#sequencetrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-zoom.svg" alt="Reference bases rendered as coloured letters" loading="lazy"></span>
    <span class="plot-card__body"><small>Reference symbols</small><strong>SequenceTrack</strong><span>Bases become letters, blocks or a zoom hint according to the actual pixel resolution.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#dynseqtrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-regulation.svg" alt="Per-base attribution drawn as scaled letters over a promoter" loading="lazy"></span>
    <span class="plot-card__body"><small>A signed number per base</small><strong>DynseqTrack</strong><span>Model attribution drawn as the bases themselves, hanging below the line where the model pulled away.</span><b>Open reference <span aria-hidden="true">&rarr;</span></b></span>
  </a>
  <a class="plot-card" href="../../tracks/#logotrack">
    <span class="plot-card__media"><img src="../../assets/figures/example-logo.svg" alt="One DNA motif under three sequence-logo scores" loading="lazy"></span>
    <span class="plot-card__body"><small>Aligned-column composition</small><strong>LogoTrack</strong><span>Probabilities, information, enrichment or depletion with explicit alphabet and background contracts.</span><b>Open reference <span aria-hidden="true">→</span></b></span>
  </a>
</div>

## Choose the baseline first

| If the value… | Start with | Why |
|:--|:--|:--|
| cannot be negative and zero is the floor | `CoverageTrack` | Area, line or bars rise from a meaningful zero. |
| is interpreted relative to zero, one or another reference | `WindowTrack` | Both sides of the baseline remain visible. |
| belongs to one strand at a named site | `MethylationTrack` | Strand and read support remain part of the datum. |
| is a symbol rather than a number | `SequenceTrack` or `LogoTrack` | The renderer changes representation with resolution or column composition. |

## Related routes

- [Reads and molecules](reads-molecules.md) for per-read methylation and raw nanopore current.
- [Variation and association](variation-association.md) when a position carries a call rather than a signal.
- [Annotation and coordinates](annotation-coordinates.md) for features and rulers under the signal.
