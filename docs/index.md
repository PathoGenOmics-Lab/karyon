---
template: home.html
title: Genome figures where every row lines up
description: >-
  karyon is a Rust library and command line program that stacks tracks such as
  coverage, genes, variants, reads and phylogenies over one shared coordinate
  axis and writes a standalone SVG.
hide:
  - navigation
  - toc
---

<section class="kh-section" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">One command, one figure</p>

## Name a place, then your files

Each file is a row, and every row is drawn on the same coordinates, so the
gene, the depth over it and the calls inside it line up without anyone placing
them.
</div>

<figure class="k-plate" markdown>
![A stack of four rows over two kilobases of the rpoB locus: a depth profile with a dropout in it, a reference row that says to zoom in to see bases, the gene with its resistance determining region marked inside it, variant lollipops coloured by consequence, and a coordinate ruler underneath](assets/figures/example.svg){ width="900" height="305" }
</figure>

```bash
karyon NC_000962.3:761,000-762,999 \
  depth.bg --aggregate min reference.fa annotation.gff3 variants.vcf \
  --title 'rpoB locus, resistance determining region' \
  -o example.svg
```

</section>

<section class="kh-section" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">What do you have?</p>

## Pick your data
</div>

<div class="k-plots kh-cards" markdown>

-   [![Reads, calls and genes](assets/start/reads.svg){ .k-light width="720" height="248" loading="lazy" }![Reads, calls and genes](assets/start/reads-dark.svg){ .k-dark width="720" height="248" loading="lazy" }](your-data/reads.md)

    **[Reads, calls and genes](your-data/reads.md)**
    A BAM, a VCF and an annotation over one place.

-   [![An association scan](assets/start/scan.svg){ .k-light width="720" height="185" loading="lazy" }![An association scan](assets/start/scan-dark.svg){ .k-dark width="720" height="185" loading="lazy" }](your-data/scan.md)

    **[An association scan](your-data/scan.md)**
    A table from PLINK, REGENIE, SAIGE or another association tool.

-   [![A tree and its samples](assets/start/tree.svg){ .k-light width="720" height="717" loading="lazy" }![A tree and its samples](assets/start/tree-dark.svg){ .k-dark width="720" height="717" loading="lazy" }](your-data/tree.md)

    **[A tree and its samples](your-data/tree.md)**
    A Newick tree and a sheet of what you know about each sample.

-   [![An alignment and its tree](assets/start/alignment.svg){ .k-light width="720" height="642" loading="lazy" }![An alignment and its tree](assets/start/alignment-dark.svg){ .k-dark width="720" height="642" loading="lazy" }](your-data/alignment.md)

    **[An alignment and its tree](your-data/alignment.md)**
    An aligned FASTA, in the order of a tree drawn beside it.

-   [![Two assemblies](assets/start/assemblies.svg){ .k-light width="720" height="234" loading="lazy" }![Two assemblies](assets/start/assemblies-dark.svg){ .k-dark width="720" height="234" loading="lazy" }](your-data/assemblies.md)

    **[Two assemblies](your-data/assemblies.md)**
    How two assemblies line up, from a PAF alignment.

-   [![A whole sequence](assets/start/genome.svg){ .k-light width="720" height="227" loading="lazy" }![A whole sequence](assets/start/genome-dark.svg){ .k-dark width="720" height="227" loading="lazy" }](your-data/genome.md)

    **[A whole sequence](your-data/genome.md)**
    The depth of one sample or several along a whole chromosome.

</div>

<p class="kh-cards--more" markdown>[Every kind of figure](plots/index.md){ .md-button }</p>

</section>

<section class="kh-section k-live" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">Running in this page</p>

## Try it here

This is karyon itself, compiled to WebAssembly. Drag the figure along the
genome, zoom with the buttons or the `+` and `-` keys, and watch the command
above it change: every frame is a new run of it.
</div>

<div class="k-stage">
  <div class="k-stage-head">
    <span class="k-stage-name" data-karyon-name>karyon, drawn in advance</span>
    <span class="k-stage-hint" data-karyon-hint></span>
    <span class="k-stage-keys">
      <button type="button" class="k-key" data-karyon-out disabled>Zoom out</button>
      <button type="button" class="k-key" data-karyon-in disabled>Zoom in</button>
      <button type="button" class="k-key" data-karyon-reset disabled>Reset</button>
    </span>
  </div>
  <pre class="k-stage-command" data-karyon-command><code>NC_000962.3:761,000-762,999 \
  --coverage depth.bg --label depth --aggregate min \
  --features genes.gff3 --label annotation \
  --variants calls.vcf --label variants \
  --title 'rpoB locus, resistance determining region'</code></pre>
  <div class="k-stage-plot" data-karyon-plot>
    <img src="assets/figures/example-live.svg" alt="A stack of three rows over two kilobases: a depth profile with a dropout in it, the rpoB gene running off both edges with its resistance determining region on the row beneath it, and variant lollipops coloured by consequence, over a coordinate ruler" width="860" height="265" loading="lazy">
  </div>
  <p class="k-stage-status" data-karyon-status aria-live="polite">Drawn in advance from the command above. It becomes interactive once the program has loaded.</p>
</div>

<details class="k-files" markdown>
<summary>The three input files</summary>

These are the files the command reads, written out in full. Move the window
somewhere none of them has data and the command answers the way it would in a
terminal, for example `no variants in the region`, instead of drawing an empty
figure.

<div class="k-file" markdown>
<p class="k-file-name">depth.bg</p>
<pre data-karyon-file="depth.bg"><code>NC_000962.3 756999 759999 62
NC_000962.3 759999 760999 58
NC_000962.3 760999 761899 57
NC_000962.3 761899 762029 3
NC_000962.3 762029 763999 60
NC_000962.3 763999 766999 54</code></pre>
</div>

<div class="k-file" markdown>
<p class="k-file-name">genes.gff3</p>
<pre data-karyon-file="genes.gff3"><code>##gff-version 3
NC_000962.3 . gene 759807 763325 . + . Name=rpoB
NC_000962.3 . gene 761082 761162 . + . Name=RRDR</code></pre>
</div>

<div class="k-file" markdown>
<p class="k-file-name">calls.vcf</p>
<pre data-karyon-file="calls.vcf"><code>NC_000962.3 760106 . C T . . AF=0.09;ANN=T|synonymous_variant|LOW|rpoB
NC_000962.3 761052 . C T . . AF=0.12;ANN=T|synonymous_variant|LOW|rpoB
NC_000962.3 761109 . G T . . AF=0.98;ANN=T|missense_variant|MODERATE|rpoB
NC_000962.3 761139 . C T . . AF=0.55;ANN=T|missense_variant|MODERATE|rpoB
NC_000962.3 761155 . T C . . AF=1.00;ANN=C|missense_variant|MODERATE|rpoB
NC_000962.3 761156 . C T . . AF=0.21;ANN=T|synonymous_variant|LOW|rpoB
NC_000962.3 761606 . G A . . AF=0.07;ANN=A|synonymous_variant|LOW|rpoB
NC_000962.3 762206 . C T . . AF=0.15;ANN=T|synonymous_variant|LOW|rpoB</code></pre>
</div>

</details>


</section>



<section class="kh-cite" markdown>

Using karyon in a paper? Cite the version you used: the
[citation page](about/citation.md) has the reference. Built at I²SysBio,
University of Valencia-CSIC, and the FISABIO Joint Research Unit Infection and
Public Health, Valencia, Spain.

</section>

<script src="assets/karyon-wasm.js" defer></script>
<script src="assets/karyon-live.js" defer></script>
