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
<p class="kh-eyebrow">One region, every track</p>

## Name the region once. Every track lines up on it.

A karyon figure is a stack of tracks over one coordinate axis. The gene, the
depth over it and the variants inside it land on the same position because they
are drawn on the same scale, not because anyone placed them.
</div>

<figure class="k-plate" markdown>
![A stack of four rows over two kilobases of the rpoB locus: a depth profile with a dropout in it, a reference row that says to zoom in to see bases, the gene with its resistance determining region marked inside it, variant lollipops coloured by consequence, and a coordinate ruler underneath](assets/figures/example.svg){ width="900" height="305" }
</figure>

<div class="kh-code" markdown>

=== "Rust"

    ```rust
    use karyon::{plot, Aggregate};

    plot("NC_000962.3:761000-762999")?
        .title("rpoB locus, resistance determining region")
        .add_coverage(depth).label("depth")
        .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
        .add_sequence(bases).label("reference")
        .add_features(genes).label("annotation")
        .add_variants(variants).label("variants")
        .save("example.svg")?;
    ```

=== "Command line"

    ```bash
    karyon NC_000962.3:761,000-762,999 \
      depth.bg --aggregate min reference.fa annotation.gff3 variants.vcf \
      --title 'rpoB locus, resistance determining region' \
      -o example.svg
    ```

</div>

Each file is drawn as what its name says it holds and labelled after itself.
The ruler, the depth axis, the colour key and the accessible title and
description are added for you. The whole program is
[`examples/locus.rs`](https://github.com/PathoGenOmics-Lab/karyon/blob/main/examples/locus.rs),
and the figure it writes is built in
[`examples/figures/locus.rs`](https://github.com/PathoGenOmics-Lab/karyon/blob/main/examples/figures/locus.rs).

</section>

<section class="kh-section" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">36 track types</p>

## What it draws

From a read pileup to a phylogeny around a globe, sorted by the question you
are asking of your data.
</div>

<div class="k-plots kh-cards" markdown>

-   [![A read pileup with mismatches, insertions, deletions and spliced alignments](assets/figures/example-pileup.svg){ width="920" height="473" loading="lazy" }](plots/reads-molecules.md)

    **[Reads and molecules](plots/reads-molecules.md)**
    Pileups, split reads, single molecules.

-   [![A chromosome ideogram with its bands and a highlighted region](assets/figures/example-ideogram.svg){ width="900" height="275" loading="lazy" }](plots/annotation-coordinates.md)

    **[Annotation](plots/annotation-coordinates.md)**
    Ideograms, genes, transcripts.

-   [![Variable sites across isolates, ordered by the phylogeny beside them](assets/figures/example-snps.svg){ width="900" height="388" loading="lazy" }](plots/variation-association.md)

    **[Variation](plots/variation-association.md)**
    Variants, variable sites, scans.

-   [![A dotplot with synteny ribbons between two genomes](assets/figures/example-synteny.svg){ width="900" height="437" loading="lazy" }](plots/comparisons-alignments.md)

    **[Comparisons](plots/comparisons-alignments.md)**
    Dotplots, synteny, alignments.

-   [![Two phylogenies of the same isolates face to face, with the tips that moved linked across](assets/figures/example-tanglegram.svg){ width="760" height="234" loading="lazy" }](plots/phylogeny-clades.md)

    **[Phylogeny](plots/phylogeny-clades.md)**
    Trees, clades, sample traits.

-   [![Sequence logos scored three ways](assets/figures/example-logo.svg){ width="900" height="379" loading="lazy" }](plots/signal-sequence.md)

    **[Signal and sequence](plots/signal-sequence.md)**
    Coverage, logos, methylation.

-   [![Branch rate mixtures, recurrence links and genomic site-wise selection evidence](assets/figures/example-selection-atlas.svg){ width="1508" height="1053" loading="lazy" }](plots/evolution-surveillance.md)

    **[Evolution and surveillance](plots/evolution-surveillance.md)**
    Selection, lineages over time.

-   [![A dated phylogeny drawn around a globe, each tip linked to where it was sampled](assets/figures/example-phylo-map.svg){ width="1916" height="868" loading="lazy" }](plots/whole-genomes-geography.md)

    **[Whole genomes and maps](plots/whole-genomes-geography.md)**
    Circular plots, assemblies, geography.

</div>

<p class="kh-cards--more" markdown>[Open the gallery](plots/index.md){ .md-button } [Every track type](tracks/index.md){ .md-button }</p>

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

<section class="kh-section" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">Why karyon</p>

## Small, strict and honest about the data
</div>

<div class="kh-points" markdown>

<div class="kh-point" markdown>
### Aligned by construction
One scale is handed to every track, so nothing is lined up by hand and nothing
drifts when the region changes.
</div>

<div class="kh-point" markdown>
### Honest about gaps
When the data cannot support a picture, karyon stops with a clear message
instead of drawing something that only looks right.
</div>

<div class="kh-point" markdown>
### Nothing else to install
No runtime dependencies. A Rust toolchain is all it needs, and the same code
runs in this page as WebAssembly.
</div>

<div class="kh-point" markdown>
### Ready to publish
Plain SVG 1.1 with a title and a description a screen reader can use, and the
same input always gives the same file.
</div>

</div>

</section>

<section class="kh-section" markdown>

<div class="kh-head" markdown>
<p class="kh-eyebrow">Three ways in</p>

## From Rust, from the shell or from this browser
</div>

<div class="kh-ways" markdown>

<div class="kh-way" markdown>
### Rust library

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

[The Rust API](guide/plot.md)
</div>

<div class="kh-way" markdown>
### Command line

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
karyon chr1:1-60 --coverage depth.txt -o first.svg
```

[Every flag](guide/cli.md)
</div>

<div class="kh-way" markdown>
### In your browser

The [playground](playground.md) runs the command line on your own files, and
the [tree viewer](tree.md) opens a Newick tree of up to a million tips. Nothing
leaves your machine.
</div>

</div>

</section>

<section class="kh-cite" markdown>

Using karyon in a paper? Cite the version you used: the
[citation page](about/citation.md) has the reference. Built at I²SysBio,
University of Valencia-CSIC, and the FISABIO Joint Research Unit Infection and
Public Health, Valencia, Spain.

</section>

<script src="assets/karyon-wasm.js" defer></script>
<script src="assets/karyon-live.js" defer></script>
