<div align="center">
  <h1>karyon</h1>
  <p><strong>Genomic track plots for Rust. Composable tracks over a shared coordinate axis, rendered to standalone SVG.</strong></p>

  <p>
    <a href="https://github.com/PathoGenOmics-Lab/karyon/actions/workflows/ci.yml"><img src="https://github.com/PathoGenOmics-Lab/karyon/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
    <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
    <img src="https://img.shields.io/badge/Rust-1.74%2B-orange.svg" alt="Rust 1.74+">
    <img src="https://img.shields.io/badge/dependencies-0-brightgreen.svg" alt="Zero dependencies">
    <a href="https://github.com/PathoGenOmics-Lab"><img src="https://img.shields.io/badge/PathoGenOmics-lab-red?" alt="PathoGenOmics Lab"></a>
  </p>
</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>

General plotting libraries know about points and lines. They do not know that a
position is a base, that a gene has a strand, that a pixel at genome scale
covers two thousand bases, or that a figure is worthless if its tracks do not
line up. `karyon` is the small amount of code that does know those things.

It draws what a genome browser draws: a stack of tracks over one shared
coordinate axis, so read depth, the reference bases, the gene models and the
variant calls all agree on where position 761,410 is.

<img src="assets/example.svg" alt="A coverage profile with a dropout, the reference sequence, two gene models and variants coloured by consequence, all over one coordinate axis" width="100%">

Zoom in and the same tracks show individual bases. Nothing about the tracks
changes, only the region:

<img src="assets/example-zoom.svg" alt="The same locus at base resolution, with the reference sequence drawn as coloured letters" width="100%">

## Quick start

Not on crates.io yet, so point Cargo at the repository:

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

```rust
use karyon::{AxisTrack, CoverageTrack, Feature, FeatureTrack, Figure, Region, Strand};

let region = Region::parse("NC_000962.3:761000-763000")?;

Figure::new(region)
    .title("rpoB locus")
    .push(CoverageTrack::new(760_999, depth).label("depth"))
    .push(FeatureTrack::new(vec![
        Feature::new(761_050, 762_100).name("rpoB").strand(Strand::Forward),
    ]).label("genes"))
    .push(AxisTrack::new())
    .save_svg("rpoB.svg")?;
```

Run the figures above yourself:

```bash
cargo run --example locus -- assets
```

## Tracks

| Track | What it draws | Notes |
|:------|:--------------|:------|
| `CoverageTrack` | Per-base signal: depth, GC content, mappability | Area, line or bars. Bins to one point per pixel with max, mean or min. Optional log scale |
| `SequenceTrack` | The reference bases | Letters when zoomed in, coloured blocks when not, a hint when the bases are thinner than a pixel |
| `FeatureTrack` | Genes, exons, repeats, primers | Strand arrows, automatic packing into rows so nothing overlaps, labels inside or beside |
| `VariantTrack` | SNPs, indels, any point event | Lollipops scaled by value, or ticks when dense. Coloured and legended by category |
| `LogoTrack` | Sequence logos | Probability, information content, or enrichment and depletion around a baseline. Arbitrary alphabets |
| `AxisTrack` | The coordinate ruler | Round tick positions, one unit for the whole ruler, bp, kb or Mb as the zoom demands |

Each is an implementation of one trait with no privileged access to the figure.
A track type that is not here is about thirty lines: see the example on
[`Track`](src/track/mod.rs).

## Sequence logos

A classic logo can only say "this symbol is common here". It measures
conservation, so a column that is nearly uniform comes out flat, and the
biology hiding in that column stays hidden.

`LogoTrack` draws the classic logo and the alternative. In
`EnrichmentDepletion` scaling, each symbol scores `log2(p / q)` against a
background and the column is recentred on its own median, so enriched symbols
stack above the line and depleted ones hang below it. This is the plot
[Logolas](https://github.com/kkdey/Logolas) calls an EDLogo.

<img src="assets/example-logo.svg" alt="The same eight column motif drawn three ways: as probabilities, as information content in bits, and as enrichment above a line with depletion below it" width="100%">

Look at position 4. It is near uniform, so the bits panel says there is nothing
there, and it has no T at all, which only the third panel can tell you.

```rust
use karyon::{Figure, LogoScaling, LogoTrack, Region};

let logo = LogoTrack::from_sequences(0, &alignment)
    .alphabet_size(4)                              // count the bases that never appear
    .scaling(LogoScaling::EnrichmentDepletion)
    .background([("A", 0.35), ("C", 0.15), ("G", 0.15), ("T", 0.35)])
    .label("motif");

Figure::new(Region::new("motif", 0, 8)?)
    .push(logo)
    .save_svg("motif.svg")?;
```

Symbols are arbitrary strings, so an alphabet is not limited to four letters or
to one character each:

<img src="assets/example-logo-protein.svg" alt="A sequence logo whose symbols are three letter amino acid codes" width="70%">

Two things are worth knowing before you read one of these. Absent symbols
dominate the scale, because absence really is a strong signal, and how far they
fall is set by `LogoTrack::smoothing`. That smoothing is a fraction of the
column mass rather than a pseudocount, so the same motif plots identically
whether you pass counts out of 500 or a probability matrix.

## Coordinates

Positions are **0-based and half-open** everywhere, the BED convention. The two
exceptions are the ones a reader sees: `Region::parse` accepts the 1-based
inclusive locus strings that samtools and IGV use, and tick labels are printed
in that same form.

```rust
let region = Region::parse("chr1:101-200")?; // what you type in IGV
assert_eq!(region.start(), 100);             // 0-based
assert_eq!(region.end(), 200);               // exclusive
assert_eq!(region.len(), 100);
```

A VCF `POS` or a GFF `start` is therefore `pos - 1` on the way in. This is the
single most common source of off-by-one bugs in genomic plotting, so the crate
states the convention on every constructor rather than leaving it to be
guessed.

## Design

**Zero runtime dependencies.** The crate compiles in under a second and adds
nothing to your dependency tree. The SVG writer is a few hundred lines because
that is all this needs.

**Scale aware by construction.** A coverage track over four megabases does not
emit four million points; it bins to one value per pixel column, choosing max,
mean or min. A sequence track past single-base resolution draws a hint instead
of a million rectangles. A 4 Mb genome-wide figure comes out under 100 KB, which
is checked by a test.

**Output that survives publication.** Plain SVG 1.1, no scripts, no external
references, no embedded fonts. It opens unchanged in a browser, in Inkscape and
in Illustrator, and every element stays selectable when a reviewer asks for the
gene labels to be bigger.

**Deterministic.** The same input renders byte-identical output, and category
colours follow first appearance rather than hash order. A figure that recolours
itself when a sample is added is not one you can put in a paper.

**No file parsing.** `karyon` takes vectors of numbers and structs, not paths to
BAM files. Reading genomic formats is a solved problem that
[noodles](https://github.com/zaeleus/noodles) and
[rust-bio](https://github.com/rust-bio/rust-bio) solve better, and keeping them
out is what makes the dependency count zero.

## Roadmap

Not implemented yet, in the order they are likely to arrive:

- Read pileup track, with mismatches coloured against the reference
- Ideogram and karyogram, for whole-chromosome context
- Dotplot and synteny ribbons between two sequences
- Manhattan plot, with its own y axis and significance line
- More logo scorings: KL divergence and the unscaled log odds Logolas offers
- PNG output, likely behind a feature flag so the default stays dependency-free

## Installation

```bash
git clone https://github.com/PathoGenOmics-Lab/karyon
cd karyon
cargo test
```

Nothing is published to crates.io yet. `cargo add karyon` will work once 0.1.0
is released there.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
