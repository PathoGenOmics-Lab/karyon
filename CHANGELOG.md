# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-07-31

Sequence logos gain the rest of the scoring schemes `Logolas` offers, and the
score is separated from where the baseline sits.

### Added

- Four more scores on top of the existing three: `LogoScore::KullbackLeibler`
  (`p log2(p/q)`, the symbol's contribution to the divergence),
  `LogoScore::Difference` (`p - q`), `LogoScore::Ratio` (`p / q`) and
  `LogoScore::OddsRatio`. Together with `Probability`, `InformationContent` and
  `LogOdds` these cover the `Logolas` set.
- `Centering`, which decides where the baseline of a background-relative logo
  sits, as a quantile of the column rather than a fixed median. `Quantile(0.5)`
  is the default and matches the `quant` argument of `Logolas`;
  `Centering::None` leaves the baseline at "exactly the background".
- `LogoTrack::edlogo`, shorthand for log odds centred on the median.
- `LogoTrack::max_extent`, pinning how far the tallest stack reaches so that
  two panels of a figure can be compared honestly.
- `LogoScore::uses_background` and `LogoScore::unit`.
- A third example figure, `assets/example-logo-scores.svg`, showing the same
  four columns under all five background-relative scores.

### Changed

- `LogoScaling` is now `LogoScore`, and `LogoTrack::scaling` is
  `LogoTrack::score`. The old `EnrichmentDepletion` variant was the pairing of
  one score with one centring, which is now spelled `LogoTrack::edlogo()` or
  `.score(LogoScore::LogOdds).centering(Centering::median())`.
- Smoothing now applies to every score measured against a background rather
  than to the one that used to exist.

### Fixed

- An information content logo labelled its axis `2 bits` while drawing it to
  the tallest column on screen, so the annotation and the picture disagreed.
  The two absolute scores now get the fixed axis their convention asks for:
  probabilities to one, information content to `log2(K)`, the way WebLogo draws
  DNA from zero to two bits. Two figures of two motifs are comparable again.
- A fitted axis now measures only the columns inside the region on display,
  instead of letting an off-screen column set the scale.
- The scale annotation reports the value at the bottom of a two sided logo as
  well as the one at the top.

## [0.2.0] - 2026-07-31

### Added

- `LogoTrack`, a sequence logo over consecutive positions, with three scalings:
  `Probability`, `InformationContent` (the classic `seqLogo` plot in bits) and
  `EnrichmentDepletion`, which scores each symbol as `log2(p / q)` against a
  background, recentres the column on its median, and stacks enriched symbols
  above a baseline with depleted ones below it. That last one is the plot
  `Logolas` calls an EDLogo, and it is the only one of the three that can show
  a symbol which is absent rather than merely uncommon.
- `LogoColumn`, whose symbols are arbitrary strings rather than single letters,
  so amino acid codes, codons and k-mers all plot.
- `LogoTrack::from_sequences` to count an alignment, and
  `LogoTrack::from_matrix` to take a position weight matrix.
- `LogoTrack::background` for a non-uniform background, `alphabet_size` for
  counting symbols that never appear, `smoothing` for how loudly absence
  speaks, `symbol_color` for per-symbol colours, and `order` for which end of
  the stack the tallest symbol goes.
- `SvgWriter::glyph`, which stretches one glyph to an exact box using
  `textLength`, so the renderer does the fitting rather than this crate
  guessing font metrics.
- `AxisTrack::center_on_bases`, putting each tick in the middle of its base
  instead of on its left edge. A ruler marks boundaries, but a logo column
  wants its number underneath it.
- `Theme::cap_height_ratio`, the cap height of the font as a fraction of the
  font size, which is what turns a box height into a font size for a logo.

### Changed

- `Theme` gained a field, so any code building one with a struct literal rather
  than from `Theme::light()` will need updating.

## [0.1.0] - 2026-07-31

First release. Everything below is new.

### Added

- `Figure`, which stacks tracks over one shared coordinate axis, reserves a
  label gutter only when a track asks for one, clips each track to its own
  band, and computes the image height from its contents.
- `Region`, 0-based half-open, with `Region::parse` for the 1-based inclusive
  locus strings used by samtools and IGV, and a `Display` that round-trips
  through it.
- `Scale`, the shared mapping between base positions and pixels, with
  `x`, `x_center`, `pos_at_x` and the two resolution helpers.
- `CoverageTrack` for per-base signal, in area, line or bar style, binned to one
  point per pixel column with max, mean or min, with an optional log scale and
  a pinnable maximum.
- `SequenceTrack` for the reference bases, drawn as letters, as coloured blocks
  or as a zoom hint depending on how many pixels a base gets.
- `FeatureTrack` for annotated intervals, with strand arrows, per-feature
  colours, and packing into as many rows as the current zoom needs so that
  neither features nor their labels overlap.
- `VariantTrack` for point events, as lollipops scaled by value or as ticks,
  coloured by category in order of first appearance, with a legend.
- `AxisTrack`, whose ticks land on round 1-based coordinates and whose labels
  share one unit across the whole ruler.
- `Theme`, with light and dark presets, the Okabe-Ito categorical palette and
  IGV-style nucleotide colours, all of it public and replaceable.
- `Track`, the trait behind all of the above, so a track type the crate does not
  ship can live outside it.
- `examples/locus.rs`, which renders the two figures in the README from a fixed
  seed.

[0.3.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.3.0
[0.2.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.2.0
[0.1.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.1.0
