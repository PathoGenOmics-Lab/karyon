# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.12.0] - 2026-07-31

### Added

- `SnpTrack` and `SnpSite`, a panel of the columns that vary and nothing else.
  An alignment of closely related genomes is almost entirely agreement: thirty
  kilobases carrying thirty-four differences would spend 99.9% of a plot on the
  part that says nothing. `SnpTrack::from_alignment` finds the variable columns
  against a reference row and drops the rest.

  The idea is the one [snipit](https://github.com/aineniamh/snipit) is built
  around. The implementation and the drawing are this crate's own.

- `SvgWriter::text_rotated`, since a five digit position has to stand on end to
  fit under a column narrower than it is.

### Notes

- The x axis of a `SnpTrack` is not linear in the genome. Two adjacent columns
  may be nine bases apart or nine kilobases apart, and the spacing says nothing
  about which. Every column therefore carries its own position, and an
  `AxisTrack` does not belong underneath one: a ruler there would be a lie.
- A cell matching the reference is drawn as a quiet bar rather than a letter,
  on the same reasoning as the alignment track: the matches are the noise.

## [0.11.0] - 2026-07-31

### Added

- `MsaTrack` and `MsaSequence`, a multiple sequence alignment drawn row by row.
  The coordinates are alignment columns rather than genomic positions, which
  are two different things, so the figure's region is the column space and the
  ruler counts columns. Ungapping a row back to reference coordinates is a real
  operation with real decisions in it and is not done silently.
- `MsaDisplay`, defaulting to `Differences`. In a real alignment most cells
  agree and the agreement is the noise, so rows are drawn as a quiet bar and
  only what disagrees with the comparison row is painted. `MsaTrack::compare_to`
  names a reference row; without one the consensus is used.
- `MsaTrack::consensus`, whose ties break by which row came first so that the
  same alignment always gives the same answer.
- `MsaColoring` and `ResidueClass`, six physicochemical classes for proteins,
  six being how many hues the validated palette has. Cysteine goes with the
  hydrophobics, histidine with the positives, tyrosine with the polars, and
  glycine and proline keep their own.
- `MsaTrack::max_rows`, with the sequences that do not fit counted and reported
  on the figure rather than dropped quietly.
- Two more example figures, one nucleotide alignment under a conservation logo
  built from the same sequences, and one protein alignment by residue class.

### Notes

- Neighbouring cells of the same colour are merged into one rectangle. Twelve
  rows of a hundred and twenty columns come out as a hundred and twelve
  rectangles rather than one thousand four hundred and forty, which is the
  difference between a figure and a file no viewer will open.

## [0.10.0] - 2026-07-31

The last of the comparison plots on the roadmap: one set of alignment blocks,
drawn two ways.

### Added

- `AlignmentBlock`, one alignment between a stretch of the query and a stretch
  of the target. Both spans ascend and the strand is a flag, which is how PAF
  records it and saves a conversion on the way in.
- `DotplotTrack`, which gives the target the vertical axis. Forward blocks
  climb, reversed ones descend, and a translocation sits off the main diagonal,
  so a rearrangement has a shape rather than a description.
- `SyntenyTrack`, which gives the target a second bar and joins the two with
  ribbons. An inversion becomes a twist. Ribbons are translucent so two
  crossing ones read as two, and every block is also drawn solid on both bars,
  so a thin ribbon still shows exactly what it connects.
- Both take the whole target sequence with `target_length` or a window of it
  with `target_range`, and infer the span from the blocks when given neither.
- Two more example figures: two chromosomes disagreeing in three places, and
  the inversion close up.

## [0.9.0] - 2026-07-31

Two more kinds of plot, chosen because they answer consecutive questions.

### Added

- `ManhattanTrack` and `Association`, association statistics along the sequence
  with a threshold line and hits coloured and ringed above it.
  `Association::from_p_value` converts a p-value to `-log10(p)` and floors a
  zero at the smallest positive double rather than letting it become infinity.
  `genome_wide_threshold` is the conventional `-log10(5e-8)`, documented as a
  Bonferroni correction for a million independent tests and therefore the wrong
  number for most organisms.
- `MatrixTrack`, `MatrixRow` and `CellScale`, a sample by site matrix: a
  genotype matrix out of a VCF, or a pangenome presence matrix once its genes
  have coordinates. Row names go in the axis strip, so they size themselves.
  `CellScale::Sequential` for quantities and presence, `CellScale::Categorical`
  for genotypes, where the numbers name things rather than measure them.
- `assets/example-association.svg`, a Manhattan tower with the gene under it
  and the isolates carrying the haplotype under that, all on one x axis.

### Notes

- A sequential matrix ramp starts a step off the page rather than on it. Three
  things have to look different in a genotype matrix: a sample that does not
  carry the allele, a sample that was never typed, and empty page. A pure
  surface at the bottom of the ramp collapses the first two into the third.
- `MatrixTrack` gives every cell a minimum width, since variant sites are
  points and would otherwise be sub-pixel slivers. A cell's width therefore
  says nothing about how much sequence it covers.

## [0.8.1] - 2026-07-31

### Added

- A value axis on `VariantTrack`, so an allele frequency can be read off the
  plot rather than guessed from the height of a stem. `VariantTrack::show_scale`
  turns it off.

  It appears only when it would mean something. `VariantStyle::Tick` draws every
  mark at full height by design, so it gets no axis; and a track whose variants
  carry no values has no ceiling to label, only a convention that a valueless
  stem is drawn full height. Putting a "1" on that would be inventing a
  measurement. A pinned `max` counts as a ceiling even with no values.

  Frequencies keep their decimals, since a scale that rounded 0.64 to 1 would be
  worse than no scale, while large values still shorten to `1.5k` and `2M`.

## [0.8.0] - 2026-07-31

Quantitative tracks get a real value axis. Until now a coverage track showed
only its maximum, crammed into the top left corner on a translucent patch, and a
reader could not get a number off the plot at any other height.

### Added

- `Track::y_axis_width`, which a track implements to ask for a strip between
  the labels and the plotting area, and `DrawContext::axis`, the rectangle it
  gets. The figure reserves the widest request across every track, so the
  plotting areas still line up and the shared x axis survives. A track that
  wants no axis returns zero, gets a zero-width strip, and cannot draw outside
  itself: the clip covers the band and the strip together.
- A value axis on `CoverageTrack`: the ceiling and zero, labelled outside the
  plot, with a hairline across it at the top of the scale. Two ticks rather
  than a ladder, because a coverage profile is read for its shape and its order
  of magnitude and six gridlines would be more ink than the thing they measure.
- The same for `LogoTrack`, whose top and bottom values move out of the band.

### Changed

- Track labels sit to the left of the value axes, so a track with an axis and
  one without still line their names up.

## [0.7.0] - 2026-07-31

### Added

- `IdeogramTrack`, a chromosome drawn end to end with a marker showing which
  part of it the rest of the figure is showing. It is the one track that does
  not map its data through the shared `Scale`, because a track that only showed
  the region on display could not say where that region is.
- `Band` and `Stain`, with `Stain::from_name` reading the `gieStain` column of
  the UCSC `cytoBand` table so a row of that file converts without a lookup
  table of your own. An unrecognised stain becomes the palest shade rather than
  an invented dark one.
- `IdeogramTrack::bare` for a chromosome with no banding, which is what a
  bacterial genome wants, and `IdeogramTrack::highlight` for marking a span
  other than the one on display.
- `theme::mix`, which blends two colours. The cytogenetic grey ladder is mixed
  from the theme's own ink and page, so a dark figure gets a dark ladder rather
  than a light one left uninverted.
- `SvgWriter::rect_outline`, `SvgWriter::path_stroked` and
  `SvgWriter::begin_clip_path`, the last so bands can be painted inside an
  arm's silhouette and keep its rounded end.
- Two more example figures: a banded chromosome above its detail tracks, and
  the H37Rv chromosome as a bare outline with rpoB marked on it.

## [0.6.0] - 2026-07-31

A pass over how the figures look, with the colour half decided by measurement
rather than taste.

### Fixed

- The categorical palette had two slots, a brown and a grey, that were **1.8
  apart under protanopia** and 7.8 apart even with full colour vision. Any
  figure with seven or eight categories was drawing two of them the same. The
  palette is now six hues, every pair of which clears the separation floor.
- The dark palette failed on every count: half of it sat outside the lightness
  band a dark surface needs, four entries were desaturated enough to read as
  grey, and two were 10.4 apart. It has been replaced with steps chosen against
  the dark surface and validated there.

### Added

- `BaseColors::colorblind_safe`, four nucleotide colours whose closest pair is
  11.0 apart under deuteranopia. The conventional IGV colours remain the
  default and are now named `BaseColors::conventional`, with their measured
  weakness stated: adenine and guanine are 1.7 apart under protanopia, which is
  the transition pair.
- `Theme::corner_radius` and `Theme::insertion`, the latter so a pileup's
  insertion marks stop being a hardcoded hex.
- `SvgWriter::rect_rounded` and `SvgWriter::circle_ringed`.
- `assets/example-dark.svg`, because a theme nobody has looked at is a theme
  nobody has checked.

### Changed

- A coverage area is now a wash under a drawn line rather than a saturated
  block: the line carries the shape and the fill only says which side of it is
  under the curve. Bars stay solid, since there the bar is the mark.
- Coverage fits its axis with a little headroom, so the tallest point reads as
  a peak instead of something that ran out of band. A pinned maximum is still
  taken literally.
- Variant markers are large enough to see and carry a ring in the page colour,
  so two variants a base apart read as two variants rather than one blob.
- Features and reads have rounded ends.
- Lines are 2 pixels rather than 1.2.

## [0.5.0] - 2026-07-31

### Added

- `PileupTrack`, aligned reads stacked into rows the way a genome browser
  stacks them, with mismatches painted against the reference.
- `CigarOp` and `Read`, which model an alignment rather than an interval. Every
  SAM operation has somewhere to go and consumes what the specification says it
  consumes, so converting a BAM record is one arm per letter. `M`, `=` and `X`
  all become `Match`: the track compares the sequences itself rather than
  trusting the operation.
- `Read::segments`, which walks the CIGAR and places every piece on the
  reference, and `Read::base_at`, which finds the read base at a reference
  position through whatever insertions and deletions lie in between.
- `PileupTrack::max_rows`, defaulting to forty, with the reads that do not fit
  counted and reported on the figure instead of dropped quietly.
  `PileupTrack::layout` returns that count.
- `ReadColoring::Strand` for strand bias at a glance, and
  `PileupTrack::fade_by_quality`, because a read with mapping quality zero
  could have come from anywhere and should not be drawn as solidly as one that
  could not.
- A fifth example figure, `assets/example-pileup.svg`, whose depth profile is
  computed from the same reads the pileup draws.

## [0.4.0] - 2026-07-31

Empirical Bayes stabilisation, the other half of the Logolas paper. A logo can
now tell four sequences from four thousand.

### Added

- `karyon::dash`, an implementation of Dirichlet adaptive shrinkage. Each
  column of counts is modelled as multinomial with a prior that is a mixture of
  Dirichlet distributions, all centred on a background and differing in
  concentration. The mixture weights are fitted by EM across every column at
  once, so a well sampled column overrules the prior while a thin one is pulled
  towards the background. From Dey, Xie and Stephens, BMC Bioinformatics 19:473
  (2018).
- `LogoTrack::stabilize` and `LogoTrack::stabilize_with` to apply it, and
  `LogoTrack::sample_size` for the case where the input is a probability matrix
  and carries no record of how many sequences it came from.
- `LogoTrack::dash_fit`, returning the fitted mixture weights and, through
  `DashFit::shrinkage`, how far each column actually moved. A figure that
  quietly moved the numbers should be able to say by how much.
- `karyon::dash::ln_gamma`, since the standard library has none and this crate
  has no dependencies.
- A fourth example figure, `assets/example-logo-stability.svg`, showing one set
  of proportions at three sample sizes drawn raw and shrunk.

### Changed

- `LogoTrack::smoothing` is not applied to a stabilised track. The shrunk
  composition is already strictly positive, so the smoothing would be a second
  repair on top of a first.

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

[0.12.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.12.0
[0.11.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.11.0
[0.10.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.10.0
[0.9.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.9.0
[0.8.1]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.8.1
[0.8.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.8.0
[0.7.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.7.0
[0.6.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.6.0
[0.5.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.5.0
[0.4.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.4.0
[0.3.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.3.0
[0.2.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.2.0
[0.1.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.1.0
