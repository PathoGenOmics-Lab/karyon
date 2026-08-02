# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Every figure now says what it is. The document carries a `<title>`, a `<desc>`,
  `role="img"` and `aria-labelledby`, so a screen reader is given a sentence
  instead of several thousand unnamed rectangles, and a reader whose image did
  not load is given the alt text. `Figure::description` sets that text; without
  it the description is composed only from what the figure knows for certain,
  the region and the labels of its tracks in order, which says what is there
  without claiming anything about what it shows.
- `SvgWriter::begin_titled` and `SvgWriter::describe`, and a tooltip on every
  glyph worth pointing at: a gene gives its name, span and strand, a variant its
  position, consequence and allele fraction, a read its span and CIGAR. A
  `<title>` as the first child of a group is SVG 1.1, so this needs no script,
  fetches nothing, and survives Inkscape and Illustrator.
- What gets named follows the rule the crate already binned by: the thing that
  decides whether a datum is drawn on its own decides whether it is named. A
  coverage bin is the resolution of the output rather than of the input and gets
  nothing; a dense row-major track names the row and not the cell; a glyph under
  a pixel is nothing a pointer can rest on. The worst a figure grew was a third,
  and the genome-wide figure did not grow at all.
- `SvgWriter::begin_titled_inert` for decoration drawn over data. The ideogram's
  region marker is a translucent rectangle written after the bands, so its own
  tooltip took every hover inside the highlighted span and the bands it points
  at were the one part of the chromosome a reader could not interrogate.
- `plot` and `Plot`, a short way to write a figure down. `Figure` takes
  tracks that already exist, which is right when a track is built elsewhere and
  wrong for the ordinary case: a handful of tracks written out once, in place.
  There the plumbing was most of the code. A `Region` to unwrap, a `use` line
  naming every track type, a `push` around each one, an `AxisTrack` to remember,
  and the window start repeated on every track that takes one even though the
  figure already knew it. `plot("NC_000962.3:761,000-763,000")?` starts a stack
  and one `add_` call per track fills it, in the order they stack.
- `Plot` holds the track added last in its type parameter, so `Plot::label`
  names it and `Plot::adjust` hands it to a closure that sees the concrete type.
  Every builder method on the track is therefore in reach, and a name that does
  not exist on it is a compile error rather than a call that goes nowhere. This
  is the reason the type has a parameter at all.
- Two things a plot fills in, both of which can be undone. An `AxisTrack` goes
  on the bottom, since a figure without coordinates along it is rarely what
  anyone meant; `Plot::add_axis` puts it elsewhere and `Plot::remove_axis`
  leaves it out. Tracks carrying an array rather than coordinates start at the
  left edge of the region, with an `_at` form for the ones that do not.
- `Plot::done`, which puts the pending track away without adding another. A
  stack built in a loop or behind a condition needs every arm to have one type,
  and the type parameter that makes `adjust` work is what takes that away.
- `Plot::add_track` and `Plot::add_boxed`, for a track that is already built:
  one an alternative constructor made, one read back before it is drawn, or one
  from outside the crate. This is not a fallback but the ordinary way to write
  a heavily configured track, since the builder chain says it in one line and
  `adjust` would wrap the same calls in a closure.
- A command line front end, `karyon`, in its own binary target. The library
  still reads no files and still has no dependencies: everything the command
  reads is line based text, and the reading lives in `src/bin/karyon/`. The
  grammar is the facade with spaces instead of dots, because `argv` is already
  an ordered list whose later words describe the earlier ones: each `--<track>`
  flag starts a track and the flags after it describe that one, so the order of
  the flags is the order of the stack.
- Twelve tracks have a standard text format to read and are the ones the
  command has. BED, bedGraph, cytoBand, GFF3, VCF, SAM, FASTA, Newick and two
  tables. BAM, CRAM and BCF are not read at all and are not meant to be: they
  come in through a pipe, since `samtools depth`, `samtools view` and
  `bcftools view` already write what these readers take.
- Coordinates are read as each format defines them, which is the one thing in
  the command that would fail silently. BED, bedGraph and cytoBand are 0-based
  half-open; GFF3, VCF, SAM and `samtools depth` are 1-based. Every reader has
  a test that pins a known base through the conversion, and the end coordinate
  is pinned separately from the start.
- `Tree`, `Clade` and `Placement` at the crate root. `TreeTrack`, `CladeTrack`
  and `TanglegramTrack` were re-exported there and all three take a `Tree`,
  which was not.
- `From<Error> for std::io::Error`. A program that draws figures spends its
  errors on writing files, so its functions return `io::Result`, and the one
  call that parses a region was the only thing in such a function that could
  not use `?`. That is a poor reason to reach for `unwrap`, and it was the
  second line of all fourteen chains in the examples.
- `Empty`, `Slot` and `Named`, the plumbing that gives `Plot` its type
  parameter, under `karyon::plot` rather than at the crate root. `Empty` is
  deliberately not a `Track`, so `label` and `adjust` before any track has been
  added are compile errors rather than calls that go nowhere, and a figure
  cannot be asked to draw the empty slot. `Slot` and `Named` are sealed.
- `CodonTrack`, a ruler in codons. A variant in a coding sequence is named by
  residue (BRAF V600E, TP53 R175H, rpoB S450L) and the crate could only speak
  bases, so no figure could be pointed at with the name its result has. It
  partitions a coding sequence into codons, numbers them, and translates them
  where a letter fits. The partition is itself the claim: two changes in one
  codon are competing alleles at one residue, not a double mutant. On the
  reverse strand codon 1 sits at the highest coordinate, which is the whole
  reason this is a track and not a division by three: roughly half the coding
  sequences in any annotation run backwards, and getting their numbering wrong
  is silent.
- `CodonTrack::genetic_code`, since the translation table was hard-coded to NCBI
  table 1. Table 11 gives the same residues so bacteria, archaea and plastids
  were fine, but a mitochondrial or ciliate sequence was quietly translated into
  a plausible protein that was wrong.
- `SplitReadTrack`, `SplitRead` and `SplitSegment`, reads that align in pieces.
  A `PileupTrack` read is one start, one CIGAR and one strand, so a molecule
  that visits three places cannot be written down in it, and a `StructuralTrack`
  arc starts from a finished call, so by the time there is an arc the evidence
  has been summarised away. This is the evidence: one row per molecule, one bar
  per alignment, connectors saying in what order and orientation it went. A
  backward hop dips under the row and a forward one arches over, so a read
  crossing an inversion is a different picture from one crossing a deletion. A
  transposition is three segments and not an arc.
- `CladeTrack` and `CladeBlock`, genomic intervals painted onto a phylogeny. A
  `MatrixTrack` cell is one base wide and cells never merge, so a matrix can
  only say that six samples each carry something here, which is six
  observations. A block whose height covers a clade says one. The track can lie
  in exactly one way and refuses to: rows inside a block that do not carry it
  are cut out, so a paraphyletic set can never pass for a clade.
- `TranscriptionUnitTrack`, `TranscriptionUnit` and `Terminator`. Nothing in the
  crate touched transcription. The span from the bent arrow to the hairpin is
  one RNA molecule, so the genes under it are co-transcribed and a promoter
  mutation upstream changes all of them at once. The 5' leader is drawn hollow,
  which makes its length a distance to measure rather than a number in a
  caption, and a leaderless transcript is a different picture rather than a
  different label.
- `examples/reading.rs` and four assets, plus gallery panels S to V. The sheet is
  twenty-two panels.
- The README track table now lists all twenty-nine track types. It had stopped at
  `AxisTrack` and was fourteen tracks out of date.

### Added

- `tests/properties.rs`, which asserts invariants rather than examples. A test
  written beside the code inherits the blind spot of whoever wrote it, because
  the person writing the assertion is the person who believed the code: a
  reverse-strand span was wrong here from the first commit and every test of it
  agreed with it. An invariant can be checked without knowing the right answer,
  which is the only way to find something nobody thought of. Seven of them so
  far: every generated figure renders a valid document, rendering is
  deterministic, the reported dimensions are the ones in the document, a locus
  string survives a round trip, the scale never goes backwards, changing the
  theme changes no geometry, and a track added below never shortens the figure.
- The figures are generated by a seeded LCG rather than by a property testing
  crate, so the dependency count stays at zero even for the tests, and a
  failure is reproducible from the seed in the message. What it feeds in is the
  point: not plausible depths but `NaN`, both infinities, `f64::MAX`,
  `u64::MAX`, empty vectors, one-base regions and names carrying markup. Ten
  thousand figures take about three seconds.

### Fixed

- Six coordinate overflows, every one found by the property suite on its first
  run and none reachable with ordinary data: `Scale::bounds`, `Read::end` and
  the CIGAR walk, the codon walk in `OrfTrack`, and the sequence track's right
  edge. They are the same defect in six places. A coordinate near the top of
  the range does not survive a trip through `f64`, since `u64::MAX - 10` is not
  representable, rounds to 2^64 and comes back as `u64::MAX`, so anything added
  to it overflows. All of them saturate now, which is what the arithmetic in
  the rest of the crate already did: a panic is what happens in a debug build
  and a wrapped coordinate is what happens in a release one, and the wrapped
  coordinate is the worse of the two because it draws.

Thirty-five defects, found by hunting the crate along seven failure modes and
then handing each claim to an independent reviewer whose job was to refute it.
Seven were refuted. Every one below has a reproduction that runs and a test
that pins the answer.

**Wrong figures.**

- Reverse-frame open reading frames were a whole codon out at both ends. They
  swallowed the stop codon that is supposed to bound them and dropped their
  last codon, so `OrfTrack` reported 94 codons where there were 93 and painted
  the bar over the very tick that ends it, while the three forward lanes beside
  it were right. Frame -k of a sequence walks the codons frame +k of its
  reverse complement walks, and over 5,520 random sequences 5,519 disagreed
  with their own mirror. They now all agree, and the property is a test.
- `MatrixTrack` put its cell half a floor-width to the right of the site it
  stands for, so it did not line up with the Manhattan point above it in the
  same figure.
- `contrast_ink` chose the less legible ink for adenine and guanine. The
  contrast ratios say so, and now the ink follows them.
- `AxisRing` placed its labels from the outer radius alone, so a thin ruler
  printed its coordinates on top of the ring inside it.
- Every track was clipped to the widest value axis in the figure, so a track
  that asked for no axis could paint into the strip another track had reserved.
- A track label longer than the gutter was drawn off the left edge of the
  canvas and lost its first characters; it is ellipsised now.
- `Panels` sized the sheet from the panel widths alone, so a caption could be
  painted over the next column.
- `Panels::share_out` did not return the most level partition among those tying
  on the tallest column.
- `FeatureTrack` packed features that were off screen, so the band height and
  the row of every visible feature depended on data nobody could see.
- `AxisTrack`'s tick loop overflowed at the top of the coordinate range: a
  panic in a debug build, and 503 wrapped-around tick labels in a release one.
- A doubled quote inside a quoted Newick label started a new taxon, inventing
  tips that are not in the tree.
- `Genome::at` did not bound-check, so a coordinate past the end of one
  sequence was drawn inside the next one and counted as mapped. It now accepts
  exactly the positions `Genome::locate` inverts, which is a narrowing: the end
  of a half-open interval is `Genome::offset`, not `at`.

**Data dropped without a word.**

- `CodonTrack` drew nothing at all when the coding sequence length was not a
  multiple of three, or when the trailing partial codon was in view.
- A structural variant whose two breakpoints straddle the whole window was
  dropped, so the largest events were the ones that vanished.
- `Tree::parse_newick` silently discarded every unnamed leaf, so a legal Newick
  tree lost tips and changed topology.
- A UTF-8 byte order mark made the first data row of every column reader
  disappear.
- Track heights were computed against `Theme::default` while the rows were
  packed with the figure's own theme, so any other font size could cut a row
  off the bottom of a band.
- `Rings::link` collapsed the origin-spanning end of a chord to a hairline.
- `--matrix` whose sites all fall outside the region was a silently empty lane
  and exit 0.

**Output that was not valid SVG.**

- `--color` was pasted into an attribute unescaped. That is attribute
  injection, not a malformed file: a value carrying `<` put elements of the
  caller's choosing into the document. Metacharacters are refused now.
- `svg::num` wrote the literal `inf` for finite inputs above `f64::MAX / 1000`,
  which was the one guard the writer had against a non-finite number reaching
  the page.
- `escape` let XML-forbidden control characters through, so a name read from a
  user's file could produce a document no parser accepts.
- A theme font size below 1 made tracks emit a negative `font-size`.
- A non-finite figure height was written as `height="0"`, and on a sheet it
  stacked every later panel on top of it.

**Panics on data.**

- `u64` coordinate arithmetic overflowed across eight constructors.
- A large finite `Figure::width` panicked with a capacity overflow.
- A named origin-spanning feature panicked `FeatureRing` with a subtraction
  overflow.
- A number on the command line unwound instead of returning an error.

**And the gate that should have caught some of this.** The CI step that fails
when the committed figures disagree with the code rendered ten of the sixteen
examples, so six figures could drift unnoticed, and it never checked the copies
the documentation site serves. It renders all sixteen now and checks both.

- A nested figure named itself, which left every panel of a sheet unnameable.
  `<title>` resolves to the innermost element under the pointer, so the title a
  figure writes for itself shadowed the one `Panels` puts on the panel over the
  panel's whole area: all twenty-two panel tooltips in the gallery existed in
  the file and none could be reached. A document rendered with an id prefix is
  being nested, so it no longer names itself; `role="img"` inside another
  `role="img"` also hid its contents from a screen reader rather than describing
  them.
- `SnpTrack` hung its tooltip on the rotated column label rather than on the
  data, so the panel itself answered nothing, and `show_positions(false)` took
  every tooltip in the track with it. The cell is what a pointer lands on, so
  the cell now carries the site, the reference base and what that sample has.
  Only a cell that differs is named: agreement is the state the panel exists to
  see past, the colour already says it, and naming it cost a third of the file.
- `MatrixTrack` counted every site in the track while the drawing loop skipped
  the ones outside the region, so a zoomed figure described cells that were not
  on the page. Both numbers now count over the predicate the drawing uses.
- Four tracks spelled the coordinate conversion themselves as `start + 1`,
  which dropped the guard on a degenerate interval, so a zero-length span came
  out backwards: a clade block at 100 to 100 read `101 to 100`. There is one
  span formatter in the crate now and all of them call it.

### Changed

- Every module doc is at one level. The crate had been written by many hands
  and read like two crates: eighteen modules carried a single line of `//!`,
  among them `figure.rs` at 462 lines of code, `theme.rs`, `codon.rs` and all
  four command line readers, while `bisulfite.rs` carried three named sections
  and the crate root carried eight. A module doc now opens with one summary
  line and an orienting paragraph, and takes named sections in twos or more,
  since exactly one heading is worse than none. Headings state a claim or ask a
  question, the way "The frames are not the strands" does, rather than saying
  Overview or Usage. The spread went from 1 to 106 lines to 17 to 105, with the
  middle half of the crate between 28 and 35.
- A module doc states the decision and the item doc keeps the reasoning, which
  is a rule the levelling pass broke before it held it: seven modules ended up
  restating their own `///` docs, and in two of them the same sentence stood in
  three places counting the documentation site. Those are back to one place
  each.
- The documentation site is one shape too. Every page title is a noun phrase,
  which matters because the title is the navigation label, and thirteen of the
  fifteen pages now close with the `## Next` block the other tools in the lab
  use. The two that do not are the home page, which has the card grid, and the
  changelog, which is an include.
- The examples are no longer all one organism, and three of them were showing
  an assay the organism does not have. The bisulfite panel drew a CpG island on
  a genome whose methylome is 6mA and which has no 5mC at all; it is now the
  human H19/IGF2 imprinting control region, which is where a per-molecule
  bisulfite plot with two populations comes from. The per-strand methylation
  panel claimed Dam and GATC in a genome with no `dam`; it is now *E. coli*
  K-12 at oriC, with the eleven real GATC sites SeqA holds hemimethylated. The
  pangenome matrix drew an accessory genome with gained and lost islands in a
  species that has neither, and is now *Klebsiella pneumoniae*. The clade panel
  needed a block that was lost more than once, which a clonal organism does not
  supply, so it is SARS-CoV-2 lineage-defining deletions with the nsp6 SGF
  deletion as the recurrent one.
- Coordinates in the remaining figures now match the annotation they name.
  `rpoB` was drawn 1,050 bases long instead of 3,519 with an invented
  neighbouring gene beside it; the RRDR box sat on codons 525 to 572 rather
  than 426 to 452, contradicting the frame `CodonTrack` pins in its own
  doctest; `RD1` was drawn at 8 kb, on top of *gyrA* and *gyrB*, rather than at
  4.35 Mb; `eccCa1` was drawn on the reverse strand when its locus tag has no
  `c` suffix. The split-read panel asserted an empty donor site and twice the
  depth over that same donor, which no single sample can show, and put the
  element in two orientations for one insertion event.
- Documentation throughout now states the general case rather than the bacterial
  or tuberculosis one. Twenty-three sentences across the API docs and the README
  said what the library is for in terms of one clade: circular plotting existed
  because "a bacterial chromosome has no ends", protein coordinates because of
  what is "clinically interesting about a bacterial genome", split reads for
  `IS6110`. The crate is organism-agnostic and the prose now reads that way.
  Example and doctest data are untouched, since a figure has to be of something
  and rpoB at H37Rv coordinates is as good as anything.

### Removed

- `AccumulationTrack`, `AccumulationCurve`, `DistanceTrack`, `FrequencyTrack`
  and `Frequency`, with their examples, assets and gallery panels. All three
  tracks failed the only entry test this crate has: none of their `draw`
  methods ever reads `ctx.scale`. Their x axis was a count of genomes, a count
  of genomes and a list of sample names, so what they actually drew was a bar
  chart, a line chart with a quantile ribbon and a clustered heatmap. A general
  plotting library draws those better, and `src/track/distance.rs` had said so
  in its own header since the day it landed: "Nothing genomic. Both axes of
  this track are the sample list, so it does not use the figure's shared
  coordinate system."
- Breaking, as the names were public re-exports. The analyses they carried are
  worth keeping somewhere, but a rarefaction over a presence matrix is a
  statistic, not a plot type, and it does not need a `Track` to compute it.

### Fixed

- `MethylationTrack::hemimethylated` only paired the two strands when both
  calls sat on the same coordinate, which is never true of a palindrome: the
  two 6mA of a `GATC` are one base apart, and so are the two 5mC of a `CpG`.
  The API forced the caller to lie about one of the two positions, and the
  *E. coli* oriC panel found it the moment its coordinates became real. It now
  pairs the nearest partner within `MethylationTrack::pair_within`, one base by
  default.
- `Panels::columns` now cuts the panels where the tallest column comes out
  shortest, and levels the rest, instead of filling each column until it has
  passed its share. The old rule compounded: one column overshooting left every
  column after it short, and taking three panels off the gallery sheet turned a
  159 pixel difference between the columns into 317. The cuts are now chosen by
  weighing all of them, which for a sheet of eighteen panels is nothing, and the
  same sheet comes out 66 pixels shorter with a 111 pixel difference.

## [0.14.0] - 2026-07-31

### Added

- `Panels`, several figures in one document, which is what a multi-panel paper
  figure is and what an overview of a whole library needs. Each figure is
  embedded as a nested `<svg>` inside a translated group, so nothing is
  reparsed and nothing is rewritten: a panel on a sheet is the same picture as
  the panel on its own, and there is a test that asserts exactly that.
- `assets/gallery.svg` and `examples/gallery.rs`, one panel per kind of plot.
  **When a new track type is added it gets a panel here**, because an overview
  that quietly stops covering everything is worse than no overview: it looks
  complete.
- `MatrixTrack::tree`, the same phylogeny ordering `SnpTrack` already had. For
  a pangenome it is the difference between a figure and a shrug: sorted by
  sample name an accessory region is a speckle, sorted by descent it is a
  rectangle, and a rectangle is a claim about the biology.
- `assets/example-pangenome.svg`, a presence and absence matrix ordered by its
  tree.
- `track::tree::leaf_order`, the row permutation both tracks now share, so they
  cannot drift apart on how an unnamed row is handled.

### Fixed

- `MsaTrack::consensus` broke ties by the **last** residue seen rather than the
  first, contradicting its own documentation. It used `max_by_key`, which
  returns the last maximum. Two residues tied in a column therefore picked the
  wrong one, though the result was at least stable between runs.
- Six clippy lints that had been failing CI: three functions over the argument
  limit, two byte-string literals and a single-character `push_str`.

### Changed

- The CI workflow runs on `workflow_dispatch` only. Pushing no longer spends
  Actions minutes; run it from the Actions tab or with
  `gh workflow run ci.yml --repo PathoGenOmics-Lab/karyon --ref main`.

## [0.13.0] - 2026-07-31

Phylogenies, and the reason to have them here rather than in a tree viewer.

### Added

- `karyon::tree`, with `Tree::parse_newick` reading nested clades, branch
  lengths, quoted names, scientific notation and internal labels, taking an
  internal label as a support value when it parses as a number and as a name
  when it does not. Nodes live in one flat list and refer to each other by
  index, so a ten thousand deep ladder lays out without touching the stack.
- `Tree::layout`, placing every node at a depth and a row, as a phylogram from
  the branch lengths or as a cladogram counting branches instead.
- `TreeTrack`, a phylogeny as a track of its own, with `TreeShape` choosing
  between the two.
- `SnpTrack::tree`, which is the pairing the tree is here for. It draws the
  phylogeny in the strip beside the panel and **sorts the rows by descent**, so
  a clade's shared substitutions line up into a block instead of scattering
  down the panel in sample-name order. Rows match leaves by name, and a sample
  the tree does not mention keeps its place at the bottom rather than
  vanishing: a row silently dropped from a figure is worse than a row out of
  order.
- `Error::InvalidNewick`, so a malformed tree is rejected with a reason rather
  than guessed at.

### Notes

- Branches are drawn rectangular rather than as diagonals. A diagonal would
  imply the tree says something about the space between two rows, and it does
  not.

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

[0.14.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.14.0
[0.13.0]: https://github.com/PathoGenOmics-Lab/karyon/releases/tag/v0.13.0
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
