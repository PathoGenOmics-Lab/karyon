# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Metadata columns beside the rows of six tracks. A track drawn as a row per
  named thing answers which ones, and the question after it is almost always
  what they were, which is in a sample sheet rather than in the file the track
  was drawn from. `MatrixTrack`, `MsaTrack`, `SnpTrack`, `CladeTrack`,
  `DomainTrack` and `LocusTrack` take `Traits` now, and draw one narrow strip
  per attribute in the strip they already reserve for the row names and the
  dendrogram. A pileup has rows too and they are reads rather than samples, so
  it does not take one.
- These are not a track and are not going to be. A sample's lineage is not at a
  base, so a track of them would need an x nobody has and the first pan would
  slide a sample's lineage off the end of that sample's own row. They are drawn
  in the gutter and survive every pan and zoom untouched, which is the whole
  reason they are attached to a track rather than stacked beside one.
- `read::sheet`, which reads a sample sheet: a header, names in column one, and
  a column per thing known about them. A field is a number when it parses as
  one, a flag when it spells one, and text otherwise; an empty field, a `.`, an
  `NA` and anything that parses to a number that is not a number are all
  absent, which is the same three spellings `read::table` already reads as
  nothing.
- `--traits <FILE>` on the command line, with `--columns <A,B,C>` to pick the
  columns and their order. `--columns` may be written before the `--traits` it
  picks from, since two modifiers of one track are not in an order. A sheet
  naming none of the rows is refused, and so is a column the sheet has not got,
  which is named alongside the ones it has.
- The `snps` example carries the metadata it was always about. Lineage,
  resistance and collection year sit between the isolate names and the panel,
  read out of a sample sheet the example prints. Lineage runs in three blocks
  because the phylogeny put the rows in that order and a lineage is inherited;
  resistance does not, which is the whole reason two strips are worth reading
  side by side. The reference and the one isolate the tree does not carry are
  absent from the sheet, and their cells say so rather than guessing.
- Twenty-one examples in the playground rather than six, which between them
  use every one of the twenty-five flags the command has. Each was run through
  the binary before it was written down, so none of them is an example that
  does not draw.
- The playground's examples open in a panel rather than a menu, because they
  are pictures and pictures need room. The measurements are the ones the tool
  it is modelled on uses, since they are good ones: a panel of `min(68rem,
  100vw)` at full height over a dimmed and slightly blurred page, a pill to
  close, a search across the top, and a grid of cards that lift when the
  pointer is on them. Escape and a click on the page behind close it.
- Each card's preview is drawn by the program when the panel opens, out of the
  same files the example loads. The tool this follows ships a screenshot per
  example; a picture of a figure is a different claim from the figure, and at
  a millisecond each there is no reason to make it.
- The site is laid out the way a tool's site is laid out: a hero, a band of
  prose, a figure the reader can take hold of, the concepts, what it will not
  do, and where to go. Five destinations rather than eight, and no file moved,
  because a page's URL comes from its path and every link that exists points at
  one of those. Headings are set in their own face and section headings carry a
  rule, so a page of eight hundred lines can be scanned rather than read.
- The home page draws its own figure. The picture it ships with is what the
  printed command actually draws, and the program replaces it on arrival, so
  the page is correct with JavaScript off, with the wasm missing, and for the
  moment the fetch takes. It can then be dragged, zoomed with the wheel once it
  has focus, and moved with the arrow keys, and each of those re-runs the
  program over the three files printed underneath it, which the script reads
  out of the page so that what is shown and what is drawn cannot drift.
- `assets/karyon-wasm.js`, one bridge for the two pages that run the program.
  The home page and the playground had a copy each of the protocol, the command
  line splitting and the region arithmetic, and two copies of a protocol are two
  things that drift.
- The playground is a workbench rather than a page: a toolbar, an editor with
  one tab per file, and the figure filling its own pane, with a splitter
  between them and a stacked layout for a narrow screen.
- Figures in the playground can be taken hold of. **Interactive** makes a drag
  pan and the wheel zoom, and each of those re-runs the whole program: the
  region in the command box is what is being changed, so every frame is a
  figure `karyon` drew at that region from those files rather than a picture
  being scaled. A figure of twenty thousand rows over two megabases redraws in
  about sixteen milliseconds, which is what makes it possible at all.
- A playground on the documentation site: the command line running in the
  reader's own browser, over files typed into the page. It is the same code the
  terminal runs, because it is the same code. `cli::stack::build` takes a
  closure that answers with a source's text rather than opening a path, which
  is what makes it possible: a browser has no disk, so a shell hands the
  grammar one and the page hands it a textarea. Every flag, every reader and
  every refusal is what the terminal gives, including the counts, and the
  figure is a real `<svg>` element rather than an image of one.
- `playground/`, the bridge that runs over. It is a crate of its own and not a
  member of anything, because it is the one place in the repository that needs
  raw pointers and the library forbids them, and because doing it with
  `wasm-bindgen` would mean a dependency tree for a crate that has none. One
  buffer in and one buffer out, both length prefixed, and the protocol is
  written down in the module. The wasm is built by the documentation workflow
  rather than committed, so the page cannot serve a program older than the
  library it came from.
- `read::bisulfite` reads a `bismark_methylation_extractor` file, which is one
  row per cytosine per read. The track is a matrix, one shared list of sites and
  one call per site per molecule indexed by the site's place in that list, and
  nothing downstream checks that a caller got it right, so a row built by
  pushing calls in the order the file listed them puts every call after the
  first gap one column to the left: a methylation pattern that never existed,
  drawn as cleanly as one that did. Every row is built at the full width and
  written into by position. Both mates of a pair are one fragment and so one
  row, and where they overlap and disagree neither call is kept.
- `read::domain` reads an `InterProScan` table, and its axis is residues rather
  than bases, which makes it the only track here that is not drawn over a
  genome. Column one names the row rather than selecting it, as in the locus
  reader, since the figure is the comparison. The protein's length is column
  three and is never the furthest domain: a protein whose last annotated domain
  ends at 300 may run to 800, and a backbone drawn to 300 says the domain
  reaches the C terminus. `hmmscan --domtblout` is deliberately not read; which
  of its columns holds the protein depends on which program wrote it, and the
  only thing that says so is a comment line that goes the moment output is
  piped.
- `--bisulfite` and `--domains`, which is twenty-five of the thirty-three track
  types the command reaches, and eight left in the library.
- `--context` and `--analysis`, which join `--modification` as the flags that
  say which of the several things a file holds to draw. All three are one
  mechanism now, `Kind::selector`, spelled by what is being chosen and
  exhaustive so that a track added without one cannot silently stack whatever
  it found.
- `read::methyl` reads bedMethyl, as `modkit pileup` writes it. Coordinates
  pass straight through, since it counts from nought like BED. Two things a
  pileup writes are not measurements and are refused as such: a position with no
  valid coverage, which the file marks by writing nought in every count and
  which read as a fraction is a mark on the baseline saying the cytosine is
  unmodified; and column eleven, the rendered percentage, which is the one
  field a pileup can write the word `NaN` into. The fraction is the two counts
  instead, so a value that is not a number cannot be built. A file counting more
  than one modification is refused until `--modification` says which, because
  `m` and `h` at one cytosine stacked on one axis are two marks naming neither.
- `read::structural` reads structural calls out of a VCF. `POS` is the base
  before a symbolic event, so it becomes the 0-based start unchanged and `END`
  the exclusive end unchanged: both conversions are the identity, for two
  different reasons, and copying the point reader's `POS - 1` moves every call
  two bases left. `SVLEN` is taken absolute, since 4.3 writes a deletion's
  length negative and 4.4 positive. A call stating no length at all is refused
  rather than made one base wide, `<CNV>` is counted rather than drawn as a
  duplication, and a breakend finds its mate in its own ALT so that the
  reciprocal record does not draw the same arc again.
- `read::split` reads molecules that aligned in pieces, from SAM and its `SA`
  tag. The order of the pieces is the whole claim the figure makes, so it is
  computed from where each alignment sat on the read rather than assumed. A
  CIGAR runs along the reference, so a reverse-strand alignment's clips are the
  far end of the molecule and its place on the read is measured from the other
  side; without that step every read crossing an inversion comes out in the
  opposite order and draws the mirror of what happened. Only primary alignments
  are read, since a supplementary one is already an entry in its primary's tag
  and reading both counts every piece twice.
- `--methylation`, `--structural` and `--split-reads`, which is twenty-three of
  the thirty-three track types the command reaches. `--modification` names the
  modification a pileup counted.
- `read::clade` reads the GFF3 that recombination detection writes, `Gubbins`
  above all, whose ninth column carries the taxa an interval is carried by. It
  is the one file shape in ordinary use holding a span and its carriers in one
  record. The list is split before it is decoded, because the ninth column
  spends the comma on its own syntax and decoding first turns one name into two;
  a repeated taxon is dropped and counted, because a block claiming more
  carriers than the tree has rows leaves what it covers to a subtraction with no
  answer; and a file naming exactly one sequence is read whatever that sequence
  is called, because `Gubbins` writes the literal `SEQUENCE` and filtering on
  the region's own name would draw a phylogeny with nothing on it, which reads
  as a statement that there was no recombination here.
- `read::locus` reads gene neighbourhoods and the homologies between them.
  Column one names the genome rather than selecting it, which is the one place
  an interval file is read differently here, so the file is the concatenation a
  shell already produces. `links` takes BLAST tabular output, which DIAMOND and
  others write too, or two or three columns of names, and does the join a
  `Homology` needs: it refers to genes by position and every file that produces
  one names them instead. The join is exact and every way it fails is counted
  rather than resolved, since a name pointed at gene nought draws a well formed
  ribbon into the leftmost gene and takes that gene's unmatched outline off.
- `--clades` and `--loci`, taking their second file by name as `--tanglegram`
  does: `--with-tree` for the phylogeny a clade track is painted onto and
  `--links` for what joins one locus to the next. Both are required. Each of
  these tracks draws a finished-looking figure without its second file and each
  of those figures says something strong and false, a locus track most of all,
  which marks every gene no homology reaches and so turns a forgotten flag into
  the loudest positive finding it can make. That is twenty of the thirty-three
  track types the command reaches.
- `--identity` says whether a homology file's third column is a percentage, as
  BLAST and DIAMOND write it, or a fraction, as others do. Left out it is worked
  out from the values and refused where they cannot say, because read the wrong
  way round every ribbon becomes a perfect match and nothing fails.
- The command line moved into the library as `karyon::cli`, and stopped
  opening files. `cli::stack::build` now takes a closure that answers with a
  source's text, so the same grammar that a shell drives from disk can be
  driven from a browser holding the files in memory. `cli::args::Kind::ALL`
  lists the tracks it can draw, for a caller offering a choice of them.
- Five more track flags. `--synteny` and `--dotplot` draw a PAF from
  `minimap2` as ribbons or as a dot plot; `--orfs` and `--logo` compute off the
  same FASTA and aligned FASTA that `--sequence` and `--msa` already take; and
  `--tanglegram` puts two phylogenies face to face. That is eighteen of the
  thirty-three track types the command reaches, up from thirteen.
- `read::align_pairs` reads PAF. Both of its coordinate pairs are already
  0-based and half-open, so it is the one reader in the directory that moves
  nothing, and a whole-genome file says how many rows belonged to a pair the
  figure is not about rather than stacking them on one axis in silence.
- A track flag can now name a second file, for a track whose data is not one
  file: `--tanglegram left.nwk --against right.nwk`. The second is spelled by
  what it means rather than by where it sits, as every other modifier is, and a
  track that takes one is refused without it, because a tanglegram of one tree
  against itself has no crossings and no crossings is what a correct answer
  looks like.
- An integrated evolutionary-genomics and surveillance vocabulary. Rectangular
  `TreeTrack` branches can now be orthogonal, diagonal or curved without
  changing topology; `AncestralStateLayer`, `BranchEventLayer` and
  `BranchIntervalLayer` keep state posteriors, direct events and uncertainty
  on their correct nodes or edges in rectangular, circular and unrooted
  projections. `PhylodynamicTrack` adds linear or log trajectories with
  uncertainty ribbons, while `SurveillanceTrack` adds observed lineage counts,
  denominators, stacked or line views, sampling floors and explicit frequency
  or growth alerts. The `evolutionary_surveillance` example joins these with
  molecular selection, genomic layers and a core-versus-accessory tanglegram
  on one eight-panel sheet.
- Molecular-selection figures can now preserve the structure of branch-site
  and site-wise results rather than compressing every analysis into one branch
  colour. `BranchRateMixture` draws weighted ω-class capsules on rectangular,
  circular and unrooted trees; `HomoplasyLayer` connects equal direct branch
  events as recurrence candidates; and `SelectionTrack` aligns p-value or
  posterior evidence above a signed, neutral-centred site effect. Missing
  estimates remain missing, source weights and rates stay exact in SVG
  tooltips, and evidence is never allowed to impersonate effect size. The
  `selection_atlas` example combines all three with protein architecture and
  iTOL-style terminal metadata.

### Changed

- `TraitColumn`, `TraitScale` and `TraitStyle` moved from `track::tree` to
  `track::traits`, unchanged, and are still exported from `karyon`,
  `karyon::track` and `karyon::track::tree`. They stopped being a phylogeny's
  business when six other tracks started drawing them, and the tree now draws
  its own through the same function the others do, so a lineage cannot be one
  colour beside a tree and another beside the matrix under it. Every committed
  figure is byte for byte what it was.
- Two counts in the prose were wrong before this and are fixed with it: the
  command line reaches twenty-five track types, not eighteen, and `--help`
  prints twenty-four track flags and the ruler, not twelve.
- Release builds now check arithmetic for overflow. A wrapped add is a figure
  that is quietly wrong, which is the one thing this crate is not allowed to
  be: `Genome` added lengths without saturating, and a release build turned
  that into boundaries running backwards and every point drawn on the wrong
  sequence, with nothing said. A branch per arithmetic operation is nothing in
  a crate whose work is string formatting, and it converts the whole class from
  silent to loud. A hundred and fifty thousand generated figures and every
  example render clean under it.
- Nine functions turned an `f64` into text; there are now three, and the two
  that remain outside `svg` are four-line wrappers over them. Two of the copies
  were literally identical and a third was the same plus a check for `NaN` that
  the others never received, which is what a copy costs: a fix reaches the copy
  the bug was found in and stops there. The same mistake had by then appeared
  three times in three modules. `svg::text_rounded` and `svg::text_exact` are
  the two ways a number reaches a reader, and both carry the guards; `svg::num`
  is unchanged and remains the way a number reaches an attribute.
  `svg::text_number`, added earlier in this same unreleased version, is gone in
  favour of `text_rounded`.
- `svg::finite_within` replaces the two private helpers, one in the tree track
  and one in the map, that differed only in whether the upper end was named.
  Almost every setting a drawing takes is a length, an angle or a fraction and
  none of them has a reading for `NaN`, so there is now one answer to what
  happens then. Four pixel thresholds that had no answer at all, in
  `LogoTrack`, `PileupTrack` and `SequenceTrack`, use it.
- The three largest files became directory modules, split along what they
  actually do rather than at a line count. `src/track/tree.rs` was 4,884 lines,
  two and a half times the next largest, and is now seven files: the three
  projections apart (`rectangular`, `radial`, `unrooted`), the layout they all
  read (`scene`), what is drawn on top of any of them (`decorate`), the two
  ways a branch length is given a unit (`scale`), and the type and its
  vocabulary in `mod.rs`. `src/map.rs` became the plain map, `phylo`, the
  marks both compositions share (`draw`) and the projection maths (`project`),
  which is worth its own file because every function in it can answer that it
  will not place a point. `src/tree.rs` gave up its Newick, annotated Newick
  and NEXUS parsing to `tree::parse`, which is where all the format's ambiguity
  lives. Nothing else moved: the remaining files are one subject each, and
  splitting a coherent module into thirds makes it harder to read, not easier.
  Every public path is unchanged, and all forty-one committed figures render
  byte for byte the same, which is what makes this a move rather than an edit.
- The format readers moved from the command line binary into the library, as
  `karyon::read`. Four thousand lines that parse BED, bedGraph, GFF3, VCF, SAM,
  cytoBand, `samtools depth`, FASTA and Newick were reachable only by running
  the command, so writing Rust against this crate meant writing a VCF parser
  first, which is the opposite of what a crate that ships twelve file-backed
  tracks should ask for. The rule that made them live there is unchanged and
  now says something sharper: every reader takes a `&str`, so nothing in the
  library opens a path, and the binary keeps the one part that is genuinely
  the command line's. `Format` moved with them, since it names file formats
  rather than flags.

### Fixed

- Closing a file tab wrote the closed file's text into whichever file took its
  place. The editor remembered which file it held by position, and removing one
  shifts every file after it down, so the save that runs on the way out landed
  on the wrong object: closing `depth.bg` left `genes.gff3` holding a bedGraph.
  It remembers the file itself now, so a file that has gone is not written to
  and one that has moved still gets its own text.
- The splitter announced itself as a separator with a value and answered to no
  key at all. Arrows move it, Shift moves it further, Home and End take it to
  its ends, and it reports where it is. Its orientation follows the layout
  rather than being frozen at whichever one the page opened in.
- Interactive turned on a pan and zoom no keyboard could reach: only the home
  page's figure had keys. The playground's figure joins the tab order while it
  is interactive and leaves it when it is not, and answers to the same keys.
  Its wheel now waits for focus too, so a wheel over it scrolls the page.
- A refusal in the playground was announced to nobody. The status line is a
  live region, as the home page's already was.
- A file could not be removed without a pointer, and the close control being a
  span inside the tab made every tab's name read as "depth.bg times". It is a
  button with its own name, the tab carries its own, and F2 renames a file and
  Delete removes one.
- The grabbing cursor vanished a pixel into a drag, because a redraw rewrote
  the whole class attribute and took `pg-dragging` with it.
- A region the pan and zoom arithmetic was happy to reach stopped the program
  dead, and the page said nothing. `MAX_SPAN` was `1 << 32`, which the comment
  beside it called "below the point where the buffer cannot be allocated";
  thirty-two gigabytes is above that point on every target there is, and on the
  32-bit one the documentation site runs on, a `Vec` cannot exceed `isize::MAX`
  bytes, so a per-base track of `1 << 28` bases is a capacity overflow, which
  is a panic, which in a build that aborts on one is a trap with no `Result` to
  carry it. Measured: a span of 268,435,455 draws and 268,435,456 stops the
  program. The limit is that number now, which leaves every chromosome anyone
  sequences inside it, and says why.
- Nothing caught that trap, and both pages write the new command down before
  they draw, so twenty-two turns of the wheel left the command box saying
  `chr1:1-271,050,503` over the previous figure, the previous region readout
  and the previous timing, with nothing marked as failed. The page asserted a
  region it had never drawn, which is this crate's own bug class in a browser.
  A trap is an ordinary refusal now, and the wheel cannot reach a span the
  program will not draw.
- Without the program, the example panel came up as twenty-one empty boxes,
  which reads as twenty-one figures that failed. The boxes go and one sentence
  says why; the cards keep their command and their files, which are still worth
  reading and still run in a terminal.
- Dragging a figure towards the first base of a sequence shrank the window
  instead of moving it. The start was clamped to 1 and the end left where it
  was, so from `chr1:1-1,000` five drags gave 700, 490, 343, 240 and 168 bases
  and the figure zoomed itself in while the reader was only moving it sideways.
  The span asked for is kept and the window slid. The same clamp stopped a zoom
  out from ever reaching its own limit.
- A `--label` that looked like a locus, written before the real one, was taken
  for the region, so a drag rewrote the label and left the figure where it was.
  The four flags that stand alone are known now, and the word after every other
  flag is skipped, as the parser itself does.
- The example panel said `aria-modal` and did not mean it: with it open the page
  behind was still reachable by Tab. It is a real `dialog` opened with
  `showModal`, which contains the focus and makes the background inert.
- The Examples button went on saying it was expanded over a panel that had
  closed, because the bookkeeping hung off the dialog's `close` event, which a
  listener added directly in front of a `close()` call was measured never
  receiving. It is done where the closing happens instead.
- On a 375px screen the toolbar wrapped but the groups inside it did not, so
  `Full screen` was drawn from 326 to 426 with the window ending at 375 and
  nothing scrolling: the button could not be reached at all.
- The example previews were drawn on the next animation frame, which a tab
  that is not on screen never gets, so a panel opened in a background tab came
  up with every card empty. They need no layout, so they are drawn as the card
  is made.
- The example panel came up cream in the dark scheme, with light text on it.
  Material writes the palette onto `body` and the rule was written against
  `:root`, which is the same mistake this stylesheet had made once before.
- The playground's first example refused on the fourth drag. Its three rows
  covered the region exactly, so moving the window took it off the data, which
  is the program being right and the example being too small to show anything
  else. It is six bedGraph rows and eight calls over ten kilobases now, and the
  home page says out loud that a window taken far enough will always find
  nothing, and that what comes back then is a sentence rather than an empty
  figure.
- A wheel over the figure no longer eats the page's scroll. It zooms only once
  the figure has focus, which is also what makes the arrow keys reachable
  without a pointer.
- A GFF3 span whose end came before its start was read as a one-base feature at
  the start rather than refused. The module doc had said all along that an
  inverted span stops the read, and only BED and cytoBand did it, so a gene
  written `400 100` was drawn confidently three hundred bases from where either
  number put it.
- A clade block naming one taxon twice counted two carriers of one row, so the
  rows it cuts out came out of a subtraction with no answer: a panic where this
  crate checks its arithmetic, and a number near the top of a `u64` in a release
  build that does not. Carriers are distinct rows now, which is what the figure
  draws either way.
- A track whose second file was missing was refused by a match with a fallback
  arm, so a track added without a spelling would have drawn without a file it
  cannot do without. It is exhaustive now, as the flag spellings next to it are.
- A value that is not a number was drawn as the strongest one on the page.
  `mix` blends two colours by an amount, `clamp` propagates a NaN, and `NaN as
  u8` saturates to nought, so an amount that was not a number left the ramp
  entirely and came back `#000000`: darker than the dark end, on scales where
  darker means more. Every ramp in the crate goes through that one function.
  A homology of unknown identity was the visible case, drawn blacker than a
  perfect match.
- `Homology::identity` is `Option<f64>` rather than `f64`, matching
  `AlignmentBlock`. Two genes matching with nobody stating how closely is an
  ordinary case, and `Homology::new` now records a NaN as the absence it is
  instead of clamping it into the range. `Homology::unstated` says it outright.
  An unstated identity draws at the pale end of the ramp with a dashed edge, so
  the figure distinguishes "nobody said" from "barely alike" without anyone
  having to point at it.
- `AlignmentBlock::identity` had the same hole and is closed the same way.
- A tree annotation whose number is not a number printed as `NaN` in every
  tooltip, trait strip and node label, through the one function the crate had
  already centralised annotation formatting in. It is spelled with the same
  glyph an absent annotation already used, since a number field that is not a
  number and a field that is not there are the same fact about the same field.
  The infinities are untouched: a file can carry a number too large for an
  `f64`, it parses to one, and reporting it verbatim says what the file said.
- `tests/properties.rs` checks tooltips as well as attributes. The writer's
  guard covers numbers written into attributes, and a tooltip is a caller's
  number put through `format!` into text a person reads, which is the gap all
  of the above came through. The property found the tree annotation case on
  its own, in a track nobody was looking at.
- A track flag written without its path named the wrong flag. The spelling was
  kept in step with the list of tracks by a fallback arm rather than by the
  compiler, so `--synteny` with nothing after it reported that `--axis` needed
  a value. It is an exhaustive match now, and the help text is checked against
  the parser's own list of tracks instead of a copy of it.
- A tanglegram summarising one crossing said `1 crossings`.
- `SquiggleTrack` never finished drawing over a wide region. The envelope
  branch walked the axis positions a pixel column covers rather than the
  samples in it, so over a chromosome each column asked thousands of bases
  whether they happened to be a sample, and over the whole coordinate range it
  asked billions and the figure never came back. The walk is now cut down to
  where samples can be before it starts.
- `Genome` added lengths and gaps without saturating, in all seven places it
  accumulates. A debug build panicked; a release build wrapped, which is worse,
  because `boundaries` then came back decreasing and every point in a
  genome-wide plot landed on the wrong sequence with nothing said.
- A figure whose tracks each reported a height that was a number could still
  report a total that was not. The height check reads one track at a time, so
  two finite heights adding to infinity passed it, and the document then said
  `height="0"` while `Figure::dimensions` said infinity: the file and the API
  disagreed about the same figure. The total now has a ceiling at the largest
  integer an f64 holds exactly.
- Tree scale bars, support labels and clade support tooltips printed a large
  value as three hundred digits. They were using `num`, which writes a
  coordinate and stays plain on purpose. `svg::text_number` is the sibling of
  it for numbers a person reads, agreeing with `num` to the byte below two to
  the fifty-third and switching to an exponent above it. This was the third
  place the same mistake appeared, after `VariantTrack` and `Map`, which is why
  it now has a function of its own rather than a third patch.
- A map tooltip printed `f64::MAX` as three hundred and nine digits. The same
  bug as the one fixed in `VariantTrack`, in the other writer, found because
  the properties had never been pointed at a map: `data_number` handed the
  value straight to a formatter, which expands it in full. An f64 stops
  holding consecutive integers at two to the fifty-third, so past that all but
  the first sixteen digits are the formatter filling in a gap rather than
  anything measured, and the same happens under a thousandth as a run of
  zeros. Outside that range the value is now written as an exponent, which
  says the same number and says how much of it is known.
- Six committed figures disagreed with the code that draws them: the MSA row
  labels had moved to full ink, four phylogenetic figures had shifted
  geometry, and the gallery had grown. A seventh, the node faces figure,
  reached the documentation without its `assets/` original ever being added.
  The staleness gate could not have caught that last one, because `git diff`
  compares against the index and a file that was never added is not stale to
  it, it is invisible; the gate now stages the renders first.

### Added

- A branch-wise dN/dS visual grammar for `TreeTrack`. `dnds` uses a
  colour-vision-safe logarithmic scale centred on ω = 1 instead of stretching
  the observed minimum and maximum across a generic ramp. Values are direct
  branch annotations and never inherited; missing or invalid fits are dotted.
  `dnds_significance` can use a p- or q-value to carry evidence through stroke
  weight while colour continues to carry effect size. Rectangular, circular,
  unrooted and cladogram examples share one explicit legend and exact SVG
  tooltips.
- The generated figures now reach all thirty-three track types, not sixteen, and
  build the tree ones through their whole builder rather than their
  constructor. The fourteen original gaps had never been given a hostile input and eleven
  numeric setters had never been given a hostile number; four bugs came out of
  it, each one now pinned by a test that was watched failing first.
- A property over every reader at once: the same interval written as BED,
  GFF3, cytoBand, SAM, bedGraph, `samtools depth`, VCF, an association table
  and a genotype matrix has to come back as the same two 0-based numbers. The
  audit beside the readers already pins every format to one known base; this
  asks the same question without picking the base, because a reader can be
  right at position 100 and wrong at 0, and because two readers disagreeing is
  what an off-by-one becomes once a figure stacks them. Every one of the six
  places the crate converts from 1-based was shifted on purpose to confirm the
  property notices. It also checks that naming a format and letting the reader
  work it out give the same answer, which is where these readers went wrong in
  practice rather than in the arithmetic.
- Twelve properties over phylogenies and maps, the two subsystems the suite
  had never covered, and each one checked by breaking the library on purpose
  to confirm it notices. Trees: every node placed exactly once with the tips
  filling the rows, a parent between its children, the rows identical whether
  depth is measured or counted, a cladogram depth equal to the number of
  branches, and a parser that decides rather than panics on a file that is not
  a tree. Rerooting is held to the only thing it must not do, by an oracle
  that walks the undirected edges itself: every tip-to-tip distance survives
  rooting on an internal node, on an outgroup and at the midpoint, and a
  refused rerooting leaves the tree untouched. Maps: a valid document, the
  same document twice, no location dropped without the notice saying how many,
  coordinates kept rather than clamped, and nothing projected off the page.
  The last of those covers path data as well as attributes, since that is
  where a projection actually lands.
- Row-major `Panels` grids for balanced comparison sheets, plus quieter clade
  fields, halo-separated node compositions, compact annotation legend chips,
  aligned-tree dividers and washed domain architectures.
- Node-attached bubbles, pies, donuts and stacked bars from exact numeric
  annotations; clade highlights that become rectangular bands, radial sectors
  or unrooted fields; tree-aligned MSA rows; and a new `DomainTrack` for domain,
  motif, exon or repeat architectures ordered by phylogenetic descent.
- Crossing-aware tanglegram untangling, equal-angle unrooted trees and layered
  iTOL-style annotation strips, bars, binary markers and symbols.
- Projection-independent branch evidence for phylograms: support can use scaled
  symbols, exact labels and a threshold in either common convention;
  branch-specific annotations can follow their edge without ancestral
  inheritance; and automatic or exact scale bars retain evolutionary units in
  rectangular, circular and unrooted coordinates.
- Validated rerooting around an internal node, a strictly monophyletic outgroup
  or the weighted midpoint. Tip distances and edge support are preserved,
  invalid requests do not mutate the tree, and rooted projections can mark the
  selected root explicitly without leaking that marker into unrooted views.
- A visual plot catalogue for the documentation, separating discovery from the
  exhaustive track reference. Eight category pages group all thirty-three
  tracks and the three standalone drawings by biological question, data shape
  and coordinate system, with responsive image-led cards, shorter section-only
  tables of contents and direct routes to the exact API contract.
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
- Typed phylogenetic annotations through `parse_annotated_newick` and
  `parse_nexus`, including BEAST, NHX, rootedness markers and Nexus translation
  tables. `Tree` now exposes MRCA, ancestor and descendant queries, rotation,
  ladderising, internal-node rerooting, subtree extraction and collapse.
- Dated `TreeTrack` layouts, metadata-driven branch colours, visible internal
  nodes, non-destructive clade triangles and categorical or continuous trait
  columns aligned to the terminal taxa. The synthetic surveillance example is
  a committed visual regression target for the complete composition.
- Circular `TreeTrack` projections with complete or partial sweeps, controllable
  start angle and central gap, outward or inward branch direction, concentric
  time guides, annular categorical or continuous trait rings and
  non-destructive collapsed-clade wedges. The four-view radial gallery is a
  committed visual regression target.
- `Map`, `GeoLocation`, `GeoFlow` and three explicit geographic projections.
  Locations carry exact coordinates, counts and categories; categories use
  shape as well as colour, weighted links retain direction, and invalid,
  ambiguous or hidden inputs are counted rather than clamped or dropped. The
  Natural Earth 1:110m land geometry is public-domain data compiled into the
  crate, so rendering stays offline and deterministic.
- `PhyloMap` and `PhyloConnector`, composing an inward circular time tree,
  phylogram or cladogram around a central map. Terminal annotations match a
  unique supplied coordinate table, with aggregated, individual or absent
  connectors; the renderer does not infer coordinates, ancestral states or
  migration. A two-view circular phylogeography gallery pins the composition
  as a visual regression target.
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
- The README track table now lists all thirty-three track types. It had stopped
  at `AxisTrack`.

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

- A variant tooltip printed `f64::MAX` as three hundred and ten digits. The
  guard against that existed and its fallback expanded the number anyway, and
  its threshold was above the point where the precision it was protecting had
  already gone: an `f64` holds consecutive integers to 2^53 and therefore
  consecutive hundredths only to a hundredth of that, so 999,999,999,999,999
  came back reading `.04`. Past that the magnitude is all there is to say, and
  a new invariant says no tooltip may carry more than twenty digits in a row,
  which is what `u64::MAX` needs and more than any grouped number reaches.
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
