# Contributing

Issues and pull requests are welcome at
[PathoGenOmics-Lab/karyon](https://github.com/PathoGenOmics-Lab/karyon).

## Reporting a problem

A plotting library fails differently from a program that prints numbers: the
figure still renders. A variant one base to the left, a gene numbered from the
wrong end, a track clipped to somebody else's band, all of them come out as a
perfectly valid SVG that is wrong. So a report needs three things:

- **The code or the command**, complete, with every flag or every `add_` call.
  A snippet that someone else can paste into `examples/` and run is the fastest
  route to a fix.
- **The version.** `karyon --version` for the command, the `Cargo.lock` entry
  or the git commit for the library.
- **The SVG itself.** It is text, so it goes in an issue as an attachment
  without being converted to anything. If it is large, the few elements that
  are wrong plus the figure's dimensions are usually enough.

Say what you expected to see and what you saw. "The lollipop for S450L sits one
base to the left of the codon the ruler numbers 450" is a report; "the codon
track is broken" is not.

If the data cannot be shared, synthetic input of the same shape almost always
reproduces the problem. The positions, the lengths, the strand and the CIGAR
are what the drawing depends on; the bases themselves rarely are. Every example
in `examples/` generates its own data from a fixed seed for exactly this
reason, so there is a pattern to copy.

## What a change has to pass

The same four gates CI runs:

```bash
cargo test
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
```

`cargo test` covers more than it looks like. The suite is in four places:

- **Unit tests** next to the code, which is most of them: the arithmetic of
  every track, the scale, the parsers of `Region` and `Tree`, the panel layout.
- **The command line's own tests** in `src/bin/karyon/`, over the grammar in
  `args.rs` and the readers in `read/`. Every reader has a test that pins a
  known position through the 0-based or 1-based conversion, because that is the
  failure nobody sees.
- **Integration tests** in `tests/render.rs`, which check the thing a user
  actually gets: a well-formed document, no non-finite number anywhere in it,
  byte-identical output between two runs, unique clip ids, a variant landing on
  the centre of its base, and a 4 Mb genome-wide figure staying under 100 KB.
- **Doc tests**, since every ```` ```rust ```` block in the crate documentation
  is compiled and run. A doc example that stops compiling is a failing test.

`cargo test --release` is worth running too, and CI does: a few of the checks
are about floating point and layout arithmetic that the optimiser is allowed to
rearrange.

## The figures are part of the build

Everything under `assets/` is rendered by an example, and rendering is
deterministic, so a stale figure is a diff. After any change that could touch
the drawing, re-render and look at what moved:

```bash
cargo run --example locus -- assets
cargo run --example logo -- assets
cargo run --example pileup -- assets
cargo run --example ideogram -- assets
cargo run --example association -- assets
cargo run --example synteny -- assets
cargo run --example msa -- assets
cargo run --example snps -- assets
cargo run --example pangenome -- assets
cargo run --example gallery -- assets
```

CI runs exactly those ten and then `git diff --exit-code -- assets`, so a
rendering change that was not committed fails the build. The other examples
write into the same directory and are worth running for the same reason:

```bash
cargo run --example genomewide -- assets
cargo run --example circular -- assets
cargo run --example selection -- assets
cargo run --example shapes -- assets
cargo run --example reading -- assets
cargo run --example niche -- assets
```

That is every example there is, and between them they produce every file in
`assets/`, byte for byte.

The documentation site reads its figures from `docs/assets/figures/`, which is
a copy. Refresh it after re-rendering, or the site keeps showing the old
drawing:

```bash
cp assets/*.svg docs/assets/figures/
```

A diff in `assets/` is not a problem in itself. It is the review: open the old
and the new figure side by side and check the change is the one you meant.

## The documentation

The site is MkDocs Material and the pages are in `docs/`:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

```bash
mkdocs build --strict
```

`--strict` turns a broken internal link, a dead anchor or an image whose file
is not there into a failed build, so a page that points at a figure nobody
committed does not get published. Every image needs alt text that says what the
figure shows, not what number it is.

Two house rules. **English throughout**, in the prose, the code and the
comments. And **no em-dash characters anywhere**, U+2014, which is easy to
check without typing one:

```bash
grep -rn "$(printf '\xe2\x80\x94')" docs src examples README.md CHANGELOG.md
```

That has to find nothing. A comma, a colon or a full stop says the same thing
and survives every editor and every font.

## Adding a track type

The entry test first, because it is the reason this crate exists rather than a
general plotting library: **does the track live on the genomic coordinate
axis?** If its `draw` never reads `ctx.scale`, its x is a sample list or a
category and the plot is a bar chart, a line chart or a heatmap that happened to
be handed genomic data, which matplotlib already draws better. Three track
types were removed under this rule rather than kept for the sake of a longer
list.

A track that passes it is a small amount of code. The trait is in
`src/track/mod.rs` and asks for two things, a height for a given scale and a
`draw`, with `label` and `y_axis_width` on top of them; the figure decides where
the band goes and clips the output to it. A complete track type is on the
`Track` documentation, and it is about a dozen lines.

What a new track shipped by the crate then needs:

- **Its own file** in `src/track/`, exported from `src/track/mod.rs` and
  re-exported from `src/lib.rs`.
- **An entry in the `tracks!` list in `src/plot.rs`**, which is what makes
  `Plot::label` and `Plot::adjust` work on it, and an `add_` method beside the
  others. Tracks that lay an array along the axis also get the `_at` form.
- **A panel in `examples/gallery.rs`.** The gallery is meant to cover
  everything, and an overview that quietly stops being complete is worse than
  no overview, because it still looks complete.
- **Tests**, including one that pins a known coordinate to a known pixel. Track
  arithmetic is where the silent errors live.
- **A row in the README table** and a line in `CHANGELOG.md` saying why the
  track exists, not just that it was added.

A reader for it belongs in the command line front end only if the format is
line-based text. The library reads no files and has no dependencies, and both
of those are on purpose.

## Conventions the code follows

- **Coordinates are 0-based and half-open**, everywhere, with two exceptions
  a reader sees: `Region::parse` takes 1-based inclusive locus strings, and tick
  labels are printed the same way. Every constructor says which it takes.
- **Scale awareness is not optional.** A track handed four million points must
  not emit four million elements. Bin to one value per pixel column, or draw a
  hint, the way `CoverageTrack` and `SequenceTrack` do.
- **Deterministic output.** The same input renders byte-identical output, which
  is what makes `git diff -- assets` a test at all. Category colours follow
  first appearance rather than hash order, because a figure that recolours
  itself when a sample is added cannot go in a paper.
- **Plain SVG 1.1.** No scripts, no external references, no embedded fonts, so
  the file opens unchanged in a browser, in Inkscape and in Illustrator.
- **Errors are values.** A bad locus string is an `Error`, not a panic, and it
  converts into `io::Error` so a region and the file it renders to can share
  one `?`.

## CI

The workflow is `.github/workflows/ci.yml` and it is **`workflow_dispatch`
only**. Pushing does not start it. Run it when you want it, from the Actions
tab or with:

```bash
gh workflow run ci.yml --repo PathoGenOmics-Lab/karyon --ref main
```

It has two jobs. The first runs format, clippy, `cargo test`,
`cargo test --release` and `cargo doc` with `RUSTDOCFLAGS: -D warnings`, on
Ubuntu and macOS. The second renders the ten README figures and fails if
`git diff -- assets` is not empty.

Running it on purpose rather than on every push is a deliberate choice about
Actions minutes, not an oversight. It also means the gates above are yours to
run locally before you open a pull request.
