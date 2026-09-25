# Contributing

How to report a figure that came out wrong, what a change has to pass before
it is merged, and the conventions the code keeps. Issues and pull requests are
welcome at [PathoGenOmics-Lab/karyon](https://github.com/PathoGenOmics-Lab/karyon).
{ .k-lead }

## Reporting a problem

| You have | Where it goes |
|:--|:--|
| a figure that is wrong, a crash, or a refusal that should not happen | an [issue](https://github.com/PathoGenOmics-Lab/karyon/issues) |
| a question about how to draw something, or why it came out as it did | [Q&A](https://github.com/PathoGenOmics-Lab/karyon/discussions/categories/q-a) |
| a hang, a panic or broken escaping on input someone else supplied | a private report, as [SECURITY.md](https://github.com/PathoGenOmics-Lab/karyon/blob/main/.github/SECURITY.md) describes |

A plotting library fails differently from a program that prints numbers: the
figure still renders. A variant one base to the left, a gene numbered from the
wrong end and a track clipped into its neighbour's band all come out as a
perfectly valid SVG that is wrong. So a report needs three things:

- **The code or the command, complete**, with every `add_` call or every flag.
  A `fn main` someone else can drop into `examples/` and run is the fastest
  route to a fix.
- **The version.** `karyon --version` for the command, or the `Cargo.lock`
  entry or git commit for the library.
- **The SVG itself.** It is text, so it attaches to an issue as it is. If it is
  large, the few elements that are wrong and the figure's dimensions are
  usually enough.

Say what you expected to see and what you saw. "The lollipop for S450L sits
one base to the left of the codon the ruler numbers 450" is a report; "the
codon track is broken" is not.

!!! tip "Before reporting an empty track"
    A row naming another sequence, or lying outside the region on display, is
    skipped without a word. That is how a whole-genome file can be handed over
    to draw one locus, and it is the commonest reason for a track that comes
    out empty. See
    [Coordinates](../how-it-works/coordinates.md#what-a-files-numbers-become).

If the data cannot be shared, synthetic input of the same shape almost always
reproduces the problem: the positions, the lengths, the strand and the CIGAR
are what the drawing depends on, and the bases rarely are. Every example in
`examples/` generates its data from a fixed seed, so there is a pattern to
copy.

## What a change has to pass

CI runs on every pull request, and again on `main` after a merge. The core of
it is five commands, on Ubuntu and on macOS:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo test --release
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
```

The rest of the workflow, `.github/workflows/ci.yml`, checks what those five do
not reach:

| Check | CI job | To run it yourself |
|:--|:--|:--|
| The oldest supported compiler, 1.74, still builds everything | Oldest supported compiler | `cargo +1.74 check --all-targets --locked`, with that toolchain installed through rustup |
| The playground's WebAssembly bridge, a crate of its own in `playground/` | Check & Test | the same format, lint and test commands with `--manifest-path playground/Cargo.toml`, and a `--target wasm32-unknown-unknown` release build |
| The tree viewer's scripts | Check & Test | `node tests/tree-canvas.test.js` and `node tests/tree-radio.test.js` |
| Every figure is current | Example renders | see [the figures are part of the build](#the-figures-are-part-of-the-build) |
| Every tracked file has a label | Every file has a label | `python3 .github/scripts/labeler-coverage.py` |

A change to the documentation, or to anything the site is built from, also has
to pass `mkdocs build --strict` in the Docs workflow.

`cargo test` covers more than it looks like. The suite is in four places:

- **Unit tests beside the code**, which is most of them: the arithmetic of every
  track, the scale, the parsers of `Region` and `Tree`, the command line's
  grammar in `src/cli/`, and the readers in `src/read/`. Every reader pins a
  known base through its conversion, and `src/read/audit.rs` brings every
  format to the same base, because an off-by-one is the failure nobody sees.
- **`tests/properties.rs`**, which asserts what must hold for every figure and
  then generates figures from a seeded generator until one breaks it: every
  document is valid, rendering is deterministic, every format puts the same
  interval on the same bases, and hostile names never break the document.
- **`tests/render.rs`**, which checks what a user actually gets: a well-formed
  document, no number that is not finite, byte-identical output from two runs,
  unique clip ids, a variant on the middle of its base, and a four-megabase
  figure under 100 KB.
- **Doc tests.** The examples in the crate's documentation are compiled and
  run, so an example that stops compiling is a failing test.

`cargo test --release` runs the same suite optimised, the way the command is
built for use. The release profile keeps overflow checks on
(`overflow-checks = true` in `Cargo.toml`), so arithmetic that would wrap
fails loudly there too, rather than drawing a figure that is quietly wrong.

The pull request template asks for the rest: what you ran and what it printed,
a test for new behaviour, and for a fix a test that you watched fail before
the change.

## The figures are part of the build

Everything under `assets/` is drawn by an example, and rendering is
deterministic, so a figure that was not drawn again after a change is a diff.
After any change that could touch the drawing, render every example and copy
the result into the site's own copies:

```bash
for example in examples/*.rs; do
  cargo run --example "$(basename "$example" .rs)" -- assets
done
cp assets/*.svg docs/assets/figures/
```

CI does the same and fails when a committed figure, under `assets/` or under
`docs/assets/figures/`, disagrees with the code that draws it. Its list of
examples in `ci.yml` has to match `examples/`, so a new example is added there
too, or the job fails and names it. The run also turns every figure into a
PNG and uploads them as the `visual-gallery` artifact, which shows a rendering
change without checking out the branch.

An example builds its figures in a file of its own under `examples/figures/`,
one function per SVG named after the file, so `example-genomewide.svg` is
`example_genomewide`, taking a theme, a width and a region. The example writes
them, and the playground compiles the same file to draw them live on the site,
so a new figure is listed in `playground/src/committed.rs` as well. A test
there compares that list with `assets/` and names any figure one has and the
other does not.

A diff in `assets/` is not a problem in itself. It is the review: open the old
and the new figure side by side, and check the change is the one you meant.

## The documentation

The site is MkDocs Material. The pages are in `docs/`, the site's stylesheets
live under `docs/stylesheets/`, and its figures are the copies in
`docs/assets/figures/`.

```bash
pip install -r requirements-docs.txt
mkdocs serve
mkdocs build --strict
```

`--strict` turns a broken internal link, a dead anchor or a missing image into
a failed build, so a page that points at a figure nobody committed is never
published. Link between pages with relative links to `.md` files, and give
every figure alt text that says what it shows, not what it is called.

The [Playground](../playground.md) and the [Tree viewer](../tree.md) run the
command line compiled to WebAssembly, and that build is not committed. To try
either page locally, build it into `docs/assets/` first:

```bash
rustup target add wasm32-unknown-unknown
cargo build --release --target wasm32-unknown-unknown --manifest-path playground/Cargo.toml
cp playground/target/wasm32-unknown-unknown/release/karyon_playground.wasm docs/assets/
```

Two house rules hold for the documentation, the code and its comments alike.
**English throughout, with British spelling**, as the code and its `--help`
already use: colour, behaviour, organise. And **no em-dash characters
anywhere** (U+2014), which is easy to check without typing one:

```bash
grep -rn "$(printf '\xe2\x80\x94')" docs src examples README.md CHANGELOG.md
```

That has to find nothing. A comma, a colon or a full stop says the same thing.

## Adding a track type

[Writing a track](../how-it-works/extending.md) has the trait, a complete
example and the test a track has to pass first: whether it lives on the
coordinate axis at all. A track that the crate itself ships then needs:

- **Its own file** in `src/track/`, exported from `src/track/mod.rs` and
  re-exported from `src/lib.rs`.
- **An entry in the `tracks!` list in `src/plot.rs`**, which is what makes
  `Plot::label` and `Plot::adjust` work on it, and an `add_` method beside the
  others. A track that lays an array along the axis also gets an `_at` form.
- **Tests**, including one that pins a known coordinate to a known pixel. Track
  arithmetic is where the silent errors live.
- **A figure**: a panel in `examples/figures/gallery.rs`, or an example of its
  own, rendered into `assets/` and copied to `docs/assets/figures/`.
- **Documentation**: an entry on its family's page of the
  [track catalogue](../tracks/index.md), and a line in `CHANGELOG.md` that
  says why the track exists, not only that it was added.

A reader for its format belongs in `src/read/` only if the format is
line-based text, and it takes a `&str` rather than a path. It wants a fixture
in `src/read/audit.rs` on the same base as every other format. For the command
line to reach the track, it needs a `Kind` in `src/cli/args.rs` and a builder
in `src/cli/stack.rs`.

## Conventions the code follows

- **Coordinates are 0-based and half-open**, everywhere, with the two
  exceptions a reader sees: locus strings and the numbers printed on a figure.
  A constructor that takes anything else says so in its documentation. See
  [Coordinates](../how-it-works/coordinates.md).
- **The work follows the pixels.** Where the data are denser than the output,
  a track reduces them to what the pixels can show: one value per pixel
  column, merged runs, or a hint, the way `CoverageTrack`, `MsaTrack` and
  `SequenceTrack` do. [Scale](../how-it-works/scale.md) shows how each track
  does it, and which point tracks still draw a mark per datum.
- **Output is deterministic.** The same input renders byte-identical output:
  there is no clock, no hash iteration order and no random tie-break, and the
  examples generate their data from fixed seeds. That is what makes
  `git diff -- assets` a test at all.
- **Categories are coloured in order of first appearance**, not by hashing
  their names, so adding a sample adds a colour at the end rather than
  reshuffling the ones already there.
- **Plain SVG 1.1.** No scripts, no external references and no embedded fonts,
  so a file opens unchanged in a browser, in Inkscape and in Illustrator.
- **Errors are values.** A bad locus string is an `Error`, not a panic, and it
  converts into `std::io::Error`; a file a reader cannot parse stops on the
  line that failed, and the `ReadError` says which.
- **No dependencies and no unsafe code, and the library reads no files.** Both
  dependency tables in `Cargo.toml` are empty, `unsafe` is forbidden, and every
  reader takes a `&str`, so where the text comes from is the caller's
  decision.
- **Every public item is documented.** `missing_docs` is a warning, and the
  checks deny warnings, so an undocumented item fails the build.
- **The oldest supported compiler is 1.74**, the `rust-version` in
  `Cargo.toml`.

??? info "Why colour follows first appearance"
    A figure that recolours itself when a sample is added cannot go in a paper:
    the caption written against last week's version is now wrong, and nothing
    announces it. `VariantTrack::categories` walks the variants and collects the
    categories in the order it meets them, and the palette is indexed by that
    position. The cost is that two figures meant to agree on what a colour
    means have to be given their data in one order.

## Where next

<div class="grid cards" markdown>

-   **[Writing a track](../how-it-works/extending.md)**

    The trait a new track implements, and what the figure does for it.

-   **[Changelog](changelog.md)**

    What changed in each release, and why.

-   **[Citation](citation.md)**

    What to record about the version you used.

</div>
