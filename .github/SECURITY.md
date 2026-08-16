# Security Policy for karyon

## What a vulnerability means here

`karyon` reads text and writes an SVG document. It has no runtime dependencies
and no development dependencies, it opens no sockets, it runs no subprocess,
and the library opens no files at all: every reader in `karyon::read` takes a
`&str`, so the program embedding the crate is the one that touches the disk.
That leaves a small and specific attack surface, and it is worth naming it
precisely, because a report that lands in one of these classes gets handled as
a vulnerability and a report that does not is a plain bug that belongs in a
public issue.

The inputs that come from someone else are VCF, SAM, GFF3, BED, bedGraph,
cytoBand, `samtools depth`, FASTA and Newick files, and the strings inside them:
sequence names, feature names, sample names, tip labels. The output goes into a
document that a browser will render. So there are three realistic classes.

### 1. A hang or unbounded work on hostile input

A figure that never comes back is a denial of service against whatever renders
it, and this is not hypothetical. `SquiggleTrack` never finished drawing over a
wide region: the envelope branch walked the axis positions a pixel column
covers rather than the samples in it, so over a chromosome each column asked
thousands of bases whether they happened to be a sample, and over the whole
coordinate range it asked billions. Ordinary data plus a wide region was enough
to reach it. Chosen data would have been worse. It is fixed and recorded in
`CHANGELOG.md`.

The rule the crate holds itself to is that a track handed four million points
must not emit four million elements, and the work must scale with the pixels or
with the samples, never with the coordinate range. An input that breaks that
rule is in scope.

### 2. A panic on hostile input

A panic in a library is a crash in the caller. Coordinate arithmetic is where
these live, because a position near the top of the `u64` range does not survive
a trip through `f64`, and everything downstream of that overflows. Six of them
came out of the property suite on its first run, in `Scale::bounds`, in
`Read::end` and the CIGAR walk, in the codon walk of `OrfTrack` and at the right
edge of the sequence track. `Genome` accumulated lengths without saturating in
seven places. A large finite `Figure::width` panicked with a capacity overflow.
All of those are fixed and each has a test that pins it.

Note that release builds set `overflow-checks = true` on purpose. A wrapped
coordinate is worse than a panic, because it draws: `Genome::boundaries` came
back decreasing and every point in a genome-wide plot landed on the wrong
sequence with nothing said. So the intended behaviour of arithmetic that goes
out of range is a saturating value where the drawing has a sensible answer, and
a loud abort where it does not. An input that produces a wrong figure instead is
the more serious finding of the two.

### 3. An escaping failure

Names come from a user's file and go into the document, as `<title>` text and
inside attributes. `svg::escape` replaces the five XML metacharacters and drops
the characters XML 1.0 has no way to write at all, and nothing reaches the
output without passing through it. If a name can close an element, open one of
its own, or break out of an attribute, that is a vulnerability rather than a
formatting problem, because the file is a document and an attacker choosing its
elements chooses what a viewer executes.

Both halves of that have already gone wrong once. `--color` was pasted into an
attribute unescaped, which is attribute injection: a value carrying `<` put
elements of the caller's choosing into the document, and metacharacters are
refused there now. Separately, `escape` let XML-forbidden control characters
through, so a name read from a user's file could produce a document no parser
accepts.

`tests/properties.rs` is the standing check on this class. It generates ten
thousand figures per run from a seeded LCG, and three names in ten of the ones
it feeds in are hostile on purpose: `a<b&c"d'e`, `<script>alert(1)</script>`,
and a bidi override next to a raw tab. Every rendered document is then read back and must
have no raw `<` in any title, no `&` that does not open one of the five entities
the writer emits, no control byte anywhere, no non-finite or negative number in
a numeric attribute, and no `url(#...)` that resolves to nothing. A report in
this class is most useful when it can be turned into another one of those
invariants, so send the input rather than only the conclusion.

### What is not a vulnerability

A figure that is drawn wrong, a track that clips, a parse error on a malformed
file that is correctly reported as an error, a locus string that is rejected: all
bugs, all welcome, all belong in a public issue with the `bug` label. The
project treats a wrong figure as seriously as a crash, so nothing is lost by
reporting it in the open.

## Reporting a vulnerability

Please do not open a public issue for anything in the three classes above.
Email <paula.ruiz.rodriguez@csic.es> instead. When private vulnerability
reporting is enabled on the repository, the **Report a vulnerability** button on
the [Security tab](https://github.com/PathoGenOmics-Lab/karyon/security) reaches
the same people; email is the route that always works.

A report is actionable when it carries:

- **The input, or a way to make it.** Every example under `examples/` generates
  its own data from a fixed seed, and `tests/properties.rs` prints the seed of
  any figure that fails, so a seed is often the whole reproduction. If real data
  triggered it, synthetic input of the same shape almost always reproduces it,
  since the positions, lengths, strands and CIGARs are what the drawing depends
  on and the bases rarely are.
- **The version.** `karyon --version` for the command, the `Cargo.lock` entry or
  the git commit for the library.
- **What happens.** The panic message and backtrace, the region and the wall
  clock time for a hang, or the fragment of the SVG that escaped.

We ask that you give us the chance to ship a fix before publishing the details.
You will get an acknowledgement that a human read the report, and a fix for
anything in class 1 or 3 is treated as urgent.

## Supported versions

`karyon` is pre-1.0 and is not published on crates.io yet, so there are no
long-lived release branches to backport into and no packaged builds in anyone
else's distribution channel to coordinate with. Fixes land on `main` and go out
in the next release. Use the latest release, or `main`.

## What a fix looks like

Every fix in these classes ships with a test that pins the answer, and where the
failure can be stated without knowing the right answer it becomes a property in
`tests/properties.rs` rather than an example beside the code. That distinction is
the reason the suite exists: a test written next to the code inherits the blind
spot of whoever wrote it, because the person writing the assertion is the person
who believed the code. The escaping rules, the finiteness rules and the
coordinate arithmetic are all pinned that way.

Fixes and hardening changes carry the `security` label, which also keeps them
out of the stale sweep.

## Notes for anyone embedding karyon

- **There is no dependency tree to audit.** `[dependencies]` and
  `[dev-dependencies]` are both empty, and that is a property of the crate
  rather than a coincidence of the current version. No advisory database has
  anything to say about a build of this library.
- **The trust boundary is yours.** The library never opens a path. You read the
  bytes, you decide how many of them to read, and you hand a `&str` to a reader
  in `karyon::read`. A file too large to hold in memory is your call to make,
  not the crate's.
- **The output is plain SVG 1.1**, with no scripts, no external references and
  no embedded fonts, so nothing in a rendered figure fetches anything when it is
  opened.
- **Serve a rendered figure as an image, not as inline markup**, if the names in
  it came from somebody else. `image/svg+xml` in an `<img>` or as a separate
  request keeps the document out of your page's origin. The escaping is meant to
  make inlining safe, and this is the layer that holds if it ever does not.

## Contact

- Email: <paula.ruiz.rodriguez@csic.es>
- The repository's
  [security overview](https://github.com/PathoGenOmics-Lab/karyon/security) for
  the current posture and any published advisories.
