# Contributing to karyon

The contributing guide lives on the documentation site:
**[Contributing](https://pathogenomics-lab.github.io/karyon/about/contributing/)**,
whose source is [`docs/about/contributing.md`](../docs/about/contributing.md).
It covers what a report of a wrong figure has to carry, the four gates a change
has to pass, how to add a track type, and the conventions the code follows.
Read it before you write code. This file is only the pointer GitHub looks for,
and the three things that catch people out.

## 1. Run the gates yourself

CI is `workflow_dispatch` only. Pushing does not start it, so nothing checks
your branch until somebody runs the workflow by hand. Run the same four
commands locally before you open a pull request:

```bash
cargo test
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
```

`cargo test --release` is worth a fifth run: some of the checks are about
floating point and layout arithmetic that the optimiser is allowed to rearrange,
and release builds also turn on overflow checks.

## 2. The figures are part of the build, not illustrations

Everything under `assets/` is rendered by an example, and rendering is
deterministic, so a figure that was not re-rendered is a diff. CI renders every
example and then fails if `git diff -- assets` is not empty. After any change
that could touch the drawing, re-render (the guide lists the commands), copy
the result into the site's own copy, and commit both:

```bash
cp assets/*.svg docs/assets/figures/
```

A diff in `assets/` is not a problem in itself. It is the review: open the old
figure and the new one side by side and check that what moved is what you meant
to move.

## 3. English throughout, and no em-dash characters

Prose, code and comments are all in English. The em-dash, U+2014, is not used
anywhere in this repository, because a comma, a colon or a full stop says the
same thing and survives every editor and every font. You can check without
typing one:

```bash
grep -rn "$(printf '\xe2\x80\x94')" docs src examples README.md CHANGELOG.md
```

That has to find nothing.

## If you are adding a track type

There is an entry test before any of the above, and it is the reason this crate
exists rather than a general plotting library: **does the track live on the
genomic coordinate axis?** If its `draw` never reads `ctx.scale`, its x is a
sample list or a category, and the plot is a bar chart, a line chart or a
heatmap that was handed genomic data, which other tools already draw better.
Three track types were removed under that rule. The guide has the rest of what
a new track needs.

## Reporting a problem

Open an issue at
[PathoGenOmics-Lab/karyon/issues](https://github.com/PathoGenOmics-Lab/karyon/issues).
A plotting library fails differently from a program that prints numbers, since
the figure still renders and looks fine, so a report needs the code or the
command in full, the version, and the SVG itself, which is text and attaches as
it is. Say what you expected to see and what you saw.

A hang, a panic or an escaping failure on input somebody else supplied is a
security report rather than an issue. See [SECURITY.md](SECURITY.md).

By participating you agree to the
[Code of Conduct](CODE_OF_CONDUCT.md).
