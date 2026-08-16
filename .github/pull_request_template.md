<!--
Fill in what applies and delete what does not. A short pull request needs a
short description; nobody is asking for an essay to fix a typo.
-->

## What this changes

<!-- One or two sentences. What is different after this is merged? -->

## Why

<!-- The problem it solves. Link an issue with "Closes #123" if there is one. -->

## How it was verified

<!--
The important part, and the one a reviewer cannot reconstruct.

Not "it should work", but what you ran and what it printed. For a fix, the test
that pins it and the fact that you watched it fail before the change: a test
written after the fix and never seen red is a test of nothing, and this crate
has had one of those. For a change to the drawing, which figures under assets/
moved and what you saw when you opened the old and the new one side by side.

If CI already covers it, say which of the two jobs does.
-->

## What moved in the public API

<!--
Nothing, or the items added, renamed or removed. karyon is a library before it
is a command, it is pre-1.0 with no crates.io release, and callers pin a commit,
so a rename is free today and expensive the week after publication. Say it here
so the decision is made on purpose.

New public items need documentation: missing_docs is a warning and cargo doc
runs with warnings denied, so an undocumented one fails the build rather than
being noticed later.
-->

## Checklist

- [ ] `cargo test`, `cargo test --release`, `cargo fmt --all -- --check` and `cargo clippy --all-targets -- -D warnings` all pass.
- [ ] `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` is clean.
- [ ] New behaviour has a test, and a fix has one that I watched fail without the change, or I have said above why it does not.
- [ ] The figures were re-rendered and I looked at the diff, and `docs/assets/figures/` was refreshed with `cp assets/*.svg docs/assets/figures/`, or nothing here can reach the drawing.
- [ ] Documentation under `docs/` is updated if the change is user-visible, and `mkdocs build --strict` passes.
- [ ] `CHANGELOG.md` gained a line saying why the change exists, not only that it was made.
- [ ] Both dependency tables in `Cargo.toml` are still empty, and nothing added needs a toolchain newer than the declared 1.74.
- [ ] English throughout, and `grep -rn "$(printf '\xe2\x80\x94')" docs src examples README.md CHANGELOG.md` finds nothing.
