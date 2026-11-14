---
description: >-
  Install the karyon program or add the library to a Rust project, set up tab
  completion and the manual page, and pin a commit for reproducible figures.
---

# Installation

Install the `karyon` program with one command, or add the library to a Rust
project with one line. Either way, a Rust toolchain is all you need. This page
has the details; [Start here](../start.md) has the short way.
{ .k-lead }

## The fast path

The command line, built and put on your `PATH` by Cargo:

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
```

The library, as one line in your project's `Cargo.toml`:

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

!!! note "Not on crates.io yet"
    Until it is published, `cargo install karyon` and `cargo add karyon` will
    not find it, which is why both lines point Cargo at the repository.

!!! tip "Try it without installing anything"
    The [Playground](../playground.md) runs the real command line in your
    browser, compiled to WebAssembly. Nothing you type is uploaded.

## Check that it works

`cargo install` puts the program in `~/.cargo/bin`, which rustup adds to your
`PATH`:

```bash
karyon --version   # prints karyon and the version number
karyon --help      # what to type, on one screen
karyon help coverage   # one track: what it reads and the options it takes
```

For the library, this small program builds a figure and reports its size:

```rust
fn main() -> std::io::Result<()> {
    let svg = karyon::plot("chr1:1-1000")?
        .add_coverage(vec![30.0; 1000])
        .to_svg();
    println!("{} bytes of SVG", svg.len());
    Ok(())
}
```

A project that depends on the library never builds the command line program,
so adding karyon adds nothing else to your dependency tree.

## Tab completion and the manual page

The program writes its own completion script for bash, zsh and fish, and its
own manual page, from the tables its parser and its help are built from, so
each matches the version that wrote it. A script completes the flags, the
words a flag takes, as `--theme light`, and the names of files after a flag
that reads one. No release has been published yet; each release archive
will carry all four already written, the scripts in `completions/` and the
page in `man/`.

=== "bash"

    ```bash
    mkdir -p ~/.local/share/bash-completion/completions
    karyon completions bash > ~/.local/share/bash-completion/completions/karyon
    ```

    That is where the bash-completion package looks. Without it, load the
    script from `~/.bashrc` instead:

    ```bash
    eval "$(karyon completions bash)"
    ```

=== "zsh"

    ```zsh
    mkdir -p ~/.zfunc
    karyon completions zsh > ~/.zfunc/_karyon
    ```

    Then in `~/.zshrc`, before `compinit` runs, or in place of the line that
    runs it:

    ```zsh
    fpath=(~/.zfunc $fpath)
    autoload -Uz compinit && compinit
    ```

=== "fish"

    ```fish
    mkdir -p ~/.config/fish/completions
    karyon completions fish > ~/.config/fish/completions/karyon.fish
    ```

    fish loads it from there by itself, the first time a shell completes
    karyon.

The manual page goes where `man` looks for one: beside each folder on your
`PATH`, in `../share/man`. For the program Cargo put in `~/.cargo/bin`:

```bash
mkdir -p ~/.cargo/share/man/man1
karyon man > ~/.cargo/share/man/man1/karyon.1
man karyon
```

A command is only ever the first word, where the place goes, so a gene or a
sequence called `help`, `completions` or `man` is still drawn: write it after
a file, as `karyon genes.gff3 man`.

## Get a Rust toolchain

Skip this if `rustc --version` already prints 1.74 or newer.

=== "rustup (recommended)"

    ```bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    ```

=== "conda"

    ```bash
    conda install -c conda-forge rust
    ```

## Pin a commit for reproducible figures

A git dependency follows the repository's default branch, so `cargo update`
can move you to newer code. When a figure has to come out the same next year,
pin the commit you tested:

```toml
[dependencies]
# <sha> is the commit you tested against, full or short.
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon", rev = "<sha>" }
```

The command line takes the same pin:

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon --rev <sha>
```

Rendering is deterministic: the same input and the same version give a
byte-identical SVG or PDF, so a pinned figure can be checked with a plain
`diff`.

## Build from a clone

```bash
git clone https://github.com/PathoGenOmics-Lab/karyon
cd karyon
cargo build --release          # the program is target/release/karyon
cargo test                     # the full test suite
cargo run --quiet -- --help    # run the command line without installing it
```

Nearly every figure on this site comes from a program in `examples/`, and each
one takes the directory to write into:

```bash
cargo run --example locus -- assets
```

The code that builds each of those figures is in `examples/figures/`, and the
site compiles the same files into the page, so the figures you see here are
drawn in your browser rather than shown: in the page's light or dark, at the
width of the column. A figure over a stretch of genome can be dragged along it
and zoomed into, down to the bases, with the buttons under it, a pinch, or the
wheel with ctrl or cmd held; the others zoom as a picture. **SVG** under a
figure saves the view you are looking at.

## Requirements

| | |
|:--|:--|
| Rust | 1.74 or newer (the MSRV), edition 2021. An older toolchain stops with a message naming the version it needs. |
| Dependencies | None. The `[dependencies]` table in `Cargo.toml` is empty. |
| System libraries | None: no cairo, fontconfig, OpenSSL, Python or headless browser. |
| Platforms | Any target with Rust's standard library, WebAssembly included. Linux, macOS and Windows are built and tested on every pull request. No release has been published yet; when one is, its Windows archive will carry its own C runtime, so there will be no Visual C++ redistributable to install. |
| Input | Text formats, gzipped or not, and BAM, BCF, bigWig, bigBed, 2bit and Juicer's `.hic` read directly through their indexes. CRAM comes in through `samtools`, piped. |
| Output | Standalone SVG, which names its fonts rather than embedding them, or PDF, which embeds a subset of each [Liberation face](../how-it-works/pdf.md#text) it sets its text in, under the SIL Open Font License. A name ending in `.pdf` writes PDF. |

## Next

<div class="grid cards" markdown>

-   **[Start here](../start.md)**

    A figure from the shell in three steps.

-   **[A first program](quickstart.md)**

    The same kind of figure from Rust.

-   **[Core ideas](concepts.md)**

    Regions, tracks, the shared scale and the builder that stacks them.

-   **[Command line](../guide/cli.md)**

    Every flag the `karyon` command takes.

-   **[Gallery](../plots/index.md)**

    Every kind of plot karyon draws, found by what you want to show.

</div>
