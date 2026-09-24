# Installation

karyon is one Rust package that contains both a library and a command line
program. The only thing you need is a Rust toolchain: there is nothing else to
install, no C library to link against and no dependency to download.

!!! note "Not on crates.io yet"
    `cargo add karyon` and `cargo install karyon` will not find it yet, so the
    commands on this page point Cargo at the GitHub repository instead. Once it
    is published, those two commands will work and nothing else here changes.

## 1. Get a Rust toolchain

Skip this if `rustc --version` already prints 1.74 or newer.

=== "rustup (recommended)"

    ```bash
    curl https://sh.rustup.rs -sSf | sh
    ```

=== "conda"

    ```bash
    conda install -c conda-forge rust
    ```

## 2. Install the command line

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
```

This compiles the `karyon` program and puts it in `~/.cargo/bin`, which rustup
adds to your `PATH`. Check that it works:

```bash
karyon --version
karyon --help
```

`--help` prints the whole command grammar. The [Quickstart](quickstart.md) walks
through a first figure, and [Command line](../guide/cli.md) documents every flag.

## 3. Or use it as a library

Add karyon to your project's `Cargo.toml`:

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

A small program to check that it builds:

```rust
fn main() -> std::io::Result<()> {
    let svg = karyon::plot("chr1:1-1000")?
        .add_coverage(vec![30.0; 1000])
        .to_svg();
    println!("{} bytes of SVG", svg.len());
    Ok(())
}
```

Projects that use the library never build the command line program, so adding
karyon adds nothing else to your dependency tree.

### Pin a version for reproducible figures

A git dependency follows the default branch, so `cargo update` can move you to
newer code. When a figure has to come out exactly the same later, pin the commit
you tested:

```toml
[dependencies]
# <sha> is the commit you tested against, full or short.
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon", rev = "<sha>" }
```

Rendering is deterministic: the same input and the same version produce a
byte-identical SVG.

## Build from a clone

```bash
git clone https://github.com/PathoGenOmics-Lab/karyon
cd karyon
cargo build --release      # the program ends up in target/release/karyon
cargo test                 # the full test suite
```

The examples draw the figures used on this site. Each one takes an output
directory:

```bash
cargo run --example locus -- assets
```

## Requirements in detail

| | |
|:--|:--|
| **Rust** | 1.74 or newer, edition 2021. An older toolchain gets a clear message naming the version it needs. |
| **Runtime dependencies** | None. Both dependency tables in `Cargo.toml` are empty. |
| **System libraries** | None: no cairo, fontconfig, OpenSSL, Python or headless browser. |
| **Input formats** | Line-based text only. Convert BAM, CRAM and BCF with `samtools` or `bcftools` and pipe the text in. |
| **Output** | Standalone SVG 1.1, which names its fonts rather than embedding them. |

## Next

- [Quickstart](quickstart.md): a first figure from the shell and from Rust.
- [Core concepts](concepts.md): regions, tracks and the shared scale.
- [Gallery](../plots/index.md): every kind of plot karyon draws.
