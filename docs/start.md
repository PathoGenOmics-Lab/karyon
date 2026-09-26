---
title: Start here
description: Install karyon and draw a figure of reads, variant calls and genes in three steps.
---

# Start here

Three steps, and all you need is a terminal.

## 1. Install karyon

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
```

`cargo` comes with Rust. If you do not have it, the one command at
[rustup.rs](https://rustup.rs) installs it.

## 2. Get the example files

Reads, variant calls, genes and the files of every other page, about 40
kilobytes in all:

```bash
curl -O https://pathogenomics-lab.github.io/karyon/data/examples.zip
unzip examples.zip
```

## 3. Draw them

```bash
karyon rpoB reads.bam genes.gff3 calls.vcf.gz -o rpoB.svg
```

<figure class="k-start" markdown>
![The depth of the reads over the gene rpoB, with a stretch no read covers, the gene as an arrow, and seven variant calls as lollipops as tall as their allele frequency](assets/start/reads.svg){ .k-light width="720" height="248" }
![The same figure on the dark page](assets/start/reads-dark.svg){ .k-dark width="720" height="248" }
</figure>

Open `rpoB.svg` in a web browser. Point at a call or a gene, and the browser
says what it is and where.

## How to read a command

1. **The place comes first:** a gene, a sequence, or a region such as
   `NC_000962.3:761,000-763,000`.
2. **Then your files.** Each one is a row, top to bottom in the order you
   write them.
3. **Then `-o` and the file to write.**

Options for one row go right after its file: `reads.bam --height 100`.

## Next

[What do you have?](your-data/index.md) has a page for each kind of data, with
the command and the figure it draws.
