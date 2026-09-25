# Your first figure

Point karyon at the files you have and a place, and it draws them on one axis.
The first two steps use your own files; the ones after make their own, so all
you need for them is karyon itself.
{ .k-lead }

<div class="k-steps" markdown>

### Install the command line

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
karyon --version
```

This needs a Rust toolchain. [Installation](installation.md) covers getting
one, and adding karyon to a Rust project instead.

### Draw the files you have

Name a place, then the files, then where the figure goes:

```bash
karyon rpoB reads.bam genes.gff3 calls.vcf.gz -o rpoB.svg
```

Open `rpoB.svg` in a browser: the depth of the reads, the genes and the calls,
over the gene rpoB with a margin either side, and a ruler underneath.

| Part | What it does |
|:--|:--|
| `rpoB` | The place, always first. A gene the annotation names, as here; a sequence's name, drawn whole; or a region such as `NC_000962.3:761,000-763,000`, 1-based and inclusive as `samtools` and IGV write it. |
| `reads.bam` | A track for each file, of the kind its name says. A BAM is drawn as its depth, read through the `.bai` beside it, so only the reads over the place are read. |
| `genes.gff3` | The genes, each drawn once and named. GTF and BED work the same way. |
| `calls.vcf.gz` | The calls. A file compressed with gzip or bgzip is read as the file inside. |
| `-o rpoB.svg` | The output file. Without it the SVG goes to standard output, ready for a pipe. |

Each track is called after its file in the left margin, and the depth of a
BAM after the file and `depth`, as `reads depth`. To choose the kind
yourself, put its flag in front: `--pileup reads.bam` draws the reads rather
than their depth. The options for a track go after it: `reads.bam --label depth
--height 60`. `karyon help pileup` lists what one track takes, and `karyon
--help` fits on one screen.

A file whose name does not say what it holds, a `.tsv` or a `.txt`, needs its
track's flag in front. When a place is not in any of the files, karyon says so
and names what each file does hold.

### Draw a scan or a tree

An association tool's table, and a tree with a sheet of what is known about
its samples, are drawn the same way:

```bash
karyon 1 gwas.assoc --threshold genome-wide -o scan.svg
karyon tree.nwk --traits samples.tsv --columns lineage -o tree.svg
```

- **A place is named as the files name it.** PLINK writes a chromosome as `1`
  where a FASTA may call it `NC_000962.3`, so a PLINK table is placed on `1`.
  A sequence no file gives the length of is drawn as far as the rows reach;
  write the span, as `1:1-4,411,532`, to draw all of it.
- **The table is read by its header.** PLINK, PLINK 2, REGENIE, BOLT-LMM,
  GEMMA, SAIGE and the GWAS Catalog each name their columns, and a p-value is
  drawn as -log10 of itself. `--threshold genome-wide` draws the line at
  p = 5e-8, and `--threshold 1e-5` wherever you say.
- **A figure of trees alone needs no place.** `--traits` puts the sheet's
  columns beside the tips, every column unless `--columns` names the ones to
  draw, and a key under the figure names each colour.
- **One position is not a place.** Around a variant, write a span, as
  `NC_000962.3:761,000-761,400`; given `NC_000962.3:761,200`, karyon answers
  with that span.

### Or make a figure from nothing

The first line writes a small depth file, so there is nothing to download. The
second draws it:

```bash
seq 1 60 | awk '{print "chr1\t" $1 "\t" (20 + $1 % 7)}' > depth.txt
karyon chr1:1-60 --coverage depth.txt -o first.svg
```

`--coverage` says what `depth.txt` is, since a `.txt` could hold anything: three
columns as `samtools depth` writes them. The track is called `depth`, after its
file.

### Stack several tracks

Write a sixty-base reference and a VCF with three calls, then draw all three
files on the one axis:

```bash
printf '>chr1\n%s\n' GATTACAGGCTTACCGATCGATGCAAGCTTGGCCATTAGCGATCCGTAAGCTTACGGATC > ref.fa
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > calls.vcf
printf 'chr1\t%s\t.\t%s\t%s\t.\tPASS\tAF=%s\n' 18 T C 0.95 31 G A 0.4 47 T TA 0.6 >> calls.vcf

karyon chr1:1-60 \
  --coverage depth.txt --height 45 \
  ref.fa --label reference \
  calls.vcf --height 40 \
  --title 'Three tracks over sixty bases' \
  -o stack.svg
```

**The order of the words is the order of the stack.** Each file, or each
track flag such as `--coverage`, starts a new track, and the options after it
(`--label`, `--height`) describe that track until the next one starts. Move
the `calls.vcf` line above `--coverage` and the calls are drawn on top. Figure
options such as `--title` and `-o` belong to no track and can go anywhere.

<figure class="k-plate" markdown>
![Sixty bases of a locus: a depth profile, the reference drawn as coloured letters, three variant lollipops over the bases they change, and a ruler counting single bases](../assets/figures/example-zoom.svg){ width="900" height="222" loading="lazy" }
<figcaption>The same three tracks over sixty bases of the rpoB gene, drawn by the crate's <code>locus</code> example.</figcaption>
</figure>

- **Every file keeps its own coordinates.** VCF, GFF3, SAM and `samtools depth`
  count from 1; BED and bedGraph count from 0. Each is converted on the way in,
  so all of them land in the same place.
- **Variants are coloured by what they are.** The colour comes from the
  consequence an annotator wrote in the VCF (`ANN` or `BCSQ`), or, without one,
  from the call itself: here two substitutions and an insertion. The height of
  each lollipop follows its `AF`, and a call without one stands full height.
- **Letters follow the zoom.** A base is printed as a letter while it is at
  least 7 pixels wide and as a coloured block down to 0.6 of a pixel, with a
  key under the figure naming the colours. Below that the track prints a hint
  to zoom in rather than a smear.
- **Reads are compared with the reference the figure draws.** A `--pileup`
  under a `--sequence` paints every base that differs from it.

!!! note "Give `--sequence` the whole reference, or a slice that says where it is"
    `--sequence` takes the only record of the FASTA, or the one named like the
    region's sequence when it holds several, and cuts the region out of it by
    position. A record is read from base 1 unless its header says otherwise,
    as the `>chr1:101-160` that `samtools faidx` writes does, so give it the
    whole sequence or such a slice. A record with no base in the window is
    refused, saying which bases it holds. Here the two are the same, because
    the region starts at base 1.

CRAM, BCF and bigWig are not read directly: hand a track one and it says which
`samtools`, `bcftools` or UCSC command writes what it reads, to put in place of
the file's name.

### Draw the same kind of figure from Rust

Make a new project and add the library:

```bash
cargo new first-figure
cd first-figure
```

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

Replace `src/main.rs` with this whole program:

```rust
use karyon::{plot, Aggregate, Feature, Strand, Variant};

fn main() -> std::io::Result<()> {
    // One depth value per base of the region, with a dropout in the middle:
    // the shape a deletion leaves behind.
    let depth: Vec<f64> = (0..2_000)
        .map(|i| if (900..1_030).contains(&i) { 3.0 } else { 55.0 - (i % 23) as f64 })
        .collect();
    let bases: Vec<u8> = b"ACGT".iter().cycle().take(2_000).copied().collect();

    plot("NC_000962.3:761000-762999")?
        .title("rpoB locus, resistance determining region")
        .add_coverage(depth)
        .label("depth")
        .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
        .add_sequence(bases)
        .label("reference")
        .add_features(vec![
            Feature::new(759_806, 763_325)
                .name("rpoB")
                .strand(Strand::Forward),
            Feature::new(761_081, 761_162)
                .name("RRDR")
                .strand(Strand::Forward)
                .color("#d55e00"),
        ])
        .label("annotation")
        .add_variants(vec![
            Variant::new(761_051).value(0.12).category("synonymous"),
            Variant::new(761_108).value(0.98).category("missense"),
            Variant::new(761_138).value(0.55).category("missense"),
            Variant::new(761_154).value(1.00).category("missense"),
            Variant::new(761_155).value(0.21).category("synonymous"),
        ])
        .label("variants")
        .save("rpoB.svg")?;

    Ok(())
}
```

`cargo run` writes `rpoB.svg`:

<figure class="k-plate" markdown>
![A stack of tracks over two kilobases of the rpoB locus: a depth profile with a dropout in it, a reference sequence too zoomed out to show its letters, the gene with the resistance determining region marked inside it, variant lollipops coloured by consequence, and a coordinate ruler underneath](../assets/figures/example.svg){ width="900" height="305" loading="lazy" }
<figcaption>The <code>locus</code> example is this program with seeded, more realistic depth and bases: <code>cargo run --example locus -- assets</code>.</figcaption>
</figure>

- **One call per track, in the order they stack.** `add_coverage` comes first,
  so depth is on top. Every track gets the same horizontal scale from the
  figure, which is what keeps them aligned.
- **The region is written once.** `add_coverage` and `add_sequence` lay their
  arrays from the left edge of the region. When an array starts somewhere
  else, use `add_coverage_at(start, values)` or `add_sequence_at(start, seq)`.
- **`label` and `adjust` apply to the track just added.** `label` names it in
  the left margin; `adjust` hands you the track itself to configure. Here it
  asks for `Aggregate::Min`, so the dropout survives when several bases share
  a pixel.
- **The ruler is added for you** at the bottom. `add_axis` puts it somewhere
  else and `remove_axis` leaves it out.
- **A bad region is an error, not a panic.** `plot` reports it as a
  `karyon::Error`, which converts into `std::io::Error`, so one `?` covers
  both the region and the file write.
- **`to_svg()` returns the document as a `String`** instead of writing a file,
  which is what a web service or a test wants.

Over two thousand bases the reference is too narrow for letters and prints a
hint instead. [Zooming to base resolution](../recipes.md#zooming-to-base-resolution)
draws the same stack over sixty bases, where every letter shows.

!!! warning "Two coordinate conventions"
    The region string and the tick labels are 1-based and inclusive, like
    `samtools` and IGV. Everything else in the API is 0-based and half-open,
    like BED. rpoB is 759,807 to 763,325 in the annotation, so it is
    `Feature::new(759_806, 763_325)` here, and the call at VCF `POS` 761,155
    is `Variant::new(761_154)`. The file readers subtract the one for you;
    [Coordinates](../how-it-works/coordinates.md) has the full story.

### Where to go next

<div class="grid cards" markdown>

-   **[Core ideas](concepts.md)**

    Regions, tracks, the shared scale and what karyon refuses to draw.

-   **[Gallery](../plots/index.md)**

    Find the plot for your data by what you want to show.

-   **[Recipes](../recipes.md)**

    Short, complete answers to common tasks, from the shell and from Rust.

-   **[The Rust API](../guide/plot.md)**

    Everything `plot()` and `Figure` can do.

-   **[Command line](../guide/cli.md)**

    Every flag the `karyon` command takes.

</div>

</div>
