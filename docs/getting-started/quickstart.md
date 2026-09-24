# Quickstart

Two short routes to a first figure: one from the shell, one from Rust. Both
draw a stack of tracks over one shared coordinate axis and write it as a
standalone SVG. Pick the one that suits you; the shell needs no Rust code at
all.

## From the shell

Install the command line first (see [Installation](installation.md)):

```bash
cargo install --git https://github.com/PathoGenOmics-Lab/karyon
```

### A figure in two lines

This makes its own input file, so you can run it right now:

```bash
seq 1 60 | awk '{print "chr1\t" $1 "\t" (20 + $1 % 7)}' > depth.txt
karyon chr1:1-60 --coverage depth.txt --label depth -o first.svg
```

Open `first.svg` in a browser. The first argument is the region, `--coverage`
adds a track that reads `depth.txt`, `--label` names it, and `-o` is the output
file. Without `-o` the SVG goes to standard output.

### Several tracks

Each track flag, such as `--coverage` or `--variants`, starts a new track, and
the flags after it describe that track. **The order of the flags is the order of
the tracks, from top to bottom:**

```bash
karyon NC_000962.3:761,121-761,180 \
  --coverage depth.txt  --label depth     --height 45 \
  --sequence H37Rv.fa   --label reference \
  --variants calls.vcf  --label variants  --height 40 \
  --title 'The same locus at base resolution' \
  -o rpoB-zoom.svg
```

![The same locus over sixty bases: a depth profile, the reference sequence drawn as coloured letters, three variant lollipops over the bases they change, and a ruler counting single bases](../assets/figures/example-zoom.svg){ width="900" height="223" loading="lazy" }

- `depth.txt` is what `samtools depth` writes, `H37Rv.fa` is the reference
  FASTA and `calls.vcf` is a VCF.
- Each file is read in its own format's coordinates (BED and bedGraph are
  0-based, GFF3, VCF, SAM and `samtools depth` are 1-based), and all of them end
  up at the same place in the figure.
- The variant colours come from the `ANN` or `BCSQ` consequence when the VCF has
  one, and otherwise from the call itself: a substitution, an insertion or a
  deletion.
- `karyon --help` lists every track flag and option.

!!! note "Give `--sequence` the whole reference"
    `--sequence` wants the FASTA the window is cut out of, not just the bases on
    display. A FASTA holding only those sixty bases would put them at positions 0
    to 59, far from the window, and the track would come out empty.

### BAM, CRAM and BCF

karyon reads text formats only. For binary files, let `samtools` or `bcftools`
write the text and pipe it in; any track file can be `-` for standard input:

```bash
samtools depth -a -r NC_000962.3:761000-763000 aln.bam \
  | karyon NC_000962.3:761,000-763,000 --coverage - --label depth -o rpoB.svg
```

## From Rust

Add the dependency (see [Installation](installation.md)):

```toml
[dependencies]
karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon" }
```

This is a whole program:

```rust
use karyon::{plot, Aggregate, Feature, Strand, Variant};

fn main() -> std::io::Result<()> {
    // One value per base of the window, however you got it. The dip is the
    // shape a deletion leaves behind.
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
            Variant::new(761_108).value(0.98).category("missense"),
            Variant::new(761_154).value(1.00).category("missense"),
            Variant::new(761_155).value(0.21).category("synonymous"),
        ])
        .label("variants")
        .save("rpoB.svg")?;

    Ok(())
}
```

`cargo run` writes `rpoB.svg`:

![A stack of tracks over two kilobases of the rpoB locus: a depth profile with a dropout in it, a reference sequence too zoomed out to show its letters, the gene with the resistance determining region marked inside it, variant lollipops coloured by consequence, and a coordinate ruler underneath](../assets/figures/example.svg){ width="900" height="306" loading="lazy" }

What the program does:

- **One call per track, in the order they stack.** `add_coverage` is at the top
  because it comes first. Tracks only say how tall they are; the figure stacks
  them and gives every one the same horizontal scale, which is what keeps them
  aligned.
- **The region is written once.** `add_coverage` lays its array from the left
  edge of the region. When your array starts somewhere else, use
  `add_coverage_at(start, values)`.
- **`label` and `adjust` apply to the track just added.** `label` names it in
  the left margin; `adjust` gives you the track itself to configure. Here it asks
  for `Aggregate::Min`, so the dropout survives when several bases share a
  pixel.
- **The ruler is added for you** at the bottom. `add_axis` puts one elsewhere
  and `remove_axis` leaves it out.
- **Errors are `io::Error`.** A bad region string is an error, not a panic, so
  `?` works in a function that returns `io::Result`.
- **Detail follows the zoom.** Over two thousand bases the sequence track shows
  a hint instead of unreadable letters; zoom in to sixty bases and the letters
  appear, with no other change to the code.

`to_svg()` returns the SVG as a string instead of writing a file, which is what
a web service or a test wants.

!!! warning "Coordinates"
    The region string and the tick labels are 1-based and inclusive, like
    `samtools` and IGV. Everything else in the API is **0-based and
    half-open**, like BED: `rpoB` is 759,807 to 763,325 in the 1-based
    annotation, so it is `Feature::new(759_806, 763_325)` here, and the call at
    VCF `POS` 761,155 is `Variant::new(761_154)`. The readers subtract the one
    for you.

The two figures on this page come from the `locus` example, which is this
program with seeded random data: `cargo run --example locus -- assets`.

## Next

- [Core concepts](concepts.md): regions, tracks and the shared scale.
- [Gallery](../plots/index.md): find the plot for your data.
- [Writing a figure in Rust](../guide/plot.md): everything `plot()` can do.
- [Command line](../guide/cli.md): every flag the `karyon` command takes.
