# Core concepts

Six ideas explain almost everything karyon does. Read them once and the rest of
the documentation will make sense faster.

<div class="grid cards" markdown>

-   :material-map-marker-radius: **Region**

    ---

    Every figure is drawn over one region: a sequence name and a span, such as
    `chr1:1,000-2,000`. You name it once and every track is drawn on it.

    [:octicons-arrow-right-24: Coordinates](../how-it-works/coordinates.md)

-   :material-view-sequential: **Track**

    ---

    A track is one horizontal band: coverage, genes, variants, reads, a tree.
    Tracks stack top to bottom in the order you add them.

    [:octicons-arrow-right-24: Track reference](../tracks.md)

-   :material-ruler: **Scale**

    ---

    One scale turns positions into pixels, and every track uses the same one.
    That is why the rows line up without being aligned by hand.

    [:octicons-arrow-right-24: Scale](../how-it-works/scale.md)

-   :material-code-braces: **Plot and Figure**

    ---

    `plot()` is the short way to build a figure, one call per track. `Figure` is
    the lower level it builds, for when you need more control.

    [:octicons-arrow-right-24: Writing a figure in Rust](../guide/plot.md)

-   :material-earth: **Genome, Rings and Panels**

    ---

    Containers for what one region cannot hold: several sequences end to end, a
    circular chromosome, and several drawings on one sheet.

    [:octicons-arrow-right-24: Whole genomes and geography](../plots/whole-genomes-geography.md)

-   :material-file-document-outline: **Readers**

    ---

    `karyon::read` turns text formats such as BED, VCF and GFF3 into the values
    tracks take. The command line uses the same readers.

    [:octicons-arrow-right-24: File formats](../guide/formats.md)

</div>

## Region

A region is a sequence name and a span. It is the first thing a figure is built
from, and it is never empty.

```rust
let figure = karyon::plot("NC_000962.3:761,000-762,999")?;
```

The region string is written the way `samtools` and IGV write it, 1-based and
inclusive, and commas in the numbers are ignored. The name is only a label:
karyon never looks it up or compares it with a file, so a region can just as
well count the columns of an alignment (`alignment:1-320`) or the samples of a
nanopore signal.

## Track

A track owns one horizontal band. It knows how tall it wants to be and how to
draw itself, and nothing else: it is never told where it sits or what the other
tracks are. The figure stacks the bands in the order you add them and gives every
one the same horizontal axis.

That is what lets you read a figure by running a finger down it. Everything a
finger crosses at one position is the same position in every band.

There are 36 track types. The [gallery](../plots/index.md) sorts them by what
you want to show, and the [track reference](../tracks.md) lists each one with its
options.

## Scale

A figure builds one scale and hands it to every track. The scale knows where the
region starts, how long it is, and where and how wide the plotting area is.

Two consequences are worth knowing:

- **Nothing is clamped.** A gene that starts before the window is drawn as the
  whole shape it is, and the edge of the figure crops it. Clamping would move its
  end to the border and make it look shorter than it is.
- **A pixel can cover many bases.** Across a whole genome one pixel may stand for
  thousands of bases, so a track has to choose what to show for them. Coverage,
  for example, can show the highest, lowest or mean value in each pixel. The
  lowest is what keeps a dropout visible:

    ```rust
    .add_coverage(depth)
    .adjust(|track| track.aggregate(Aggregate::Min))
    ```

## Plot and Figure

`plot()` is the quickest way to write a figure: one `add_` call per track, in the
order they stack, with the ruler added at the bottom for you.

```rust
use karyon::plot;

plot("chr1:1-1000")?
    .add_coverage(vec![30.0; 1000]).label("depth")
    .save("depth.svg")?;
```

It builds a `Figure`, which you can also use directly when a track is made by
another constructor, or when you want to build the stack in a loop. The figure's
width is a setting; its height is simply the sum of its tracks.

[Writing a figure in Rust](../guide/plot.md) and [the Figure
API](../guide/figure.md) cover both.

## Genome, Rings and Panels

Some pictures do not fit on one region of one sequence:

| Container | Use it for |
|:----------|:-----------|
| `Genome` | Several sequences laid end to end, such as the contigs of an assembly, drawn as one figure. |
| `Rings` | A circular chromosome or plasmid, where position becomes an angle. |
| `Panels` | Several drawings on one sheet, labelled A, B, C in the order you add them. |

## Readers

`karyon::read` parses line-based text formats into the values the tracks take.
Every reader takes the text as a string rather than a path, so the library does
no file I/O of its own. That is also why the same code runs at a shell and inside
a web page.

Binary formats such as BAM, CRAM and BCF are not read directly. Convert them
with the tools that already do it well, and pipe the text in:

```bash
samtools depth -a -r NC_000962.3:761000-763000 aln.bam \
  | karyon NC_000962.3:761,000-763,000 --coverage - --label depth -o rpoB.svg
```

[File formats](../guide/formats.md) lists every format, which columns are read,
and which coordinate convention each one uses.

## Coordinates in one paragraph

The region string and the tick labels are **1-based and inclusive**, like
`samtools` and IGV, because those are the numbers a person reads. Everything
else in the API, such as `Feature::new` and `Variant::new`, is **0-based and
half-open**, like BED. So a VCF position goes in as `POS - 1`. The readers do
that subtraction for you, and [Coordinates](../how-it-works/coordinates.md) has
the full story.

## What karyon will not do

Some tracks make a claim rather than just a picture, and when the data cannot
support the claim, karyon refuses with a clear message instead of drawing
something that looks right but is not. For example:

| Situation | What karyon does |
|:----------|:-----------------|
| A tanglegram given only one tree | Refuses: it needs a second tree with `--against`. One tree drawn against itself would show no crossings, which looks like a perfect result. |
| A homology join where no names matched | Refuses and names the first gene it could not find, instead of marking every gene as having no match. |
| Identity values that could be a fraction or a percentage | Asks you to say which with `--identity`, instead of guessing and turning every ribbon into a perfect match. |
| Methylation positions with no reads | Skips and counts them, instead of drawing them as 0 % methylated. |

All four are the same mistake: a value standing in for the absence of a value.

A few plainer limits:

- It does not resample, smooth or interpolate. A track draws the numbers it was
  given, or says it could not.
- It does not read BAM, CRAM or BCF directly; pipe them through `samtools` or
  `bcftools` as shown above.
- It is not on crates.io yet, so it [installs from the
  repository](installation.md).

!!! warning "A common trap: where coverage starts"
    `add_coverage` puts its first value at the left edge of the region. If your
    array starts somewhere else and you change the region, the values silently
    shift and the figure still looks plausible. Use
    `add_coverage_at(start, values)` whenever the array and the region do not
    start at the same position.
