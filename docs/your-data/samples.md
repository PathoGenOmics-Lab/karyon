---
title: Many samples in windows
description: Draw the depth, copy number or methylation of many samples in windows along a genome, in the order of a tree.
---

# Many samples in windows

You have a value for many samples in windows along a genome: the depth of each
sample, a copy number, a methylation level.

```bash
karyon NC_000962.3 --heatmap depths.tsv --relative \
  --with-tree tree.nwk --label depth -o samples.svg
```

<figure class="k-start" markdown>
![Forty samples in the order of a tree, each a row of cells along the chromosome, pale where the depth is usual: one clade has a blue cell where it lost a stretch, and three samples a pink one where they carry a stretch twice](../assets/start/heatmap.svg){ .k-light width="720" height="602" }
![The same figure on the dark page](../assets/start/heatmap-dark.svg){ .k-dark width="720" height="602" }
</figure>

One row per sample and one cell per window, in the order of the tree's tips.
`--relative` reads each sample against its own median, so 1× is its usual
depth, drawn pale: a sample sequenced deeper is not a darker row, and what
changed along the genome is what stands out, a loss in blue and a gain in
pink. Here one clade lost a stretch, and three samples carry another twice. A table of windows does not say how long the
sequence is, so karyon draws as far as its last window and says so; write the
span, as `NC_000962.3:1-4,411,532`, to set it yourself.

## Your table

As `bedtools unionbedg` writes it: a sequence, a start and an end, counted from
0 as BED is, then one column per sample, named in the header.

```bash
bedtools unionbedg -header -names S01 S02 S03 -i S01.bg S02.bg S03.bg > depths.tsv
```

deepTools' `multiBigwigSummary --outRawCounts` writes the same shape, and is
read too, as is the long form, a sequence, a start, an end, a sample and its
value to a row.

## Change it

| To | Write |
|:--|:--|
| The values as they are | leave out `--relative` |
| Log ratios, read either side of nought | `--center 0` in place of `--relative` |
| What is known about each sample, beside the rows | `--traits samples.tsv` |
| Thinner rows | `--row-height 6` |
| One stretch of the genome | `NC_000962.3:1,300,000-1,600,000` as the place |

The example file: [depths.tsv](../data/depths.tsv), with
[tree.nwk](../data/tree.nwk). Every option: `karyon help heatmap`, or the
[command line reference](../guide/cli.md).
