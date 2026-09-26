---
title: A whole sequence
description: Draw the depth of reads along a whole chromosome in windows, for one sample or several, to see what stands out at a glance.
---

# A whole sequence

You have the depth of your reads in windows along a whole chromosome, as
bedGraph, for one sample or several.

```bash
karyon NC_000962.3 sampleA.bedgraph sampleB.bedgraph -o genome.svg
```

<figure class="k-start" markdown>
![The depth of two samples along a whole chromosome in 10 kb windows: the first drops to nothing over one stretch and doubles over another, the second is level throughout](../assets/start/genome.svg){ .k-light width="720" height="227" }
![The same figure on the dark page](../assets/start/genome-dark.svg){ .k-dark width="720" height="227" }
</figure>

The place is the whole sequence, by its name, and each file is a row. The first
sample drops to nothing where it has lost a stretch and doubles where it
carries one twice. A bedGraph does not say how long the sequence is, so karyon
draws as far as the files reach and says so; write the span, as
`NC_000962.3:1-4,411,532`, to set it yourself.

## Change it

| To | Write |
|:--|:--|
| Zoom into the stretch with no reads | `NC_000962.3:1,400,000-1,560,000` in place of `NC_000962.3` |
| A log scale for the depth | `sampleA.bedgraph --log` |
| Make depth windows from a BAM | `mosdepth --by 10000 sample reads.bam`, then draw `sample.regions.bed.gz` |

The example files: [sampleA.bedgraph](../data/sampleA.bedgraph) and
[sampleB.bedgraph](../data/sampleB.bedgraph). Every option:
`karyon help coverage`, or the [command line reference](../guide/cli.md).
