---
title: A nanopore signal
description: Draw the raw current of a nanopore read, from a SLOW5 file.
---

# A nanopore signal

You have the raw current of nanopore reads, as a SLOW5 file.

```bash
karyon reads.slow5 -o signal.svg
```

<figure class="k-start" markdown>
![The current of one read over 2,400 samples, stepping between levels from about 60 to 125 picoamperes](../assets/start/signal.svg){ .k-light width="720" height="156" }
![The same figure on the dark page](../assets/start/signal-dark.svg){ .k-dark width="720" height="156" }
</figure>

The ruler counts samples, and the current is in picoamperes, put there with the
read's own digitisation, offset and range. The first read of the file is drawn
and named in the margin, and the command says how many others the file holds.

## Your file

SLOW5 is text. A BLOW5 file is the same in binary, and `slow5tools view`
writes it as text. A file of plain numbers, one sample after another, is read
too, as picoamperes already.

The bases the basecaller called come from its SAM or BAM, as Dorado writes it
with `--emit-moves`: its move table puts each base over the stretch of
current it was called from.

## Change it

| To | Write |
|:--|:--|
| Another read | `--read read_2` |
| The bases the basecaller called, over the current | `--with-moves moves.sam` |
| Part of the read | `sample:1-500` as the place, before the file |
| Another colour | `--color '#d55e00'` |

The example files: [reads.slow5](../data/reads.slow5), and
[moves.sam](../data/moves.sam) for the bases. Every option:
`karyon help squiggle`, or the [command line reference](../guide/cli.md).
