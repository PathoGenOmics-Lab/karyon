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

## Change it

| To | Write |
|:--|:--|
| Another read | `--read read_2` |
| Part of the read | `sample:1-500` as the place, before the file |
| Another colour | `--color '#d55e00'` |

The example file: [reads.slow5](../data/reads.slow5). Every option:
`karyon help squiggle`, or the [command line reference](../guide/cli.md).
