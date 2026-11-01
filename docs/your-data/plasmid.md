---
title: A plasmid or circular genome
description: Draw a plasmid, an organelle genome or another sequence with no ends as a circle, its genes as arrows on a ring, from its annotation alone.
---

# A plasmid or circular genome

You have a plasmid, or another sequence with no ends, and its annotation as
GFF3, as an annotation tool such as Bakta or Prokka writes it.

```bash
karyon pEX1 --circular pEX1.gff3 -o plasmid.svg
```

<figure class="k-start" markdown>
![The plasmid drawn as a circle, its name and length in the middle and a ruler round the outside from 0 at the top: five genes as arcs ending in a point, repA, tetA and mobA on the outer half pointing clockwise and bla and tetR on the inner half pointing back, tetR and tetA pointing away from each other, and two short arcs with no point in a muted ink for the origins of replication and of transfer; the key underneath names the ring genes and what its colours mean](../assets/start/plasmid.svg){ .k-light width="720" height="751" }
![The same figure on the dark page](../assets/start/plasmid-dark.svg){ .k-dark width="720" height="751" }
</figure>

The place is the sequence, by its name, and `--circular` draws it round:
position zero at twelve o'clock, running clockwise, inside a ruler. A circle
closes where its sequence ends, which the annotation's `##sequence-region`
line says here; a FASTA of the sequence among the files says it too. Each gene
is an arc ending in a point at its 3' end, on the outer half of the ring for
the forward strand and the inner half for the reverse, so the divergent pair
tetR and tetA point away from each other. What the annotation calls something
other than a gene, here the origins of replication and of transfer, is an arc
with no point in a muted ink, keyed by what the annotation calls it. The ring
is named for what it holds, and the key under the circle says what each colour
means.

## Change it

| To | Write |
|:--|:--|
| Add the sequence's GC skew as a second ring, which turns over at the origin of replication | `pEX1.fa` after `pEX1.gff3` |
| One colour for every gene | `pEX1.gff3 --color '#6a3d9a'` |
| Name the ring yourself | `pEX1.gff3 --label annotation` |
| Leave the names off | `pEX1.gff3 --no-names` |
| Leave out the key | `--no-legend` |
| Another title in the middle | `--title 'pEX1, a plasmid'` |
| Draw it along a line, the genes as arrows | leave out `--circular` |
| Draw one stretch of it along a line | `pEX1:3,000-5,500 pEX1.gff3`, with no `--circular` |

Calls and depth over the plasmid are rings of their own in the same way, in
the order written, the first outermost; [a whole sequence](genome.md#as-a-circle)
draws a chromosome round with both.

The example files: [pEX1.gff3](../data/pEX1.gff3) and
[pEX1.fa](../data/pEX1.fa). Every option: `karyon help features`, or the
[command line reference](../guide/cli.md#a-whole-sequence-as-a-circle).
