---
title: Two assemblies
description: Draw how two assemblies of a genome line up, from a PAF alignment, as ribbons or as a dot plot.
---

# Two assemblies

You have an alignment of one assembly against another, in PAF, as `minimap2`
writes it.

```bash
karyon asm1_chr1 assemblies.paf -o assemblies.svg
```

<figure class="k-start" markdown>
![Two assemblies of one chromosome joined by ribbons: straight ribbons where they agree, a crossed pink pair where a stretch is turned round, and crossing ribbons where a stretch has moved](../assets/start/assemblies.svg){ .k-light width="720" height="234" }
![The same figure on the dark page](../assets/start/assemblies-dark.svg){ .k-dark width="720" height="234" }
</figure>

The place is the sequence of the first assembly, as the PAF names it. The
ribbons join each stretch to where it aligns in the second assembly, and the
key says which colour is the same strand and which is reversed.

## Change it

| To | Write |
|:--|:--|
| Draw it as a dot plot | `--dotplot assemblies.paf` |
| Zoom into one stretch | `asm1_chr1:100,000-300,000` in place of `asm1_chr1` |
| Make the PAF | `minimap2 -x asm5 first.fa second.fa > assemblies.paf` |

The example file: [assemblies.paf](../data/assemblies.paf). Every option:
`karyon help synteny`, or the [command line reference](../guide/cli.md).
