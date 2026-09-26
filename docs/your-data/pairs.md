---
title: Pairs of positions
description: Draw linkage between variants, contacts between the bins of a genome, or scores between sites, as a triangle or as arcs.
---

# Pairs of positions

You have a value for pairs of positions: linkage between variants from PLINK,
contacts between the bins of a Hi-C map, or scores between sites.

```bash
karyon rpoB genes.gff3 linkage.ld -o pairs.svg
```

<figure class="k-start" markdown>
![A gene with thirty-six variants under it, and under them a triangle in which each pair of variants is a cell coloured by its linkage: three dark triangles where variants are inherited together](../assets/start/pairs.svg){ .k-light width="720" height="386" }
![The same figure on the dark page](../assets/start/pairs-dark.svg){ .k-dark width="720" height="386" }
</figure>

Each pair of variants is a cell under the point half way between them, as deep
as they are far apart, and coloured by its r² from 0 to 1. Variants inherited
together are the dark triangles.

## Your file

- **Linkage**: PLINK's `--r2` writes a `.ld` table, read as it is.
  `--ld-window-r2 0` keeps the weak pairs as well, which a complete triangle
  needs.
- **Contacts or loops**: BEDPE, as `cooler dump --join` writes a contact map.
- **Anything else**: a table headed `pos1`, `pos2` and a value, its positions
  counted from 1.

## Change it

| To | Write |
|:--|:--|
| Arcs between a few pairs far apart | nothing: a few pairs are drawn as arcs by themselves, and `--style arcs` or `--style triangle` chooses |
| Only the strong pairs | `--threshold 0.5` |
| Contacts, which fall by orders of magnitude | `--log` |
| Scores between the calls, over the calls | `rpoB genes.gff3 calls.vcf.gz --pairs epistasis.tsv` |

The example files: [linkage.ld](../data/linkage.ld) and
[epistasis.tsv](../data/epistasis.tsv). Every option: `karyon help pairs`, or
the [command line reference](../guide/cli.md).
