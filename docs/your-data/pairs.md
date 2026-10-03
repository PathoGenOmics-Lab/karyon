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
- **Contacts**: a `.hic`, as Juicer's tools and hictk write it, named as it
  is: `karyon chr1:20,000,001-22,000,000 contacts.hic`. It is read at the
  finest of its resolutions that cuts the window into 250 bins or fewer, which
  a note names, and `--resolution 10000` picks another it holds. Its raw counts
  are drawn, with none of its normalisations applied. A `.hic` lists only the
  cells that hold a count, and a cell it does not list is left as the page,
  not drawn in the pale end of the key; the triangle reaches as far from the
  diagonal as the farthest cell it lists in the window.
- **Loops, or a cooler file**: BEDPE, as loop callers write their loops and
  `cooler dump --join` writes a contact map. A `.cool` is drawn with
  `--pairs <(cooler dump --join -r REGION map.cool)`, and a `.mcool` named on
  its own is answered with the steps that write it so.
- **Anything else**: a table headed `pos1`, `pos2` and a value, its positions
  counted from 1.

## Change it

| To | Write |
|:--|:--|
| Arcs between a few pairs far apart | nothing: a few pairs are drawn as arcs by themselves, and `--style arcs` or `--style triangle` chooses |
| Only the strong pairs | `--threshold 0.5` |
| Weak linkage read on a scale of its own, rather than against an r² of 1 | `--max 0.5` |
| Contacts, which fall by orders of magnitude | `--log` |
| A contact map at bins of 10 kb, where the file holds them | `--resolution 10000` |
| The ruler along the diagonal, over the triangle | `--axis` before the file |
| Scores between the calls, over the calls | `rpoB genes.gff3 calls.vcf.gz --pairs epistasis.tsv` |

The example files: [linkage.ld](../data/linkage.ld) and
[epistasis.tsv](../data/epistasis.tsv). Every option: `karyon help pairs`, or
the [command line reference](../guide/cli.md).
