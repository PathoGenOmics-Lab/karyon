---
title: A tree and its samples
description: Draw a Newick tree with strips of what you know about each sample beside its tips.
---

# A tree and its samples

You have a tree (Newick) and a sample sheet: a table with the sample names in
its first column and one column for each thing you know about them.

```bash
karyon tree.nwk --traits samples.tsv --columns lineage,country -o tree.svg
```

<figure class="k-start" markdown>
![A tree of forty samples with a strip of lineage colours and a column of country symbols beside the tips, and a key underneath](../assets/start/tree.svg){ .k-light width="720" height="717" }
![The same figure on the dark page](../assets/start/tree-dark.svg){ .k-dark width="720" height="717" }
</figure>

Each column you name is drawn beside the tips, and the key underneath says
which colour is which. A column with more values than colours uses shapes as
well. A tree needs no place.

## Change it

| To | Write |
|:--|:--|
| Colour the branches too | `--color-by lineage` |
| Draw it as a circle | `--projection circular` |
| Show the support values of 70 and over | `--support-style labels --threshold 70` |
| Fold a large tree to fit | `--max-rows 60` |
| Draw one clade | `--focus S26,S27` |

The example files: [tree.nwk](../data/tree.nwk) and
[samples.tsv](../data/samples.tsv). Every option: `karyon help tree`, or the
[command line reference](../guide/cli.md).
