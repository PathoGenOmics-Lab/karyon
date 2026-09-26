---
title: An alignment and its tree
description: Draw a multiple sequence alignment with its rows in the order of a tree, and the tree beside them.
---

# An alignment and its tree

You have an alignment (FASTA with every record the same length) and, if you
like, a tree of the same samples (Newick).

```bash
karyon --msa aln.fasta --with-tree tree.nwk -o alignment.svg
```

<figure class="k-start" markdown>
![Forty rows of a 300-column alignment in the order of the tree drawn beside them, with a coloured mark wherever a sample differs from the consensus, and a key to the colours of the four bases](../assets/start/alignment.svg){ .k-light width="720" height="642" }
![The same figure on the dark page](../assets/start/alignment-dark.svg){ .k-dark width="720" height="642" }
</figure>

Each row is a sample, in the order of the tree's tips. A cell is coloured
where the sample differs from the consensus, and the key names the colours.
An alignment needs no place: it is drawn over all its columns.

## Change it

| To | Write |
|:--|:--|
| Only the columns that vary | `--snps aln.fasta --with-tree tree.nwk` |
| Every base, not only the differences | `--style all` after the alignment |
| Compare with one sample | `--compare-to S01` |
| Add what you know about each sample | `--traits samples.tsv --columns lineage` |
| Zoom into columns 100 to 160 | `aln.fasta:100-160 --msa aln.fasta --with-tree tree.nwk` |

The example files: [aln.fasta](../data/aln.fasta), [tree.nwk](../data/tree.nwk)
and [samples.tsv](../data/samples.tsv). Every option: `karyon help msa`, or
the [command line reference](../guide/cli.md).
In Rust, it is [a few lines](../guide/phylogenetics.md#put-an-alignment-in-the-order-of-the-tree).
