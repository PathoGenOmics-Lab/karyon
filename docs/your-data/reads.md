---
title: Reads, calls and genes
description: Draw the depth of a BAM, the calls of a VCF and the genes of a GFF3 over one place.
---

# Reads, calls and genes

You have aligned reads (a BAM with its `.bai`), variant calls (a VCF, gzipped
or not) and an annotation (GFF3 or GTF).

```bash
karyon rpoB reads.bam genes.gff3 calls.vcf.gz -o rpoB.svg
```

<figure class="k-start" markdown>
![The depth of the reads over the gene rpoB, the gene as an arrow, and seven calls as lollipops](../assets/start/reads.svg){ .k-light width="720" height="248" }
![The same figure on the dark page](../assets/start/reads-dark.svg){ .k-dark width="720" height="248" }
</figure>

- **reads depth**: how many reads cover each base. The dip is a stretch no
  read covers.
- **genes**: each gene once, with its name and its direction.
- **calls**: each variant at its position, as tall as its allele frequency,
  coloured by what the annotation says it does.

## Change it

| To | Write |
|:--|:--|
| See the reads, not their depth | `--pileup reads.bam` |
| See the bases, zoomed in | `NC_000962.3:761,100-761,200 --pileup reads.bam ref.fa` |
| Draw a wider stretch | `NC_000962.3:755,000-770,000` in place of `rpoB` |
| Make a row taller | `reads.bam --height 100` |
| Rename a row | `calls.vcf.gz --label "variant calls"` |

The example files:
[reads.bam](../data/reads.bam), [reads.bam.bai](../data/reads.bam.bai),
[genes.gff3](../data/genes.gff3), [calls.vcf.gz](../data/calls.vcf.gz) and
[ref.fa](../data/ref.fa). Every option: `karyon help coverage`, or the
[command line reference](../guide/cli.md).
