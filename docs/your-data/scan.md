---
title: An association scan
description: Draw a genome-wide association table as a Manhattan plot, with its significance line.
---

# An association scan

You have a table from an association tool: PLINK, PLINK 2, REGENIE, SAIGE,
GEMMA, BOLT-LMM or the GWAS Catalog.

```bash
karyon 1 gwas.assoc --threshold genome-wide -o scan.svg
```

<figure class="k-start" markdown>
![A scan along chromosome 1, flat except for one peak whose markers rise above the dashed line at p = 5e-8](../assets/start/scan.svg){ .k-light width="720" height="185" }
![The same figure on the dark page](../assets/start/scan-dark.svg){ .k-dark width="720" height="185" }
</figure>

The place, `1`, is the chromosome as the table names it. Each marker is as
high as `-log10` of its p-value, and the dashed line is p = 5e-8.

## Change it

| To | Write |
|:--|:--|
| Zoom into the peak | `1:700,000-820,000` in place of `1` |
| Draw the line elsewhere | `--threshold 1e-5` |
| Add the genes under the peak | `NC_000962.3:700,000-820,000 gwas.assoc --threshold genome-wide genes.gff3 --rename 1=NC_000962.3` |

`--rename` is for a table that names the chromosome `1` where the annotation
calls it `NC_000962.3`.

## Colour a peak by linkage

Zoomed into the peak, each marker coloured by its linkage with the strongest,
as LocusZoom draws one, with the recombination rate under it:

```bash
karyon 1:661,000-861,000 gwas.assoc --ld lead.ld --threshold genome-wide \
  --coverage recombination.bedgraph --style line --label cM/Mb -o locus.svg
```

<figure class="k-start" markdown>
![The markers of the peak coloured from grey to blue by their r-squared with the strongest, which is a diamond labelled with its position, and under them a recombination rate with two hotspots](../assets/start/locus.svg){ .k-light width="720" height="286" }
![The same figure on the dark page](../assets/start/locus-dark.svg){ .k-dark width="720" height="286" }
</figure>

`lead.ld` is PLINK's linkage of the lead with its neighbours, as
`plink --r2 --ld-snp snp00342 --ld-window-kb 100 --ld-window 99999
--ld-window-r2 0` writes it. The recombination rate is any bedGraph.

The example files: [gwas.assoc](../data/gwas.assoc), [lead.ld](../data/lead.ld)
and [recombination.bedgraph](../data/recombination.bedgraph). Every option:
`karyon help manhattan`, or the [command line reference](../guide/cli.md).
