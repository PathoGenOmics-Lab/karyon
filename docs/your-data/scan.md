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

The example file: [gwas.assoc](../data/gwas.assoc). Every option:
`karyon help manhattan`, or the [command line reference](../guide/cli.md).
