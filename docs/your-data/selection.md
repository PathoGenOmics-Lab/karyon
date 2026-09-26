---
title: Selection along a gene
description: Draw a test of selection at each site of a gene, from HyPhy or a table of your own.
---

# Selection along a gene

You have a test of selection at each codon of a gene: HyPhy's FEL or MEME, a
site model from PAML, or a table of your own.

```bash
karyon --selection fel.csv -o selection.svg
```

<figure class="k-start" markdown>
![Three hundred sites of a gene: above, the evidence at each, with two stretches of sites rising past p = 0.05 as diamonds; below, each site's omega, most below 1 and those same sites well above it](../assets/start/selection.svg){ .k-light width="720" height="216" }
![The same figure on the dark page](../assets/start/selection-dark.svg){ .k-dark width="720" height="216" }
</figure>

Above, the evidence at each site, as `-log10` of its p-value, and the sites at
p ≤ 0.05 as diamonds. Below, its direction: ω, which is dN/dS, above 1 upwards
and below 1 downwards. The ruler counts sites from 1, as HyPhy does.

## Your table

HyPhy's CSV is read as it is: `alpha`, `beta` and `p-value`, one row per site
in order. Any other table names its columns in a header: a `site` or `codon`
where the rows are not in order from 1, the rates as `alpha` and `beta`, `dS`
and `dN`, or their ratio as `omega`, and the evidence as a `p-value` or a
`posterior`. MEME's `beta-`, `beta+` and `p+` are drawn as its two rate
classes.

## Change it

| To | Write |
|:--|:--|
| A stricter line | `--threshold 0.01` |
| A table of posteriors, from FUBAR or a Bayes empirical Bayes | nothing: the header says so, and `--threshold 0.95` moves the line |
| Some of the sites | `site:50-200` as the place, before the file |

The example file: [fel.csv](../data/fel.csv). Every option:
`karyon help selection`, or the [command line reference](../guide/cli.md).
