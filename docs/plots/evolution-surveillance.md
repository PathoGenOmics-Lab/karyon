# Evolution and surveillance

Plots for results fitted upstream: ancestral states, selection on branches and
at sites, population trajectories and lineage counts through time. karyon
draws each one with its uncertainty in view; it fits none of the models.
{ .k-lead }

<figure class="k-plate" markdown>
![An eight-panel synthetic atlas: orthogonal and diagonal phylograms with host strips, a curved tree with ancestral-state donuts, mutation symbols and concordance whiskers, a circular tree coloured by branch omega, an unrooted mutation map, a core-versus-accessory tanglegram, protein domains over site-wise selection and observed variants, and an effective population size trajectory above stacked lineage frequencies](../assets/figures/example-evolutionary-surveillance.svg){ width="1408" height="2050" loading="lazy" }
<figcaption>A and B, tree geometries. C, ancestral states and branch events. D, selection on branches. E, the same branch evidence unrooted. F, a tanglegram. G, site-wise selection. H, a phylodynamic trajectory over lineage frequencies.</figcaption>
</figure>

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| What state was each ancestor in, and where did it change? | [Ancestral states and branch events](../tracks/phylogeny.md#treetrack) | `TreeTrack` with `AncestralStateLayer`, `BranchEventLayer` and `BranchIntervalLayer` |
| On which branches is ω above or below one? | [Selection on branches](../tracks/phylogeny.md#treetrack) | `TreeTrack` with `.dnds()`, `BranchRateMixture` and `HomoplasyLayer` |
| Which codons are under selection, and in which direction? | [Site-wise selection](../tracks/variation.md#selectiontrack) | `SelectionTrack` |
| How did the effective population size or R change through time? | [Phylodynamic trajectory](../tracks/evolution-surveillance.md#phylodynamictrack) | `PhylodynamicTrack` |
| Which lineages are rising, and out of how many samples? | [Lineage surveillance](../tracks/evolution-surveillance.md#surveillancetrack) | `SurveillanceTrack` |

All five are built in Rust. The command line draws the tree itself with
`--tree`, but none of these layers or tracks. To compare two trees, see the
[tanglegram](phylogeny-clades.md).

??? info "Why karyon draws these results but fits none of them"
    Each of these results is an estimate that another program made, and each
    plot keeps the estimate apart from its uncertainty and from what was
    observed. A missing value stays missing rather than becoming zero.

    | You supply | karyon draws | karyon does not |
    |:--|:--|:--|
    | Ancestral state probabilities | Donuts on internal nodes, and a mark where a confident state changes | Reconstruct ancestral states |
    | Events on a branch | Ordered symbols on that branch, and dashed curves between branches carrying the same event | Infer events, or claim convergence |
    | An estimate with bounds | A whisker on a branch, or a ribbon through time | Estimate the interval |
    | Branch or site ω | Colours centred on ω = 1, capsules for rate classes, evidence in a tier of its own | Fit a codon model |
    | Effective size, R or growth through time | A line on a linear or log scale | Fit a coalescent model or a clock |
    | Lineage counts and totals | Stacked composition or lines, with alerts that state their reason | Smooth, extrapolate, or fill a missing count with zero |

## Plots

<div class="k-plots" markdown>

-   [![An eight-panel synthetic atlas whose third panel is a curved tree with ancestral-state donuts on internal nodes, mutation symbols on branches and concordance whiskers](../assets/figures/example-evolutionary-surveillance.svg){ width="1408" height="2050" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Ancestral states and branch events](../tracks/phylogeny.md#treetrack)**
    State probabilities as donuts on internal nodes, events on the branch that owns them, and a branch estimate with whiskers (panels C and E).

-   [![One tree drawn four ways with branches coloured by dN/dS on a scale centred at one: rectangular with amino acid changes and host and resistance columns, circular with metadata rings, unrooted, and as a cladogram](../assets/figures/example-phylo-dnds.svg){ width="1508" height="1390" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Selection on branches](../tracks/phylogeny.md#treetrack)**
    Branches coloured by ω, cool below one and warm above it, with the branches that pass a significance cut drawn heavier; fitted rate classes and recurrent changes can go on top.

-   [![A molecular selection atlas: rate classes and recurrent changes on a rectangular tree, mean branch omega on a circular tree, and two site-wise scans over protein domains with evidence above signed omega effects](../assets/figures/example-selection-atlas.svg){ width="1508" height="1053" loading="lazy" }](../tracks/variation.md#selectiontrack)

    **[Site-wise selection](../tracks/variation.md#selectiontrack)**
    Evidence, as a p-value or a posterior, in one tier and the signed log2(ω) effect in another, so a significant purifying site still reads as purifying.

-   [![An eight-panel synthetic atlas whose last panel is an effective population size trajectory on a log scale with its uncertainty ribbon, above stacked lineage frequencies](../assets/figures/example-evolutionary-surveillance.svg){ width="1408" height="2050" loading="lazy" }](../tracks/evolution-surveillance.md#phylodynamictrack)

    **[Phylodynamic trajectory](../tracks/evolution-surveillance.md#phylodynamictrack)**
    An estimate through time as a line with its interval as a ribbon, on a linear or log scale, with a reference such as R = 1 (panel H).

-   [![An eight-panel synthetic atlas whose last panel ends in stacked lineage frequencies by month, with markers where a lineage passed an alert](../assets/figures/example-evolutionary-surveillance.svg){ width="1408" height="2050" loading="lazy" }](../tracks/evolution-surveillance.md#surveillancetrack)

    **[Lineage surveillance](../tracks/evolution-surveillance.md#surveillancetrack)**
    Lineage counts over their totals as stacked composition or as lines, with alerts for frequency and growth that never hide the counts (panel H).

</div>

From a clone of the repository, `cargo run --example evolutionary_surveillance`
writes the atlas to the current directory, so the code behind every panel is
there to read.

## Related

<div class="grid cards" markdown>

-   **[Phylogeny and clades](phylogeny-clades.md)**

    Trees in every layout, two trees face to face, and spans painted onto
    clades.

-   **[Variation and association](variation-association.md)**

    Site-wise selection beside the variants, genes and codons it concerns.

-   **[Phylogenetics](../guide/phylogenetics.md)**

    Reading annotated trees, and the builders behind every layer on this page.

-   **[Recipes](../recipes.md)**

    Complete programs that stack several tracks into one figure.

</div>
