# Phylogeny and clades

Plots that put a tree beside the evidence: annotated trees in several layouts,
two trees face to face, genomic spans painted onto clades, and a tree around a
map. karyon draws the tree you give it; it does not infer one.
{ .k-lead }

## How to choose

| Your question | Plot | Build it with |
|:--|:--|:--|
| What does the tree say, with its support and metadata beside it? | [Annotated tree](../tracks/phylogeny.md#treetrack) | `TreeTrack`, `--tree` |
| Where do two trees over the same taxa disagree? | [Tanglegram](../tracks/phylogeny.md#tanglegramtrack) | `TanglegramTrack`, `--tanglegram` with `--against` |
| Which clades carry a genomic span, and is it one event? | [Clade blocks](../tracks/phylogeny.md#cladetrack) | `CladeTrack`, `--clades` with `--with-tree` |
| Where were the tips sampled? | [Tree around a map](../guide/maps.md) | `PhyloMap`, Rust only |

To sort rows of sites, genes, residues or domains by a tree, hand the tree to
`SnpTrack`, `MatrixTrack`, `MsaTrack` or `DomainTrack` with `.tree()`: their
plots are under [variation](variation-association.md) and
[comparisons](comparisons-alignments.md).

## Plots

<div class="k-plots" markdown>

-   [![A dated outbreak tree with branches coloured by country and aligned country and depth columns, beside the same tree with two named clades collapsed into triangles](../assets/figures/example-phylogenetics.svg){ width="1540" height="368" loading="lazy" .k-wide }](../tracks/phylogeny.md#treetrack)

    **[Annotated tree](../tracks/phylogeny.md#treetrack)**
    A phylogram, cladogram or time tree with metadata aligned to its tips, branches coloured by a trait and clades collapsed without changing the tree.

-   [![Four radial views of one outbreak tree: a circular time tree with country and depth rings, a 250-degree fan with a collapsed clade, time radiating inwards, and a circular cladogram](../assets/figures/example-phylo-layouts.svg){ width="1400" height="1228" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Circular and radial trees](../tracks/phylogeny.md#treetrack)**
    The same tree as a full circle, a partial fan or radiating inwards, with metadata as rings; the topology and the order of the tips stay the same.

-   [![An unrooted tree and a circular cladogram carrying the same four datasets: a country strip, radial depth bars, binary resistance marks and host symbols](../assets/figures/example-phylo-annotations.svg){ width="1480" height="712" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Unrooted trees and datasets](../tracks/phylogeny.md#treetrack)**
    An unrooted layout centred on the topology rather than on the root the file happened to use, with strips, bars, binary marks and symbols as datasets.

-   [![One phylogram drawn rectangular, circular and unrooted, with support values as symbols and labels, mutations written along their branches and a branch-length scale bar](../assets/figures/example-phylo-evidence.svg){ width="1739" height="630" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Support, events and scale bars](../tracks/phylogeny.md#treetrack)**
    Support as symbols or labels above a threshold, each branch's own events written along it, and a scale bar in branch-length units.

-   [![The same phylogram rooted three ways, at the source root, on a checked monophyletic outgroup and at the weighted midpoint, each root marked with a diamond](../assets/figures/example-phylo-reroot.svg){ width="1739" height="360" loading="lazy" .k-wide }](../tracks/phylogeny.md#treetrack)

    **[Rooting choices](../tracks/phylogeny.md#treetrack)**
    Reroot at a node, on an outgroup checked for monophyly or at the midpoint, with the new root marked.

-   [![Abundance bubbles and host bars on a rectangular tree with a highlighted transmission cluster, ancestral-host donuts and a highlighted sector on a radial tree, and an alignment and domain architectures sorted by the same tree](../assets/figures/example-phylo-faces.svg){ width="1386" height="660" loading="lazy" }](../tracks/phylogeny.md#treetrack)

    **[Node glyphs and clade highlights](../tracks/phylogeny.md#treetrack)**
    Bubbles, pies, donuts or stacked bars on nodes and a band or sector behind a clade; a missing value draws no glyph rather than a zero.

-   [![Core and accessory genome trees of eight isolates face to face, matching tips joined across the middle, with a header reporting seven crossings reduced to five](../assets/figures/example-tanglegram.svg){ width="760" height="234" loading="lazy" }](../tracks/phylogeny.md#tanglegramtrack)

    **[Tanglegram](../tracks/phylogeny.md#tanglegramtrack)**
    Two trees facing each other with every shared tip joined; untangling rotates clades to remove crossings but never changes either tree.

-   [![Lineage-defining deletions painted across a SARS-CoV-2 phylogeny, each block spanning the lineages that carry it, with the recurrent one cut out where a lineage between the carriers lacks it](../assets/figures/example-clades.svg){ width="880" height="241" loading="lazy" .k-wide }](../tracks/phylogeny.md#cladetrack)

    **[Clade blocks](../tracks/phylogeny.md#cladetrack)**
    A genomic span painted across the clade that carries it, with any row inside the span that does not carry it cut out, so a scattered set never passes for a clade.

-   [![Two circular phylogenies around an orthographic globe: a calendar time tree with one link per reported location, and a partial cladogram with a link from every sample](../assets/figures/example-phylo-map.svg){ width="1916" height="868" loading="lazy" }](../guide/maps.md)

    **[Tree around a map](../guide/maps.md)**
    A circular tree framing a map, its tips linked to the locations you supply, place by place or sample by sample; nothing about migration is inferred.

</div>

## Related

<div class="grid cards" markdown>

-   **[Evolution and surveillance](evolution-surveillance.md)**

    Ancestral states and selection on the same trees, and trajectories and
    lineage counts through time.

-   **[Variation and association](variation-association.md)**

    Variable sites and genotype matrices sorted by a tree.

-   **[Phylogenetics](../guide/phylogenetics.md)**

    Reading trees and their annotations, dates, layouts, rooting and
    collapsing.

-   **[Tree viewer](../tree.md)**

    Open a tree in the browser and explore it, from a handful of tips to a
    million.

</div>
