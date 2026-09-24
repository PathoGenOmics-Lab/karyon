# Citation

How to cite karyon and the version you used, who made it, and the methods and
file formats it builds on, which deserve citing in their own right.
{ .k-lead }

## How to cite karyon

karyon does not have a DOI yet, so cite the repository together with the
version you used:

> Ruiz-Rodriguez P, Coscolla M. *karyon: genomic track plots for Rust.*
> PathoGenOmics Lab. <https://github.com/PathoGenOmics-Lab/karyon>

```text
@misc{karyon,
  author       = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  title        = {karyon: genomic track plots for Rust},
  howpublished = {PathoGenOmics Lab},
  url          = {https://github.com/PathoGenOmics-Lab/karyon},
  note         = {Version 0.14.0},
  year         = {2026}
}
```

### Record the version

=== "Rust"

    The library comes from git until it is on crates.io, so `Cargo.lock`
    records the exact commit your figures were drawn with:

    ```toml
    [[package]]
    name = "karyon"
    version = "0.14.0"
    source = "git+https://github.com/PathoGenOmics-Lab/karyon#<commit>"
    ```

    Pin that commit in `Cargo.toml` to draw with it again later:

    ```toml
    [dependencies]
    karyon = { git = "https://github.com/PathoGenOmics-Lab/karyon", rev = "<commit>" }
    ```

=== "Command line"

    ```bash
    karyon --version
    ```

    ```text
    karyon 0.14.0
    ```

The version matters because rendering is deterministic: the same input draws a
byte-identical figure, so a figure can be regenerated exactly, but only with
the version that drew it. A new default or a changed layout is a different
figure from the same data.

!!! tip "Say what the axis counts"
    When the figure is the finding rather than an illustration, say in the
    caption what its horizontal axis counts. A stack of tracks counts positions
    on a sequence, a variable-site panel counts sites, an alignment counts
    columns, a domain panel counts residues and a squiggle counts samples, and
    "drawn with karyon" does not say which. See
    [Coordinates](../how-it-works/coordinates.md#the-region-is-a-coordinate-system).

## Authors

**Paula Ruiz-Rodriguez** and **Mireia Coscolla**

I²SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection
and Public Health, Valencia, Spain.

## Methods it builds on

Most of what karyon draws is a standard representation, and the ones that are
not carry somebody else's idea. When a figure leans on one of these, cite it
as well.

| What karyon draws | Whose idea | Reference |
|:--|:--|:--|
| Sequence logos scaled by information content, `LogoScore::InformationContent` | Schneider and Stephens; the axis fixed from zero to log2 of the alphabet size and the tallest symbol on top follow WebLogo | [1], [2] |
| Enrichment and depletion logos, `LogoTrack::edlogo`, the other `LogoScore` schemes, and the shrinkage behind `LogoTrack::stabilize` and the `dash` module | Logolas | [3] |
| Only the columns that vary, `SnpTrack` | the idea snipit is built around; the implementation and the drawing are karyon's own | [4] |
| Model attribution drawn as the bases themselves, `DynseqTrack` | the convention of BPNet and the dynseq browser track | [5], [6] |
| The default base colours, `BaseColors::conventional` | the IGV-style convention readers of genome figures already know: A green, C blue, G orange, T red | [7] |
| The land under a `Map` | Natural Earth's 1:110m land polygons, which are in the public domain | [Natural Earth](https://www.naturalearthdata.com/) |

The genome-wide significance line that `ManhattanTrack::genome_wide_threshold`
and `--threshold genome-wide` draw is `-log10(5e-8)`, a Bonferroni correction
for a million independent tests. It is a convention from human association
studies rather than a property of any one study, and it is often the wrong
number elsewhere: cite whatever fixed the threshold you actually used.

## Formats it reads

The readers open no files: `karyon::read` parses the line-based text a
genomics shell already writes, and the caller decides where the text came
from. The formats themselves are defined elsewhere. How each one's coordinates
are converted is on
[Coordinates](../how-it-works/coordinates.md#what-a-files-numbers-become).

| Format | Read by | Reference |
|:--|:--|:--|
| SAM | `--pileup`, `--split-reads` | [8] |
| VCF | `--variants`, `--structural` | [9] |
| BED, bedGraph, cytoBand | `--features`, `--loci`, `--coverage`, `--windows`, `--dynseq`, `--ideogram` | [10] |
| GFF3 | `--features`, `--loci`, `--clades` | [the GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) |
| PAF, from minimap2 | `--synteny`, `--dotplot` | [11] |
| bedMethyl, from modkit | `--methylation` | [modkit](https://github.com/nanoporetech/modkit) |
| Bismark methylation extractor output | `--bisulfite` | [12] |
| STAR `SJ.out.tab` | `--junctions` | [13] |
| InterProScan TSV | `--domains` | [14] |
| Gubbins recombination GFF | `--clades` | [15] |
| CNVkit `.cns` | `--copy-number` | [16] |
| ASCAT segments | `--copy-number` | [17] |
| NEXUS | `Tree::parse_nexus` | [18] |
| BLAST tabular output | `--links` | [19] |

## References

1. Schneider TD, Stephens RM. Sequence logos: a new way to display consensus
   sequences. *Nucleic Acids Research*. 1990;18(20):6097-6100.
   [doi:10.1093/nar/18.20.6097](https://doi.org/10.1093/nar/18.20.6097)
2. Crooks GE, Hon G, Chandonia JM, Brenner SE. WebLogo: a sequence logo
   generator. *Genome Research*. 2004;14(6):1188-1190.
   [doi:10.1101/gr.849004](https://doi.org/10.1101/gr.849004)
3. Dey KK, Xie D, Stephens M. A new sequence logo plot to highlight enrichment
   and depletion. *BMC Bioinformatics*. 2018;19:473.
   [doi:10.1186/s12859-018-2489-3](https://doi.org/10.1186/s12859-018-2489-3)
4. O'Toole Á, Aziz A, Maloney D. Publication-ready single nucleotide
   polymorphism visualization with snipit. *Bioinformatics*.
   2024;40(8):btae510.
   [doi:10.1093/bioinformatics/btae510](https://doi.org/10.1093/bioinformatics/btae510)
5. Avsec Ž, Weilert M, Shrikumar A, Krueger S, Alexandari A, Dalal K, et al.
   Base-resolution models of transcription-factor binding reveal soft motif
   syntax. *Nature Genetics*. 2021;53(3):354-366.
   [doi:10.1038/s41588-021-00782-6](https://doi.org/10.1038/s41588-021-00782-6)
6. Nair S, Barrett A, Li D, Raney BJ, Lee BT, Kerpedjiev P, et al. The dynseq
   browser track shows context-specific features at nucleotide resolution.
   *Nature Genetics*. 2022;54(11):1581-1583.
   [doi:10.1038/s41588-022-01194-w](https://doi.org/10.1038/s41588-022-01194-w)
7. Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G,
   et al. Integrative genomics viewer. *Nature Biotechnology*.
   2011;29(1):24-26. [doi:10.1038/nbt.1754](https://doi.org/10.1038/nbt.1754)
8. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The
   Sequence Alignment/Map format and SAMtools. *Bioinformatics*.
   2009;25(16):2078-2079.
   [doi:10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)
9. Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. The
   variant call format and VCFtools. *Bioinformatics*. 2011;27(15):2156-2158.
   [doi:10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330)
10. Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, et al. The
    human genome browser at UCSC. *Genome Research*. 2002;12(6):996-1006.
    [doi:10.1101/gr.229102](https://doi.org/10.1101/gr.229102)
11. Li H. Minimap2: pairwise alignment for nucleotide sequences.
    *Bioinformatics*. 2018;34(18):3094-3100.
    [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
12. Krueger F, Andrews SR. Bismark: a flexible aligner and methylation caller
    for Bisulfite-Seq applications. *Bioinformatics*. 2011;27(11):1571-1572.
    [doi:10.1093/bioinformatics/btr167](https://doi.org/10.1093/bioinformatics/btr167)
13. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. STAR:
    ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15-21.
    [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
14. Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, et al.
    InterProScan 5: genome-scale protein function classification.
    *Bioinformatics*. 2014;30(9):1236-1240.
    [doi:10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031)
15. Croucher NJ, Page AJ, Connor TR, Delaney AJ, Keane JA, Bentley SD, et al.
    Rapid phylogenetic analysis of large samples of recombinant bacterial whole
    genome sequences using Gubbins. *Nucleic Acids Research*. 2015;43(3):e15.
    [doi:10.1093/nar/gku1196](https://doi.org/10.1093/nar/gku1196)
16. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: genome-wide copy number
    detection and visualization from targeted DNA sequencing. *PLoS
    Computational Biology*. 2016;12(4):e1004873.
    [doi:10.1371/journal.pcbi.1004873](https://doi.org/10.1371/journal.pcbi.1004873)
17. Van Loo P, Nordgard SH, Lingjærde OC, Russnes HG, Rye IH, Sun W, et al.
    Allele-specific copy number analysis of tumors. *Proceedings of the
    National Academy of Sciences of the USA*. 2010;107(39):16910-16915.
    [doi:10.1073/pnas.1009843107](https://doi.org/10.1073/pnas.1009843107)
18. Maddison DR, Swofford DL, Maddison WP. NEXUS: an extensible file format for
    systematic information. *Systematic Biology*. 1997;46(4):590-621.
    [doi:10.1093/sysbio/46.4.590](https://doi.org/10.1093/sysbio/46.4.590)
19. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al.
    BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421.
    [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

## Licence

karyon is released under the
[MIT licence](https://github.com/PathoGenOmics-Lab/karyon/blob/main/LICENSE),
copyright 2026 Paula Ruiz-Rodriguez and Mireia Coscolla. A plotting library is
meant to be a dependency, and a permissive licence lets any tool depend on it
whatever that tool's own licence is.

## Where next

<div class="grid cards" markdown>

-   **[Changelog](changelog.md)**

    What moved between the version you cited and this one.

-   **[Contributing](contributing.md)**

    How to report a figure that came out wrong, and what a change has to pass.

</div>
