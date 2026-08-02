# Citation

karyon is not archived anywhere yet. There is no Zenodo DOI and no
`CITATION.cff` in the repository, so until there is, cite the repository and
the version you used:

> Ruiz-Rodriguez P, Coscolla M. *karyon: genomic track plots for Rust.*
> PathoGenOmics Lab. <https://github.com/PathoGenOmics-Lab/karyon>

Record the version as well. `karyon --version` prints it for the command, and
the `Cargo.lock` entry carries it for the library, together with the exact
commit when the dependency comes from git, which it does until karyon is on
crates.io. The version is worth recording because rendering is deterministic:
the same input produces byte-identical output, so a figure can be regenerated
exactly, but only against the version that drew it. A new default or a changed
layout is a different figure from the same data.

If the figure is the finding rather than the illustration, say what the axis
is. A stack of tracks is in genomic coordinates, a variable-site panel counts
sites, an alignment counts columns and a squiggle counts samples, and "drawn
with karyon" does not say which.

## Authors

**Paula Ruiz-Rodriguez** and **Mireia Coscolla**, I²SysBio, University of
Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health,
Valencia, Spain.

## Background

Most of what karyon draws is a standard representation, and the ones that are
not have somebody else's idea in them. Those are worth citing in their own
right when a figure leans on them.

**Sequence logos** scaled by information content are Schneider and Stephens's.
`LogoScore::InformationContent` is that plot, and its fixed axis from zero to
`log2(K)` is WebLogo's convention.

> Schneider TD, Stephens RM. *Sequence logos: a new way to display consensus
> sequences.* Nucleic Acids Research. 1990;18(20):6097-6100.

> Crooks GE, Hon G, Chandonia JM, Brenner SE. *WebLogo: a sequence logo
> generator.* Genome Research. 2004;14(6):1188-1190.

**The enrichment and depletion logo**, `LogoTrack::edlogo`, and the Dirichlet
adaptive shrinkage behind `LogoTrack::stabilize` and the `dash` module, are
both from the Logolas paper. The scoring schemes `LogoScore` offers follow the
same source.

> Dey KK, Xie D, Stephens M. *A new sequence logo plot to highlight enrichment
> and depletion.* BMC Bioinformatics. 2018;19:473.

**Drawing only the columns that vary** is the idea
[snipit](https://github.com/aineniamh/snipit) is built around. The
implementation in `SnpTrack` and the drawing are this crate's own.

**The genome-wide significance threshold** that `genome_wide_threshold` draws
is `-log10(5e-8)`, a Bonferroni correction for a million independent tests. It
is a convention from human GWAS rather than a property of any particular study,
and it is frequently the wrong number elsewhere. Cite whatever fixed the
threshold you actually used.

## Formats

karyon the crate reads no files. karyon the command reads the line-based text
formats a genomics shell already writes, and those are defined elsewhere:

> Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis
> G, Durbin R. *The Sequence Alignment/Map format and SAMtools.*
> Bioinformatics. 2009;25(16):2078-2079.

> Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker
> RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R. *The variant call
> format and VCFtools.* Bioinformatics. 2011;27(15):2156-2158.

> Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D.
> *The Human Genome Browser at UCSC.* Genome Research. 2002;12(6):996-1006.

BED, bedGraph and the `cytoBand` table come from the UCSC browser and are
0-based and half-open. GFF3, VCF, SAM and `samtools depth` are 1-based and
inclusive. Both come out at the same place in the figure, and every reader has
a test that pins a known base through the conversion.

The default nucleotide colours are IGV's, because a figure that recolours the
bases surprises every reader.

> Robinson JT, Thorvaldsdottir H, Winckler W, Guttman M, Lander ES, Getz G,
> Mesirov JP. *Integrative Genomics Viewer.* Nature Biotechnology.
> 2011;29(1):24-26.

## License

[GPL-3.0-or-later](https://github.com/PathoGenOmics-Lab/karyon/blob/main/LICENSE).

## Next

- [Changelog](changelog.md), for what moved between the version cited and this
  one.
- [Contributing](contributing.md), for reporting a figure that came out wrong.
