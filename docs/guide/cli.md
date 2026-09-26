# Command line

Draw a figure straight from your files without writing any Rust: name a region,
then stack one track per file, in the order you want them drawn.
{ .k-lead }

## The grammar

```text
karyon <PLACE> [FILE | TRACK FILE]... [OPTIONS]
```

Four rules cover every command:

1. **The place comes first.** A region, a gene or a whole sequence:
    - a region is a locus string, 1-based and inclusive as samtools and IGV
      write it, a sequence name, a colon and a span, with commas and
      underscores in the numbers ignored, so `NC_000962.3:761,000-762,999` and
      `NC_000962.3:761000-762999` are the same 2,000 bases;
    - a gene is a name the figure's annotation gives, and the figure is that
      gene with a tenth of its length either side, titled with its name;
    - a sequence's name is the whole sequence.

    A figure made only of `--tree`, `--tanglegram` and `--snps` tracks takes no
    place, since none of them is drawn in a window.
2. **Each file, or each track flag and its file, starts a track.** A file
   named on its own is the kind of track its name says, as the table below
   has it; a flag in front chooses the kind, as `--pileup reads.bam`. Tracks
   stack from top to bottom in the order you write them.
3. **The options after a track describe that track**, up to the next one. Each
   is said once per track: a second `--label` for the same track is refused,
   since it was almost always meant for the next one. A track given no
   `--label` is called after its file.
4. **Figure options** such as `--title` and `-o` belong to no track and can go
   anywhere on the line.

| A file named | Is drawn as |
|:--|:--|
| `.bam`, `.cram` | the depth of its reads (`--coverage`); `--pileup` draws the reads |
| `.sam` | its reads (`--pileup`) |
| `.vcf`, `.bcf` | its calls (`--variants`) |
| `.gff3`, `.gff`, `.gtf`, `.bed` | features (`--features`); a `.bed` that is modkit's bedMethyl, as methylation |
| `.bedgraph`, `.bg`, `.bdg` | a signal (`--coverage`) |
| `.fa`, `.fasta`, `.fna` | the reference (`--sequence`) |
| `.aln`, `.afa` | an alignment (`--msa`) |
| `.nwk`, `.newick`, `.tree`, `.treefile` | a phylogeny (`--tree`) |
| `.paf` | synteny (`--synteny`) |
| `.assoc`, `.glm.linear`, `.regenie` | an association scan (`--manhattan`) |
| `.ld`, `.bedpe` | pairs of positions (`--pairs`) |
| `.slow5` | a read's signal (`--squiggle`) |
| a name holding `genetic_map` | a recombination map (`--recombination`) |

Any of these may end in `.gz`. A name that could be several things, `.tsv` or
`.txt`, needs its track's flag, except the genetic maps HapMap and the
imputation panels ship, whose names say what they are.

An option and its value may be written as two words or joined by `=`, as
`--label depth` or `--label=depth`.

### Asking the program

`karyon --help` fits on one screen: the grammar, every track by what it draws,
and the figure options. `karyon help coverage`, or `--help` written after a
track flag, prints one track's entry and only the options that track takes,
with a link to its page here. `karyon help all` prints everything, every track
and every option in full. A flag typed wrong is answered with the one it was
probably meant to be:

```text
$ karyon chr1:1-10 --coverge depth.bedgraph
karyon: unknown flag --coverge; did you mean --coverage?
```

A flag another tool spells, `--vcf`, `--region`, `--metadata` or `--legend`,
is answered with how karyon says the same thing:

```text
$ karyon rpoB --vcf calls.vcf.gz
karyon: unknown flag --vcf; variant calls are --variants FILE, or the VCF named on its own
```

### A worked example

```bash
karyon NC_000962.3:761,000-762,999 \
  --coverage depth.bedgraph --label depth --aggregate min --height 70 \
  --sequence H37Rv.fa --label reference \
  --features genes.gff3 --label annotation \
  --variants calls.vcf --label variants \
  --title 'rpoB locus, resistance determining region' -o rpoB.svg
```

<figure class="k-plate" markdown>
![Four bands over one kilobase ruler: read depth, the reference sequence, the rpoB gene with its resistance determining region marked, and five variants coloured missense or synonymous](../assets/figures/example.svg){ width="900" height="305" loading="lazy" }
<figcaption>One band per track flag, in the order the flags were written.</figcaption>
</figure>

Read the command one flag at a time:

- `NC_000962.3:761,000-762,999` is the window every track is drawn over.
- `--coverage depth.bedgraph` opens the first band. The three options after it
  belong to it: `--label` names it in the left gutter, `--aggregate min` keeps
  a dropout visible where one pixel covers several bases, and `--height` sets
  how tall it is.
- `--sequence`, `--features` and `--variants` open the next three bands, each
  named by the `--label` that follows it.
- `--title` and `-o` describe the figure, so where they sit does not matter.
- A coordinate ruler is added along the bottom without being asked for.

Every reader skips rows on another sequence and rows outside the window, so you
can hand over a genome-wide file and only the window is drawn. The
[playground](../playground.md) runs this same command line in your browser, on
example files you can edit.

### Spaces instead of dots

The grammar is the `plot()` builder of [the Rust API](plot.md) written with
spaces instead of dots. A track flag is an `add_` call, and each option after it
is a call on that track:

```text
--coverage depth.bedgraph --label depth --aggregate min
.add_coverage(depth).label("depth").adjust(|track| track.aggregate(Aggregate::Min))
```

The worked example, as Rust:

```rust
use karyon::{plot, Aggregate};

plot("NC_000962.3:761,000-762,999")?
    .title("rpoB locus, resistance determining region")
    .add_coverage(depth)
    .label("depth")
    .adjust(|track| track.aggregate(Aggregate::Min).height(70.0))
    .add_sequence(bases)
    .label("reference")
    .add_features(genes)
    .label("annotation")
    .add_variants(calls)
    .label("variants")
    .save("rpoB.svg")?;
```

The library takes values rather than paths, so `depth`, `bases`, `genes` and
`calls` are already in hand here. Reading them from files is the one thing the
command adds, and it does that with the readers in `karyon::read`, described in
[File formats](formats.md).

## Track flags

Thirty-five flags, one per track the command can draw. Each takes one
file, or `-` for [standard input](#standard-input), except `--axis`, which
reads nothing.

| Flag | Draws | Reads | Track |
|:--|:--|:--|:--|
| `--coverage <FILE>` | per-base signal | [bedGraph](formats.md#bedgraph), [`samtools depth`](formats.md#samtools-depth), [a bare column of values](formats.md#a-bare-column-of-values) or [a BAM](#binary-formats), whose depth it counts | [CoverageTrack](../tracks/signal-sequence.md#coveragetrack) |
| `--copy-number <FILE>` | segmented copy number | [a segment table](formats.md#the-segment-table): CNVkit `.cns`, ASCAT or `.seg` | [CopyNumberTrack](../tracks/variation.md#copynumbertrack) |
| `--dynseq <FILE>` | per-base model attribution, drawn as the bases themselves | [bedGraph](formats.md#bedgraph), with the reference from `--with-sequence` | [DynseqTrack](../tracks/signal-sequence.md#dynseqtrack) |
| `--junctions <FILE>` | splice junctions as arcs weighted by their reads | [`SJ.out.tab`](formats.md#sj-out-tab) | [JunctionTrack](../tracks/reads-molecules.md#junctiontrack) |
| `--sequence <FILE>` | the reference bases | [FASTA](formats.md#fasta) | [SequenceTrack](../tracks/signal-sequence.md#sequencetrack) |
| `--features <FILE>` | genes and other intervals, a gene drawn once | [BED](formats.md#bed), [GFF3 or GTF](formats.md#gff3) | [FeatureTrack](../tracks/annotation.md#featuretrack) |
| `--variants <FILE>` | point calls | [VCF](formats.md#vcf) | [VariantTrack](../tracks/variation.md#varianttrack) |
| `--windows <FILE>` | a statistic in windows | [bedGraph](formats.md#bedgraph) | [WindowTrack](../tracks/signal-sequence.md#windowtrack) |
| `--manhattan <FILE>` | association statistics | [a table of position and value](formats.md#the-association-table) | [ManhattanTrack](../tracks/variation.md#manhattantrack) |
| `--recombination <FILE>` | recombination rates, as a line in cM/Mb | [a genetic map, or a bedGraph of rates](formats.md#a-recombination-map) | [CoverageTrack](../tracks/signal-sequence.md#coveragetrack) |
| `--tree <FILE>` | a phylogeny | [Newick](formats.md#newick) | [TreeTrack](../tracks/phylogeny.md#treetrack) |
| `--msa <FILE>` | a multiple sequence alignment | [aligned FASTA](formats.md#aligned-fasta) | [MsaTrack](../tracks/comparison.md#msatrack) |
| `--snps <FILE>` | the variable sites of an alignment | [aligned FASTA](formats.md#aligned-fasta) | [SnpTrack](../tracks/variation.md#snptrack) |
| `--ideogram <FILE>` | cytogenetic bands | [a cytoBand table](formats.md#cytoband) | [IdeogramTrack](../tracks/whole-genome.md#ideogramtrack) |
| `--matrix <FILE>` | a value per sample per site | [a matrix table](formats.md#the-matrix-table) | [MatrixTrack](../tracks/variation.md#matrixtrack) |
| `--heatmap <FILE>` | a value per sample per window, as a heatmap | [a table of windows](formats.md#the-table-of-windows), as `bedtools unionbedg` writes it | [MatrixTrack](../tracks/variation.md#matrixtrack) |
| `--pileup <FILE>` | aligned reads | [SAM text](formats.md#sam) from `samtools view`, or [a BAM](#binary-formats) | [PileupTrack](../tracks/reads-molecules.md#pileuptrack) |
| `--synteny <FILE>` | alignment ribbons between two sequences | [PAF](formats.md#paf) from minimap2 | [SyntenyTrack](../tracks/comparison.md#syntenytrack) |
| `--dotplot <FILE>` | the same alignments as a dot plot | [PAF](formats.md#paf) | [DotplotTrack](../tracks/comparison.md#dotplottrack) |
| `--orfs <FILE>` | open reading frames in six frames | [FASTA](formats.md#fasta), the file `--sequence` takes | [OrfTrack](../tracks/annotation.md#orftrack) |
| `--logo <FILE>` | a sequence logo | [aligned FASTA](formats.md#aligned-fasta), the file `--msa` takes | [LogoTrack](../tracks/signal-sequence.md#logotrack) |
| `--tanglegram <FILE>` | two phylogenies face to face | [Newick](formats.md#newick), the left tree; `--against` names the right | [TanglegramTrack](../tracks/phylogeny.md#tanglegramtrack) |
| `--clades <FILE>` | spans carried by named taxa, painted onto a phylogeny | [GFF3 with a `taxa` attribute](formats.md#gubbins-clade-blocks), as Gubbins writes it; `--with-tree` names the tree | [CladeTrack](../tracks/phylogeny.md#cladetrack) |
| `--loci <FILE>` | gene neighbourhoods from several genomes | [BED or GFF3 whose first column names the genome](formats.md#gene-neighbourhoods); `--links` names the homologies | [LocusTrack](../tracks/comparison.md#locustrack) |
| `--methylation <FILE>` | modified bases per strand | [bedMethyl](formats.md#bedmethyl) from modkit | [MethylationTrack](../tracks/signal-sequence.md#methylationtrack) |
| `--structural <FILE>` | structural calls as arcs between their breakpoints | [VCF with symbolic alleles or `SVTYPE`](formats.md#structural-vcf) | [StructuralTrack](../tracks/variation.md#structuraltrack) |
| `--pairs <FILE>` | pairs of places and a value, as a triangle or as arcs | [PLINK's `.ld`, BEDPE or a table of pairs](formats.md#pairs-of-positions) | [PairTrack](../tracks/variation.md#pairtrack) |
| `--split-reads <FILE>` | molecules that aligned in pieces | [SAM carrying an `SA` tag](formats.md#sam-with-sa-tags), or a BAM holding them | [SplitReadTrack](../tracks/reads-molecules.md#splitreadtrack) |
| `--bisulfite <FILE>` | methylation one molecule at a time | [a Bismark methylation extractor file](formats.md#the-bismark-extractor-file) | [BisulfiteTrack](../tracks/reads-molecules.md#bisulfitetrack) |
| `--domains <FILE>` | protein domains on an axis of residues | [an InterProScan table](formats.md#the-interproscan-table) | [DomainTrack](../tracks/comparison.md#domaintrack) |
| `--frequencies <FILE>` | how often each group was seen at each time | [a table of counts over time](formats.md#counts-over-time) | [SurveillanceTrack](../tracks/evolution-surveillance.md#surveillancetrack) |
| `--phylodynamics <FILE>` | an estimate over time with its interval | [a table of estimates over time](formats.md#estimates-over-time) | [PhylodynamicTrack](../tracks/evolution-surveillance.md#phylodynamictrack) |
| `--selection <FILE>` | a test of selection at each site of a gene | [HyPhy's FEL or MEME, or a table of sites](formats.md#selection-by-site) | [SelectionTrack](../tracks/variation.md#selectiontrack) |
| `--squiggle <FILE>` | the current of a nanopore read | [SLOW5, or a column of samples](formats.md#slow5) | [SquiggleTrack](../tracks/reads-molecules.md#squiggletrack) |
| `--axis` | the coordinate ruler, where the flag sits | nothing | [AxisTrack](../tracks/scales-keys.md#axistrack) |

`--matrix` and `--heatmap` draw one track type from two shapes of table, and
`--coverage` and `--recombination` another, so thirty-three types are drawn
here. The other four are reached from Rust
only, and the [track catalogue](../tracks/index.md) lists all thirty-seven.

A few things about track flags are worth knowing before they surprise you:

- **`--orfs` and `--logo` work out their track from a file another flag also
  takes**: reading frames from the FASTA `--sequence` takes, a logo from the
  alignment `--msa` takes. Either can sit under the track it was derived from,
  reading the same file.
- **The ruler goes wherever something is measured against it.** It is added at
  the bottom of any figure holding a track laid on the coordinates. `--axis`
  puts it where the flag sits instead, and `--no-axis` leaves the automatic one
  out. A phylogeny is not laid on the coordinates (its x is branch length), nor
  is a panel of variable sites (its x is a site index) or an ideogram (its x is
  the whole chromosome), so a figure of nothing but `--tree`, `--tanglegram`,
  `--snps` and `--ideogram` gets no ruler unless you write `--axis`.
- **Some files are their own place.** An alignment, a table over time, a table
  of the sites of a gene and a read's signal need no place named: the figure is
  drawn over all of it, and its ruler counts columns, weeks, sites or samples
  rather than bases, numbering them from 1 as the file does. A place narrows
  it, as `week:10-30` or `site:50-200`. A table whose times have fractions, as
  a skyline in decimal years has, is a continuous time: its ruler and its
  tooltips write the times as the file does, to a thousandth, from nought.
- **A track flag takes the next word as its file, whatever it is.** A forgotten
  path swallows the flag after it, and the error arrives a word late:

    ```text
    $ karyon NC_000962.3:761,000-763,000 --coverage --label depth
    karyon: one region per figure, and "depth" is a second one
    ```

## Track options

Each option describes the track before it. An option given to a track that has
no use for it is refused by name rather than ignored, as in
`--aggregate means nothing to a features track`, and where an earlier file
takes it, the refusal names that file: `it is an option of reads.bam, so write
it right after reads.bam`. The one exception is `--label`, which every track
takes.

| Option | Takes | Applies to | When left out |
|:--|:--|:--|:--|
| `--label <TEXT>` | any text | every track, `--axis` included | no name in the gutter |
| `--against <FILE>` | a Newick file, or `-` | `--tanglegram` | required |
| `--with-sequence <FILE>` | a FASTA file, or `-` | `--dynseq`, `--pileup` | required by `--dynseq`; a pileup reads against the figure's `--sequence`, and with neither draws every read agreeing |
| `--with-tree <FILE>` | a Newick file, or `-` | `--clades`, `--msa`, `--snps`, `--matrix`, `--heatmap`, `--domains` | required by `--clades`; for the others the rows stay in the order of their file, and with it they take the order of its tips and the tree is drawn beside them |
| `--links <FILE>` | BLAST tabular, or two or three columns of names, or `-` | `--loci` | required |
| `--ld <FILE>` | [PLINK's `.ld`](formats.md#pairs-of-positions) of the lead against its neighbours, or `-` | `--manhattan` | every point in one colour |
| `--with-recombination <FILE>` | a [genetic map](formats.md#a-recombination-map), or `-` | `--manhattan` | no rate laid over the scan |
| `--with-moves <FILE>` | SAM or BAM as Dorado writes it with `--emit-moves`, or `-` | `--squiggle` | the current alone, with no bases over it |
| `--identity <UNIT>` | `percent` or `fraction` | `--loci` | worked out from the values, and refused when they cannot say |
| `--modification <CODE>` | `m`, `h`, `a` or another modkit code | `--methylation` | the one code in the file; refused when it holds several |
| `--context <NAME>` | `CpG`, `CHG` or `CHH` | `--bisulfite` | the one context in the file; refused when it holds several |
| `--analysis <NAME>` | `Pfam`, `PANTHER` or another member database | `--domains` | the one analysis in the file; refused when it holds several |
| `--read <NAME>` | a read the SLOW5 file names | `--squiggle` | the first read, and the command says how many the file holds |
| `--ploidy <COPIES>` | a number of copies above 0, as in `2` | `--copy-number` | required |
| `--sample <NAME>` | a sample the table names | `--copy-number` | the one sample; refused when the table holds several |
| `--traits <FILE>` | a [sample sheet](formats.md#the-sample-sheet), or `-` | `--matrix`, `--heatmap`, `--msa`, `--snps`, `--clades`, `--domains`, `--loci`, `--tree` | no strips |
| `--columns <A,B,C>` | column names, comma separated | the tracks `--traits` applies to, and only with a sheet | every column, in the sheet's order |
| `--height <PX>` | pixels | `--coverage`, `--copy-number`, `--dynseq`, `--sequence`, `--variants`, `--windows`, `--manhattan`, `--recombination`, `--ideogram`, `--synteny`, `--dotplot`, `--methylation`, `--structural`, `--pairs`, `--junctions`, `--frequencies`, `--phylodynamics`, `--selection`, `--squiggle`, `--axis` | the track's own |
| `--threshold <V|genome-wide>` | a number in the file's units, so a p-value for a file of p-values, or `genome-wide` for -log10(5e-8) on a scan | `--manhattan`; `--tree`, as the least support worth showing; `--phylodynamics`, as a dashed reference; `--selection`, as the p-value or posterior a site needs; `--pairs`, as the least value drawn; `--frequencies`, as the frequency a lineage is flagged at | no line on a scan; every support value on a tree; no reference; p = 0.05, or a posterior of 0.9; every pair; no lineage flagged |
| `--projection <HOW>` | `rectangular`, `circular` or `unrooted` | `--tree` | `rectangular` |
| `--color-by <KEY>` | a column of the `--traits` sheet, or an annotation in the file | `--tree` | one colour for every branch |
| `--support-style <HOW>` | `none`, `symbols`, `labels` or `both` | `--tree` | `none`: support is in the tooltips only |
| `--mutations <KEY>` | the annotation the changes are kept under; needs `--carrying` | `--tree` | no changes read |
| `--highlight <NAMES>` | clade names, comma separated | `--tree` | nothing highlighted |
| `--carrying <CHANGE>` | a change, as the file spells it; needs `--mutations` | `--tree` | nothing marked |
| `--shape <HOW>` | `phylogram` or `cladogram` | `--tree` | `phylogram` |
| `--no-scale-bar` | nothing | `--tree` | a scale bar on a phylogram with branch lengths |
| `--focus <NAME[,N]>` | a clade label, a tip, or two tips | `--tree` | the whole tree |
| `--compare-to <NAME>` | a row, named as its FASTA header names it | `--msa`, `--snps` | the consensus for `--msa`; the first record for `--snps` |
| `--no-counts` | nothing | `--snps`, `--junctions` | counts printed |
| `--min-reads <COUNT>` | a whole number of reads | `--methylation`, `--junctions` | 5 behind a methylation site; 1 across a junction |
| `--fade-by-mapq` | nothing | `--pileup` | every read at full strength |
| `--relative` | nothing | `--heatmap` | the values as they are |
| `--center <V>` | a number, as in `0` for a log ratio | `--heatmap` | one hue from nought up; `1` with `--relative` |
| `--growth <RISE>` | a rise in frequency from one time to the next, above 0 and at most 1, as in `0.15` | `--frequencies` | no rise flagged |
| `--min-total <N>` | a whole number of samples from 1 | `--frequencies` | every time drawn |
| `--counts` | nothing | `--frequencies` | frequencies |
| `--row-height <PX>` | pixels above 0 | `--features`, `--msa`, `--snps`, `--matrix`, `--heatmap`, `--pileup`, `--orfs`, `--tree`, `--tanglegram`, `--clades`, `--split-reads`, `--bisulfite`, `--domains` | the track's own |
| `--max-rows <N|all>` | a number of rows from 1, or `all` | `--pileup`, `--msa`, `--snps`, `--bisulfite`, `--tree` | 40 for the first four; no cap on a tree |
| `--no-names` | nothing | `--features`, `--msa`, `--snps`, `--matrix`, `--heatmap`, `--split-reads`, `--structural`, `--bisulfite`, `--domains`, `--loci`, `--clades` | names drawn |
| `--aggregate <HOW>` | `max`, `mean` or `min` | `--coverage` | `max` |
| `--style <HOW>` | `area`, `line` or `bars` for coverage; `steps` or `line` for windows; `tick` or `lollipop` for variants; `differences` or `all` for an alignment; `stacked` or `line` for frequencies; `triangle` or `arcs` for pairs | `--coverage`, `--windows`, `--variants`, `--msa`, `--frequencies`, `--pairs` | `area`, `steps`, `lollipop`, `differences` and `stacked`; for pairs, a triangle where most places were measured against the next one, and linkage always |
| `--log` | nothing | `--coverage`, `--phylodynamics`, `--pairs` | a linear scale |
| `--color <HEX>` | a colour, as in `'#d55e00'` | `--coverage`, `--features`, `--junctions`, `--phylodynamics`, `--squiggle`, `--pairs`, `--recombination` | the theme's colours |
| `--format <NAME>` | `bedgraph`, `depth` or `values` for coverage; `bed` or `gff3` for features and loci | `--coverage`, `--features`, `--loci` | told from the file |

`--height` and `--row-height` never apply to the same track. A track sized by
its rows (a feature track, a pileup, an alignment, a tree) takes `--row-height`
and grows with its data, down to a minimum row height of its own. Most other
tracks have a fixed height and take `--height`. `--logo` and `--loci` take
neither.

### Tracks drawn from two files

A track flag takes one path, so the tracks whose data is two files name the
second with an option, spelled by what the file is:

| Track | First file | Second file |
|:--|:--|:--|
| `--tanglegram` | the left-hand tree | `--against`, the right-hand tree |
| `--clades` | the blocks and the taxa carrying them | `--with-tree`, the phylogeny they are painted onto |
| `--loci` | the genes of each genome | `--links`, the homologies between neighbouring rows |
| `--dynseq` | one score per base | `--with-sequence`, the reference the letters are drawn from |
| `--pileup` | the aligned reads | `--with-sequence`, optional: the reference mismatches are read against, the figure's `--sequence` when not given |
| `--manhattan` | the scan | `--ld`, optional: the linkage of each variant with the lead, which colours the points as LocusZoom does; `--with-recombination`, optional: a genetic map laid over it, read off a scale on the right |
| `--msa`, `--snps`, `--matrix`, `--heatmap`, `--domains` | the rows | `--with-tree`, optional: the tree the rows are ordered by and drawn beside |
| `--squiggle` | the read's current | `--with-moves`, optional: the basecaller's record of the read, whose move table puts each base over its stretch of current |

The first four are refused without their second file:

```text
$ karyon chr1:1-1000 --tanglegram before.nwk
karyon: a tanglegram track is drawn from two files, and --against names the second
```

??? info "Why a missing second file is an error"
    Each of these tracks would still draw something, and what it drew would look
    finished and say something false. A tanglegram of one tree against itself
    has no crossings, which is what a perfect result looks like. A locus track
    with no homologies outlines every gene as having no counterpart, which reads
    as a discovery.

A pileup is the one that can do without. It colours the bases that disagree
with a reference, and when `--with-sequence` gives it none it reads against the
reference the figure draws:

```bash
karyon NC_000962.3:761,000-763,000 \
  --sequence H37Rv.fa --label reference \
  --pileup aln.bam --label reads -o pileup.svg
```

`--with-sequence` is for a pileup read against a reference the figure does not
draw. With neither, every read is drawn agreeing, since a mismatch is a base
that differs from something. A FASTA
given to `--sequence`, `--orfs` or `--with-sequence` that holds one record is
used whatever its header says; one that holds several is a genome, and the
record named like the region's sequence is used, or the command is refused with
the names the file does hold.

A tanglegram names its two trees after their files, so the figure says which is
which:

```bash
karyon --tanglegram before.nwk --against after.nwk --label topology -o tangle.svg
```

A locus track joins the names in `--links` to the gene names of the loci file
exactly. The names in a search result come from the FASTA it ran against and
the names in an annotation from its own columns, and they are often not the
same strings. A join that finds nothing is refused, with the first names that
failed:

```text
$ karyon locus:1-2,500 --loci loci.bed --links hits.tsv
karyon: --loci hits.tsv: no gene name in this file names anything in the loci, starting with lcl|NC_000962.3_cds_NP_218391.1_3874, lcl|NC_002755.2_cds_MT3978
```

`--identity` says whether the third column of the links is a percentage, as
BLAST and DIAMOND write it, or a fraction. Left out, a file with any value above
1 is read as percentages, and a file whose values are all at or below 1 is
refused, because either reading would draw without failing. The
[homology table](formats.md#the-homology-table) and
[gene neighbourhoods](formats.md#gene-neighbourhoods) sections have the columns.

### Choosing one of several things in a file

Four formats can hold several datasets in one file: a modkit pileup counting two
modifications, a Bismark file with three contexts, an InterProScan table from a
dozen member databases, and a segment table with several samples. A file holding
one is drawn. A file holding several, with none chosen, is refused with what it
holds:

```text
$ karyon NC_000913.3:900-1,100 --methylation dual.bed
karyon: --methylation dual.bed holds h, m, and --modification says which to draw
```

The option that chooses is `--modification`, `--context`, `--analysis` or
`--sample`. A SLOW5 file holds many reads and is the exception: its first read
is drawn, the command says how many others there are, and `--read` names
another. A methylation, bisulfite or domain band is named after what it shows
unless `--label` names it.

`--copy-number` also needs `--ploidy`, the copy number that counts as balanced.
The segment file does not say it, and a rule in the wrong place turns every gain
into a loss.

### How deep a stack of rows goes

`--pileup`, `--msa`, `--snps` and `--bisulfite` stack their data in rows: a
pileup packs its reads into rows, and the other three give each sequence or
molecule a row of its own. Each stops at 40 rows and prints how many it left
out. `--max-rows` moves the cap, and `all` removes it:

```bash
karyon chr1:1-400 --pileup reads.sam --max-rows 10 --label reads
karyon chr1:1-400 --pileup reads.sam --max-rows all --label reads
```

A tree has no cap unless you give one, and it answers differently: it folds its
smallest clades into triangles until it fits, so every tip is still on the
figure inside a triangle that says how many it holds.

```bash
karyon --tree big.nwk --max-rows 200 --label phylogeny
```

### The row others are read against

An alignment draws only the cells that disagree with the consensus, and a
variable-site panel keeps only the columns where some row disagrees with the
first record. `--compare-to` names the row to compare against instead, spelled
as its FASTA header spells it:

```bash
karyon aln:1-900 --msa aln.fa --compare-to H37Rv --label alignment
```

A name the file does not hold is refused with the names it does, and so is a
name two records share. `--style all` draws every cell of an alignment rather
than only the disagreements.

### Sample sheets beside the rows { #what-is-known-about-the-rows }

`--traits` takes a [sample sheet](formats.md#the-sample-sheet), joins it to a
track's rows by name, and draws one narrow strip per column beside them.
`--columns` picks the columns and their order:

```bash
karyon NC_000962.3:1-4,411,532 \
  --matrix genotypes.tsv --label 'resistance alleles' \
  --traits samples.tsv --columns lineage,drug,depth
```

Eight tracks take a sheet: `--matrix`, `--heatmap`, `--msa`, `--snps`,
`--clades`, `--domains`, `--loci` and `--tree`. A pileup has rows too, but they are reads,
so `--traits` is refused there. On a tree the strips sit beside the tips, or
become rings on a circular or unrooted one, and a folded clade shows what its
tips agree on and nothing where they differ. The strips are not placed at
coordinates, so they stay put when the region changes.

Two refusals to expect. A sheet whose names match none of the rows is refused
with the first names it holds. A column that is not in the sheet is refused with
the columns it has:

```text
$ karyon NC_000962.3:1-4,411,532 --matrix genotypes.tsv --traits samples.tsv --columns linage
karyon: --matrix samples.tsv has no column called linage; it has lineage, host, depth, drug
```

### Phylogenies

A `--tree` track has the most options of any track. A typical figure:

```bash
karyon --tree big.nwk --max-rows 60 \
  --traits samples.tsv --color-by lineage --support-style symbols
```

- `--projection` lays the tree out as `rectangular`, `circular` or `unrooted`.
  A circle sizes itself so its tip labels clear each other, up to the width of
  the figure, so a big tree wants a wider figure or fewer rows.
- `--shape cladogram` makes every branch one step, the shape to read a topology
  by when the branch lengths are noise or absent.
- `--color-by` colours each branch by a column of the `--traits` sheet or by an
  annotation the Newick carries. A clade whose tips all agree takes the colour
  too, so a lineage comes out as a coloured clade.
- `--support-style` makes support values readable without hovering, and
  `--threshold` hides the ones below it.
- A phylogram draws a scale bar, a rule in its own branch-length units, and
  `--no-scale-bar` leaves it out. It is not the ruler at the bottom, which
  measures the region, and a cladogram or a tree with no branch lengths draws
  none.
- `--focus` draws one clade and nothing else, named by its own label, by a tip
  inside it, or by two tips it spans. A folded triangle gives the pair in its
  tooltip, and the [tree viewer](../tree.md) opens a clade the same way.
- `--highlight` draws a band behind each named clade.

`--mutations` and `--carrying` ask who carries a change. An annotated Newick
keeps the changes on its branches, under a key the writing tool chose:

```text
(a[&mutations="A123T,S:D614G"]:0.1,b:0.2)[&mutations="C241T"]:0.3;
```

```bash
karyon --tree tree.nwk --mutations mutations --carrying S:D614G
```

`--mutations` names the key, and each of the two is refused without the other.
`A123T`, `S:D614G`, and either with `nt:` or `aa:` in front are all read.
`--carrying` then marks everything at or below a branch where that change
happened, which is every tip that carries it, and colours by the answer unless
`--color-by` says otherwise. A change that arose twice marks both clades.

A clade, tip, change or `--color-by` key the tree does not hold is refused with
what it does hold, rather than drawing the whole tree:

```text
$ karyon --tree tree.nwk --focus ERR9
karyon: --tree tree.nwk has no tip or clade called ERR9; it has ERR01, ERR02, ERR03
```

[Phylogenetics](phylogenetics.md) covers the same tree from Rust, with the
layouts and annotations in more depth.

### Dense variant calls

A lollipop is a stem with a ringed head, and it reads well up to a few hundred
calls. Past that the heads overlap, and every call still costs its share of the
file. `--style tick` draws a plain vertical mark and folds the calls of one
colour on one pixel column into a single mark:

```bash
karyon chr1:1-4,000,000 --variants calls.vcf --style tick --label variants
```

A tick does not show the value, so it gives up the height of each call and its
tooltip. Colours are kept: two categories on one pixel column are two ticks.

??? info "How much smaller"
    Two hundred thousand calls in four categories over 4 Mb came out at 61 MB
    as lollipops and 0.28 MB as ticks. A figure 900 pixels wide has fewer than
    900 pixel columns to draw in, and ticks draw at most one mark per column
    and colour.

### Telling a file's format

`--coverage` reads three formats and `--features` and `--loci` two each. The
format is told from the file (by its column count, a `##gff-version` line or
its seventh column), and `--format` overrides the guess where it would be
wrong. The words are `bedgraph`, `depth`, `values`, `bed` and `gff3`, and
[Telling formats apart](formats.md#telling-formats-apart) explains each guess.

The case that needs it most is `samtools depth` over two or more files, which
writes four columns, the same count as a bedGraph. karyon notices and says so,
and `--format depth` reads the first sample:

```bash
samtools depth -a -r NC_000962.3:761000-763000 sample1.bam sample2.bam \
  | karyon NC_000962.3:761,000-763,000 --coverage - --format depth --label sample1
```

## Figure options

| Option | What it does | Default |
|:--|:--|:--|
| `--title <TEXT>` | a title above the stack | no title |
| `--width <PX>` | the width of the figure, at most 100,000 pixels | 900 |
| `--theme <NAME>` | `light` or `dark` | `light` |
| `--no-axis` | leaves out the automatic ruler; an `--axis` track stays | a ruler at the bottom |
| `--no-region-label` | leaves out the locus printed at the top right | printed |
| `--no-legend` | leaves out the key to the colours of a tree's branches, of `--traits` strips, and of bases drawn as blocks too narrow for their letters | drawn under the figure |
| `--rename <FROM=TO>` | reads a sequence a file calls `FROM` as the figure's `TO`, as `--rename 1=NC_000962.3` for a PLINK table beside a FASTA; several joined by commas, or the flag again | each file's own names |
| `-o`, `--output <FILE>` | writes the figure to a file | standard output |
| `-h`, `--help` | prints the help that fits on a screen, or after a track flag that track's; `karyon help all` prints all of it | |
| `-V`, `--version` | prints the version | |

`-h` and `-V` are answered before anything else on the line is looked at, so
`karyon --version` prints the version whatever else you wrote.

The dark theme is a set of colours chosen for a dark background, not an
inversion of the light one:

```bash
karyon NC_000962.3:761,000-762,999 \
  --coverage depth.bedgraph --label depth --aggregate min --height 70 \
  --sequence H37Rv.fa --label reference \
  --features genes.gff3 --label annotation \
  --variants calls.vcf --label variants \
  --title 'The same locus, dark theme' --theme dark -o rpoB-dark.svg
```

<figure class="k-plate" markdown>
![The rpoB locus with depth, reference, annotation and variant bands](../assets/figures/example.svg){ width="900" height="305" loading="lazy" }
<figcaption>Drawn in this page's theme; switch the page to dark with the button at the top to see what `--theme dark` draws.</figcaption>
</figure>

## Standard input

Any file a command names can be `-`, which reads standard input: a track's own
file, and the files named by `--against`, `--with-sequence`, `--with-tree`,
`--links`, `--ld` and `--traits`. There is one standard input, so only one of them can
take it:

```text
$ karyon NC_000962.3:761,000-763,000 --coverage - --variants -
karyon: only one track can read from standard input
```

A pipe is read once and kept, so what comes through one is read as a file
would be: an alignment, a table over time or over the sites of a gene, and a
read's signal are their own place from standard input too, and a table whose
times have fractions is read as a continuous time.

The readers drop blank lines, `#` comment and header lines, and SAM `@`
headers, so a tool's output pipes in as it comes, with nothing to strip first.
Tab and space separated files both read, with the exceptions listed in
[How a file is read](formats.md#how-a-file-is-read).

## Output

The figure goes to standard output unless `-o` names a file, so it can go into a
pipe or a redirect:

```bash
karyon Chr1:1-50,000 --manhattan gwas.tsv --label association > scan.svg
```

`-o` always names a file: `-o -` writes a file called `-`. Leave `-o` out to
write to standard output.

The figure is SVG whatever the file is called, so a name that promises another
format, such as `fig.png` or `fig.pdf`, is refused rather than written as SVG
under it. Write `fig.svg` and convert it with `rsvg-convert`, Inkscape or a
browser.

The whole figure is built before any of it is written. A command that fails
writes nothing, so a figure left from an earlier run under the same name is not
replaced.

## Compressed and binary files { #binary-formats }

A file compressed with gzip or bgzip is read as the text inside it, by every
track and from standard input, so a `.vcf.gz`, a `.gff3.gz` or a `.bed.gz` is
given as it is.

A BAM is read by `--coverage`, which draws the depth of its reads, and by
`--pileup` and `--split-reads`, which draw the reads. The `.bai` beside it,
`reads.bam.bai` or `reads.bai`, says which blocks hold the reads over the
region, so a figure of one gene reads that gene's blocks and no more; without
one the file is read from its start. Depth is counted as `samtools depth -a`
counts it: every base, with reads that are unmapped, secondary, failing
quality checks or duplicates left out, and a deletion not counted as covered.

```bash
karyon NC_000962.3:761,000-763,000 --coverage aln.bam --pileup aln.bam -o rpoB.svg
```

CRAM, BCF and bigWig are not read. Hand a track one and it says what the file
is and what to write in place of its name:

```text
$ karyon chr1:1-5,000 --variants calls.bcf
karyon: --variants calls.bcf: the file is BCF, and karyon reads text; write <(bcftools view calls.bcf) where its name is, or turn it into text first
```

`<(command)` is the shell handing the command's output over as though it were
a file, in bash and zsh, so it works for any track and for a second file as
well, where `-` can be given to only one:

```bash
karyon NC_000913.3:3,423,000-3,424,000 genes.gff3 \
  --pileup <(samtools view -h -T ecoli.fa aln.cram NC_000913.3:3423000-3424000) \
  --variants <(bcftools view calls.bcf) -o reads.svg
```

!!! tip "Secondary and supplementary alignments"
    `--pileup` draws every mapped record it is given, secondary and
    supplementary ones included. Leave them out on the way in with
    `samtools view -F 0x900` if you do not want them.

## When something fails

A failing command prints one line to standard error, starting with `karyon:`,
and exits with status 1. Success exits with 0, and so do `--help` and
`--version`. The first problem stops the command, and nothing is written.

A figure drawn with something its reader should know is still drawn and still
exits with 0, and the something goes to standard error, after `karyon:` too:
a sequence no file gives the length of, drawn only as far as its rows reach; a
BAM named on its own over a window of reads, drawn as its depth; and bases too
narrow for their letters, with the `--width` that would letter them:

```text
$ karyon NC_000962.3:761,100-761,500 aln.bam H37Rv.fa -o reads.svg
karyon: aln.bam is drawn as its depth; --pileup aln.bam draws its reads
karyon: the bases are blocks of colour at this width, too narrow for their letters; --width 3000 draws the letters
```

The command line is checked before any file is opened:

```text
$ karyon --coverage depth.bedgraph
karyon: the first argument is the region, as in NC_000962.3:761,000-763,000; only a figure of --tree, --tanglegram, --snps, --msa and --logo tracks goes without one

$ karyon NC_000962.3:0-1000 --coverage depth.bedgraph
karyon: invalid locus "NC_000962.3:0-1000": 1-based coordinates start at 1, not 0

$ karyon NC_000962.3:761,000-763,000 --label depth
karyon: --label describes the track before it, and no track has been given yet

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --label depth --label reads
karyon: --label is given twice to one coverage track, which takes one; a flag describes the track written before it

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --height NaN
karyon: --height does not take "NaN", only a number of pixels above nought, as in 80

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph -o rpoB.png
karyon: rpoB.png names a PNG file, and karyon writes SVG: write the figure to a file ending in .svg and convert it, with rsvg-convert, Inkscape or a browser

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --aggregate median
karyon: --aggregate does not take "median", only max, mean or min

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --style steps
karyon: --style does not take "steps", only area, line or bars for a coverage track

$ karyon NC_000962.3:761,000-763,000 --features genes.gff3 --aggregate min
karyon: --aggregate means nothing to a features track

$ karyon NC_000962.3:761,000-763,000 depth.bedgraph genes.gff3 --aggregate min
karyon: --aggregate means nothing to a features track; it is an option of depth.bedgraph, so write it right after depth.bedgraph

$ karyon chr8:1-1000 --copy-number segments.cns
karyon: --copy-number needs --ploidy, since where balanced sits is not in the file
```

Then the files, in the order the tracks are stacked. Every message names the
flag and the file, and the line number when a line did not say what it should.
The line number counts comment and blank lines, so it is the line number in
your editor:

```text
$ karyon NC_000962.3:761,000-763,000 --features nowhere.bed
karyon: --features nowhere.bed: No such file or directory (os error 2)

$ karyon chr1:1-1000 --features broken.bed
karyon: --features broken.bed: line 2: end is not a number: "three-hundred-and-fifty"

$ karyon 1:1-5,000 --variants calls.vcf
karyon: --variants calls.vcf: no variants in 1:1-5000, though the file holds 12 on chr1

$ karyon aln:1-9 --msa aln.fa
karyon: --msa aln.fa: line 3: an alignment has every record the same length, and "sample_02" is 1 shorter than "sample_01", which is 9 columns
```

A file that parsed and held nothing for the window is an error too, because an
empty band is almost always a wrong region or a wrong sequence name rather than
a fact worth drawing. The message takes one of two shapes, the second when the
file held records and none of them reached the window:

```text
$ karyon chr7:1-1000 --features genes.gff3
karyon: --features genes.gff3: no features in the region

$ karyon NC_011900.1:5,000-9,000 --clades gubbins.gff --with-tree tree.nwk
karyon: --clades gubbins.gff: no clade blocks in NC_011900.1:5000-9000, though the file holds 1 on SEQUENCE
```

Check the sequence name first. It has to match the region's exactly: `chr1`,
`1` and `NC_000001.11` are three different sequences to every reader.

## Where next

<div class="grid cards" markdown>

-   **[File formats](formats.md)**

    What each reader takes from a file, column by column, and where its
    coordinates land.

-   **[Recipes](../recipes.md)**

    Whole pipelines that end in one of these commands.

-   **[The Rust API](plot.md)**

    The same figures built with the `plot()` builder.

-   **[Track catalogue](../tracks/index.md)**

    Every track type, including the ones only Rust can reach.

</div>
