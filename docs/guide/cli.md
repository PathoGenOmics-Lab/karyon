# Command line

`karyon` is a second front end onto the same library, and it reaches twenty-eight
of the thirty-six track types: twenty-seven that have a file to read, and the
coordinate ruler, which needs none. Trees drawn with metadata, maps and the
selection views are library only.

This page is the grammar: which flags start a track, which describe the one
before them, and what the command says when one of them is wrong.

Every command on this page runs in the [playground](../playground.md), which is
this same code compiled to WebAssembly: type the files into the box and the
figure is drawn by your own browser.

## The grammar is the stack

A figure is a list of tracks in the order they are drawn, and `argv` is already
an ordered list whose later words describe the earlier ones. That is exactly
what `Plot` is, so the grammar is the obvious one:

- each `--<track>` flag starts a track,
- the flags after it describe that track until the next one starts,
- **the order of the flags is the order of the stack**.

Figure flags such as `--title` are attached to nothing and may sit anywhere.
The correspondence with the Rust API is one to one, spaces instead of dots:

```text
--coverage depth.bg --label depth --aggregate min

    .add_coverage(..).label("depth").adjust(|t| t.aggregate(Aggregate::Min))
```

So the same figure is the same list of words either way:

=== "Shell"

    ```bash
    karyon NC_000962.3:761,000-763,000 \
      --coverage depth.bedgraph --label depth --aggregate min \
      --sequence H37Rv.fa \
      --features genes.gff3     --label annotation \
      --variants calls.vcf      --label variants \
      --title 'rpoB locus' -o rpoB.svg
    ```

=== "Rust"

    ```rust
    use karyon::{plot, Aggregate};

    plot("NC_000962.3:761,000-763,000")?
        .title("rpoB locus")
        .add_coverage(depth)
        .label("depth")
        .adjust(|track| track.aggregate(Aggregate::Min))
        .add_sequence(bases)
        .add_features(genes)
        .label("annotation")
        .add_variants(calls)
        .label("variants")
        .save("rpoB.svg")?;
    ```

    The library takes values rather than paths, so `depth`, `bases`, `genes` and
    `calls` are already in hand. Reading them from a bedGraph, a FASTA, a GFF3
    and a VCF is the only thing the command does that the library does not.

![A coverage profile with a dropout, the reference sequence, two gene models and variants coloured by consequence, all over one coordinate axis](../assets/figures/example.svg){ width="900" height="306" loading="lazy" }

## The region

The first argument is the region every track is drawn over, written as a locus
string: a sequence name, a colon, and a 1-based inclusive span. Commas,
underscores and spaces inside the numbers are dropped, so
`NC_000962.3:761,000-763,000` and `NC_000962.3:761000-763000` are the same
2,000 bases.

One region per figure. It is the only argument that is not a flag or a flag's
value, so it may sit anywhere on the line, but it reads better first and the
help text puts it there. A second bare word is an error rather than a silent
choice between two regions.

Every reader skips rows on another sequence and rows outside the window, so a
genome-wide file can be handed over whole and only the window comes back.

!!! warning "Coordinates"
    The locus string and the tick labels are 1-based and inclusive, the form
    samtools and IGV use. Everything else in the crate is 0-based and half-open.
    Files are read as each format defines them: BED, bedGraph and cytoBand are
    0-based half-open, and GFF3, VCF, SAM and `samtools depth` are 1-based
    inclusive. Both come out at the same place in the figure, and every reader
    has a test that pins a known base through the conversion.

## Track flags

Twenty-eight flags, twenty-seven of which take a file. `-` in place of a path means
standard input.

The parsing is not in this binary. It lives in the library as
[`karyon::read`](formats.md), and every reader there takes a `&str` rather than
a path, so nothing in the crate opens a file to read one and the dependency
count stays at zero. What the binary keeps is opening the path, and every format
it opens is line based text.

| Flag | The track | What it reads |
|:-----|:----------|:--------------|
| `--coverage <FILE>` | Per-base signal | bedGraph, `samtools depth`, or a bare column of values |
| `--copy-number <FILE>` | Segmented copy number | a segment table: CNVkit `.cns`, ASCAT, or `.seg` |
| `--dynseq <FILE>` | Per-base model attribution | bedGraph, with `--with-sequence` |
| `--junctions <FILE>` | Splice junctions as arcs | an aligner's `SJ.out.tab` |
| `--sequence <FILE>` | The reference bases | FASTA |
| `--features <FILE>` | Genes and other intervals | BED or GFF3 |
| `--variants <FILE>` | Point calls | VCF |
| `--windows <FILE>` | A statistic in windows | bedGraph |
| `--manhattan <FILE>` | Association statistics | a table of position and value |
| `--tree <FILE>` | A phylogeny | Newick |
| `--msa <FILE>` | A multiple sequence alignment | aligned FASTA |
| `--snps <FILE>` | The variable sites of an alignment | aligned FASTA |
| `--ideogram <FILE>` | Cytogenetic bands | a cytoBand table |
| `--matrix <FILE>` | A value per sample per site | a table |
| `--pileup <FILE>` | Aligned reads | SAM text, as `samtools view` writes it |
| `--synteny <FILE>` | Alignment ribbons between two sequences | PAF, as `minimap2` writes it |
| `--dotplot <FILE>` | The same alignments as a dot plot | the same PAF |
| `--orfs <FILE>` | Open reading frames in six frames | FASTA, the same file `--sequence` takes |
| `--logo <FILE>` | A sequence logo | aligned FASTA, the same file `--msa` takes |
| `--tanglegram <FILE>` | Two phylogenies face to face | Newick, and a second one named by `--against` |
| `--clades <FILE>` | Spans carried by named taxa, painted onto a phylogeny | GFF3 with a `taxa` attribute, as `Gubbins` writes it, and a tree named by `--with-tree` |
| `--loci <FILE>` | Gene neighbourhoods from several genomes | BED or GFF3 whose first column names the genome, and the homologies named by `--links` |
| `--methylation <FILE>` | Modified bases per strand | bedMethyl, as `modkit pileup` writes it |
| `--structural <FILE>` | Structural calls as arcs between their breakpoints | VCF carrying symbolic alleles or `SVTYPE` |
| `--split-reads <FILE>` | Molecules that aligned in pieces | SAM text carrying an `SA` tag |
| `--bisulfite <FILE>` | Methylation one molecule at a time | a `bismark_methylation_extractor` file |
| `--domains <FILE>` | Protein domains, on an axis of residues | an `InterProScan` table |
| `--axis` | The coordinate ruler | nothing |

Twenty-eight of the thirty-six track types are reachable here, and those are
the ones it has: twenty-seven have a file the command can put in front of
them, and the coordinate ruler needs none. The rest are library only, either
because their format is binary, because no single standard exists for what they
draw, or because nobody has written the reader yet;
[Track catalogue](../tracks.md) says which is which for each of them. What each
reader accepts, column by column, is in [Formats](formats.md).

`--orfs` and `--logo` compute rather than read: reading frames off the same
FASTA `--sequence` takes, and a logo off the same aligned FASTA `--msa` takes.
So either can be stacked under the track it was derived from without naming a
second file.

A ruler is added at the bottom without being asked for. `--axis` puts one where
the flag sits and cancels the automatic one, so writing `--axis` first is a
ruler on top and nothing at the bottom. `--no-axis` leaves it out entirely.

A track flag takes the word after it as its file, whatever that word is, so a
forgotten path swallows the next flag and the error arrives a word or two later
than the mistake.

## Track modifiers

Each of these describes the track before it. They are not all universal,
because most of them are a setting only some tracks have.

| Flag | Value | Applies to | Default |
|:-----|:------|:-----------|:--------|
| `--label <TEXT>` | any text | every track, `--axis` included | no name in the gutter |
| `--against <FILE>` | a path, or `-` | `--tanglegram` | none, and it is required |
| `--with-sequence <FILE>` | a path, or `-` | `--dynseq`, `--pileup` | required by `--dynseq`, optional for `--pileup` |
| `--with-tree <FILE>` | a path, or `-` | `--clades` | none, and it is required |
| `--links <FILE>` | a path, or `-` | `--loci` | none, and it is required |
| `--identity <UNIT>` | `percent`, `fraction` | `--loci` | worked out from the values, and refused where they cannot say |
| `--modification <CODE>` | `m`, `h`, `a`, any modkit code | `--methylation` | the one the file holds, and refused where it holds several |
| `--context <NAME>` | `CpG`, `CHG`, `CHH` | `--bisulfite` | the one the file holds, and refused where it holds several |
| `--analysis <NAME>` | `Pfam`, `PANTHER`, any member database | `--domains` | the one the file holds, and refused where it holds several |
| `--ploidy <COPIES>` | a number of copies, as in `2` | `--copy-number` | none, and it is required |
| `--sample <NAME>` | a sample the table names | `--copy-number` | the one the file holds, and refused where it holds several |
| `--traits <FILE>` | a path, or `-` | `--matrix`, `--msa`, `--snps`, `--clades`, `--domains`, `--loci`, `--tree` | no strip |
| `--columns <A,B,C>` | column names, comma separated | every track `--traits` applies to, and only after one | every column the sheet has, in its own order |
| `--max-rows <N\|all>` | a number of rows, or `all` | `--pileup`, `--msa`, `--snps`, `--bisulfite`, `--tree` | 40 for the first four, and each says how many it left out; no cap on a tree |
| `--row-height <PX>` | a number of pixels above nought | `--features`, `--msa`, `--snps`, `--matrix`, `--pileup`, `--orfs`, `--tree`, `--tanglegram`, `--clades`, `--split-reads`, `--bisulfite`, `--domains` | the track's own |
| `--projection <HOW>` | `rectangular`, `circular` or `unrooted` | `--tree` | rectangular |
| `--focus <NAME[,NAME]>` | a clade label, a tip name, or two tip names | `--tree` | the whole tree |
| `--compare-to <NAME>` | a row name, as its FASTA header spells it | `--msa`, `--snps` | the consensus for an alignment, the first record for a variable-site panel |
| `--no-counts` | none | `--snps`, `--junctions` | the counts are printed |
| `--min-reads <COUNT>` | a whole number of reads | `--methylation`, `--junctions` | each track's own floor |
| `--fade-by-mapq` | none | `--pileup` | every read at full strength |
| `--no-names` | none | `--features`, `--msa`, `--snps`, `--matrix`, `--split-reads`, `--structural`, `--bisulfite`, `--domains`, `--loci`, `--clades` | the names are drawn |
| `--threshold <V\|genome-wide>` | a number, or `genome-wide` | `--manhattan` | no line, on purpose |
| `--height <PX>` | a number of pixels | `--coverage`, `--sequence`, `--variants`, `--windows`, `--manhattan`, `--ideogram`, `--synteny`, `--dotplot`, `--axis` | the track's own |
| `--aggregate <HOW>` | `max`, `mean`, `min` | `--coverage` | `max` |
| `--style <HOW>` | `area`, `line`, `bars` for `--coverage`; `steps`, `line` for `--windows`; `lollipop`, `tick` for `--variants`; `differences`, `all` for `--msa` | `--coverage`, `--windows`, `--variants`, `--msa` | `area`, `steps`, `lollipop` and `differences` |
| `--log` | none | `--coverage` | linear |
| `--color <HEX>` | as in `'#d55e00'` | `--coverage`, `--features` | the theme's colours |
| `--format <NAME>` | `bedgraph`, `depth`, `values`, `bed`, `gff3` | `--coverage`, `--features` | told from the file |

### When the calls are too dense to tell apart

A lollipop is a stem with a ringed head, and it reads well up to a few hundred
calls. Past that the heads overlap into a smear, and the file pays for every one
of them: two hundred thousand calls over four megabases is fifty megabytes of
document for a picture eight hundred and sixty-seven pixel columns wide.

`--style tick` is the other answer. A tick is a plain vertical mark that ignores
the value, and two ticks of one colour on one pixel column are the same ink, so
only the first is drawn:

```bash
karyon chr1:1-4,000,000 --variants calls.vcf --style tick --label variants
```

The same two hundred thousand calls come out at seventy-four kilobytes, in a
twelfth of the time and a fifth of the memory. What you give up is the value,
which a tick does not show, and the tooltip, which a mark nobody can point at
alone was never going to carry. Colours are kept: two categories on one column
are two ticks.

### How deep a stack is drawn

Four tracks lay one row per record: `--pileup` a row per read, `--msa` and
`--snps` a row per sequence, `--bisulfite` a row per molecule. A pile a thousand
reads deep is a figure a thousand rows tall, which is not a figure, so each of
them stops at forty rows and says how many it left out.

Forty is a guess about the reader's screen rather than about their data, and
`--max-rows` is how it is moved:

```bash
karyon chr1:1-400 --pileup reads.sam --max-rows 10 --label reads
karyon chr1:1-400 --pileup reads.sam --max-rows all --label reads
```

On a three hundred read pileup over four hundred bases those come out 204 and
820 pixels tall against the 534 the default gives, and only the first two carry
the line saying what was dropped, because the third dropped nothing.

`all` is the word for no cap at all, and it is a word rather than a very large
number because those are not the same thing: a cap that happens not to bite
still writes the line saying nothing was dropped, and no cap at all does not.

The other row tracks are not here because they have no cap to move. A feature
track packs into as many rows as the features need and no more, so there is
nothing to raise or lower.

A tree takes the same flag and answers it differently. The four above stop
opening rows and count what they left out, which a tree cannot do: a tip is not
interchangeable with the tip below it, and cutting the list would cut a clade in
half. So a tree collapses instead, smallest clade first, until it fits, and
every tip is still on the figure inside a triangle that says how many it holds:

```bash
karyon phylo:1-1 --no-axis --tree big.nwk --max-rows 200 --label phylogeny
```

Sixty thousand tips draw a figure 900,058 pixels tall without it, which is a
thousand screens, and there is no other way down: `--row-height` floors at two
pixels, so twenty thousand tips cannot be brought under forty thousand pixels by
any setting. With `--max-rows 200` the same tree is 3,058 pixels and 102
kilobytes, and all sixty thousand tips are accounted for on it.

A tree has no cap unless one is asked for, because a phylogeny of three hundred
tips is an ordinary figure and folding it by default would fold figures nobody
asked to fold.

### One clade of a tree

`--focus` draws one clade and nothing else. It takes a clade's own label, a tip
inside it, or the two tips a clade spans:

```bash
karyon phylo:1-1 --no-axis --tree big.nwk --focus outbreak
karyon phylo:1-1 --no-axis --tree big.nwk --focus L4_D001,L4_H148
```

The pair is not something to work out by hand. A folded triangle says it: rest
a pointer on one and it reads `clade (13 tips), L4_D001 to L4_H148`, and the two
names at the end are what opens it. The tips under a node are a run, so the
first and the last of them pick out one clade and no other.

A name the tree has not got is refused against the names it has, rather than
quietly drawing the whole tree.

The [tree viewer](../tree.md) is this flag with a hand on it. Dragging and
rolling the wheel move a transform, which costs a composite and never a render,
so the picture keeps up on a million tips exactly as it does on a hundred. When
the hand stops, the rows the view is looking at are read back, the clade holding
them is worked out, and the program is asked for that clade: a wheel gesture
takes twenty thousand tips to eight thousand, then four, then two, each drawn in
full. Pulling back out returns a level at a time. Clicking a triangle is the
same thing done in one step, `--focus` on the pair that triangle prints.

### The shape a phylogeny is laid out in

`--projection` takes `rectangular`, `circular` or `unrooted`.

```bash
karyon phylo:1-1 --no-axis --tree big.nwk --projection circular
```

A circle sizes itself so its tip labels clear each other, up to the width of the
figure, so a big tree wants a wider one or fewer rows. An unrooted drawing puts
each name at the end of its own branch while they fit, and gathers them onto a
ring with a leader each once they would touch.

### Which row the others are read against

An alignment draws only the cells that disagree, and a variable-site panel keeps
only the columns where something does. Both need something to disagree with. An
alignment uses the consensus, which is a decision its own documentation calls
deliberate; a panel uses whichever record the file happened to hold first, which
is not a decision anyone made.

`--compare-to` names it, spelled as the FASTA header spells it:

```bash
karyon aln:1-900 --msa aln.fa --compare-to H37Rv --label alignment
```

A name that is not in the file is refused, and so is one that two records share,
because both builders take a row number and neither complains about a wrong one:
an alignment quietly falls back to the consensus and a panel comes out empty, and
both of those look like figures rather than mistakes.

`--style all` is the other half of the same idea. It draws every cell rather than
only the disagreements, which is what to reach for when the question is what the
sequence is rather than where it differs.

### What is known about the rows

A track drawn as rows answers which ones. `--traits` answers what they were.
It takes [a sample sheet](formats.md#the-sample-sheet), joins it to the rows by
name, and draws one narrow strip per column between the row names and the
figure:

```bash
karyon NC_000962.3:1-4,411,532 \
  --matrix genotypes.tsv --label "resistance alleles" \
  --traits samples.tsv --columns lineage,drug,depth
```

Every strip is beside the row it belongs to, whatever order the rows are in, so
a phylogeny attached to the track reorders the strips with it. Nothing in the
sheet is at a coordinate, so panning and zooming leave the strips where they
are.

Seven tracks have rows a sheet can name: `--matrix`, `--msa`, `--snps`,
`--clades`, `--domains`, `--loci` and `--tree`. A pileup has rows too, and they
are reads rather than samples, so `--traits` is refused there rather than
accepted and ignored.

A phylogeny takes one the same way, and two things about it are worth knowing.
The strips sit to the right of the tip names on a rectangular tree and become
rings outside the tips on a circular or unrooted one, and asking for them on an
unrooted tree gathers its names onto a ring, since rings need every tip at one
radius to line up against. And a folded clade says what its tips agree on: a
clade whose samples are all one lineage is drawn as that lineage, one holding
two is left empty rather than shown as either, and one tip with nothing recorded
is enough to withhold it, because a clade cannot be called uniform on the
strength of the members that happen to have been typed.

```bash
karyon phylo:1-1 --no-axis --tree big.nwk --max-rows 60 \
  --traits samples.tsv --columns lineage,country
```

Two refusals are worth knowing before they happen. A sheet whose names match
none of the rows is refused, naming the first few it did hold, because the
figure it would otherwise draw is a strip of empty outlines beside every row
and reads as "nothing is known about any of these". And a name in `--columns`
that the sheet has not got is refused with the columns it has, since that is
nearly always a spelling:

```text
karyon: --matrix samples.tsv has no column called linage; it has lineage, drug, depth
```

### A track whose data is not one file

Four modifiers carry a file rather than a setting, for the tracks whose data is
two files. A `--<track>` flag takes one path, so the second is named, and it is
named by what it means rather than by where it sits, which is the rule every
other modifier follows:

| Track | The first file | The second |
|:------|:---------------|:-----------|
| `--tanglegram` | the left-hand tree | `--against`, the right-hand tree |
| `--clades` | the blocks and their taxa | `--with-tree`, the phylogeny they are painted onto |
| `--loci` | the genes of each genome | `--links`, what joins one row to the next |
| `--dynseq` | one score per base | `--with-sequence`, the reference the letters are drawn from |
| `--pileup` | the aligned reads | `--with-sequence`, the reference they are compared against |

The last of those is the one that is optional. The first four tracks cannot be
drawn without their second file and are refused without it. A pileup can: given
no reference it draws every read agreeing, because a mismatch is a base that
differs from something, and without a reference there is nothing for a base to
differ from. Give it one and it colours what disagrees:

```bash
samtools view aln.bam NC_000962.3:761000-763000 \
  | karyon NC_000962.3:761,000-763,000 \
      --sequence H37Rv.fa --label reference \
      --pileup - --with-sequence H37Rv.fa --label reads -o pileup.svg
```

The same FASTA twice is not a mistake. The first draws the reference as a track
of its own, so the letters are on the page; the second hands the same letters to
the pileup, which is what lets it tell a mismatch from a match. Either is useful
without the other.

A tanglegram is two phylogenies:

```bash
karyon chr1:1-1000 --no-axis \
  --tanglegram before.nwk --against after.nwk --label topology -o tangle.svg
```

All three are required, and for one reason: each of these tracks draws a
finished-looking figure without its second file, and each of those figures says
something strong and false. A tanglegram given one tree twice has no crossings,
and no crossings is what a perfect result looks like:

```console
$ karyon chr1:1-1000 --tanglegram before.nwk
karyon: a tanglegram track is drawn from two files, and --against names the second
```

The two trees are named in the figure after the files they came from, since two
phylogenies side by side with nothing over them do not say which is which.
`--no-axis` is worth adding: a tanglegram has no genomic coordinates, so the
ruler underneath it measures nothing.

A clade track is the same shape, and the tree is what fixes the rows:

```bash
karyon NC_011900.1:1-2,221,315 \
  --clades gubbins.recombination_predictions.gff --with-tree tree.nwk \
  --label recombination -o clades.svg
```

`Gubbins` writes the literal `SEQUENCE` in its first column whatever the
reference was called, so a clade file naming exactly one sequence is read
whatever that sequence is named. A file naming several is a whole genome, and
then the region picks among them and the rest are counted. Either way, a file
that does hold blocks and none of them in the window says so rather than drawing
a bare tree:

```console
$ karyon NC_011900.1:2,000,000-2,100,000 --clades gubbins.gff --with-tree tree.nwk
karyon: --clades gubbins.gff: no clade blocks in NC_011900.1:2000000-2100000, though the file holds 47 on SEQUENCE
```

A locus track stacks gene neighbourhoods from several genomes. Its first column
names the genome rather than selecting it, which is the one place an interval
file is read differently here, so the file is the concatenation a shell already
produces:

```bash
cat H37Rv.bed CDC1551.bed Erdman.bed > loci.bed
karyon ESX-1:1-4,000 --loci loci.bed --links hits.tsv --label 'ESX-1' -o esx.svg
```

`--links` takes BLAST tabular output, `-outfmt 6` or `-outfmt 7`, which DIAMOND
and others write too, or two or three columns of names. It joins its names to
the genes by exact match, and that join is the thing worth watching: the names
in a search result come from the FASTA it ran against and the names in an
annotation come from its ninth column, and those are routinely not the same
strings. A join that found nothing draws every gene in every genome outlined as
having no counterpart, which reads as a discovery, so it is refused instead:

```console
$ karyon ESX-1:1-4,000 --loci loci.bed --links hits.tsv
karyon: --loci hits.tsv: no gene name in this file names anything in the loci, starting with lcl|NC_000962.3_cds_NP_215181.1_667
```

`--identity` says whether the third column is a percentage, as BLAST and DIAMOND
write it, or a fraction, as some others do. Left out it is worked out from the
values, and a file whose values are all at or below one could be either, so it
is refused by name rather than guessed at: read the wrong way round, every
ribbon in the figure becomes a perfect match and nothing fails.

`--height` is on the tracks that do not size themselves. A feature track, a
pileup, an alignment, a matrix and a tree are as tall as the number of rows
their data needs, so a height would be a number fighting the layout; a coverage
band or a ruler has a height of its own and takes one.

`--aggregate` says what a pixel column does when it covers more than one base.
It and `--log` are coverage settings, since coverage is the track that bins.
`min` is worth reaching for when a dropout is the thing not to smooth away.

`--style` spells two vocabularies with one flag, and the error names the one
that fits the track it landed on:

```console
$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --style steps
karyon: --style does not take "steps", only area, line or bars for a coverage track
```

`--color` is written into the SVG as it stands, so any colour an SVG understands
works; a hex is the spelling that survives every renderer. Every modifier except
`--label` and `--format` says so when it lands on a track that has no use for
it, rather than being accepted and doing nothing:

```console
$ karyon NC_000962.3:761,000-763,000 --features genes.gff3 --aggregate min
karyon: --aggregate means nothing to a features track
```

`--format` is the exception: it is accepted after any track and consulted only
by `--coverage` and `--features`, the two readers that have more than one file
shape to tell apart. Its words are `bedgraph` (or `bg`), `depth`, `values`,
`bed` and `gff3` (or `gff`, or `gtf`, both read as GFF3). Naming an interval
format for a coverage track is refused rather than guessed at, and naming a
signal format for a feature track leaves the guess to run.

!!! warning "`--sequence` wants the whole sequence"
    The FASTA reader takes the first record of the file and the region indexes
    into it from its first base, so a `samtools faidx ref.fa chr:start-end`
    slice is read as if it began at base 1 of the chromosome. That draws an
    empty sequence track, or the wrong bases, without an error. Hand it the
    reference, and let the region do the cutting.

## Figure options

| Flag | Effect |
|:-----|:-------|
| `--title <TEXT>` | a title above the stack |
| `--width <PX>` | the figure width, 900 by default |
| `--theme <NAME>` | `light` or `dark` |
| `--no-axis` | leave the ruler out |
| `--no-region-label` | leave out the locus printed at the top right |
| `-o`, `--output <FILE>` | write to a file rather than standard output |
| `-h`, `--help` | the whole grammar, which is the specification |
| `-V`, `--version` | the version and nothing else |

`-h` and `-V` are answered before anything else is looked at, so
`karyon --version` prints the version whatever else is on the line.

The dark theme is a selected set of colours rather than an inversion of the
light one:

```bash
karyon NC_000962.3:761,000-763,000 \
  --coverage depth.bedgraph --label depth \
  --features genes.gff3 --label annotation \
  --theme dark -o rpoB-dark.svg
```

![The same locus figure on a dark background](../assets/figures/example-dark.svg){ width="900" height="306" loading="lazy" }

## Standard input

Any track file may be `-`, and one track may take it, since there is only one
standard input to go around:

```console
$ karyon NC_000962.3:761,000-763,000 --coverage - --variants -
karyon: only one track can read from standard input
```

Header and comment lines are dropped by every reader: `#` for BED, GFF3 and VCF,
`@` for a SAM header, and blank lines throughout. A tool's output therefore
pipes in as it comes, with nothing to strip first. Tab separated and space
separated files both read.

## Output

Standard output by default, so a figure goes into a pipe or a redirect:

```bash
karyon Chr1:1-50,000 --manhattan gwas.tsv --label association > scan.svg
```

`-o` writes a file instead. It always names a file: `-o -` writes a file called
`-`, and leaving `-o` out is the way to ask for standard output.

## Binary formats

BAM, CRAM and BCF are not read here, and are not meant to be. `samtools` and
`bcftools` already write exactly what these readers take, so the pipeline is the
parser and the library keeps its zero dependencies.

Coverage from an alignment:

```bash
samtools depth -a -r chr20:1,000,000-1,050,000 sample.bam \
  | karyon chr20:1,000,000-1,050,000 --coverage - --label depth -o depth.svg
```

Reads from a CRAM, under the annotation they fall in:

```bash
samtools view -T ecoli.fa aln.cram NC_000913.3:3,423,000-3,424,000 \
  | karyon NC_000913.3:3,423,000-3,424,000 \
      --features genes.gff3 --label genes \
      --pileup - --label reads -o reads.svg
```

Calls from a BCF, which comes out of `bcftools` as VCF text with its `##`
header, dropped on the way in:

```bash
bcftools view -r Chr1:1,000,000-1,001,000 calls.bcf \
  | karyon Chr1:1,000,000-1,001,000 --variants - --label calls -o calls.svg
```

!!! note "`samtools depth` over more than one file"
    It writes one depth column per file, so four columns for two alignments,
    which is also the shape of a bedGraph. The reader tells the two apart by the
    column count and would take the first depth for an interval end. Pass
    `--format depth` and the first sample is read.

## When it fails

Every failure exits non-zero with one line on standard error, naming the flag
and the file it was reading, and the line number when a file did not say what it
claimed to. The whole figure is built before a byte of it is written, so a file
that would not parse leaves no output behind to be mistaken for a result.

The command line itself is checked before anything is opened:

```console
$ karyon --coverage depth.bedgraph
karyon: the first argument is the region, as in NC_000962.3:761,000-763,000

$ karyon NC_000962.3:0-1000 --coverage depth.bedgraph
karyon: invalid locus "NC_000962.3:0-1000": 1-based coordinates start at 1, not 0

$ karyon NC_000962.3:761,000-763,000 --label depth
karyon: --label describes the track before it, and no track has been given yet

$ karyon NC_000962.3:761,000-763,000 --coverage depth.bedgraph --aggregate median
karyon: --aggregate does not take "median", only max, mean or min
```

Then the files:

```console
$ karyon NC_000962.3:761,000-763,000 --features nowhere.bed
karyon: --features nowhere.bed: No such file or directory (os error 2)

$ karyon chr1:1-1000 --features broken.bed
karyon: --features broken.bed: line 2: end is not a number: "three-hundred-and-fifty"

$ karyon NC_000962.3:761,000-761,500 --coverage two-samples.depth
karyon: --coverage two-samples.depth: line 1: end is before start
```

A file that opened and parsed and held nothing for this window is also an error,
because an empty track is almost always the wrong region or the wrong sequence
name rather than a fact worth drawing:

```console
$ karyon chr7:1-1000 --features genes.gff3
karyon: --features genes.gff3: no features in the region
```

Rows on another sequence and rows outside the window are not errors on their
own. They are skipped, which is what lets a whole genome annotation be handed to
a two kilobase figure.

## Next

- [File formats](formats.md), for what each reader accepts, column by column.
- [Recipes](../recipes.md), for whole pipelines that end in one of these
  commands.
- [Plot API](plot.md), for the same figure written in Rust.
