# Coordinates

karyon counts every position one way: 0-based and half-open, as BED does. This
page gives that rule, the two places that count from one instead, and what each
file format's numbers turn into on the way in.
{ .k-lead }

## In short

| You have | In karyon it is |
|:--|:--|
| a BED, bedGraph, cytoBand or PAF start and end | both as they are |
| a GFF3 start and end | `start - 1`, and the end as it is |
| a VCF, SAM or `samtools depth` position | `pos - 1` |
| a locus string such as `chr1:101-200` | `Region::parse("chr1:101-200")`, which converts it |
| a codon number, such as the 450 of S450L | `CodonTrack::span_of(450)`, which counts from 1, and from the highest coordinate on the reverse strand |

The readers in `karyon::read` make these conversions for you, and the command
line reads its files with them. You make a conversion yourself only when you
build a `Feature`, a `Variant` or another value by hand.

## One convention

```rust
use karyon::Region;

let region = Region::new("chr1", 100, 200)?; // bases 100 to 199
assert_eq!(region.len(), 100);               // end - start, with no + 1
assert!(region.contains(199));               // the last base
assert!(!region.contains(200));              // the end is not in the region
```

The first base of a sequence is `0`, and `end` is one past the last base
included. Two useful things follow. The length of an interval is `end - start`,
with no correction anywhere. And two intervals touch exactly when one's `end`
equals the other's `start`, so a gap of no bases is written as zero rather than
as minus one.

Every constructor that takes a position counts this way unless its
documentation says otherwise: `Region::new`, `Feature::new`, `Variant::new`,
`Window::new`, `Read::new`, `Band::new`, `Association::new`, `MethylSite::new`,
`Junction::new`, `StructuralVariant::new`, `CladeBlock::new`,
`SplitSegment::new`, `DomainFeature::new`, `CodonTrack::new`, and the `_at`
forms on the `plot()` builder, such as `add_coverage_at` and `add_sequence_at`.

??? info "Why one convention matters"
    An off-by-one in a genomic figure does not crash anything. The figure still
    draws, the tracks still line up with each other, and the only thing wrong
    with it is that every mark sits one base from where the data put it. That
    is why the conversion happens in as few places as possible, and why the
    test suite pushes a known base through every reader and checks where it
    lands.

## The two exceptions

The 1-based, inclusive counting of samtools, IGV, GFF3 and VCF appears in two
places, and in both a person reads the number rather than a program.

**The locus string.** `Region::parse` takes the string you would paste into a
genome browser, and `to_string()` prints it back the same way. The region you
give the command line is the same string.

```rust
let region = Region::parse("chr1:101-200")?;    // what you type into IGV
assert_eq!(region.start(), 100);                // 0-based
assert_eq!(region.end(), 200);                  // exclusive
assert_eq!(region.len(), 100);
assert_eq!(region.to_string(), "chr1:101-200");
```

**The numbers printed on a figure.** The tick labels of `AxisTrack`, the codon
numbers of `CodonTrack`, the locus in the top right corner and the positions in
tooltips all count from one, so any of them can go straight into a browser's
search box. `AxisTrack` picks its ticks in 1-based space, so that the labels
come out round, and draws each one at `scale.x(pos - 1)`.

!!! warning "Converting by hand"
    A VCF `POS` or a GFF3 `start` is `pos - 1` on the way in. A GFF3 `end` goes
    in **unchanged**: a 1-based inclusive end is already one past the last base
    once counting starts at zero. A BED `start` and `end` both go in as they
    are.

```rust
use karyon::{Feature, Variant};

// A VCF line at POS 761,155.
let call = Variant::new(761_155 - 1);
assert_eq!(call.pos, 761_154);

// A GFF3 record from 759,807 to 763,325, 1-based and inclusive.
let gene = Feature::new(759_807 - 1, 763_325);
assert_eq!((gene.start, gene.end), (759_806, 763_325));
assert_eq!(gene.len(), 3_519);

// The same gene as a BED line, 759806 763325, goes in unchanged.
assert_eq!(Feature::new(759_806, 763_325), gene);
```

## The Region API

A `Region` is a sequence name and a half-open interval on it. It is all a
`Figure` needs to know about what it draws: the one `Scale` that every track is
mapped through is built from it.

| Method | Returns | Counts |
|:--|:--|:--|
| `Region::new(seq, start, end)` | `Result<Region, Error>` | 0-based, half-open |
| `Region::parse(locus)` | `Result<Region, Error>` | 1-based, inclusive |
| `seq()` | `&str` | |
| `start()` | `u64` | 0-based, the first base included |
| `end()` | `u64` | exclusive, one past the last base |
| `len()` | `u64` | `end - start`, always at least 1 |
| `contains(pos)` | `bool` | 0-based, and `false` at `end` |
| `display_start()` | `u64` | 1-based, `start + 1` |
| `display_end()` | `u64` | 1-based inclusive, equal to `end` |
| `to_string()` | `String` | `seq:display_start-display_end` |

```rust
let region = Region::parse("chr1:101-200")?;
assert_eq!(region.display_start(), 101);            // what a reader is shown
assert_eq!(region.display_end(), 200);
assert_eq!(Region::new("chr1", 100, 200)?, region); // the same window, 0-based
```

`parse` forgives what a locus string picks up in transit and refuses anything
that would change the answer. Thousands separators (`,` and `_`) and spaces
inside the numbers are ignored, and the string is split at its **last** colon,
so a sequence name that contains colons survives. `to_string()` writes the
numbers without separators, and what it writes parses back to the same region.

```rust
let region = Region::parse("NC_000962.3:761,100-761,200")?;
assert_eq!(region, Region::new("NC_000962.3", 761_099, 761_200)?);
assert_eq!(region.to_string(), "NC_000962.3:761100-761200");
assert_eq!(Region::parse("gi|123|ref|NC_1.1:5-8")?.seq(), "gi|123|ref|NC_1.1");
```

Anything else is an error rather than a guess:

```rust
assert!(Region::parse("chr1:0-200").is_err());   // 1-based coordinates start at 1
assert!(Region::parse("chr1:200-100").is_err()); // end is before start
assert!(Region::parse("chr1:100").is_err());     // expected start-end after the colon
assert!(Region::new("chr1", 100, 100).is_err()); // a figure needs at least one base
```

The first three are `Error::InvalidLocus`, which carries the string and the
reason, and the last is `Error::EmptyRegion`. Both convert into
`std::io::Error`, so parsing a region and writing its figure to a file can
share one `?`.

## A base is a span, not a point

```rust
use karyon::{Region, Scale};

let scale = Scale::new(&Region::new("chr1", 100, 200)?, 50.0, 500.0);
assert_eq!(scale.x(100), 50.0);        // left edge of the first base
assert_eq!(scale.x_center(100), 52.5); // and its middle
assert_eq!(scale.x(200), 550.0);       // the end of the region is the right edge
```

`Scale` maps a 0-based position to an x in the image. The base at position `p`
covers the pixels from `x(p)` to `x(p + 1)`, so `Scale::x` gives the left edge
of a base and `Scale::x_center` gives its middle. Which of the two a track uses
is a real choice:

| Mark | Drawn at | Because |
|:--|:--|:--|
| a variant | `x_center(pos)` | a call is an event at one base |
| a feature | from `x(start)` to `x(end)` | it is an interval, and its half-open end lands on the edge of the next base |
| a ruler tick | the left edge of the base it numbers | a ruler marks boundaries |

`AxisTrack::center_on_bases(true)` moves the ticks to the middle of the bases
they number. That suits a sequence logo or a short motif, where a base is a
column you can see rather than a fraction of a pixel.

<figure class="k-plate" markdown>
![Sixty bases of the rpoB locus at base resolution: a depth profile, the reference drawn as coloured letters, and variant lollipops standing over the middle of the bases they call](../assets/figures/example-zoom.svg){ width="900" height="221" loading="lazy" }
</figure>

## What a file's numbers become

`karyon::read` is the one place where the convention is not uniform, because
each file format defines its own. Every reader says which one it reads and has
a test that pins a known base through the conversion, and a property test
checks that nine of the formats put the same interval on the same bases,
wherever it starts. What each format should look like is on
[File formats](../guide/formats.md).

| Format | Flags | The file counts | On the way in |
|:--|:--|:--|:--|
| BED | `--features`, `--loci` | 0-based, half-open | unchanged |
| bedGraph | `--coverage`, `--windows`, `--dynseq` | 0-based, half-open | unchanged |
| cytoBand | `--ideogram` | 0-based, half-open | unchanged |
| bedMethyl, from modkit | `--methylation` | 0-based | unchanged |
| CNVkit `.cns` | `--copy-number` | 0-based, half-open | unchanged |
| PAF | `--synteny`, `--dotplot` | 0-based, half-open, on both sequences | unchanged |
| GFF3 | `--features`, `--loci`, `--clades` | 1-based, inclusive | `start - 1`, end unchanged |
| ASCAT segments, `.seg` | `--copy-number` | 1-based, inclusive | `start - 1`, end unchanged |
| STAR `SJ.out.tab` | `--junctions` | 1-based, inclusive, on the intron | `start - 1`, end unchanged |
| InterProScan TSV | `--domains` | 1-based, inclusive, in residues | `start - 1`, stop unchanged |
| VCF | `--variants` | 1-based | `POS - 1` |
| SAM | `--pileup`, `--split-reads` | 1-based | `POS - 1`, and the same for each `pos` in an `SA` tag; the end is walked from the CIGAR |
| `samtools depth` | `--coverage` | 1-based | `pos - 1` |
| Bismark methylation extractor | `--bisulfite` | 1-based | `pos - 1` |
| a table of positions and values | `--manhattan` | 1-based | `pos - 1` |
| positions along a table's header row | `--matrix` | 1-based | `pos - 1` |
| VCF of structural calls | `--structural` | `POS` is the base before the event | `POS` and `END` unchanged |
| a bare column of values | `--coverage` | no coordinates | the first value is the first base of the region |
| FASTA | `--sequence`, `--orfs`, `--with-sequence` | no coordinates | byte `n` of the record is the base at 0-based `n`, cut down to the region |
| aligned FASTA | `--msa`, `--snps`, `--logo` | no coordinates | the alignment column, counted from 0, is the position |
| Newick | `--tree`, `--tanglegram` | no coordinates | none: a phylogeny is not drawn on the axis |

Two rows need a second look. A structural call is the exception among the
1-based formats: a symbolic allele such as `<DEL>` cannot be written without a
reference base, so `POS` is the base *before* the event, which is already the
0-based start, and `END` is already the half-open end. Taking one off `POS`, as
for an ordinary VCF, would move every call a base to the left. And
`--sequence`, `--orfs` and `--with-sequence` take the only record of a FASTA
file whatever its name, or the one named like the region when the file holds
several.

Three more things follow from how the readers filter:

- **A whole-genome file is fine to hand over.** Rows naming another sequence
  are skipped, and rows outside the region are dropped rather than carried into
  a track that would not draw them. A bedGraph interval is clipped to the
  region, so a genome-wide file is never widened into memory one base at a
  time. A cytoBand table is the exception: all the bands of the named sequence
  are kept, because an ideogram draws the whole of it.
- **Overlap counts, containment does not.** An interval is kept when
  `end > region.start()` and `start < region.end()`. A gene that runs off both
  edges of the window is drawn running off both edges, which is what rpoB does
  in most figures of its locus.
- **A deletion belongs to the window even when its anchor does not.** VCF
  writes a deletion one base to the left of the bases it removes, so the reader
  tests the span that `REF` spells rather than the anchor alone. A `POS` one
  base before the window with a five-base `REF` is still a call about the
  window.

The same gene, as BED and as GFF3, comes out as the same pair of numbers:

=== "Rust"

    ```rust
    use karyon::{read, Region};

    let region = Region::parse("NC_000962.3:759,000-764,000")?;
    let bed = "NC_000962.3\t759806\t763325\trpoB\t0\t+\n";
    let gff3 = "NC_000962.3\tRefSeq\tgene\t759807\t763325\t.\t+\t.\tName=rpoB\n";
    let vcf = "NC_000962.3\t761155\t.\tC\tT\t900\tPASS\tAF=0.98\n";

    let from_bed = read::interval::features(bed, &region, None)?;
    let from_gff3 = read::interval::features(gff3, &region, None)?;
    assert_eq!(from_bed, from_gff3); // one gene, whichever file it came from
    assert_eq!(read::point::variants(vcf, &region)?[0].pos, 761_154);
    ```

=== "Command line"

    ```bash
    # rpoB.gff3, rpoB.bed and calls.vcf hold the rows above.
    # The two gene tracks line up to the base.
    karyon NC_000962.3:761,100-761,200 \
      --features rpoB.gff3 --label GFF3 \
      --features rpoB.bed --label BED \
      --variants calls.vcf -o rpoB.svg
    ```

## The region is a coordinate system

A region's name is never looked up: nothing checks that the sequence exists or
how long it is. The readers use the name to pick their rows out of a file, and
the figure prints it in its corner. Beyond that, the axis counts whatever the
data counts, and several tracks put something other than bases on it. The
region is then written in that unit.

<figure class="k-plate" markdown>
![Thirty-four variable sites from 30 kb of twelve isolates and a reference, spaced evenly as columns with each site's position printed on end beneath it, beside a tree and three trait strips](../assets/figures/example-snps.svg){ width="900" height="385" loading="lazy" }
</figure>

| Track | The axis counts | A region for it |
|:--|:--|:--|
| [`MsaTrack`](../tracks/comparison.md#msatrack) | alignment columns, gaps included | `Region::new("alignment", 0, 900)` for 900 columns |
| [`SnpTrack`](../tracks/variation.md#snptrack) | variable sites, evenly spaced | `Region::new("sites", 0, 20)` for 20 sites |
| [`SquiggleTrack`](../tracks/reads-molecules.md#squiggletrack) | samples of the raw signal, which is time | `Region::new("read", 0, samples)` |
| [`DomainTrack`](../tracks/comparison.md#domaintrack) | residues of a protein | `Region::parse("P00533:1-1210")` |

- **An alignment is not ungapped for you.** Mapping a row back to reference
  coordinates is a real operation with real decisions in it, and karyon does
  not make them behind your back. The ruler under an alignment counts columns.
- **A variable-site panel is not linear in the genome.** It drops the columns
  that agree and spaces the rest evenly, so two neighbouring columns may be
  nine bases or nine kilobases apart. Each column prints its own position
  underneath instead, counted from one like the ruler's, and an `AxisTrack`
  does not belong under the panel: `SnpTrack` answers `false` to
  `Track::on_coordinates`, so a plot of the panel alone gets no ruler.
  `SnpTrack::from_alignment` numbers alignment columns from 0, so the first
  column is labelled 1. Move every position together with `SnpTrack::offset`,
  or build the `SnpSite` values yourself, 0-based like every other position.
- **An ideogram shows where the region is, not what is in it.**
  [`IdeogramTrack`](../tracks/whole-genome.md#ideogramtrack) draws the whole
  sequence across the plot and marks the part in view, so its x is not the
  ruler's, and it answers `false` to `Track::on_coordinates` too.
- **A phylogeny is not on the axis at all.** A tree's x is a branch length or a
  depth, so `TreeTrack` and `TanglegramTrack` answer `false` to
  `Track::on_coordinates`, and a plot holding nothing but trees gets no ruler.
- **A comparison measures only its first sequence.** The region is on the
  query of `SyntenyTrack` and `DotplotTrack`; the target is measured against
  the height of the dot plot or the width of the lower bar.

### Several sequences on one axis

`Genome` goes the other way: it lays several sequences end to end and hands
back the one region that covers them, so every track works across all of them
at once.

```rust
use karyon::Genome;

let genome = Genome::new([("chr1", 1_000_000u64), ("chr2", 600_000)]);
let at = genome.at("chr2", 1_000).unwrap();
assert_eq!(at, 1_001_000);
assert_eq!(genome.locate(at), Some(("chr2", 1_000)));
assert_eq!(genome.region().seq(), "genome");
```

The cost is that the axis is a concatenation, so a distance across a join is
not a distance. [`GenomeTrack`](../tracks/whole-genome.md#genometrack) draws
where the joins are and names each sequence, because a ruler of global
coordinates would be a ruler of a coordinate system nothing else uses.
`Genome::at` answers `None` for a name the genome does not have and for a
position past the end of its sequence, which is what a file called against
another reference, or 1-based positions handed in unconverted, look like. Names
are matched exactly: `chr1` and `1` never meet.

## Codons, a third numbering

<figure class="k-plate" markdown>
![rpoB codons 439 to 465 drawn as numbered cells carrying their residues, the variants H445Y and S450L standing over the codons they change, and a ruler of bases underneath](../assets/figures/example-codons.svg){ width="880" height="169" loading="lazy" }
</figure>

A change in a coding sequence is named by residue: BRAF V600E, TP53 R175H, rpoB
S450L. [`CodonTrack`](../tracks/scales-keys.md#codontrack) is the ruler that
can be pointed at with those names. It takes the coding sequence as a 0-based
half-open span, like everything else, and numbers the codons **from 1**, like
every protein coordinate. A GFF3 or GenBank CDS therefore goes in as
`start - 1`, with its end unchanged.

```rust
use karyon::{CodonTrack, Strand};

// rpoB, forward strand, 0-based half-open.
let ruler = CodonTrack::new(759_806, 763_325, Strand::Forward);
assert_eq!(ruler.codons(), 1_173);
assert_eq!(ruler.codon_of(761_154), Some(450));
assert_eq!(ruler.span_of(450), Some((761_153, 761_156)));
```

On the reverse strand codon 1 sits at the **highest** coordinate and the
numbering runs right to left. `span_of` still reports a span low coordinate
first, on both strands, so it can go straight to a `Scale`.

```rust
let ruler = CodonTrack::new(1_000, 1_030, Strand::Reverse);
assert_eq!(ruler.span_of(1), Some((1_027, 1_030)));  // codon 1 at the top
assert_eq!(ruler.span_of(10), Some((1_000, 1_003))); // the last at the bottom
assert_eq!(ruler.codon_of(1_029), Some(1));
assert_eq!(ruler.codon_of(1_000), Some(10));
```

Hand in the reference as it is, never its reverse complement. On the reverse
strand the track complements the bases and reads them backwards itself.

```rust
// The reference reads TTACCACAT. Complemented and read backwards that is
// ATGTGGTAA: methionine, tryptophan, stop.
let cds = CodonTrack::new(0, 9, Strand::Reverse).sequence(0, b"TTACCACAT".to_vec());
assert_eq!(cds.residue_of(1), Some(b'M'));
assert_eq!(cds.residue_of(2), Some(b'W'));
assert_eq!(cds.residue_of(3), Some(b'*'));
```

`residue_of` answers `Some(b'*')` for a stop, and `None` when no sequence is
attached, when the slice does not reach that codon, or when one of its bases is
ambiguous. A codon with no letter is drawn without one rather than guessed at.

??? info "Why a track, and not a division by three"
    The partition is itself the claim. Two changes at different bases of one
    codon are competing alleles at one residue rather than a double mutant, and
    two changes in neighbouring codons are two substitutions however few bases
    apart they are. Neither statement can be made on a ruler of bases.

    Roughly half the coding sequences in any annotation run backwards, and
    numbering one of them from the wrong end is silent: the figure still draws
    and merely names the wrong residue. The chevron the track draws on codon 1
    (`CodonTrack::show_start`, on by default) is the one thing on a
    reverse-strand ruler that says where the count starts at a glance.

??? info "A coding sequence whose length is not a multiple of three"
    The partial codon left over is not counted and not drawn, since a third of
    a residue is not a residue. It sits at the end the count finishes on, which
    depends on the strand:

    ```rust
    let forward = CodonTrack::new(0, 11, Strand::Forward);
    assert_eq!(forward.codons(), 3);
    assert_eq!(forward.codon_of(9), None); // the leftover is at the high end

    let reverse = CodonTrack::new(0, 11, Strand::Reverse);
    assert_eq!(reverse.codon_of(10), Some(1));
    assert_eq!(reverse.codon_of(1), None); // and here at the low end
    ```

??? info "A genetic code that reassigns a residue"
    Translation uses NCBI table 1. Table 11, which bacteria, archaea and
    plastids use, differs from it only in which codons may start a protein, so
    the residues are identical. A table that reassigns a residue has to be
    passed in, because translating with the wrong one gives a plausible protein
    that is wrong.

    `genetic_code` takes the sixty-four residues in NCBI order, `TTT`, `TTC`,
    `TTA`, `TTG`, `TCT` and on to `GGG`, which is the `AAs` line of an NCBI
    translation table copied across unchanged.

    ```rust
    // Vertebrate mitochondrial, NCBI table 2: TGA is tryptophan, not a stop.
    let mito = CodonTrack::new(0, 3, Strand::Forward)
        .sequence(0, b"TGA".to_vec())
        .genetic_code(b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG");
    assert_eq!(mito.residue_of(1), Some(b'W'));
    ```

## Where next

<div class="grid cards" markdown>

-   **[Scale](scale.md)**

    What happens to these positions once one pixel covers more than one of them.

-   **[File formats](../guide/formats.md)**

    What each reader expects a file to look like, and what it refuses.

-   **[Track catalogue](../tracks/index.md)**

    Every track, including the ones whose axis is not counted in bases.

</div>
