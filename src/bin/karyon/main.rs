//! `karyon`, the command line front end.
//!
//! The parsing lives in the library, as [`karyon::read`], because a Rust
//! caller wanting a VCF has the same problem this binary does and should not
//! have to write the reader twice. What stays here is the part that is
//! genuinely the command line's: opening the path. Every reader takes a
//! `&str`, so nothing in `src/lib.rs` touches a disk to read, and `cargo add
//! karyon` still brings in no dependencies, because all nine formats are lines
//! of text.
//!
//! The grammar is in [`args`] and the walk from it to a figure in [`stack`],
//! which is also the only place a file is opened.
//!
//! Errors stay values until `main`, which is the only place that prints one and
//! the only place that picks an exit code. Everything under it returns a
//! `Result`, so a test can drive a whole command line without spawning a
//! process, and a command line that makes no sense is an error message rather
//! than a panic in front of the person who wrote it.
//!
//! The figure goes to standard output unless `-o` names a file, since a track
//! plot is usually one step of a pipeline rather than the end of one, and the
//! same is true on the way in: any track may read `-`, and one of them may.

use karyon::cli::{args, stack};

use std::fs;
use std::io::{self, Write};
use std::process::ExitCode;

/// What `--help` prints: enough to write a first command, on one screen.
///
/// The whole of it is `karyon help all`, [`HELP`], which is two hundred lines
/// and was what `--help` printed. Everything a reader needs to choose is here,
/// and everything they need to use one track is one `karyon help <track>`
/// away, which prints only the options that track takes.
const SHORT: &str = "\
karyon, genomic tracks on one shared coordinate axis

USAGE
    karyon <PLACE> <FILE>... [options] -o fig.svg

    karyon rpoB reads.bam genes.gff3 calls.vcf.gz -o rpoB.svg
    karyon 1 gwas.assoc --threshold genome-wide -o scan.svg
    karyon tree.nwk --traits samples.tsv --columns lineage -o tree.svg

The place comes first: a region as chr1:10,000-20,000, a gene the annotation
names, or a sequence drawn whole, named as the files name it; a figure of
trees alone needs none. Each file is a track of the kind its name says, BAM,
VCF, GFF3, GTF, BED, bedGraph, FASTA, Newick, PAF or PLINK, and .gz or not;
its options come after it. A BAM is drawn as its depth, and a track flag
chooses another kind, as --pileup reads.bam. The figure is SVG, on standard
output unless -o names one.

TRACKS, by what they draw
    signal and sequence   --coverage --windows --methylation --sequence
                          --logo --dynseq
    annotation            --features --orfs
    variation             --variants --structural --copy-number --snps
                          --matrix --heatmap --manhattan --pairs
                          --recombination
    reads and molecules   --pileup --split-reads --bisulfite --junctions
    comparison            --msa --domains --dotplot --synteny --loci
    phylogeny             --tree --tanglegram --clades
    over time and sites   --frequencies --phylodynamics --selection
                          --squiggle
    whole genome          --ideogram
    scales                --axis

FIGURE OPTIONS, anywhere on the line
    --title <TEXT>   --width <PX>   --theme light|dark   -o <FILE>

MORE
    karyon help <track>   what one track reads and the options it takes,
                          as in karyon help coverage
    karyon help all       every track and every option, in full
    https://pathogenomics-lab.github.io/karyon/
";

/// Where the guide describes each track, under the site's address.
fn guide_page(kind: args::Kind) -> &'static str {
    use args::Kind;
    match kind {
        Kind::Coverage => "tracks/signal-sequence/#coveragetrack",
        Kind::Windows => "tracks/signal-sequence/#windowtrack",
        Kind::Methylation => "tracks/signal-sequence/#methylationtrack",
        Kind::Sequence => "tracks/signal-sequence/#sequencetrack",
        Kind::Logo => "tracks/signal-sequence/#logotrack",
        Kind::Dynseq => "tracks/signal-sequence/#dynseqtrack",
        Kind::Features => "tracks/annotation/#featuretrack",
        Kind::Orfs => "tracks/annotation/#orftrack",
        Kind::Variants => "tracks/variation/#varianttrack",
        Kind::Structural => "tracks/variation/#structuraltrack",
        Kind::Pairs => "tracks/variation/#pairtrack",
        Kind::CopyNumber => "tracks/variation/#copynumbertrack",
        Kind::Snps => "tracks/variation/#snptrack",
        Kind::Matrix => "tracks/variation/#matrixtrack",
        Kind::Heatmap => "tracks/variation/#matrixtrack",
        Kind::Manhattan => "tracks/variation/#manhattantrack",
        Kind::Recombination => "tracks/signal-sequence/#coveragetrack",
        Kind::Pileup => "tracks/reads-molecules/#pileuptrack",
        Kind::SplitReads => "tracks/reads-molecules/#splitreadtrack",
        Kind::Bisulfite => "tracks/reads-molecules/#bisulfitetrack",
        Kind::Junctions => "tracks/reads-molecules/#junctiontrack",
        Kind::Msa => "tracks/comparison/#msatrack",
        Kind::Domains => "tracks/comparison/#domaintrack",
        Kind::Dotplot => "tracks/comparison/#dotplottrack",
        Kind::Synteny => "tracks/comparison/#syntenytrack",
        Kind::Loci => "tracks/comparison/#locustrack",
        Kind::Tree => "tracks/phylogeny/#treetrack",
        Kind::Tanglegram => "tracks/phylogeny/#tanglegramtrack",
        Kind::Clades => "tracks/phylogeny/#cladetrack",
        Kind::Ideogram => "tracks/whole-genome/#ideogramtrack",
        Kind::Frequencies => "tracks/evolution-surveillance/#surveillancetrack",
        Kind::Phylodynamics => "tracks/evolution-surveillance/#phylodynamictrack",
        Kind::Selection => "tracks/variation/#selectiontrack",
        Kind::Squiggle => "tracks/reads-molecules/#squiggletrack",
        Kind::Axis => "tracks/scales-keys/#axistrack",
    }
}

/// The site's address, which every page of help ends by pointing at.
const GUIDE: &str = "https://pathogenomics-lab.github.io/karyon/";

/// A value each option that takes one accepts, for asking the parser whether
/// a track takes the option at all.
///
/// The help on a track lists the options the parser answers to after that
/// track and no others, and asks the parser rather than keeping a list, so
/// the two cannot come apart. The question needs a value the option would
/// take, since an option reads its value before it looks at the track, and a
/// value it would refuse says nothing about the track. A test checks there is
/// one here for every option that takes a value.
const SAMPLES: &[(&str, &str)] = &[
    ("--label", "depth"),
    ("--against", "b.nwk"),
    ("--with-sequence", "ref.fa"),
    ("--with-tree", "t.nwk"),
    ("--ld", "lead.ld"),
    ("--links", "l.tsv"),
    ("--identity", "percent"),
    ("--modification", "m"),
    ("--context", "CpG"),
    ("--analysis", "Pfam"),
    ("--read", "r1"),
    ("--ploidy", "2"),
    ("--sample", "s1"),
    ("--traits", "s.tsv"),
    ("--columns", "a"),
    ("--height", "50"),
    ("--threshold", "7"),
    ("--projection", "circular"),
    ("--color-by", "lineage"),
    ("--support-from", "posterior"),
    ("--support-style", "both"),
    ("--mutations", "muts"),
    ("--highlight", "a"),
    ("--carrying", "A1T"),
    ("--shape", "cladogram"),
    ("--focus", "a"),
    ("--compare-to", "r"),
    ("--min-reads", "2"),
    ("--row-height", "10"),
    ("--max-rows", "10"),
    ("--aggregate", "max"),
    ("--style", "line"),
    ("--color", "#d55e00"),
    ("--format", "bedgraph"),
];

/// The lines of one section of [`HELP`], from under its heading to the next.
fn section(heading: &str) -> &'static str {
    let from = HELP
        .find(&format!("\n{heading}"))
        .and_then(|at| HELP[at + 1..].find('\n').map(|end| at + 1 + end + 1))
        .unwrap_or(HELP.len());
    let rest = &HELP[from..];
    // A heading is a line that starts at the margin with a capital.
    let to = rest
        .match_indices('\n')
        .map(|(at, _)| at + 1)
        .find(|at| rest[*at..].starts_with(|c: char| c.is_ascii_uppercase()))
        .unwrap_or(rest.len());
    &rest[..to]
}

/// The entries of a section, each a flag and every line written about it.
fn entries(section: &str) -> Vec<(&str, String)> {
    let mut found: Vec<(&str, String)> = Vec::new();
    for line in section.lines() {
        let continues = line.starts_with("     ") && !found.is_empty();
        if continues {
            let text = &mut found.last_mut().expect("an entry to continue").1;
            text.push_str(line);
            text.push('\n');
            continue;
        }
        let Some(entry) = line
            .strip_prefix("    ")
            .filter(|entry| entry.starts_with('-'))
        else {
            continue;
        };
        let flag = entry.split_whitespace().next().unwrap_or(entry);
        found.push((flag.trim_end_matches(','), format!("{line}\n")));
    }
    found
}

/// One track's own entry for an option whose shared entry speaks for several.
///
/// `karyon help all` keeps the shared entry, since it is read as a whole; the
/// help on one track gives only what that track does with the option. The
/// tree's help was seventy lines, and a reader after its two options scanned
/// past what a scan does with a threshold and a pileup with a cap.
fn said_for(kind: args::Kind, flag: &str) -> Option<&'static str> {
    use args::Kind;
    Some(match (kind, flag) {
        (Kind::Tree, "--traits") => {
            "    --traits <FILE>      a sample sheet drawn as strips beside the tips, or as
                         rings around them in a circle: a header, names in
                         column one, one column per thing known about them. A
                         folded clade shows what its tips agree on
"
        }
        (Kind::Msa | Kind::Snps | Kind::Matrix | Kind::Heatmap | Kind::Domains, "--with-tree") => {
            "    --with-tree <FILE>   a Newick tree: the rows are drawn in the order of
                         its tips, with the tree beside them
"
        }
        (Kind::Tree, "--threshold") => {
            "    --threshold <V>      the least support worth showing: the weaker values
                         --support-style would draw are left out
"
        }
        (Kind::Tree, "--row-height") => {
            "    --row-height <PX>    how tall one tip's row is; a row too short for a
                         name shrinks the name with it
"
        }
        (Kind::Tree, "--max-rows") => {
            "    --max-rows <N|all>   fold the smallest clades until the tree fits in N
                         rows; every tip stays on the figure, inside a
                         triangle saying how many it holds. No cap by default
"
        }
        (Kind::Phylodynamics, "--threshold") => {
            "    --threshold <V>      a dashed reference line at this value, as 1 for a
                         reproductive number
"
        }
        (Kind::Selection, "--threshold") => {
            "    --threshold <V>      the evidence a site needs to count as selected: a
                         p-value, 0.05 by default, or a posterior, 0.9 by
                         default, where the table holds posteriors only
"
        }
        (Kind::Pairs, "--threshold") => {
            "    --threshold <V>      the least value a pair is drawn with, as 0.2 for an r²
"
        }
        (Kind::Manhattan, "--threshold") => {
            "    --threshold <V|genome-wide>
                         the line the scan is read against, labelled with its
                         value, in the units the file is in: a p-value where it
                         holds p-values. genome-wide is p = 5e-8, a correction
                         for a million tests and wrong wherever a million were
                         not run
"
        }
        (_, "--traits") => {
            "    --traits <FILE>      a sample sheet drawn as strips beside the rows: a
                         header, names in column one, one column per thing
                         known about them
"
        }
        (_, "--max-rows") => {
            "    --max-rows <N|all>   how deep the track is drawn before it stops and
                         counts the rest; 40 by default, and all lifts it
"
        }
        _ => return None,
    })
}

/// Whether the parser takes `flag` after a `kind` track, asked with a value
/// of the kind the flag reads.
fn takes(kind: args::Kind, flag: &str) -> bool {
    let value = SAMPLES
        .iter()
        .find(|(name, _)| *name == flag)
        .map(|(_, value)| *value);
    args::takes(kind, flag, value)
}

/// The help on one topic: everything for `all`, and for a track its entry,
/// the options it takes and the page of the guide about it.
fn help_on(topic: &str) -> Result<String, String> {
    if topic == "all" {
        return Ok(HELP.to_string());
    }
    let dashed = if topic.starts_with("--") {
        topic.to_string()
    } else {
        format!("--{topic}")
    };
    let Some(kind) = args::Kind::ALL
        .iter()
        .copied()
        .find(|kind| kind.dashed() == dashed)
    else {
        let near = args::nearest_flag(&dashed)
            .filter(|flag| args::Kind::ALL.iter().any(|kind| kind.dashed() == *flag));
        return Err(match near {
            Some(flag) => format!(
                "no track is called {topic}; did you mean {}? karyon --help lists them",
                &flag[2..]
            ),
            None => format!(
                "no track is called {topic}; karyon --help lists them, and karyon help all has everything"
            ),
        });
    };

    let mut out = String::new();
    for (flag, text) in entries(section("TRACKS")) {
        if flag == kind.dashed() {
            out.push_str(&text);
        }
    }
    let options: Vec<String> = entries(section("TRACK OPTIONS"))
        .into_iter()
        .filter(|(flag, _)| takes(kind, flag))
        .map(|(flag, text)| said_for(kind, flag).map_or(text, str::to_string))
        .collect();
    if !options.is_empty() {
        out.push_str(&format!(
            "\nOPTIONS, each written after {} and before the next track\n",
            kind.dashed()
        ));
        for text in options {
            out.push_str(&text);
        }
    }
    out.push_str(
        "\nFIGURE OPTIONS, anywhere on the line: --title, --width, --theme, --no-axis,\n\
         --no-region-label, --no-legend, --rename and -o.\n",
    );
    out.push_str(&format!(
        "\nMore, with examples: {GUIDE}{}\n",
        guide_page(kind)
    ));
    Ok(out)
}

/// What `karyon help all` prints: every track and every option, in full.
const HELP: &str = "\
karyon, genomic track plots on one shared coordinate axis

USAGE
    karyon <PLACE> [FILE | TRACK FILE]... [OPTIONS]

The place comes first: a 1-based inclusive locus string, a gene the figure's
annotation names, drawn with a margin, or a sequence's name, drawn whole. A
file named on its own is a track of the kind its name says: BAM and CRAM draw
their depth, SAM its reads, VCF its calls, GFF3, GTF and BED features,
bedGraph a signal, FASTA the reference, Newick a tree, PAF synteny, and a
PLINK or REGENIE table a scan; a .gz is read as the file inside. Each track
flag starts a track of its own kind, and the flags after a track describe
that one, so the order of the words is the order of the stack. A coordinate ruler is added at the bottom
unless --axis puts one elsewhere or --no-axis leaves it out, and unless nothing
in the figure is laid on the coordinates: a phylogeny is not, so a stack of
trees gets no ruler measuring a window it is not drawn in, and a figure made
only of --tree, --tanglegram and --snps tracks takes no region at all. An
alignment, a table over time or over the sites of a gene, and a read's signal
are their own place: named nowhere, the figure is laid over all of it, and
its ruler counts columns, weeks, sites or samples rather than bases. Any
track file may be - for standard input, and one track may take it.

TRACKS
    --coverage <FILE>    per-base signal: bedGraph, samtools depth, values, or
                         the depth of a BAM's reads
    --copy-number <FILE> segmented copy number, a caller's segment table;
                         the ploidy is a track option and is required
    --dynseq <FILE>      per-base model attribution, bedGraph, drawn as the
                         bases themselves; takes --with-sequence
    --junctions <FILE>   splice junctions as arcs weighted by the reads that
                         crossed them, an aligner's SJ.out.tab
    --sequence <FILE>    the reference bases, FASTA
    --features <FILE>    genes and other intervals, BED, GFF3 or GTF; a gene
                         is drawn once, not again for its transcripts, exons
                         and CDS
    --variants <FILE>    point calls, VCF
    --windows <FILE>     a statistic in windows, bedGraph
    --manhattan <FILE>   association statistics, a table of position and value;
                         a column headed P, pvalue or p_wald is drawn as -log10
    --recombination <FILE>
                         recombination rates along a chromosome, drawn as a line
                         in cM/Mb: a genetic map as HapMap or IMPUTE writes one,
                         or a bedGraph of rates
    --tree <FILE>        a phylogeny, Newick or NEXUS
    --msa <FILE>         a multiple sequence alignment, aligned FASTA
    --snps <FILE>        the variable sites of an alignment, aligned FASTA
    --ideogram <FILE>    cytogenetic bands, a cytoBand table
    --matrix <FILE>      a value per sample per site, a table
    --heatmap <FILE>     a value per sample per window along the sequence: a
                         depth, a copy number or a methylation level, as
                         bedtools unionbedg writes it, a sequence, a start and
                         an end, then a column per sample under a header that
                         names them
    --pileup <FILE>      aligned reads, a BAM or SAM text; takes
                         --with-sequence, and colours what disagrees with it
    --synteny <FILE>     alignment ribbons between two sequences, PAF from
                         minimap2; the most-aligned target is drawn and named
    --dotplot <FILE>     the same PAF as a dot plot
    --orfs <FILE>        open reading frames in six frames, from the same FASTA
                         --sequence takes
    --logo <FILE>        a sequence logo counted from aligned FASTA
    --tanglegram <FILE>  two phylogenies face to face, Newick; this is the
                         left one and --against names the right
    --clades <FILE>      spans carried by named taxa painted onto a phylogeny,
                         GFF3 with a taxa attribute as Gubbins writes it; the
                         tree is named by --with-tree
    --loci <FILE>        gene neighbourhoods from several genomes, BED or GFF3
                         whose first column names the genome; what joins them
                         is named by --links
    --methylation <FILE> modified bases per strand, bedMethyl from modkit;
                         --modification says which one when a file holds several
    --structural <FILE>  structural calls as arcs between their breakpoints, a
                         VCF carrying symbolic alleles or SVTYPE
    --pairs <FILE>       pairs of places and a value between them: linkage as
                         PLINK writes it in a .ld table, contacts or loops as
                         BEDPE, or a table headed pos1, pos2 and a value.
                         Linkage and contacts are a triangle under the axis,
                         and a few pairs far apart are arcs
    --split-reads <FILE> molecules that aligned in pieces, a BAM or SAM carrying
                         an SA tag; only primary alignments are read
    --bisulfite <FILE>   methylation one molecule at a time, a Bismark
                         methylation extractor file; --context says which
    --domains <FILE>     protein domains on an axis of residues, an
                         InterProScan table; --analysis says which
    --frequencies <FILE> how many of each group were seen at each time, out
                         of how many, drawn as frequencies: lineages over
                         weeks, or mutations over the passages of an
                         experiment. A table headed time (or week, day, year),
                         lineage (or mutation, group), count and total
    --phylodynamics <FILE>
                         an estimate over time with its interval, as a
                         skyline, a population size or a reproductive number:
                         a table headed time, estimate (or mean, median) and
                         lower and upper where there is an interval
    --selection <FILE>   a test of selection at each site of a gene, HyPhy's
                         FEL or MEME as CSV, or a table of site, alpha and beta
                         (or dS and dN, or omega) and a p-value or a posterior
    --squiggle <FILE>    the current of a nanopore read, a SLOW5 file as
                         slow5tools view writes it, or one number per sample;
                         --read says which read
    --axis               the coordinate ruler, put where this flag sits

TRACK OPTIONS, each describing the track before it, once
    --label <TEXT>       the name in the left gutter
    --against <FILE>     the right-hand tree of a tanglegram
    --with-sequence <FILE> the reference, FASTA. A dynseq track draws its
                         letters from it and cannot do without one; a pileup
                         colours the bases that disagree with it, and draws
                         every read agreeing when it is not given
    --with-tree <FILE>   the phylogeny a clade track paints onto, Newick
    --ld <FILE>          the linkage of each variant of a scan with its lead, a
                         PLINK .ld table of the lead against its neighbours:
                         each point is coloured by its r² with the lead, and
                         the lead is a diamond with its position over it, as
                         in LocusZoom
    --links <FILE>       the homologies between the rows of a locus track,
                         BLAST tabular, or two or three columns of names
    --identity <UNIT>    percent or fraction, for a homology file whose third
                         column could be either
    --modification <CODE> m, h, a or another modkit code, for a pileup that
                         counted more than one
    --context <NAME>     CpG, CHG or CHH, for an extractor file holding more
                         than one
    --analysis <NAME>    Pfam, PANTHER or another member database, for a
                         domain table holding more than one
    --read <NAME>        which read of a SLOW5 file holding several; the first
                         by default
    --ploidy <COPIES>    where balanced sits on a copy number ladder, as in 2;
                         required, since it is not in the file
    --sample <NAME>      which sample of a segment table holding several
    --traits <FILE>      a sample sheet drawn as strips beside the rows, for
                         the tracks that have rows: a header, names in column
                         one, one column per thing known about them. A tree
                         takes one too, as strips beside its tips or as rings
                         around them, and a folded clade shows what its tips
                         agree on and nothing where they differ
    --columns <A,B,C>    which columns of that sheet to draw, in this order;
                         every column of it by default
    --height <PX>        for the tracks that do not size themselves by rows
    --threshold <V|genome-wide>
                         the line a scan is read against, in the units the
                         file is in, so a p-value where it holds p-values;
                         genome-wide is -log10(5e-8), which is a
                         correction for a million tests and wrong wherever a
                         million were not run. On a phylogeny it is the least
                         support worth showing, and hides the weaker values
                         --support-style would otherwise draw. On an estimate
                         over time it is a dashed reference line, as 1 for a
                         reproductive number; on a selection test, the p-value
                         (0.05 by default) or posterior (0.9) a site needs;
                         on pairs, the least value a pair is drawn with
    --projection <HOW>   rectangular, circular or unrooted, for a phylogeny.
                         A circle sizes itself so its tip labels clear each
                         other, up to the width of the figure, so a big tree
                         wants a wider one or fewer rows
    --color-by <KEY>     colour each branch of a phylogeny by one column of its
                         --traits sheet, or by an annotation the file carries.
                         A clade whose tips all agree takes that colour too, so
                         a lineage is a coloured clade; a key the tree does not
                         carry is refused with the keys it does
    --support-style <HOW>
                         none, symbols, labels or both, for the support values
                         a phylogeny carries; support is always in the tooltips
                         and this is what makes it readable without hovering
    --support-from <KEY> the annotation a clade's support is kept under, as
                         posterior in a BEAST tree or prob in a MrBayes one
    --mutations <KEY>    the annotation each branch of a phylogeny keeps its
                         changes under, as an annotated Newick writes them:
                         A123T, S:D614G, or either with an nt: or aa: in front;
                         needs --carrying, which is what reads them
    --highlight <NAMES>  draw a band behind each named clade, comma separated
    --carrying <CHANGE>  mark everything at or below a branch where that change
                         happened, which is everything that carries it; needs
                         --mutations to say where the changes are, and colours
                         by the answer unless --color-by says otherwise
    --shape <HOW>        phylogram or cladogram, for a phylogeny: whether a
                         branch is as long as its length says or every branch
                         is one step, which is the shape to read a topology by
    --no-scale-bar       leave out the rule in the tree's own branch-length
                         units, which a phylogram draws by default. It is not
                         the ruler along the bottom: that one measures the
                         region, and is left out of a figure holding nothing
                         but phylogenies
    --focus <NAME[,N]>   draw one clade of a phylogeny and nothing else, named
                         by its own label, by a tip inside it, or by two tips
                         it spans. A folded triangle says the pair it spans in
                         its tooltip, so what a figure calls a clade is what
                         opens it
    --compare-to <NAME>  the row every other row is read against, named as its
                         FASTA header names it. An alignment compares against
                         the consensus without it and a variable-site panel
                         against whichever record came first, which is not a
                         decision anyone made
    --no-counts          leave out the count printed beside each thing counted:
                         the differences per sample on a variable-site panel,
                         the reads over each junction arc. A count already too
                         wide for its arc was never drawn, so a zoomed out
                         figure does not change
    --min-reads <COUNT>  the fewest reads behind a methylation site or across a
                         junction for it to be drawn; how many were left out is
                         printed on the track
    --fade-by-mapq       draw a read at less than full strength the lower its
                         mapping quality, since a read that could have come
                         from anywhere should not look as solid as one that
                         could not. A faded read is drawn square, without the
                         arrowhead that says which way it ran, and its
                         mismatches stay at full strength
    --relative           read each sample of a heatmap against its own median,
                         so 1× is its usual value: a sample sequenced deeper
                         is no longer a darker row from end to end, and a
                         deletion or a duplication is what stands out
    --row-height <PX>    how tall one row is, for the tracks that size
                         themselves by rows rather than by --height; each has
                         a minimum of its own and will not be drawn under it,
                         and a row too short for a name shrinks the name with it
    --max-rows <N|all>   how deep a pileup, alignment, variable-site panel or
                         molecule grid is drawn before it stops and counts the
                         rest; 40 by default, and all lifts it. A tree takes it
                         too and answers differently: it collapses the smallest
                         clades until it fits, so every tip is still on the
                         figure inside a triangle saying how many it holds, and
                         it has no cap unless one is asked for
    --no-names           leave out the name written on or beside each thing a
                         track draws, which is not the track's own name in the
                         gutter: that one is --label
    --aggregate <HOW>    max, mean or min, when a pixel covers many bases
    --style <HOW>        area, line or bars for coverage, steps or line for
                         windows, tick or lollipop for variants, differences
                         or all for an alignment, stacked or line for
                         frequencies, triangle or arcs for pairs
    --log                a log scale, for coverage, for an estimate over time
                         and for the colours of pairs, as a contact map is read
    --color <HEX>        as in '#d55e00'
    --format <NAME>      bedgraph, depth or values for coverage, bed or gff3
                         for features and loci, when the file cannot be told
                         by looking at it

FIGURE OPTIONS
    --title <TEXT>
    --width <PX>         900 by default
    --theme <NAME>       light or dark
    --no-axis            leave out the ruler
    --no-region-label    leave out the locus printed at the top right
    --no-legend          leave out the key to the colours a tree's branches and
                         the strips of a --traits sheet are painted in, and
                         to the bases where they are blocks too narrow for
                         their letters
    --rename <FROM=TO>   read a sequence a file calls FROM as the figure's TO,
                         as --rename 1=NC_000962.3 for a PLINK table beside a
                         FASTA; several joined by commas, or the flag again
    -o, --output <FILE>  standard output by default. The figure is SVG, so a
                         name ending in .png, .pdf or another format is refused
    -h, --help
    -V, --version

COORDINATES
    BED, bedGraph and cytoBand are read 0-based and half-open. GFF3, VCF, SAM
    and samtools depth are read 1-based and inclusive. Both come out at the
    same place in the figure.

COMPRESSED AND BINARY FILES
    A file compressed with gzip or bgzip is read as the text inside it. A BAM
    is read by --coverage, --pileup and --split-reads, through the .bai beside
    it when there is one, so only the reads over the region are read. CRAM,
    BCF and bigWig are not read here; a track handed one says which command
    writes what it reads, and a pipe brings that in:

    samtools depth -a -r NC_000962.3:761000-763000 aln.bam \\
      | karyon NC_000962.3:761,000-763,000 --coverage - --label depth -o rpoB.svg

EXAMPLE
    karyon NC_000962.3:761,000-763,000 \\
      --coverage depth.bedgraph --label depth --aggregate min \\
      --sequence H37Rv.fa \\
      --features genes.gff3 --label annotation \\
      --variants calls.vcf --label variants \\
      --title 'rpoB locus' -o rpoB.svg
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("karyon: {message}");
            ExitCode::FAILURE
        }
    }
}

/// Everything `main` does, with the errors still values.
fn run(args: &[String]) -> Result<(), String> {
    // Nothing at all is someone finding out what this is, and the answer is
    // the help rather than an error telling them to ask for it.
    if args.is_empty() {
        print!("{SHORT}");
        return Ok(());
    }
    let request = args::parse(args).map_err(|error| error.to_string())?;
    let invocation = match request {
        args::Request::Help => {
            print!("{SHORT}");
            return Ok(());
        }
        args::Request::HelpOn(topic) => {
            print!("{}", help_on(&topic)?);
            return Ok(());
        }
        args::Request::Version => {
            println!("karyon {}", karyon::VERSION);
            return Ok(());
        }
        args::Request::Draw(invocation) => invocation,
    };

    // Read through `Disk`, which takes compressed files out of their wrapper
    // and reads a BAM a window at a time through its index.
    let mut disk = stack::Disk::default();
    let svg = stack::build_files(&invocation, &mut disk, |_, _| None)
        .map_err(|error| error.to_string())?;
    // A figure drawn with something its reader should know is still a figure,
    // so what it says goes to standard error and the exit stays nought.
    for note in &disk.notes {
        eprintln!("karyon: {note}");
    }
    match &invocation.output {
        Some(path) => {
            fs::write(path, svg).map_err(|error| format!("{}: {error}", path.display()))?
        }
        None => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            out.write_all(svg.as_bytes())
                .and_then(|()| out.flush())
                .map_err(|error| format!("standard output: {error}"))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn help_is_printed_rather_than_a_figure() {
        assert!(run(&["--help".to_string()]).is_ok());
    }

    #[test]
    fn a_command_line_that_makes_no_sense_is_an_error_and_not_a_panic() {
        let error = run(&["--nonsense".to_string()]).unwrap_err();
        assert!(error.contains("unknown flag"), "{error}");
    }

    /// Against the parser's own list rather than a copy of it, so a track
    /// wired up and left out of the help text is a failing test instead of a
    /// flag nobody can find.
    #[test]
    fn the_help_text_names_every_track_flag() {
        for kind in args::Kind::ALL {
            assert!(
                HELP.contains(kind.dashed()),
                "the help text does not mention {}",
                kind.dashed()
            );
        }
        // And the modifiers that name data rather than a setting: the unit a
        // column is in, and every flag that picks one of the several things a
        // file holds, taken from the parser rather than written out here.
        assert!(HELP.contains("--identity"), "--identity is undocumented");
        for kind in args::Kind::ALL {
            if let Some(flag) = kind.selector() {
                assert!(HELP.contains(flag), "{flag} is undocumented");
            }
        }
        // And every flag that carries a second file, which is not a track.
        for kind in args::Kind::ALL {
            if let Some(flag) = kind.second_flag() {
                assert!(HELP.contains(flag), "{flag} is undocumented");
            }
        }
    }

    /// Every flag the parser answers to, against the help text and the guide.
    ///
    /// The track flags above are taken from a list the parser exports, and the
    /// modifiers are not: they are match arms, and nothing but this reads them.
    /// So they are read here, out of the source, and both places that describe
    /// them are checked against the one place that implements them. Without
    /// this, `--projection` was added to the program and reached neither the
    /// help text nor the guide, and nothing said so.
    #[test]
    fn the_help_text_and_the_guide_name_every_flag_the_parser_answers_to() {
        const PARSER: &str = include_str!("../../cli/args.rs");
        const GUIDE: &str = include_str!("../../../docs/guide/cli.md");
        // A match arm on a flag, at the indentation the parse loop is written
        // at, so a flag named in a comment or a message is not mistaken for one
        // the parser answers to.
        //
        // Three shapes of arm: `"--label" =>`, `"-o" | "--output" =>`, and
        // `flag @ ("--against" | "--with-tree" | ...) =>`, which binds the
        // spelling for the one mechanism several flags share. Reading only the
        // first shape left the second and third out of this check, and so out
        // of the list a mistyped flag is matched against.
        let mut flags: Vec<String> = Vec::new();
        for line in PARSER.lines() {
            let Some(rest) = line.strip_prefix("            ") else {
                continue;
            };
            let rest = rest.strip_prefix("flag @ (").unwrap_or(rest);
            if !rest.starts_with("\"-") {
                continue;
            }
            let Some(arms) = rest.split(" =>").next() else {
                continue;
            };
            if arms.len() == rest.len() {
                continue;
            }
            for piece in arms.trim_end_matches(')').split(" | ") {
                let flag = piece.trim().trim_matches('"').to_string();
                if flag.starts_with('-') && flag != "-" && !flags.contains(&flag) {
                    flags.push(flag);
                }
            }
        }
        assert!(
            flags.len() > 40,
            "only {} flags found; the parse loop has been rewritten and this no longer reads it",
            flags.len()
        );
        // The list a mistyped flag is matched against is the same list.
        let mut listed: Vec<&str> = args::FLAGS
            .iter()
            .copied()
            .filter(|flag| !matches!(*flag, "--help" | "--version"))
            .collect();
        let mut answered: Vec<&str> = flags.iter().map(String::as_str).collect();
        listed.sort_unstable();
        answered.sort_unstable();
        assert_eq!(listed, answered, "args::FLAGS and the parse loop disagree");
        for flag in &flags {
            assert!(
                HELP.contains(flag.as_str()),
                "{flag} is not in the help text"
            );
            assert!(
                GUIDE.contains(flag.as_str()),
                "{flag} is not in docs/guide/cli.md"
            );
        }
    }

    /// Every option that takes a value takes one: given twice, to one track or
    /// to the figure, it is refused rather than the last one winning.
    ///
    /// Driven by the help text, so every option written there is checked, on
    /// a track the parser says takes it. `--highlight` is the one that adds to
    /// a list, so a second one is a second clade and not a contradiction.
    #[test]
    fn every_option_that_takes_a_value_refuses_a_second_one() {
        let mut checked = 0;
        for (flag, text) in entries(section("TRACK OPTIONS")) {
            let first = text.lines().next().unwrap_or_default();
            if flag == "--highlight" || !first.contains(" <") {
                continue;
            }
            let kind = args::Kind::ALL
                .iter()
                .copied()
                .find(|kind| takes(*kind, flag))
                .unwrap_or_else(|| panic!("no track takes {flag}"));
            let (_, value) = SAMPLES
                .iter()
                .find(|(name, _)| *name == flag)
                .unwrap_or_else(|| panic!("{flag} takes a value and has no sample"));
            let file = if kind == args::Kind::Axis {
                ""
            } else {
                "x.txt"
            };
            let line = format!(
                "chr1:1-10 {} {file} {flag} {value} {flag} {value}",
                kind.dashed()
            );
            let error = run(&line
                .split_whitespace()
                .map(String::from)
                .collect::<Vec<_>>())
            .expect_err(&line);
            assert!(error.contains("given twice"), "{line}: {error}");
            checked += 1;
        }
        for (flag, value) in [
            ("--title", "a"),
            ("--width", "500"),
            ("--theme", "dark"),
            ("-o", "a.svg"),
        ] {
            let line = format!("chr1:1-10 {flag} {value} {flag} {value}");
            let error = run(&line
                .split_whitespace()
                .map(String::from)
                .collect::<Vec<_>>())
            .expect_err(&line);
            assert!(error.contains("given twice"), "{line}: {error}");
            checked += 1;
        }
        assert!(checked > 30, "only {checked} options were checked");
    }

    /// The help a reader meets first, which has to fit where they meet it.
    /// A track's own help gives what that track does with an option, and
    /// only for an option the track takes. The tree's ran to seventy lines,
    /// through what a scan does with a threshold and a pileup with a cap.
    #[test]
    fn a_track_s_help_says_what_it_does_with_a_shared_option() {
        for kind in args::Kind::ALL {
            for (flag, _) in entries(section("TRACK OPTIONS")) {
                let Some(own) = said_for(kind, flag) else {
                    continue;
                };
                if !takes(kind, flag) {
                    continue;
                }
                assert!(
                    own.starts_with(&format!("    {flag}")),
                    "{}'s {flag} entry does not start with the flag",
                    kind.dashed()
                );
                for line in own.lines() {
                    assert!(line.chars().count() <= 80, "wider than a terminal: {line}");
                }
            }
        }
        let tree = help_on("tree").unwrap();
        for elsewhere in ["genome-wide", "pileup", "scan"] {
            assert!(
                !tree.contains(elsewhere),
                "the tree's help speaks of {elsewhere}"
            );
        }
        assert!(tree.lines().count() < 60, "{} lines", tree.lines().count());
        let scan = help_on("manhattan").unwrap();
        assert!(!scan.contains("phylogeny"), "{scan}");
        assert!(scan.contains("p = 5e-8"), "{scan}");
    }

    #[test]
    fn the_short_help_fits_on_a_screen_and_names_every_track() {
        let lines: Vec<&str> = SHORT.lines().collect();
        assert!(lines.len() <= 40, "{} lines", lines.len());
        for line in &lines {
            assert!(line.chars().count() <= 80, "wider than a terminal: {line}");
        }
        for kind in args::Kind::ALL {
            let named = SHORT.split_whitespace().any(|word| word == kind.dashed());
            assert!(named, "the short help does not name {}", kind.dashed());
        }
        assert!(SHORT.contains("karyon help <track>"));
        assert!(SHORT.contains("karyon help all"));
    }

    /// Every track has help of its own, which points at a page of the guide
    /// that exists.
    #[test]
    fn every_track_has_help_and_a_page_of_the_guide() {
        let root = env!("CARGO_MANIFEST_DIR");
        for kind in args::Kind::ALL {
            let help = help_on(&kind.dashed()[2..]).expect("a track has help");
            assert!(
                help.starts_with(&format!("    {}", kind.dashed())),
                "{help}"
            );
            let (page, anchor) = guide_page(kind)
                .split_once("/#")
                .expect("a page and an anchor");
            let path = format!("{root}/docs/{page}.md");
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
            assert!(
                text.contains(&format!("{{ #{anchor} }}")),
                "{path} has no #{anchor}"
            );
            // With the dashes or without them, and after a track flag.
            assert_eq!(help_on(kind.dashed()).unwrap(), help);
        }
    }

    /// What the help on a track lists is what the parser takes after it: a
    /// few facts, each of which a hand-kept list could get wrong.
    #[test]
    fn the_help_on_a_track_lists_the_options_it_takes_and_no_others() {
        let options = |track: &str| -> Vec<String> {
            help_on(track)
                .unwrap()
                .lines()
                .skip(1)
                .filter_map(|line| line.strip_prefix("    --"))
                .map(|rest| format!("--{}", rest.split_whitespace().next().unwrap_or(rest)))
                .collect()
        };
        let coverage = options("coverage");
        for flag in [
            "--label",
            "--height",
            "--aggregate",
            "--style",
            "--log",
            "--format",
        ] {
            assert!(
                coverage.contains(&flag.to_string()),
                "coverage lacks {flag}"
            );
        }
        for flag in ["--ploidy", "--projection", "--row-height", "--traits"] {
            assert!(
                !coverage.contains(&flag.to_string()),
                "coverage lists {flag}"
            );
        }
        assert!(options("copy-number").contains(&"--ploidy".to_string()));
        assert!(options("tree").contains(&"--projection".to_string()));
        assert!(!options("tree").contains(&"--aggregate".to_string()));
        assert!(options("pileup").contains(&"--with-sequence".to_string()));
    }

    /// A sample that an option would refuse says nothing about which tracks
    /// take it, since the value is read before the track is looked at: every
    /// sample has to be one some track accepts.
    #[test]
    fn every_sample_is_a_value_its_option_accepts_somewhere() {
        for (flag, text) in entries(section("TRACK OPTIONS")) {
            if !text.lines().next().unwrap_or_default().contains(" <") {
                continue;
            }
            let (_, value) = SAMPLES
                .iter()
                .find(|(name, _)| *name == flag)
                .unwrap_or_else(|| panic!("{flag} takes a value and has no sample"));
            let accepted = args::Kind::ALL.iter().any(|kind| {
                let mut line = vec!["chr1:1-10".to_string(), kind.dashed().to_string()];
                if *kind != args::Kind::Axis {
                    line.push("x.txt".to_string());
                }
                line.extend([flag.to_string(), value.to_string()]);
                !matches!(
                    args::parse(&line),
                    Err(args::ArgError::WrongTrack { .. } | args::ArgError::BadValue { .. })
                )
            });
            assert!(accepted, "no track takes {flag} {value}");
        }
    }

    #[test]
    fn a_topic_that_is_no_track_says_which_one_was_meant() {
        let error = help_on("coverge").unwrap_err();
        assert!(error.contains("did you mean coverage?"), "{error}");
        let error = help_on("genes").unwrap_err();
        assert!(error.contains("karyon --help lists them"), "{error}");
        assert_eq!(help_on("all").unwrap(), HELP);
    }

    #[test]
    fn nothing_at_all_and_help_after_a_track_print_help_and_not_an_error() {
        assert!(run(&[]).is_ok());
        let line: Vec<String> = "chr1:1-10 --coverage d.bg --help"
            .split_whitespace()
            .map(String::from)
            .collect();
        assert!(run(&line).is_ok());
        assert!(run(&["help".to_string(), "tree".to_string()]).is_ok());
        assert!(run(&["help".to_string(), "nothing".to_string()]).is_err());
    }

    /// The same loop for the words `--style` takes, against the parser's own
    /// list. A style wired up and left out of the help is a value a reader has
    /// no way to discover.
    #[test]
    fn the_help_text_names_every_style() {
        // Against the --style entry rather than the whole text. "all" occurs
        // nine times in this help and "line" and "bars" occur elsewhere too,
        // so a search of the whole thing passes for a style that is wired up
        // and written nowhere, which is the one failure this test exists for.
        let entry = HELP
            .split_once("    --style <HOW>")
            .expect("the help text has no --style entry")
            .1
            .split_once("\n    --")
            .expect("the --style entry runs to the end of the help")
            .0;
        for (_, word) in args::Style::ALL {
            assert!(
                entry.contains(word),
                "the --style entry does not mention {word}"
            );
        }
    }

    /// The other direction, which closes the loop. `Kind::dashed` is
    /// exhaustive so the compiler keeps the spellings honest, but `Kind::ALL`
    /// is a list, and a list is the kind of thing that falls behind. Checking
    /// the help text against it and it against the help text makes the two a
    /// pair: a track wired into the parser and written up here but left out of
    /// the list fails, which is the way round the compiler cannot see.
    #[test]
    fn the_help_text_names_no_track_the_list_has_not_got() {
        let tracks = HELP
            .split_once("\nTRACKS\n")
            .expect("the help text has no TRACKS section")
            .1
            .split_once("\nTRACK OPTIONS")
            .expect("the TRACKS section has no end")
            .0;

        let known: Vec<&str> = args::Kind::ALL.iter().map(|k| k.dashed()).collect();
        // A second-path flag is named inside the track it belongs to, since
        // that is where a reader looks for it. Taken from the parser rather
        // than written out, so a third spelling does not have to be added here.
        // Both halves of it: the file a track cannot be drawn without, and the
        // one it will take and can do without. A pileup names the second kind.
        let seconds: Vec<&str> = args::Kind::ALL
            .iter()
            .filter_map(|kind| kind.second_flag())
            .chain(
                args::Kind::ALL
                    .iter()
                    .filter_map(|kind| kind.optional_second()),
            )
            .collect();
        let selectors: Vec<&str> = args::Kind::ALL
            .iter()
            .filter_map(|kind| kind.selector())
            .collect();
        let mut found = 0;
        for word in tracks.split_whitespace() {
            // A flag named mid sentence carries the sentence's punctuation, and
            // the flag is the word rather than the comma after it.
            let word = word.trim_end_matches([',', ';', '.']);
            if !word.starts_with("--") {
                continue;
            }
            found += 1;
            assert!(
                known.contains(&word) || seconds.contains(&word) || selectors.contains(&word),
                "{word} is documented as a track and is not in Kind::ALL"
            );
        }
        assert!(
            found >= known.len(),
            "the TRACKS section found only {found}"
        );
    }
}
