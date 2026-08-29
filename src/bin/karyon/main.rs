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

/// What `--help` prints.
const HELP: &str = "\
karyon, genomic track plots on one shared coordinate axis

USAGE
    karyon <REGION> [TRACK...] [OPTIONS]

The region comes first, as a 1-based inclusive locus string. Each track flag
starts a track and the flags after it describe that one, so the order of the
flags is the order of the stack. A coordinate ruler is added at the bottom
unless --axis puts one elsewhere or --no-axis leaves it out. Any track file
may be - for standard input, and one track may take it.

TRACKS
    --coverage <FILE>    per-base signal: bedGraph, samtools depth, or values
    --copy-number <FILE> segmented copy number, a caller's segment table;
                         the ploidy is a track option and is required
    --dynseq <FILE>      per-base model attribution, bedGraph, drawn as the
                         bases themselves; takes --with-sequence
    --junctions <FILE>   splice junctions as arcs weighted by the reads that
                         crossed them, an aligner's SJ.out.tab
    --sequence <FILE>    the reference bases, FASTA
    --features <FILE>    genes and other intervals, BED or GFF3
    --variants <FILE>    point calls, VCF
    --windows <FILE>     a statistic in windows, bedGraph
    --manhattan <FILE>   association statistics, a table of position and value
    --tree <FILE>        a phylogeny, Newick
    --msa <FILE>         a multiple sequence alignment, aligned FASTA
    --snps <FILE>        the variable sites of an alignment, aligned FASTA
    --ideogram <FILE>    cytogenetic bands, a cytoBand table
    --matrix <FILE>      a value per sample per site, a table
    --pileup <FILE>      aligned reads, SAM text from samtools view; takes
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
    --split-reads <FILE> molecules that aligned in pieces, SAM carrying an SA
                         tag; only primary alignments are read
    --bisulfite <FILE>   methylation one molecule at a time, a Bismark
                         methylation extractor file; --context says which
    --domains <FILE>     protein domains on an axis of residues, an
                         InterProScan table; --analysis says which
    --axis               the coordinate ruler, put where this flag sits

TRACK OPTIONS, each describing the track before it
    --label <TEXT>       the name in the left gutter
    --against <FILE>     the right-hand tree of a tanglegram
    --with-sequence <FILE> the reference, FASTA. A dynseq track draws its
                         letters from it and cannot do without one; a pileup
                         colours the bases that disagree with it, and draws
                         every read agreeing when it is not given
    --with-tree <FILE>   the phylogeny a clade track paints onto, Newick
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
                         file is in; genome-wide is -log10(5e-8), which is a
                         correction for a million tests and wrong wherever a
                         million were not run
    --projection <HOW>   rectangular, circular or unrooted, for a phylogeny.
                         A circle sizes itself so its tip labels clear each
                         other, up to the width of the figure, so a big tree
                         wants a wider one or fewer rows
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
                         or all for an alignment
    --log                a log scale
    --color <HEX>        as in '#d55e00'
    --format <NAME>      bedgraph, depth, values, bed or gff3, when the file
                         cannot be told by looking at it

FIGURE OPTIONS
    --title <TEXT>
    --width <PX>         900 by default
    --theme <NAME>       light or dark
    --no-axis            leave out the ruler
    --no-region-label    leave out the locus printed at the top right
    -o, --output <FILE>  standard output by default
    -h, --help
    -V, --version

COORDINATES
    BED, bedGraph and cytoBand are read 0-based and half-open. GFF3, VCF, SAM
    and samtools depth are read 1-based and inclusive. Both come out at the
    same place in the figure.

BINARY FORMATS
    BAM, CRAM and BCF are not read here. They come in through a pipe, since
    samtools and bcftools already write what these readers take:

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
    let request = args::parse(args).map_err(|error| error.to_string())?;
    let invocation = match request {
        args::Request::Help => {
            print!("{HELP}");
            return Ok(());
        }
        args::Request::Version => {
            println!("karyon {}", karyon::VERSION);
            return Ok(());
        }
        args::Request::Draw(invocation) => invocation,
    };

    let svg =
        stack::build(&invocation, stack::open_from_disk).map_err(|error| error.to_string())?;
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
        let mut flags: Vec<String> = Vec::new();
        for line in PARSER.lines() {
            let Some(rest) = line.strip_prefix("            \"--") else {
                continue;
            };
            let Some(arms) = rest.split(" =>").next() else {
                continue;
            };
            if arms.len() == rest.len() {
                continue;
            }
            let whole = format!("\"--{arms}");
            for piece in whole.split(" | ") {
                let flag = piece.trim().trim_matches('"').to_string();
                if flag.starts_with("--") && !flags.contains(&flag) {
                    flags.push(flag);
                }
            }
        }
        assert!(
            flags.len() > 40,
            "only {} flags found; the parse loop has been rewritten and this no longer reads it",
            flags.len()
        );
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
