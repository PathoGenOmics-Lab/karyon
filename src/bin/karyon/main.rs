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
    --pileup <FILE>      aligned reads, SAM text from samtools view
    --synteny <FILE>     alignment ribbons between two sequences, PAF from
                         minimap2; the most-aligned target is drawn and named
    --dotplot <FILE>     the same PAF as a dot plot
    --orfs <FILE>        open reading frames in six frames, from the same FASTA
                         --sequence takes
    --logo <FILE>        a sequence logo counted from aligned FASTA
    --tanglegram <FILE>  two phylogenies face to face, Newick; this is the
                         left one and --against names the right
    --axis               the coordinate ruler, put where this flag sits

TRACK OPTIONS, each describing the track before it
    --label <TEXT>       the name in the left gutter
    --against <FILE>     the other file, for a track drawn from two of them:
                         the right-hand tree of a tanglegram
    --height <PX>        for the tracks that do not size themselves by rows
    --aggregate <HOW>    max, mean or min, when a pixel covers many bases
    --style <HOW>        area, line, bars, or steps for a window track
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
            println!("karyon {}", env!("CARGO_PKG_VERSION"));
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
        // And the flag that carries a second file, which is not a track.
        assert!(HELP.contains("--against"), "--against is undocumented");
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
        let mut found = 0;
        for word in tracks.split_whitespace() {
            if !word.starts_with("--") {
                continue;
            }
            found += 1;
            assert!(
                known.contains(&word) || word == "--against",
                "{word} is documented as a track and is not in Kind::ALL"
            );
        }
        assert!(
            found >= known.len(),
            "the TRACKS section found only {found}"
        );
    }
}
