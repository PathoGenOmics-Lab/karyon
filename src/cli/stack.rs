//! Turning a parsed command line into a figure.
//!
//! Every arm here is one track of the stack, built in the order the flags were
//! written, which is the whole reason the grammar in [`args`](crate::cli::args) looks the
//! way it does.
//!
//! The tracks are built and then handed over with
//! [`Plot::add_track`](crate::Plot::add_track) rather than through the `add_`
//! methods, because a command line always has the settings in hand before the
//! track exists: `--label` and `--height` have already been read by the time
//! this runs, so there is nothing left for `Plot::label` to do afterwards.
//!
//! A file that opened, parsed, and held nothing on the sequence in the window
//! is an error here rather than a track with nothing in it. An empty lane reads
//! as a stretch of the genome where there is no data, which is a different
//! claim from a file that names another sequence or a locus typed one digit
//! out, and the figure has no way to tell the reader which of them happened.
//! Every error names the flag that asked for the file and the file it was,
//! since a stack is as many files as it has tracks.

use std::fmt;
use std::fs;
use std::io::{self, Read as _, Seek as _};
use std::path::Path;

use crate::{
    Aggregate, BisulfiteTrack, CladeTrack, CopyNumberTrack, CoverageTrack, DomainTrack,
    DotplotTrack, DynseqTrack, FeatureTrack, Figure, GenotypeTrack, IdeogramTrack, JunctionTrack,
    LocusTrack, LogoTrack, ManhattanTrack, MatrixTrack, MethylationTrack, MsaSequence, MsaTrack,
    OrfTrack, PairStyle, PairTrack, PhylodynamicScale, PhylodynamicTrack, PileupTrack, Plot,
    Region, SelectionEvidence, SelectionTrack, SequenceTrack, SnpTrack, SplitReadTrack,
    SquiggleTrack, StructuralTrack, SurveillanceTrack, SyntenyTrack, TanglegramTrack, Theme, Track,
    Tree, TreeTrack, VariantTrack, WindowStyle, WindowTrack,
};

use crate::cli::args::{
    Invocation, Kind, Palette, Place, ShadePlace, Source, Style, Threshold, TrackSpec, TreeSupport,
};
use crate::read;
use crate::track::traits::Traits;
use crate::Mutations;
use crate::TreeShape;

/// What went wrong once the command line was understood.
#[derive(Debug)]
pub enum BuildError {
    /// A file would not open.
    Open {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What the operating system said.
        cause: io::Error,
    },
    /// A file opened and did not say what it claimed to.
    Parse {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called, or `standard input`.
        path: String,
        /// Which line, and what was wrong with it.
        cause: read::ReadError,
    },
    /// A file held nothing the figure could use.
    Empty {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What was expected in it.
        wanted: &'static str,
    },
    /// A bigWig or a bigBed whose index names no sequence the window is on.
    ///
    /// kent's writers index only the sequences that hold data, so a file
    /// written against a whole genome's lengths names only those, and holds
    /// no rows on the rest, as the text it was written from holds none
    /// there. One place of several is drawn as a band saying so, as it is
    /// from that text, where the whole figure was refused; a figure of one
    /// place is refused naming the sequences the file has, with the
    /// `--rename` that would draw it where one plainly would.
    Absent {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What was expected in it.
        wanted: &'static str,
        /// What the file said, naming the sequences it has.
        said: String,
    },
    /// A file held what was wanted, somewhere the figure is not.
    ///
    /// Separate from [`BuildError::Empty`] because the two are different
    /// mistakes with the same symptom. An empty file is a wrong path; a full
    /// file with nothing in the window is a wrong locus, or a file whose first
    /// column names its sequence something the region does not.
    Elsewhere {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What was expected in it.
        wanted: &'static str,
        /// How many of them the file did hold.
        held: usize,
        /// The sequences it named, where naming them helps.
        named: String,
        /// The locus that was asked for.
        region: String,
        /// The `--rename` that would draw it, where one plainly would.
        rename: Option<Box<(String, String)>>,
    },
    /// A reference whose bases stop before the region starts, or start after
    /// it ends.
    ///
    /// The record was the right one and none of it is in the window, so the
    /// track would be drawn empty: letters, frames and mismatches all need a
    /// base to stand on. It was drawn empty, and the command exited nought.
    Beyond {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// The record, by the name the file gives it.
        record: String,
        /// The first base it holds, 1-based.
        first: u64,
        /// The last base it holds, 1-based.
        last: u64,
        /// The locus that was asked for.
        region: String,
    },
    /// A place named by a word that no file of the figure names.
    Nowhere {
        /// The word.
        name: String,
        /// The files of the figure that name sequences, in a header or on
        /// their rows, gathered by the sequences they name.
        held: Vec<Naming>,
        /// Genes named nearly the same, which may be what was meant.
        near: Vec<String>,
        /// Whether any file of the figure is an annotation to look genes up in.
        annotated: bool,
        /// The `--rename` that would find it, where one plainly would.
        rename: Option<Box<(String, String)>>,
    },
    /// A gene named at more than one place.
    Several {
        /// The name.
        name: String,
        /// Where each one is.
        places: Vec<String>,
    },
    /// A bigWig, a bigBed, a 2bit, a BCF or a `.hic`, which karyon reads as
    /// it is, handed to a track that does not draw what it holds, or given a
    /// `--format`, which says what the columns of a text file are. And a
    /// `.hic` named for a file a track reads as text, such as `--ld`'s, which
    /// no one command writes as text.
    ///
    /// Read as the text its tool writes, a bigWig handed to `--pileup` was
    /// refused with the command that turns it into bedGraph, which `--pileup`
    /// would have refused next.
    OtherTrack {
        /// Which track was handed it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What its first bytes say it is.
        binary: Binary,
        /// Whether `--format` was given, which is the fault, rather than the
        /// track.
        format: bool,
    },
    /// A file that is not text, with the command that writes the text.
    NotText {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// What its first bytes say it is.
        binary: Binary,
        /// The command to write in place of its name, where there is a name.
        instead: Option<String>,
        /// Whether the file was named on its own, its track read off its
        /// name, so that what goes in its place has to bring the flag: a
        /// `<(...)` or a `-` has no name to read a track off.
        alone: bool,
    },
    /// A `--resolution` the track's file cannot be drawn at: a file that is
    /// not a `.hic`, whose bins are the ones it was written in, or a `.hic`
    /// that does not hold that resolution, which is told the ones it does.
    Unresolved {
        /// Which track wanted it.
        track: &'static str,
        /// What the file was called.
        path: String,
        /// The resolution asked for, in bases.
        asked: u32,
        /// The resolutions the `.hic` holds, or `None` for a file that is not
        /// one.
        held: Option<Vec<u32>>,
    },
    /// A threshold given as a number no p-value can be, for a scan whose file
    /// held p-values.
    ///
    /// The threshold is in the file's units, and 7.3 was the right number for
    /// a file of `-log10(p)` and is no p-value at all.
    NotAPValue {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// The number given.
        given: f64,
    },
    /// A file holds several of a thing and the command asked for none of them.
    ///
    /// The `--format` case turned round: there the shape is ambiguous and the
    /// file cannot say, and here the file says several things and only one of
    /// them is a track. Picking the first would draw one of them under a label
    /// that names none.
    Ambiguous {
        /// Which track wanted it.
        track: &'static str,
        /// What it was called.
        path: String,
        /// The flag that settles it.
        flag: &'static str,
        /// What the file holds, so the choice can be made without opening it.
        choices: Vec<String>,
    },
    /// Two files were read and nothing in one names anything in the other.
    ///
    /// The join is names, and names from two tools are routinely not the same
    /// strings. A figure drawn from a join that found nothing is not blank: it
    /// is every gene marked as having no counterpart, or a phylogeny with no
    /// block on it, and both of those read as a finding.
    Unjoined {
        /// Which track wanted it.
        track: &'static str,
        /// The file whose names found nothing.
        path: String,
        /// What kind of name did not join.
        what: &'static str,
        /// What it was matched against.
        against: &'static str,
        /// A few of the names, so the mismatch can be seen at a glance.
        examples: Vec<String>,
    },
    /// Something was asked for by name and the file has nothing of that name:
    /// a column of a sheet, a clade, tip, change or annotation of a tree, a
    /// record of a FASTA, or the row `--compare-to` names.
    ///
    /// Not [`BuildError::Ambiguous`], which is a file holding several things
    /// and a command naming none of them. Here the command named one and the
    /// file has not got it, which is nearly always a spelling, so the names
    /// it does have are worth printing beside it.
    Unnamed {
        /// Which track wanted it.
        track: &'static str,
        /// What the file was called.
        path: String,
        /// What kind of thing was being named, as a reader would say it.
        what: &'static str,
        /// The one that is not in the file.
        wanted: String,
        /// The ones that are.
        held: Vec<String>,
    },
    /// A name meant to pick one thing out of a file that holds several of it.
    ///
    /// Picking the first would be a figure drawn against something the caller
    /// did not choose, and looking like the one they asked for.
    Repeated {
        /// Which track wanted it.
        track: &'static str,
        /// What the file was called.
        path: String,
        /// What kind of thing was being named.
        what: &'static str,
        /// The name that picks more than one.
        name: String,
        /// How many it picks.
        held: usize,
    },
    /// A copy number track with no ploidy to read its levels against.
    ///
    /// The parser refuses this, so it reaches here only from an [`Invocation`]
    /// built by hand, whose fields are all public. Defaulting it would draw a
    /// confident ladder whose rule came from nowhere, and the rule is what
    /// separates a gain from a loss.
    MissingPloidy {
        /// Which track is short of it.
        track: &'static str,
    },
    /// A track drawn from two files was handed one.
    ///
    /// The parser refuses this, so it reaches here only from an
    /// [`Invocation`] built by hand, whose fields are all public.
    MissingSecond {
        /// Which track is short of a file.
        track: &'static str,
    },
    /// A `--colors` naming a column no `--traits` sheet of the figure has,
    /// or a value no row of one holds in that column.
    ///
    /// Refused rather than passed over, as a misspelt `--columns` is: the
    /// colours would paint nothing, and the figure would come out in the
    /// palette looking as though they had been read.
    NotColored {
        /// The column `--colors` names.
        column: String,
        /// The value it names, or `None` where no sheet has the column.
        value: Option<String>,
        /// The sheets of the figure, by what they were called.
        sheets: Vec<String>,
        /// What they hold instead: their columns, or the values of this one.
        held: Vec<String>,
    },
    /// A `--colors` for a column of a sheet that no track draws: left out
    /// of every `--columns`, and colouring no tree's branches.
    ///
    /// The colours would reach it and show nowhere, which is a figure that
    /// looks as though they had been read and were wrong.
    ColorsUndrawn {
        /// The column.
        column: String,
    },
    /// A `--colors` for a column whose every value is a number.
    ///
    /// A strip draws such a column on a ramp, and so do the branches
    /// coloured by it, and a colour chosen for a value reaches neither.
    ColorsOfNumbers {
        /// The column.
        column: String,
    },
    /// The tree would not parse.
    Tree {
        /// The flag that asked for it, since more than one takes a phylogeny
        /// and a tanglegram takes two by different names.
        flag: &'static str,
        /// What it was called, or `standard input`.
        path: String,
        /// Why it is not a tree.
        cause: crate::Error,
    },
    /// A `--shade` the figure has nowhere to draw.
    ///
    /// Refused rather than left out, as a file with nothing in the window is:
    /// a figure drawn without the stretch it was asked to mark looks like one
    /// where nothing is there to mark. Only a stretch on the right sequence
    /// and outside the window is let go with a note, since a page that moves
    /// the window runs the same command again at every step.
    Unshaded {
        /// The value as it was written.
        given: String,
        /// Why there is nowhere to draw it.
        why: ShadeRefusal,
    },
    /// A track drawn over a place, in a figure that names none, found once
    /// its file is read: a `.bed` named on its own that holds features, where
    /// mosdepth's depth in windows, a bedGraph by another name, would have
    /// been drawn across the genome. Said as the parser says it of a track
    /// whose name tells it, so `karyon genes.bed` is answered as it always was.
    Placeless(crate::cli::args::ArgError),
    /// A `--codons` with no one unbroken coding sequence to count.
    ///
    /// Refused rather than drawn over whatever came nearest: a ruler that
    /// numbers the wrong stretch names the wrong residue at every codon, and
    /// looks exactly as right as one that numbers the right stretch.
    Uncounted(CodonRefusal),
    /// A `--circular` with no whole sequence to draw round, or a track that
    /// turned out, once read, to have no ring.
    Uncircled(CircleRefusal),
}

/// Why a `--circular` has no circle to draw.
#[derive(Debug, Clone, PartialEq)]
pub enum CircleRefusal {
    /// A circle handed to a builder of figures along a sequence, which
    /// [`build_circle`] draws instead.
    NotAFigure,
    /// The place is a gene, which is part of a sequence.
    Gene {
        /// The gene, as its annotation spells it.
        name: String,
        /// The sequence it is on.
        sequence: String,
    },
    /// No file says how long the sequence is, so nothing says where the
    /// circle closes.
    ///
    /// A figure along a sequence is drawn as far as its rows reach and says
    /// so. A circle cannot be: a ring that closes where the last row of a
    /// bedGraph happens to end puts the origin's neighbour in the wrong
    /// place, and every angle round it with it.
    NoLength {
        /// The sequence.
        sequence: String,
        /// The file whose rows reach furthest on it, and how far.
        reached: String,
        /// How far they reach.
        reach: u64,
    },
    /// A span written from base 1 that ends where no file says the sequence
    /// does.
    Length {
        /// The sequence.
        sequence: String,
        /// The end written.
        written: u64,
        /// The length the files state.
        stated: u64,
    },
    /// A sequence the files give different lengths.
    ///
    /// A figure along it takes the first file's word. A circle closed at
    /// either length draws the other file's rows round the wrong angles, and
    /// which one it closed at was the order the files were written in.
    Lengths {
        /// The sequence.
        sequence: String,
        /// Each file that states a length, and the length, in the order of
        /// the tracks.
        said: Vec<(String, u64)>,
    },
    /// A file named on its own that turned out, once read, to be a kind with
    /// no ring: a `.bed` that is modkit's bedMethyl.
    NoRing {
        /// The file.
        path: String,
        /// What it turned out to be.
        kind: Kind,
    },
}

impl fmt::Display for CircleRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::track::axis::group_thousands;
        match self {
            CircleRefusal::NotAFigure => write!(
                f,
                "--circular draws a circle, which is no figure along a sequence: draw it \
                 with build_circle"
            ),
            CircleRefusal::Gene { name, sequence } => write!(
                f,
                "{name} is a gene, and a circle is a whole sequence: name the sequence it \
                 is on, as karyon {sequence} --circular"
            ),
            CircleRefusal::NoLength {
                sequence,
                reached,
                reach,
            } => write!(
                f,
                "a circle closes where {sequence} ends, and no file says where that is, \
                 only that {reached} reaches {}: add its FASTA, a BAM, a VCF with \
                 ##contig or a GFF3 with ##sequence-region, or write the span from 1, as \
                 {sequence}:1-LENGTH",
                group_thousands(*reach)
            ),
            CircleRefusal::Length {
                sequence,
                written,
                stated,
            } => write!(
                f,
                "{sequence}:1-{} closes the circle at {}, and the files say {sequence} is \
                 {} bases long: write {sequence} alone to draw all of it",
                group_thousands(*written),
                group_thousands(*written),
                group_thousands(*stated)
            ),
            CircleRefusal::Lengths { sequence, said } => {
                let said: Vec<String> = said
                    .iter()
                    .map(|(file, length)| format!("{file} says {} bases", group_thousands(*length)))
                    .collect();
                write!(
                    f,
                    "a circle closes where {sequence} ends, and the files disagree on where \
                     that is: {}; draw it from files that agree on how long {sequence} is",
                    joined(&said)
                )
            }
            CircleRefusal::NoRing { path, kind } => write!(
                f,
                "{path} holds what {} draws, once read, and --circular draws {} tracks as \
                 rings: draw it along the sequence, without --circular",
                kind.dashed(),
                crate::cli::args::on_a_circle()
            ),
        }
    }
}

/// Why a `--codons` has no coding sequence to count.
#[derive(Debug, Clone, PartialEq)]
pub enum CodonRefusal {
    /// The figure has no annotation to find a gene in.
    NoAnnotation,
    /// The figure is placed on a sequence drawn whole, which is no gene.
    WholeSequence {
        /// The sequence.
        sequence: String,
    },
    /// The annotation writes no CDS for the gene the figure is placed on.
    NoCds {
        /// The gene, as the annotation spells it.
        gene: String,
        /// The annotations, by what they were called.
        files: Vec<String>,
    },
    /// No gene of the annotation codes in the place written.
    NoneCodes {
        /// The place, as it was written.
        place: String,
        /// The annotations, by what they were called.
        files: Vec<String>,
    },
    /// More than one gene codes in the place.
    Several {
        /// The genes.
        genes: Vec<String>,
        /// The place, as it was written or named.
        place: String,
    },
    /// The transcripts of one gene code different stretches.
    Isoforms {
        /// The gene.
        gene: String,
        /// Its transcripts, by their names.
        transcripts: Vec<String>,
    },
    /// The coding sequence is in pieces with introns between them.
    Spliced {
        /// The gene.
        gene: String,
        /// How many pieces.
        pieces: usize,
    },
    /// The gene is on no strand, so neither end is where it starts.
    NoStrand {
        /// The gene.
        gene: String,
    },
    /// The rows of the coding sequence leave the frame they began in, as a
    /// ribosomal slippage is written: two rows that share bases, or that
    /// meet with a phase that does not carry the frame on.
    Frameshift {
        /// The gene.
        gene: String,
        /// The first base read in the new frame, 0-based.
        at: u64,
    },
    /// The coding sequence does not begin on its start codon: its 5'-most
    /// row has a phase of 1 or 2, or says the CDS goes on past it, as NCBI
    /// writes a CDS at the edge of a contig.
    Partial {
        /// The gene.
        gene: String,
        /// The phase of that row, where that is what says so.
        phase: Option<u8>,
    },
    /// The CDS names a translation table NCBI does not list.
    UnknownTable {
        /// The gene.
        gene: String,
        /// The number its CDS gives.
        table: u8,
    },
}

/// Why a `--shade` has nowhere to go.
#[derive(Debug, Clone, PartialEq)]
pub enum ShadeRefusal {
    /// Nothing in the figure is laid on the coordinates.
    NothingToShade,
    /// The stretch is on a sequence no place of the figure is on.
    Elsewhere {
        /// The sequence it is on.
        sequence: String,
        /// The places the figure is drawn over.
        places: Vec<String>,
    },
    /// A gene to shade, and no annotation to look it up in.
    NoAnnotation,
    /// A gene to shade that the annotation does not name.
    NoSuchGene {
        /// Names it does give that are nearly the same.
        near: Vec<String>,
    },
    /// A span with no sequence, on a figure whose axis is every sequence.
    WholeGenome,
    /// A gene, on a figure across the whole genome, which reads no
    /// annotation to look one up in: a GFF3 beside it needs a place.
    GenomeGene,
}

/// The kind of system a message that names a command is written for.
///
/// It is the program's target and not the shell that decides, because what
/// differs is what the program can open. `<(command)` hands the command's
/// output over as a path to a pipe, `/dev/fd/63`, which a program built for
/// Linux or macOS opens like a file and a program built for Windows cannot,
/// whether cmd, PowerShell or Git Bash started it. A pipe into standard input
/// works in every one of those shells, and `-` is where karyon reads it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shell {
    /// bash or zsh on Linux, macOS or WSL, which runs the Linux build.
    Posix,
    /// A Windows build, in any Windows shell.
    Windows,
}

/// The system this build of the program runs on.
///
/// `cfg!` and not `#[cfg]`, so that both wordings are compiled, linted and
/// tested on every system rather than one of them only where it ships. A
/// page in a browser is wasm32, not Windows, and keeps the `<(...)` wording.
const HOST: Shell = if cfg!(windows) {
    Shell::Windows
} else {
    Shell::Posix
};

/// What to write so a command's output is read in place of a file, as a
/// message says it to `shell`, with the track's flag where the file was named
/// on its own.
///
/// Every message that offers a command in place of a file's name says it
/// through here, so none can offer a Windows user a `<(...)` that their
/// shell cannot hand over. A file named on its own was given its track by its
/// name, and what replaces it has none to go by: a bare `<(...)` is read as a
/// place called `/dev/fd/63`, and a bare `-` is refused for want of a track.
/// So `flag`, the track the file was given, is written in front of either.
fn in_place(shell: Shell, command: &str, flag: Option<&str>) -> String {
    let flag = flag.map(|track| format!("--{track} ")).unwrap_or_default();
    match shell {
        Shell::Posix => format!("write {flag}<({command}) where its name is"),
        Shell::Windows => {
            format!("pipe what {command} writes into karyon, with {flag}- where its name is")
        }
    }
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::Open { track, path, cause } => write!(f, "--{track} {path}: {cause}"),
            BuildError::Parse { track, path, cause } => write!(f, "--{track} {path}: {cause}"),
            BuildError::Empty {
                track,
                path,
                wanted,
            } => write!(f, "--{track} {path}: no {wanted} in the region"),
            BuildError::Absent {
                track, path, said, ..
            } => write!(f, "--{track} {path}: {said}"),
            BuildError::Elsewhere {
                track,
                path,
                wanted,
                held,
                named,
                region,
                rename,
            } => {
                write!(f, "--{track} {path}: no {wanted} in {region}")?;
                write!(f, ", though the file holds {held}")?;
                if !named.is_empty() {
                    write!(f, " on {named}")?;
                }
                if let Some(pair) = rename {
                    let (from, to) = pair.as_ref();
                    write!(f, "; if {from} is {to}, add --rename {from}={to}")?;
                }
                Ok(())
            }
            BuildError::Nowhere {
                name,
                held,
                near,
                annotated,
                rename,
            } => {
                write!(f, "no gene and no sequence is called {name}")?;
                match near.as_slice() {
                    [] => write!(f, ".")?,
                    [one] => write!(f, "; did you mean {one}?")?,
                    [many @ .., last] => {
                        write!(f, "; did you mean {} or {last}?", many.join(", "))?
                    }
                }
                // What each file calls its sequences, since a name that is
                // nowhere is most often one file's name for what another
                // calls otherwise, as PLINK writes 1 for NC_000962.3.
                for (files, sequences) in held {
                    let files: Vec<(String, usize)> =
                        files.iter().map(|file| (file.clone(), 0)).collect();
                    let verb = if files.len() == 1 { "names" } else { "name" };
                    let noun = if sequences.len() == 1 {
                        "sequence"
                    } else {
                        "sequences"
                    };
                    write!(
                        f,
                        " {} {verb} the {noun} {}.",
                        listed_as(&files, "files"),
                        listed(sequences)
                    )?;
                }
                if let Some(pair) = rename {
                    let (from, to) = pair.as_ref();
                    write!(f, " If {from} is {to}, add --rename {from}={to}.")?;
                }
                if near.is_empty() && !annotated {
                    write!(
                        f,
                        " To find a gene by name, add the annotation that names it, as --features genes.gff3."
                    )?;
                }
                Ok(())
            }
            BuildError::Several { name, places } => write!(
                f,
                "{name} is named at {} places, {}; give the one you mean as a region",
                places.len(),
                places.join(" and ")
            ),
            BuildError::NotText {
                track,
                path,
                binary,
                instead,
                alone,
            } => match instead {
                Some(command) => write!(
                    f,
                    "--{track} {path}: {binary}; {}, or turn it into text first",
                    in_place(HOST, command, alone.then_some(*track))
                ),
                // A BCF is read from its bytes, which a pipe of text gives
                // none of, and the text bcftools writes is what a pipe takes.
                None if *binary == Binary::Bcf => write!(
                    f,
                    "--{track} {path}: the file is BCF, which karyon reads from a file \
                     named on the command line and not from a pipe; name the file instead, \
                     or pipe in the text bcftools view writes"
                ),
                // A pipe is read once, front to back, and these are read by
                // going to where their index says.
                None if binary.drawn_by().is_some() => write!(
                    f,
                    "--{track} {path}: the file is {}, which is read through the index it \
                     holds, and a pipe cannot be read that way; name the file instead",
                    binary.called()
                ),
                None => match binary.advice() {
                    Some(advice) => write!(f, "--{track} {path}: {binary}; {advice}"),
                    None => write!(
                        f,
                        "--{track} {path}: {binary}; pipe it through the tool that writes \
                         it as text first"
                    ),
                },
            },
            BuildError::Placeless(said) => write!(f, "{said}"),
            BuildError::Uncounted(why) => write!(f, "{why}"),
            BuildError::Uncircled(why) => write!(f, "{why}"),
            BuildError::OtherTrack {
                track,
                path,
                binary,
                format,
            } => {
                let (holds, drawn_by) = binary.drawn_by().unwrap_or(("text", "a track"));
                if *format {
                    write!(
                        f,
                        "--{track} {path}: --format says what the columns of a text file are, \
                         and the file is {}, which says what it holds itself",
                        binary.called()
                    )
                } else {
                    write!(
                        f,
                        "--{track} {path}: the file is {}, {holds}, which {drawn_by}",
                        binary.called()
                    )
                }
            }
            BuildError::Unresolved {
                track,
                path,
                asked,
                held,
            } => match held {
                None => write!(
                    f,
                    "--{track} {path}: --resolution picks one of the resolutions a .hic holds, \
                     and {path} is drawn at the bins it was written in"
                ),
                Some(held) => write!(
                    f,
                    "--{track} {path} has no {}-base resolution; it holds {}",
                    crate::track::axis::group_thousands(u64::from(*asked)),
                    read::hic::listed(held)
                ),
            },
            BuildError::NotAPValue { track, path, given } => write!(
                f,
                "--{track} {path} holds p-values, so --threshold is a p-value too, between 0 \
                 and 1, and {given} is not one; give it as 5e-8, or as genome-wide"
            ),
            BuildError::Beyond {
                track,
                path,
                record,
                first,
                last,
                region,
            } => write!(
                f,
                "--{track} {path}: {record} holds bases {} to {}, and none of them is in {region}",
                crate::track::axis::group_thousands(*first),
                crate::track::axis::group_thousands(*last),
            ),
            BuildError::Ambiguous {
                track,
                path,
                flag,
                choices,
            } => write!(
                f,
                "--{track} {path} holds {}, and {flag} says which to draw",
                choices.join(", ")
            ),
            BuildError::Unjoined {
                track,
                path,
                what,
                against,
                examples,
            } => {
                write!(
                    f,
                    "--{track} {path}: no {what} in this file names anything in {against}"
                )?;
                if !examples.is_empty() {
                    write!(f, ", starting with {}", examples.join(", "))?;
                }
                Ok(())
            }
            BuildError::Unnamed {
                track,
                path,
                what,
                wanted,
                held,
            } => {
                write!(f, "--{track} {path} has no {what} called {wanted}; it has ")?;
                // A plain Newick carries no annotation at all, which is the
                // ordinary way to ask a tree for a colour key it has not got,
                // and a list of nothing would leave the sentence unfinished.
                if held.is_empty() {
                    write!(f, "none")
                } else {
                    write!(f, "{}", held.join(", "))
                }
            }
            BuildError::Repeated {
                track,
                path,
                what,
                name,
                held,
            } => write!(
                f,
                "--{track} {path} has {held} {what}s called {name}, so the name does not pick one"
            ),
            BuildError::MissingPloidy { track } => write!(
                f,
                "a {track} track is drawn against a ploidy, and none was given"
            ),
            BuildError::MissingSecond { track } => write!(
                f,
                "a {track} track is drawn from two files, and only one was given"
            ),
            BuildError::NotColored {
                column,
                value,
                sheets,
                held,
            } => {
                let several = sheets.len() > 1;
                let named: Vec<(String, usize)> =
                    sheets.iter().map(|sheet| (sheet.clone(), 0)).collect();
                let sheets = listed_as(&named, "sheets");
                match value {
                    None => write!(
                        f,
                        "--colors names a column called {column}, and {sheets} {} none; {} {}",
                        if several { "have" } else { "has" },
                        if several { "they have" } else { "it has" },
                        held.join(", ")
                    ),
                    // A column of nothing but gaps, or a sheet of a header
                    // alone, holds no value to list, and the sentence was
                    // left unfinished at "holds".
                    Some(value) if held.is_empty() => write!(
                        f,
                        "--colors names {value} in {column}, and no row of {sheets} holds a \
                         value in {column}"
                    ),
                    Some(value) => write!(
                        f,
                        "--colors names {value} in {column}, and no row of {sheets} holds it; \
                         {column} holds {}",
                        held.join(", ")
                    ),
                }
            }
            BuildError::ColorsUndrawn { column } => write!(
                f,
                "--colors paints {column}, and no track draws it: name it in --columns, or \
                 colour a tree's branches by it with --color-by {column}"
            ),
            BuildError::ColorsOfNumbers { column } => write!(
                f,
                "--colors names {column}, a column of numbers, which is drawn as a ramp; \
                 --colors paints a column of words"
            ),
            BuildError::Tree { flag, path, cause } => write!(f, "{flag} {path}: {cause}"),
            BuildError::Unshaded { given, why } => match why {
                ShadeRefusal::NothingToShade => write!(
                    f,
                    "--shade {given}: nothing in this figure is laid on the coordinates to \
                     shade; a phylogeny's x is a branch length, and an ideogram's the whole \
                     chromosome"
                ),
                ShadeRefusal::Elsewhere { sequence, places } => write!(
                    f,
                    "--shade {given} is on {sequence}, and the figure is drawn over {}",
                    joined(places)
                ),
                ShadeRefusal::NoAnnotation => write!(
                    f,
                    "--shade {given}: a gene is looked up in the figure's annotation, and it \
                     has none; add the GFF3, GTF or BED, or write the span"
                ),
                ShadeRefusal::NoSuchGene { near } => match near.as_slice() {
                    [] => write!(f, "--shade {given}: no gene of that name in the annotation"),
                    _ => write!(
                        f,
                        "--shade {given}: no gene of that name; did you mean {}?",
                        joined_or(near)
                    ),
                },
                ShadeRefusal::WholeGenome => write!(
                    f,
                    "--shade {given}: a figure across the whole genome is shaded on one of \
                     its sequences, as 7:1,001-2,000"
                ),
                ShadeRefusal::GenomeGene => write!(
                    f,
                    "--shade {given}: a figure across the whole genome reads no annotation \
                     to look a gene up in; write the gene's span on its sequence, as \
                     7:1,001-2,000"
                ),
            },
        }
    }
}

/// Places, as a sentence lists them: `a`, `a and b`, `a, b and c`.
fn joined(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [many @ .., last] => format!("{} and {last}", many.join(", ")),
    }
}

/// Names, as a question offers them: `a`, `a or b`, `a, b or c`.
fn joined_or(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [many @ .., last] => format!("{} or {last}", many.join(", ")),
    }
}

impl std::error::Error for BuildError {}

impl fmt::Display for CodonRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Five names and how many more, since a place as wide as a
        // chromosome holds thousands of genes and a refusal is one line.
        let names = |names: &[String], things: &str| {
            listed_as(
                &names
                    .iter()
                    .map(|name| (name.clone(), 0))
                    .collect::<Vec<_>>(),
                things,
            )
        };
        let files = |files: &[String]| names(files, "files");
        match self {
            CodonRefusal::NoAnnotation => write!(
                f,
                "--codons counts the codons of a gene in the figure's annotation, and the \
                 figure has none: add the GFF3, GTF or BED that writes its CDS"
            ),
            CodonRefusal::WholeSequence { sequence } => write!(
                f,
                "--codons counts the codons of one gene, and {sequence} is a whole sequence: \
                 place the figure on a gene by its name, or on a stretch one gene codes in"
            ),
            CodonRefusal::NoCds { gene, files: held } => write!(
                f,
                "{gene} has no CDS in {}, and --codons counts from the first base of a CDS, \
                 which a gene's own span does not say: add the annotation's CDS rows, or a \
                 BED12 whose thick span is the part that codes",
                files(held)
            ),
            CodonRefusal::NoneCodes { place, files: held } => write!(
                f,
                "no gene of {} codes in {place}, and --codons counts the codons of one",
                files(held)
            ),
            CodonRefusal::Several { genes, place } => write!(
                f,
                "{} code in {place}, and --codons counts the codons of one gene: place the \
                 figure on one of them by its name, as karyon {} in place of {place}",
                names(genes, "genes"),
                genes.first().map_or("GENE", String::as_str)
            ),
            CodonRefusal::Isoforms { gene, transcripts } => write!(
                f,
                "the transcripts of {gene}, {}, code different stretches, and --codons counts \
                 one: place the figure on a transcript by its name, as karyon {} in place of \
                 {gene}",
                names(transcripts, "transcripts"),
                transcripts.first().map_or("TRANSCRIPT", String::as_str)
            ),
            CodonRefusal::Spliced { gene, pieces } => write!(
                f,
                "{gene} codes in {pieces} pieces with introns between them, and --codons \
                 counts one unbroken coding sequence: across the introns it would number \
                 bases that are never translated"
            ),
            CodonRefusal::NoStrand { gene } => write!(
                f,
                "{gene} is on no strand, and --codons counts from its start codon, which the \
                 strand puts at one end or the other: write + or - in the strand column of \
                 its rows"
            ),
            CodonRefusal::Frameshift { gene, at } => write!(
                f,
                "the CDS of {gene} changes frame at {}, where its rows overlap or meet out \
                 of frame, as a ribosomal slippage is written, and --codons counts codons in \
                 one frame: past {0} it would number every codon in the frame the ribosome \
                 left",
                crate::track::axis::group_thousands(at.saturating_add(1))
            ),
            CodonRefusal::Partial {
                gene,
                phase: Some(phase),
            } => write!(
                f,
                "the CDS of {gene} has phase {phase}, so its first {} {} a codon that begins \
                 before it and the annotation does not hold its start codon, and --codons \
                 counts from the start codon: every codon would be numbered from one that is \
                 not the first",
                if *phase == 1 { "base" } else { "2 bases" },
                if *phase == 1 { "ends" } else { "end" },
            ),
            CodonRefusal::Partial { gene, phase: None } => write!(
                f,
                "the CDS of {gene} is written as partial at its 5' end, so the annotation \
                 does not hold its start codon, and --codons counts from the start codon: \
                 every codon would be numbered from one that is not the first"
            ),
            CodonRefusal::UnknownTable { gene, table } => write!(
                f,
                "the CDS of {gene} names translation table {table}, which NCBI does not list; \
                 --genetic-code N says which table to read it with"
            ),
        }
    }
}

/// Builds the figure the command line asked for and renders it.
///
/// # Errors
///
/// Returns the first file that would not open, would not parse, or held
/// nothing inside the region. The command line itself has already been checked
/// by [`crate::cli::args::parse`], so everything here is about the data.
pub fn build(
    invocation: &Invocation,
    open: impl FnMut(&Source) -> io::Result<String>,
) -> Result<String, BuildError> {
    build_with(invocation, open, |_, _| None)
}

/// The same, for a caller that has already read one of these trees.
///
/// `parsed` is offered the name and the text of every phylogeny before it is
/// read, and may answer with one it made earlier. A shell never can: it runs
/// once and reads each file once. A page driving the program on every move
/// does, and the difference is most of the work: reading a million tip tree
/// takes 361 ms of a 578 ms figure in a browser, where folding it and drawing
/// sixty rows of it take 189 between them.
///
/// # Errors
///
/// The same as [`build`].
pub fn build_with(
    invocation: &Invocation,
    mut open: impl FnMut(&Source) -> io::Result<String>,
    parsed: impl FnMut(&str, &str) -> Option<Tree>,
) -> Result<String, BuildError> {
    build_files(invocation, &mut open, parsed)
}

/// The same, reading through [`Files`], which may do more than hand over
/// text: [`Disk`] reads a BAM a window at a time.
///
/// # Errors
///
/// The same as [`build`].
pub fn build_files(
    invocation: &Invocation,
    files: &mut dyn Files,
    parsed: impl FnMut(&str, &str) -> Option<Tree>,
) -> Result<String, BuildError> {
    let mut theme = match invocation.theme {
        Palette::Dark => Theme::dark(),
        Palette::Light => Theme::light(),
    };
    if let Some(ground) = &invocation.background {
        theme.background = ground.clone();
    }
    if invocation.circular {
        return build_circle(invocation, files, theme).map(|circle| circle.to_svg());
    }
    if invocation.more.is_empty() {
        return build_figure(invocation, files, parsed, theme, None)
            .map(|built| built.figure.to_svg());
    }
    build_sheet(invocation, files, parsed, theme).map(|sheet| sheet.to_svg())
}

/// A figure of several places, as `karyon rpoB katG inhA reads.bam
/// genes.gff3`: one panel a place, each the same tracks over its own place,
/// one under the other with their plotting areas aligned, the title over
/// them all and the key once under them.
///
/// Each panel is the figure its place would be on its own, so a gene is
/// titled with its name and a locus says itself at the top right, and what
/// is piped in is read once for all of them.
///
/// # Errors
///
/// The first place that would not draw, as [`build`] says it.
pub fn build_sheet(
    invocation: &Invocation,
    files: &mut dyn Files,
    mut parsed: impl FnMut(&str, &str) -> Option<Tree>,
    theme: Theme,
) -> Result<crate::Panels, BuildError> {
    if invocation.circular {
        return Err(BuildError::Uncircled(CircleRefusal::NotAFigure));
    }
    let mut kept = KeptStdin { files, stdin: None };
    let first = match (&invocation.region, &invocation.named) {
        (Some(region), _) => Some(Place::Locus(region.clone())),
        (None, Some(name)) => Some(Place::Named(name.clone())),
        (None, None) => None,
    };
    let places: Vec<Place> = first
        .into_iter()
        .chain(invocation.more.iter().cloned())
        .collect();
    let mut sheet = crate::Panels::new().theme(theme.clone());
    if let Some(title) = &invocation.title {
        sheet = sheet.title(title);
    }
    let mut legend = crate::track::legend::Legend::new();
    // Each shade is settled across every panel, once they are all drawn.
    let mut shading = Shading::new(invocation);
    let mut names = Vec::with_capacity(places.len());
    let mut figures = Vec::with_capacity(places.len());
    // The width that letters the bases of every panel, said once under them
    // all. Said a panel at a time, each panel named the width its own span
    // needed, and a reader who took the first width still had the panels of
    // longer places in blocks, which said so again.
    let mut letters = None;
    for place in &places {
        let mut one = invocation.clone();
        one.more = Vec::new();
        one.title = None;
        one.legend = false;
        match place {
            Place::Locus(region) => {
                one.region = Some(region.clone());
                one.named = None;
                names.push(region.to_string());
            }
            Place::Named(name) => {
                one.region = None;
                one.named = Some(name.clone());
                names.push(name.clone());
            }
        }
        let (built, wanted) = build_one(
            &one,
            &mut kept,
            &mut parsed,
            theme.clone(),
            None,
            Some(&mut shading),
        )?;
        letters = letters.max(wanted);
        gather(&mut legend, &built.legend);
        figures.push(built.figure);
    }
    settle_shades(invocation, &shading.fates, &mut kept)?;
    note_letters(&mut kept, letters);
    // One scale across the panels as well as down each: the depth over rpoB
    // and the depth over katG read off one ceiling, or the eye compares two.
    if invocation.same_scale {
        let extents = crate::Extent::join(figures.iter().flat_map(Figure::extents));
        figures = figures
            .into_iter()
            .map(|figure| figure.share_extents(extents.clone()))
            .collect();
    }
    for figure in &figures {
        sheet = sheet.push_bare(figure);
    }
    if invocation.legend && !legend.is_empty() {
        let key = crate::Figure::new(Region::new("key", 0, 1).expect("a one-base window"))
            .width(invocation.width.unwrap_or(900.0))
            .theme(theme)
            .show_region_label(false)
            .push(crate::track::legend::LegendTrack::new(legend));
        sheet = sheet.push_bare(&key);
    }
    Ok(sheet.description(format!(
        "A karyon figure of {} places, one panel each over the same tracks: {}.",
        places.len(),
        names.join(", ")
    )))
}

/// A figure a command line builds, before it is written out.
pub struct Built {
    /// The figure, for [`Figure::to_svg`], or for
    /// [`Figure::to_svg_with_id_prefix`] where a page holds several.
    pub figure: Figure,
    /// The stretch of a genome it is drawn over, where it is drawn over one:
    /// the place the command line wrote or named, which [`build_figure`] can
    /// draw the same command over another stretch of. `None` for a figure
    /// that is its own place, as an alignment's columns or a table's weeks,
    /// and for one with no place, as a tree.
    pub along: Option<Region>,
    /// The key to the colours the tracks paint, drawn under the figure where
    /// the command line asks for one, and kept here either way, for a sheet of
    /// several places that draws it once under them all.
    pub legend: crate::track::legend::Legend,
}

/// What [`build_files`] draws, as the figure rather than its text: in
/// `theme`, and over `window` where one is given, in place of the place the
/// command line wrote or named.
///
/// A page that runs a command line draws it in the page's own colours, and
/// moves it along the genome under its reader's hand, which is the same
/// command over another stretch. A figure placed by a gene's name keeps the
/// gene as its title wherever it is moved to.
///
/// # Errors
///
/// The same as [`build`].
pub fn build_figure(
    invocation: &Invocation,
    files: &mut dyn Files,
    parsed: impl FnMut(&str, &str) -> Option<Tree>,
    theme: Theme,
    window: Option<&Region>,
) -> Result<Built, BuildError> {
    // A circle is not a figure, and drawn as one it would be the same tracks
    // along the sequence, which is not what was asked for and looks right.
    if invocation.circular {
        return Err(BuildError::Uncircled(CircleRefusal::NotAFigure));
    }
    let (built, letters) = build_one(invocation, files, parsed, theme, window, None)?;
    note_letters(files, letters);
    Ok(built)
}

/// Says how wide a figure would have to be for its bases to be letters, where
/// they are blocks of colour too narrow for them and a figure that wide is one
/// worth drawing.
fn note_letters(files: &mut dyn Files, letters: Option<u64>) {
    if let Some(width) = letters.filter(|width| *width <= 100_000) {
        files.note(&format!(
            "the bases are blocks of colour at this width, too narrow for their \
             letters; --width {width} draws the letters"
        ));
    }
}

/// The arcs a signal is cut into on a circle.
///
/// About one for every two pixels round the outside of a circle the default
/// size, so no arc is narrower than the eye follows, and a ring is at most a
/// thousand sectors however long the sequence is or however many lines its
/// file has: `samtools depth` over a chromosome of four megabases is four
/// million of them.
const ARCS: usize = 1_000;

/// The whitespace round a circle, which `--width` takes off its side.
const CIRCLE_MARGIN: f64 = 14.0;

/// The most named features a ring writes the names of. Past this a ring of
/// annotation drawn with names is a wheel of unreadable text, which is why
/// [`crate::FeatureRing::show_names`] is off by default.
const NAMED_ON_A_RING: usize = 20;

/// Draws a command line given `--circular`: its place, one whole sequence,
/// as a circle, each track a ring of it in the order written, the first
/// outermost, inside the ruler unless `--axis` puts it elsewhere or
/// `--no-axis` leaves it out. A breakend join on the sequence is a chord
/// across the middle, and the key under the circle names each ring, outside
/// in, and what its colours mean.
///
/// The parser has refused every track with no ring and every option a ring
/// would leave unsaid. What is left to refuse is about the files: a place
/// that turns out to be a gene, and a sequence no file gives the length of,
/// since a circle has to close where its sequence ends.
///
/// # Errors
///
/// The same as [`build`], and [`BuildError::Uncircled`] for those two.
pub fn build_circle(
    invocation: &Invocation,
    files: &mut dyn Files,
    theme: Theme,
) -> Result<crate::Rings, BuildError> {
    let mut kept = KeptStdin { files, stdin: None };
    let files: &mut dyn Files = &mut kept;
    let region = whole_sequence(invocation, files)?;
    let length = region.end();

    let mut rings = crate::Rings::new(length)
        .theme(theme.clone())
        .margin(CIRCLE_MARGIN);
    // `--width` is the side of the square, as it is the width of a figure.
    if let Some(width) = invocation.width {
        rings = rings.diameter(width - 2.0 * CIRCLE_MARGIN);
    }
    // The middle says which molecule this is and how long, which is what a
    // figure's locus says at its top right.
    let sequence = region.seq();
    let bases = format!("{} bases", crate::track::axis::group_thousands(length));
    rings = match (&invocation.title, invocation.region_label) {
        (Some(title), true) => rings.title(title).subtitle(format!("{sequence}, {bases}")),
        (Some(title), false) => rings.title(title),
        (None, true) => rings.title(sequence).subtitle(bases),
        (None, false) => rings,
    };

    for ringed in circle_rings(invocation, &region, &theme, files)? {
        if invocation.legend {
            if let Some(name) = ringed.name.filter(|name| !name.is_empty()) {
                rings = rings.key(name, ringed.legend);
            }
        }
        rings = match ringed.ring {
            RingOf::Signal(ring) => rings.push(ring),
            RingOf::Other(ring) => rings.push_boxed(ring),
        };
        // Nearly opaque, where a link between two spans is a wash: a join
        // between two breakends is a ribbon a pixel and a half wide at
        // each end, and at the wash's 0.35 it all but vanished.
        for (from, to, color) in ringed.chords {
            rings = rings.link_colored(from, to, Some(color), 0.8);
        }
    }
    Ok(rings)
}

/// The rings of a circle round `region`, outside in, each track's from its
/// file, on one reach per kind where `--same-scale` asks.
fn circle_rings(
    invocation: &Invocation,
    region: &Region,
    theme: &Theme,
    files: &mut dyn Files,
) -> Result<Vec<Ringed>, BuildError> {
    let sequence = region.seq();
    let mut built: Vec<Ringed> = Vec::new();
    if invocation.axis {
        built.push(Ringed::ruler(None));
    }
    for spec in &invocation.tracks {
        if spec.kind == Kind::Axis {
            built.push(Ringed::ruler(Some(spec)));
            continue;
        }
        let ringed = match ringed(spec, region, theme, files) {
            Ok(ringed) => ringed,
            // A file that calls the sequence by a name --rename gives it is
            // read by that name, as a figure along it reads one.
            Err(error) => {
                let mut again = None;
                for alias in called_by(invocation, sequence).into_iter().skip(1) {
                    let Ok(renamed) = Region::new(alias, region.start(), region.end()) else {
                        continue;
                    };
                    if let Ok(ringed) = ringed(spec, &renamed, theme, files) {
                        again = Some(ringed);
                        break;
                    }
                }
                match again {
                    Some(ringed) => ringed,
                    None => return Err(explained(error, spec, Some(region), files)),
                }
            }
        };
        built.push(ringed);
    }

    // One reach for every ring of one kind, so two depths read off one
    // scale, the way `--same-scale` puts two bands on one ceiling.
    if invocation.same_scale {
        for kind in [Kind::Coverage, Kind::Windows, Kind::Sequence] {
            let reach = built
                .iter()
                .filter(|ringed| ringed.kind == kind)
                .filter_map(|ringed| match &ringed.ring {
                    RingOf::Signal(ring) => Some(ring.reach()),
                    RingOf::Other(_) => None,
                })
                .fold(0.0f64, f64::max);
            if reach > 0.0 {
                for ringed in built.iter_mut().filter(|ringed| ringed.kind == kind) {
                    if let RingOf::Signal(ring) = &mut ringed.ring {
                        *ring = ring.clone().extent(reach);
                    }
                }
            }
        }
    }
    Ok(built)
}

/// The sequence a circle is drawn round, whole: named, as long as the files
/// say it is, or written from base 1, as long as written.
fn whole_sequence(invocation: &Invocation, files: &mut dyn Files) -> Result<Region, BuildError> {
    match (&invocation.region, &invocation.named) {
        (Some(written), _) => {
            // Written from base 1, which the parser has checked, it is the
            // whole of a sequence that long, unless a file says the sequence
            // runs elsewhere: a circle closed at the wrong base draws every
            // ring round the wrong angle, and looks right. Every length the
            // files state is checked, read from their headers and indexes
            // alone; a file that will not open is its track's to report.
            if let Ok(surveyed) = survey(written.seq(), invocation, files, Asked::Lengths) {
                if let Some(stated) = agreed(written.seq(), &surveyed)? {
                    if stated != written.end() {
                        return Err(BuildError::Uncircled(CircleRefusal::Length {
                            sequence: written.seq().to_string(),
                            written: written.end(),
                            stated,
                        }));
                    }
                }
            }
            Ok(written.clone())
        }
        (None, Some(name)) => {
            let surveyed = survey(name, invocation, files, Asked::Place)?;
            agreed(name, &surveyed)?;
            let placed = located(name, invocation, files, surveyed)?;
            if let Some(gene) = placed.gene {
                return Err(BuildError::Uncircled(CircleRefusal::Gene {
                    name: gene,
                    sequence: placed.region.seq().to_string(),
                }));
            }
            if let Some(reached) = placed.reached {
                return Err(BuildError::Uncircled(CircleRefusal::NoLength {
                    sequence: placed.region.seq().to_string(),
                    reached,
                    reach: placed.region.end(),
                }));
            }
            Ok(placed.region)
        }
        (None, None) => Err(BuildError::Placeless(
            crate::cli::args::ArgError::CircleWithoutPlace,
        )),
    }
}

/// The one length the files state for `sequence`, `None` where none states
/// any, and refused where two disagree.
///
/// A figure along a sequence takes the first file's word for it. A circle
/// cannot: closed at the first file's length, a VCF called on a sequence half
/// as long as the FASTA beside it cut every ring at its end, and with the
/// FASTA named first the same files closed it at the other, each without a
/// word, so the circle drawn depended on the order the files were written in.
fn agreed(sequence: &str, surveyed: &Survey<'_>) -> Result<Option<u64>, BuildError> {
    let mut said: Vec<(String, u64)> = Vec::new();
    for stated in surveyed
        .lengths
        .iter()
        .filter(|stated| stated.sequence == sequence)
    {
        let pair = (stated.file.clone(), stated.length);
        if !said.contains(&pair) {
            said.push(pair);
        }
    }
    match said.split_first() {
        None => Ok(None),
        Some(((_, first), rest)) if rest.iter().all(|(_, length)| length == first) => {
            Ok(Some(*first))
        }
        Some(_) => Err(BuildError::Uncircled(CircleRefusal::Lengths {
            sequence: sequence.to_string(),
            said,
        })),
    }
}

/// A track as a ring, before it is put on the circle.
struct Ringed {
    /// What it was drawn as, for `--same-scale` to put the rings of one kind
    /// on one reach.
    kind: Kind,
    ring: RingOf,
    /// Its name, for the key; the ring carries it too, as its tooltip.
    name: Option<String>,
    /// What its colours mean, for the key.
    legend: crate::Legend,
    /// Chords across the middle: a breakend join's two ends, each a span,
    /// and its colour.
    chords: Vec<Join>,
}

/// A chord across a circle: its two ends, each a span, and its colour.
type Join = ((u64, u64), (u64, u64), String);

/// A ring, kept as a signal where it is one, since `--same-scale` sets its
/// reach once every ring is read.
enum RingOf {
    Signal(crate::SignalRing),
    Other(Box<dyn crate::Ring>),
}

impl Ringed {
    /// The ruler, outside the first ring or where `--axis` was written.
    fn ruler(spec: Option<&TrackSpec>) -> Self {
        let mut ring = crate::AxisRing::new();
        let label = spec.and_then(|spec| spec.label.clone());
        if let Some(label) = &label {
            ring = ring.label(label);
        }
        if let Some(height) = spec.and_then(|spec| spec.height) {
            ring = ring.thickness(height);
        }
        Ringed {
            kind: Kind::Axis,
            ring: RingOf::Other(Box::new(ring)),
            // A ruler's colours mean nothing, and its marks say what it is.
            name: None,
            legend: crate::Legend::new(),
            chords: Vec::new(),
        }
    }
}

/// [`ring_of`], read again whole where a window through an index was
/// refused, as [`track`] does for a band.
fn ringed(
    spec: &TrackSpec,
    region: &Region,
    theme: &Theme,
    files: &mut dyn Files,
) -> Result<Ringed, BuildError> {
    let slurped = slurp(spec, region, ARCS as f64, files, false)?;
    if slurped.origin != Origin::Window {
        return ring_of(spec, region, theme, slurped);
    }
    match ring_of(spec, region, theme, slurped) {
        Err(BuildError::Parse { .. }) => {
            let slurped = slurp(spec, region, ARCS as f64, files, true)?;
            ring_of(spec, region, theme, slurped)
        }
        other => other,
    }
}

/// The ring a track becomes on a circle round `region`, from its file once
/// it is read: what [`built`] is to a band.
fn ring_of(
    spec: &TrackSpec,
    region: &Region,
    theme: &Theme,
    slurped: Slurped,
) -> Result<Ringed, BuildError> {
    let Slurped {
        text,
        path,
        origin,
        probe,
        held: _,
        reference: native_reference,
        most: _,
        absent,
    } = slurped;
    // A file named on its own may hold another kind than its name says, as
    // a figure along it finds, and that kind may have no ring.
    let told = match origin {
        Origin::Text | Origin::Window => spec
            .guessed
            .then(|| refine(spec.kind, probe.as_deref().unwrap_or(&text)))
            .flatten(),
        Origin::Bam | Origin::Native(_) => None,
    };
    let refined;
    let spec = match told {
        Some(kind) if kind != spec.kind => {
            if !kind.ring() {
                return Err(BuildError::Uncircled(CircleRefusal::NoRing { path, kind }));
            }
            refined = TrackSpec {
                kind,
                ..spec.clone()
            };
            &refined
        }
        _ => spec,
    };
    let name = spec.kind.flag();
    let empty = |wanted: &'static str| match &absent {
        Some(said) => BuildError::Absent {
            track: name,
            path: path.clone(),
            wanted,
            said: said.clone(),
        },
        None => BuildError::Empty {
            track: name,
            path: path.clone(),
            wanted,
        },
    };
    // The file's name, as a band is called in its gutter. A reference is
    // drawn as its GC skew, which is what its ring says.
    let label = spec.label.clone().or_else(|| match spec.kind {
        Kind::Sequence => default_label(spec).map(|stem| format!("{stem} GC skew")),
        _ => default_label(spec),
    });
    let (above, below) = match &spec.color {
        Some(color) => (color.clone(), color.clone()),
        None => (theme.color(0).to_string(), theme.color(1).to_string()),
    };
    let mut legend = crate::Legend::new();
    let mut chords = Vec::new();

    let ring = match spec.kind {
        Kind::Coverage => {
            // Painted span by span, as a band of depth is, and cut into arcs
            // with the aggregate a band reduces a pixel column with.
            let mut painted = CoverageTrack::from_spans(region, std::iter::empty());
            let spans = wrap(
                name,
                &path,
                read::signal::fold_spans(
                    &text,
                    region,
                    match origin {
                        Origin::Bam | Origin::Native(_) => Some(crate::Format::BedGraph),
                        Origin::Text | Origin::Window => spec.format,
                    },
                    |start, end, value| painted.paint(start, end, value),
                ),
            )?;
            drop(text);
            if spans == 0 {
                return Err(empty("values"));
            }
            let painted = painted.aggregate(spec.aggregate.unwrap_or(Aggregate::Max));
            // Read either side of its median, so a stretch lost dips inside
            // the line and one carried twice stands outside it: against
            // nought every arc of a sequenced genome stands outside, and a
            // loss is only a shorter one.
            let ring = crate::SignalRing::new(painted.binned(ARCS))
                .baseline_at_median()
                .colors(above.clone(), below.clone());
            let median = crate::svg::text_rounded(ring.baseline_value(), 2);
            legend = legend
                .key(format!("above its median, {median}"), above)
                .key("below it", below);
            RingOf::Signal(ring)
        }
        Kind::Windows => {
            let windows = wrap(name, &path, read::signal::windows(&text, region))?;
            if windows.is_empty() {
                return Err(empty("windows"));
            }
            // A stretch no window covers is no answer, not nought, and an
            // arc over several windows is their mean, so a value either side
            // of the line is not pulled to one side of it by a maximum.
            let mut painted = CoverageTrack::blank(region).aggregate(Aggregate::Mean);
            for window in &windows {
                painted.paint(window.start, window.end, window.value);
            }
            let ring =
                crate::SignalRing::new(painted.binned(ARCS)).colors(above.clone(), below.clone());
            legend = legend.key("above 0", above).key("below 0", below);
            RingOf::Signal(ring)
        }
        Kind::Sequence => {
            let reference = match native_reference {
                Some(reference) => reference,
                None => sequence(name, &path, &text, region)?,
            };
            let (from, bases) = reference.clip(region)?;
            let window = (region.len() / ARCS as u64).max(1);
            let skew = WindowTrack::gc_skew(from, &bases, window);
            let ring = crate::SignalRing::new(skew.windows().to_vec())
                .colors(above.clone(), below.clone());
            legend = legend
                .key("more G than C", above)
                .key("more C than G", below);
            RingOf::Signal(ring)
        }
        Kind::Features => {
            let format = match origin {
                Origin::Native(_) => Some(crate::Format::Bed),
                _ => spec.format,
            };
            let features = wrap(name, &path, read::interval::features(&text, region, format))?;
            if features.is_empty() {
                return Err(empty("features"));
            }
            let named = features
                .iter()
                .filter(|feature| feature.name.is_some())
                .count();
            let reverse = features
                .iter()
                .any(|feature| feature.strand == crate::Strand::Reverse);
            let unknown = features
                .iter()
                .any(|feature| feature.strand == crate::Strand::Unknown);
            let forward = features
                .iter()
                .any(|feature| feature.strand == crate::Strand::Forward);
            // The colours a band of features paints, a strand each: the
            // forward strand on the outer half of the ring and the reverse on
            // the inner. One colour asked for paints both, as it does a band.
            if spec.color.is_none() {
                if forward || unknown {
                    let said = match (forward, unknown) {
                        (true, false) => "forward strand, outer half",
                        (false, true) => "no strand, outer half",
                        _ => "forward or no strand, outer half",
                    };
                    legend = legend.key(said, above.clone());
                }
                if reverse {
                    legend = legend.key("reverse strand, inner half", below.clone());
                }
            }
            let mut ring = crate::FeatureRing::new(features)
                .show_names(!spec.no_names && named <= NAMED_ON_A_RING);
            if spec.color.is_some() {
                ring = ring.colors(above, below);
            }
            RingOf::Other(Box::new(named_ring(ring, label.as_deref())))
        }
        Kind::Variants => {
            let variants = wrap(
                name,
                &path,
                recorded(origin, &text, read::point::variants(&text, region)),
            )?;
            if variants.is_empty() {
                return Err(empty("variants"));
            }
            // A colour each consequence, from the most damaging down, as a
            // band of calls deals them. The reader names one for every call,
            // from its shape where nothing annotated it, so the colour after
            // the last is for a caller that hands over none.
            let ranked = read::point::ranked(&variants);
            let slot = |category: Option<&String>| -> usize {
                category
                    .and_then(|category| ranked.iter().position(|known| known == category))
                    .unwrap_or(ranked.len())
            };
            let marks: Vec<(u64, usize)> = variants
                .iter()
                .map(|variant| (variant.pos, slot(variant.category.as_ref())))
                .collect();
            for (index, category) in ranked.iter().enumerate() {
                legend = legend.line(category.clone(), theme.color(index));
            }
            let mut ring = crate::MarkerRing::categorised(marks);
            if let Some(height) = spec.height {
                ring = ring.thickness(height);
            }
            if let Some(label) = &label {
                ring = ring.label(label);
            }
            RingOf::Other(Box::new(ring))
        }
        Kind::Structural => {
            let found = wrap(
                name,
                &path,
                recorded(origin, &text, read::structural::variants(&text, region)),
            )?;
            if found.records == 0 {
                return Err(empty("variant calls"));
            }
            if found.variants.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "structural calls",
                    held: found.records,
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }
            // The colour a band of calls gives each class. A call that covers
            // reference is its footprint on the ring, an insertion one base
            // held open by the ring's floor, and a breakend join on the
            // sequence is a chord between its two ends.
            let colors = StructuralTrack::new(Vec::new());
            let mut footprints = Vec::new();
            let mut kinds: Vec<crate::SvKind> = Vec::new();
            for call in &found.variants {
                let color = colors.color_of(call.kind, theme);
                if !kinds.contains(&call.kind) {
                    kinds.push(call.kind);
                }
                if call.kind == crate::SvKind::Translocation {
                    // A join is read as the span from the base after its POS
                    // to its mate's base, and a band names the first and the
                    // last of it, so the chord joins those two. Taken from
                    // the span's end, the chord's target was the base after
                    // the mate, which neither the file nor the band names.
                    let last = call.end.max(call.start + 1);
                    chords.push(((call.start, call.start + 1), (last - 1, last), color));
                    continue;
                }
                let called = match &call.name {
                    Some(id) => format!("{} {id}", call.kind.name()),
                    None => call.kind.name().to_string(),
                };
                footprints.push(
                    crate::Feature::new(call.start, call.end.max(call.start + 1))
                        .name(called)
                        .color(color),
                );
            }
            for kind in &kinds {
                let color = colors.color_of(*kind, theme);
                legend = if *kind == crate::SvKind::Translocation {
                    legend.line(format!("{}, across the middle", kind.name()), color)
                } else {
                    legend.key(kind.name(), color)
                };
            }
            let shown = !spec.no_names && footprints.len() <= NAMED_ON_A_RING;
            let mut ring = crate::FeatureRing::new(footprints)
                .split_strands(false)
                .show_names(shown);
            if let Some(height) = spec.height {
                ring = ring.thickness(height);
            }
            RingOf::Other(Box::new(named_ring(ring, label.as_deref())))
        }
        // The parser lets no other kind onto a circle, and `refine` turns a
        // file into one with no ring only where the refusal above says so.
        _ => unreachable!("the parser lets only kinds with a ring onto a circle"),
    };
    let ring = match ring {
        RingOf::Signal(mut ring) => {
            if let Some(height) = spec.height {
                ring = ring.thickness(height);
            }
            if let Some(label) = &label {
                ring = ring.label(label);
            }
            RingOf::Signal(ring)
        }
        other => other,
    };
    Ok(Ringed {
        kind: spec.kind,
        ring,
        name: label,
        legend,
        chords,
    })
}

/// A ring of features under its name, where it has one.
fn named_ring(ring: crate::FeatureRing, label: Option<&str>) -> crate::FeatureRing {
    match label {
        Some(label) => ring.label(label),
        None => ring,
    }
}

/// The same, for a panel of a sheet of several places where `sheet` is
/// given: a track with nothing in the place is drawn as a band that says so
/// rather than refusing the figure, and what became of each shade is put in
/// `sheet` for the sheet to settle once every panel has had it, beside where
/// each gene to shade is, which the first panel looks up for them all.
///
/// Beside the figure, the width it would letter its bases at, where it draws
/// a reference as blocks too narrow for their letters, for the caller to say
/// once: a sheet says the widest any of its panels needs.
fn build_one(
    invocation: &Invocation,
    files: &mut dyn Files,
    mut parsed: impl FnMut(&str, &str) -> Option<Tree>,
    theme: Theme,
    window: Option<&Region>,
    sheet: Option<&mut Shading>,
) -> Result<(Built, Option<u64>), BuildError> {
    let tolerant = sheet.is_some();
    let mut kept = KeptStdin { files, stdin: None };
    let files: &mut dyn Files = &mut kept;
    if invocation.genome_wide() && window.is_none() {
        return build_genome(invocation, files, theme).map(|built| (built, None));
    }
    // A figure of phylogenies and variable-site panels names no region, and
    // none of its tracks asks the window anything. The figure still wants
    // one to lay its width out over, so it is given one that nothing prints:
    // a figure that shows no window draws no locus and no ruler.
    let unnamed = Region::new("phylogeny", 0, 1).expect("a one-base window is a window");
    // A place named by a word is found in the figure's own files.
    let placed = match (&invocation.region, &invocation.named) {
        (None, Some(name)) => Some(place(name, invocation, files)?),
        _ => None,
    };
    let known = window
        .or(invocation.region.as_ref())
        .or(placed.as_ref().map(|placed| &placed.region));
    // An alignment is its own place: a figure of one, named nowhere, is laid
    // over all its columns. It was refused until a region was made up for it,
    // and the one to make up was the name of one of its rows. A table of
    // counts over time, of estimates, of sites and a read's signal are their
    // own places the same way.
    let decimals = time_decimals(invocation, files);
    let columns = match known {
        None => own_place(invocation, files, decimals)?,
        Some(_) => None,
    };
    // A place written for a continuous time names it in its own units, as
    // `year:2010-2016`, and the tables are read to a thousandth of one.
    let rescaled = match (&invocation.region, decimals) {
        (Some(written), 1..) if window.is_none() && all_times(invocation) => {
            let scale = 10u64.pow(decimals);
            Region::new(
                written.seq(),
                (written.start() + 1).saturating_mul(scale),
                written.end().saturating_mul(scale).saturating_add(1),
            )
            .ok()
        }
        _ => None,
    };
    let region = rescaled
        .as_ref()
        .or(known)
        .or(columns.as_ref())
        .unwrap_or(&unnamed);
    let counting = counted(invocation, region, decimals);
    // Along a genome where the place is one, and not where the ruler counts
    // weeks, sites, samples or columns.
    let along = known
        .filter(|_| counting.is_none() && rescaled.is_none())
        .cloned();
    // A sequence no file gives the length of ends where its rows do, which
    // is a figure worth drawing and an end worth saying where it came from:
    // three simulated users read it as the end of the chromosome.
    let reached = placed
        .as_ref()
        .and_then(|placed| placed.reached.as_ref())
        .filter(|_| window.is_none());
    if let Some(file) = reached {
        let sequence = region.seq();
        // A FASTA that calls the sequence otherwise draws it only once the
        // figure is placed on its name, which --rename then reads the table
        // under; where a --rename already did that, the FASTA is all it needs.
        let renamed = invocation.renames.iter().any(|(_, to)| to == sequence);
        let otherwise = if renamed {
            String::new()
        } else {
            format!(
                "; where they give it another name, place the figure on that name and add \
                 --rename {sequence}=THAT_NAME"
            )
        };
        files.note(&format!(
            "{sequence} is drawn to {}, as far as {file} reaches, since no file says how \
             long it is. To draw all of it, write the span, as {sequence}:1-LENGTH, or add \
             its FASTA or BAM{otherwise}",
            crate::track::axis::group_thousands(region.end())
        ));
    }
    let mut plot = Plot::over(region.clone());
    // A figure placed by a gene's name is about that gene, and says so.
    let gene = placed.as_ref().and_then(|placed| placed.gene.as_ref());
    if let Some(title) = invocation.title.as_ref().or(gene) {
        plot = plot.title(title);
    }
    // The key to every colour a track paints by a category, gathered as the
    // tracks are built and drawn under the ruler.
    let mut legend = crate::track::legend::Legend::new();
    // The reference the figure draws, which a pileup given none of its own
    // reads its reads against: naming the same FASTA twice was the only way to
    // see a mismatch, and a pileup without it drew every read as agreeing.
    let reference = invocation
        .tracks
        .iter()
        .find(|track| track.kind == Kind::Sequence)
        .and_then(|track| track.source.clone());
    if let Some(width) = invocation.width {
        plot = plot.width(width);
    }
    plot = plot.theme(theme.clone());
    // A week, a site or a sample says no more at the top right than the
    // ruler says under it.
    if !invocation.region_label || counting.as_ref().is_some_and(|counting| !counting.columns) {
        plot = plot.remove_region_label();
    }
    // A ruler of bases would print a year as 2,015 and a count of samples in
    // kilobases, so a figure whose ruler counts something else puts in its
    // own, once the tracks are in.
    if !invocation.axis || counting.is_some() {
        plot = plot.remove_axis();
    }

    colored_as_asked(invocation, files)?;
    for spec in &invocation.tracks {
        // The ruler is the one track that reads nothing, and the one the plot
        // has to be told about so it does not append a second.
        if spec.kind == Kind::Axis {
            let mut axis = plot.add_axis();
            if let Some(counting) = &counting {
                let places = counting.decimals;
                axis = axis
                    .adjust(|axis| axis.counting().decimals(places))
                    .label(&counting.unit);
            }
            if let Some(label) = &spec.label {
                axis = axis.label(label);
            }
            if let Some(height) = spec.height {
                axis = axis.adjust(|track| track.height(height));
            }
            plot = axis.done();
            continue;
        }
        // The codon ruler reads no file either: its gene is the annotation's
        // and its letters the reference's, so it is built from the figure.
        if spec.kind == Kind::Codons {
            let ruler = codons(
                spec,
                invocation,
                placed.as_ref(),
                region,
                reference.as_ref(),
                files,
            )?;
            plot = plot.add_boxed(ruler);
            continue;
        }
        let genotyped = invocation
            .tracks
            .iter()
            .any(|other| other.kind == Kind::Genotypes && other.source == spec.source);
        let context = Context {
            region,
            theme: &theme,
            reference: reference.as_ref(),
            decimals,
            genotyped,
            colors: &invocation.colors,
            width: invocation.width.unwrap_or(900.0),
        };
        let built = match track(spec, &context, files, &mut parsed, &mut legend) {
            Ok(built) => built,
            // A file that calls the figure's sequence by a name --rename
            // gives it is read by that name, and drawn under the figure's.
            Err(error) => {
                let mut again = None;
                for alias in called_by(invocation, region.seq()).into_iter().skip(1) {
                    let Ok(renamed) = Region::new(alias, region.start(), region.end()) else {
                        continue;
                    };
                    let context = Context {
                        region: &renamed,
                        theme: &theme,
                        reference: reference.as_ref(),
                        decimals,
                        genotyped,
                        colors: &invocation.colors,
                        width: invocation.width.unwrap_or(900.0),
                    };
                    if let Ok(built) = track(spec, &context, files, &mut parsed, &mut legend) {
                        again = Some(built);
                        break;
                    }
                }
                match (again, error) {
                    (Some(built), _) => built,
                    // One place of several with nothing on it is a finding,
                    // a gene with no calls, and refused it took every other
                    // panel of the figure with it.
                    (
                        None,
                        BuildError::Empty { wanted, .. } | BuildError::Absent { wanted, .. },
                    ) if tolerant => Box::new(Nothing {
                        label: spec.label.clone().or_else(|| default_label(spec)),
                        said: format!("no {wanted} here"),
                    }),
                    (None, error) => return Err(explained(error, spec, known, files)),
                }
            }
        };
        plot = plot.add_boxed(built);
    }
    // After the ruler, which closing the plot puts in, so the key is not taken
    // for a track measured against it.
    let mut figure = plot.into_figure();
    if invocation.same_scale {
        figure = figure.same_scale();
    }
    if let Some(counting) = counting.as_ref().filter(|_| invocation.axis) {
        figure = figure.push_ruler(
            crate::AxisTrack::new()
                .counting()
                .decimals(counting.decimals)
                .label(&counting.unit),
        );
    }
    // The stretches asked to be shaded, on the panel's own sequence, settled
    // here for a figure of one place and by the sheet for one of several.
    let over = known.unwrap_or(region);
    figure = match sheet {
        Some(shading) => shaded(figure, invocation, files, region, over, decimals, shading)?,
        None => {
            let mut shading = Shading::new(invocation);
            let figure = shaded(
                figure,
                invocation,
                files,
                region,
                over,
                decimals,
                &mut shading,
            )?;
            settle_shades(invocation, &shading.fates, files)?;
            figure
        }
    };
    // What the tracks need explained at the zoom the figure is drawn at, as
    // the colours of bases too narrow for their letters.
    let key = figure.key();
    let bases = theme.bases.legend();
    // Only for a reference drawn as blocks: an alignment keys its colours the
    // same way, and its reader came for the pattern down the rows, not for
    // the letters a width of several thousand pixels would write.
    let sequence = invocation
        .tracks
        .iter()
        .any(|spec| spec.kind == Kind::Sequence);
    // And how wide the figure would have to be for the letters, which a
    // reader asked to show a sequence came for.
    let letters = (sequence && bases.items().iter().all(|item| key.items().contains(item)))
        .then(|| {
            let px = figure.px_per_bp();
            let span = region.len() as f64;
            let now = figure.dimensions().0;
            ((now + (crate::track::sequence::LETTER_PX - px) * span) / 100.0).ceil() * 100.0
        })
        .map(|wanted| wanted as u64);
    gather(&mut legend, &key);
    if invocation.legend && !legend.is_empty() {
        figure = figure.push(crate::track::legend::LegendTrack::new(legend.clone()));
    }
    Ok((
        Built {
            figure,
            along,
            legend,
        },
        letters,
    ))
}

/// A track of one place of several that has nothing there: the band it would
/// have been, under the label it would have had, saying so.
struct Nothing {
    label: Option<String>,
    said: String,
}

impl Track for Nothing {
    fn height(&self, _scale: &crate::scale::Scale) -> f64 {
        26.0
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn noun(&self) -> &str {
        "an empty track"
    }

    fn draw(&self, ctx: &mut crate::track::DrawContext<'_>) {
        let size = (ctx.theme.font_size - 1.0).max(6.0);
        ctx.svg.text(
            ctx.band.x + ctx.band.w / 2.0,
            ctx.band.y + ctx.band.h / 2.0 + size * 0.35,
            &self.said,
            &ctx.theme.muted,
            size,
            crate::svg::Anchor::Middle,
        );
    }
}

/// A scan's threshold line, where the command line draws one: as a p-value
/// wherever the file held p-values, and in the file's own units otherwise.
/// And the top of its scale, where `--max` pins one.
fn thresholded(
    track: ManhattanTrack,
    spec: &TrackSpec,
    p_values: bool,
    path: &str,
) -> Result<ManhattanTrack, BuildError> {
    let track = match spec.max {
        Some(max) => track.max(max),
        None => track,
    };
    Ok(match spec.threshold {
        None => track,
        Some(Threshold::GenomeWide) => track.genome_wide_threshold(),
        // In the units the file is in, so a p-value where the file held
        // p-values, drawn where its points are.
        Some(Threshold::At(value)) if p_values => {
            if !(value > 0.0 && value <= 1.0) {
                return Err(BuildError::NotAPValue {
                    track: spec.kind.flag(),
                    path: path.to_string(),
                    given: value,
                });
            }
            track.p_value_threshold(value)
        }
        Some(Threshold::At(value)) => track.threshold(value),
    })
}

/// The order a reader counts chromosomes in: the numbered ones by number,
/// with or without `chr` in front, then X, Y and the mitochondrion, then any
/// other sequence, contigs and unplaced scaffolds, as their names sort.
fn chromosome_order(a: &str, b: &str) -> std::cmp::Ordering {
    fn rank(name: &str) -> (u8, u64) {
        let lower = name.to_ascii_lowercase();
        let bare = lower.strip_prefix("chr").unwrap_or(&lower);
        if let Ok(number) = bare.parse::<u64>() {
            return (0, number);
        }
        match bare {
            "x" => (1, 0),
            "y" => (2, 0),
            "m" | "mt" => (3, 0),
            _ => (4, 0),
        }
    }
    rank(a)
        .cmp(&rank(b))
        .then_with(|| crate::track::traits::natural(a, b))
}

/// A figure with no place: every sequence its files name, end to end, in the
/// order a reader counts chromosomes, with the bands a genome-wide plot is
/// read by and each sequence named under the tracks in place of a ruler of
/// positions no file uses.
///
/// A scan, a signal, windows and a segment table are each read whole, every
/// row on the sequence it names, under the name `--rename` gives it, so a
/// table that calls a chromosome `1` and a bedGraph that calls it `chr1` are
/// laid on one sequence where `--rename 1=chr1` says they are the same.
///
/// Each sequence is as long as the furthest any file reaches on it, which is
/// how a scan is drawn: an association table says where its markers are and
/// not how long the chromosomes they are on run, and a bedGraph says where
/// its depth was counted. A bigWig says how long each sequence it names is,
/// and is taken at its word. A figure whose files name one sequence is that
/// sequence drawn whole, ruler and all, as though it had been written.
fn build_genome(
    invocation: &Invocation,
    files: &mut dyn Files,
    theme: Theme,
) -> Result<Built, BuildError> {
    let width = invocation.width.unwrap_or(900.0);
    let mut spread = Vec::with_capacity(invocation.tracks.len());
    for spec in &invocation.tracks {
        if spec.source.is_some() {
            spread.push(spread_of(spec, invocation, files)?);
        }
    }
    let mut lengths: Vec<(String, u64)> = Vec::new();
    let mut found: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for track in &spread {
        for (sequence, reach) in track.reaches() {
            match found.get(&sequence) {
                Some(at) => lengths[*at].1 = lengths[*at].1.max(reach),
                None => {
                    found.insert(sequence.clone(), lengths.len());
                    lengths.push((sequence, reach));
                }
            }
        }
    }
    // One sequence is no genome to lay out, and is the figure of it, as long
    // as the files reach on it: written as a place, it has a ruler in its own
    // bases and a page can move along it.
    if let [(sequence, length)] = &lengths[..] {
        if let Ok(whole) = Region::new(sequence, 0, (*length).max(1)) {
            let mut one = invocation.clone();
            one.region = Some(whole);
            return build_one(&one, files, |_, _| None, theme, None, None).map(|(built, _)| built);
        }
    }
    lengths.sort_by(|a, b| chromosome_order(&a.0, &b.0));
    let genome = crate::Genome::new(lengths);

    let mut plot = Plot::over(genome.region())
        .remove_region_label()
        .remove_axis()
        .theme(theme);
    if let Some(title) = &invocation.title {
        plot = plot.title(title);
    }
    if let Some(width) = invocation.width {
        plot = plot.width(width);
    }
    for track in spread {
        plot = plot.add_boxed(track.laid(&genome, width)?);
    }
    let mut figure = plot.add_genome(genome.clone()).into_figure();
    figure = shaded_genome(figure, invocation, &genome, files)?;
    if invocation.same_scale {
        figure = figure.same_scale();
    }
    let legend = figure.key();
    if invocation.legend && !legend.is_empty() {
        figure = figure.push(crate::track::legend::LegendTrack::new(legend.clone()));
    }
    Ok(Built {
        figure,
        along: None,
        legend,
    })
}

/// One track of a figure across a whole genome, read: what it draws on each
/// sequence its file names, each under the name the figure gives it.
struct Spread {
    /// The track, as its file turned out to be where its name was a guess.
    spec: TrackSpec,
    /// What its file is called in a message.
    path: String,
    /// What it holds.
    held: Spreading,
}

/// What a track across a whole genome holds, by its kind.
enum Spreading {
    /// A scan's points on each sequence, 0-based on it.
    Scan {
        sequences: Vec<(String, Vec<crate::Association>)>,
        p_values: bool,
    },
    /// A signal painted on each sequence from its base nought, and how far
    /// its rows reach there.
    Signal(Vec<(String, CoverageTrack, u64)>),
    /// A bigWig, read once the genome is laid out, at the zoom level its
    /// length over the figure's width wants: each sequence it names under the
    /// figure's name and its own, and the length its index gives it.
    BigWig {
        file: Box<dyn Seekable>,
        sequences: Vec<(String, String, u64)>,
    },
    /// Windows on each sequence.
    Windows(Vec<(String, Vec<crate::Window>)>),
    /// The segments called on each sequence, at the track's ploidy.
    Copies(Vec<(String, Vec<crate::CopyNumberSegment>)>),
}

impl Spread {
    /// Each sequence the track names, under the figure's name for it, and
    /// how far along it the track reaches.
    fn reaches(&self) -> Vec<(String, u64)> {
        match &self.held {
            Spreading::Scan { sequences, .. } => sequences
                .iter()
                .map(|(name, points)| {
                    let end = points
                        .iter()
                        .map(|point| point.pos.saturating_add(1))
                        .max()
                        .unwrap_or(1);
                    (name.clone(), end)
                })
                .collect(),
            Spreading::Signal(sequences) => sequences
                .iter()
                .map(|(name, _, reach)| (name.clone(), *reach))
                .collect(),
            Spreading::BigWig { sequences, .. } => sequences
                .iter()
                .map(|(name, _, length)| (name.clone(), *length))
                .collect(),
            Spreading::Windows(sequences) => sequences
                .iter()
                .map(|(name, windows)| {
                    (
                        name.clone(),
                        windows.iter().map(|window| window.end).max().unwrap_or(1),
                    )
                })
                .collect(),
            Spreading::Copies(sequences) => sequences
                .iter()
                .map(|(name, segments)| {
                    (
                        name.clone(),
                        segments
                            .iter()
                            .map(|segment| segment.end)
                            .max()
                            .unwrap_or(1),
                    )
                })
                .collect(),
        }
    }

    /// The track, laid over `genome`, in a figure `width` pixels wide.
    fn laid(self, genome: &crate::Genome, width: f64) -> Result<Box<dyn Track>, BuildError> {
        let Spread { spec, path, held } = self;
        let spec = &spec;
        let name = spec.kind.flag();
        let label = spec.label.clone().or_else(|| default_label(spec));
        // Where each sequence starts, by name: an assembly can be a hundred
        // thousand scaffolds, and a walk along the genome from its first to
        // find each of them is billions of comparisons, which with the same
        // walk in `gathered` took fourteen seconds to draw their windows.
        let mut starts: std::collections::HashMap<&str, u64> =
            std::collections::HashMap::with_capacity(genome.len());
        for (sequence, start, _) in genome.spans() {
            starts.entry(sequence).or_insert(start);
        }
        let offset = |sequence: &str| starts.get(sequence).copied().unwrap_or(0);
        // Every position is laid with a saturating sum, as the genome's own
        // lengths are: a length is whatever a file said it was, and two of
        // them can sum past what a u64 holds, where a plain sum panics.
        Ok(match held {
            Spreading::Scan {
                sequences,
                p_values,
            } => {
                let points: Vec<crate::Association> = sequences
                    .iter()
                    .flat_map(|(sequence, points)| {
                        let offset = offset(sequence);
                        points.iter().map(move |point| {
                            crate::Association::new(offset.saturating_add(point.pos), point.value)
                        })
                    })
                    .collect();
                let mut track = ManhattanTrack::new(points).bands(genome.boundaries());
                if p_values {
                    track = track.axis_title("-log10 p");
                }
                let mut track = thresholded(track, spec, p_values, &path)?;
                if let Some(height) = spec.height {
                    track = track.height(height);
                }
                Box::new(named(track, label, ManhattanTrack::label))
            }
            Spreading::Signal(sequences) => {
                let mut parts: std::collections::HashMap<String, CoverageTrack> = sequences
                    .into_iter()
                    .map(|(sequence, track, _)| (sequence, track))
                    .collect();
                let track = CoverageTrack::end_to_end(
                    genome
                        .sequences()
                        .iter()
                        .map(|sequence| (parts.remove(&sequence.name), sequence.length)),
                );
                Box::new(named(
                    dressed_coverage(track, spec, None),
                    label,
                    CoverageTrack::label,
                ))
            }
            Spreading::BigWig {
                mut file,
                sequences,
            } => {
                // As many bases to a pixel as the whole figure is wide, as a
                // window of one sequence is read: a zoom level is never
                // coarser than the drawing for it. Windows and the scores
                // under bases are read as written, as they are over a window.
                let (per_pixel, aggregate) = match spec.kind {
                    Kind::Coverage => (
                        genome.total() as f64 / width.max(1.0),
                        spec.aggregate.unwrap_or(Aggregate::Max),
                    ),
                    _ => (1.0, Aggregate::Max),
                };
                let read =
                    read::bigwig::genome(&mut file, per_pixel, aggregate).map_err(|cause| {
                        BuildError::Open {
                            track: name,
                            path: path.clone(),
                            cause: unreadable(cause),
                        }
                    })?;
                // The figure's name for each of the file's, by the file's.
                let names: std::collections::HashMap<&str, &str> = sequences
                    .iter()
                    .map(|(figure, own, _)| (own.as_str(), figure.as_str()))
                    .collect();
                let figured = |own: &str| names.get(own).map(|figure| figure.to_string());
                if read.iter().all(|(_, _, signal)| signal.spans.is_empty()) {
                    return Err(BuildError::Empty {
                        track: name,
                        path,
                        wanted: "values",
                    });
                }
                if spec.kind == Kind::Windows {
                    let mut windows = Vec::new();
                    for (own, _, signal) in &read {
                        let Some(sequence) = figured(own) else {
                            continue;
                        };
                        let offset = offset(&sequence);
                        windows.extend(signal.spans.iter().map(|(start, end, value)| {
                            crate::Window::new(
                                offset.saturating_add(*start),
                                offset.saturating_add(*end),
                                *value,
                            )
                        }));
                    }
                    return Ok(Box::new(named(
                        dressed_windows(windows, spec),
                        label,
                        WindowTrack::label,
                    )));
                }
                let mut most: Option<f64> = None;
                let mut parts: std::collections::HashMap<String, CoverageTrack> =
                    std::collections::HashMap::new();
                for (own, _, signal) in read {
                    let Some(sequence) = figured(&own) else {
                        continue;
                    };
                    if let Some(high) = signal.most {
                        most = Some(most.map_or(high, |most| most.max(high)));
                    }
                    let mut painted = CoverageTrack::unbounded(0);
                    for (start, end, value) in signal.spans {
                        painted.paint(start, end, value);
                    }
                    parts.insert(sequence, painted);
                }
                let track = CoverageTrack::end_to_end(
                    genome
                        .sequences()
                        .iter()
                        .map(|sequence| (parts.remove(&sequence.name), sequence.length)),
                );
                Box::new(named(
                    dressed_coverage(track, spec, most),
                    label,
                    CoverageTrack::label,
                ))
            }
            Spreading::Windows(sequences) => {
                let windows: Vec<crate::Window> = sequences
                    .iter()
                    .flat_map(|(sequence, windows)| {
                        let offset = offset(sequence);
                        windows.iter().map(move |window| {
                            crate::Window::new(
                                offset.saturating_add(window.start),
                                offset.saturating_add(window.end),
                                window.value,
                            )
                        })
                    })
                    .collect();
                Box::new(named(
                    dressed_windows(windows, spec),
                    label,
                    WindowTrack::label,
                ))
            }
            Spreading::Copies(sequences) => {
                // Checked by the parser, so this is reached only from an
                // Invocation built by hand.
                let Some(ploidy) = spec.ploidy else {
                    return Err(BuildError::MissingPloidy { track: name });
                };
                let segments: Vec<crate::CopyNumberSegment> = sequences
                    .iter()
                    .flat_map(|(sequence, segments)| {
                        let offset = offset(sequence);
                        segments
                            .iter()
                            .map(move |segment| crate::CopyNumberSegment {
                                start: offset.saturating_add(segment.start),
                                end: offset.saturating_add(segment.end),
                                copy: segment.copy,
                            })
                    })
                    .collect();
                let mut track = CopyNumberTrack::at_ploidy(segments, ploidy).across(genome);
                if let Some(height) = spec.height {
                    track = track.height(height);
                }
                Box::new(named(track, label, CopyNumberTrack::label))
            }
        })
    }
}

/// Reads one track of a figure across a whole genome: every row of its file,
/// each on the sequence it names, under the name the figure gives it.
fn spread_of(
    spec: &TrackSpec,
    invocation: &Invocation,
    files: &mut dyn Files,
) -> Result<Spread, BuildError> {
    let name = spec.kind.flag();
    let source = spec.source.as_ref().expect("a track with a file");
    let opened = native(files, source).map_err(|cause| BuildError::Open {
        track: name,
        path: called(source),
        cause,
    })?;
    if let Some((binary, mut file)) = opened {
        let path = called(source);
        let refused = |format: bool| BuildError::OtherTrack {
            track: name,
            path: path.clone(),
            binary,
            format,
        };
        // `--format` says what the columns of a text file are, and a bigWig
        // says it itself.
        if spec.format.is_some() {
            return Err(refused(true));
        }
        if binary != Binary::BigWig || !matches!(spec.kind, Kind::Coverage | Kind::Windows) {
            return Err(refused(false));
        }
        let sequences = read::bigwig::sequences(&mut file)
            .map_err(|cause| BuildError::Open {
                track: name,
                path: path.clone(),
                cause: unreadable(cause),
            })?
            .into_iter()
            .map(|(own, length)| (renamed(invocation, own.clone()), own, length))
            .collect();
        return Ok(Spread {
            spec: spec.clone(),
            path,
            held: Spreading::BigWig { file, sequences },
        });
    }
    let (text, path) =
        fetch(name, source, files).map_err(|error| explained(error, spec, None, files))?;
    // A `.bed` named on its own is told once it is read: mosdepth's depth in
    // windows is a bedGraph by another name, and features need a place.
    let mut spec = spec.clone();
    if spec.guessed {
        if let Some(kind) = refine(spec.kind, &text) {
            spec.kind = kind;
        }
    }
    if !spec.kind.genome_wide() {
        return Err(BuildError::Placeless(if invocation.tracks.len() == 1 {
            crate::cli::args::ArgError::NoRegion
        } else {
            crate::cli::args::ArgError::NotGenomeWide {
                track: spec.kind.flag(),
                file: Some(path),
                tied: None,
            }
        }));
    }
    let name = spec.kind.flag();
    let held = match spec.kind {
        Kind::Manhattan => {
            let read = wrap(name, &path, read::point::genome_associations(&text))?;
            if read.sequences.iter().all(|(_, points)| points.is_empty()) {
                return Err(BuildError::Empty {
                    track: name,
                    path,
                    wanted: "association statistics",
                });
            }
            Spreading::Scan {
                sequences: gathered(invocation, read.sequences),
                p_values: read.p_values,
            }
        }
        Kind::Coverage => {
            // Painted as each row is read, a sequence at a time, each on a
            // profile of its own from its base nought: the rows of one
            // sequence come together, and a profile across the whole genome
            // painted in the file's order, where it is not the order the
            // genome is laid out in, splices runs into the middle of millions.
            let mut sequences: Vec<(String, CoverageTrack, u64)> = Vec::new();
            let mut at: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
            let mut last: Option<(String, usize)> = None;
            let painted = wrap(
                name,
                &path,
                read::signal::fold_genome_spans(&text, spec.format, |own, start, end, value| {
                    let index = match &last {
                        Some((called, index)) if called == own => *index,
                        _ => {
                            let figure = renamed(invocation, own.to_string());
                            let index = *at.entry(figure.clone()).or_insert_with(|| {
                                sequences.push((figure, CoverageTrack::unbounded(0), 0));
                                sequences.len() - 1
                            });
                            last = Some((own.to_string(), index));
                            index
                        }
                    };
                    let (_, track, reach) = &mut sequences[index];
                    track.paint(start, end, value);
                    *reach = (*reach).max(end);
                }),
            )?;
            if painted == 0 {
                return Err(BuildError::Empty {
                    track: name,
                    path,
                    wanted: "values",
                });
            }
            Spreading::Signal(sequences)
        }
        Kind::Windows => {
            let read = wrap(name, &path, read::signal::genome_windows(&text))?;
            if read.iter().all(|(_, windows)| windows.is_empty()) {
                return Err(BuildError::Empty {
                    track: name,
                    path,
                    wanted: "windows",
                });
            }
            Spreading::Windows(gathered(invocation, read))
        }
        Kind::CopyNumber => {
            let Some(ploidy) = spec.ploidy else {
                return Err(BuildError::MissingPloidy { track: name });
            };
            one_sample(&spec, &path, &text)?;
            let read = wrap(
                name,
                &path,
                read::segments::genome_copy_numbers(&text, ploidy, spec.sample.as_deref()),
            )?;
            if read.records == 0 {
                return Err(BuildError::Empty {
                    track: name,
                    path,
                    wanted: "segments",
                });
            }
            if read.sequences.is_empty() {
                return Err(BuildError::Empty {
                    track: name,
                    path,
                    wanted: "called segments",
                });
            }
            Spreading::Copies(gathered(invocation, read.sequences))
        }
        _ => unreachable!("only the kinds drawn across a genome are read across one"),
    };
    Ok(Spread { spec, path, held })
}

/// Each sequence of a file read whole, under the name the figure gives it,
/// with what two of the file's names `--rename` makes one sequence hold
/// joined, in the order the file first names them.
///
/// Each name is looked up where it has been put, rather than looked for
/// among the names before it, which across an assembly of a hundred thousand
/// scaffolds is five billion comparisons.
fn gathered<T>(invocation: &Invocation, sequences: Vec<(String, Vec<T>)>) -> Vec<(String, Vec<T>)> {
    let mut joined: Vec<(String, Vec<T>)> = Vec::with_capacity(sequences.len());
    let mut at: std::collections::HashMap<String, usize> =
        std::collections::HashMap::with_capacity(sequences.len());
    for (own, rows) in sequences {
        let figure = renamed(invocation, own);
        match at.get(&figure) {
            Some(index) => joined[*index].1.extend(rows),
            None => {
                at.insert(figure.clone(), joined.len());
                joined.push((figure, rows));
            }
        }
    }
    joined
}

/// Refuses a segment table that names no sample where one was asked for, or
/// several where none was, before its segments are read.
fn one_sample(spec: &TrackSpec, path: &str, text: &str) -> Result<(), BuildError> {
    let name = spec.kind.flag();
    let held = wrap(name, path, read::segments::samples(text))?;
    if spec.sample.is_some() && held.is_empty() {
        // A flag accepted and then ignored gives a figure that is not the one
        // asked for and does not look wrong: this table names no samples, so
        // the whole of it would be drawn under a name the command asked to
        // pick out of it.
        return Err(BuildError::Ambiguous {
            track: name,
            path: path.to_string(),
            flag: "--sample",
            choices: vec!["no sample column".to_string()],
        });
    }
    if held.len() > 1 && spec.sample.is_none() {
        return Err(BuildError::Ambiguous {
            track: name,
            path: path.to_string(),
            flag: "--sample",
            choices: held,
        });
    }
    Ok(())
}

/// A coverage track as the options after it ask: how a column of many bases
/// is summed up, its style, its scale and its colour. `most` lifts the scale
/// to the most of the values a bigWig's zoom level summarises.
fn dressed_coverage(painted: CoverageTrack, spec: &TrackSpec, most: Option<f64>) -> CoverageTrack {
    let mut track = painted.aggregate(spec.aggregate.unwrap_or(Aggregate::Max));
    // A bigWig's zoom level painted with its bins' means or least, scaled to
    // the most of the values under them, as the values as written scale it.
    if let Some(most) = most {
        track = track.reaching(most);
    }
    if let Some(style) = spec.style.and_then(|style| style.coverage()) {
        track = track.style(style);
    }
    if spec.log {
        track = track.log_scale(true);
    }
    if let Some(max) = spec.max {
        track = track.max(max);
    }
    if let Some(color) = &spec.color {
        track = track.color(color);
    }
    if let Some(height) = spec.height {
        track = track.height(height);
    }
    track
}

/// A window track as the options after it ask: its style, its scale and its
/// height.
fn dressed_windows(windows: Vec<crate::Window>, spec: &TrackSpec) -> WindowTrack {
    let mut track = WindowTrack::new(windows).style(
        spec.style
            .and_then(|s| s.window())
            .unwrap_or(WindowStyle::Steps),
    );
    // The band is symmetric about its line, so the top is the one number to
    // pin, and the bottom is as far below.
    if let Some(max) = spec.max {
        let baseline = track.baseline_value();
        track = track.extent(max - baseline);
    }
    if let Some(height) = spec.height {
        track = track.height(height);
    }
    track
}

/// What became of one `--shade` across the panels of a figure, which is one
/// panel unless several places were written.
///
/// Settled once every panel is drawn, since a stretch on one place of a
/// sheet is no fault of the others: `--shade c2:1-100` over `c1` and `c2`
/// shades the second panel and is refused only where no panel is on `c2`.
#[derive(Debug, Default)]
struct ShadeFate {
    /// Whether some panel drew it.
    drawn: bool,
    /// The windows on its sequence it fell outside of.
    outside: Vec<String>,
    /// Every window it was looked for in.
    windows: Vec<String>,
    /// The sequence it is on, for the refusal of one no panel is on.
    sequence: String,
}

/// Where an annotation puts a gene a `--shade` names: the sequence, under the
/// name the figure gives it, and the span. The sequence is `None` for a gene
/// of `--loci`, whose first column names a genome and not a sequence, and
/// which is drawn on the figure's own axis whatever genome it names.
type GenePlace = (Option<String>, u64, u64);

/// The `--shade`s of a command line, as the panels of its figure take them.
#[derive(Debug)]
struct Shading {
    /// What became of each, in the order they were written.
    fates: Vec<ShadeFate>,
    /// Where the annotations put the gene each names, in the same order, and
    /// nothing for one that names no gene. Looked up the first time a panel
    /// shades anything and kept for the rest, since every panel has the same
    /// answer: looked up again for each shade and each panel, four genes over
    /// a sheet of four places went through a 63 MB GFF3 sixteen more times.
    genes: Option<Vec<Vec<GenePlace>>>,
}

impl Shading {
    fn new(invocation: &Invocation) -> Shading {
        Shading {
            fates: invocation
                .shades
                .iter()
                .map(|_| ShadeFate::default())
                .collect(),
            genes: None,
        }
    }
}

/// A window as a reader writes one, `chr1:1-4,000`, for the messages a shade
/// is answered with.
fn written(region: &Region) -> String {
    use crate::track::axis::group_thousands;
    format!(
        "{}:{}-{}",
        region.seq(),
        group_thousands(region.display_start()),
        group_thousands(region.display_end())
    )
}

/// A window of a time axis as its ruler and its tooltips write a time,
/// `year:2010.25-2015.75`: in the units of the table, from the first time to
/// the last, rather than in the thousandths it is drawn at, and never grouped
/// as a count of bases is.
fn written_in_time(region: &Region, decimals: u32) -> String {
    use crate::track::axis::time_text;
    format!(
        "{}:{}-{}",
        region.seq(),
        time_text(region.start(), decimals),
        time_text(region.end().saturating_sub(1), decimals)
    )
}

/// Shades on `figure`, drawn over `region` and written as `over`, every
/// stretch `--shade` asks for that is on its sequence, and says in `shading`
/// what became of each.
///
/// A place on the sequence under any name `--rename` gives it, and a span with
/// no sequence on whatever the figure is drawn over. A gene over its own span
/// as the annotation gives it, without the margin a figure placed on it gets,
/// and over each of its places where it has several; a gene of `--loci` where
/// its row draws it, whichever genome that is. On a continuous time, each is
/// moved into the thousandths the tables are read in, as the place is.
fn shaded(
    figure: Figure,
    invocation: &Invocation,
    files: &mut dyn Files,
    region: &Region,
    over: &Region,
    decimals: u32,
    shading: &mut Shading,
) -> Result<Figure, BuildError> {
    let Some(first) = invocation.shades.first() else {
        return Ok(figure);
    };
    if !figure.takes_shades() {
        return Err(BuildError::Unshaded {
            given: first.given.clone(),
            why: ShadeRefusal::NothingToShade,
        });
    }
    let Shading { fates, genes } = shading;
    if genes.is_none() {
        *genes = Some(gene_places(invocation, files)?);
    }
    let genes = genes.as_deref().unwrap_or_default();
    let aliases = called_by(invocation, region.seq());
    // Moved into thousandths wherever the tables are read in them, as the
    // place is.
    let scale = (decimals > 0 && all_times(invocation)).then(|| 10u64.pow(decimals));
    // A time is said in its own units, as its ruler and its tooltips write
    // it, rather than in the thousandths it is drawn at, and a year is never
    // grouped: a table of whole years said `2,012 to 2,013` under a ruler
    // reading 2012. Only where the ruler counts time, though: a skyline beside
    // a depth is drawn under a ruler of bases, which groups.
    let times = all_times(invocation) && counted(invocation, region, decimals).is_some();
    let timed = |(start, end): (u64, u64)| {
        let (start, end) = match scale {
            Some(scale) => (
                start.saturating_add(1).saturating_mul(scale),
                end.saturating_mul(scale).saturating_add(1),
            ),
            None => (start, end),
        };
        let said = times.then(|| {
            use crate::track::axis::time_text;
            format!(
                "{} to {}",
                time_text(start, decimals),
                time_text(end.saturating_sub(1).max(start), decimals)
            )
        });
        (start, end, said)
    };
    // The window each shade is said to be outside of, or refused against, in
    // the units the ruler counts: on a time, the thousandths a table with
    // fractions is drawn at said `year:2,010,251-2,015,751` under a ruler
    // reading 2011 to 2015.
    let window = if times {
        written_in_time(region, decimals)
    } else {
        written(over)
    };
    let mut figure = figure;
    for ((shading, fate), places) in invocation.shades.iter().zip(fates.iter_mut()).zip(genes) {
        let (spans, sequence): (Vec<(u64, u64, Option<String>)>, String) = match &shading.place {
            ShadePlace::Locus(at) => {
                let here = aliases.contains(&at.seq());
                let spans = if here {
                    vec![timed((at.start(), at.end()))]
                } else {
                    Vec::new()
                };
                (spans, at.seq().to_string())
            }
            ShadePlace::Along(start, end) => {
                (vec![timed((*start, *end))], region.seq().to_string())
            }
            ShadePlace::Gene(_) => {
                let mut sequences: Vec<String> = Vec::new();
                for sequence in places
                    .iter()
                    .filter_map(|(sequence, _, _)| sequence.as_ref())
                {
                    if !sequences.contains(sequence) {
                        sequences.push(sequence.clone());
                    }
                }
                let here = places
                    .iter()
                    .filter(|(sequence, _, _)| {
                        sequence
                            .as_ref()
                            .map_or(true, |sequence| sequence == region.seq())
                    })
                    .map(|(_, start, end)| (*start, *end, None))
                    .collect();
                (here, joined(&sequences))
            }
        };
        fate.windows.push(window.clone());
        if spans.is_empty() {
            fate.sequence = sequence;
            continue;
        }
        let mut seen = false;
        for (start, end, said) in spans {
            let mut shade = crate::Shade::new(start, end);
            if let Some(name) = &shading.name {
                shade = shade.name(name);
            }
            if let Some(said) = said {
                shade = shade.described(said);
            }
            seen |= shade.touches(region);
            figure = figure.shade(shade);
        }
        if seen {
            fate.drawn = true;
        } else {
            fate.outside.push(window.clone());
        }
    }
    Ok(figure)
}

/// Every place the figure's annotations put each gene a `--shade` names, in
/// the order the shades are written, and nothing for one that names no gene.
///
/// Each place once: one gene written as a gene, a transcript and a CDS
/// overlaps itself and is one place, merged as [`place`] merges the place a
/// figure is drawn over. Each annotation is read once, and every name looked
/// up in the one pass over it.
fn gene_places(
    invocation: &Invocation,
    files: &mut dyn Files,
) -> Result<Vec<Vec<GenePlace>>, BuildError> {
    let wanted: Vec<&str> = invocation
        .shades
        .iter()
        .filter_map(|shading| match &shading.place {
            ShadePlace::Gene(name) => Some(name.as_str()),
            _ => None,
        })
        .collect();
    let mut places: Vec<Vec<GenePlace>> = invocation.shades.iter().map(|_| Vec::new()).collect();
    if wanted.is_empty() {
        return Ok(places);
    }
    let mut annotated = false;
    let mut spans: Vec<Vec<GenePlace>> = vec![Vec::new(); wanted.len()];
    let mut names: Vec<String> = Vec::new();
    for spec in &invocation.tracks {
        if !matches!(spec.kind, Kind::Features | Kind::Loci) {
            continue;
        }
        let Some(source) = spec.source.as_ref() else {
            continue;
        };
        // A file that will not open has already been refused by its track.
        let Some(text) = annotation(files, source) else {
            continue;
        };
        annotated = true;
        let found = read::interval::named_each(&text, &wanted);
        // A row of --loci is drawn on the figure's axis whatever genome its
        // first column names, so its gene is shaded there too: held to the
        // sequence of that name, it was refused as elsewhere in the very
        // figure that drew it.
        let loci = spec.kind == Kind::Loci;
        for (spans, found) in spans.iter_mut().zip(found.spans) {
            spans.extend(found.into_iter().map(|(sequence, start, end)| {
                let sequence = (!loci).then(|| renamed(invocation, sequence));
                (sequence, start, end)
            }));
        }
        names.extend(found.names);
    }
    let mut spans = spans.into_iter();
    for (shading, places) in invocation.shades.iter().zip(places.iter_mut()) {
        let ShadePlace::Gene(name) = &shading.place else {
            continue;
        };
        let refuse = |why: ShadeRefusal| BuildError::Unshaded {
            given: shading.given.clone(),
            why,
        };
        if !annotated {
            return Err(refuse(ShadeRefusal::NoAnnotation));
        }
        let mut found = spans.next().unwrap_or_default();
        found.sort();
        for (sequence, start, end) in found {
            match places.last_mut() {
                Some(last) if last.0 == sequence && start <= last.2 => last.2 = last.2.max(end),
                _ => places.push((sequence, start, end)),
            }
        }
        if places.is_empty() {
            return Err(refuse(ShadeRefusal::NoSuchGene {
                near: near_names(&names, name),
            }));
        }
    }
    Ok(places)
}

/// Says what became of each shade once every panel has had it: nothing for
/// one some panel drew, a note for one only outside the windows, and a
/// refusal for one on a sequence no panel is on.
fn settle_shades(
    invocation: &Invocation,
    fates: &[ShadeFate],
    files: &mut dyn Files,
) -> Result<(), BuildError> {
    for (shading, fate) in invocation.shades.iter().zip(fates) {
        if fate.drawn {
            continue;
        }
        if !fate.outside.is_empty() {
            files.note(&format!(
                "--shade {} is outside {}, so it is not drawn",
                shading.given,
                joined(&fate.outside)
            ));
            continue;
        }
        return Err(BuildError::Unshaded {
            given: shading.given.clone(),
            why: ShadeRefusal::Elsewhere {
                sequence: fate.sequence.clone(),
                places: fate.windows.clone(),
            },
        });
    }
    Ok(())
}

/// Shades on a figure across the whole genome, through the offsets its
/// sequences are laid end to end at.
///
/// Its own function because the layout of a whole genome is its own: a place
/// is on one sequence, found by the name the figure gives it, which is the one
/// a `--rename` makes it, or by the name its files give it, and held to the
/// furthest any file reaches on it, which is where that sequence ends in the
/// figure, so a shade never runs into the next one.
/// The tooltip and the alt text name the place as it was written, since the
/// shared axis counts through every sequence before it.
fn shaded_genome(
    figure: Figure,
    invocation: &Invocation,
    genome: &crate::Genome,
    files: &mut dyn Files,
) -> Result<Figure, BuildError> {
    let mut figure = figure;
    for shading in &invocation.shades {
        let refuse = |why: ShadeRefusal| BuildError::Unshaded {
            given: shading.given.clone(),
            why,
        };
        let at = match &shading.place {
            ShadePlace::Locus(at) => at,
            ShadePlace::Along(..) => return Err(refuse(ShadeRefusal::WholeGenome)),
            ShadePlace::Gene(_) => return Err(refuse(ShadeRefusal::GenomeGene)),
        };
        let figured = renamed(invocation, at.seq().to_string());
        let found = std::iter::once(at.seq())
            .chain(std::iter::once(figured.as_str()))
            .chain(
                invocation
                    .renames
                    .iter()
                    .filter(|(_, to)| to == at.seq())
                    .map(|(from, _)| from.as_str()),
            )
            .find_map(|name| {
                genome
                    .sequences()
                    .iter()
                    .find(|sequence| sequence.name == name)
            });
        let Some(sequence) = found else {
            return Err(refuse(ShadeRefusal::Elsewhere {
                sequence: at.seq().to_string(),
                places: genome
                    .sequences()
                    .iter()
                    .map(|sequence| sequence.name.clone())
                    .collect(),
            }));
        };
        let Some(start) = genome.at(&sequence.name, at.start()) else {
            files.note(&format!(
                "--shade {} is past the furthest any file reaches on {}, so it is not drawn",
                shading.given, sequence.name
            ));
            continue;
        };
        let offset = start - at.start();
        let end = offset.saturating_add(at.end().min(sequence.length));
        let mut shade = crate::Shade::new(start, end).described(written(at));
        if let Some(name) = &shading.name {
            shade = shade.name(name);
        }
        figure = figure.shade(shade);
    }
    Ok(figure)
}

/// Files that name the same sequences, and those sequences, each with how
/// many rows the files give it.
type Naming = (Vec<String>, Vec<(String, usize)>);

/// Where a place named by a word is, and whether it was a gene.
struct Placed {
    region: Region,
    /// The gene's name as its annotation spells it, where the place is a gene.
    gene: Option<String>,
    /// Where the gene is, without the margin the figure is drawn with,
    /// 0-based and half-open: every row going by its name, merged.
    extent: Option<(u64, u64)>,
    /// The file whose rows reach furthest, where no file says how long the
    /// sequence is and the figure ends where the rows do.
    reached: Option<String>,
}

/// The name the figure gives a sequence a file calls `name`: the one
/// `--rename` gives it, or its own.
fn renamed(invocation: &Invocation, name: String) -> String {
    invocation
        .renames
        .iter()
        .find(|(from, _)| *from == name)
        .map_or(name, |(_, to)| to.clone())
}

/// The names the files may call the figure's sequence `name` by: its own, and
/// every one `--rename` makes it.
fn called_by<'a>(invocation: &'a Invocation, name: &'a str) -> Vec<&'a str> {
    std::iter::once(name)
        .chain(
            invocation
                .renames
                .iter()
                .filter(|(_, to)| to == name)
                .map(|(from, _)| from.as_str()),
        )
        .collect()
}

/// Finds a place named by a word in the figure's own files.
///
/// A sequence first: a word that names one is that sequence, whole, as long
/// as a file says it is (a FASTA record, a BAM header, a GFF3's
/// `##sequence-region`, a VCF's `##contig`), or as far as any file reaches on
/// it. Then a gene, looked up in every annotation of the figure by the names
/// GFF3, GTF and BED give it, and drawn with a margin of a tenth of its length
/// either side, a hundred bases at least, so its ends are on the page. One
/// gene written at several levels, gene, transcript and CDS, is one place; a
/// name at two places is refused with both, and a name at none with the
/// nearest names the annotation has.
fn place(name: &str, invocation: &Invocation, files: &mut dyn Files) -> Result<Placed, BuildError> {
    let surveyed = survey(name, invocation, files, Asked::Place)?;
    located(name, invocation, files, surveyed)
}

/// What [`survey`] reads a figure's files for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Asked {
    /// Where a word is: a sequence as long as a file says, a gene, or a
    /// sequence as far as the rows reach.
    Place,
    /// How long the files say a sequence is, and nothing more, which is all
    /// a span written from base 1 is checked against. No gene is looked up
    /// and no file is read whole for how far its rows reach: a VCF read
    /// through its index with no `##contig` was read whole for its reach,
    /// which was then thrown away, so a circle of one sequence cost as much
    /// as every row the file holds on the others.
    Lengths,
}

/// A sequence's length as one file states it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Stated {
    /// The sequence, under the name the figure gives it.
    sequence: String,
    length: u64,
    /// The file that says so, as the figure calls it.
    file: String,
}

/// What a figure's files say of a word, each read once, for [`located`].
struct Survey<'a> {
    /// Every length a file states, in the order of the tracks.
    lengths: Vec<Stated>,
    /// Each row that goes by the word as a gene's name: its sequence and
    /// where it is.
    spans: Vec<(String, u64, u64)>,
    /// The names the annotations give, for a word that is none of them.
    names: Vec<String>,
    /// How far each file's rows reach, in the order of the tracks, or the
    /// file to read for it once nothing else has placed the figure.
    reaches: Vec<Reach<'a>>,
    /// Whether any annotation was read for a gene.
    annotated: bool,
    /// The gene's name as its annotation spells it.
    spelled: Option<String>,
}

/// Each sequence `source` gives the length of, under the name the figure
/// gives it, with the file.
fn stated_by(
    invocation: &Invocation,
    source: &Source,
    held: impl IntoIterator<Item = (String, u64)>,
) -> Vec<Stated> {
    let file = called(source);
    held.into_iter()
        .map(|(sequence, length)| Stated {
            sequence: renamed(invocation, sequence),
            length,
            file: file.clone(),
        })
        .collect()
}

/// Reads the figure's files for what they say of `name`, as far as `asked`
/// needs.
fn survey<'a>(
    name: &str,
    invocation: &'a Invocation,
    files: &mut dyn Files,
    asked: Asked,
) -> Result<Survey<'a>, BuildError> {
    let placing = asked == Asked::Place;
    let mut found = Survey {
        lengths: Vec::new(),
        spans: Vec::new(),
        names: Vec::new(),
        reaches: Vec::new(),
        annotated: false,
        spelled: None,
    };
    let aliases = called_by(invocation, name);
    let open_error = |spec: &TrackSpec, source: &Source, cause: io::Error| BuildError::Open {
        track: spec.kind.flag(),
        path: called(source),
        cause,
    };
    for spec in &invocation.tracks {
        let sources = [spec.source.as_ref(), spec.second.as_ref()];
        for source in sources.into_iter().flatten() {
            let own = spec.source.as_ref() == Some(source);
            // A BCF names its sequences, and the lengths of some, in its
            // header, and its records are read for how far they reach only
            // where nothing else places the figure, as a file read through
            // its index is.
            if let Some((header, _)) = own.then(|| bcf(files, source)).flatten() {
                let held = header
                    .contigs
                    .into_iter()
                    .filter_map(|(sequence, length)| Some((sequence, length?)));
                found.lengths.extend(stated_by(invocation, source, held));
                if placing {
                    found.reaches.push(Reach::Later(spec, source));
                }
                continue;
            }
            let text = match files
                .sequences(source)
                .map_err(|cause| open_error(spec, source, cause))?
            {
                Some(held) => {
                    found.lengths.extend(stated_by(invocation, source, held));
                    // A bigBed names its sequences in its index and its genes
                    // in its rows, which are read whole only where the name
                    // is no sequence any file has named so far.
                    let gene = placing
                        && matches!(spec.kind, Kind::Features | Kind::Loci)
                        && own
                        && !found.lengths.iter().any(|stated| stated.sequence == name);
                    match gene.then(|| annotation(files, source)).flatten() {
                        Some(rows) => rows,
                        None => continue,
                    }
                }
                // A file read through its index is read here for its header
                // alone, which says how long its sequences are where it says
                // at all, and its rows are read whole, for how far they reach,
                // only where nothing else places the figure. Read whole here
                // and again by its track, the calls beside the annotation a
                // gene is found in cost 8.9 s and 1.3 GB for one gene of a
                // VCF of 825 MB, which now takes 7 ms. An index that is not
                // trusted is said by the track, as `slurp` reads it. An
                // annotation is read whole for a gene's name, and for its
                // header alone where no gene is looked up.
                None if own
                    && through_index(spec.kind)
                    && (spec.kind != Kind::Features || !placing) =>
                {
                    match indexed(files, source) {
                        Ok(opened) => {
                            found.lengths.extend(stated_by(
                                invocation,
                                source,
                                sequence_lengths(&opened.head.text),
                            ));
                            if placing {
                                found.reaches.push(Reach::Later(spec, source));
                            }
                            continue;
                        }
                        Err(_) => match files.text(source) {
                            Ok(text) => text,
                            Err(_) => continue,
                        },
                    }
                }
                // A file that will not open is its track's to report.
                None => match files.text(source) {
                    Ok(text) => text,
                    Err(_) => continue,
                },
            };
            found
                .lengths
                .extend(stated_by(invocation, source, sequence_lengths(&text)));
            // A PAF writes the length of every query it aligns, which is the
            // sequence a synteny figure or a dot plot is drawn along. It was
            // not asked, so a figure placed on its own query was refused.
            if matches!(spec.kind, Kind::Synteny | Kind::Dotplot) && own {
                found
                    .lengths
                    .extend(stated_by(invocation, source, paf_query_lengths(&text)));
            }
            if !placing {
                continue;
            }
            if matches!(spec.kind, Kind::Features | Kind::Loci) && own {
                found.annotated = true;
                let named = read::interval::named(&text, name);
                found.spelled = found.spelled.or(named.spelled);
                found.spans.extend(
                    named
                        .spans
                        .into_iter()
                        .map(|(sequence, start, end)| (renamed(invocation, sequence), start, end)),
                );
                found.names.extend(named.names);
            }
            found.reaches.push(Reach::Read(
                reached_by(spec.kind, &text, &aliases),
                called(source),
            ));
        }
    }
    Ok(found)
}

/// Where a word is, from what the figure's files say of it: [`place`], once
/// they are read.
fn located(
    name: &str,
    invocation: &Invocation,
    files: &mut dyn Files,
    surveyed: Survey<'_>,
) -> Result<Placed, BuildError> {
    let Survey {
        lengths,
        mut spans,
        names,
        reaches,
        annotated,
        spelled,
    } = surveyed;
    let aliases = called_by(invocation, name);
    if let Some(stated) = lengths.iter().find(|stated| stated.sequence == name) {
        let region = Region::new(name, 0, stated.length.max(1))
            .map_err(|_| nowhere(name, invocation, files, &names, annotated))?;
        return Ok(Placed {
            region,
            gene: None,
            extent: None,
            reached: None,
        });
    }

    // One gene written at several levels overlaps itself: merged, it is one
    // place. What is left over is as many places as the name has.
    spans.sort();
    let mut places: Vec<(String, u64, u64)> = Vec::new();
    for (sequence, start, end) in spans {
        match places.last_mut() {
            Some(last) if last.0 == sequence && start <= last.2 => last.2 = last.2.max(end),
            _ => places.push((sequence, start, end)),
        }
    }
    match places.as_slice() {
        [(sequence, start, end)] => {
            let margin = ((end - start) / 10).max(100);
            let length = lengths
                .iter()
                .find(|stated| stated.sequence == *sequence)
                .map(|stated| stated.length);
            let stop = (end + margin).min(length.unwrap_or(u64::MAX));
            let region = Region::new(sequence, start.saturating_sub(margin), stop.max(start + 1))
                .map_err(|_| nowhere(name, invocation, files, &names, annotated))?;
            return Ok(Placed {
                region,
                reached: None,
                gene: Some(spelled.unwrap_or_else(|| name.to_string())),
                extent: Some((*start, *end)),
            });
        }
        [] => {}
        many => {
            return Err(BuildError::Several {
                name: name.to_string(),
                places: many
                    .iter()
                    .map(|(sequence, start, end)| {
                        Region::new(sequence, *start, *end)
                            .map_or_else(|_| sequence.clone(), |region| region.to_string())
                    })
                    .collect(),
            })
        }
    }

    // The file whose rows reach furthest, the earlier of two that reach as
    // far.
    let mut furthest: Option<(u64, String)> = None;
    for reach in reaches {
        let (reach, file) = match reach {
            Reach::Read(reach, file) => (reach, file),
            Reach::Later(spec, source) => match whole_text(files, source) {
                Ok(text) => (reached_by(spec.kind, &text, &aliases), called(source)),
                Err(_) => continue,
            },
        };
        if let Some(reach) = reach {
            if furthest.as_ref().map_or(true, |(known, _)| reach > *known) {
                furthest = Some((reach, file));
            }
        }
    }
    if let Some((end, file)) = furthest {
        if let Ok(region) = Region::new(name, 0, end.max(1)) {
            return Ok(Placed {
                region,
                gene: None,
                extent: None,
                reached: Some(file),
            });
        }
    }
    Err(nowhere(name, invocation, files, &names, annotated))
}

/// How far a file's rows reach on a sequence, or where to find out.
enum Reach<'a> {
    /// As far as this, in the file called this, where they reach it at all.
    Read(Option<u64>, String),
    /// A file read through its index, whose rows are read whole for this only
    /// where nothing else places the figure.
    Later(&'a TrackSpec, &'a Source),
}

/// How far the rows of a track's file reach on the sequence, under any name
/// the figure gives it.
fn reached_by(kind: Kind, text: &str, aliases: &[&str]) -> Option<u64> {
    aliases
        .iter()
        .filter_map(|alias| {
            if kind == Kind::Manhattan {
                read::point::association_table(
                    text,
                    &Region::new(*alias, 0, 1 << 28)
                        .unwrap_or_else(|_| Region::new("x", 0, 1).expect("a window")),
                )
                .ok()
                .and_then(|table| table.points.iter().map(|point| point.pos + 1).max())
            } else {
                sequence_column(kind).and_then(|column| reach_on(text, &column, alias))
            }
        })
        .max()
}

/// The largest bigBed read whole to look a gene up by its name. A bigBed is
/// read a window at a time, and only a name has to be found in every row: one
/// of 29.6 MB holding 4.7 million rows of BED6 was read whole in 1.9 s, as
/// 187 MB of BED, so this bound is about four seconds and 400 MB of rows.
/// Many sequences cost no more, since each row's is looked up by its number:
/// one of 20.9 MB naming 200,000 sequences places a figure in 0.9 s, where
/// looking each row's sequence up among them all took 16 s.
const WHOLE_BIGBED: u64 = 64 << 20;

/// The rows of an annotation a gene is looked up in by its name: a file's
/// text, or every row of a bigBed as BED. `None` for a file that will not
/// open, which its track reports, and for a bigBed larger than
/// [`WHOLE_BIGBED`], which is said.
fn annotation(files: &mut dyn Files, source: &Source) -> Option<String> {
    annotation_within(files, source, WHOLE_BIGBED)
}

/// [`annotation`], reading a bigBed whole only up to `most` bytes of it.
fn annotation_within(files: &mut dyn Files, source: &Source, most: u64) -> Option<String> {
    match native(files, source) {
        Ok(Some((Binary::BigBed, mut file))) => {
            let size = file.seek(io::SeekFrom::End(0)).ok()?;
            if size > most {
                files.note(&format!(
                    "{} is a bigBed of {} MB, which is not read whole to find a gene by its \
                     name; write the gene's place, as chr1:1,001-2,000",
                    called(source),
                    size.div_ceil(1 << 20)
                ));
                return None;
            }
            read::bigbed::bed(file, None).ok()
        }
        Ok(Some(_)) => None,
        _ => files.text(source).ok(),
    }
}

/// The codon ruler `--codons` asks for.
///
/// Over the coding sequence of the gene the figure is placed on, or of the one
/// gene that codes in the place written, as the figure's annotations write it:
/// each transcript read for what it codes, and the one stretch every
/// transcript of the gene agrees on. The place written and not the window
/// drawn, so a page that moves the window keeps counting the gene it started
/// on, as a figure placed by a gene's name keeps its title wherever it is
/// moved to. Translated with the table `--genetic-code` names, or the one the
/// CDS names, or the standard one, from the figure's reference where it has
/// one, and numbers alone where it has none.
fn codons(
    spec: &TrackSpec,
    invocation: &Invocation,
    placed: Option<&Placed>,
    region: &Region,
    reference: Option<&Source>,
    files: &mut dyn Files,
) -> Result<Box<dyn Track>, BuildError> {
    /// One transcript that codes: its rows where they fall on its exons, the
    /// stretches they make, the annotation it is in, the text it was read
    /// from, and the name of the sequence it is on there.
    struct Model {
        feature: crate::Feature,
        rows: Vec<read::interval::CdsRow>,
        pieces: Vec<(u64, u64)>,
        file: usize,
        text: usize,
        sequence: String,
    }
    let refuse = BuildError::Uncounted;
    // A gene placed by its name is looked for inside its own rows, and a
    // place written is looked through for the one gene that codes in it.
    let (sequence, start, end, gene, place) = match (placed, &invocation.region) {
        (Some(placed), _) => match (&placed.gene, placed.extent) {
            (Some(gene), Some((start, end))) => (
                placed.region.seq().to_string(),
                start,
                end,
                Some(gene.clone()),
                gene.clone(),
            ),
            _ => {
                return Err(refuse(CodonRefusal::WholeSequence {
                    sequence: placed.region.seq().to_string(),
                }))
            }
        },
        (None, Some(written)) => (
            written.seq().to_string(),
            written.start(),
            written.end(),
            None,
            written.to_string(),
        ),
        (None, None) => {
            return Err(BuildError::Placeless(
                crate::cli::args::ArgError::CodonsWithoutPlace,
            ))
        }
    };

    // The texts the models were read from: an annotation's window, read
    // through the index beside it, or the whole of it.
    let mut texts: Vec<String> = Vec::new();
    let mut annotations: Vec<String> = Vec::new();
    let mut models: Vec<Model> = Vec::new();
    let mut unread: Option<&Source> = None;
    for track in &invocation.tracks {
        if track.kind != Kind::Features {
            continue;
        }
        let Some(source) = track.source.as_ref() else {
            continue;
        };
        let path = called(source);
        // A `.bed` named on its own that turns out to be a signal is no
        // annotation, and a file that will not open is said below, where the
        // ruler would otherwise say the figure had no annotation at all.
        let signal = |text: &str| track.guessed && refine(track.kind, text).is_some();
        let mut whole: Option<usize> = None;
        let mut annotated = false;
        for alias in called_by(invocation, &sequence) {
            let Ok(over) = Region::new(alias, start, end.max(start + 1)) else {
                continue;
            };
            // Through the index beside the file where its own track reads its
            // window through one: a place of a few hundred bases is a few
            // blocks of an annotation of a whole genome, and read whole, it
            // cost the time and the memory of every row of every gene to
            // number a hundred codons. A row the window's reader refuses is
            // read again whole, so the refusal names its line in the file,
            // and so is a window holding a CDS under nothing, whose other
            // pieces need not be over it.
            let mut found = None;
            if let Ok(window) = windowed(track, &over, files, source) {
                if signal(window.probe.as_deref().unwrap_or(&window.text)) {
                    break;
                }
                let read = read::interval::coding(&window.text, &over, track.format);
                if let (Ok(read), false) = (read, read::interval::pieced(&window.text, alias)) {
                    texts.push(window.text);
                    found = Some((read, texts.len() - 1));
                }
            }
            let (read, text) = match found {
                Some(found) => found,
                None => {
                    let text = match whole {
                        Some(text) => text,
                        None => {
                            let Some(text) = annotation(files, source) else {
                                unread = unread.or(Some(source));
                                break;
                            };
                            if signal(&text) {
                                break;
                            }
                            texts.push(text);
                            *whole.insert(texts.len() - 1)
                        }
                    };
                    // A broken row is said as its own track says it, since
                    // the ruler may be built first and a gene with no CDS is
                    // not what is wrong.
                    let read = wrap(
                        Kind::Features.flag(),
                        &path,
                        read::interval::coding(&texts[text], &over, track.format),
                    )?;
                    (read, text)
                }
            };
            annotated = true;
            models.extend(read.into_iter().map(|cds| {
                let rows = coded(&cds);
                Model {
                    pieces: stretches(&rows),
                    rows,
                    feature: cds.feature,
                    file: annotations.len(),
                    text,
                    sequence: alias.to_string(),
                }
            }));
        }
        if annotated {
            annotations.push(path);
        }
    }
    if annotations.is_empty() {
        if let Some(source) = unread {
            fetch(Kind::Features.flag(), source, files)?;
        }
        return Err(refuse(CodonRefusal::NoAnnotation));
    }

    // A gene's own rows hold every transcript of it, so a transcript coding
    // outside them is another gene's; a place written holds whatever codes
    // in it, if only a codon.
    let mut chosen: Vec<&Model> = models
        .iter()
        .filter(|model| match gene {
            Some(_) => {
                !model.pieces.is_empty()
                    && model
                        .pieces
                        .iter()
                        .all(|&(from, to)| from >= start && to <= end)
            }
            None => model
                .pieces
                .iter()
                .any(|&(from, to)| from < end && to > start),
        })
        .collect();
    let agree = |models: &[&Model]| {
        models.iter().all(|model| {
            model.pieces == models[0].pieces && model.feature.strand == models[0].feature.strand
        })
    };
    // A gene in the intron of the one asked for is inside its rows too, and
    // the name tells the two apart.
    if let Some(gene) = &gene {
        if !agree(&chosen) {
            let named: Vec<&Model> = chosen
                .iter()
                .copied()
                .filter(|model| {
                    [&model.feature.name, &model.feature.gene]
                        .into_iter()
                        .flatten()
                        .any(|name| name.eq_ignore_ascii_case(gene))
                })
                .collect();
            if !named.is_empty() {
                chosen = named;
            }
        }
    }
    let called_as = |model: &Model| {
        model
            .feature
            .gene
            .clone()
            .or_else(|| model.feature.name.clone())
            .unwrap_or_else(|| {
                Region::new(&model.sequence, model.feature.start, model.feature.end)
                    .map_or_else(|_| model.sequence.clone(), |at| at.to_string())
            })
    };
    let Some(first) = chosen.first().copied() else {
        // A place written over one gene that the annotation writes no CDS
        // for is that gene's want of one, which is the thing to say. Every
        // model read touches the place; one that codes elsewhere, as a gene
        // whose untranslated end is all the place holds, has a CDS, and the
        // place is merely not on it.
        let mut uncoded: Vec<String> = Vec::new();
        for model in &models {
            let own = called_as(model);
            if !uncoded.contains(&own) {
                uncoded.push(own);
            }
        }
        if models.iter().any(|model| !model.pieces.is_empty()) {
            uncoded.clear();
        }
        return Err(refuse(match (gene, uncoded.as_slice()) {
            (Some(gene), _) => CodonRefusal::NoCds {
                gene,
                files: annotations,
            },
            (None, [only]) => CodonRefusal::NoCds {
                gene: only.clone(),
                files: annotations,
            },
            (None, _) => CodonRefusal::NoneCodes {
                place,
                files: annotations,
            },
        }));
    };
    if !agree(&chosen) {
        let mut genes: Vec<String> = Vec::new();
        let mut transcripts: Vec<String> = Vec::new();
        for model in &chosen {
            let own = called_as(model);
            if !genes.contains(&own) {
                genes.push(own);
            }
            if let Some(name) = &model.feature.name {
                if !transcripts.contains(name) {
                    transcripts.push(name.clone());
                }
            }
        }
        return Err(refuse(if genes.len() > 1 {
            CodonRefusal::Several { genes, place }
        } else {
            CodonRefusal::Isoforms {
                gene: gene.unwrap_or_else(|| genes.swap_remove(0)),
                transcripts,
            }
        }));
    }
    let name = gene.unwrap_or_else(|| called_as(first));
    if first.pieces.len() > 1 {
        return Err(refuse(CodonRefusal::Spliced {
            gene: name,
            pieces: first.pieces.len(),
        }));
    }
    let strand = first.feature.strand;
    if !matches!(strand, crate::Strand::Forward | crate::Strand::Reverse) {
        return Err(refuse(CodonRefusal::NoStrand { gene: name }));
    }
    if let Some(at) = frameshift(&first.rows, strand) {
        return Err(refuse(CodonRefusal::Frameshift { gene: name, at }));
    }
    // Codon 1 is the start codon only where the CDS begins on one, which its
    // 5'-most row says: the first of them forwards and the last backwards.
    let five = if strand == crate::Strand::Reverse {
        first.rows.last()
    } else {
        first.rows.first()
    };
    if let Some(five) = five {
        let phase = five.phase.filter(|phase| *phase > 0);
        if phase.is_some() || five.partial {
            return Err(refuse(CodonRefusal::Partial { gene: name, phase }));
        }
    }
    let (from, to) = first.pieces[0];

    // The table the flag names, or else the one the CDS names, and the
    // standard one where neither names any.
    let named = read::interval::translation_table(&texts[first.text], &first.sequence, from, to);
    let table = match (spec.genetic_code, named) {
        (Some(asked), Some(named)) if asked != named => {
            files.note(&format!(
                "--genetic-code {asked} reads {name}, whose CDS in {} names table {named}",
                annotations[first.file]
            ));
            asked
        }
        (Some(asked), _) => asked,
        (None, Some(named)) => named,
        (None, None) => 1,
    };
    let Some(residues) = crate::track::codon::ncbi_table(table) else {
        return Err(refuse(CodonRefusal::UnknownTable { gene: name, table }));
    };
    let left = (to - from) % 3;
    if left > 0 {
        files.note(&format!(
            "the last {} of the CDS of {name} {} not a whole codon, and {} left off the ruler",
            if left == 1 { "base" } else { "2 bases" },
            if left == 1 { "is" } else { "are" },
            if left == 1 { "is" } else { "are" },
        ));
    }

    let mut ruler = crate::CodonTrack::new(from, to, strand)
        .genetic_code(residues)
        .label(spec.label.clone().unwrap_or_else(|| name.clone()));
    if let Some(color) = &spec.color {
        ruler = ruler.color(color);
    }
    // The letters of a codon half in the window are worth reading, so the
    // bases are read two past either edge. A reference that will not give
    // them leaves the ruler numbers alone, and its own track says why.
    if let Some(source) = reference {
        for alias in called_by(invocation, region.seq()) {
            let Ok(wide) = Region::new(
                alias,
                region.start().saturating_sub(2),
                region.end().saturating_add(2),
            ) else {
                continue;
            };
            let read = second_sequence(Kind::Sequence.flag(), source, &wide, files)
                .and_then(|reference| reference.clip(&wide));
            if let Ok((at, bases)) = read {
                ruler = ruler.sequence(at, bases);
                break;
            }
        }
    }
    Ok(Box::new(ruler))
}

/// The rows of a transcript's CDS where they fall on its exons, as they are
/// translated, which for a BED12 is its thick span cut by its blocks, in the
/// order of their starts. A piece whose 5' end is not its row's, cut out of
/// the row by an intron, keeps neither the row's phase nor its word that the
/// CDS goes on past it, which are said of the row's 5' end.
fn coded(cds: &read::interval::Cds) -> Vec<read::interval::CdsRow> {
    let exons = &cds.feature.exons;
    if exons.is_empty() {
        return cds.rows.clone();
    }
    let reverse = cds.feature.strand == crate::Strand::Reverse;
    let mut pieces = Vec::new();
    for row in &cds.rows {
        for &(start, end) in exons {
            let (lo, hi) = (row.start.max(start), row.end.min(end));
            if lo >= hi {
                continue;
            }
            let five = if reverse {
                hi == row.end
            } else {
                lo == row.start
            };
            pieces.push(read::interval::CdsRow {
                start: lo,
                end: hi,
                phase: row.phase.filter(|_| five),
                partial: row.partial && five,
            });
        }
    }
    pieces.sort_unstable_by_key(|row| (row.start, row.end));
    pieces
}

/// The stretches the rows of a CDS cover, in order and apart: the rows
/// joined where they touch or overlap, which is how far they reach and not
/// how they are read.
fn stretches(rows: &[read::interval::CdsRow]) -> Vec<(u64, u64)> {
    let mut joined: Vec<(u64, u64)> = Vec::new();
    for row in rows {
        match joined.last_mut() {
            Some(last) if row.start <= last.1 => last.1 = last.1.max(row.end),
            _ => joined.push((row.start, row.end)),
        }
    }
    joined
}

/// Where the rows of a CDS, read from its 5' end, leave the frame they
/// began in, 0-based: the first base read of a row that overlaps the one
/// before it, as NCBI writes a ribosomal slippage, two rows sharing a base,
/// or that meets it with a phase that does not carry its frame on. `None`
/// for rows read in one frame from end to end.
///
/// Joined into one stretch, as a drawing joins them, two rows that share a
/// base count it once where the ribosome reads it twice, and every codon
/// past the slip was numbered in the frame the ribosome left.
fn frameshift(rows: &[read::interval::CdsRow], strand: crate::Strand) -> Option<u64> {
    let reverse = strand == crate::Strand::Reverse;
    let mut order: Vec<&read::interval::CdsRow> = rows.iter().collect();
    if reverse {
        order.reverse();
    }
    let mut read = 0u64;
    let mut last: Option<&read::interval::CdsRow> = None;
    for row in order {
        if let Some(last) = last {
            let (overlaps, meets) = if reverse {
                (row.end > last.start, row.end == last.start)
            } else {
                (row.start < last.end, row.start == last.end)
            };
            // The bases of this row that finish the codon the rows before it
            // left open.
            let carried = (3 - read % 3) % 3;
            let shifts = meets && row.phase.is_some_and(|phase| u64::from(phase) != carried);
            if overlaps || shifts {
                return Some(if reverse { row.end - 1 } else { row.start });
            }
        }
        read += row.end - row.start;
        last = Some(row);
    }
    None
}

/// The refusal of a name no file of the figure has, with the genes named
/// nearly the same and the sequences each file does name.
fn nowhere(
    name: &str,
    invocation: &Invocation,
    files: &mut dyn Files,
    names: &[String],
    annotated: bool,
) -> BuildError {
    // Read again, and only here, as the empty window's message is: a figure
    // that finds its place has no need to know.
    let mut held: Vec<Naming> = Vec::new();
    for spec in &invocation.tracks {
        let Some(source) = spec.source.as_ref() else {
            continue;
        };
        // A file read through its index names its sequences in the index,
        // and the lengths of some in its header, so neither is read whole.
        let named = through_index(spec.kind)
            .then(|| indexed(files, source).ok())
            .flatten()
            .map(|found| {
                let mut sequences: Vec<(String, usize)> = sequence_lengths(&found.head.text)
                    .into_iter()
                    .map(|(name, _)| (name, 0))
                    .collect();
                sequences.extend(found.index.names().iter().map(|name| (name.clone(), 0)));
                sequences
            });
        // A BCF names every sequence in its header, whatever its index says.
        let named = match bcf(files, source) {
            Some((header, _)) => Some(
                header
                    .contigs
                    .into_iter()
                    .map(|(name, _)| (name, 0))
                    .collect(),
            ),
            None => named,
        };
        let mut sequences: Vec<(String, usize)> = match (files.sequences(source), named) {
            (Ok(Some(lengths)), _) => lengths.into_iter().map(|(name, _)| (name, 0)).collect(),
            (_, Some(named)) => named,
            _ => match files.text(source) {
                Ok(text) => {
                    let mut sequences: Vec<(String, usize)> = sequence_lengths(&text)
                        .into_iter()
                        .map(|(name, _)| (name, 0))
                        .collect();
                    sequences.extend(rows_on(spec.kind, &text));
                    sequences
                }
                Err(_) => continue,
            },
        };
        let mut seen: Vec<String> = Vec::new();
        sequences.retain(|(name, _)| {
            let first = !seen.contains(name);
            seen.push(name.clone());
            first
        });
        if sequences.is_empty() {
            continue;
        }
        let same = held.iter().position(|(_, known)| {
            known.len() == sequences.len() && known.iter().zip(&sequences).all(|(a, b)| a.0 == b.0)
        });
        match same {
            Some(at) => held[at].0.push(called(source)),
            None if held.len() < 4 => held.push((vec![called(source)], sequences)),
            None => {}
        }
    }
    let near = near_names(names, name);
    // A name the files hold under another spelling is most likely that one.
    // Where no annotation is there to look genes up in, a word that is no
    // sequence is most likely the one sequence the files do name.
    let mut every: Vec<(String, usize)> = Vec::new();
    for (_, sequences) in &held {
        for (sequence, count) in sequences {
            if !every.iter().any(|(known, _)| known == sequence) {
                every.push((sequence.clone(), *count));
            }
        }
    }
    let rename = if near.is_empty() {
        rename_for(&every, name, !annotated).map(Box::new)
    } else {
        None
    };
    BuildError::Nowhere {
        name: name.to_string(),
        rename,
        held,
        near,
        annotated,
    }
}

/// The three of `names` nearest `name`: within two edits, or the same
/// letters in another case, closest first.
fn near_names(names: &[String], name: &str) -> Vec<String> {
    let mut near: Vec<(usize, &String)> = names
        .iter()
        .filter_map(|candidate| {
            let distance = if candidate.eq_ignore_ascii_case(name) {
                0
            } else {
                crate::cli::args::edits(&candidate.to_ascii_lowercase(), &name.to_ascii_lowercase())
            };
            (distance <= 2).then_some((distance, candidate))
        })
        .collect();
    near.sort();
    near.dedup_by(|a, b| a.1 == b.1);
    near.into_iter()
        .take(3)
        .map(|(_, candidate)| candidate.clone())
        .collect()
}

/// The sequences a text file says the lengths of: FASTA records, a SAM
/// header's `@SQ` lines, a GFF3's `##sequence-region` and a VCF's `##contig`.
fn sequence_lengths(text: &str) -> Vec<(String, u64)> {
    let mut found = Vec::new();
    if text.starts_with('>') {
        if let Ok(records) = read::seq::fasta(text) {
            for (name, bases) in records {
                found.push((name, bases.len() as u64));
            }
        }
        return found;
    }
    for line in text
        .lines()
        .take_while(|line| line.starts_with('#') || line.starts_with('@'))
    {
        if let Some(rest) = line.strip_prefix("##sequence-region") {
            let fields: Vec<&str> = rest.split_whitespace().collect();
            if let [name, _, end] = fields.as_slice() {
                if let Ok(end) = end.parse() {
                    found.push((name.to_string(), end));
                }
            }
        } else if let Some(rest) = line.strip_prefix("##contig=<") {
            let field = |key: &str| {
                rest.trim_end_matches('>')
                    .split(',')
                    .find_map(|pair| pair.strip_prefix(key))
                    .map(str::to_string)
            };
            if let (Some(name), Some(length)) = (field("ID="), field("length=")) {
                if let Ok(length) = length.parse() {
                    found.push((name, length));
                }
            }
        } else if line.starts_with("@SQ") {
            let field = |key: &str| line.split('\t').find_map(|pair| pair.strip_prefix(key));
            if let (Some(name), Some(length)) = (field("SN:"), field("LN:")) {
                if let Ok(length) = length.parse() {
                    found.push((name.to_string(), length));
                }
            }
        }
    }
    found
}

/// The query sequences a PAF aligns, each with the length its second column
/// gives it, once each.
fn paf_query_lengths(text: &str) -> Vec<(String, u64)> {
    let mut found: Vec<(String, u64)> = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let (Some(name), Some(length)) = (fields.next(), fields.next()) else {
            continue;
        };
        let Ok(length) = length.parse() else {
            continue;
        };
        if !name.is_empty() && !found.iter().any(|(held, _)| held == name) {
            found.push((name.to_string(), length));
        }
    }
    found
}

/// How far a file's rows on `name` reach, by the largest number in the
/// columns that hold positions.
fn reach_on(text: &str, column: &SequenceColumn, name: &str) -> Option<u64> {
    let mut furthest: Option<u64> = None;
    for (_, line) in read::lines(text) {
        let cols = read::columns(line);
        if cols.len() < column.width || cols[column.name].trim() != name {
            continue;
        }
        for at in column
            .position
            .iter()
            .chain(std::iter::once(&(column.position[0] + 1)))
        {
            if let Some(value) = cols
                .get(*at)
                .and_then(|field| field.trim().parse::<u64>().ok())
            {
                furthest = furthest.max(Some(value));
            }
        }
    }
    furthest
}

/// What every track of one figure is built against.
struct Context<'a> {
    region: &'a Region,
    theme: &'a Theme,
    /// The FASTA a `--sequence` track reads, if the figure has one.
    reference: Option<&'a Source>,
    /// The places the figure's times are read to: nought for whole units,
    /// or `read::series::DECIMALS` where a table has fractions of one.
    decimals: u32,
    /// Whether a `--genotypes` track of the figure reads the file this track
    /// does, so a VCF drawn as its calls need not say its samples can be.
    genotyped: bool,
    /// The colours `--colors` chose for the values of a sheet's columns,
    /// which every sheet of the figure paints alike.
    colors: &'a [(String, Vec<(String, String)>)],
    /// How wide the figure is drawn, in pixels, which says how many bases a
    /// pixel holds and so which zoom level of a bigWig to read.
    width: f64,
}

/// Adds a track's keys to the figure's, each once: a lineage coloured beside
/// a tree and beside a matrix is one key, since it is one colour.
fn gather(into: &mut crate::track::legend::Legend, from: &crate::track::legend::Legend) {
    *into = std::mem::take(into).and(from);
}

/// Asks the caller for a source's text, and names it for any error message.
///
/// The whole of this module's contact with the outside world used to be here,
/// as an `fs::read_to_string` and a read of standard input. It is a closure now
/// because the two callers that matter want different answers: a shell opens
/// the path, and a browser looks the name up in whatever the editor is holding.
/// Neither is more correct, and the grammar does not care, so the grammar stops
/// deciding.
/// Reads the sheet a track's `--traits` names, or nothing where it named none.
fn sheet(
    spec: &TrackSpec,
    files: &mut dyn Files,
) -> Result<Option<(read::sheet::Sheet, String)>, BuildError> {
    let Some(source) = spec.traits.as_ref() else {
        return Ok(None);
    };
    let name = spec.kind.flag();
    let (text, path) = fetch(name, source, files)?;
    let held = wrap(name, &path, read::sheet::sheet(&text))?;
    Ok(Some((held, path)))
}

/// The metadata columns a sheet becomes once it is joined to a track's rows.
///
/// The join is names and it is checked here rather than left to the drawing,
/// because a sheet that names none of these rows draws a strip of empty
/// outlines beside every one of them, and a figure that says "nothing is known
/// about any of these" looks exactly like a figure that read the wrong file.
///
/// The colours `--colors` chose go to every column of the sheet they name,
/// drawn as a strip or not, so a phylogeny coloured by one with no strip of
/// it beside the tree paints its branches in them too.
fn strip(
    spec: &TrackSpec,
    sheet: Option<&(read::sheet::Sheet, String)>,
    rows: &[String],
    colors: &[(String, Vec<(String, String)>)],
) -> Result<Option<Traits>, BuildError> {
    let Some((held, path)) = sheet else {
        return Ok(None);
    };
    let track = spec.kind.flag();

    let wanted: Vec<String> = match &spec.columns {
        Some(named) => {
            for column in named {
                if !held.columns.contains(column) {
                    return Err(BuildError::Unnamed {
                        what: "column",
                        track,
                        path: path.clone(),
                        wanted: column.clone(),
                        held: held.columns.clone(),
                    });
                }
            }
            named.clone()
        }
        None => held.columns.clone(),
    };

    if held.covers(rows.iter().map(String::as_str)) == 0 {
        return Err(BuildError::Unjoined {
            track,
            path: path.clone(),
            what: "name",
            against: "the rows this track drew",
            examples: held.names().take(3).map(str::to_string).collect(),
        });
    }

    // From the sheet rather than from its rows, so every column deals its
    // levels the palette in the order the file lists them. A phylogeny is
    // handed these same columns, and that shared order is what makes a
    // lineage one colour beside the tree and beside the matrix under it.
    let traits = colors.iter().fold(
        Traits::from_sheet(held).strips(wanted),
        |traits, (key, chosen)| traits.colors(key, chosen.iter().cloned()),
    );
    Ok(Some(traits))
}

/// Refuses a `--colors` that would paint nothing: a column no `--traits`
/// sheet of the figure has, a column of numbers, which is drawn as a ramp,
/// a value no row holds in the column, or a column no track draws.
///
/// Over every sheet of the figure at once, before any track is built,
/// because the colours are the figure's and each sheet takes what it holds
/// of them: a tree's sheet may name a country the matrix's does not, and
/// the country is painted where it is named. One sheet at a time, the
/// matrix would refuse a colour the tree was asked for.
fn colored_as_asked(invocation: &Invocation, files: &mut dyn Files) -> Result<(), BuildError> {
    if invocation.colors.is_empty() {
        return Ok(());
    }
    // Each sheet once, and each track with the sheet it was given.
    let mut sheets: Vec<(read::sheet::Sheet, String)> = Vec::new();
    let mut given: Vec<(&TrackSpec, usize)> = Vec::new();
    for spec in &invocation.tracks {
        if let Some(held) = sheet(spec, files)? {
            let at = match sheets.iter().position(|(_, path)| *path == held.1) {
                Some(at) => at,
                None => {
                    sheets.push(held);
                    sheets.len() - 1
                }
            };
            given.push((spec, at));
        }
    }
    // Capped, since a column can hold as many values as there are rows.
    let shown = |names: Vec<String>| -> Vec<String> {
        let mut shown: Vec<String> = names.iter().take(24).cloned().collect();
        if names.len() > shown.len() {
            shown.push(format!("and {} more", names.len() - shown.len()));
        }
        shown
    };
    for (column, chosen) in &invocation.colors {
        let having: Vec<&read::sheet::Sheet> = sheets
            .iter()
            .filter(|(held, _)| held.columns.contains(column))
            .map(|(held, _)| held)
            .collect();
        if having.is_empty() {
            let mut held: Vec<String> = Vec::new();
            for (sheet, _) in &sheets {
                for name in &sheet.columns {
                    if !held.contains(name) {
                        held.push(name.clone());
                    }
                }
            }
            return Err(BuildError::NotColored {
                column: column.clone(),
                value: None,
                sheets: sheets.iter().map(|(_, path)| path.clone()).collect(),
                held: shown(held),
            });
        }
        // As `Traits::strips` decides it: a column whose every stated value
        // is a number is a ramp.
        let numeric = |sheet: &read::sheet::Sheet| {
            let mut stated = sheet
                .order
                .iter()
                .filter_map(|name| sheet.rows.get(name))
                .filter_map(|row| row.get(column))
                .peekable();
            stated.peek().is_some() && stated.all(|value| value.as_number().is_some())
        };
        if having.iter().all(|sheet| numeric(sheet)) {
            return Err(BuildError::ColorsOfNumbers {
                column: column.clone(),
            });
        }
        // A set, because a column can hold a value per row, and a search of
        // a list per value took eight seconds at 120,000 rows.
        let levels: std::collections::BTreeSet<String> = having
            .iter()
            .flat_map(|sheet| sheet.levels(column))
            .collect();
        if let Some((value, _)) = chosen.iter().find(|(value, _)| !levels.contains(value)) {
            // In the order a key lists them, which is the order a reader
            // looks a misspelt value up in.
            let mut levels: Vec<String> = levels.into_iter().collect();
            levels.sort_by(|a, b| crate::track::traits::natural(a, b));
            return Err(BuildError::NotColored {
                column: column.clone(),
                value: Some(value.clone()),
                sheets: sheets
                    .iter()
                    .filter(|(held, _)| held.columns.contains(column))
                    .map(|(_, path)| path.clone())
                    .collect(),
                held: shown(levels),
            });
        }
        // A column the sheet has and no track draws, as a strip or as the
        // colour of a tree's branches, takes the colours and shows none.
        let drawn = given.iter().any(|(spec, at)| {
            sheets[*at].0.columns.contains(column)
                && (spec
                    .columns
                    .as_ref()
                    .map_or(true, |named| named.contains(column))
                    || spec.color_by.as_ref() == Some(column))
        });
        if !drawn {
            return Err(BuildError::ColorsUndrawn {
                column: column.clone(),
            });
        }
    }
    Ok(())
}

/// Where the text a track is handed came from, which says how much of its
/// file it is and what can be told from it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Origin {
    /// The file as it is, its own lines, all of them.
    Text,
    /// Some of the file's own lines, as they are: its header and the rows
    /// over the window that the tabix index beside it holds. They say what
    /// the file is as its whole text does, but for the counts of its rows,
    /// which [`Slurped::held`] has from the index, and the line a row is on,
    /// which a refusal is read whole again to name.
    Window,
    /// Lines written from a BAM over the window: its reads as SAM, or the
    /// depth they add up to as bedGraph. They say what they were written as
    /// and nothing of the file's own format, and they hold the window alone.
    Bam,
    /// Lines written from a binary format read as it is, over the window: a
    /// bigWig's values as bedGraph and a bigBed's rows as BED, or nothing
    /// for a 2bit, whose bases are [`Slurped::reference`]. Like a BAM's, they
    /// are what they were written as.
    Native(Binary),
}

/// A track's own file, read.
struct Slurped {
    /// What the track's reader takes.
    text: String,
    /// What the file is called in a message.
    path: String,
    /// Where `text` came from.
    origin: Origin,
    /// The lines that tell what the file is, for [`refine`], where `text`
    /// holds only some of the file's own lines and may hold none of them:
    /// its header and its first row, or a BCF's whole header, which names the
    /// samples its sites leave out. `None` where `text` is the file, which
    /// says it itself.
    probe: Option<String>,
    /// How many rows the file holds, where something other than `text`
    /// counted them, as an index counts a file it hands over a window of.
    /// `None` where `text` is all there is to count.
    held: Option<usize>,
    /// The bases over the window of a 2bit, read straight into the reference
    /// a FASTA is read into, where `text` holds nothing.
    reference: Option<Reference>,
    /// The most of a bigWig's values under the window, which the bins of a
    /// zoom level painted with their mean or their least do not reach, for
    /// the track's scale to reach instead.
    most: Option<f64>,
    /// What a bigWig or a bigBed said of a window on a sequence its index
    /// does not name, where `text` holds no rows: a track that finds none
    /// says this, rather than that the region holds none.
    absent: Option<String>,
}

/// A track's own file, read as much of it as the track needs.
///
/// The files are asked in this order, and the first to answer gives the text:
///
/// 1. A BAM, for a track that draws reads or the depth they add up to: its
///    reads over the window as SAM, or their depth as bedGraph, through
///    [`Files::reads`] and [`Files::depth`].
/// 2. A binary format read natively, turned over the window into text a
///    reader already takes.
/// 3. The rows over the window of a bgzipped text file, through the index
///    beside it.
/// 4. The whole text, through [`Files::text`].
///
/// Each comes before the next because the next would misread its file or
/// read more of it. A BAM is not text, and is read through the `.csi` or the
/// `.bai` beside it. A BCF has a `.csi` beside it as a bgzipped VCF does, and
/// only its own reader can read the blocks it points into. A window through an index reads
/// a few blocks of a file the whole text reads all of.
///
/// The second is a bigWig, a bigBed or a 2bit, told by its first bytes and
/// read through [`Files::seekable`] by [`natively`], a BCF, read by
/// [`calls`] through the `.csi` beside it, or a `.hic`, read by
/// [`contact_map`] through the index it holds. The third is a file
/// compressed with bgzip with a `.csi` or a `.tbi` beside it, read by
/// [`windowed`] for a track [`through_index`] says draws only the rows over
/// the window, and only where `whole` is false: a track whose window was
/// refused is read again whole, so the refusal names its line in the file.
fn slurp(
    spec: &TrackSpec,
    region: &Region,
    width: f64,
    files: &mut dyn Files,
    whole: bool,
) -> Result<Slurped, BuildError> {
    let Some(source) = spec.source.as_ref() else {
        return Ok(Slurped {
            text: String::new(),
            path: String::new(),
            origin: Origin::Text,
            probe: None,
            held: None,
            reference: None,
            most: None,
            absent: None,
        });
    };
    let bam = match spec.kind {
        Kind::Coverage => files.depth(source, region),
        Kind::Pileup | Kind::SplitReads => files.reads(source, region),
        _ => Ok(None),
    }
    .map_err(|cause| BuildError::Open {
        track: spec.kind.flag(),
        path: called(source),
        cause,
    })?;
    if let Some(text) = bam {
        return Ok(Slurped {
            text,
            path: called(source),
            origin: Origin::Bam,
            probe: None,
            held: None,
            reference: None,
            most: None,
            absent: None,
        });
    }
    let opened = native(files, source).map_err(|cause| BuildError::Open {
        track: spec.kind.flag(),
        path: called(source),
        cause,
    })?;
    if let Some((binary, file)) = opened {
        if binary == Binary::Bcf {
            return calls(spec, region, files, source, file);
        }
        if binary == Binary::Hic {
            return contact_map(spec, region, files, file, called(source));
        }
        // As many bases to a pixel as the whole figure is wide, which is a
        // few more than the plot inside its gutter holds: a zoom level is
        // never coarser than the drawing for it.
        let per_pixel = region.len() as f64 / width.max(1.0);
        return natively(spec, region, per_pixel, binary, file, called(source));
    }
    let mut untrusted = None;
    if !whole && through_index(spec.kind) {
        match windowed(spec, region, files, source) {
            Ok(window) => return Ok(window),
            Err(why) => untrusted = why,
        }
    }
    let (text, path) = fetch(spec.kind.flag(), source, files)?;
    if let Some(untrusted) = untrusted {
        untrusted.say(spec, region, &text, files);
    }
    Ok(Slurped {
        text,
        path,
        origin: Origin::Text,
        probe: None,
        held: None,
        reference: None,
        most: None,
        absent: None,
    })
}

/// A source that is a binary format a track reads as it is, a bigWig, a
/// bigBed, a 2bit, a BCF or a `.hic`, by its first bytes, opened at its
/// start. `None`
/// for any other source, and for one the files cannot give as bytes, which
/// [`Files::text`] then reads or says why not.
fn native<F: Files + ?Sized>(
    files: &mut F,
    source: &Source,
) -> io::Result<Option<(Binary, Box<dyn Seekable>)>> {
    let Some(mut file) = files.seekable(source)? else {
        return Ok(None);
    };
    Ok(match sniffed(&mut file)? {
        Some(
            binary @ (Binary::BigWig | Binary::BigBed | Binary::TwoBit | Binary::Bcf | Binary::Hic),
        ) => Some((binary, file)),
        _ => None,
    })
}

/// What a file's first bytes say it is, and for a file bgzip wrote, what its
/// first block says is inside it: a BAM and a BCF are BGZF on the outside, as
/// a bgzipped VCF is, and their names need not say so. Left at its start.
///
/// The first block is at most 64 KiB, inflated once for each file a figure
/// asks this of. A BCF written bare, as gzip leaves `-Ou`'s BGZF, is told by
/// its magic alone.
fn sniffed<S: io::Read + io::Seek + ?Sized>(file: &mut S) -> io::Result<Option<Binary>> {
    let mut first = Vec::with_capacity(18);
    io::Read::take(&mut *file, 18).read_to_end(&mut first)?;
    file.seek(io::SeekFrom::Start(0))?;
    let mut binary = Binary::of(&first, None);
    if read::bgzf::is_bgzf(&first) {
        let mut inside = [0u8; 4];
        let got = read::bgzf::Bgzf::new(&mut *file)
            .fill(&mut inside)
            .unwrap_or(0);
        file.seek(io::SeekFrom::Start(0))?;
        match &inside[..got] {
            b"BAM\x01" => binary = Some(Binary::Bam),
            b"BCF\x02" => binary = Some(Binary::Bcf),
            _ => {}
        }
    }
    Ok(binary)
}

/// A bigWig, a bigBed or a 2bit, read over the window as the track that
/// draws what it holds takes it, or refused, naming the track that does.
///
/// A bigWig drawn as a signal is read from the zoom level `per_pixel` bases
/// to a pixel holds two bins of, painted as `--aggregate` takes a pixel;
/// read for windows or for scores under bases, it is read as written. A
/// bigBed is its rows as BED, and a 2bit its bases.
fn natively(
    spec: &TrackSpec,
    region: &Region,
    per_pixel: f64,
    binary: Binary,
    mut file: Box<dyn Seekable>,
    path: String,
) -> Result<Slurped, BuildError> {
    let track = spec.kind.flag();
    // `--format` says what the columns of a text file are, and these say it
    // themselves.
    if spec.format.is_some() {
        return Err(BuildError::OtherTrack {
            track,
            path,
            binary,
            format: true,
        });
    }
    let mut absent = None;
    let (text, reference, most) = match (binary, spec.kind) {
        (Binary::BigWig, Kind::Coverage | Kind::Windows | Kind::Dynseq)
        | (Binary::BigBed, Kind::Features) => {
            let written = match binary {
                Binary::BigWig => {
                    // Windows are drawn whole, and a score under each base,
                    // so those two read the values as written.
                    let (per_pixel, aggregate) = match spec.kind {
                        Kind::Coverage => (per_pixel, spec.aggregate.unwrap_or(Aggregate::Max)),
                        _ => (1.0, Aggregate::Max),
                    };
                    read::bigwig::window(&mut file, region, per_pixel, aggregate).map(|signal| {
                        (
                            read::bigwig::bedgraph(region.seq(), &signal.spans),
                            signal.most,
                        )
                    })
                }
                _ => read::bigbed::bed(&mut file, Some(region)).map(|text| (text, None)),
            };
            match written {
                Ok((text, most)) => (text, None, most),
                Err(error) => {
                    let (cause, unnamed) = refused_window(error, binary, file.as_mut(), region);
                    if !unnamed {
                        return Err(BuildError::Open { track, path, cause });
                    }
                    // A sequence the index does not name holds no rows in
                    // the file, and the track says it found none.
                    absent = Some(cause.to_string());
                    (String::new(), None, None)
                }
            }
        }
        (Binary::TwoBit, Kind::Sequence | Kind::Orfs) => (
            String::new(),
            Some(two_bit(track, &path, file.as_mut(), region)?),
            None,
        ),
        _ => {
            return Err(BuildError::OtherTrack {
                track,
                path,
                binary,
                format: false,
            })
        }
    };
    Ok(Slurped {
        text,
        path,
        origin: Origin::Native(binary),
        probe: None,
        held: None,
        reference,
        most,
        absent,
    })
}

/// A Juicer `.hic`, read over the window as the BEDPE `hictk dump --join`
/// writes for it, which the reader of pairs takes as it takes that, or
/// refused, naming the track that draws it.
///
/// At the resolution `--resolution` names, refused where the file does not
/// hold it, or at the finest that cuts the window into no more than
/// [`read::hic::BINS`] bins, the coarsest where none does, which a note says
/// where the file holds finer: a contact map drawn at its finest over a
/// chromosome is millions of cells, and a figure too large to open.
fn contact_map(
    spec: &TrackSpec,
    region: &Region,
    files: &mut dyn Files,
    mut file: Box<dyn Seekable>,
    path: String,
) -> Result<Slurped, BuildError> {
    let track = spec.kind.flag();
    if spec.format.is_some() || spec.kind != Kind::Pairs {
        return Err(BuildError::OtherTrack {
            track,
            path,
            binary: Binary::Hic,
            format: spec.format.is_some(),
        });
    }
    let open = |path: &str, cause: read::ReadError| BuildError::Open {
        track,
        path: path.to_string(),
        cause: unreadable(cause),
    };
    let header = read::hic::header_of(&mut file).map_err(|error| open(&path, error))?;
    let resolution = match spec.resolution {
        Some(asked) if header.resolutions.contains(&asked) => asked,
        Some(asked) => {
            return Err(BuildError::Unresolved {
                track,
                path,
                asked,
                held: Some(header.resolutions),
            })
        }
        None => read::hic::resolution_for(&header.resolutions, region, read::hic::BINS)
            .ok_or_else(|| {
                open(
                    &path,
                    read::ReadError::whole(
                        "the .hic holds no resolution in bases, only in restriction fragments, \
                         which are not drawn",
                    ),
                )
            })?,
    };
    let mut absent = None;
    let cells = match header.contacts(&mut file, region, resolution) {
        Ok(cells) => cells,
        Err(error) => {
            let (cause, unnamed) = refused_window(error, Binary::Hic, file.as_mut(), region);
            if !unnamed {
                return Err(BuildError::Open { track, path, cause });
            }
            // A sequence the file does not name holds no contacts in it, and
            // the track says it found none.
            absent = Some(cause.to_string());
            Vec::new()
        }
    };
    // Said once the window is read, so a sequence the file does not have is
    // what a figure refused for it says, and nothing about its bins.
    if spec.resolution.is_none() && absent.is_none() {
        let finer = header
            .resolutions
            .iter()
            .copied()
            .filter(|size| *size < resolution)
            .max();
        let size = u64::from(resolution);
        let bins = region.end().saturating_sub(1) / size - region.start() / size + 1;
        let grouped = crate::track::axis::group_thousands;
        if bins > read::hic::BINS {
            let which = if header.resolutions.len() > 1 {
                "its coarsest"
            } else {
                "the one it holds"
            };
            files.note(&format!(
                "{path} is drawn at {}-base bins, {which}, which cut the window into {}",
                grouped(size),
                grouped(bins)
            ));
        } else if let Some(finer) = finer {
            files.note(&format!(
                "{path} is drawn at {}-base bins, the finest of its {} resolutions that keeps \
                 the window to {} bins; --resolution {finer} draws finer",
                grouped(size),
                header.resolutions.len(),
                read::hic::BINS
            ));
        }
    }
    Ok(Slurped {
        text: read::hic::bedpe(region.seq(), &cells),
        path,
        origin: Origin::Native(Binary::Hic),
        probe: None,
        held: None,
        reference: None,
        most: None,
        absent,
    })
}

/// A BCF, read as the VCF text `bcftools view` prints for it, which the
/// track's reader takes as it takes a VCF's, or refused, naming the tracks
/// that draw it.
///
/// `--variants` and `--structural` read the sites alone, and leave the
/// samples' columns undecoded, which are most of a cohort's file and nothing
/// either track draws; `--genotypes` reads each sample's `GT` alone. The rows
/// over the window are read through the `.csi` beside the file where there is
/// one it trusts, as a bgzipped VCF's are through its index, and from every
/// record of the file where there is not. Structural calls are read whole, as
/// they are from a VCF, since an arc is drawn from the lower of its two
/// breakends.
///
/// The header goes with the text as [`Slurped::probe`], whole: the sites
/// alone name no sample, and a cohort's calls named on their own say how many
/// samples `--genotypes` would draw.
fn calls(
    spec: &TrackSpec,
    region: &Region,
    files: &mut dyn Files,
    source: &Source,
    mut file: Box<dyn Seekable>,
) -> Result<Slurped, BuildError> {
    use read::bcf::Fields;
    let track = spec.kind.flag();
    let path = called(source);
    let other = |format: bool| BuildError::OtherTrack {
        track,
        path: path.clone(),
        binary: Binary::Bcf,
        format,
    };
    if spec.format.is_some() {
        return Err(other(true));
    }
    let fields = match spec.kind {
        Kind::Variants | Kind::Structural => Fields::Sites,
        Kind::Genotypes => Fields::Genotypes,
        _ => return Err(other(false)),
    };
    let open = |error: read::ReadError| BuildError::Open {
        track,
        path: path.clone(),
        cause: unreadable(error),
    };
    let header = read::bcf::header_of(&mut file).map_err(open)?;
    let text = if spec.kind == Kind::Structural {
        read::bcf::whole(&mut file, fields)
    } else {
        match trusted_csi(files, source, &path) {
            Some(index) => match read::bcf::window(&mut file, Some(&index), region, fields) {
                Ok(text) => Ok(text),
                // A file that reads whole, read past an index that does not
                // fit it, says so; a damaged file is refused either way.
                Err(error) => {
                    let text = read::bcf::window(&mut file, None, region, fields);
                    if text.is_ok() {
                        let name = csi_name(files, source);
                        files.note(&format!(
                            "{name} does not describe {path}: {error}, so the file was read \
                             whole; bcftools index -f {path} writes it again"
                        ));
                    }
                    text
                }
            },
            None => read::bcf::window(&mut file, None, region, fields),
        }
    }
    .map_err(open)?;
    Ok(Slurped {
        text,
        path,
        origin: Origin::Native(Binary::Bcf),
        probe: Some(header.text),
        held: None,
        reference: None,
        most: None,
        absent: None,
    })
}

/// What the `.csi` beside a source is called, as a note names it.
fn csi_name(files: &mut dyn Files, source: &Source) -> String {
    files.beside(source, ".csi").ok().flatten().map_or_else(
        || format!("the index beside {}", called(source)),
        |found| found.name,
    )
}

/// The `.csi` beside a BCF, where there is one and it is trusted, looked for as
/// htslib looks: `calls.bcf.csi`, then `calls.csi`. A BCF has no other index.
///
/// Not trusted, each said in a note, and the file read whole, which draws the
/// same figure: an index older than its file, in whole seconds, as a
/// bgzipped VCF's is not trusted; one that does not read as an index; and one
/// that is not a CSI written for a BCF, as a tabix index renamed would be.
fn trusted_csi(files: &mut dyn Files, source: &Source, path: &str) -> Option<read::index::Index> {
    match csi_beside(files, source, path) {
        Ok(index) => index,
        Err(note) => {
            files.note(&note);
            None
        }
    }
}

/// The `.csi` beside a BCF as [`trusted_csi`] takes it, without a word: the
/// index where there is one and it is trusted, `None` where there is none,
/// and the note that says why where it is not trusted.
fn csi_beside(
    files: &mut dyn Files,
    source: &Source,
    path: &str,
) -> Result<Option<read::index::Index>, String> {
    let again = format!("bcftools index -f {path} writes it again");
    let beside = match files.beside(source, ".csi") {
        Ok(Some(beside)) => beside,
        Ok(None) => return Ok(None),
        Err(error) => {
            return Err(format!(
                "the index beside {path} could not be read ({error}), so the file was read whole"
            ))
        }
    };
    let name = beside.name;
    if beside.older {
        return Err(format!(
            "{name} is older than {path}, so it was not trusted and the file was read whole; \
             {again}"
        ));
    }
    match read::index::parse(&beside.bytes) {
        Ok(index) if index.kind() == read::index::Kind::Csi && index.columns().is_none() => {
            Ok(Some(index))
        }
        Ok(_) => Err(format!(
            "{name} is not the CSI bcftools index writes for a BCF, so {path} was read whole; \
             {again}"
        )),
        Err(error) => Err(format!(
            "{name} cannot be read as an index ({error}), so {path} was read whole; {again}"
        )),
    }
}

/// A source that is a BCF, with its header, opened at its start. `None` for
/// any other source, and for a BCF whose header will not read, which its
/// track reports.
fn bcf(files: &mut dyn Files, source: &Source) -> Option<(read::bcf::Header, Box<dyn Seekable>)> {
    match native(files, source) {
        Ok(Some((Binary::Bcf, mut file))) => {
            let header = read::bcf::header_of(&mut file).ok()?;
            file.seek(io::SeekFrom::Start(0)).ok()?;
            Some((header, file))
        }
        _ => None,
    }
}

/// How many records a BCF holds on each of its sequences, for the message of
/// a window that held none: from the counts the `.csi` beside it keeps, where
/// it is the index its track reads through and keeps them, and otherwise from
/// the file's records, each read as far as the number of its sequence. `None`
/// for a source that is not a BCF.
///
/// The index is the one [`calls`] trusts and [`read::bcf::window`] reads
/// through: one either passes over is passed over here as well, without a
/// word, since a track that reads through an index has said why already.
/// Counted from, the index of another file put that file's records under
/// this file's names.
fn bcf_rows(files: &mut dyn Files, source: &Source) -> Option<Vec<(String, usize)>> {
    let (_, mut file) = bcf(files, source)?;
    let path = called(source);
    if let Ok(Some(index)) = csi_beside(files, source, &path) {
        if let Ok(counted) = read::bcf::counted(&mut file, Some(&index)) {
            return Some(counted);
        }
    }
    read::bcf::counted(file, None).ok()
}

/// The whole text of a source: [`Files::text`], or for a BCF, every record's
/// site as `bcftools view -G` prints it, which is what a VCF of the same calls
/// would be read for when a figure is placed and its places named.
fn whole_text(files: &mut dyn Files, source: &Source) -> io::Result<String> {
    if let Some((_, file)) = bcf(files, source) {
        return read::bcf::whole(file, read::bcf::Fields::Sites).map_err(unreadable);
    }
    files.text(source)
}

/// The bases of a 2bit over the window, as the reference a FASTA's record is
/// read into, or the refusal of a window that holds none of them, as a
/// FASTA's is refused.
fn two_bit(
    track: &'static str,
    path: &str,
    file: &mut dyn Seekable,
    region: &Region,
) -> Result<Reference, BuildError> {
    let read = read::twobit::bases(&mut *file, region)
        .map_err(|error| refused_window(error, Binary::TwoBit, file, region).0)
        .map_err(|cause| BuildError::Open {
            track,
            path: path.to_string(),
            cause,
        })?;
    if read.bases.is_empty() {
        return Err(BuildError::Beyond {
            track,
            path: path.to_string(),
            record: read.name,
            first: 1,
            last: read.length,
            region: region.to_string(),
        });
    }
    Ok(Reference {
        track,
        path: path.to_string(),
        sequence: read.name.clone(),
        name: read.name,
        offset: read.start,
        bases: read.bases,
    })
}

/// A window a binary file would not give, as the error a source gives, with
/// the `--rename` that would draw it where the file plainly calls the
/// figure's sequence otherwise: `1` for `chr1`, or the one sequence it has.
/// Beside it, whether the file's index names sequences and none of them is
/// the window's, which is why the window was not given.
fn refused_window(
    error: read::ReadError,
    binary: Binary,
    file: &mut dyn Seekable,
    region: &Region,
) -> (io::Error, bool) {
    let names = file.seek(io::SeekFrom::Start(0)).ok().and_then(|_| {
        match binary {
            Binary::BigWig => read::bigwig::sequences(&mut *file),
            Binary::BigBed => read::bigbed::sequences(&mut *file),
            Binary::Hic => read::hic::sequences(&mut *file),
            _ => read::twobit::sequences(&mut *file),
        }
        .ok()
    });
    let held: Vec<(String, usize)> = names
        .unwrap_or_default()
        .into_iter()
        .map(|(name, _)| (name, 0))
        .collect();
    let named = held.iter().any(|(name, _)| name == region.seq());
    let unnamed = !named && !held.is_empty();
    let hint = unnamed
        .then(|| rename_for(&held, region.seq(), true))
        .flatten()
        .map(|(from, to)| format!("; if {from} is {to}, add --rename {from}={to}"))
        .unwrap_or_default();
    (
        io::Error::new(io::ErrorKind::InvalidData, format!("{error}{hint}")),
        unnamed,
    )
}

/// Whether a track of `kind` can be drawn from the rows a tabix index finds
/// over its window, drawing what it draws from the whole file.
///
/// This is a second copy of what the readers do with a row, and has to be
/// kept with them: a kind is here when every row it draws lies at one place,
/// is drawn only where it lies, and is read without the rows around it, so
/// the rows over the window are every row it draws. A reader that comes to
/// draw a row away from where it lies, or to read the rows of a whole file
/// for something, takes its kind off this list. [`windowable`] asks a few
/// more of a file whose kind is here: a GFF3 has to say it is one, since its
/// window is widened to its genes, and [`windowed`] has it write a row for
/// every transcript its exons name; a bedMethyl is read through the index only
/// with `--modification`, since the codes a file holds are offered from all
/// its rows; and a table of windows has to be the wide one, since the long
/// one names its samples on its rows. A scan with no header is read as
/// p-values when every value it holds lies between nought and one, and a
/// window of such values in a file that holds others is refused for it,
/// which [`track`] answers by reading the file whole, as it does any window
/// a reader refuses: the figure is the whole file's either way.
///
/// Read whole, with an index beside them, each for what was measured on 300
/// to 400 random windows of synthetic files against their whole text:
///
/// - Structural calls, `--structural`: an arc is drawn from the lower of its
///   two rows, which lies outside a window it crosses, and 22 windows of 300
///   drew differently. A `<DEL>` with an `END` is indexed over its span.
/// - Genetic maps, `--recombination` and `--with-recombination`: a rate runs
///   from each row to the next, so the rows either side of the window belong
///   to it, and 290 of 300 drew differently.
/// - A GTF, and a GFF3 that does not say it is one: UCSC's GTF has only exon
///   and CDS rows, so a window inside an intron holds no row of its gene even
///   widened, and 328 of 400 drew differently. GENCODE and Ensembl ship the
///   same genes as GFF3, which with its gene rows drew none differently.
/// - A GFF3 whose exons name a transcript it has no row for, for the same
///   reason: the transcript reaches from its first exon to its last, and no
///   row says so.
/// - Pairs, LD tables and PAF, which name a second place the index was not
///   built on; SAM text, which has the BAM's own index; the Bismark extractor
///   file, which is not sorted by position; cytoBand, whose sequence is as
///   long as all its bands; segment tables, read for their samples over the
///   whole file; gene neighbourhoods, joined across genomes; and InterProScan,
///   whose coordinates are residues.
fn through_index(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Variants
            | Kind::Genotypes
            | Kind::Coverage
            | Kind::Windows
            | Kind::Dynseq
            | Kind::Junctions
            | Kind::Methylation
            | Kind::Manhattan
            | Kind::Heatmap
            | Kind::Features
    )
}

/// A bgzipped text file with a tabix index beside it that describes it.
struct Indexed {
    /// The file, as bytes to go back and forth in.
    file: Box<dyn Seekable>,
    /// Its index.
    index: read::index::Index,
    /// Its header and its first row, read and checked against the index.
    head: read::tabix::Head,
    /// What the index is called, as a note names it.
    name: String,
}

impl Indexed {
    /// The sequences the file has rows on, each with how many, in the order
    /// the file first names them, as the index counts them. `None` where the
    /// index leaves a count out, which the format allows, and where the
    /// counts add up to more than a count can hold, which only a damaged
    /// index says. Added up as they came, such counts stopped the program
    /// with an overflow; an index that counts nothing has the file read
    /// whole where its count is wanted.
    fn counted(&self) -> Option<Vec<(String, usize)>> {
        let names = self.index.names();
        let counted: Vec<(String, usize)> = names
            .iter()
            .enumerate()
            .map(|(at, name)| {
                let rows = self.index.summary(at)?.placed;
                Some((name.clone(), usize::try_from(rows).ok()?))
            })
            .collect::<Option<_>>()?;
        counted
            .iter()
            .try_fold(0usize, |all, (_, rows)| all.checked_add(*rows))?;
        Some(counted.into_iter().filter(|(_, rows)| *rows > 0).collect())
    }
}

/// An index beside a file that was not trusted, and why, for the note
/// [`Untrusted::say`] writes once the file is read whole.
struct Untrusted {
    /// The note, which ends with what writes the index again.
    note: String,
    /// The columns the index says the file has, where it reads as an index.
    columns: Option<read::index::Columns>,
}

impl Untrusted {
    /// Why an index was not trusted, as [`indexed`] gives it, with the
    /// columns of `index` where it read as one.
    fn new(note: String, index: Option<&read::index::Index>) -> Option<Self> {
        Some(Untrusted {
            note,
            columns: index.and_then(read::index::Index::columns),
        })
    }

    /// Says why the index was not trusted, where `text`, the whole file read
    /// in its place, shows that a sound one would have been read over
    /// `region`: by the file's header and first row, as [`windowable`] tells
    /// them, and for a GFF3 by every transcript its exons on the region's
    /// sequence name having a row. A GTF, a bedMethyl with no code named and
    /// the other files read whole with an index beside them are read whole
    /// however sound their index, so the note is left out: written again as
    /// it asked, the index changed nothing.
    fn say(self, spec: &TrackSpec, region: &Region, text: &str, files: &mut dyn Files) {
        let header = text
            .split_inclusive('\n')
            .take_while(|line| line.starts_with('#'))
            .map(str::len)
            .sum();
        let zero_based = self.columns.map_or(true, |columns| columns.zero_based);
        let sequence = Region::new(region.seq(), 0, u64::MAX).ok();
        let through = match windowable(spec, &text[..header], text, zero_based) {
            Some((_, models)) => !models || !read::interval::parentless(text, sequence.as_ref()),
            None => false,
        };
        if through {
            files.note(&self.note);
        }
    }
}

/// The tabix index beside a bgzipped text file, read and checked against the
/// file. `Err(None)` for a file that is not gzip or has no index beside it,
/// and `Err` with why for one whose index is not trusted: the file is read
/// whole as though there were none, and [`Untrusted::say`] says why where
/// the index would have been read.
///
/// The index is looked for as htslib looks for it, a `.csi` before a `.tbi`,
/// each after the whole name and then in place of its last extension:
/// `calls.vcf.gz.csi`, `calls.vcf.csi`, `calls.vcf.gz.tbi`, `calls.vcf.tbi`.
/// tabix 1.24 reads a `.csi` beside a damaged `.tbi`, and refuses a damaged
/// `.csi` beside a sound `.tbi`.
///
/// Not trusted, each with its note: an index beside a file compressed with
/// gzip rather than bgzip, which has no blocks for an index to point to; one
/// older than its file, in whole seconds, which may be for an earlier version
/// of it; one that does not read as an index; and one that does not describe
/// the file, by where its rows begin and what its first row is on. htslib
/// warns of an older index and reads through it anyway, and a row the earlier
/// version did not have is then left out of the figure with nothing to show
/// for it. Here the cost is a whole read, which draws the same figure, and a
/// note asking for the index again; a copy that did not keep the times, `cp`
/// without `-p`, `rsync` without `-t` or an archive unpacked without them, is
/// reason enough for it.
///
/// The note waits for the whole read because the header and the first row
/// that say whether the file is read through its index at all are read
/// through the index: written here, an index older than a GTF asked to be
/// written again, and the GTF was read whole as before once it was.
fn indexed(files: &mut dyn Files, source: &Source) -> Result<Indexed, Option<Untrusted>> {
    let Source::Path(path) = source else {
        return Err(None);
    };
    let mut file = files.seekable(source).ok().flatten().ok_or(None)?;
    // The gzip header, and the extra field bgzip writes the size of each
    // block in, which is what an index needs to point into the file.
    let mut first = Vec::with_capacity(18);
    file.by_ref()
        .take(12)
        .read_to_end(&mut first)
        .map_err(|_| None)?;
    if !read::gzip::is_gzip(&first) {
        return Err(None);
    }
    if first.len() == 12 && first[3] & 0x04 != 0 {
        let extra = u64::from(u16::from_le_bytes([first[10], first[11]]));
        file.by_ref()
            .take(extra)
            .read_to_end(&mut first)
            .map_err(|_| None)?;
    }
    // A BAM and a BCF are bgzip on the outside too, and are named for what
    // they are by their own reader or by the whole text's refusal, by their
    // names or by what their first block holds.
    if matches!(
        Binary::of(&first, Some(path)),
        Some(Binary::Bam | Binary::Bcf)
    ) || matches!(sniffed(&mut file), Ok(Some(Binary::Bam | Binary::Bcf)))
    {
        return Err(None);
    }
    let beside = match files.beside(source, ".csi") {
        Ok(None) => files.beside(source, ".tbi"),
        found => found,
    };
    let file_name = called(source);
    let beside = match beside {
        Ok(Some(beside)) => beside,
        Ok(None) => return Err(None),
        Err(error) => {
            return Err(Untrusted::new(
                format!(
                    "the index beside {file_name} could not be read ({error}), so the file \
                     was read whole"
                ),
                None,
            ))
        }
    };
    let name = beside.name.clone();
    let parsed = read::index::parse(&beside.bytes);
    let again = rewritten(&name, parsed.as_ref().ok(), &file_name);
    if !read::bgzf::is_bgzf(&first) {
        return Err(Untrusted::new(
            format!(
                "{file_name} is compressed with gzip rather than bgzip, so {name} beside it \
                 has no blocks to point to, and the file was read whole; gunzip it, bgzip it \
                 and index it again to read it a window at a time"
            ),
            parsed.as_ref().ok(),
        ));
    }
    if beside.older {
        return Err(Untrusted::new(
            format!(
                "{name} is older than {file_name}, so it was not trusted and the file was \
                 read whole; {again}"
            ),
            parsed.as_ref().ok(),
        ));
    }
    let index = match parsed {
        Ok(index) => index,
        Err(error) => {
            return Err(Untrusted::new(
                format!(
                    "{name} cannot be read as an index ({error}), so {file_name} was read \
                     whole; {again}"
                ),
                None,
            ))
        }
    };
    file.seek(io::SeekFrom::Start(0)).map_err(|_| None)?;
    let head = match read::tabix::head(&mut file, &index) {
        Ok(head) => head,
        Err(error) => {
            return Err(Untrusted::new(
                not_described(&name, &file_name, &error, &again),
                Some(&index),
            ))
        }
    };
    Ok(Indexed {
        file,
        index,
        head,
        name,
    })
}

/// The note for an index that does not describe the file beside it.
fn not_described(name: &str, file: &str, error: &read::ReadError, again: &str) -> String {
    format!("{name} does not describe {file}: {error}, so the file was read whole; {again}")
}

/// The command that writes an index again, as the one beside the file was
/// written: `tabix -f -p vcf calls.vcf.gz`, with `-C` for a `.csi`, and the
/// columns written out for a table tabix has no word for.
fn rewritten(name: &str, index: Option<&read::index::Index>, file: &str) -> String {
    use read::index::Preset;
    let csi = if name.ends_with(".csi") { " -C" } else { "" };
    let Some(columns) = index.and_then(read::index::Index::columns) else {
        return format!(
            "tabix -f{csi} {file}, with the columns it was written for, writes it again"
        );
    };
    let options = match (columns.preset, columns.comment, columns.skip) {
        (Preset::Vcf, '#', 0) => "-p vcf".to_string(),
        (Preset::Sam, '@', 0) => "-p sam".to_string(),
        (Preset::Generic, '#', 0)
            if columns.zero_based
                && (columns.sequence, columns.start) == (1, 2)
                && columns.end == Some(3) =>
        {
            "-p bed".to_string()
        }
        (Preset::Generic, '#', 0)
            if !columns.zero_based
                && (columns.sequence, columns.start) == (1, 4)
                && columns.end == Some(5) =>
        {
            "-p gff".to_string()
        }
        _ => {
            let mut options = format!("-s{} -b{}", columns.sequence, columns.start);
            if let Some(end) = columns.end {
                options.push_str(&format!(" -e{end}"));
            }
            if columns.zero_based {
                options.push_str(" -0");
            }
            if columns.skip > 0 {
                options.push_str(&format!(" -S{}", columns.skip));
            }
            if columns.comment != '#' {
                options.push_str(&format!(" -c'{}'", columns.comment));
            }
            options
        }
    };
    format!("tabix -f{csi} {options} {file} writes it again")
}

/// What a track's file is read as through the index beside it, and whether
/// its window is widened to the gene models over it, by the file's header and
/// `probe`, which holds its first row; `None` for a file read whole whatever
/// its index. `zero_based` is whether the index counts from nought, as a
/// BED's does.
///
/// A file named on its own is the kind its first row says. A feature file is
/// a GFF3 that says it is one, whose window is widened, or a BED indexed as
/// one; a GTF, and a GFF3 that does not say it is one, are read whole. So is
/// a bedMethyl with no `--modification` and a long table of windows, as
/// [`through_index`] says why.
fn windowable(
    spec: &TrackSpec,
    header: &str,
    probe: &str,
    zero_based: bool,
) -> Option<(Kind, bool)> {
    let kind = spec
        .guessed
        .then(|| refine(spec.kind, probe))
        .flatten()
        .unwrap_or(spec.kind);
    let gff3 = header.lines().any(|line| {
        line.strip_prefix("##gff-version")
            .is_some_and(|version| version.trim_start().starts_with('3'))
    });
    let models = match kind {
        Kind::Features => match (gff3, spec.format) {
            (true, None | Some(crate::Format::Gff3)) => true,
            (false, None | Some(crate::Format::Bed)) if zero_based => false,
            _ => return None,
        },
        Kind::Methylation if spec.selects.is_none() => return None,
        Kind::Heatmap if read::table::is_long(probe) => return None,
        kind if through_index(kind) => false,
        _ => return None,
    };
    Some((kind, models))
}

/// How far past its first row a GFF3 read through its index is read for
/// exons naming a transcript it never writes, which [`opens_parentless`]
/// looks for.
const OPENING: u64 = 100_000;

/// Whether the rows a GFF3 opens with, over [`OPENING`] bases from its first
/// row, hold a model with no row of its own, as
/// [`read::interval::parentless`] tells one: every exon of a file written
/// without its transcripts' rows does, wherever the window is. `true` where
/// they cannot be read, for the file to be read whole, and `false` for a
/// file of no rows.
fn opens_parentless(found: &mut Indexed, files: &mut dyn Files, path: &str) -> bool {
    let Some(first) = found.head.first_row.clone() else {
        return false;
    };
    let opening = found
        .index
        .columns()
        .and_then(|columns| read::tabix::placed(&first, &columns))
        .and_then(|(sequence, start)| {
            Region::new(sequence, start, start.saturating_add(OPENING)).ok()
        });
    let Some(opening) = opening else {
        return true;
    };
    rows_over(found, &opening, files, path).map_or(true, |text| {
        read::interval::parentless(&text, Some(&opening))
    })
}

/// A track's own file read over the window through the tabix index beside
/// it, where [`through_index`] lists its kind, the index is trusted and the
/// file is one the kind reads the same from its window. `Err(None)` for the
/// whole text to be read, and `Err` with why the index was not trusted, for
/// [`Untrusted::say`] to say once the whole text shows whether it would have
/// been read.
///
/// The text is the header and the rows over the window, which the readers
/// take as they take the whole file. Beside it, the header and the first row
/// for [`refine`] to tell a file named on its own by, and the count of the
/// file's rows the index keeps, for the tracks that say how many rows a file
/// held elsewhere. Those tracks need the count, and are read whole from an
/// index that left it out. A GFF3 is read twice: over the window, and then
/// over as far as the gene models over it reach, so a gene whose intron
/// covers the window comes with its exons.
///
/// A GFF3 whose exons name a transcript it never writes is read whole, as a
/// GTF is: the transcript is drawn from its first exon to its last, and a
/// window between two of them holds neither, so it was left out of the
/// figure, and a window over one drew that exon alone, where the whole file
/// draws the transcript across it. Every row over the window that names its
/// transcript has it among the rows read, since a transcript's row spans its
/// exons, so one that does not is such a file; and the rows the file opens
/// with say the same of a file that writes no transcript's row at all, which
/// is how such files are written. A file that writes its transcripts' rows
/// at its start and leaves them out further on is drawn through its index
/// without the transcripts that leave no row over the window.
fn windowed(
    spec: &TrackSpec,
    region: &Region,
    files: &mut dyn Files,
    source: &Source,
) -> Result<Slurped, Option<Untrusted>> {
    let mut found = indexed(files, source)?;
    let columns = found.index.columns().ok_or(None)?;
    let head = &found.head;
    let probe = format!(
        "{}{}\n",
        head.text,
        head.first_row.as_deref().unwrap_or_default()
    );
    let (kind, models) = windowable(spec, &head.text, &probe, columns.zero_based).ok_or(None)?;
    let held = found
        .counted()
        .map(|counted| counted.iter().map(|(_, rows)| rows).sum());
    if held.is_none() && matches!(kind, Kind::Junctions | Kind::Dynseq | Kind::Methylation) {
        return Err(None);
    }
    let path = called(source);
    let mut text = rows_over(&mut found, region, files, &path).ok_or(None)?;
    if models {
        let (from, to) = read::interval::reach(&text, region, spec.format);
        let mut over = region.clone();
        if from < region.start() || to > region.end() {
            over = Region::new(region.seq(), from, to).map_err(|_| None)?;
            text = rows_over(&mut found, &over, files, &path).ok_or(None)?;
        }
        if read::interval::parentless(&text, Some(&over))
            || opens_parentless(&mut found, files, &path)
        {
            return Err(None);
        }
    } else if kind == Kind::Features
        && read::interval::flavour(&text, spec.format)
            != read::interval::flavour(&probe, spec.format)
    {
        // A BED tells itself from a GFF3 by its first row, and the window's
        // first row has to say what the file's does.
        return Err(None);
    }
    Ok(Slurped {
        text,
        path,
        origin: Origin::Window,
        probe: Some(probe),
        held,
        reference: None,
        most: None,
        absent: None,
    })
}

/// The header of a file read through its index and the rows over `region`,
/// or `None` where a row is not where the index says, which is noted, and
/// the file is read whole.
fn rows_over(
    found: &mut Indexed,
    region: &Region,
    files: &mut dyn Files,
    path: &str,
) -> Option<String> {
    match read::tabix::rows(&mut found.file, &found.index, region) {
        Ok(rows) => Some(found.head.text.clone() + &rows),
        Err(error) => {
            let again = rewritten(&found.name, Some(&found.index), path);
            files.note(&not_described(&found.name, path, &error, &again));
            None
        }
    }
}

/// The tree a file holds, Newick or NEXUS, read once and kept by `parsed`
/// where the caller keeps trees, and the first of several, which is said.
///
/// NEXUS, as BEAST, MrBayes and FigTree write it, was refused with `more than
/// one root`, since it was read as Newick.
fn read_tree(
    flag: &'static str,
    path: &str,
    text: &str,
    parsed: &mut dyn FnMut(&str, &str) -> Option<Tree>,
    files: &mut dyn Files,
) -> Result<Tree, BuildError> {
    let tree = match parsed(path, text.trim()) {
        Some(tree) => tree,
        None => Tree::parse(text.trim()).map_err(|cause| BuildError::Tree {
            flag,
            path: path.to_string(),
            cause,
        })?,
    };
    let held = Tree::count_trees(text);
    if held > 1 {
        files.note(&format!(
            "{flag} {path} holds {} trees, and the first is drawn",
            crate::track::axis::group_thousands(held as u64)
        ));
    }
    Ok(tree)
}

/// What a source is called in a message.
fn called(source: &Source) -> String {
    match source {
        Source::Path(path) => path.display().to_string(),
        Source::Stdin => "standard input".to_string(),
    }
}

/// The tree `--with-tree` names for a track drawn as rows of samples, which
/// orders the rows as its tips and is drawn beside them, or `None` where no
/// tree was named.
///
/// A tree none of whose tips names a row is refused: drawn, it would order
/// nothing and hang beside rows it says nothing about.
fn row_tree(
    spec: &TrackSpec,
    names: &[String],
    files: &mut dyn Files,
    parsed: &mut dyn FnMut(&str, &str) -> Option<Tree>,
) -> Result<Option<Tree>, BuildError> {
    let Some(source) = spec.second.as_ref() else {
        return Ok(None);
    };
    let name = spec.kind.flag();
    let (newick, path) = fetch(name, source, files)?;
    let tree = read_tree("--with-tree", &path, &newick, parsed, files)?;
    let tips = tree.leaf_names();
    if !names.iter().any(|row| tips.contains(row)) {
        return Err(BuildError::Unjoined {
            track: name,
            path,
            what: "tip",
            against: "the rows",
            examples: tips
                .into_iter()
                .filter(|tip| !tip.is_empty())
                .take(3)
                .collect(),
        });
    }
    Ok(Some(tree))
}

/// The place a figure that names no region is drawn over, where its tracks
/// have one of their own, or `None` where none has.
///
/// The columns of an alignment, named after its file. The times of a table of
/// counts or of estimates, from the first to the last, named for the unit its
/// header counts them in: `week`, `year`. The sites of a gene, named `site`,
/// and the samples of a read, `sample`. Every such file is read for its
/// extent and the place is all of them together, so a skyline and the counts
/// it was estimated from share one axis.
fn own_place(
    invocation: &Invocation,
    files: &mut dyn Files,
    decimals: u32,
) -> Result<Option<Region>, BuildError> {
    let mut place: Option<(String, u64, u64)> = None;
    for spec in invocation
        .tracks
        .iter()
        .filter(|spec| spec.kind.own_place())
    {
        let name = spec.kind.flag();
        let Some(source) = spec.source.as_ref() else {
            continue;
        };
        // Read twice, here for its extent and then for its rows, which a
        // pipe can be only because `build_files` keeps what it gave.
        let (text, path) = fetch(name, source, files)?;
        let (unit, start, end) = match spec.kind {
            Kind::Msa | Kind::Logo => {
                let rows = wrap(name, &path, read::seq::alignment(&text))?;
                let width = rows.iter().map(|(_, row)| row.len()).max().unwrap_or(0);
                if width == 0 {
                    return Err(BuildError::Empty {
                        track: name,
                        path,
                        wanted: "sequences",
                    });
                }
                let file = Path::new(&path)
                    .file_name()
                    .map_or(path.clone(), |file| file.to_string_lossy().to_string());
                (file, 0, width as u64)
            }
            Kind::Frequencies => {
                let (rows, unit) = wrap(name, &path, read::series::counts(&text, decimals))?;
                let times = rows.iter().map(|row| row.time);
                let (first, last) = (times.clone().min(), times.max());
                (unit, first.unwrap_or(0), last.map_or(1, |last| last + 1))
            }
            Kind::Phylodynamics => {
                let (points, unit) = wrap(name, &path, read::series::estimates(&text, decimals))?;
                let times = points.iter().map(|point| point.time);
                let (first, last) = (times.clone().min(), times.max());
                (unit, first.unwrap_or(0), last.map_or(1, |last| last + 1))
            }
            Kind::Selection => {
                let sites = wrap(name, &path, read::series::selection(&text))?;
                let at = sites.iter().map(|site| site.pos);
                let (first, last) = (at.clone().min(), at.max());
                (
                    "site".to_string(),
                    first.unwrap_or(0),
                    last.map_or(1, |last| last + 1),
                )
            }
            Kind::Squiggle => {
                let signal = wrap(
                    name,
                    &path,
                    read::series::squiggle(&text, spec.selects.as_deref()),
                )?;
                ("sample".to_string(), 0, signal.samples.len() as u64)
            }
            _ => continue,
        };
        place = Some(match place {
            None => (unit, start, end),
            Some((first, from, to)) => (first, from.min(start), to.max(end)),
        });
    }
    let Some((unit, start, end)) = place else {
        return Ok(None);
    };
    Region::new(unit, start, end)
        .map(Some)
        .map_err(|error| BuildError::Open {
            track: "figure",
            path: String::new(),
            cause: io::Error::new(io::ErrorKind::InvalidData, error.to_string()),
        })
}

/// A ruler that counts something other than bases.
struct Counting {
    /// What it counts, which it says in its margin: `column`, `week`, `site`.
    unit: String,
    /// Whether it counts an alignment's columns. The locus at the top right
    /// names the alignment's file and stays; a week, a site or a sample
    /// names nothing the ruler does not.
    columns: bool,
    /// The places a continuous time is kept to, or nought for whole units.
    decimals: u32,
}

/// What the ruler counts, where every track measured against it has a place
/// of its own, or `None` for a ruler of bases.
///
/// An alignment counts columns. Anything else is counted in the unit its
/// place is named for, which is the unit its header gave or the word a
/// region named it by, as `week:10-30`.
fn counted(invocation: &Invocation, region: &Region, decimals: u32) -> Option<Counting> {
    let first = invocation
        .tracks
        .iter()
        .find(|spec| spec.kind.own_place())?;
    let measured = |kind: Kind| {
        !matches!(
            kind,
            Kind::Tree | Kind::Tanglegram | Kind::Snps | Kind::Axis
        )
    };
    if invocation
        .tracks
        .iter()
        .any(|spec| measured(spec.kind) && !spec.kind.own_place())
    {
        return None;
    }
    let columns = first.kind.over_columns();
    Some(Counting {
        unit: if columns {
            "column".to_string()
        } else {
            region.seq().to_string()
        },
        columns,
        decimals: if all_times(invocation) { decimals } else { 0 },
    })
}

/// Whether every track with a place of its own is a table over time, so the
/// place is a time and may be a continuous one.
fn all_times(invocation: &Invocation) -> bool {
    invocation
        .tracks
        .iter()
        .filter(|spec| spec.kind.own_place())
        .all(|spec| matches!(spec.kind, Kind::Frequencies | Kind::Phylodynamics))
}

/// The places the figure's times are read to: `read::series::DECIMALS` where a
/// table of counts or of estimates named as a file has a fraction of its unit
/// in it, and nought, whole units counted from one, where none has. A table on
/// standard input is not looked at, since not every way of reading one can
/// read it twice; it says so if it turns out to have fractions.
fn time_decimals(invocation: &Invocation, files: &mut dyn Files) -> u32 {
    let fractional = invocation
        .tracks
        .iter()
        .filter(|spec| matches!(spec.kind, Kind::Frequencies | Kind::Phylodynamics))
        .filter_map(|spec| spec.source.as_ref())
        .any(|source| {
            files
                .text(source)
                .is_ok_and(|text| read::series::fractional_times(&text))
        });
    if fractional {
        read::series::DECIMALS
    } else {
        0
    }
}

/// Reads one source, and says what it was called.
///
/// Split out from [`slurp`] because a track drawn from two files opens the
/// second the same way it opened the first, and both belong to the same track
/// as far as an error message is concerned.
fn fetch(
    track: &'static str,
    source: &Source,
    files: &mut dyn Files,
) -> Result<(String, String), BuildError> {
    let name = called(source);
    let text = files.text(source).map_err(|cause| BuildError::Open {
        track,
        path: name.clone(),
        cause,
    })?;
    Ok((text, name))
}

/// Which of the several things a file holds this figure is about.
///
/// One is the answer, whatever the flag says. Several with none named is
/// refused rather than taken the first of, since drawing one of them under a
/// label that names none is the whole of what the flag is for.
/// The one thing in the file, or the one the reader asked for.
///
/// Enumerating what a file holds means reading all of it, and where the reader
/// has already named which one they want, that pass answers a question nobody
/// asked: `chosen` hands a named value straight back without checking it
/// against the list, and the only other use of the list is an emptiness check
/// the second pass makes again on its own. Measured on a bedMethyl of four
/// hundred thousand rows, the enumerating pass was 1.353 seconds against 1.297
/// for the pass that reads the calls, so skipping it halves the work.
fn selected<F>(
    track: &'static str,
    path: &str,
    flag: &'static str,
    asked: &Option<String>,
    enumerate: F,
) -> Result<Option<String>, BuildError>
where
    F: FnOnce() -> Result<std::collections::BTreeMap<String, usize>, BuildError>,
{
    if let Some(name) = asked {
        return Ok(Some(name.clone()));
    }
    let held = enumerate()?;
    if held.is_empty() {
        return Ok(None);
    }
    chosen(track, path, flag, &held, asked).map(Some)
}

fn chosen(
    track: &'static str,
    path: &str,
    flag: &'static str,
    held: &std::collections::BTreeMap<String, usize>,
    asked: &Option<String>,
) -> Result<String, BuildError> {
    match (asked, held.len()) {
        (Some(name), _) => Ok(name.clone()),
        (None, 1) => Ok(held.keys().next().cloned().unwrap_or_default()),
        (None, _) => Err(BuildError::Ambiguous {
            track,
            path: path.to_string(),
            flag,
            choices: held.keys().cloned().collect(),
        }),
    }
}

/// The part of a path a figure has room to print.
///
/// A tanglegram names its two trees after the files they came from, and the
/// name is drawn over the tree, so an absolute path would run across the
/// figure. The last component is what distinguishes two trees in practice and
/// is still exactly what was typed, rather than something made up for the
/// caption.
///
/// The cut is at the separator the system writes, which on Windows is `\` as
/// well as `/`. Cut at `/` alone, `C:\runs\before.nwk` was printed whole over
/// its tree, folders and all, while on Linux a `\` is a letter a file's name
/// may hold and stays in it.
fn shortened(path: &str) -> &str {
    Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path)
}

/// Opens what a command line names, the way a shell would.
///
/// This is what the `karyon` binary hands to [`build`], and it is the only
/// thing in the crate that reads a path. A caller with no filesystem, which is
/// every caller in a browser, passes its own closure instead.
pub fn open_from_disk(source: &Source) -> io::Result<String> {
    let (bytes, path) = match source {
        Source::Path(path) => (fs::read(path)?, Some(path.as_path())),
        Source::Stdin => {
            let mut bytes = Vec::new();
            io::stdin().read_to_end(&mut bytes)?;
            (bytes, None)
        }
    };
    decoded(bytes, path)
}

/// The text a file's bytes hold, whether they were read from a path or held
/// in memory, with `path` the name they go by, for saying what they are when
/// they are not text.
fn decoded(bytes: Vec<u8>, path: Option<&Path>) -> io::Result<String> {
    // Compressed text is text: a `.vcf.gz`, a `.bed.gz` or anything bgzip
    // wrote is taken out of its wrapper here, so every reader takes it as it
    // takes the plain file. BAM and BCF are compressed the same way and are
    // not text inside, so they are named rather than unwrapped, by their name
    // before any work is done and by their first bytes after.
    let bytes = if read::gzip::is_gzip(&bytes) {
        if let Some(binary @ (Binary::Bam | Binary::Bcf)) = Binary::of(&bytes, path) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, binary));
        }
        let inside = read::gzip::decompress(&bytes)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        for (magic, binary) in [
            (&b"BAM\x01"[..], Binary::Bam),
            (&b"BCF\x02"[..], Binary::Bcf),
        ] {
            if inside.starts_with(magic) {
                return Err(io::Error::new(io::ErrorKind::InvalidData, binary));
            }
        }
        if let Some(binary @ (Binary::BigWig | Binary::BigBed | Binary::TwoBit | Binary::Hic)) =
            Binary::of(&inside, None)
        {
            let gzipped = Gzipped {
                binary,
                path: path.map(|path| path.display().to_string()),
            };
            return Err(io::Error::new(io::ErrorKind::InvalidData, gzipped));
        }
        inside
    } else {
        bytes
    };
    let mut text = String::from_utf8(bytes).map_err(|error| {
        // What a genomics file is when it is not text is nearly always one of
        // a few formats, and naming it is what turns "stream did not contain
        // valid UTF-8" into the command that reads it.
        match Binary::of(error.as_bytes(), path) {
            Some(binary) => io::Error::new(io::ErrorKind::InvalidData, binary),
            None => io::Error::new(io::ErrorKind::InvalidData, error),
        }
    })?;
    // A byte order mark is a mark on the file and not its first character.
    // Windows PowerShell's `Out-File -Encoding utf8` and a spreadsheet saved
    // as UTF-8 CSV both write one, and it went on to the reader: the line
    // readers drop it, but a Newick tree was "more than one root" at
    // character 2 and a SLOW5 file had "a raw sample is not a number" on
    // line 1. Dropped here, after the bytes are decoded, it is dropped once
    // for every reader the command line and the playground reach, for text
    // out of a gzip wrapper as for plain text, and for any text a later
    // reader takes through this function.
    if text.starts_with('\u{feff}') {
        text.drain(..'\u{feff}'.len_utf8());
    }
    Ok(text)
}

/// A bigWig, a bigBed, a 2bit or a `.hic` compressed with gzip, which karyon
/// reads as it is once it is not.
///
/// Each is read where its index says a block is, and compressed, every byte
/// of it is somewhere else. Answered as a file that is not text, a gzipped
/// bigWig was told that karyon reads text and given the `bigWigToBedGraph`
/// that writes it as text, which refuses a gzipped file as well. Carried
/// inside the [`io::Error`] [`open_from_disk`] returns, as [`Binary`] is.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Gzipped {
    /// What the gzip holds.
    binary: Binary,
    /// What the file is called, where it has a name.
    path: Option<String>,
}

impl fmt::Display for Gzipped {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let called = self.binary.called();
        // A contact map is called by what it is the first time, and by its
        // ending the second, which reads as one sentence.
        let again = match self.binary {
            Binary::Hic => ".hic",
            _ => called,
        };
        write!(
            f,
            "the file is {called} compressed with gzip, and a {again} is read through the \
             index it holds, which the compression hides; "
        )?;
        // gunzip takes the `.gz` off the name and keeps the file it was
        // given, and with any other name it refuses the file.
        match self.path.as_deref().and_then(|path| {
            path.strip_suffix(".gz")
                .filter(|plain| !plain.is_empty())
                .map(|plain| (path, plain))
        }) {
            Some((path, plain)) => write!(
                f,
                "gunzip -k {path} writes {plain}, which karyon reads as it is"
            ),
            None => write!(
                f,
                "decompress it into a file and name that file, which karyon reads as it is"
            ),
        }
    }
}

impl std::error::Error for Gzipped {}

/// Bytes a reader can go back and forth in, as a file on disk or a buffer can
/// be read and a pipe cannot: what a reader of a binary format, or of a
/// compressed file through the index beside it, takes.
///
/// Everything that reads and seeks is one, so whatever a caller holds will do.
/// [`Files::seekable`] hands one over boxed and owned, which leaves the files
/// free to be asked for the index beside it while it is read.
pub trait Seekable: io::Read + io::Seek {}

impl<T: io::Read + io::Seek> Seekable for T {}

/// A file found beside the one a command line names, as an index sits beside
/// the file it indexes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Beside {
    /// What it is called, as a message would name it.
    pub name: String,
    /// What it holds.
    pub bytes: Vec<u8>,
    /// Whether it was last written before the file it sits beside, where the
    /// files say when they were written, as on a disk; false where they say
    /// nothing, as files held in memory do not. An index older than its file
    /// may have been made for an earlier version of it.
    pub older: bool,
}

/// Where a figure's files come from.
///
/// Text is the one question every source answers: a shell reads a path, a
/// page looks a name up among its buffers, a test hands over a literal, and
/// any closure from a source to its text is a `Files`.
///
/// Two more are what every reader that does not read a file whole stands on:
/// [`Files::seekable`], the file as bytes to go back and forth in, and
/// [`Files::beside`], the file its index would be. A source that cannot
/// answer one says `None`, which leaves the track to [`Files::text`]. [`Disk`]
/// and [`Held`] answer both. A type that wraps another `Files` has to pass
/// both on, since one that does not turns every reader built on them off
/// behind it, and the file is read whole, or refused, with nothing said.
///
/// The questions about a BAM are answered from those two unless a type
/// answers them itself: the depth and the reads over a window, through the
/// `.csi` or the `.bai` beside it, the sequences its header names, and one
/// read by its name.
pub trait Files {
    /// The text a source holds.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn text(&mut self, source: &Source) -> io::Result<String>;

    /// The source as bytes to read out of order, as they are and not taken
    /// out of a gzip wrapper, or `None` for a source this cannot give that
    /// way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being opened.
    fn seekable(&mut self, _source: &Source) -> io::Result<Option<Box<dyn Seekable>>> {
        Ok(None)
    }

    /// The file beside a source whose name is the source's with `ending`
    /// after it, or in place of its last extension, looked for in that order
    /// as samtools and tabix look for an index: `.bai` finds `reads.bam.bai`,
    /// and then `reads.bai`. `None` where there is neither, and for a source
    /// with no name to put an ending on.
    ///
    /// # Errors
    ///
    /// Whatever stopped one that is there being read.
    fn beside(&mut self, _source: &Source, _ending: &str) -> io::Result<Option<Beside>> {
        Ok(None)
    }

    /// The depth of the reads a binary alignment file holds over `region`, as
    /// bedGraph, or `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn depth(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(bam_window(self, source, region)?.map(|(_, records)| {
            read::bam::bedgraph(region.seq(), region, &read::bam::depth(&records, region))
        }))
    }

    /// The reads a binary alignment file holds over `region`, as SAM text, or
    /// `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn reads(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(bam_window(self, source, region)?
            .map(|(header, records)| read::bam::sam(&header, &records)))
    }

    /// The sequences a binary file names, each with its length, or `None`
    /// for a source that is not one this can read: a BAM's header, the index
    /// a bigWig, a bigBed or a 2bit holds, and a `.hic`'s header. A BCF's header is read where
    /// a figure is placed, since it may name a sequence with no length.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn sequences(&mut self, source: &Source) -> io::Result<Option<Vec<(String, u64)>>> {
        if let Some((binary, file)) = native(self, source)? {
            return match binary {
                Binary::BigWig => read::bigwig::sequences(file),
                Binary::BigBed => read::bigbed::sequences(file),
                Binary::TwoBit => read::twobit::sequences(file),
                Binary::Hic => read::hic::sequences(file),
                // A BCF's header may name a sequence with no length, and its
                // records are read for how far they reach, as a VCF's rows
                // are, where the figure is placed.
                _ => return Ok(None),
            }
            .map(Some)
            .map_err(unreadable);
        }
        let Some(file) = opened_bam(self, source)? else {
            return Ok(None);
        };
        read::bam::header_of(file)
            .map(|header| Some(header.references))
            .map_err(unreadable)
    }

    /// The records of one read in a binary alignment file, by its name, as
    /// SAM text, or `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn named_read(&mut self, source: &Source, name: &str) -> io::Result<Option<String>> {
        let Some(file) = opened_bam(self, source)? else {
            return Ok(None);
        };
        read::bam::named(file, name)
            .map(|(header, records)| Some(read::bam::sam(&header, &records)))
            .map_err(unreadable)
    }

    /// Something about a figure drawn anyway that whoever asked for it should
    /// know, such as where the end of a sequence was taken from. The default
    /// keeps it, as a page with nowhere to print it does.
    fn note(&mut self, _message: &str) {}
}

impl<F: FnMut(&Source) -> io::Result<String>> Files for F {
    fn text(&mut self, source: &Source) -> io::Result<String> {
        self(source)
    }
}

/// A file that could not be read, as the error a source gives.
fn unreadable(error: read::ReadError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

/// A source that is a BAM, opened at its start: bgzip on the outside and a
/// BAM's header inside. `None` for any other source, and for one the files
/// cannot give as bytes, which [`Files::text`] then reads or says why not.
fn opened_bam<F: Files + ?Sized>(
    files: &mut F,
    source: &Source,
) -> io::Result<Option<Box<dyn Seekable>>> {
    let Some(mut file) = files.seekable(source)? else {
        return Ok(None);
    };
    let mut first = [0u8; 3];
    if file.read(&mut first)? < first.len() || !read::gzip::is_gzip(&first) {
        return Ok(None);
    }
    file.seek(io::SeekFrom::Start(0))?;
    if read::bam::header_of(&mut file).is_err() {
        return Ok(None);
    }
    file.seek(io::SeekFrom::Start(0))?;
    Ok(Some(file))
}

/// The header and the reads over `region`, for a source that is a BAM, read
/// through the index beside it where there is one and from its start where
/// there is not.
///
/// The index is looked for as htslib looks for it, a `.csi` before a `.bai`,
/// each after the whole name and then in place of its last extension:
/// `reads.bam.csi`, `reads.csi`, `reads.bam.bai`, `reads.bai`. `samtools
/// index -c` writes the first, for a sequence longer than a BAI has room for,
/// and samtools reads it before a `.bai` beside it, so a BAM with both is read
/// through the index samtools reads. Either one, damaged, is refused, as
/// samtools refuses it.
fn bam_window<F: Files + ?Sized>(
    files: &mut F,
    source: &Source,
    region: &Region,
) -> io::Result<Option<(read::bam::Header, Vec<read::bam::Record>)>> {
    let Some(file) = opened_bam(files, source)? else {
        return Ok(None);
    };
    let index = match files.beside(source, ".csi")? {
        Some(found) => match read::index::parse(&found.bytes).map_err(unreadable)? {
            index if index.kind() == read::index::Kind::Csi => Some(index),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("{} is not the CSI samtools index -c writes", found.name),
                ))
            }
        },
        None => match files.beside(source, ".bai")? {
            Some(found) => Some(read::bam::index(&found.bytes).map_err(unreadable)?),
            None => None,
        },
    };
    read::bam::window(file, index.as_ref(), region)
        .map(Some)
        .map_err(unreadable)
}

/// The names a file beside `path` is looked for under, in the order htslib
/// looks: `ending` after the whole name, then in place of its last extension.
fn beside_names(path: &Path, ending: &str) -> [std::path::PathBuf; 2] {
    let mut after = path.as_os_str().to_owned();
    after.push(ending);
    let extension = ending.strip_prefix('.').unwrap_or(ending);
    [
        std::path::PathBuf::from(after),
        path.with_extension(extension),
    ]
}

/// Whether a file last written at `this` was written before one last written
/// at `that`.
///
/// In whole seconds: an archive, or a copy to another kind of disk, keeps no
/// finer a time than that, and a file and its index written by one command in
/// the same second are neither older than the other. A time that is not known
/// is not older.
fn older(this: Option<std::time::SystemTime>, that: Option<std::time::SystemTime>) -> bool {
    let seconds = |time: Option<std::time::SystemTime>| {
        time.and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|since| since.as_secs())
    };
    matches!((seconds(this), seconds(that)), (Some(this), Some(that)) if this < that)
}

/// Bytes held once and handed out as often as they are asked for, without a
/// copy each time: a page's files, and what a pipe gave.
#[derive(Debug, Clone)]
struct Shared(std::sync::Arc<Vec<u8>>);

impl From<Vec<u8>> for Shared {
    fn from(bytes: Vec<u8>) -> Self {
        Shared(std::sync::Arc::new(bytes))
    }
}

impl AsRef<[u8]> for Shared {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

/// The files a figure is drawn from, with what standard input held kept once
/// it is read.
///
/// A figure reads some sources twice: an alignment, a table over time or a
/// signal for its extent and then for its rows, and every table of times to
/// see whether any of them has fractions. A pipe can be read once, so what it
/// gave is kept here, whatever the files underneath keep, and a table piped in
/// is its own place as it is from a file. It was refused, and a table with
/// fractions piped into a figure of whole units was refused too.
struct KeptStdin<'a> {
    files: &'a mut dyn Files,
    stdin: Option<String>,
}

impl Files for KeptStdin<'_> {
    fn text(&mut self, source: &Source) -> io::Result<String> {
        if !matches!(source, Source::Stdin) {
            return self.files.text(source);
        }
        if let Some(text) = &self.stdin {
            return Ok(text.clone());
        }
        let text = self.files.text(source)?;
        self.stdin = Some(text.clone());
        Ok(text)
    }

    fn seekable(&mut self, source: &Source) -> io::Result<Option<Box<dyn Seekable>>> {
        // Standard input is kept here as the text it held, which is all a
        // reader of it is given: asked for its bytes underneath, it would be
        // read past, and the text asked for after found empty.
        if matches!(source, Source::Stdin) {
            return Ok(None);
        }
        self.files.seekable(source)
    }

    fn beside(&mut self, source: &Source, ending: &str) -> io::Result<Option<Beside>> {
        self.files.beside(source, ending)
    }

    fn depth(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        self.files.depth(source, region)
    }

    fn reads(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        self.files.reads(source, region)
    }

    fn sequences(&mut self, source: &Source) -> io::Result<Option<Vec<(String, u64)>>> {
        self.files.sequences(source)
    }

    fn named_read(&mut self, source: &Source, name: &str) -> io::Result<Option<String>> {
        self.files.named_read(source, name)
    }

    fn note(&mut self, message: &str) {
        self.files.note(message);
    }
}

/// Keeps `message` unless it is kept already.
///
/// A figure tells the files the same thing more than once: a figure of
/// several places builds its tracks once a panel, and a file read under the
/// name `--rename` gives its sequence is read again under that name. Said
/// each time, a BAM over two genes printed that it was drawn as its depth
/// twice, word for word.
fn keep_once(notes: &mut Vec<String>, message: &str) {
    if !notes.iter().any(|kept| kept == message) {
        notes.push(message.to_string());
    }
}

/// The files a command line names, read from disk.
///
/// Text is read whole, and compressed text is taken out of its wrapper, by
/// [`open_from_disk`]. A BAM is read a window at a time, through the `.csi`
/// or the `.bai` beside it where there is one, so a figure of one gene reads
/// the blocks that gene is in, and so is a BCF through its `.csi`, and a
/// bgzipped text file through the `.csi` or `.tbi` beside it, for a track that
/// draws only the rows over its window.
///
/// A pipe the shell named, as `<(zcat genes.gff3.gz)`, cannot be read twice,
/// so it is read whole the first time anything asks and its bytes are kept,
/// since placing a figure by a gene's name reads the annotation before the
/// track does. Looking for a BAM's magic in one read its first three bytes
/// and dropped them, and the text read after began at the fourth: a bedGraph
/// of `chr1` rows was said to hold rows on `1`. Standard input is read once and
/// kept by [`build_files`], which every figure is drawn through. A file on
/// disk is read again instead: kept, every file was held twice while its
/// track was built, once here and once in the text handed out, and a depth
/// file of 176 MB took 435 MB to draw where it now takes 259.
#[derive(Debug, Default)]
pub struct Disk {
    kept: std::collections::HashMap<Source, Shared>,
    /// What [`Files::note`] was told, each thing once, for the command line
    /// to print.
    pub notes: Vec<String>,
}

impl Disk {
    /// What a pipe the command line names held, read the first time it is
    /// asked for and kept after. `None` for a file on disk, which is read
    /// again, and for one that is not there, which reading it says.
    fn pipe(&mut self, source: &Source) -> io::Result<Option<Shared>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        if let Some(bytes) = self.kept.get(source) {
            return Ok(Some(bytes.clone()));
        }
        match fs::metadata(path) {
            Ok(meta) if !meta.is_file() => {
                let bytes = Shared::from(fs::read(path)?);
                self.kept.insert(source.clone(), bytes.clone());
                Ok(Some(bytes))
            }
            _ => Ok(None),
        }
    }
}

impl Files for Disk {
    fn note(&mut self, message: &str) {
        keep_once(&mut self.notes, message);
    }

    fn text(&mut self, source: &Source) -> io::Result<String> {
        match (self.pipe(source)?, source) {
            (Some(bytes), Source::Path(path)) => decoded(bytes.as_ref().to_vec(), Some(path)),
            _ => open_from_disk(source),
        }
    }

    fn seekable(&mut self, source: &Source) -> io::Result<Option<Box<dyn Seekable>>> {
        if let Some(bytes) = self.pipe(source)? {
            return Ok(Some(Box::new(io::Cursor::new(bytes))));
        }
        match source {
            Source::Path(path) => Ok(Some(Box::new(io::BufReader::new(fs::File::open(path)?)))),
            Source::Stdin => Ok(None),
        }
    }

    fn beside(&mut self, source: &Source, ending: &str) -> io::Result<Option<Beside>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        let written = |path: &Path| fs::metadata(path).and_then(|meta| meta.modified()).ok();
        for candidate in beside_names(path, ending) {
            if candidate.is_file() {
                return Ok(Some(Beside {
                    name: candidate.display().to_string(),
                    bytes: fs::read(&candidate)?,
                    older: older(written(&candidate), written(path)),
                }));
            }
        }
        Ok(None)
    }
}

/// The files a command line names, held in memory by name.
///
/// What [`Disk`] reads from a path, read from bytes a caller already has: a
/// page that fetched them, a service that was sent them, a test. Compressed
/// text is taken out of its wrapper, a BAM is read a window at a time through
/// the `.csi` or the `.bai` held beside it, a BCF through its `.csi`, and a
/// bgzipped text file through the `.csi` or `.tbi` held beside it, and one read is found by its name, so a figure
/// drawn from these files is the figure a shell draws from the same files on
/// disk. Bytes held say nothing of when they were written, so an index held
/// for an earlier version of its file is told only by not describing it. A
/// file is found by its name as the command line writes it, and the file
/// beside it, as an index, by that name with its ending.
///
/// ```
/// use karyon::cli::{args, stack};
///
/// let mut files = stack::Held::new();
/// files.insert("depth.bg", "chr1\t0\t100\t12\n");
/// let argv = ["chr1:1-100", "depth.bg"].map(String::from);
/// let args::Request::Draw(invocation) = args::parse(&argv).unwrap() else {
///     unreachable!("a command line that draws")
/// };
/// let svg = stack::build_files(&invocation, &mut files, |_, _| None).unwrap();
/// assert!(svg.starts_with("<svg"));
/// ```
#[derive(Debug, Default)]
pub struct Held {
    files: std::collections::BTreeMap<String, Shared>,
    /// What [`Files::note`] was told, each thing once.
    pub notes: Vec<String>,
}

impl Held {
    /// Holds no files yet.
    pub fn new() -> Self {
        Held::default()
    }

    /// Holds `bytes` as the file called `name`, in place of any held before.
    pub fn insert(&mut self, name: impl Into<String>, bytes: impl Into<Vec<u8>>) {
        self.files.insert(name.into(), Shared::from(bytes.into()));
    }

    /// Whether a file called `name` is held.
    pub fn contains(&self, name: &str) -> bool {
        self.files.contains_key(name)
    }

    /// The names of the files held, in order.
    pub fn names(&self) -> impl Iterator<Item = &str> + '_ {
        self.files.keys().map(String::as_str)
    }

    /// The bytes held for a source, and the name they are held under.
    fn held<'s>(&self, source: &'s Source) -> io::Result<(&'s Path, &[u8])> {
        let Source::Path(path) = source else {
            return Err(io::Error::other(
                "nothing is piped into files held in memory",
            ));
        };
        let name = path.display().to_string();
        match self.files.get(&name) {
            Some(bytes) => Ok((path.as_path(), bytes.as_ref())),
            None => {
                let names: Vec<&str> = self.names().collect();
                let held = if names.is_empty() {
                    ", and no files are".to_string()
                } else {
                    format!("; the files held are {}", names.join(", "))
                };
                Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("no file called {name} is held{held}"),
                ))
            }
        }
    }
}

impl Files for Held {
    fn note(&mut self, message: &str) {
        keep_once(&mut self.notes, message);
    }

    fn text(&mut self, source: &Source) -> io::Result<String> {
        let (path, bytes) = self.held(source)?;
        decoded(bytes.to_vec(), Some(path))
    }

    /// The bytes held for a source, or `None` for one not held, which
    /// [`Files::text`] then says, naming the files that are.
    fn seekable(&mut self, source: &Source) -> io::Result<Option<Box<dyn Seekable>>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        Ok(self
            .files
            .get(&path.display().to_string())
            .map(|bytes| Box::new(io::Cursor::new(bytes.clone())) as Box<dyn Seekable>))
    }

    fn beside(&mut self, source: &Source, ending: &str) -> io::Result<Option<Beside>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        Ok(beside_names(path, ending)
            .into_iter()
            .find_map(|candidate| {
                let name = candidate.display().to_string();
                self.files.get(&name).map(|bytes| Beside {
                    name,
                    bytes: bytes.as_ref().to_vec(),
                    older: false,
                })
            }))
    }
}

/// A file that is not text, by the format its first bytes say it is.
///
/// Carried inside the [`io::Error`] [`open_from_disk`] returns, so that
/// [`build`] can turn it into a message naming the tool that reads it, and a
/// caller opening files some other way is not obliged to know about it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Binary {
    /// gzip or bgzip, and not BAM or BCF by its name.
    Gzip,
    /// bzip2.
    Bzip2,
    /// xz.
    Xz,
    /// Zstandard.
    Zstd,
    /// BAM, which is bgzip by its bytes and BAM by its name.
    Bam,
    /// CRAM.
    Cram,
    /// BCF, which is bgzip by its first bytes and BCF by its name or by what
    /// its first block holds, or BCF by its first bytes where it is written
    /// bare.
    Bcf,
    /// bigWig.
    BigWig,
    /// bigBed.
    BigBed,
    /// UCSC's 2bit.
    TwoBit,
    /// HDF5, as cooler writes a contact map at one resolution, a `.cool`.
    Cool,
    /// HDF5, as cooler writes a contact map at several resolutions, a
    /// `.mcool`, which is one by its bytes and the other by its name.
    Mcool,
    /// Juicer's `.hic`.
    Hic,
}

impl Binary {
    /// The format a file's first bytes, and for bgzip its name, say it is.
    pub fn of(bytes: &[u8], path: Option<&Path>) -> Option<Binary> {
        let named = |ending: &str| {
            path.and_then(|path| path.extension())
                .and_then(|extension| extension.to_str())
                .is_some_and(|extension| extension.eq_ignore_ascii_case(ending))
        };
        let magic = |mark: &[u8]| bytes.starts_with(mark);
        // The UCSC formats write a four byte number in the machine's order, so
        // either order is the format.
        let either = |mark: [u8; 4]| {
            let mut turned = mark;
            turned.reverse();
            magic(&mark) || magic(&turned)
        };
        Some(if magic(&[0x1f, 0x8b]) {
            if named("bam") {
                Binary::Bam
            } else if named("bcf") {
                Binary::Bcf
            } else {
                Binary::Gzip
            }
        } else if magic(b"BCF\x02") {
            Binary::Bcf
        } else if magic(b"CRAM") {
            Binary::Cram
        } else if magic(b"BZh") {
            Binary::Bzip2
        } else if magic(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]) {
            Binary::Xz
        } else if magic(&[0x28, 0xb5, 0x2f, 0xfd]) {
            Binary::Zstd
        } else if either([0x26, 0xfc, 0x8f, 0x88]) {
            Binary::BigWig
        } else if either([0xeb, 0xf2, 0x89, 0x87]) {
            Binary::BigBed
        } else if either([0x43, 0x27, 0x41, 0x1a]) {
            Binary::TwoBit
        } else if magic(b"\x89HDF\r\n\x1a\n") {
            if named("mcool") {
                Binary::Mcool
            } else {
                Binary::Cool
            }
        } else if magic(b"HIC\0") {
            Binary::Hic
        } else {
            return None;
        })
    }

    /// What to do with a file no single command turns into text, where there
    /// is more to say than to pipe it through the tool that writes it.
    fn advice(self) -> Option<&'static str> {
        Some(match self {
            Binary::Mcool => {
                "cooler ls lists its resolutions, and cooler dump --join -r REGION \
                 FILE::/resolutions/N writes one of them as the BEDPE --pairs reads"
            }
            _ => return None,
        })
    }

    /// What a format karyon reads as it is holds, and which tracks draw it,
    /// as a sentence says them; `None` for a format read through a pipe.
    fn drawn_by(self) -> Option<(&'static str, &'static str)> {
        Some(match self {
            Binary::BigWig => (
                "a signal along the sequence",
                "--coverage, --windows and --dynseq draw",
            ),
            Binary::BigBed => ("intervals along the sequence", "--features draws"),
            Binary::TwoBit => (
                "the bases of a reference",
                "--sequence and --orfs draw and --with-sequence reads",
            ),
            Binary::Bcf => (
                "calls along the sequence",
                "--variants, --genotypes and --structural draw",
            ),
            Binary::Hic => ("contacts between the bins of a sequence", "--pairs draws"),
            _ => return None,
        })
    }

    /// What the format is called, as a sentence would say it.
    fn called(self) -> &'static str {
        match self {
            Binary::Gzip => "compressed with gzip",
            Binary::Bzip2 => "compressed with bzip2",
            Binary::Xz => "compressed with xz",
            Binary::Zstd => "compressed with zstd",
            Binary::Bam => "BAM",
            Binary::Cram => "CRAM",
            Binary::Bcf => "BCF",
            Binary::BigWig => "bigWig",
            Binary::BigBed => "bigBed",
            Binary::TwoBit => "2bit",
            Binary::Cool => "a contact map in cooler's HDF5",
            Binary::Mcool => "a contact map at several resolutions in cooler's HDF5",
            Binary::Hic => "a contact map in Juicer's .hic",
        }
    }

    /// The command that writes the text a `kind` track reads out of `path`,
    /// cut to the window where the tool can do that, or `None` for a format
    /// no one command turns into text, which [`Binary::advice`] speaks for.
    fn reader(self, kind: Kind, path: &str, region: Option<&Region>) -> Option<String> {
        let window = region.map(|region| region.to_string());
        let near = |region: &Region| {
            format!(
                "-chrom={} -start={} -end={}",
                region.seq(),
                region.start(),
                region.end()
            )
        };
        Some(match self {
            Binary::Gzip => format!("gzip -dc {path}"),
            Binary::Bzip2 => format!("bzip2 -dc {path}"),
            Binary::Xz => format!("xz -dc {path}"),
            Binary::Zstd => format!("zstd -dc {path}"),
            Binary::Bam | Binary::Cram => match (kind, window) {
                (Kind::Coverage, Some(window)) => format!("samtools depth -a -r {window} {path}"),
                (Kind::Coverage, None) => format!("samtools depth -a {path}"),
                (_, Some(window)) => format!("samtools view -h {path} {window}"),
                (_, None) => format!("samtools view -h {path}"),
            },
            Binary::Bcf => format!("bcftools view {path}"),
            Binary::BigWig => match region {
                Some(region) => format!("bigWigToBedGraph {} {path} /dev/stdout", near(region)),
                None => format!("bigWigToBedGraph {path} /dev/stdout"),
            },
            Binary::BigBed => match region {
                Some(region) => format!("bigBedToBed {} {path} /dev/stdout", near(region)),
                None => format!("bigBedToBed {path} /dev/stdout"),
            },
            Binary::TwoBit => match region {
                Some(region) => format!("twoBitToFa -seq={} {path} /dev/stdout", region.seq()),
                None => format!("twoBitToFa {path} /dev/stdout"),
            },
            Binary::Cool => match window {
                Some(window) => format!("cooler dump --join -r {window} {path}"),
                None => format!("cooler dump --join {path}"),
            },
            Binary::Mcool | Binary::Hic => return None,
        })
    }
}

impl fmt::Display for Binary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "the file is {}, and karyon reads text", self.called())
    }
}

impl std::error::Error for Binary {}

/// What a track is called when `--label` does not say: its file's name, less
/// the folder, the compression and the format, as `reads` for
/// `data/reads.bam` and `calls` for `calls.vcf.gz`.
///
/// A track in the gutter with no name is a band a reader has to work out, and
/// the file's name is the name its owner gave it. Standard input and a pipe
/// the shell names, `/dev/fd/63`, name nothing anyone chose, and a tanglegram
/// is two files, so neither gets one.
fn default_label(spec: &TrackSpec) -> Option<String> {
    // A phylogeny is plain to see, and named after its file it put a stray
    // word in the margin of every tree: three simulated users asked what it was.
    if matches!(spec.kind, Kind::Tanglegram | Kind::Axis | Kind::Tree) {
        return None;
    }
    let Some(Source::Path(path)) = &spec.source else {
        return None;
    };
    if path.starts_with("/dev") || path.starts_with("/proc") {
        return None;
    }
    let file = path.file_name()?.to_str()?;
    let file = file
        .strip_suffix(".gz")
        .or_else(|| file.strip_suffix(".bgz"))
        .unwrap_or(file);
    let (stem, extension) = file.rsplit_once('.').unwrap_or((file, ""));
    if stem.is_empty() {
        return None;
    }
    // A BAM drawn as a line is its depth, and a name that says reads would
    // have the line read as a count of them.
    if spec.kind == Kind::Coverage && extension.eq_ignore_ascii_case("bam") {
        return Some(format!("{stem} depth"));
    }
    // A VCF's genotypes are drawn under its calls as often as not, and two
    // bands both called `calls` would leave the reader to tell them apart.
    if spec.kind == Kind::Genotypes {
        return Some(format!("{stem} genotypes"));
    }
    Some(stem.to_string())
}

/// The kind a file named on its own turns out to be once it is read, where
/// that differs from what its name said.
///
/// Two cases, both `.bed` by name: modkit's bedMethyl, eighteen columns with
/// a modification code in the fourth, is methylation; and four columns whose
/// fourth is a number is a bedGraph, a signal.
fn refine(kind: Kind, text: &str) -> Option<Kind> {
    if kind != Kind::Features {
        return None;
    }
    let (_, first) = read::lines(text).next()?;
    let cols = read::columns(first);
    let number = |at: usize| {
        cols.get(at)
            .is_some_and(|field| field.trim().parse::<f64>().is_ok())
    };
    if cols.len() >= 18 && cols[3].trim().len() <= 8 && number(9) && number(10) {
        return Some(Kind::Methylation);
    }
    if cols.len() == 4 && number(1) && number(2) && number(3) {
        return Some(Kind::Coverage);
    }
    None
}

/// A refusal made plainer where the program can see what went wrong.
///
/// The two a first attempt meets most. A file that is not text failed with
/// "stream did not contain valid UTF-8", which says nothing about what the
/// file is or what reads it; it now names the format and the command to write
/// in place of the file's name. And a window that held nothing said only that,
/// though the usual reason is a file that names its sequences `1` where the
/// region says `chr1`; it now says what the file does hold, and where.
fn explained(
    error: BuildError,
    spec: &TrackSpec,
    region: Option<&Region>,
    files: &mut dyn Files,
) -> BuildError {
    match error {
        BuildError::Open { track, path, cause } => {
            let binary = cause
                .get_ref()
                .and_then(|inner| inner.downcast_ref::<Binary>())
                .copied();
            match binary {
                Some(binary) => {
                    // A pipe has no name to write a command in place of. The
                    // file is told to be one by its own source, of the several
                    // a track reads, its own, a second file, a sheet of traits
                    // or a genetic map, and not by the track's: a bigWig named
                    // after --ld, beside a scan piped in, was asked for the
                    // name it had been given.
                    let piped = spec
                        .sources()
                        .any(|source| *source == Source::Stdin && called(source) == path);
                    let instead = (!piped)
                        .then(|| binary.reader(spec.kind, &path, region))
                        .flatten();
                    // The track's own file, and not another file it reads,
                    // which came with a flag of its own.
                    let own = spec.source.as_ref().is_some_and(|own| called(own) == path);
                    // A format karyon reads as it is, named for another file
                    // of the track, which is read as text, such as the
                    // linkage --ld gives a scan, is told which track draws
                    // it where no one command writes it as text: it is no
                    // pipe, and its name is the one thing it was not missing.
                    if instead.is_none() && !piped && !own && binary.drawn_by().is_some() {
                        return BuildError::OtherTrack {
                            track,
                            path,
                            binary,
                            format: false,
                        };
                    }
                    // Named on its own, with no flag at all.
                    let alone = spec.guessed && own;
                    BuildError::NotText {
                        track,
                        path,
                        binary,
                        instead,
                        alone,
                    }
                }
                None => BuildError::Open { track, path, cause },
            }
        }
        BuildError::Empty {
            track,
            path,
            wanted,
        } => {
            // Read again, and only here: a pipe cannot be, and a figure that
            // drew has no need to know. A file with an index is not read
            // again at all where the index counts its rows, which is the
            // usual way this fails, a file that says 1 where the figure says
            // chr1, answered without inflating a whole genome's calls.
            let held = match (region, spec.source.as_ref()) {
                (Some(_), Some(source @ Source::Path(_))) => bcf_rows(files, source)
                    .or_else(|| {
                        through_index(spec.kind)
                            .then(|| indexed(files, source).ok())
                            .flatten()
                            .and_then(|found| found.counted())
                    })
                    .unwrap_or_else(|| {
                        files
                            .text(source)
                            .map(|text| rows_on(spec.kind, &text))
                            .unwrap_or_default()
                    }),
                _ => Vec::new(),
            };
            match region {
                Some(region) if !held.is_empty() => BuildError::Elsewhere {
                    track,
                    path,
                    wanted,
                    held: held.iter().map(|(_, count)| count).sum(),
                    named: listed(&held),
                    region: region.to_string(),
                    rename: rename_for(&held, region.seq(), true).map(Box::new),
                },
                _ => BuildError::Empty {
                    track,
                    path,
                    wanted,
                },
            }
        }
        other => other,
    }
}

/// Where a track's file writes the sequence a row is on.
struct SequenceColumn {
    /// The column holding the sequence's name.
    name: usize,
    /// Columns one of which holds a position, which is what tells a row from
    /// a header: a header has a word there.
    position: &'static [usize],
    /// The fewest columns a row has.
    width: usize,
}

/// The column each format keeps its sequence in, for the tracks that read
/// one.
fn sequence_column(kind: Kind) -> Option<SequenceColumn> {
    let at = |name: usize, position: &'static [usize], width: usize| SequenceColumn {
        name,
        position,
        width,
    };
    Some(match kind {
        Kind::Coverage
        | Kind::Windows
        | Kind::Dynseq
        | Kind::Junctions
        | Kind::Methylation
        | Kind::Variants
        | Kind::Genotypes
        | Kind::Structural
        | Kind::Ideogram
        | Kind::CopyNumber => at(0, &[1], 2),
        Kind::Heatmap => at(0, &[1], 4),
        Kind::Recombination => at(0, &[1], 3),
        Kind::Pairs => at(0, &[1], 3),
        // BED puts the start in column two and GFF3 in column four.
        Kind::Features => at(0, &[1, 3], 3),
        Kind::Clades => at(0, &[3], 4),
        // A table of two columns is a position and a value, and names none.
        Kind::Manhattan => at(0, &[1], 3),
        Kind::Pileup | Kind::SplitReads | Kind::Bisulfite => at(2, &[3], 4),
        Kind::Domains => at(0, &[6], 7),
        _ => return None,
    })
}

/// The sequences a track's file has rows on, each once with how many rows,
/// in the order the file first names them.
///
/// An association table says in its header which column that is, and BOLT-LMM
/// writes it second, behind each variant's own name.
fn rows_on(kind: Kind, text: &str) -> Vec<(String, usize)> {
    if kind == Kind::Manhattan {
        return read::point::association_sequences(text);
    }
    sequence_column(kind).map_or_else(Vec::new, |column| sequences_named(text, column))
}

/// The sequences a file's rows are on, each once with how many rows it has,
/// in the order the file first names them.
fn sequences_named(text: &str, column: SequenceColumn) -> Vec<(String, usize)> {
    let mut held: Vec<(String, usize)> = Vec::new();
    for (_, line) in read::lines(text) {
        let cols = read::columns(line);
        let placed = column.position.iter().any(|at| {
            cols.get(*at)
                .is_some_and(|field| field.trim().parse::<u64>().is_ok())
        });
        if cols.len() < column.width || !placed {
            continue;
        }
        let name = cols[column.name].trim();
        match held.iter_mut().find(|(seen, _)| seen == name) {
            Some((_, count)) => *count += 1,
            None => held.push((name.to_string(), 1)),
        }
    }
    held
}

/// The widest window over which a BAM named on its own says it could be
/// drawn as reads: a thousand bases, where a read is a few dozen pixels long
/// and its mismatches can be told apart.
const READS_WINDOW: u64 = 1_000;

/// The `--rename` that would draw a file's rows under the figure's name
/// `figure`, where one plainly would: the name among `held` that is the
/// figure's own with a `chr` more or less, or, where `alone` allows, the one
/// sequence the file names. None where the file already names the figure's.
fn rename_for(held: &[(String, usize)], figure: &str, alone: bool) -> Option<(String, String)> {
    if held.iter().any(|(name, _)| name == figure) {
        return None;
    }
    let bare = |name: &str| {
        let lower = name.to_ascii_lowercase();
        lower.strip_prefix("chr").unwrap_or(&lower).to_string()
    };
    let spelled: Vec<&String> = held
        .iter()
        .map(|(name, _)| name)
        .filter(|name| bare(name) == bare(figure))
        .collect();
    let from = match (spelled.as_slice(), held) {
        ([one], _) => (*one).clone(),
        ([], [(one, _)]) if alone => one.clone(),
        _ => return None,
    };
    Some((from, figure.to_string()))
}

/// Names for a sentence: `chr1`, `chr1 and chr2`, `1, 2 and 3`, and past
/// five, how many more.
fn listed(held: &[(String, usize)]) -> String {
    listed_as(held, "sequences")
}

/// [`listed`], for names of things other than sequences.
fn listed_as(held: &[(String, usize)], things: &str) -> String {
    const SHOWN: usize = 5;
    let names: Vec<&str> = held
        .iter()
        .take(SHOWN)
        .map(|(name, _)| name.as_str())
        .collect();
    let rest = held.len().saturating_sub(SHOWN);
    match (names.as_slice(), rest) {
        ([one], 0) => (*one).to_string(),
        (many, 0) => format!(
            "{} and {}",
            many[..many.len() - 1].join(", "),
            many[many.len() - 1]
        ),
        (many, rest) => format!("{} and {rest} more {things}", many.join(", ")),
    }
}

/// Builds one track from its flags and its file.
///
/// A file read through its index whose window a reader refuses is read again
/// whole, and the track is built from that. A refusal names the line it is
/// on, and a row's line in the window's text is not its line in the file.
/// And a reader may refuse a window for what the rest of the file answers: a
/// scan with no header whose values over the window all lie between 0 and 1
/// is taken for p-values, which the values elsewhere may say it is not. Read
/// whole, the track draws or is refused as the whole file always had it, for
/// one more read of a figure that was being refused.
fn track(
    spec: &TrackSpec,
    context: &Context<'_>,
    files: &mut dyn Files,
    parsed: &mut dyn FnMut(&str, &str) -> Option<Tree>,
    legend: &mut crate::track::legend::Legend,
) -> Result<Box<dyn Track>, BuildError> {
    let slurped = slurp(spec, context.region, context.width, files, false)?;
    if slurped.origin != Origin::Window {
        return built(spec, context, files, parsed, legend, slurped);
    }
    match built(spec, context, files, parsed, legend, slurped) {
        Err(BuildError::Parse { .. }) => {
            let slurped = slurp(spec, context.region, context.width, files, true)?;
            built(spec, context, files, parsed, legend, slurped)
        }
        other => other,
    }
}

/// [`track`], from the file once it is read.
fn built(
    spec: &TrackSpec,
    context: &Context<'_>,
    files: &mut dyn Files,
    parsed: &mut dyn FnMut(&str, &str) -> Option<Tree>,
    legend: &mut crate::track::legend::Legend,
    slurped: Slurped,
) -> Result<Box<dyn Track>, BuildError> {
    let region = context.region;
    let theme = context.theme;
    let Slurped {
        text,
        path,
        origin,
        probe,
        held: counted,
        reference: native_reference,
        most,
        absent,
    } = slurped;
    // A file named on its own was placed by its name, and a few names hide
    // another format: modkit writes its bedMethyl as `.bed`. Lines written
    // from a binary file say only what they were written as.
    let told = match origin {
        Origin::Text | Origin::Window => spec
            .guessed
            .then(|| refine(spec.kind, probe.as_deref().unwrap_or(&text)))
            .flatten(),
        Origin::Bam | Origin::Native(_) => None,
    };
    let refined;
    let spec = match told {
        Some(kind) if kind != spec.kind => {
            refined = TrackSpec {
                kind,
                ..spec.clone()
            };
            &refined
        }
        _ => spec,
    };
    let name = spec.kind.flag();
    let label = spec.label.clone().or_else(|| default_label(spec));
    let height = spec.height;

    let empty = |wanted: &'static str| match &absent {
        Some(said) => BuildError::Absent {
            track: name,
            path: path.clone(),
            wanted,
            said: said.clone(),
        },
        None => BuildError::Empty {
            track: name,
            path: path.clone(),
            wanted,
        },
    };
    // Read before the match rather than inside the arms, because two arms
    // shadow `text` with a second file of their own and a sheet fetched after
    // that would be read out of the wrong one.
    let sheet = sheet(spec, files)?;

    let built: Box<dyn Track> = match spec.kind {
        Kind::Coverage => {
            // Painted span by span rather than through `from_spans`, which
            // wants the list. `samtools depth` writes a line per base, and a
            // ten million base window cost 231 MB of spans standing beside the
            // 152 MB of text they were read from and the 76 MB track they were
            // about to become.
            let mut painted = CoverageTrack::from_spans(region, std::iter::empty());
            let spans = wrap(
                name,
                &path,
                // Depth worked out from a BAM, and a bigWig's values, arrive
                // as bedGraph, whatever `--format` said about a text file, and
                // a bigWig's spans that overlap are its own, not the sign of
                // two samples' depth the guess takes them for.
                read::signal::fold_spans(
                    &text,
                    region,
                    match origin {
                        Origin::Bam | Origin::Native(_) => Some(crate::Format::BedGraph),
                        Origin::Text | Origin::Window => spec.format,
                    },
                    |start, end, value| painted.paint(start, end, value),
                ),
            )?;
            drop(text);
            // Named on its own, a BAM is its depth; over a window a few reads
            // wide the reads are what a reader came for, and nothing said
            // they could be drawn.
            if origin == Origin::Bam && spec.guessed && region.len() <= READS_WINDOW {
                files.note(&format!(
                    "{path} is drawn as its depth; --pileup {path} draws its reads"
                ));
            }
            if spans == 0 {
                return Err(empty("values"));
            }
            Box::new(named(
                dressed_coverage(painted, spec, most),
                label,
                CoverageTrack::label,
            ))
        }
        Kind::Dynseq => {
            let Some(source) = spec.second.as_ref() else {
                return Err(BuildError::MissingSecond { track: name });
            };
            let reference = second_sequence(name, source, region, files)?;

            let found = wrap(name, &path, read::dynseq::scores(&text, region))?;
            // The file's rows, counted by the index where the text is the
            // window's: a window of none in a file of many is elsewhere.
            if counted.unwrap_or(found.records) == 0 {
                return Err(empty("scores"));
            }
            if found.spans.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "scores",
                    held: counted.unwrap_or(found.records),
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }

            // Padded as far as the scores reach, rather than cut to the
            // reference or stretched to the window. A track only as long as a
            // short FASTA drops a score the reader accepted, without a word;
            // one as long as the window allocates a byte and eight more per
            // base of it, which a sixty byte file across a chromosome should
            // not be able to ask for.
            // Indexed from the start of the window, so a reference that starts
            // later, as a slice may, is padded to it with the letter for a
            // base nobody read.
            let (from, clipped) = reference.clip(region)?;
            let mut letters = vec![b'N'; (from - region.start()) as usize];
            letters.extend(clipped);
            let reach = found
                .spans
                .iter()
                .map(|(_, to, _)| to.saturating_sub(region.start()))
                .max()
                .unwrap_or(0);
            let wanted = usize::try_from(reach).unwrap_or(letters.len());
            letters.resize(wanted.max(letters.len()), b'N');
            let mut track = DynseqTrack::from_spans(region.start(), letters, found.spans);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, DynseqTrack::label))
        }
        Kind::Junctions => {
            let found = wrap(name, &path, read::junction::junctions(&text, region))?;
            if counted.unwrap_or(found.records) == 0 {
                return Err(empty("junctions"));
            }
            if !found.junctions.iter().any(crate::Junction::is_observed) {
                // Three ways a file of junctions reaches no figure, and the
                // counts say which: another sequence, another window, or no
                // read across any of them.
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "junctions",
                    held: counted.unwrap_or(found.records),
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }

            let mut track = JunctionTrack::new(found.junctions);
            if spec.no_counts {
                track = track.show_counts(false);
            }
            if let Some(reads) = spec.min_reads {
                track = track.min_reads(reads);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            if let Some(color) = spec.color.clone() {
                track = track.color(color);
            }
            Box::new(named(track, label, JunctionTrack::label))
        }
        Kind::Sequence => {
            let reference = match native_reference {
                Some(reference) => reference,
                None => sequence(name, &path, &text, region)?,
            };
            let (from, bases) = reference.clip(region)?;
            let mut track = SequenceTrack::new(from, bases);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, SequenceTrack::label))
        }
        Kind::Features => {
            // Each gene once, with the exons all its transcripts use, unless
            // each transcript was asked for.
            let read = if spec.isoforms {
                read::interval::transcripts
            } else {
                read::interval::features
            };
            // A bigBed's rows arrive as BED, cut to the columns that are
            // BED's own, and are read as BED whatever their seventh column
            // holds: guessed at, a dot there says GFF3.
            let format = match origin {
                Origin::Native(_) => Some(crate::Format::Bed),
                _ => spec.format,
            };
            let features = wrap(name, &path, read(&text, region, format))?;
            if features.is_empty() {
                return Err(empty("features"));
            }
            let mut track = FeatureTrack::new(features);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            Box::new(named(track, label, FeatureTrack::label))
        }
        Kind::Variants => {
            let variants = wrap(
                name,
                &path,
                recorded(origin, &text, read::point::variants(&text, region)),
            )?;
            if variants.is_empty() {
                return Err(empty("variants"));
            }
            // A stem is as tall as the allele fraction the VCF's AF gives, and
            // the axis says so; a file with no AF draws no axis to title.
            // Colours are dealt from the most damaging consequence down
            // rather than by which call comes first in the window, which
            // painted one consequence two ways in a gene and a zoom into it.
            let ranked = read::point::ranked(&variants);
            let mut track = VariantTrack::new(variants)
                .axis_title("AF")
                .category_order(ranked);
            if let Some(height) = height {
                track = track.height(height);
            }
            // The library's answer to density, which the command line could not
            // reach before: a genome-wide panel of two hundred thousand calls
            // was fifty megabytes of lollipops nobody could tell apart, and the
            // same panel as ticks is seventy-four kilobytes.
            if let Some(style) = spec.style.and_then(Style::variant) {
                track = track.style(style);
            }
            // Named on its own, a VCF is its calls, which is what a VCF of one
            // sample is for. A cohort's holds a row of calls per sample too, and
            // nothing said they could be drawn. Said once the track is built, so
            // a first attempt that fails and is read again under a --rename
            // says nothing twice, and not where the figure draws them already.
            // A BCF's sites name no sample, and its header, the probe, does.
            let held = if spec.guessed && !context.genotyped {
                read::point::samples(probe.as_deref().unwrap_or(&text)).len()
            } else {
                0
            };
            if held >= 2 {
                files.note(&format!(
                    "{path} is drawn as its calls; --genotypes {path} draws its {} samples, \
                     a row each",
                    crate::track::axis::group_thousands(held as u64)
                ));
            }
            Box::new(named(track, label, VariantTrack::label))
        }
        Kind::Genotypes => {
            let held = read::point::samples(&text);
            // A list, after --genotypes: the rows to draw, in this order. A
            // name the header has not got is nearly always a spelling, so the
            // names it does have are given, the first five of them, since a
            // cohort's header can name thousands.
            let wanted: Option<Vec<String>> = spec.sample.as_ref().map(|list| {
                list.split(',')
                    .map(str::trim)
                    .filter(|name| !name.is_empty())
                    .map(str::to_string)
                    .collect()
            });
            if let Some(missing) = wanted
                .iter()
                .flatten()
                .find(|name| !held.is_empty() && !held.contains(name))
            {
                let named: Vec<(String, usize)> =
                    held.iter().map(|sample| (sample.clone(), 0)).collect();
                return Err(BuildError::Unnamed {
                    track: name,
                    path: path.clone(),
                    what: "sample",
                    wanted: missing.clone(),
                    held: vec![listed_as(&named, "samples")],
                });
            }
            let found = wrap(
                name,
                &path,
                recorded(
                    origin,
                    &text,
                    read::point::genotypes(&text, region, wanted.as_deref()),
                ),
            )?;
            // A cohort's text is the larger of the two by far, and the calls
            // are all of it the track needs.
            drop(text);
            if found.sites.is_empty() {
                return Err(empty("genotypes"));
            }
            let names = found.samples.clone();
            let mut track = GenotypeTrack::new(found.samples, found.sites);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(named(track, label, GenotypeTrack::label))
        }
        Kind::Windows => {
            let windows = wrap(name, &path, read::signal::windows(&text, region))?;
            if windows.is_empty() {
                return Err(empty("windows"));
            }
            Box::new(named(
                dressed_windows(windows, spec),
                label,
                WindowTrack::label,
            ))
        }
        Kind::Manhattan => {
            let table = wrap(name, &path, read::point::association_table(&text, region))?;
            if table.points.is_empty() {
                return Err(empty("association statistics"));
            }
            let points = table.points.clone();
            let mut track = ManhattanTrack::new(points.clone());
            if let Some(source) = spec.second.as_ref() {
                let (text, ld_path) = fetch(name, source, files)?;
                let (pairs, _) = wrap(name, &ld_path, read::pairs::pairs(&text, region))?;
                let Some((lead, linkage)) = lead_linkage(&pairs, &points) else {
                    return Err(BuildError::Empty {
                        track: name,
                        path: ld_path,
                        wanted: "pairs with a tested variant",
                    });
                };
                if !points.iter().any(|point| point.pos == lead) {
                    files.note(&format!(
                        "the lead variant in {ld_path}, {}, is not one {path} tested, so \
                         no diamond marks it",
                        crate::track::axis::group_thousands(lead + 1)
                    ));
                }
                track = track.linkage(lead, linkage);
                // Called by the name the scan gives it, or the one the table
                // of linkage does, as PLINK writes the lead on every row.
                let named = table.name_at(lead).map(str::to_string).or_else(|| {
                    read::pairs::names(&text)
                        .into_iter()
                        .find(|(at, _)| *at == lead)
                        .map(|(_, name)| name)
                });
                if let Some(named) = named {
                    track = track.lead_name(named);
                }
            }
            if let Some(source) = spec.recombination.as_ref() {
                let (text, map_path) = fetch(name, source, files)?;
                let rates = wrap(name, &map_path, read::recombination::rates(&text, region))?;
                if rates.is_empty() {
                    return Err(BuildError::Empty {
                        track: name,
                        path: map_path,
                        wanted: "recombination rates",
                    });
                }
                track = track.recombination(rates);
            }
            // Drawn as -log10, and the axis says so, since the file said p.
            if table.p_values {
                track = track.axis_title("-log10 p");
            }
            let mut track = thresholded(track, spec, table.p_values, &path)?;
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, ManhattanTrack::label))
        }
        Kind::Tree => {
            let tree = read_tree("--tree", &path, &text, parsed, files)?;
            let tree = match &spec.focus {
                None => tree,
                Some(names) => {
                    let mut found = Vec::with_capacity(names.len());
                    for name in names {
                        let Some(node) = tree.node_named(name) else {
                            // Every tip and every labelled clade, so a
                            // misspelling can be seen against what is there
                            // rather than guessed at. Capped, because a
                            // million tip tree would otherwise answer a typo
                            // with a million names.
                            let mut held: Vec<String> =
                                tree.leaf_names().into_iter().take(24).collect();
                            if tree.leaf_count() > held.len() {
                                held.push(format!("and {} more", tree.leaf_count() - held.len()));
                            }
                            return Err(BuildError::Unnamed {
                                track: "tree",
                                path: path.clone(),
                                what: "tip or clade",
                                wanted: name.clone(),
                                held,
                            });
                        };
                        found.push(node);
                    }
                    // One name is that clade, or the clade a tip sits in,
                    // since a tip on its own is not a subtree anyone can read.
                    // Two names are the smallest clade holding both, which is
                    // how a folded triangle names itself in its tooltip.
                    let at = if found.len() == 1 {
                        let node = found[0];
                        if tree.nodes()[node].is_leaf() {
                            tree.nodes()[node].parent.unwrap_or(node)
                        } else {
                            node
                        }
                    } else {
                        tree.mrca(&found).unwrap_or_else(|| tree.root())
                    };
                    tree.subtree(at).unwrap_or(tree)
                }
            };
            // A sheet has to be joined against the tree that will be drawn,
            // which is the one `--focus` left behind and not the one the file
            // held, or the join would be checked against tips this figure does
            // not have.
            let mut tree = tree;
            let leaves = tree.leaf_names();
            // Joined onto the tips by `TreeTrack::traits`, as a library caller
            // joins one, once the track is made.
            let held = strip(spec, sheet.as_ref(), &leaves, context.colors)?;

            // Mutations are branch data the file keeps under a key, and
            // asking who carries one is a question about the shape of the tree.
            // The answer is written back onto the tree as an ordinary
            // annotation, so the colouring, the folding and the strips all read
            // it the way they read anything else and none of them needs to know
            // what a mutation is.
            if let Some(key) = &spec.mutations {
                let found = Mutations::read(&tree, key);
                if found.is_empty() {
                    return Err(BuildError::Empty {
                        track: name,
                        path: path.clone(),
                        wanted: "mutations under that key",
                    });
                }
                if let Some(spelling) = &spec.carrying {
                    let carriers = found.carriers(&tree, spelling);
                    if carriers.is_empty() {
                        let mut held: Vec<String> =
                            found.spellings().take(24).map(str::to_string).collect();
                        if found.distinct() > held.len() {
                            held.push(format!("and {} more", found.distinct() - held.len()));
                        }
                        return Err(BuildError::Unnamed {
                            track: "tree",
                            path: path.clone(),
                            what: "change",
                            wanted: spelling.clone(),
                            held,
                        });
                    }
                    // Only the carriers are marked. Marking both sides makes
                    // two categories, and the colours go in the order they are
                    // first seen: the root does not carry anything, so "does
                    // not" was always seen first and the majority of the tree
                    // came out in the first colour while the answer to the
                    // question sat in the second.
                    for node in carriers {
                        if let Some(into) = tree.annotations_mut(node) {
                            into.insert(
                                spelling.clone(),
                                crate::AnnotationValue::Text("carries".to_string()),
                            );
                        }
                    }
                }
            }

            let mut track = TreeTrack::new(tree);
            // Named before anything else is asked of the track, so a name that
            // is not in the file stops the figure rather than quietly drawing
            // one clade fewer than the command asked for. `highlight_named`
            // does nothing when it cannot find the name, which for a library
            // call is a choice and for a command line is a lie.
            for wanted in &spec.highlight {
                if track.tree().node_named(wanted).is_none() {
                    let mut held: Vec<String> = track
                        .tree()
                        .nodes()
                        .iter()
                        .filter_map(|clade| clade.name.clone())
                        .take(24)
                        .collect();
                    let all = track
                        .tree()
                        .nodes()
                        .iter()
                        .filter(|clade| clade.name.is_some())
                        .count();
                    if all > held.len() {
                        held.push(format!("and {} more", all - held.len()));
                    }
                    return Err(BuildError::Unnamed {
                        track: "tree",
                        path: path.clone(),
                        what: "clade",
                        wanted: wanted.clone(),
                        held,
                    });
                }
                track = track.highlight_named(wanted);
            }
            if let Some(held) = held {
                // Where the palette runs out over a column of the sheet, the
                // line under the tree that says so names the way out.
                track = track
                    .traits(held)
                    .recolour("--colors gives them colours of their own");
            }
            if let Some(projection) = spec.projection {
                track = track.projection(projection);
            }
            // Asking who carries a change is asking for it to be visible, so
            // the answer is coloured by unless the command said to colour by
            // something else.
            let coloured = spec.color_by.clone().or_else(|| spec.carrying.clone());
            if let Some(key) = coloured {
                // A key no node carries colours no branch, and the tree comes
                // out as it would have without the flag. Looked for on every
                // node rather than on the tips, since a branch takes its
                // nearest annotated ancestor's value, and on the tree as it
                // will be drawn: after `--focus` has cut it, and after the
                // sheet and the carriers have been written onto it.
                let tree = track.tree();
                let nodes = 0..tree.nodes().len();
                if !nodes
                    .clone()
                    .any(|node| tree.annotation(node, &key).is_some())
                {
                    let keys: std::collections::BTreeSet<&str> = nodes
                        .filter_map(|node| tree.annotations(node))
                        .flat_map(|held| held.keys().map(String::as_str))
                        .collect();
                    let mut held: Vec<String> =
                        keys.iter().take(24).map(|key| key.to_string()).collect();
                    if keys.len() > held.len() {
                        held.push(format!("and {} more", keys.len() - held.len()));
                    }
                    return Err(BuildError::Unnamed {
                        track: "tree",
                        path: path.clone(),
                        what: "annotation",
                        wanted: key,
                        held,
                    });
                }
                track = track.color_by(key);
            }
            // Refused rather than said under the tree, as a colour key is: a
            // support style drawing nothing on a tree whose support the file
            // keeps under another name is the figure this is asked to avoid.
            if let Some(key) = &spec.support_from {
                let tree = track.tree();
                let carried = (0..tree.nodes().len()).any(|node| {
                    !tree.nodes()[node].is_leaf()
                        && tree
                            .annotation(node, key)
                            .and_then(crate::AnnotationValue::as_number)
                            .is_some()
                });
                if !carried {
                    let keys: std::collections::BTreeSet<&str> = (0..tree.nodes().len())
                        .filter(|node| !tree.nodes()[*node].is_leaf())
                        .filter_map(|node| tree.annotations(node))
                        .flat_map(|held| {
                            held.iter()
                                .filter(|(_, value)| value.as_number().is_some())
                                .map(|(key, _)| key.as_str())
                        })
                        .collect();
                    return Err(BuildError::Unnamed {
                        track: "tree",
                        path: path.clone(),
                        what: "annotation of numbers on a clade",
                        wanted: key.clone(),
                        held: keys.iter().take(24).map(|key| key.to_string()).collect(),
                    });
                }
                track = track.support_from(key);
            }
            if let Some(style) = spec.support_style {
                track = track.support_style(match style {
                    TreeSupport::None => crate::SupportStyle::None,
                    TreeSupport::Symbols => crate::SupportStyle::Symbols,
                    TreeSupport::Labels => crate::SupportStyle::Labels,
                    TreeSupport::Both => crate::SupportStyle::SymbolsAndLabels,
                });
            }
            if let Some(minimum) = spec.threshold {
                track = track.support_threshold(minimum.drawn());
            }
            if spec.no_scale_bar {
                track = track.show_scale_bar(false);
            }
            if spec.cladogram {
                track = track.shape(TreeShape::Cladogram);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            gather(legend, &track.legend(theme));
            // What the tree was asked for and does not draw is said under
            // it, and here too, where a command that ran is read.
            for warning in track.warnings() {
                files.note(&format!("--tree {path}: {warning}"));
            }
            Box::new(named(track, label, TreeTrack::label))
        }
        // Two trees, and the grammar gives one path per flag, so the second
        // arrives by name: --tanglegram left.nwk --against right.nwk. The
        // parser refuses the flag without its pair, and this refuses it again
        // because every field of a TrackSpec is public and the parser is not
        // the only way one is built. Neither refusal is spare: a tanglegram
        // drawn from a single tree against itself has no crossings at all,
        // which is what a perfect result looks like.
        Kind::Tanglegram => {
            let Some(source) = spec.second.as_ref() else {
                return Err(BuildError::MissingSecond { track: name });
            };
            let (other, right_path) = fetch(name, source, files)?;
            let left = read_tree("--tanglegram", &path, &text, parsed, files)?;
            let right = read_tree("--against", &right_path, &other, parsed, files)?;
            // Named, because two phylogenies side by side with nothing over
            // them do not say which is which, and which is which is the whole
            // of what a tanglegram is read for.
            let mut track =
                TanglegramTrack::new(left, right).names(shortened(&path), shortened(&right_path));
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            Box::new(named(track, label, TanglegramTrack::label))
        }
        // One bedMethyl is one track only when it counted one modification. A
        // dual-mode run writes m and h at the same cytosine, and stacked on one
        // axis those are two marks at one position with nothing naming either.
        Kind::Methylation => {
            let Some(code) = selected(name, &path, "--modification", &spec.selects, || {
                wrap(name, &path, read::methyl::codes(&text))
            })?
            else {
                return Err(empty("modified bases"));
            };

            let found = wrap(name, &path, read::methyl::sites(&text, region, &code))?;
            if counted.unwrap_or(found.records) == 0 {
                return Err(empty("modified bases"));
            }
            // Only when the window listed nothing. Positions listed with no
            // valid coverage are in the window rather than elsewhere, and the
            // band counts them in its corner, as it counts the calls a floor
            // hides; refused, they were reported as though the file held its
            // calls somewhere else.
            if found.sites.is_empty() && found.no_coverage == 0 {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "modified bases",
                    held: counted.unwrap_or(found.records),
                    named: code.clone(),
                    region: region.to_string(),
                    rename: None,
                });
            }

            // The reader skips a position nobody could call rather than making
            // a site at nought per cent of it, so its count comes over by hand
            // or the band says nothing about it.
            let mut track = MethylationTrack::new(found.sites).no_coverage(found.no_coverage);
            if let Some(reads) = spec.min_reads {
                track = track.min_coverage(reads);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            // Named after the modification it counted, since the band shows one
            // of the several a file may hold and nothing else would say which.
            Box::new(match label {
                Some(label) => track.label(label),
                None => track.label(code),
            })
        }
        // Calls as arcs between their breakpoints. Every refusal in the reader
        // stands between a broken record and an arc drawn at full confidence.
        Kind::Structural => {
            let found = wrap(
                name,
                &path,
                recorded(origin, &text, read::structural::variants(&text, region)),
            )?;
            if found.records == 0 {
                return Err(empty("variant calls"));
            }
            if found.variants.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "structural calls",
                    held: found.records,
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }
            let mut track = StructuralTrack::new(found.variants);
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, StructuralTrack::label))
        }
        // One row per molecule, its pieces in the order that molecule visited
        // them, which the reader works out rather than takes from the file.
        Kind::SplitReads => {
            let found = wrap(name, &path, read::split::reads(&text, region))?;
            if found.records == 0 {
                return Err(empty("alignments"));
            }
            if found.reads.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "split reads",
                    held: found.records,
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }
            let mut track = SplitReadTrack::new(found.reads);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            Box::new(named(track, label, SplitReadTrack::label))
        }
        // One row per molecule and one column per site. The reader builds the
        // grid by position, since a call written into the wrong column is a
        // methylation pattern that never existed, drawn as cleanly as one that
        // did, and nothing downstream could tell.
        Kind::Bisulfite => {
            let Some(context) = selected(name, &path, "--context", &spec.selects, || {
                wrap(name, &path, read::bisulfite::contexts(&text))
            })?
            else {
                return Err(empty("methylation calls"));
            };

            let found = wrap(
                name,
                &path,
                read::bisulfite::molecules(&text, region, &context),
            )?;
            if found.molecules.is_empty() || found.sites.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "methylation calls",
                    held: found.records,
                    named: context.clone(),
                    region: region.to_string(),
                    rename: None,
                });
            }

            let mut track = BisulfiteTrack::new(found.sites, found.molecules);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            Box::new(match label {
                Some(label) => track.label(label),
                None => track.label(context),
            })
        }
        // Protein domains, on an axis of residues rather than of bases. Column
        // one names the row rather than selecting it, so every protein in the
        // file is drawn and they share one axis, which is what makes a domain
        // gained or lost visible at all.
        Kind::Domains => {
            let held = wrap(name, &path, read::domain::analyses(&text))?;
            if held.is_empty() {
                return Err(empty("domain annotations"));
            }
            let analysis = chosen(name, &path, "--analysis", &held, &spec.selects)?;

            let found = wrap(
                name,
                &path,
                read::domain::architectures(&text, region, &analysis),
            )?;
            // A protein with no annotated domain is a real row, so an empty
            // architecture is not the failure here; a file with no protein in
            // it at all is.
            if found.rows.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "domain annotations",
                    held: found.records,
                    named: analysis.clone(),
                    region: region.to_string(),
                    rename: None,
                });
            }
            if found.rows.iter().all(|row| row.features.is_empty()) {
                return Err(BuildError::Unjoined {
                    track: name,
                    path: path.clone(),
                    what: "domain",
                    against: "the window",
                    examples: found.proteins.iter().take(3).cloned().collect(),
                });
            }

            let names: Vec<String> = found.rows.iter().map(|row| row.name.clone()).collect();
            let mut track = DomainTrack::new(found.rows);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(match label {
                Some(label) => track.label(label),
                None => track.label(analysis),
            })
        }
        // Spans plus the taxa carrying them, painted onto a phylogeny that
        // comes from a second file. Every refusal below stands between a
        // mistake and a figure that looks like a result: a tree with no blocks
        // on it says there was no recombination here, and a tree whose taxa the
        // file never names says the same thing at more length.
        Kind::Clades => {
            let Some(source) = spec.second.as_ref() else {
                return Err(BuildError::MissingSecond { track: name });
            };
            let (newick, tree_path) = fetch(name, source, files)?;
            let tree = read_tree("--with-tree", &tree_path, &newick, parsed, files)?;

            let found = wrap(name, &path, read::clade::blocks(&text, region))?;
            if found.records == 0 {
                return Err(empty("clade blocks"));
            }
            if found.blocks.is_empty() {
                // The file did hold blocks, so say which of the two ways they
                // failed to reach the figure rather than repeating the count.
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "clade blocks",
                    held: found.records,
                    named: found.sequences.join(", "),
                    region: region.to_string(),
                    rename: None,
                });
            }

            // The join is names, and both ways it fails are counted here
            // because the track cannot tell a caller either of them: a block
            // none of whose taxa the tree has reports no unmatched taxa at all,
            // the names inside it having been dropped along with the block.
            let leaves: std::collections::BTreeSet<String> =
                tree.leaf_names().into_iter().collect();
            let carried = found
                .blocks
                .iter()
                .filter(|block| block.taxa().iter().any(|taxon| leaves.contains(taxon)))
                .count();
            if carried == 0 {
                return Err(BuildError::Unjoined {
                    track: name,
                    path: path.clone(),
                    what: "taxon",
                    against: "the phylogeny",
                    examples: found
                        .blocks
                        .iter()
                        .flat_map(|block| block.taxa())
                        .take(3)
                        .cloned()
                        .collect(),
                });
            }

            let names: Vec<String> = leaves.iter().cloned().collect();
            let mut track = CladeTrack::new(tree, found.blocks);
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            Box::new(named(track, label, CladeTrack::label))
        }
        // Gene neighbourhoods from several genomes, and a second file saying
        // what joins one to the next. The links are required: the track marks
        // every gene no homology reaches, so the absence of the file is drawn
        // as the strongest positive finding the track can make.
        Kind::Loci => {
            let Some(source) = spec.second.as_ref() else {
                return Err(BuildError::MissingSecond { track: name });
            };
            let found = wrap(name, &path, read::locus::loci(&text, region, spec.format))?;
            if found.records == 0 {
                return Err(empty("genes"));
            }
            if found.loci.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "genes",
                    held: found.records,
                    named: String::new(),
                    region: region.to_string(),
                    rename: None,
                });
            }

            let (text, link_path) = fetch(name, source, files)?;
            let joined = wrap(
                name,
                &link_path,
                read::locus::links(&text, &found.loci, spec.identity),
            )?;
            if joined.records == 0 {
                return Err(BuildError::Empty {
                    track: name,
                    path: link_path,
                    wanted: "homologies",
                });
            }
            // Nothing joined is the figure this whole arm exists to refuse. It
            // is not an empty plot: it is every gene in every genome outlined
            // as having no counterpart, which reads as a discovery.
            if joined.links.is_empty() {
                return Err(BuildError::Unjoined {
                    track: name,
                    path: link_path,
                    what: "gene name",
                    against: "the loci",
                    examples: joined.unjoined.iter().take(3).cloned().collect(),
                });
            }

            let names: Vec<String> = found.loci.iter().map(|locus| locus.name.clone()).collect();
            let mut track = LocusTrack::new(found.loci).links(joined.links);
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            Box::new(named(track, label, LocusTrack::label))
        }
        // A PAF names both sequences on every row, and an AlignmentBlock keeps
        // neither, so something has to choose which pair the figure is about.
        // The query is the sequence the region is on, which is not a choice.
        // The target is, and guessing it silently would draw a comparison
        // nobody asked for, so the pick is the most-aligned target, it is
        // deterministic, and the ribbon prints both names so the figure says
        // which two sequences it compared rather than leaving it to be assumed.
        Kind::Synteny | Kind::Dotplot => {
            let query = region.seq();
            let found = wrap(name, &path, read::align_pairs::targets(&text, query))?;
            let target = found
                .first()
                .map(|(name, _)| name.clone())
                .ok_or_else(|| empty(spec.kind.flag()))?;
            let alignments = wrap(
                name,
                &path,
                read::align_pairs::blocks(&text, query, &target),
            )?;
            if alignments.blocks.is_empty() {
                return Err(empty(spec.kind.flag()));
            }
            if spec.kind == Kind::Dotplot {
                let mut track = DotplotTrack::new(alignments.blocks);
                if let Some(length) = alignments.target_length {
                    track = track.target_length(length);
                }
                if let Some(height) = height {
                    track = track.height(height);
                }
                Box::new(named(track, label, DotplotTrack::label))
            } else {
                let mut track = SyntenyTrack::new(alignments.blocks).names(query, &target);
                if let Some(length) = alignments.target_length {
                    track = track.target_length(length);
                }
                if let Some(height) = height {
                    track = track.height(height);
                }
                Box::new(named(track, label, SyntenyTrack::label))
            }
        }
        // The two conventions in this file are opposite, and which one a track
        // follows is not a matter of taste. `--sequence` clips its bases to the
        // window and anchors them at the window's start; `--snps` anchors at
        // nought and lets the alignment column be the coordinate. An ORF is
        // read off the reference, so it takes the first; a logo is counted down
        // alignment columns, so it takes the second. Using the window's start
        // for a logo offsets every column by it, and the figure looks fine.
        Kind::Orfs => {
            let reference = match native_reference {
                Some(reference) => reference,
                None => sequence(name, &path, &text, region)?,
            };
            let (from, bases) = reference.clip(region)?;
            let mut track = OrfTrack::new(from, bases);
            if let Some(px) = spec.row_height {
                track = track.lane_height(px);
            }
            Box::new(named(track, label, OrfTrack::label))
        }
        Kind::Logo => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            let rows: Vec<String> = sequences
                .iter()
                .map(|row| String::from_utf8_lossy(&row.residues).into_owned())
                .collect();
            let track = LogoTrack::from_sequences(0, &rows);
            Box::new(named(track, label, LogoTrack::label))
        }
        Kind::Msa => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            let names: Vec<String> = sequences.iter().map(|row| row.name.clone()).collect();
            let mut track = MsaTrack::new(sequences);
            if let Some(index) = compared_row(name, &path, &names, &spec.compare_to)? {
                track = track.compare_to(index);
            }
            if let Some(display) = spec.style.and_then(Style::msa) {
                track = track.display(display);
            }
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(named(track, label, MsaTrack::label))
        }
        Kind::Snps => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            let names: Vec<String> = sequences.iter().map(|row| row.name.clone()).collect();
            let reference = compared_row(name, &path, &names, &spec.compare_to)?.unwrap_or(0);
            let mut track = SnpTrack::from_alignment(reference, &sequences);
            if spec.no_counts {
                track = track.show_counts(false);
            }
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            if spec.no_names {
                track = track.show_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(named(track, label, SnpTrack::label))
        }
        Kind::Ideogram => {
            let (length, bands) = wrap(name, &path, read::interval::cytoband(&text, region.seq()))?;
            if bands.is_empty() {
                return Err(empty("bands"));
            }
            let mut track = IdeogramTrack::new(length, bands);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, IdeogramTrack::label))
        }
        Kind::CopyNumber => {
            // Checked by the parser, so this is reached only from an
            // Invocation built by hand, whose fields are all public.
            let Some(ploidy) = spec.ploidy else {
                return Err(BuildError::MissingPloidy { track: name });
            };
            one_sample(spec, &path, &text)?;
            let found = wrap(
                name,
                &path,
                read::segments::copy_numbers(&text, region, ploidy, spec.sample.as_deref()),
            )?;
            if found.records == 0 {
                return Err(empty("segments"));
            }
            if found.segments.is_empty() {
                // The file did hold segments, so say which of the three ways
                // they failed to reach the figure rather than repeating the
                // count: another sequence, another window, or no call at all.
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "called segments",
                    held: found.records,
                    named: found.samples.join(", "),
                    region: region.to_string(),
                    rename: None,
                });
            }

            let mut track = CopyNumberTrack::at_ploidy(found.segments, ploidy);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, CopyNumberTrack::label))
        }
        Kind::Matrix => {
            let (sites, rows) = wrap(name, &path, read::table::matrix(&text, region))?;
            // No sample lines at all, and a header whose every site lies
            // outside the window, both leave nothing to draw: the second one
            // used to give a lane of names beside no cells.
            if rows.is_empty() || sites.is_empty() {
                return Err(empty("samples"));
            }
            let names: Vec<String> = rows.iter().map(|row| row.name.clone()).collect();
            let mut track = MatrixTrack::new(sites, rows);
            if let Some(max) = spec.max {
                track = track.max(max);
            }
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_row_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(named(track, label, MatrixTrack::label))
        }
        Kind::Recombination => {
            let rates = wrap(name, &path, read::recombination::rates(&text, region))?;
            if rates.is_empty() {
                return Err(empty("rates"));
            }
            // A line, since a rate is read for where it rises rather than for
            // the area under it, and the highest rate in a pixel, so a hotspot
            // narrower than a pixel is still drawn.
            let mut track = CoverageTrack::from_spans(region, rates)
                .aggregate(Aggregate::Max)
                .style(crate::CoverageStyle::Line)
                .axis_title("cM/Mb");
            if let Some(max) = spec.max {
                track = track.max(max);
            }
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, CoverageTrack::label))
        }
        Kind::Pairs => {
            // Only a .hic holds more than one size of bin to pick from.
            let contacts = origin == Origin::Native(Binary::Hic);
            if let (Some(asked), false) = (spec.resolution, contacts) {
                return Err(BuildError::Unresolved {
                    track: name,
                    path,
                    asked,
                    held: None,
                });
            }
            let (pairs, measured) = wrap(name, &path, read::pairs::pairs(&text, region))?;
            if pairs.is_empty() {
                return Err(empty("pairs"));
            }
            // A triangle where most of the pairs the places could make were
            // measured, as linkage and contacts are, and arcs where a few were.
            // Linkage is a triangle however PLINK's window filtered it, and a
            // contact map however few of its cells hold a count: one that
            // holds none is a count of nought, not a pair left unmeasured.
            let correlation = measured.as_deref().is_some_and(read::pairs::is_correlation);
            let style = spec.style.and_then(Style::pairs).unwrap_or_else(|| {
                if correlation || contacts {
                    PairStyle::Triangle
                } else {
                    PairStyle::for_pairs(&pairs)
                }
            });
            let mut track = PairTrack::new(pairs).style(style);
            // An r² of 0.4 is weak linkage whatever else is in the window, so
            // a correlation is read against one rather than its own largest.
            if correlation {
                track = track.ceiling(1.0);
            }
            // After the default, so a ceiling asked for wins over the one an
            // r² is read against until told.
            if let Some(max) = spec.max {
                track = track.ceiling(max);
            }
            if let Some(line) = spec.threshold.map(Threshold::drawn) {
                track = track.threshold(line);
            }
            if spec.log {
                track = track.log_scale(true);
            }
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, PairTrack::label))
        }
        Kind::Heatmap => {
            let (windows, mut rows) = wrap(name, &path, read::table::windows(&text, region))?;
            if rows.is_empty() || windows.is_empty() {
                return Err(empty("windows"));
            }
            // Each sample against its own usual value, so a sample sequenced
            // deeper is not a darker row from end to end and what stands out
            // is what changed along it. Its median over the windows drawn, as
            // the stretches it lost do not move a median the way they move a
            // mean; a sample with no usual value above nought keeps its own.
            if spec.relative {
                for row in &mut rows {
                    let mut drawn: Vec<f64> = row
                        .values
                        .iter()
                        .copied()
                        .filter(|v| v.is_finite())
                        .collect();
                    drawn.sort_by(f64::total_cmp);
                    let median = match drawn.len() {
                        0 => continue,
                        n if n % 2 == 1 => drawn[n / 2],
                        n => (drawn[n / 2 - 1] + drawn[n / 2]) / 2.0,
                    };
                    if median > 0.0 {
                        row.values.iter_mut().for_each(|value| *value /= median);
                    }
                }
            }
            let names: Vec<String> = rows.iter().map(|row| row.name.clone()).collect();
            let mut track = MatrixTrack::windows(windows, rows);
            if spec.relative {
                track = track.unit("×");
            }
            // A depth read against its sample's usual one is read either side
            // of one: a loss in one hue and a gain in the other, and the usual
            // depth pale, where one hue made a deletion as pale as the page.
            if let Some(center) = spec.center.or(spec.relative.then_some(1.0)) {
                track = track.scale(crate::CellScale::Diverging {
                    center,
                    spread: None,
                });
            }
            // The end of the ramp, or of its gain where it is read either
            // side of a centre; the parser has seen to it being above that.
            if let Some(max) = spec.max {
                track = track.max(max);
            }
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_row_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names, context.colors)? {
                gather(legend, &traits.legend(theme));
                track = track.traits(traits);
            }
            if let Some(tree) = row_tree(spec, &names, files, parsed)? {
                track = track.tree(tree);
            }
            Box::new(named(track, label, MatrixTrack::label))
        }
        Kind::Pileup => {
            let reads = wrap(name, &path, read::align::sam(&text, region))?;
            if reads.is_empty() {
                return Err(empty("reads"));
            }
            let mut track = PileupTrack::new(reads);
            if spec.fade_by_mapq {
                track = track.fade_by_quality(true);
            }
            if let Some(px) = spec.row_height {
                track = track.read_height(px);
            }
            // Without this the track draws every base agreeing, because a
            // mismatch is a base that differs from a reference it was never
            // given. The letters are clipped to the window and the start is the
            // window's, which is the same arrangement the sequence and dynseq
            // tracks use, so a read hanging over the left edge is compared
            // against nothing rather than against the wrong base.
            if let Some(source) = spec.second.as_ref().or(context.reference) {
                let (from, bases) = second_sequence(name, source, region, files)?.clip(region)?;
                track = track.reference(from, bases);
            }
            if let Some(cap) = spec.max_rows {
                track = track.max_rows(cap.rows());
            }
            Box::new(named(track, label, PileupTrack::label))
        }
        Kind::Frequencies => {
            let (counts, _) = wrap(name, &path, read::series::counts(&text, context.decimals))?;
            let mut track = SurveillanceTrack::new(counts).time_decimals(context.decimals);
            if let Some(style) = spec.style.and_then(Style::frequencies) {
                track = track.style(style);
            }
            if spec.counts {
                track = track.metric(crate::SurveillanceMetric::Count);
            }
            if let Some(floor) = spec.min_total {
                track = track.minimum_total(floor);
            }
            if let Some(alert) = spec.threshold.map(Threshold::drawn) {
                track = track.frequency_alert(alert);
            }
            if let Some(rise) = spec.growth {
                track = track.growth_alert(rise);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, SurveillanceTrack::label))
        }
        Kind::Phylodynamics => {
            let (points, _) = wrap(
                name,
                &path,
                read::series::estimates(&text, context.decimals),
            )?;
            let mut track = PhylodynamicTrack::new(points).time_decimals(context.decimals);
            if spec.log {
                track = track.scale(PhylodynamicScale::Log10);
            }
            if let Some(line) = spec.threshold.map(Threshold::drawn) {
                track = track.reference(line, crate::svg::text_rounded(line, 5));
            }
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, PhylodynamicTrack::label))
        }
        Kind::Selection => {
            let sites = wrap(name, &path, read::series::selection(&text))?;
            // A p-value where the table has them, as FEL and MEME write, and
            // a posterior where it has only those, as a Bayes empirical
            // Bayes table or FUBAR does.
            let evidence = if sites.iter().all(|site| site.p().is_none())
                && sites.iter().any(|site| site.probability().is_some())
            {
                SelectionEvidence::Posterior
            } else {
                SelectionEvidence::PValue
            };
            let mut track = SelectionTrack::new(sites).evidence(evidence);
            if let Some(line) = spec.threshold.map(Threshold::drawn) {
                track = match evidence {
                    SelectionEvidence::PValue => track.p_threshold(line),
                    SelectionEvidence::Posterior => track.posterior_threshold(line),
                };
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, SelectionTrack::label))
        }
        Kind::Squiggle => {
            let signal = wrap(
                name,
                &path,
                read::series::squiggle(&text, spec.selects.as_deref()),
            )?;
            if signal.reads > 1 && spec.selects.is_none() {
                files.note(&format!(
                    "{path} holds {} reads, and this is the first, {}; --read NAME draws \
                     another",
                    signal.reads, signal.read
                ));
            }
            // The bases the basecaller called, each where its stretch of
            // current starts, from the move table of the read's record.
            let moves = match spec.second.as_ref() {
                Some(source) => {
                    let wanted = signal.named.then_some(signal.read.as_str());
                    let from_bam = files.named_read(source, &signal.read).map_err(|cause| {
                        BuildError::Open {
                            track: name,
                            path: called(source),
                            cause,
                        }
                    })?;
                    let (text, moves_path) = match from_bam {
                        Some(text) => (text, called(source)),
                        None => fetch(name, source, files)?,
                    };
                    Some(wrap(name, &moves_path, read::series::moves(&text, wanted))?)
                }
                None => None,
            };
            // Named for its read, which is what tells two squiggles apart.
            let label = spec.label.clone().or(Some(signal.read));
            let mut track = SquiggleTrack::new(0, signal.samples);
            if let Some(moves) = moves {
                track = track.moves(moves);
            }
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, SquiggleTrack::label))
        }
        Kind::Axis => unreachable!("the ruler is added by build"),
        Kind::Codons => unreachable!("the codon ruler is built by build_one"),
    };
    Ok(built)
}

/// The lead variant of a table of linkage, and the r² of every other variant
/// with it.
///
/// PLINK's `--ld-snp` writes the lead into every row, so a variant in every
/// pair is the lead. A table of every pair names none, and the lead is then the
/// strongest tested variant the table has pairs for. `None` where no pair
/// names a variant.
fn lead_linkage(
    pairs: &[crate::Pair],
    points: &[crate::Association],
) -> Option<(u64, Vec<(u64, f64)>)> {
    let place = |span: (u64, u64)| span.0;
    let first = pairs.first()?;
    let shared = [place(first.first), place(first.second)]
        .into_iter()
        .find(|candidate| {
            pairs
                .iter()
                .all(|pair| place(pair.first) == *candidate || place(pair.second) == *candidate)
        });
    let lead = shared.or_else(|| {
        let paired: std::collections::BTreeSet<u64> = pairs
            .iter()
            .flat_map(|pair| [place(pair.first), place(pair.second)])
            .collect();
        points
            .iter()
            .filter(|point| paired.contains(&point.pos) && point.value.is_finite())
            .max_by(|a, b| a.value.total_cmp(&b.value))
            .map(|point| point.pos)
    })?;
    let mut linkage: Vec<(u64, f64)> = pairs
        .iter()
        .filter_map(|pair| {
            let (a, b) = (place(pair.first), place(pair.second));
            match (a == lead, b == lead) {
                (true, false) => Some((b, pair.value)),
                (false, true) => Some((a, pair.value)),
                _ => None,
            }
        })
        .collect();
    linkage.push((lead, 1.0));
    Some((lead, linkage))
}

/// Wraps a reader error with the flag and the file that produced it.
fn wrap<T>(
    track: &'static str,
    path: &str,
    result: Result<T, read::ReadError>,
) -> Result<T, BuildError> {
    result.map_err(|cause| BuildError::Parse {
        track,
        path: path.to_string(),
        cause,
    })
}

/// A refusal of a row of the text written from a BCF, which is no line of the
/// file, as a refusal of the record it was written from: named by its
/// sequence and position, as `bcftools view -r` would find it, where a VCF's
/// is named by its line. A refusal of the header loses its line alone.
fn recorded<T>(
    origin: Origin,
    text: &str,
    result: Result<T, read::ReadError>,
) -> Result<T, read::ReadError> {
    if origin != Origin::Native(Binary::Bcf) {
        return result;
    }
    result.map_err(|error| {
        let row = error
            .line
            .checked_sub(1)
            .and_then(|at| text.lines().nth(at))
            .filter(|row| !row.starts_with('#'));
        let mut fields = row.unwrap_or_default().split('\t');
        match (
            fields.next(),
            fields.next().and_then(|pos| pos.parse::<u64>().ok()),
        ) {
            (Some(sequence), Some(pos)) => read::ReadError::whole(format!(
                "the record at {sequence}:{}: {}",
                crate::track::axis::group_thousands(pos),
                error.reason
            )),
            _ => read::ReadError::whole(error.reason),
        }
    })
}

/// Applies `--label`, or leaves the track unnamed.
fn named<T>(track: T, label: Option<String>, set: fn(T, String) -> T) -> T {
    match label {
        Some(label) => set(track, label),
        None => track,
    }
}

/// Turns parsed records into alignment rows, refusing an empty file.
fn msa(
    records: Vec<(String, Vec<u8>)>,
    empty: &dyn Fn(&'static str) -> BuildError,
) -> Result<Vec<MsaSequence>, BuildError> {
    if records.is_empty() {
        return Err(empty("sequences"));
    }
    Ok(records
        .into_iter()
        .map(|(name, residues)| MsaSequence::new(name, residues))
        .collect())
}

/// The reference letters a second FASTA holds for this region, picked out of
/// it the way [`sequence`] picks them.
fn second_sequence(
    track: &'static str,
    source: &Source,
    region: &Region,
    files: &mut dyn Files,
) -> Result<Reference, BuildError> {
    // A 2bit is read a window at a time, and any other binary format a
    // track reads as it is holds no bases to read.
    let opened = native(files, source).map_err(|cause| BuildError::Open {
        track,
        path: called(source),
        cause,
    })?;
    match opened {
        Some((Binary::TwoBit, mut file)) => {
            return two_bit(track, &called(source), file.as_mut(), region)
        }
        Some((binary, _)) => {
            return Err(BuildError::OtherTrack {
                track,
                path: called(source),
                binary,
                format: false,
            })
        }
        None => {}
    }
    let (fasta, path) = fetch(track, source, files)?;
    sequence(track, &path, &fasta, region)
}

/// The record of a FASTA that a figure over this region is drawn from.
///
/// One record in the file is the record the file is about, whatever its header
/// calls it. More than one is a genome, and then the region picks: one
/// chromosome's letters under another chromosome's scores, or beside another
/// chromosome's reads, is a figure that is wrong everywhere and looks right. A
/// name two records share picks neither, since the first of them is a sequence
/// nobody chose.
///
/// Every flag that reads a FASTA for its bases comes through here, since a
/// command may hand the same file to more than one of them and it has to be
/// the same sequence to each. `--sequence` and `--orfs` took the first record
/// instead, so a region on one chromosome drew another's bases above a pileup
/// read against the right one.
fn sequence(
    track: &'static str,
    path: &str,
    fasta: &str,
    region: &Region,
) -> Result<Reference, BuildError> {
    let mut records: Vec<Reference> = wrap(track, path, read::seq::fasta(fasta))?
        .into_iter()
        .map(|(name, bases)| Reference::new(track, path, name, bases))
        .collect();
    if records.is_empty() {
        return Err(BuildError::Empty {
            track,
            path: path.to_string(),
            wanted: "sequence",
        });
    }
    if records.len() == 1 {
        return Ok(records.swap_remove(0));
    }
    let named: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|(_, record)| record.sequence == region.seq())
        .map(|(index, _)| index)
        .collect();
    // Several slices of one sequence are one sequence in pieces, and the
    // window picks the piece: only the ones it touches are candidates. Two
    // whole records sharing a name are still refused, since the first of them
    // is a sequence nobody chose.
    let slices = named.iter().all(|at| records[*at].is_slice());
    let found: Vec<usize> = if named.len() > 1 && slices {
        named
            .iter()
            .copied()
            .filter(|at| records[*at].touches(region))
            .collect()
    } else {
        named
    };
    match found.as_slice() {
        [index] => Ok(records.swap_remove(*index)),
        [] => {
            // Every name, so a sequence spelt two ways can be seen against
            // what the file calls it. Capped, as a tree's tips are, because a
            // draft assembly would otherwise answer a misspelling with every
            // contig it has.
            let mut held: Vec<String> = records
                .iter()
                .take(24)
                .map(|record| record.name.clone())
                .collect();
            if records.len() > held.len() {
                held.push(format!("and {} more", records.len() - held.len()));
            }
            Err(BuildError::Unnamed {
                track,
                path: path.to_string(),
                what: "record",
                wanted: region.seq().to_string(),
                held,
            })
        }
        many => Err(BuildError::Repeated {
            track,
            path: path.to_string(),
            what: "record",
            name: region.seq().to_string(),
            held: many.len(),
        }),
    }
}

/// One record of a FASTA: its bases, and where the first of them sits.
struct Reference {
    /// Which flag read it and from where, for a refusal to name.
    track: &'static str,
    path: String,
    /// The header's name, as the file gives it.
    name: String,
    /// The sequence the bases belong to: the name, or the part of it before
    /// the span when the record is a slice.
    sequence: String,
    /// The 0-based position of the first base on that sequence.
    offset: u64,
    bases: Vec<u8>,
}

impl Reference {
    /// A record, read as a slice where its header says it is one.
    ///
    /// `samtools faidx ref.fa chr1:101-160` writes the sixty bases under the
    /// header `>chr1:101-160`, and the only way to know they start at 101 is
    /// to read the header. They were drawn from base 1 instead, so a window on
    /// 101 to 160 held no base of a file that held exactly that window, and the
    /// band came out empty. A header is read as a span only when the span is
    /// as long as the record, so a sequence whose own name merely looks like
    /// one keeps its name and its bases where they were.
    fn new(track: &'static str, path: &str, name: String, bases: Vec<u8>) -> Self {
        let slice = name.rsplit_once(':').and_then(|(sequence, span)| {
            let (first, last) = span.split_once('-')?;
            let first: u64 = first.parse().ok()?;
            let last: u64 = last.parse().ok()?;
            (!sequence.is_empty()
                && first >= 1
                && last >= first
                && last - first + 1 == bases.len() as u64)
                .then(|| (sequence.to_string(), first - 1))
        });
        let (sequence, offset) = slice.unwrap_or_else(|| (name.clone(), 0));
        Reference {
            track,
            path: path.to_string(),
            name,
            sequence,
            offset,
            bases,
        }
    }

    /// Whether the header named a span, as `samtools faidx` writes one.
    fn is_slice(&self) -> bool {
        self.name != self.sequence
    }

    /// Whether any base of the record is in the window.
    fn touches(&self, region: &Region) -> bool {
        let end = self.offset + self.bases.len() as u64;
        self.offset < region.end() && end > region.start()
    }

    /// The bases the window holds, and the position of the first of them.
    ///
    /// Refused when it holds none. A reference that ends before the window
    /// starts drew an empty band, with every letter, frame and mismatch the
    /// track exists for missing, and the command exited nought. A window that
    /// runs past the end draws the bases there are, which is what the
    /// sequence says.
    fn clip(&self, region: &Region) -> Result<(u64, Vec<u8>), BuildError> {
        if !self.touches(region) {
            return Err(BuildError::Beyond {
                track: self.track,
                path: self.path.clone(),
                record: self.name.clone(),
                first: self.offset + 1,
                last: self.offset + self.bases.len() as u64,
                region: region.to_string(),
            });
        }
        let from = region.start().max(self.offset);
        let to = region.end().min(self.offset + self.bases.len() as u64);
        let bases = self.bases[(from - self.offset) as usize..(to - self.offset) as usize].to_vec();
        Ok((from, bases))
    }
}

/// Turns `--compare-to` into a row index, or refuses.
///
/// Refused here rather than passed on, because both builders take a `usize` and
/// neither complains about one out of range: an alignment falls back to the
/// consensus and a variable-site panel comes out empty, and both look like
/// figures. A misspelt name would be answered with a confident wrong picture.
fn compared_row(
    track: &'static str,
    path: &str,
    names: &[String],
    wanted: &Option<String>,
) -> Result<Option<usize>, BuildError> {
    let Some(wanted) = wanted else {
        return Ok(None);
    };
    let found: Vec<usize> = names
        .iter()
        .enumerate()
        .filter(|(_, name)| *name == wanted)
        .map(|(index, _)| index)
        .collect();
    match found.as_slice() {
        [index] => Ok(Some(*index)),
        [] => Err(BuildError::Unnamed {
            track,
            path: path.to_string(),
            what: "row",
            wanted: wanted.clone(),
            held: names.to_vec(),
        }),
        many => Err(BuildError::Repeated {
            track,
            path: path.to_string(),
            what: "row",
            name: wanted.clone(),
            held: many.len(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A PAF names two sequences on every row and an AlignmentBlock keeps
    /// neither, so a whole-genome file would otherwise stack alignments
    /// against different chromosomes on one axis and say nothing. This pins
    /// the three things that stops: only the chosen pair is drawn, the choice
    /// is the most-aligned target and is deterministic, and the ribbon prints
    /// both names so the figure states which two it compared.
    #[test]
    fn a_paf_naming_several_targets_draws_one_pair_and_names_it() {
        const PAF: &str = "\
ctg1\t5000\t0\t1200\t+\tchrA\t9000\t400\t1600\t1150\t1200\t60
ctg1\t5000\t1500\t2600\t-\tchrA\t9000\t3000\t4100\t1000\t1100\t60
ctg1\t5000\t100\t400\t+\tchrB\t4000\t50\t350\t280\t300\t60
ctg2\t2000\t0\t900\t+\tchrA\t9000\t100\t1000\t880\t900\t60
";
        let svg = build(&over("ctg1:1-5000", "--synteny", "a.paf"), |_| {
            Ok(PAF.to_string())
        })
        .unwrap();

        // chrA has two rows for this query and chrB one, so chrA is drawn and
        // the figure says so. chrB is not mentioned, which is the whole point.
        assert!(svg.contains("chrA"), "the ribbon does not name its target");
        assert!(
            !svg.contains("chrB"),
            "a second target reached a figure about the first"
        );

        // And the reader agrees about what it kept and what it passed over.
        let found = crate::read::align_pairs::blocks(PAF, "ctg1", "chrA").unwrap();
        assert_eq!(found.blocks.len(), 2);
        assert_eq!(
            found.passed_over, 2,
            "rows of another pair went unmentioned"
        );
        assert_eq!(
            found.target_length,
            Some(9000),
            "the target length has to come from column seven, not from the blocks"
        );
    }

    /// Both new flags follow a coordinate convention, and they are opposite
    /// ones. A logo anchored where the sequence anchors is offset by the whole
    /// window, and the figure looks perfectly reasonable, so this pins each to
    /// the library call it must agree with rather than to a shape.
    #[test]
    fn orfs_and_logos_keep_the_conventions_their_data_has() {
        use crate::{Figure, LogoTrack, OrfTrack};

        let bases: Vec<u8> = b"ACGT".iter().cycle().take(400).copied().collect();
        let fasta = format!(">ctg1\n{}\n", String::from_utf8_lossy(&bases));
        let region = Region::parse("ctg1:101-400").unwrap();

        let from_cli = build(&over("ctg1:101-400", "--orfs", "in.fa"), |_| {
            Ok(fasta.clone())
        })
        .unwrap();
        // The reference is read off the window, so it anchors at the window.
        // The track is called after its file, as every track given no label is.
        let from_library = Figure::new(region.clone())
            .push(OrfTrack::new(100, bases[100..400].to_vec()).label("in"))
            .push(crate::AxisTrack::new())
            .to_svg();
        assert_eq!(from_cli, from_library, "an ORF track moved off its window");

        let rows = ["ACGTACGTAC", "ACGTTCGTAC", "ACGAACGTAC"];
        let aligned = rows
            .iter()
            .enumerate()
            .map(|(i, r)| format!(">s{i}\n{r}\n"))
            .collect::<String>();
        let from_cli = build(&over("aln:1-10", "--logo", "aln.fa"), |_| {
            Ok(aligned.clone())
        })
        .unwrap();
        // A column of an alignment is its own coordinate, so it anchors at
        // nought, and the ruler under it counts columns rather than bases.
        let columns = || crate::AxisTrack::new().counting().label("column");
        let right = Figure::new(Region::parse("aln:1-10").unwrap())
            .push(LogoTrack::from_sequences(0, &rows).label("aln"))
            .push(columns())
            .to_svg();
        let wrong = Figure::new(Region::parse("aln:1-10").unwrap())
            .push(LogoTrack::from_sequences(1, &rows).label("aln"))
            .push(columns())
            .to_svg();
        assert_eq!(from_cli, right, "a logo moved off its columns");
        assert_ne!(right, wrong, "the two anchors are indistinguishable here");
    }
    use crate::cli::args::{parse, Request};

    fn invocation(line: &str) -> Invocation {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    /// A file with `text` in it, named after the test that asked for it.
    ///
    /// The path is not split on whitespace, since a temporary directory may
    /// have a space in its name.
    fn written(name: &str, text: &str) -> String {
        let path = std::env::temp_dir().join(format!("karyon-{}-{}", std::process::id(), name));
        fs::write(&path, text).unwrap();
        path.display().to_string()
    }

    /// A file on disk is read again rather than kept, so it is not held twice
    /// while its track is built. What cannot be read again, a pipe or a
    /// device, is kept once read.
    #[test]
    fn a_file_on_disk_is_read_again_and_a_pipe_is_kept() {
        let path = written("again.bed", "chr1\t0\t10\n");
        let mut disk = Disk::default();
        let file = Source::Path(path.clone().into());
        assert_eq!(disk.text(&file).unwrap(), "chr1\t0\t10\n");
        assert!(disk.kept.is_empty(), "a file on disk was kept");
        fs::remove_file(&path).unwrap();
        // Windows has no /dev, and no shell there names a pipe by a path a
        // Windows program can open, so there is no device to read twice.
        // `/dev/null` there is a file on the current drive that is not there.
        if cfg!(unix) {
            let device = Source::Path("/dev/null".into());
            assert_eq!(disk.text(&device).unwrap(), "");
            assert!(disk.kept.contains_key(&device), "a device was not kept");
        }
    }

    /// A pipe the shell named is asked whether it is a BAM before its text is
    /// read: by a figure placed by a gene's name, and by every track that can
    /// draw a BAM. Answered by reading the pipe, the first bytes of its text
    /// were gone, and a cohort's VCF cut with `tabix -h` lost its header. A
    /// pipe is no BAM read a window at a time, and nothing of it is read to
    /// say so.
    #[cfg(unix)]
    #[test]
    fn asking_whether_a_pipe_is_a_bam_reads_none_of_it() {
        use std::os::unix::io::AsRawFd;
        use std::process::{Command, Stdio};
        let mut child = Command::new("printf")
            .arg("##fileformat=VCFv4.2\\n")
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let pipe = child.stdout.take().unwrap();
        let named = Source::Path(format!("/dev/fd/{}", pipe.as_raw_fd()).into());
        let region = Region::parse("chr1:1-100").unwrap();
        let mut disk = Disk::default();
        assert_eq!(disk.sequences(&named).unwrap(), None);
        assert_eq!(disk.depth(&named, &region).unwrap(), None);
        assert_eq!(disk.reads(&named, &region).unwrap(), None);
        assert_eq!(disk.named_read(&named, "read_1").unwrap(), None);
        assert_eq!(disk.text(&named).unwrap(), "##fileformat=VCFv4.2\n");
        child.wait().unwrap();
    }

    /// A dynseq draws letters from a reference, and the reference has to be the
    /// one the region names.
    ///
    /// One chromosome's letters under another chromosome's scores is a figure
    /// that is wrong at every base and looks right at all of them.
    #[test]
    fn a_reference_of_several_records_is_picked_by_the_region_and_not_by_order() {
        let genome = ">chr1\nAAAAAAAAAAAAAAAAAAAA\n>chr2\nGGGGGGGGGGGGGGGGGGGG\n";
        let scores = "chr2\t0\t8\t0.5\n";

        let args: Vec<String> = "chr2:1-20 --dynseq s.bg --with-sequence g.fa"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let svg = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("g.fa") => genome.to_string(),
                _ => scores.to_string(),
            })
        })
        .unwrap();

        assert!(svg.contains("G at "), "chr2 was not the record drawn");
        assert!(
            !svg.contains("A at "),
            "chr1's letters reached a chr2 figure"
        );
    }

    #[test]
    fn a_reference_naming_none_of_the_region_is_refused_rather_than_guessed() {
        let genome = ">chr1\nAAAAAAAAAAAAAAAAAAAA\n>chr3\nCCCCCCCCCCCCCCCCCCCC\n";
        let args: Vec<String> = "chr2:1-20 --dynseq s.bg --with-sequence g.fa"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let error = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("g.fa") => genome.to_string(),
                _ => "chr2\t0\t8\t0.5\n".to_string(),
            })
        })
        .unwrap_err();
        assert!(
            matches!(error, BuildError::Unnamed { what: "record", .. }),
            "{error}"
        );
    }

    /// `--sequence` and `--orfs` read the FASTA `--with-sequence` reads, and
    /// have to pick the same record out of it.
    ///
    /// Both took the first record in the file, so a region on chrB drew chrA's
    /// bases, and chrA's reading frames, at chrB's coordinates. Each is pinned
    /// to the library call drawing chrB's own bases, since a figure from the
    /// wrong record is a perfectly good figure.
    #[test]
    fn sequence_and_orfs_pick_the_record_the_region_names_and_not_the_first() {
        use crate::{Figure, OrfTrack, SequenceTrack};

        let first: Vec<u8> = b"ACGT".iter().cycle().take(400).copied().collect();
        let named: Vec<u8> = b"ATGGCCTAA".iter().cycle().take(400).copied().collect();
        let genome = format!(
            ">chrA\n{}\n>chrB\n{}\n",
            String::from_utf8_lossy(&first),
            String::from_utf8_lossy(&named)
        );
        let region = Region::parse("chrB:101-400").unwrap();
        let drawn = |flag: &str| {
            build(&over("chrB:101-400", flag, "genome.fa"), |_| {
                Ok(genome.clone())
            })
            .unwrap()
        };

        let bases = |record: &[u8]| {
            let figure = Figure::new(region.clone())
                .push(SequenceTrack::new(100, record[100..400].to_vec()).label("genome"))
                .push(crate::AxisTrack::new());
            // At this zoom the bases are blocks, which the command line keys.
            let key = figure.key();
            figure
                .push(crate::track::legend::LegendTrack::new(key))
                .to_svg()
        };
        assert_ne!(bases(&named), bases(&first), "the records draw alike here");
        assert_eq!(
            drawn("--sequence"),
            bases(&named),
            "--sequence drew another record than the one the region names"
        );

        let frames = |record: &[u8]| {
            Figure::new(region.clone())
                .push(OrfTrack::new(100, record[100..400].to_vec()).label("genome"))
                .push(crate::AxisTrack::new())
                .to_svg()
        };
        assert_ne!(
            frames(&named),
            frames(&first),
            "the records read alike here"
        );
        assert_eq!(
            drawn("--orfs"),
            frames(&named),
            "--orfs read another record than the one the region names"
        );
    }

    /// A genome holding no record of the region's name is refused with what
    /// it does hold, the way a column or a row that is not there is refused.
    /// It is nearly always one sequence spelt two ways, as `chr1` and `1`, and
    /// the two side by side are what the reader needs to see.
    #[test]
    fn a_fasta_naming_none_of_the_region_is_refused_with_the_records_it_holds() {
        let genome = ">chrA\nAAAAAAAAAAAAAAAAAAAA\n>chrC\nCCCCCCCCCCCCCCCCCCCC\n";
        let open = |source: &Source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("genome.fa") => genome.to_string(),
                _ => "chrB\t0\t8\t0.5\n".to_string(),
            })
        };
        for (line, flag) in [
            ("chrB:1-20 --sequence genome.fa", "--sequence"),
            ("chrB:1-20 --orfs genome.fa", "--orfs"),
            (
                "chrB:1-20 --dynseq s.bg --with-sequence genome.fa",
                "--dynseq",
            ),
        ] {
            let refused = build(&invocation(line), open).unwrap_err().to_string();
            assert_eq!(
                refused,
                format!("{flag} genome.fa has no record called chrB; it has chrA, chrC"),
                "{line}"
            );
        }

        // A draft assembly answers a misspelling with two dozen of its contigs
        // and a count of the rest, rather than with every one of them.
        let draft: String = (1..=30).map(|n| format!(">contig_{n}\nACGT\n")).collect();
        let refused = build(&over("chrB:1-20", "--sequence", "draft.fa"), |_| {
            Ok(draft.clone())
        })
        .unwrap_err()
        .to_string();
        assert!(refused.ends_with("contig_24, and 6 more"), "{refused}");
    }

    /// A name two records share picks neither of them. The first would be a
    /// sequence nobody chose, drawn looking exactly like the one asked for,
    /// which is why a row two records share is refused as well.
    #[test]
    fn a_name_two_records_share_is_refused_rather_than_taken_the_first_of() {
        let twice = ">chrB\nAAAAAAAAAAAAAAAAAAAA\n>chrB\nGGGGGGGGGGGGGGGGGGGG\n";
        let refused = build(&over("chrB:1-20", "--sequence", "genome.fa"), |_| {
            Ok(twice.to_string())
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            refused,
            "--sequence genome.fa has 2 records called chrB, so the name does not pick one"
        );
    }

    /// One record is the record the file is about, whatever its header calls
    /// it, which is the rule `--with-sequence` already kept. The three flags
    /// that read a reference read it alike, so a file written as `>contig_1`
    /// under a region on chrB still draws under each of them.
    #[test]
    fn a_fasta_of_one_record_is_drawn_whatever_its_header_calls_it() {
        let single = ">contig_1\nGGGGGGGGGGGGGGGGGGGG\n";
        let open = |source: &Source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("ref.fa") => single.to_string(),
                _ => "chrB\t0\t8\t0.5\n".to_string(),
            })
        };
        for line in [
            "chrB:1-20 --sequence ref.fa",
            "chrB:1-20 --orfs ref.fa",
            "chrB:1-20 --dynseq s.bg --with-sequence ref.fa",
        ] {
            let drawn = build(&invocation(line), open);
            assert!(drawn.is_ok(), "{line}: {drawn:?}");
        }
        let svg = build(&invocation("chrB:1-20 --sequence ref.fa"), open).unwrap();
        assert!(svg.contains(">G</text>"), "the one record was not drawn");
    }

    #[test]
    fn a_segment_table_with_no_sample_column_refuses_the_flag_that_picks_one() {
        // Accepted and ignored, the whole table would be drawn under a name the
        // command asked to pick out of it.
        let table = "chromosome\tstart\tend\tcn\nchr1\t0\t500\t2\n";
        let args: Vec<String> = "chr1:1-1000 --copy-number s.cns --ploidy 2 --sample S1"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let error = build(&invocation, |_| Ok(table.to_string())).unwrap_err();
        assert!(
            matches!(
                error,
                BuildError::Ambiguous {
                    flag: "--sample",
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn a_copy_number_track_built_by_hand_with_no_ploidy_is_refused() {
        // The parser refuses it, so this is reached only from an Invocation
        // built by hand, and defaulting it would draw a confident ladder whose
        // rule came from nowhere.
        let mut invocation = invocation("chr1:1-1000 --coverage d.bg");
        invocation.tracks[0].kind = Kind::CopyNumber;
        invocation.tracks[0].ploidy = None;
        let error = build(&invocation, |_| {
            Ok("chromosome\tstart\tend\tcn\nchr1\t0\t500\t2\n".to_string())
        })
        .unwrap_err();
        assert!(matches!(error, BuildError::MissingPloidy { .. }), "{error}");
    }

    /// A sheet of metadata is a third file, and the join is names.
    ///
    /// The refusal is what these are mostly about. A sheet whose names are
    /// nobody's draws a strip of empty outlines beside every row, and that is
    /// a figure that looks finished and says nothing about anything.
    #[test]
    fn a_sheet_of_metadata_becomes_a_strip_beside_the_rows() {
        let genotypes = "sample\t10\t20\nA\t1\t0\nB\t0\t1\n";
        let sheet = "sample\tlineage\tdepth\nA\tL4\t30\nB\tL2\t50\n";

        let args: Vec<String> = "chr1:1-40 --matrix m.tsv --traits s.tsv"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let svg = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("s.tsv") => sheet.to_string(),
                _ => genotypes.to_string(),
            })
        })
        .unwrap();

        assert!(svg.contains("A; lineage L4"), "{svg}");
        assert!(svg.contains("B; depth 50"), "{svg}");
        // The heading of a column is drawn on end, since a column is narrower
        // than its name and will stay that way.
        assert!(svg.contains("rotate(-90)"), "no heading on the strip");
    }

    #[test]
    fn a_sheet_that_names_none_of_the_rows_is_refused() {
        let genotypes = "sample\t10\t20\nA\t1\t0\nB\t0\t1\n";
        let sheet = "sample\tlineage\nERR1\tL4\nERR2\tL2\n";

        let args: Vec<String> = "chr1:1-40 --matrix m.tsv --traits s.tsv"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let error = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("s.tsv") => sheet.to_string(),
                _ => genotypes.to_string(),
            })
        })
        .unwrap_err();

        assert!(
            matches!(error, BuildError::Unjoined { what: "name", .. }),
            "{error}"
        );
        // The names it did hold, so the mismatch can be seen without opening
        // either file.
        assert!(error.to_string().contains("ERR1"), "{error}");
    }

    #[test]
    fn a_column_the_sheet_has_not_got_is_refused_and_the_ones_it_has_are_named() {
        let genotypes = "sample\t10\t20\nA\t1\t0\nB\t0\t1\n";
        let sheet = "sample\tlineage\tdepth\nA\tL4\t30\nB\tL2\t50\n";

        let args: Vec<String> = "chr1:1-40 --matrix m.tsv --traits s.tsv --columns linage"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let error = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("s.tsv") => sheet.to_string(),
                _ => genotypes.to_string(),
            })
        })
        .unwrap_err();

        let said = error.to_string();
        assert!(said.contains("no column called linage"), "{said}");
        assert!(said.contains("lineage, depth"), "{said}");
    }

    #[test]
    fn the_columns_asked_for_are_the_only_ones_drawn() {
        let genotypes = "sample\t10\t20\nA\t1\t0\nB\t0\t1\n";
        let sheet = "sample\tlineage\tdepth\nA\tL4\t30\nB\tL2\t50\n";

        let args: Vec<String> = "chr1:1-40 --matrix m.tsv --traits s.tsv --columns depth"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let svg = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("s.tsv") => sheet.to_string(),
                _ => genotypes.to_string(),
            })
        })
        .unwrap();

        assert!(svg.contains("A; depth 30"), "{svg}");
        assert!(
            !svg.contains("lineage"),
            "a column nobody asked for was drawn"
        );
    }

    /// The whole of what the second path buys. Two different trees have to
    /// reach the two sides, so the closure answers a different Newick per
    /// path: hand both sides the same text and the crossing count is nought,
    /// which is indistinguishable from a real result.
    #[test]
    fn a_tanglegram_reads_two_files_and_puts_each_on_its_own_side() {
        let args: Vec<String> = "chr1:1-1000 --tanglegram before.nwk --against after.nwk"
            .split_whitespace()
            .map(String::from)
            .collect();
        let invocation = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };

        let svg = build(&invocation, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("before.nwk") => "((a,b),(c,d));",
                Source::Path(path) if path.ends_with("after.nwk") => "((a,c),(b,d));",
                other => panic!("asked for a file nobody named: {other:?}"),
            }
            .to_string())
        })
        .unwrap();

        // Both files were opened, and the figure says which side each is.
        assert!(svg.contains("before.nwk"), "the left tree is unnamed");
        assert!(svg.contains("after.nwk"), "the right tree is unnamed");

        // b and c swap between the two topologies, so exactly one tie crosses.
        // Reading the same file twice would report none.
        assert!(
            svg.contains("1 crossing") && !svg.contains("0 crossing"),
            "the two sides are not the two trees that were given"
        );
    }

    /// Every field of a TrackSpec is public, so the parser is not the only way
    /// one is built, and the playground and any other library caller will
    /// build them directly. A tanglegram short of its second tree must be an
    /// error here too rather than an unwrap.
    #[test]
    fn a_tanglegram_built_without_its_second_tree_is_refused_and_not_a_panic() {
        let mut invocation = over("chr1:1-1000", "--tree", "t.nwk");
        invocation.tracks[0].kind = Kind::Tanglegram;
        let error = build(&invocation, |_| Ok("((a,b),(c,d));".to_string())).unwrap_err();
        assert!(matches!(error, BuildError::MissingSecond { .. }), "{error}");
    }

    /// The two refusals these arms exist for, and both are about a figure
    /// that looks finished. A clade file whose taxa the tree has never heard
    /// of draws a phylogeny with nothing on it, which reads as no
    /// recombination; a homology file whose names did not join draws every
    /// gene in every genome outlined as unique, which reads as a discovery.
    #[test]
    fn two_files_that_name_nothing_in_each_other_are_refused_rather_than_drawn() {
        const GFF: &str = "SEQUENCE\t.\tCDS\t100\t900\t.\t.\t0\ttaxa=\"s1 s2\";\n";
        const BED: &str = "\
A\t0\t400\tg1\t0\t+
B\t0\t400\tg2\t0\t+
";

        // The clade side: a tree of a, b, c against a file naming s1 and s2.
        let it = pair(
            "chr1:1-1000",
            "--clades",
            "clades.gff",
            "--with-tree",
            "tree.nwk",
        );
        let error = build(&it, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("tree.nwk") => "((a,b),c);",
                _ => GFF,
            }
            .to_string())
        })
        .unwrap_err();
        assert!(
            matches!(error, BuildError::Unjoined { what: "taxon", .. }),
            "{error}"
        );
        // And the message shows the names, which is what makes it actionable.
        assert!(error.to_string().contains("s1"), "{error}");

        // A tree that does have them draws.
        let svg = build(&it, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("tree.nwk") => "((s1,s2),s3);",
                _ => GFF,
            }
            .to_string())
        })
        .unwrap();
        assert!(svg.contains("s1"));

        // The locus side: names from a search against a different FASTA.
        let it = pair("chr1:1-1000", "--loci", "loci.bed", "--links", "links.tsv");
        let error = build(&it, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("links.tsv") => "lcl|g1\tlcl|g2\t98.0\n",
                _ => BED,
            }
            .to_string())
        })
        .unwrap_err();
        assert!(
            matches!(
                error,
                BuildError::Unjoined {
                    what: "gene name",
                    ..
                }
            ),
            "{error}"
        );

        let svg = build(&it, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("links.tsv") => "g1\tg2\t98.0\n",
                _ => BED,
            }
            .to_string())
        })
        .unwrap();
        assert!(svg.contains("homology"), "the joined pair drew no ribbon");
        assert!(
            !svg.contains(", unmatched"),
            "a joined gene was still marked as having nothing to match"
        );
    }

    /// Both new tracks need a second file, and every field of a TrackSpec is
    /// public, so the parser's refusal is not the only one that has to exist.
    #[test]
    fn a_track_drawn_from_two_files_is_refused_here_too_and_not_a_panic() {
        for (flag, kind) in [
            ("--clades", Kind::Clades),
            ("--loci", Kind::Loci),
            ("--tanglegram", Kind::Tanglegram),
        ] {
            let mut it = over("chr1:1-1000", "--tree", "t.nwk");
            it.tracks[0].kind = kind;
            let error = build(&it, |_| Ok("((a,b),c);".to_string())).unwrap_err();
            assert!(
                matches!(error, BuildError::MissingSecond { .. }),
                "{flag}: {error}"
            );
        }
    }

    /// A file that holds what was asked for, somewhere the window is not, is a
    /// different mistake from an empty file and says so. Gubbins writes the
    /// literal SEQUENCE in column one, so this is the ordinary way a correct
    /// file draws nothing.
    #[test]
    fn a_full_file_with_nothing_in_the_window_says_what_it_did_hold() {
        let it = pair(
            "chr1:9000-9500",
            "--clades",
            "clades.gff",
            "--with-tree",
            "tree.nwk",
        );
        let error = build(&it, |source| {
            Ok(match source {
                Source::Path(path) if path.ends_with("tree.nwk") => "((s1,s2),s3);",
                _ => "SEQUENCE\t.\tCDS\t100\t900\t.\t.\t0\ttaxa=\"s1 s2\";\n",
            }
            .to_string())
        })
        .unwrap_err();
        let said = error.to_string();
        assert!(said.contains("though the file holds 1"), "{said}");
        assert!(said.contains("SEQUENCE"), "{said}");
    }

    /// A track drawn from two files, with both of them named.
    fn pair(locus: &str, flag: &str, path: &str, second: &str, other: &str) -> Invocation {
        let args: Vec<String> = [locus, flag, path, second, other]
            .iter()
            .map(|word| word.to_string())
            .collect();
        match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    /// One region, one track flag and one path, which may hold spaces.
    /// A whole command line, for the tests that need more than a track and a
    /// file.
    fn sheeted(line: &str) -> Invocation {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    fn over(locus: &str, flag: &str, path: &str) -> Invocation {
        let args = vec![locus.to_string(), flag.to_string(), path.to_string()];
        match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    /// A name that picks nothing, or picks two things, is refused rather than
    /// passed on. Both builders take a `usize` and neither complains about one
    /// out of range: `MsaTrack::compare_to` falls back to the consensus and
    /// `SnpTrack::from_alignment` returns an empty track, and both of those
    /// are figures. A misspelt name would be answered with a picture.
    #[test]
    fn a_row_to_compare_against_is_refused_where_the_name_picks_wrong() {
        const ALIGNED: &str = "\
>H37Rv
ACGTACGTACGTACGTACGTACGTACGTACGT
>iso1
ACGTACGTACGTACGTACGTACGTACGTAAGT
>iso2
ACGTACGTAAGTACGTACGTACGTACGTACGT
";
        const TWICE: &str = "\
>H37Rv
ACGTACGTACGTACGTACGTACGTACGTACGT
>iso1
ACGTACGTACGTACGTACGTACGTACGTAAGT
>H37Rv
ACGTACGTAAGTACGTACGTACGTACGTACGT
";
        let line = |extra: &[&str]| {
            let mut args = vec![
                "aln:1-32".to_string(),
                "--msa".to_string(),
                "a.fa".to_string(),
            ];
            args.extend(extra.iter().map(|word| word.to_string()));
            match parse(&args).unwrap() {
                Request::Draw(invocation) => *invocation,
                other => panic!("expected a figure, got {other:?}"),
            }
        };

        // A name that is there draws.
        assert!(build(
            &line(&["--compare-to", "iso1"]),
            |_| Ok(ALIGNED.to_string())
        )
        .is_ok());

        let missing = build(
            &line(&["--compare-to", "iso9"]),
            |_| Ok(ALIGNED.to_string()),
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            missing,
            "--msa a.fa has no row called iso9; it has H37Rv, iso1, iso2"
        );

        let twice = build(&line(&["--compare-to", "H37Rv"]), |_| Ok(TWICE.to_string()))
            .unwrap_err()
            .to_string();
        assert_eq!(
            twice,
            "--msa a.fa has 2 rows called H37Rv, so the name does not pick one"
        );

        // And the same on a variable-site panel, which fails differently
        // without the check: an out of range index gives an empty track.
        let mut args: Vec<String> = "aln:1-32 --snps a.fa --compare-to iso9"
            .split_whitespace()
            .map(String::from)
            .collect();
        args.truncate(5);
        let panel = match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        };
        let refused = build(&panel, |_| Ok(ALIGNED.to_string()))
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with("--snps a.fa has no row called iso9"),
            "{refused}"
        );
    }

    /// A panel of variable sites lays out its own columns, so the ruler the
    /// command line appended numbered a region the panel had thrown most of
    /// away: over ten columns of which two vary, the ticks 1 to 5 stood under
    /// the fifth column and 6 to 10 under the ninth. Alone the panel gets no
    /// ruler, and `--no-axis` has nothing left to do there. Stacked with a
    /// track on the coordinates the ruler comes back, because that track is
    /// read against it.
    #[test]
    fn a_panel_of_variable_sites_gets_no_ruler_of_its_own() {
        use crate::{Figure, MsaSequence};

        const ALIGNED: &str = ">ref\nACGTACGTAC\n>s1\nACGTTCGTAC\n>s2\nACGTTCGTGC\n";
        let rows = vec![
            MsaSequence::new("ref", b"ACGTACGTAC".to_vec()),
            MsaSequence::new("s1", b"ACGTTCGTAC".to_vec()),
            MsaSequence::new("s2", b"ACGTTCGTGC".to_vec()),
        ];
        let region = Region::parse("aln:1-10").unwrap();

        let from_cli = build(&over("aln:1-10", "--snps", "a.fa"), |_| {
            Ok(ALIGNED.to_string())
        })
        .unwrap();
        let bare = Figure::new(region.clone())
            .push(SnpTrack::from_alignment(0, &rows).label("a"))
            .to_svg();
        let ruled = Figure::new(region)
            .push(SnpTrack::from_alignment(0, &rows).label("a"))
            .push(crate::AxisTrack::new())
            .to_svg();
        assert_ne!(bare, ruled, "the ruler cannot be told apart here");
        assert_eq!(from_cli, bare, "a ruler went under the panel");

        // Stacked under depth, the ruler is the depth's, and it stays.
        let stacked = |extra: &str| {
            let line = format!("aln:1-10 --coverage d.bg --snps a.fa {extra}");
            build(&invocation(&line), |source| {
                Ok(match source {
                    Source::Path(path) if path.ends_with("a.fa") => ALIGNED.to_string(),
                    _ => "aln\t0\t10\t30\n".to_string(),
                })
            })
            .unwrap()
        };
        assert_ne!(
            stacked(""),
            stacked("--no-axis"),
            "the depth lost the ruler it is read against"
        );
    }

    #[test]
    fn a_reference_is_cut_down_to_the_region() {
        let region = Region::parse("chr1:11-20").unwrap();
        let whole = Reference::new("sequence", "r.fa", "chr1".into(), (b'A'..=b'Z').collect());
        assert_eq!(whole.clip(&region).unwrap(), (10, b"KLMNOPQRST".to_vec()));
    }

    #[test]
    fn a_reference_shorter_than_the_region_is_not_an_index_panic() {
        let short = Reference::new("sequence", "r.fa", "chr1".into(), b"ACGT".to_vec());
        let region = Region::parse("chr1:1-1000").unwrap();
        assert_eq!(short.clip(&region).unwrap(), (0, b"ACGT".to_vec()));
        // And one that holds no base of the window says so rather than
        // drawing an empty band.
        let past = Region::parse("chr1:100-200").unwrap();
        let error = short.clip(&past).unwrap_err().to_string();
        assert_eq!(
            error,
            "--sequence r.fa: chr1 holds bases 1 to 4, and none of them is in chr1:100-200"
        );
    }

    /// `samtools faidx ref.fa pX:101-160` writes sixty bases under the header
    /// `>pX:101-160`. They were drawn from base 1, so a window on exactly the
    /// bases the file held came out as an empty band.
    #[test]
    fn a_slice_samtools_cut_out_is_drawn_where_its_header_says() {
        let plasmid: String = "ACGT".repeat(40);
        let slice = format!(">pX:101-160\n{}\n", &plasmid[100..160]);
        let draw = |locus: &str| build(&over(locus, "--sequence", "s.fa"), |_| Ok(slice.clone()));
        let letters = |svg: &str| -> String {
            svg.split("<text")
                .skip(1)
                .filter_map(|piece| piece.split('>').nth(1)?.split('<').next())
                .filter(|text| matches!(*text, "A" | "C" | "G" | "T"))
                .collect()
        };
        assert_eq!(letters(&draw("pX:101-160").unwrap()), plasmid[100..160]);
        // Part of the window is in the slice, and that part is drawn where it
        // sits on the plasmid.
        assert_eq!(letters(&draw("pX:151-170").unwrap()), plasmid[150..160]);
        // None of it is, and the refusal says what the slice holds.
        let error = draw("pX:1-60").unwrap_err().to_string();
        assert!(
            error.contains("pX:101-160 holds bases 101 to 160, and none of them is in pX:1-60"),
            "{error}"
        );

        // A header that only looks like a span keeps its name and its bases
        // where they were: the span is not as long as the record.
        let named = Reference::new("sequence", "r.fa", "chr1:1-10".into(), b"ACGT".to_vec());
        assert!(!named.is_slice());
        assert_eq!(named.offset, 0);
    }

    #[test]
    fn several_slices_of_one_sequence_are_picked_by_the_window() {
        let fasta = ">chr1:1-4\nACGT\n>chr1:11-14\nTTTT\n";
        let svg = build(&over("chr1:11-14", "--sequence", "s.fa"), |_| {
            Ok(fasta.to_string())
        })
        .unwrap();
        assert_eq!(svg.matches(">T</text>").count(), 4, "{svg}");
        // Two whole records sharing a name are still refused.
        let twice = ">chr1\nACGT\n>chr1\nTTTT\n";
        let error = build(&over("chr1:1-4", "--sequence", "s.fa"), |_| {
            Ok(twice.to_string())
        })
        .unwrap_err();
        assert!(
            error.to_string().contains("2 records called chr1"),
            "{error}"
        );
    }

    /// Every flag that reads a reference refuses one with nothing in the
    /// window, rather than drawing a track with nothing to stand on: a pileup
    /// whose reference is elsewhere drew every read as agreeing with it.
    #[test]
    fn every_track_that_reads_a_reference_refuses_one_outside_the_window() {
        let fasta = ">chr1\nACGTACGTAC\n";
        for line in [
            "chr1:101-160 --sequence r.fa",
            "chr1:101-160 --orfs r.fa",
            "chr1:101-160 --pileup reads.sam --with-sequence r.fa",
            "chr1:101-160 --dynseq scores.bg --with-sequence r.fa",
        ] {
            let args: Vec<String> = line.split_whitespace().map(String::from).collect();
            let crate::cli::args::Request::Draw(invocation) =
                crate::cli::args::parse(&args).unwrap()
            else {
                unreachable!("every line draws")
            };
            let error = build(&invocation, |source| {
                let Source::Path(path) = source else {
                    unreachable!("every source is a file")
                };
                let path = path.to_string_lossy();
                Ok(if path.ends_with(".fa") {
                    fasta.to_string()
                } else if path.ends_with(".sam") {
                    "r1\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\tIIII\n".to_string()
                } else {
                    "chr1\t100\t104\t0.5\n".to_string()
                })
            })
            .unwrap_err()
            .to_string();
            assert!(
                error.contains("chr1 holds bases 1 to 10, and none of them is in chr1:101-160"),
                "{line}: {error}"
            );
        }
    }

    #[test]
    fn an_empty_stack_is_still_a_figure_with_a_ruler() {
        let svg = build(&invocation("chr1:1-1000"), open_from_disk).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.ends_with("</svg>"));
    }

    #[test]
    fn figure_flags_reach_the_figure() {
        let svg = build(
            &invocation("chr1:1-1000 --title one --width 1200"),
            open_from_disk,
        )
        .unwrap();
        assert!(svg.contains("width=\"1200\""));
        assert!(svg.contains(">one<"));
    }

    #[test]
    fn a_file_that_is_not_there_says_which_flag_wanted_it() {
        let error = build(
            &invocation("chr1:1-1000 --features nowhere.bed"),
            open_from_disk,
        )
        .unwrap_err();
        let message = error.to_string();
        assert!(message.starts_with("--features nowhere.bed:"), "{message}");
    }

    #[test]
    fn a_matrix_whose_sites_all_lie_outside_the_region_is_an_error() {
        // Two samples typed at sites 9000 and 9100, drawn over chr1:1-100. The
        // rows survive the read and every column is filtered out, which used to
        // give a figure of sample names beside no cells and exit 0.
        let path = written(
            "outside.matrix.tsv",
            "sample\t9000\t9100\nERR1\t1\t0\nERR2\t0\t1\n",
        );
        let error = build(&over("chr1:1-100", "--matrix", &path), open_from_disk).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("--matrix {path}: no samples in the region")
        );
    }

    #[test]
    fn a_matrix_with_a_site_in_the_region_still_draws() {
        let path = written("inside.matrix.tsv", "sample\t50\t60\nERR1\t1\t0\n");
        let svg = build(&over("chr1:1-100", "--matrix", &path), open_from_disk).unwrap();
        assert!(svg.contains("ERR1"), "{svg}");
    }

    /// The enumerating pass over a whole file is skipped where the reader has
    /// already named what they want, so the four answers it used to give have
    /// to keep coming from somewhere.
    ///
    /// This pins those four answers and not the skip itself: the skip is a
    /// saving and not a behaviour, so putting the pass back makes this test
    /// pass just the same. What it guards is the risk the saving carried, which
    /// is that one of the four refusals had been leaning on the pass that no
    /// longer runs.
    #[test]
    fn naming_the_modification_skips_the_pass_that_lists_them_and_keeps_every_refusal() {
        let two = concat!(
            "chr1\t1\t2\tm\t10\t+\t1\t2\t0,0,0\t10\t50.00\t5\t5\t0\t0\t0\t0\t0\n",
            "chr1\t3\t4\th\t10\t+\t3\t4\t0,0,0\t10\t50.00\t5\t5\t0\t0\t0\t0\t0\n"
        );
        let path = written("two-codes.bed", two);

        // Named and present: it draws, without ever listing the codes.
        let mut named = over("chr1:1-100", "--methylation", &path);
        named.tracks[0].selects = Some("m".to_string());
        assert!(build(&named, open_from_disk).is_ok());

        // Named and absent: the refusal comes from the pass that reads.
        let mut absent = over("chr1:1-100", "--methylation", &path);
        absent.tracks[0].selects = Some("a".to_string());
        let error = build(&absent, open_from_disk).unwrap_err().to_string();
        assert!(error.contains("no modified bases"), "{error}");

        // Not named: the file holds two, and the reader is asked which.
        let error = build(&over("chr1:1-100", "--methylation", &path), open_from_disk)
            .unwrap_err()
            .to_string();
        assert!(error.contains("--modification"), "{error}");
        assert!(error.contains('h') && error.contains('m'), "{error}");

        // Nothing in it at all, which the reading pass has to notice on its own
        // now that nothing counted the codes first.
        let empty = written("no-codes.bed", "");
        let mut named_empty = over("chr1:1-100", "--methylation", &empty);
        named_empty.tracks[0].selects = Some("m".to_string());
        let error = build(&named_empty, open_from_disk).unwrap_err().to_string();
        assert!(error.contains("no modified bases"), "{error}");
    }

    /// A window whose every position went unmeasured was refused as holding
    /// no modified bases while the file held some, which reads as though they
    /// were elsewhere. They were here, and the band has a way to say so: it
    /// counts them in its corner, as it counts the calls under the floor, and
    /// a floor that hides every call already draws the band with that count.
    #[test]
    fn a_window_where_nothing_was_measured_is_drawn_with_the_count() {
        let unmeasured = concat!(
            "chr1\t1\t2\tm\t0\t+\t1\t2\t0,0,0\t0\t0.00\t0\t0\t0\t0\t0\t0\t0\n",
            "chr1\t3\t4\tm\t0\t+\t3\t4\t0,0,0\t0\t0.00\t0\t0\t0\t0\t0\t0\t0\n",
            "chr1\t500\t501\tm\t10\t+\t500\t501\t0,0,0\t10\t50.00\t5\t5\t0\t0\t0\t0\t0\n"
        );
        let path = written("unmeasured.bed", unmeasured);
        let svg = build(&over("chr1:1-100", "--methylation", &path), open_from_disk).unwrap();
        assert!(svg.contains(">2 with no coverage</text>"), "{svg}");

        // Nothing in the window at all is still refused, with the file's rows.
        let error = build(
            &over("chr1:200-300", "--methylation", &path),
            open_from_disk,
        )
        .unwrap_err()
        .to_string();
        assert!(
            error.contains("no modified bases in chr1:200-300"),
            "{error}"
        );
    }

    #[test]
    fn a_tree_that_will_not_parse_names_the_file_it_was() {
        // The other failures of one flag name the file, and this one did not.
        let path = written("unbalanced.nwk", "((a,b\n");
        let error = build(&over("chr1:1-100", "--tree", &path), open_from_disk).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("--tree {path}: invalid Newick tree: unbalanced parentheses")
        );
    }
    /// A page that runs the program on every move reads the same file every
    /// time, and reading it is most of the work. `build_with` lets such a
    /// caller answer with a tree it made earlier, and this is the contract:
    /// the text is offered before it is read, an answer is taken as given, and
    /// no answer means the file is read as usual.
    #[test]
    fn a_caller_may_answer_with_a_tree_it_read_earlier() {
        const ON_DISK: &str = "((ondisk_a:0.1,ondisk_b:0.1):0.1,ondisk_c:0.1);";
        let held = Tree::parse_newick("((kept_a:0.1,kept_b:0.1):0.1,kept_c:0.1);").unwrap();

        // Offered, and taken.
        let mut seen = Vec::new();
        let svg = build_with(
            &over("tree:1-1", "--tree", "t.nwk"),
            |_| Ok(ON_DISK.to_string()),
            |name, text| {
                seen.push((name.to_string(), text.to_string()));
                Some(held.clone())
            },
        )
        .unwrap();
        assert_eq!(seen.len(), 1, "offered once, for the one tree");
        assert_eq!(seen[0].0, "t.nwk", "offered under the name it was asked by");
        assert_eq!(seen[0].1, ON_DISK, "offered the text, trimmed");
        assert!(svg.contains("kept_a"), "the answer was drawn: {svg}");
        assert!(!svg.contains("ondisk_a"), "and the file was not: {svg}");

        // Declined, and the file is read.
        let svg = build_with(
            &over("tree:1-1", "--tree", "t.nwk"),
            |_| Ok(ON_DISK.to_string()),
            |_, _| None,
        )
        .unwrap();
        assert!(svg.contains("ondisk_a"), "{svg}");

        // And `build` is the same thing with nothing held, so a shell is not
        // paying for a viewer's convenience.
        let plain = build(&over("tree:1-1", "--tree", "t.nwk"), |_| {
            Ok(ON_DISK.to_string())
        })
        .unwrap();
        assert_eq!(plain, svg);
    }
    /// A sheet reaches a tree by a different road from every other track. The
    /// others hand the sheet to the track and the track draws out of it; a tree
    /// reads its strips out of its own annotations, so the sheet is copied onto
    /// the tips it names. This is the check that the road arrives.
    #[test]
    fn a_sheet_of_metadata_becomes_strips_beside_a_tree() {
        const TREE: &str = "((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);";
        const SHEET: &str = "name\tlineage\tdepth\na\tL4\t30\nb\tL4\t50\nc\tL2\t70\nd\tL1\t90\n";
        let open = |source: &Source| -> io::Result<String> {
            Ok(match source {
                Source::Path(path) if path.to_string_lossy().ends_with(".nwk") => TREE.to_string(),
                _ => SHEET.to_string(),
            })
        };

        let svg = build(&sheeted("tree:1-1 --tree t.nwk --traits s.tsv"), open).unwrap();
        // Every tip gets its own value, and the heading is written out rather
        // than cut to one letter, which is what a fourteen pixel column did.
        assert!(svg.contains("a; lineage L4"), "{svg}");
        assert!(svg.contains("d; depth 90"), "{svg}");
        assert!(
            svg.contains(">lineage</text>"),
            "the heading is legible: {svg}"
        );
    }

    /// The colour drawn under each tooltip that starts with `title`, in order.
    fn painted(svg: &str, title: &str) -> Vec<String> {
        svg.split("<title>")
            .skip(1)
            .filter(|piece| piece.starts_with(title))
            .filter_map(|piece| {
                let after = piece.split("</title>").nth(1)?;
                ["fill=\"#", "stroke=\"#"]
                    .iter()
                    .filter_map(|attribute| {
                        after.find(attribute).map(|at| at + attribute.len() - 1)
                    })
                    .min()
                    .map(|at| after[at..at + 7].to_string())
            })
            .collect()
    }

    /// One sheet, a tree and a matrix under it: a lineage is one colour in
    /// both strips. The tree dealt the palette in the order it met its tips and
    /// the matrix in the order its sheet sorted the names, so each of the three
    /// lineages here came out in a different colour on each side, in a figure
    /// whose whole point is to read one strip against the other.
    #[test]
    fn a_sheet_shared_by_a_tree_and_a_matrix_colours_each_level_once() {
        const TREE: &str = "((Zed:1,Yan:1):1,(Abe:1,Bo:1):1);";
        // Three orders that all disagree: the file lists L1 first, the names
        // sort L2 first (Abe), and the tree meets L4 first (Zed).
        const SHEET: &str = "sample\tlineage\nBo\tL1\nZed\tL4\nAbe\tL2\nYan\tL4\n";
        const MATRIX: &str = "sample\t100\t200\nZed\t1\t0\nYan\t0\t1\nAbe\t1\t1\nBo\t0\t0\n";
        let open = |source: &Source| -> io::Result<String> {
            let Source::Path(path) = source else {
                unreachable!("every source here is a file")
            };
            let path = path.to_string_lossy();
            Ok(if path.ends_with(".nwk") {
                TREE
            } else if path.ends_with("m.tsv") {
                MATRIX
            } else {
                SHEET
            }
            .to_string())
        };
        let svg = build(
            &sheeted("chr:1-300 --tree t.nwk --traits s.tsv --matrix m.tsv --traits s.tsv"),
            open,
        )
        .unwrap();

        let colours = [
            ("L2", painted(&svg, "Abe; lineage L2")),
            ("L1", painted(&svg, "Bo; lineage L1")),
            ("L4", painted(&svg, "Zed; lineage L4")),
        ];
        for (level, both) in &colours {
            assert_eq!(
                both.len(),
                2,
                "{level} is drawn beside the tree and the matrix"
            );
            assert_eq!(both[0], both[1], "{level} is two colours: {both:?}");
        }
        // And the palette is dealt in the order the file lists the levels,
        // which is neither order the two tracks met them in.
        let theme = crate::Theme::light();
        assert_eq!(
            colours[1].1[0],
            theme.color(0),
            "L1 comes first in the file"
        );
        assert_eq!(colours[2].1[0], theme.color(1));
        assert_eq!(colours[0].1[0], theme.color(2));
    }

    /// A tree drawn with no region prints none, and draws no ruler: the
    /// window the figure lays its width over is one nothing is drawn in.
    #[test]
    fn a_tree_drawn_with_no_region_prints_no_window() {
        const TREE: &str = "((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);";
        let args: Vec<String> = "--tree t.nwk --label phylogeny"
            .split_whitespace()
            .map(String::from)
            .collect();
        let crate::cli::args::Request::Draw(invocation) = crate::cli::args::parse(&args).unwrap()
        else {
            unreachable!("a tree is drawn")
        };
        let svg = build(&invocation, |_| Ok(TREE.to_string())).unwrap();
        assert!(!svg.contains("phylogeny:1-1"), "{svg}");
        assert!(!svg.contains("1-1"), "a window was printed: {svg}");
        assert!(svg.contains("<title id=\"karyon-title\">phylogeny</title>"));
        // Every tip is drawn, and nothing else is written as text but the
        // length of the scale bar, which measures branches: no locus and no
        // tick of a ruler.
        let text: Vec<&str> = svg
            .split("<text")
            .skip(1)
            .filter_map(|piece| piece.split('>').nth(1)?.split('<').next())
            .collect();
        assert_eq!(text, ["phylogeny", "a", "b", "c", "d", "0.02"], "{svg}");
    }

    /// A scan as association tools write it, a column of p-values under the
    /// header P. It was drawn as written, so the hit at 4e-12 sat on the floor
    /// of the figure and the null at 0.9 at the top, and the command exited
    /// nought.
    #[test]
    fn a_scan_of_p_values_is_drawn_with_its_strongest_hit_highest() {
        const SCAN: &str = "CHR\tBP\tP\nchr1\t100\t0.5\nchr1\t200\t4e-12\nchr1\t300\t0.9\n";
        let draw = |line: &str| {
            let args: Vec<String> = line.split_whitespace().map(String::from).collect();
            let crate::cli::args::Request::Draw(invocation) =
                crate::cli::args::parse(&args).unwrap()
            else {
                unreachable!("a scan is drawn")
            };
            build(&invocation, |_| Ok(SCAN.to_string()))
        };
        let svg = draw("chr1:1-400 --manhattan gwas.tsv").unwrap();
        // Points in the order of their positions, left to right.
        let mut points: Vec<(f64, f64)> = svg
            .split("<circle cx=\"")
            .skip(1)
            .filter_map(|piece| {
                let (cx, rest) = piece.split_once('"')?;
                let cy = rest.split_once("cy=\"")?.1.split_once('"')?.0;
                Some((cx.parse().ok()?, cy.parse().ok()?))
            })
            .collect();
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert_eq!(points.len(), 3, "{svg}");
        let heights: Vec<f64> = points.iter().map(|(_, cy)| *cy).collect();
        assert!(
            heights[1] < heights[0] && heights[0] < heights[2],
            "the hit at 4e-12 is not the highest: {heights:?}"
        );
        assert!(
            svg.contains("-log10 p"),
            "the axis does not say what it shows"
        );

        // The threshold is in the file's units: a p-value, drawn at -log10 of
        // itself, which for 5e-8 is the convention asked for by name.
        let named = draw("chr1:1-400 --manhattan gwas.tsv --threshold genome-wide").unwrap();
        assert_eq!(
            draw("chr1:1-400 --manhattan gwas.tsv --threshold 5e-8").unwrap(),
            named
        );
        // The line says the p-value it is at, and the axis's title is a line
        // of its own rather than words after the top tick.
        assert!(named.contains(">p = 5e-8</text>"), "{named}");
        assert!(named.contains(">-log10 p</text>"), "{named}");
        let own = draw("chr1:1-400 --manhattan gwas.tsv --threshold 1e-5").unwrap();
        assert!(own.contains(">p = 1e-5</text>"), "{own}");
        // The number that was right for a file of -log10 values is no p-value.
        let error = draw("chr1:1-400 --manhattan gwas.tsv --threshold 7.3").unwrap_err();
        assert!(
            error
                .to_string()
                .contains("holds p-values, so --threshold is a p-value too"),
            "{error}"
        );
    }

    #[test]
    fn a_file_that_is_not_text_is_known_by_its_first_bytes() {
        let path = |name: &str| Some(Path::new(name).to_path_buf());
        let of = |bytes: &[u8], name: &str| Binary::of(bytes, path(name).as_deref());
        let gzip = [0x1f, 0x8b, 0x08, 0x04];
        assert_eq!(of(&gzip, "calls.vcf.gz"), Some(Binary::Gzip));
        assert_eq!(of(&gzip, "reads.BAM"), Some(Binary::Bam));
        assert_eq!(of(&gzip, "calls.bcf"), Some(Binary::Bcf));
        // A BCF written bare, as gzip leaves what `bcftools view -Ou` writes,
        // whatever it is called.
        assert_eq!(of(b"BCF  ", "calls.data"), Some(Binary::Bcf));
        assert_eq!(of(b"CRAM\x03\x01", "reads.cram"), Some(Binary::Cram));
        assert_eq!(of(b"BZh91AY", "a.bz2"), Some(Binary::Bzip2));
        assert_eq!(
            of(&[0xfd, b'7', b'z', b'X', b'Z', 0], "a.xz"),
            Some(Binary::Xz)
        );
        assert_eq!(of(&[0x28, 0xb5, 0x2f, 0xfd], "a.zst"), Some(Binary::Zstd));
        // The UCSC formats in either byte order.
        assert_eq!(of(&[0x26, 0xfc, 0x8f, 0x88], "d.bw"), Some(Binary::BigWig));
        assert_eq!(of(&[0x88, 0x8f, 0xfc, 0x26], "d.bw"), Some(Binary::BigWig));
        assert_eq!(of(&[0xeb, 0xf2, 0x89, 0x87], "f.bb"), Some(Binary::BigBed));
        assert_eq!(
            of(&[0x43, 0x27, 0x41, 0x1a], "g.2bit"),
            Some(Binary::TwoBit)
        );
        // cooler's HDF5, one resolution or several by the name, and Juicer's.
        let hdf5 = b"\x89HDF\r\n\x1a\n\0\0";
        assert_eq!(of(hdf5, "contacts.cool"), Some(Binary::Cool));
        assert_eq!(of(hdf5, "contacts.mcool"), Some(Binary::Mcool));
        assert_eq!(of(b"HIC\0\x08\0\0\0", "contacts.hic"), Some(Binary::Hic));
        // Text that is not UTF-8 is not a format.
        assert_eq!(of(&[b'c', b'h', b'r', 0xff], "a.bed"), None);
    }

    /// "stream did not contain valid UTF-8" was the whole of what a first try
    /// with a compressed VCF or a BAM was told.
    ///
    /// The table holds the bare command, and the message is checked for it as
    /// the system the test runs on should word it: `<(...)` on Linux and
    /// macOS, a pipe into `-` on Windows. `cfg!` rather than `HOST`, so a
    /// `HOST` that picked the wrong wording fails here as well.
    #[test]
    fn a_file_that_is_not_text_is_answered_with_the_command_that_reads_it() {
        let line = |text: &str| -> Invocation {
            let args: Vec<String> = text.split_whitespace().map(String::from).collect();
            match parse(&args).unwrap() {
                Request::Draw(invocation) => *invocation,
                other => panic!("expected a figure, got {other:?}"),
            }
        };
        let refused = |text: &str, binary: Binary| {
            build(&line(text), |_| {
                Err(io::Error::new(io::ErrorKind::InvalidData, binary))
            })
            .unwrap_err()
            .to_string()
        };
        // The flag is the one a file named on its own has to be given in
        // front of what replaces it, and `None` where the line gave it.
        for (text, binary, command, flag) in [
            (
                "chr1:1-5000 --variants calls.vcf.gz",
                Binary::Gzip,
                "gzip -dc calls.vcf.gz",
                None,
            ),
            (
                "chr1:1-5000 --pileup reads.bam",
                Binary::Bam,
                "samtools view -h reads.bam chr1:1-5000",
                None,
            ),
            (
                "chr1:1-5000 --coverage reads.bam",
                Binary::Bam,
                "samtools depth -a -r chr1:1-5000 reads.bam",
                None,
            ),
            (
                "chr1:1-5000 --variants calls.bcf",
                Binary::Bcf,
                "bcftools view calls.bcf",
                None,
            ),
            (
                "chr1:1-5000 calls.bcf",
                Binary::Bcf,
                "bcftools view calls.bcf",
                Some("variants"),
            ),
            (
                "chr1:1-5000 --coverage depth.bw",
                Binary::BigWig,
                "bigWigToBedGraph -chrom=chr1 -start=0 -end=5000 depth.bw /dev/stdout",
                None,
            ),
            (
                "chr1:1-5000 depth.bw",
                Binary::BigWig,
                "bigWigToBedGraph -chrom=chr1 -start=0 -end=5000 depth.bw /dev/stdout",
                Some("coverage"),
            ),
            (
                "chr1:1-5000 contacts.cool",
                Binary::Cool,
                "cooler dump --join -r chr1:1-5000 contacts.cool",
                Some("pairs"),
            ),
        ] {
            let error = refused(text, binary);
            assert!(
                error.contains(&in_place(HOST, command, flag)),
                "{text}: {error}"
            );
            // Named on its own, the file's replacement carries its flag, and
            // given with one, it needs none.
            assert_eq!(
                error.contains("; write --") || error.contains(", with --"),
                flag.is_some(),
                "{text}: {error}"
            );
            assert_eq!(
                error.contains(&format!("<({command})")),
                !cfg!(windows),
                "{text}: {error}"
            );
            assert!(error.contains("karyon reads text"), "{error}");
        }
        // The messages the guide prints, to the letter, each where it is
        // printed. A BCF is read as it is where its files give its bytes,
        // and is answered so only where they give text alone.
        let cram = refused("chr1:1-5000 --pileup aln.cram", Binary::Cram);
        if cfg!(windows) {
            assert_eq!(
                cram,
                "--pileup aln.cram: the file is CRAM, and karyon reads text; pipe what \
                 samtools view -h aln.cram chr1:1-5000 writes into karyon, with - where its \
                 name is, or turn it into text first"
            );
        } else {
            assert_eq!(
                cram,
                "--pileup aln.cram: the file is CRAM, and karyon reads text; write \
                 <(samtools view -h aln.cram chr1:1-5000) where its name is, or turn it into \
                 text first"
            );
        }
        let alone = refused("chr1:1-5000 aln.cram", Binary::Cram);
        assert!(
            alone.contains(&in_place(
                HOST,
                "samtools depth -a -r chr1:1-5000 aln.cram",
                Some("coverage")
            )),
            "{alone}"
        );
        // No one command: the resolution has to be picked first. A .hic is
        // read as it is, and only a pipe of one is refused, asking its name.
        for (text, binary, advice) in [
            (
                "chr1:1-5000 contacts.mcool",
                Binary::Mcool,
                "cooler ls lists its resolutions",
            ),
            (
                "chr1:1-5000 --pairs -",
                Binary::Hic,
                "a pipe cannot be read that way; name the file instead",
            ),
        ] {
            let error = refused(text, binary);
            assert!(error.contains(advice), "{text}: {error}");
            assert!(!error.contains("<("), "{text}: {error}");
        }
        // A pipe has no name to write a command in place of.
        let error = refused("chr1:1-5000 --pileup -", Binary::Bam);
        assert!(error.contains("pipe it through the tool"), "{error}");
    }

    /// A Windows build offers a pipe into `-`, which every Windows shell has,
    /// and never `<(...)`, which hands over a path no Windows program can
    /// open. Both wordings are checked on every system, since `cfg!` compiles
    /// both, so a run on Linux catches the Windows one going wrong as well.
    #[test]
    fn a_windows_build_pipes_the_command_in_rather_than_naming_it() {
        let command = "samtools view -h aln.cram chr1:1-5000";
        let windows = in_place(Shell::Windows, command, None);
        assert_eq!(
            windows,
            "pipe what samtools view -h aln.cram chr1:1-5000 writes into karyon, with - where \
             its name is"
        );
        assert!(!windows.contains("<("), "{windows}");
        let posix = in_place(Shell::Posix, command, None);
        assert_eq!(
            posix,
            "write <(samtools view -h aln.cram chr1:1-5000) where its name is"
        );
        // A file named on its own, as the guide prints it for each.
        let depth = "samtools depth -a -r chr1:1-5000 aln.cram";
        assert_eq!(
            in_place(Shell::Windows, depth, Some("coverage")),
            "pipe what samtools depth -a -r chr1:1-5000 aln.cram writes into karyon, with \
             --coverage - where its name is"
        );
        assert_eq!(
            in_place(Shell::Posix, depth, Some("coverage")),
            "write --coverage <(samtools depth -a -r chr1:1-5000 aln.cram) where its name is"
        );
        assert_eq!(HOST == Shell::Windows, cfg!(windows));
    }

    /// The advice, done as it says, is a command line karyon draws from.
    ///
    /// A file named on its own was given its track by its name, and neither
    /// `-` nor the `/dev/fd/63` a shell hands over for `<(...)` has a name to
    /// give one: told to put a bare `-` where `calls.bcf` was, a Windows user
    /// was refused for want of a track, and a bare `<(...)` was looked for as
    /// a gene or a sequence called `/dev/fd/63`. Each line here is refused,
    /// the message's words put where the file's name was, and the line parsed
    /// again; it has to come back as the same track reading the command's
    /// output.
    #[test]
    fn the_command_a_file_that_is_not_text_is_answered_with_runs_as_written() {
        let words =
            |text: &str| -> Vec<String> { text.split_whitespace().map(String::from).collect() };
        for (text, path, binary, kind) in [
            (
                "chr1:1-5000 calls.bcf",
                "calls.bcf",
                Binary::Bcf,
                Kind::Variants,
            ),
            (
                "chr1:1-5000 --variants calls.bcf",
                "calls.bcf",
                Binary::Bcf,
                Kind::Variants,
            ),
            (
                "chr1:1-5000 contacts.cool -o map.svg",
                "contacts.cool",
                Binary::Cool,
                Kind::Pairs,
            ),
            (
                "chr1:1-5000 depth.bw",
                "depth.bw",
                Binary::BigWig,
                Kind::Coverage,
            ),
        ] {
            let Request::Draw(invocation) = parse(&words(text)).unwrap() else {
                unreachable!("{text} draws a figure")
            };
            let BuildError::NotText {
                track,
                instead: Some(command),
                alone,
                ..
            } = (match build(&invocation, |_| {
                Err(io::Error::new(io::ErrorKind::InvalidData, binary))
            }) {
                Err(error) => error,
                Ok(_) => panic!("{text} drew from a file that is not text"),
            })
            else {
                panic!("{text} was not answered as a file that is not text")
            };
            let flag = alone.then_some(track);
            for shell in [Shell::Posix, Shell::Windows] {
                let advice = in_place(shell, &command, flag);
                // What goes where the name was, as karyon is handed it: the
                // shell turns `<(...)` into a path to a pipe, and the pipe
                // into a Windows build arrives on standard input.
                let put = match shell {
                    Shell::Posix => advice
                        .strip_prefix("write ")
                        .map(|rest| rest.replace(&format!("<({command})"), "/dev/fd/63")),
                    Shell::Windows => advice
                        .strip_prefix(&format!("pipe what {command} writes into karyon, with "))
                        .map(String::from),
                }
                .and_then(|rest| rest.strip_suffix(" where its name is").map(String::from))
                .unwrap_or_else(|| panic!("{text}: {advice}"));
                let followed = text.replace(path, &put);
                let parsed = match parse(&words(&followed)) {
                    Ok(Request::Draw(invocation)) => invocation,
                    other => panic!("{text}: {advice}: {followed} gave {other:?}"),
                };
                let reads = match shell {
                    Shell::Posix => Source::Path(std::path::PathBuf::from("/dev/fd/63")),
                    Shell::Windows => Source::Stdin,
                };
                assert_eq!(parsed.tracks.len(), 1, "{followed}: {parsed:?}");
                assert_eq!(parsed.tracks[0].kind, kind, "{followed}");
                assert_eq!(parsed.tracks[0].source, Some(reads), "{followed}");
                assert_eq!(parsed.region, invocation.region, "{followed}");
            }
        }
    }

    /// A tanglegram names its trees after their files, less the folders,
    /// whichever separator the system writes. The paths are built with
    /// `join`, so on Windows they are `trees\run 7\before.nwk`, which a cut
    /// at `/` alone printed whole over the tree.
    #[test]
    fn a_tanglegram_names_its_trees_after_the_files_and_not_their_folders() {
        let folder = Path::new("trees").join("run 7");
        let before = folder.join("before.nwk").display().to_string();
        let after = folder.join("after.nwk").display().to_string();
        let mut held = Held::new();
        held.insert(before.as_str(), "((a,b),(c,d));");
        held.insert(after.as_str(), "((a,c),(b,d));");
        let args = vec![
            "--tanglegram".to_string(),
            before,
            "--against".to_string(),
            after,
        ];
        let Request::Draw(invocation) = parse(&args).unwrap() else {
            unreachable!("a figure")
        };
        let svg = build_files(&invocation, &mut held, |_, _| None).unwrap();
        assert!(svg.contains(">before.nwk<"), "{svg}");
        assert!(svg.contains(">after.nwk<"), "{svg}");
        assert!(!svg.contains("run 7"), "a folder was printed over a tree");
        // A `\` is a separator on Windows and a letter of a name elsewhere.
        let typed = r"C:\runs\before.nwk";
        let wanted = if cfg!(windows) { "before.nwk" } else { typed };
        assert_eq!(shortened(typed), wanted);
    }

    /// The usual reason a window holds nothing is a file that names its
    /// sequences one way and a region that names them another.
    #[test]
    fn a_window_that_holds_nothing_says_what_the_file_holds_and_where() {
        const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\t.\t.
chr1\t200\t.\tA\tG\t.\t.\t.
chr2\t300\t.\tA\tG\t.\t.\t.
";
        let error = build(&over("1:1-5000", "--variants", "calls.vcf"), |_| {
            Ok(VCF.to_string())
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            error,
            "--variants calls.vcf: no variants in 1:1-5000, though the file holds 3 on chr1 and chr2; \
             if chr1 is 1, add --rename chr1=1"
        );
        // A header is not a sequence: a bedGraph's column two says which rows
        // are rows.
        let bedgraph = "chrom\tstart\tend\tvalue\n1\t0\t10\t5\n";
        let error = build(&over("chr1:1-5000", "--coverage", "d.bg"), |_| {
            Ok(bedgraph.to_string())
        })
        .unwrap_err()
        .to_string();
        assert!(
            error.ends_with("though the file holds 1 on 1; if 1 is chr1, add --rename 1=chr1"),
            "{error}"
        );
        // Read from a pipe, there is nothing to read twice.
        let error = build(
            &over("1:1-5000", "--variants", "-"),
            |_| Ok(VCF.to_string()),
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            error,
            "--variants standard input: no variants in the region"
        );
    }

    #[test]
    fn a_list_of_sequences_reads_as_a_sentence() {
        let held = |names: &[&str]| -> Vec<(String, usize)> {
            names.iter().map(|name| ((*name).to_string(), 1)).collect()
        };
        assert_eq!(listed(&held(&["chr1"])), "chr1");
        assert_eq!(listed(&held(&["1", "2"])), "1 and 2");
        assert_eq!(listed(&held(&["1", "2", "3"])), "1, 2 and 3");
        assert_eq!(
            listed(&held(&["1", "2", "3", "4", "5", "6", "7"])),
            "1, 2, 3, 4, 5 and 2 more sequences"
        );
    }

    /// A directory of its own for a test that reads files from disk, emptied
    /// when the test is done.
    struct Scratch(std::path::PathBuf);

    impl Scratch {
        fn new(name: &str) -> Scratch {
            let dir = std::env::temp_dir().join(format!("karyon-{name}-{}", std::process::id()));
            fs::create_dir_all(&dir).unwrap();
            Scratch(dir)
        }

        fn write(&self, name: &str, bytes: &[u8]) -> String {
            let path = self.0.join(name);
            fs::write(&path, bytes).unwrap();
            path.display().to_string()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn drawn_from_disk(line: &str) -> Result<String, BuildError> {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        let Request::Draw(invocation) = parse(&args).unwrap() else {
            unreachable!("a figure")
        };
        build_files(&invocation, &mut Disk::default(), |_, _| None)
    }

    /// A BAM was refused, with the samtools command to pipe it through; it is
    /// read as it is now, a window at a time through the index beside it.
    #[test]
    fn a_bam_on_disk_is_drawn_as_depth_and_as_reads() {
        let dir = Scratch::new("bam");
        let bam = dir.write("tiny.bam", &crate::read::bam::fixture::BAM);
        dir.write("tiny.bam.bai", &crate::read::bam::fixture::BAI);
        let svg = drawn_from_disk(&format!("chr1:10-35 --coverage {bam} --pileup {bam}")).unwrap();
        // The pileup draws every mapped record over the window, a to f; g starts
        // after it, and u is unmapped.
        for read in ["10 to 19", "15 to 26", "18 to 25", "30 to 139"] {
            assert!(
                svg.contains(&format!("read, {read}")),
                "no read {read}: {svg}"
            );
        }
        // The depth tops out at three, which samtools depth says it does.
        assert!(svg.contains(">3</text>"), "the depth's top is not 3");

        // And without the index, from the start of the file.
        fs::remove_file(dir.0.join("tiny.bam.bai")).unwrap();
        let scanned =
            drawn_from_disk(&format!("chr1:10-35 --coverage {bam} --pileup {bam}")).unwrap();
        assert_eq!(scanned, svg);
    }

    /// The index beside a BAM is found by the path as the system writes it,
    /// with backslashes on Windows, under either of its two names.
    ///
    /// A figure cannot tell: without its index a BAM is read from its start
    /// and draws the same bytes, so a run that drew the right figure from a
    /// backslashed path says nothing of the index. This asks for it instead,
    /// and an index that is found is read, so one that is not an index is
    /// refused rather than passed over.
    #[test]
    fn the_index_beside_a_bam_is_found_by_the_path_the_system_writes() {
        use crate::read::bam::fixture::{BAI, BAM};
        let dir = Scratch::new("bai");
        let bam = dir.write("tiny.bam", &BAM);
        if cfg!(windows) {
            assert!(bam.contains('\\'), "{bam}");
        }
        let source = Source::Path(std::path::PathBuf::from(&bam));
        let region = Region::parse("chr1:10-35").unwrap();
        let mut disk = Disk::default();
        assert!(disk.beside(&source, ".bai").unwrap().is_none());
        for name in ["tiny.bam.bai", "tiny.bai"] {
            let index = dir.write(name, &BAI);
            assert!(
                disk.beside(&source, ".bai").unwrap().is_some(),
                "{index} was not found beside {bam}"
            );
            fs::remove_file(&index).unwrap();
        }
        dir.write("tiny.bam.bai", b"not an index");
        assert!(bam_window(&mut disk, &source, &region).is_err());
    }

    /// Files held in memory draw what the same files on disk draw: a BAM a
    /// window at a time through the index held beside it, or from its start
    /// without one, a sequence as long as its header says, one read by its
    /// name, and compressed text out of its wrapper.
    #[test]
    fn files_held_in_memory_draw_what_the_same_files_on_disk_draw() {
        use crate::read::bam::fixture::{BAI, BAM};
        let dir = Scratch::new("held");
        let bam = dir.write("tiny.bam", &BAM);
        dir.write("tiny.bam.bai", &BAI);
        let mut held = Held::new();
        held.insert(bam.as_str(), BAM);
        held.insert(format!("{bam}.bai"), BAI);
        let from_memory = |line: &str, held: &mut Held| {
            build_files(&invocation(line), held, |_, _| None).unwrap()
        };
        let window = format!("chr1:10-35 --coverage {bam} --pileup {bam}");
        let drawn = from_memory(&window, &mut held);
        assert_eq!(drawn, drawn_from_disk(&window).unwrap());
        let mut unindexed = Held::new();
        unindexed.insert(bam.as_str(), BAM);
        assert_eq!(from_memory(&window, &mut unindexed), drawn);
        let whole = format!("chr1 {bam}");
        assert_eq!(
            from_memory(&whole, &mut held),
            drawn_from_disk(&whole).unwrap()
        );
        let source = Source::Path(bam.clone().into());
        let read = held.named_read(&source, "a").unwrap();
        assert!(
            read.as_deref().is_some_and(|sam| sam.contains("\na\t")),
            "{read:?}"
        );
        assert_eq!(read, Disk::default().named_read(&source, "a").unwrap());

        let calls = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/data/calls.vcf.gz");
        let calls = calls.display().to_string();
        let mut held = Held::new();
        held.insert(calls.as_str(), fs::read(&calls).unwrap());
        let line = format!("NC_000962.3:759,000-764,000 {calls}");
        let drawn = from_memory(&line, &mut held);
        assert_eq!(drawn, drawn_from_disk(&line).unwrap());
        // And with the index beside them, read through it, as from disk.
        held.insert(
            format!("{calls}.tbi"),
            fs::read(format!("{calls}.tbi")).unwrap(),
        );
        assert_eq!(from_memory(&line, &mut held), drawn);
        assert!(held.notes.is_empty(), "{:?}", held.notes);

        // Held beside it, the index is read: the row of seven columns on chr2
        // is never read for a window on chr1.
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/read/fixtures/indexed");
        let mut held = Held::new();
        held.insert(
            "seven.vcf.gz",
            fs::read(fixtures.join("seven.vcf.gz")).unwrap(),
        );
        let line = "chr1:1-1,000 seven.vcf.gz";
        assert!(build_files(&invocation(line), &mut held, |_, _| None).is_err());
        held.insert(
            "seven.vcf.gz.tbi",
            fs::read(fixtures.join("seven.vcf.gz.tbi")).unwrap(),
        );
        assert!(build_files(&invocation(line), &mut held, |_, _| None).is_ok());
    }

    /// A file a command names that is not held says which are, a BAM asked
    /// for as text says what it is, and nothing is piped into files held.
    #[test]
    fn a_file_not_held_says_which_files_are() {
        let mut held = Held::new();
        let calls = Source::Path("calls.vcf".into());
        let error = held.text(&calls).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::NotFound);
        assert_eq!(
            error.to_string(),
            "no file called calls.vcf is held, and no files are"
        );
        held.insert("genes.gff3", GENES);
        held.insert("tiny.bam", crate::read::bam::fixture::BAM);
        assert_eq!(
            held.text(&calls).unwrap_err().to_string(),
            "no file called calls.vcf is held; the files held are genes.gff3, tiny.bam"
        );
        let error = held.text(&Source::Path("tiny.bam".into())).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("BAM"), "{error}");
        assert!(held.text(&Source::Stdin).is_err());
        assert!(held.contains("genes.gff3") && !held.contains("calls.vcf"));
    }

    /// A byte order mark is dropped where a file is decoded, so no reader
    /// ever sees one: held in memory or read from disk, plain or out of a
    /// gzip wrapper. Only the one at the start is a mark; a second is the
    /// file's own text and is left for the reader to refuse.
    #[test]
    fn a_byte_order_mark_is_dropped_where_a_file_is_decoded() {
        const ROW: &str = "chr1\t0\t10\n";
        // "\u{feff}chr1\t0\t10\n" as `gzip` writes it.
        const GZIPPED: [u8; 33] = [
            0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x13, 0x7b, 0xbf, 0x7b, 0x7f,
            0x72, 0x46, 0x91, 0x21, 0xa7, 0x01, 0xa7, 0xa1, 0x01, 0x17, 0x00, 0x1b, 0x49, 0x84,
            0x93, 0x0d, 0x00, 0x00, 0x00,
        ];
        let mut held = Held::new();
        held.insert("marked.bed", format!("\u{feff}{ROW}"));
        held.insert("twice.bed", format!("\u{feff}\u{feff}{ROW}"));
        held.insert("marked.bed.gz", GZIPPED);
        let text = |held: &mut Held, name: &str| held.text(&Source::Path(name.into())).unwrap();
        assert_eq!(text(&mut held, "marked.bed"), ROW);
        assert_eq!(text(&mut held, "twice.bed"), format!("\u{feff}{ROW}"));
        assert_eq!(text(&mut held, "marked.bed.gz"), ROW);
        let dir = Scratch::new("bom");
        let path = dir.write("marked.bed", format!("\u{feff}{ROW}").as_bytes());
        let mut disk = Disk::default();
        assert_eq!(disk.text(&Source::Path(path.into())).unwrap(), ROW);
    }

    /// A BAM named on its own is its depth, and over a window a few reads
    /// wide the figure says the reads could be drawn: a simulated user drew
    /// four hundred bases of depth where the task wanted the reads.
    #[test]
    fn a_bam_named_on_its_own_over_a_few_reads_says_they_could_be_drawn() {
        let dir = Scratch::new("bam-note");
        let bam = dir.write("tiny.bam", &crate::read::bam::fixture::BAM);
        dir.write("tiny.bam.bai", &crate::read::bam::fixture::BAI);
        let noted = |line: &str| {
            let args: Vec<String> = line.split_whitespace().map(String::from).collect();
            let Request::Draw(invocation) = parse(&args).unwrap() else {
                unreachable!("a figure")
            };
            let mut disk = Disk::default();
            build_files(&invocation, &mut disk, |_, _| None).unwrap();
            disk.notes
        };
        assert_eq!(
            noted(&format!("chr1:10-35 {bam}")),
            [format!(
                "{bam} is drawn as its depth; --pileup {bam} draws its reads"
            )]
        );
        assert!(
            noted(&format!("chr1:10-35 --coverage {bam}")).is_empty(),
            "the depth was asked for by name"
        );
    }

    /// A figure says a thing once however many times it is told: a BAM over
    /// two places a few reads wide each said it was drawn as its depth once a
    /// panel, word for word, and so it did where each panel found its reads
    /// only on a second try, under the name `--rename` gives the sequence.
    #[test]
    fn a_note_is_said_once_however_many_panels_say_it() {
        use crate::read::bam::fixture::{BAI, BAM};
        let dir = Scratch::new("notes");
        let bam = dir.write("tiny.bam", &BAM);
        dir.write("tiny.bam.bai", &BAI);
        let said = format!("{bam} is drawn as its depth; --pileup {bam} draws its reads");
        for line in [
            format!("chr1:10-35 chr1:12-30 {bam}"),
            format!("1:10-35 1:12-30 {bam} --rename chr1=1"),
        ] {
            let mut disk = Disk::default();
            let svg = build_files(&invocation(&line), &mut disk, |_, _| None).unwrap();
            assert_eq!(
                svg.matches("<svg").count(),
                3,
                "{line}: a sheet of two panels"
            );
            assert_eq!(disk.notes, std::slice::from_ref(&said), "{line}");
            let mut held = Held::new();
            held.insert(bam.as_str(), BAM);
            held.insert(format!("{bam}.bai"), BAI);
            build_files(&invocation(&line), &mut held, |_, _| None).unwrap();
            assert_eq!(held.notes, std::slice::from_ref(&said), "{line}");
        }
    }

    /// The bytes a source holds and the file beside it, as every `Files` the
    /// command line and a page draw through answers for them. One that does
    /// not pass both on turns off every reader built on them, and nothing
    /// says so.
    #[test]
    fn every_files_answers_for_a_file_s_bytes_and_the_file_beside_it() {
        use crate::read::bam::fixture::{BAI, BAM};
        let dir = Scratch::new("door");
        let bam = dir.write("tiny.bam", &BAM);
        let replaced = dir.write("tiny.bai", &BAI);
        let source = Source::Path(bam.clone().into());
        let mut held = Held::new();
        held.insert(bam.as_str(), BAM);
        held.insert(replaced.as_str(), BAI);
        let answers = |files: &mut dyn Files| {
            let mut bytes = Vec::new();
            files
                .seekable(&source)
                .unwrap()
                .expect("the file's bytes")
                .read_to_end(&mut bytes)
                .unwrap();
            let beside = files.beside(&source, ".bai").unwrap().expect("its index");
            (bytes, beside.name, beside.bytes)
        };
        let wanted = (BAM.to_vec(), replaced.clone(), BAI.to_vec());
        assert_eq!(answers(&mut Disk::default()), wanted);
        assert_eq!(answers(&mut held), wanted);
        let mut disk = Disk::default();
        let mut kept = KeptStdin {
            files: &mut disk,
            stdin: None,
        };
        assert_eq!(answers(&mut kept), wanted);
        let mut kept = KeptStdin {
            files: &mut held,
            stdin: None,
        };
        assert_eq!(answers(&mut kept), wanted);

        // The ending after the whole name is looked for first, as samtools
        // looks; nothing is beside it under another ending; and standard input
        // has neither bytes to go back in nor a name to put an ending on.
        let after = dir.write("tiny.bam.bai", b"first");
        held.insert(after.as_str(), &b"first"[..]);
        let mut disk = Disk::default();
        for files in [&mut disk as &mut dyn Files, &mut held] {
            let found = files.beside(&source, ".bai").unwrap().unwrap();
            assert_eq!(
                (found.name, found.bytes, found.older),
                (after.clone(), b"first".to_vec(), false)
            );
            assert!(files.beside(&source, ".csi").unwrap().is_none());
            assert!(files.seekable(&Source::Stdin).unwrap().is_none());
            assert!(files.beside(&Source::Stdin, ".bai").unwrap().is_none());
        }
        // A file not held is left to the text, which names the files that are.
        let other = Source::Path("other.bam".into());
        assert!(held.seekable(&other).unwrap().is_none());
        assert!(held
            .text(&other)
            .unwrap_err()
            .to_string()
            .contains("tiny.bam"));
    }

    /// An index older than its file may be for an earlier version of it, and
    /// the times are compared in whole seconds, since an archive keeps no
    /// finer a time; a time not known is not older.
    #[test]
    fn a_file_is_older_than_another_by_whole_seconds_and_never_by_an_unknown_time() {
        use std::time::{Duration, UNIX_EPOCH};
        let at = |seconds: u64, nanos: u32| Some(UNIX_EPOCH + Duration::new(seconds, nanos));
        assert!(older(at(9, 999_999_999), at(10, 0)));
        assert!(!older(at(10, 0), at(10, 999_999_999)), "the same second");
        assert!(!older(at(11, 0), at(10, 0)));
        assert!(!older(None, at(10, 0)) && !older(at(9, 0), None));
    }

    /// The file beside one on disk says whether it was written before it, as
    /// an index made for an earlier version of its file was. `touch` sets the
    /// time here, since `File::set_modified` is newer than the oldest Rust
    /// the crate builds with.
    #[cfg(unix)]
    #[test]
    fn a_file_beside_another_on_disk_says_whether_it_is_older() {
        use crate::read::bam::fixture::{BAI, BAM};
        let dir = Scratch::new("older");
        let bam = dir.write("tiny.bam", &BAM);
        let bai = dir.write("tiny.bam.bai", &BAI);
        let source = Source::Path(bam.clone().into());
        let found = Disk::default().beside(&source, ".bai").unwrap().unwrap();
        assert!(!found.older, "written after its file");
        let touched = std::process::Command::new("touch")
            .args(["-t", "200001010000", &bai])
            .status()
            .unwrap();
        assert!(touched.success());
        let found = Disk::default().beside(&source, ".bai").unwrap().unwrap();
        assert!(found.older, "written in 2000, before its file");
        // Bytes held in memory say nothing of when they were written.
        let mut held = Held::new();
        held.insert(bam.as_str(), BAM);
        held.insert(bai.as_str(), BAI);
        assert!(!held.beside(&source, ".bai").unwrap().unwrap().older);
    }

    /// A pipe the shell names, as `<(cat depth.bg)`, is read from its first
    /// byte by every question asked of it, and draws what the file it came
    /// from draws. Asking whether it was a BAM read its first three bytes and
    /// dropped them, so its `chr1` rows were read as rows of `1` and a figure
    /// placed on `chr1` was refused, and a BAM was refused as text that is not
    /// UTF-8.
    #[cfg(unix)]
    #[test]
    fn a_pipe_named_as_a_file_is_read_from_its_first_byte() {
        use crate::read::bam::fixture::BAM;
        use std::os::unix::io::AsRawFd;
        use std::process::{Command, Stdio};
        let dir = Scratch::new("pipe");
        let depth = dir.write("depth.bg", b"chr1\t0\t100\t5\n");
        let bam = dir.write("reads.bam", &BAM);
        for (line, file) in [
            ("chr1:1-100 --coverage", &depth),
            ("chr1 --coverage", &depth),
            ("chr1 --features", &depth),
            ("chr1:10-35 --coverage", &bam),
            ("chr1:10-35 --pileup", &bam),
        ] {
            let mut child = Command::new("cat")
                .arg(file)
                .stdout(Stdio::piped())
                .spawn()
                .unwrap();
            let pipe = child.stdout.take().unwrap();
            let fd = pipe.as_raw_fd();
            let drawn = drawn_from_disk(&format!("{line} /dev/fd/{fd} --label read"));
            drop(pipe);
            child.wait().unwrap();
            let drawn = drawn.unwrap_or_else(|error| panic!("{line}: {error}"));
            let from_file = drawn_from_disk(&format!("{line} {file} --label read")).unwrap();
            assert!(drawn == from_file, "{line}: the pipe drew another figure");
        }
    }

    /// Bases too narrow for their letters say how wide the figure would have
    /// to be for them, and at that width they are letters.
    #[test]
    fn blocks_of_bases_say_the_width_their_letters_need() {
        let fasta = format!(">chr1\n{}\n", "ACGT".repeat(500));
        let held = [("ref.fa", fasta.as_str())];
        let (_, notes) = drawn_noting("chr1:1-400 ref.fa", &held);
        let note = notes
            .iter()
            .find(|note| note.starts_with("the bases are blocks of colour"))
            .expect("a note on the width");
        let width: u64 = note
            .rsplit("--width ")
            .next()
            .and_then(|rest| rest.split_whitespace().next())
            .and_then(|number| number.parse().ok())
            .expect("a width");
        let (svg, notes) = drawn_noting(&format!("chr1:1-400 ref.fa --width {width}"), &held);
        assert!(notes.is_empty(), "{notes:?}");
        let letters = svg.unwrap().matches(">A</text>").count();
        assert_eq!(letters, 100, "a hundred As in four hundred bases of ACGT");
    }

    /// A figure of several places says once the width that letters the bases
    /// of every panel. Each panel said the width its own span needed, so two
    /// places of different lengths gave two widths, and a reader who drew the
    /// figure at the first still had the longer place in blocks.
    #[test]
    fn bases_in_several_panels_say_the_one_width_all_their_letters_need() {
        let fasta = format!(">chr1\n{}\n", "ACGT".repeat(500));
        let held = [("ref.fa", fasta.as_str())];
        let widths = |line: &str| -> Vec<u64> {
            let (svg, notes) = drawn_noting(line, &held);
            svg.unwrap_or_else(|error| panic!("{line}: {error}"));
            notes
                .iter()
                .filter(|note| note.starts_with("the bases are blocks of colour"))
                .map(|note| {
                    let rest = note.rsplit("--width ").next().unwrap_or_default();
                    let number = rest.split_whitespace().next().unwrap_or_default();
                    number.parse().unwrap_or_else(|_| panic!("{note}"))
                })
                .collect()
        };
        let shorter = widths("chr1:1-400 ref.fa");
        let longer = widths("chr1:1001-1600 ref.fa");
        assert!(
            shorter.len() == 1 && longer.len() == 1,
            "{shorter:?} {longer:?}"
        );
        assert!(longer[0] > shorter[0], "{shorter:?} {longer:?}");
        // The longer place's width, whichever panel it is.
        assert_eq!(widths("chr1:1-400 chr1:1001-1600 ref.fa"), longer);
        assert_eq!(widths("chr1:1001-1600 chr1:1-400 ref.fa"), longer);
        // And at that width both panels are letters, and nothing is said.
        let line = format!("chr1:1-400 chr1:1001-1600 ref.fa --width {}", longer[0]);
        let (svg, notes) = drawn_noting(&line, &held);
        assert!(notes.is_empty(), "{notes:?}");
        let letters = svg.unwrap().matches(">A</text>").count();
        assert_eq!(letters, 100 + 150, "the As of four and six hundred bases");
    }

    /// A compressed file is text in a wrapper, and is read as the text.
    #[test]
    fn a_compressed_file_on_disk_is_read_as_the_text_inside() {
        let dir = Scratch::new("gzip");
        // "chr1 1 9 gene0 0 +" as a GFF3 row, compressed by gzip -9 -n.
        let plain = "##gff-version 3\nchr1\t.\tgene\t10\t90\t.\t+\t.\tID=g1;Name=abcA\n";
        let gz = dir.write("genes.gff3.gz", &gzip_of(plain.as_bytes()));
        let svg = drawn_from_disk(&format!("chr1:1-100 --features {gz}")).unwrap();
        assert!(svg.contains("abcA"), "{svg}");
    }

    /// A gzip member holding `data` in one stored block, which is all a test
    /// needs: the reader takes any kind of block.
    fn gzip_of(data: &[u8]) -> Vec<u8> {
        let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3];
        out.push(1);
        let length = data.len() as u16;
        out.extend_from_slice(&length.to_le_bytes());
        out.extend_from_slice(&(!length).to_le_bytes());
        out.extend_from_slice(data);
        let mut crc = !0u32;
        for byte in data {
            crc ^= u32::from(*byte);
            for _ in 0..8 {
                crc = if crc & 1 != 0 {
                    0xedb8_8320 ^ (crc >> 1)
                } else {
                    crc >> 1
                };
            }
        }
        out.extend_from_slice(&(!crc).to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
        out
    }

    /// A command line of words, drawn from files held in memory by name.
    fn drawn_from(line: &str, held: &[(&str, &str)]) -> Result<String, BuildError> {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        let Request::Draw(invocation) = parse(&args).unwrap() else {
            unreachable!("a figure")
        };
        let held: Vec<(String, String)> = held
            .iter()
            .map(|(name, text)| ((*name).to_string(), (*text).to_string()))
            .collect();
        build(&invocation, |source| {
            let Source::Path(path) = source else {
                unreachable!("every source is a file")
            };
            let path = path.to_string_lossy();
            held.iter()
                .find(|(name, _)| *name == path)
                .map(|(_, text)| text.clone())
                .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, path.to_string()))
        })
    }

    /// Standard input as a pipe gives it: once, and nothing the second time.
    struct Piped {
        text: Option<String>,
        files: Vec<(String, String)>,
    }

    impl Files for Piped {
        fn text(&mut self, source: &Source) -> io::Result<String> {
            match source {
                Source::Stdin => self
                    .text
                    .take()
                    .ok_or_else(|| io::Error::other("standard input was read twice")),
                Source::Path(path) => {
                    let path = path.to_string_lossy();
                    self.files
                        .iter()
                        .find(|(name, _)| *name == path)
                        .map(|(_, text)| text.clone())
                        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, path.to_string()))
                }
            }
        }
    }

    /// The figure `line` draws with `piped` on standard input and `held` as
    /// its files.
    fn drawn_piped(line: &str, piped: &str, held: &[(&str, &str)]) -> Result<String, BuildError> {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        let Request::Draw(invocation) = parse(&args).unwrap() else {
            unreachable!("a figure")
        };
        let mut files = Piped {
            text: Some(piped.to_string()),
            files: held
                .iter()
                .map(|(name, text)| ((*name).to_string(), (*text).to_string()))
                .collect(),
        };
        build_files(&invocation, &mut files, |_, _| None)
    }

    /// The figure `line` draws from `held`, and the notes it made.
    fn drawn_noting(
        line: &str,
        held: &[(&str, &str)],
    ) -> (Result<String, BuildError>, Vec<String>) {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
        let Request::Draw(invocation) = parse(&args).unwrap() else {
            unreachable!("a figure")
        };
        let mut files = Held::new();
        for (name, text) in held {
            files.insert(*name, *text);
        }
        let drawn = build_files(&invocation, &mut files, |_, _| None);
        (drawn, files.notes)
    }

    /// A threshold with no style draws no support, and the command that ran
    /// said nothing of it.
    #[test]
    fn a_tree_says_what_it_was_asked_for_and_did_not_draw() {
        let held = [("t.nwk", "((A:0.1,B:0.2)0.9:0.3,(C:0.15,D:0.05)0.4:0.2);")];
        let (svg, notes) = drawn_noting("--tree t.nwk --threshold 0.7", &held);
        let svg = svg.unwrap();
        assert_eq!(
            notes,
            ["--tree t.nwk: no support is drawn: a threshold was set and no support style"]
        );
        assert!(svg.contains("no support is drawn"), "{svg}");
        let (_, notes) = drawn_noting("--tree t.nwk --threshold 0.7 --support-style labels", &held);
        assert!(notes.is_empty(), "{notes:?}");
    }

    /// The ruler of columns an alignment is read against goes under it, and a
    /// tree below the alignment stays below the ruler, which does not measure
    /// it.
    #[test]
    fn the_ruler_of_columns_goes_under_the_alignment_and_not_under_a_tree() {
        let held = [("aln.fa", ALIGNMENT), ("t.nwk", ROWS_TREE)];
        let (svg, _) = drawn_noting("--msa aln.fa --tree t.nwk", &held);
        let svg = svg.unwrap();
        assert!(
            svg.contains("drawn top to bottom: aln, column and a phylogeny"),
            "{}",
            &svg[svg.find("<desc").unwrap_or(0)..][..200]
        );
    }

    /// A NEXUS file, as BEAST, MrBayes and FigTree write one, is read as a
    /// tree, named on its own too, and a file of several trees draws its
    /// first and says so. It was read as Newick and refused.
    #[test]
    fn a_nexus_tree_is_drawn_and_a_file_of_several_says_which() {
        let nexus = "#NEXUS\nbegin trees;\n translate 1 A, 2 B, 3 C;\n\
                     tree first = ((1:1,2:1)[&posterior=0.97]:1,3:2);\n\
                     tree second = ((1:1,3:1):1,2:2);\nend;";
        let (svg, notes) = drawn_noting("--tree run.trees", &[("run.trees", nexus)]);
        let svg = svg.unwrap();
        assert!(
            svg.contains(">A</text>") && svg.contains(">C</text>"),
            "{svg}"
        );
        assert_eq!(
            notes,
            ["--tree run.trees holds 2 trees, and the first is drawn"]
        );
        for alone in ["mcc.nex", "run.trees", "run.nexus", "run.nxs"] {
            let (drawn, _) = drawn_noting(alone, &[(alone, nexus)]);
            assert_eq!(drawn.unwrap(), svg, "{alone} named on its own");
        }
    }

    /// A BEAST tree's support, kept as its posterior, is drawn when asked
    /// for by that name, and a name no clade carries is refused with the ones
    /// that are.
    #[test]
    fn support_is_drawn_from_the_annotation_the_command_names() {
        let held = [(
            "mcc.tree",
            "((A:1,B:1)[&posterior=0.97,height=2]:1,(C:1,D:1)[&posterior=0.42]:1);",
        )];
        let (svg, _) = drawn_noting(
            "--tree mcc.tree --support-from posterior --support-style labels",
            &held,
        );
        assert!(svg.unwrap().contains("0.97"));
        let (refused, _) = drawn_noting("--tree mcc.tree --support-from prob", &held);
        let error = refused.unwrap_err().to_string();
        assert!(
            error.contains("prob") && error.contains("posterior") && error.contains("height"),
            "{error}"
        );
    }

    /// Four rows of eight columns, and a tree of the same four samples in
    /// another order.
    const ALIGNMENT: &str = ">one\nACGTACGT\n>two\nACGTACGA\n>three\nACGAACGT\n>four\nTCGTACGT\n";
    const ROWS_TREE: &str = "((four:1,two:1):1,(three:1,one:1):1);";

    /// The rows in the order drawn, top to bottom, read off the names written
    /// beside them.
    fn rows_drawn(svg: &str) -> Vec<String> {
        let mut rows: Vec<(f64, String)> = Vec::new();
        for piece in svg.split("<text ").skip(1) {
            let Some(content) = piece
                .split('>')
                .nth(1)
                .and_then(|text| text.strip_suffix("</text"))
            else {
                continue;
            };
            if !["one", "two", "three", "four"].contains(&content)
                || rows.iter().any(|(_, held)| held == content)
            {
                continue;
            }
            let y: f64 = piece
                .split("y=\"")
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .and_then(|y| y.parse().ok())
                .unwrap_or(f64::NAN);
            rows.push((y, content.to_string()));
        }
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        rows.into_iter().map(|(_, name)| name).collect()
    }

    /// An alignment is its own place, as a tree is: it was refused until a
    /// region was made up for it, and the one to make up was a row's name.
    #[test]
    fn an_alignment_is_drawn_over_its_own_columns_when_no_place_is_named() {
        let held = [("aln.fa", ALIGNMENT)];
        let (svg, notes) = drawn_noting("--msa aln.fa", &held);
        let svg = svg.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        assert!(
            svg.contains("aln.fa:1-8"),
            "the columns as the place: {svg}"
        );
        assert_eq!(rows_drawn(&svg), ["one", "two", "three", "four"]);
        // The same as naming the columns by hand.
        let (named, _) = drawn_noting("aln.fa:1-8 --msa aln.fa", &held);
        assert_eq!(svg, named.unwrap());
        // A pipe is read once, and its width is needed before its rows, so
        // the figure keeps what it read.
        let piped = drawn_piped("--msa -", ALIGNMENT, &[]).unwrap();
        assert_eq!(rows_drawn(&piped), ["one", "two", "three", "four"]);
        assert!(piped.contains(":1-8"), "the columns as the place");
    }

    const COUNTS: &str = "week\tlineage\tcount\ttotal\n1\tA\t9\t10\n1\tB\t1\t10\n\
                          2\tA\t6\t10\n2\tB\t4\t10\n3\tA\t2\t12\n3\tB\t10\t12\n";

    /// A table of counts takes its alerts, its sampling floor and its metric
    /// from the command line.
    #[test]
    fn a_table_of_counts_takes_its_alerts_floor_and_metric() {
        let held = [("f.tsv", COUNTS)];
        let (svg, _) = drawn_noting(
            "--frequencies f.tsv --style line --threshold 0.8 --growth 0.3",
            &held,
        );
        let svg = svg.unwrap();
        assert!(svg.contains(">≥ 80% or up 30 points</text>"), "{svg}");
        assert!(svg.contains("frequency alert &gt;= 0.8"), "{svg}");
        let (counted, _) = drawn_noting("--frequencies f.tsv --counts", &held);
        let counted = counted.unwrap();
        assert!(!counted.contains(">100%</text>"), "{counted}");
        assert!(
            counted.contains(">12</text>"),
            "a ceiling of twelve samples"
        );
        // Weeks one and two had ten samples each and week three twelve.
        let (floored, _) = drawn_noting("--frequencies f.tsv --min-total 11", &held);
        let floored = floored.unwrap();
        assert!(!floored.contains("A | time 1 |"), "{floored}");
        assert!(floored.contains("A | time 3 |"), "{floored}");
    }

    /// A table over time is its own place, as an alignment is, and the ruler
    /// under it counts its weeks: week 1 is printed 1, not 0 and not `1 bp`.
    #[test]
    fn counts_over_time_are_drawn_over_their_own_weeks() {
        let held = [("f.tsv", COUNTS)];
        let (svg, notes) = drawn_noting("--frequencies f.tsv", &held);
        let svg = svg.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        assert!(
            svg.contains(">week</text>"),
            "the ruler says what it counts"
        );
        for week in ["1", "2", "3"] {
            assert!(
                svg.contains(&format!(">{week}</text>")),
                "week {week}: {svg}"
            );
        }
        // No bases, and no locus at the top right saying what the ruler says.
        assert!(
            !svg.contains(" bp<") && !svg.contains(">week:1-3</text>"),
            "{svg}"
        );
        // A tooltip calls a week what the ruler under it calls it.
        assert!(svg.contains("A | time 1 | count 9 of 10"), "{svg}");
        // The same as naming the weeks by hand.
        let (named, _) = drawn_noting("week:1-3 --frequencies f.tsv", &held);
        assert_eq!(svg, named.unwrap());
        // An estimate beside the counts widens the place to both of them.
        let held = [
            ("f.tsv", COUNTS),
            (
                "r.tsv",
                "week\tmean\tlower\tupper\n2\t1.2\t0.9\t1.5\n5\t0.8\t0.6\t1.1\n",
            ),
        ];
        let (both, _) = drawn_noting("--frequencies f.tsv --phylodynamics r.tsv", &held);
        assert!(
            both.unwrap().contains(">5</text>"),
            "week 5 is on the ruler"
        );
        // A track measured in bases puts the ruler back in bases.
        let held = [("f.tsv", COUNTS), ("d.bg", "chr1\t0\t10\t5\n")];
        let (mixed, _) = drawn_noting("chr1:1-10 --coverage d.bg --frequencies f.tsv", &held);
        assert!(!mixed.unwrap().contains(">week</text>"));
    }

    /// A table of its own place read from a pipe is drawn over its own extent,
    /// as it is from a file, though that extent is read before its rows: the
    /// figure is the same one.
    #[test]
    fn a_table_of_its_own_place_from_a_pipe_is_its_own_place() {
        let skyline = "week\tmean\tlower\tupper\n2\t1.2\t0.9\t1.5\n5\t0.8\t0.6\t1.1\n";
        let fel = "site,alpha,beta,p-value\n1,1.0,0.2,0.8\n2,0.5,3.0,0.01\n";
        let signal = "90\n95\n120\n118\n101\n";
        for (flag, table) in [
            ("--frequencies", COUNTS),
            ("--phylodynamics", skyline),
            ("--selection", fel),
            ("--squiggle", signal),
        ] {
            let piped = drawn_piped(&format!("{flag} - --label x"), table, &[]).unwrap();
            let (named, _) = drawn_noting(&format!("{flag} t.tsv --label x"), &[("t.tsv", table)]);
            assert_eq!(piped, named.unwrap(), "{flag}");
        }
    }

    /// FEL writes p-values and a Bayes empirical Bayes table posteriors, and
    /// the threshold is read in whichever the table holds.
    #[test]
    fn a_selection_table_is_read_by_the_evidence_it_holds() {
        let held = [(
            "fel.csv",
            "site,alpha,beta,p-value\n1,1.0,0.2,0.8\n2,0.5,3.0,0.01\n",
        )];
        let (svg, _) = drawn_noting("--selection fel.csv", &held);
        let svg = svg.unwrap();
        assert!(
            svg.contains(">site</text>") && svg.contains("p ≤ 0.05"),
            "{svg}"
        );
        let (svg, _) = drawn_noting("--selection fel.csv --threshold 0.001", &held);
        assert!(svg.unwrap().contains("p ≤ 0.001"));
        let held = [("b.csv", "site,omega,posterior\n1,0.3,0.1\n2,4.0,0.97\n")];
        let (svg, _) = drawn_noting("--selection b.csv --threshold 0.95", &held);
        assert!(svg.unwrap().contains("PP ≥ 0.95"));
    }

    /// The basecaller's record of the read drawn puts each base it called
    /// over the stretch of current it was called from, the record of that
    /// read and not of another.
    #[test]
    fn a_move_table_puts_the_bases_over_the_signal_of_the_read_drawn() {
        let slow5 = "#slow5_version\t0.2.0\n\
                     #read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\t\
                     len_raw_signal\traw_signal\n\
                     r1\t0\t2048\t0\t2048\t4000\t20\t80,81,80,79,95,96,94,95,70,71,70,69,70,90,91,90,89,90,91,90\n\
                     r2\t0\t2048\t0\t2048\t4000\t4\t70,75,72,71\n";
        let sam = "r2\t4\t*\t0\t0\t*\t*\t0\t0\tG\t*\tmv:B:c,4,1\n\
                   r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\tmv:B:c,4,1,1,1,0,1\n";
        let held = [("reads.slow5", slow5), ("calls.sam", sam)];
        let (svg, _) = drawn_noting("reads.slow5 --with-moves calls.sam", &held);
        let svg = svg.unwrap();
        for base in ["A", "C", "G", "T"] {
            assert!(svg.contains(&format!(">{base}</text>")), "no {base}: {svg}");
        }
        // A record for another read is not taken for this one.
        let (other, _) = drawn_noting(
            "reads.slow5 --with-moves calls.sam",
            &[
                ("reads.slow5", slow5),
                (
                    "calls.sam",
                    "r2\t4\t*\t0\t0\t*\t*\t0\t0\tG\t*\tmv:B:c,4,1\n",
                ),
            ],
        );
        let error = other.unwrap_err().to_string();
        assert!(error.contains("no read called r1; it holds r2"), "{error}");
    }

    /// A SLOW5 holds many reads and the figure draws one, named for it, and
    /// says which when it had to choose.
    #[test]
    fn a_read_s_signal_is_drawn_over_its_samples_and_named_for_its_read() {
        let slow5 = "#slow5_version\t0.2.0\n\
                     #read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\t\
                     len_raw_signal\traw_signal\n\
                     r1\t0\t2048\t0\t2048\t4000\t4\t80,90,85,95\n\
                     r2\t0\t2048\t0\t2048\t4000\t3\t70,75,72\n";
        let held = [("reads.slow5", slow5)];
        // Named on its own, a SLOW5 is a squiggle.
        let (svg, notes) = drawn_noting("reads.slow5", &held);
        let svg = svg.unwrap();
        assert!(
            svg.contains(">r1</text>") && svg.contains(">sample</text>"),
            "{svg}"
        );
        assert!(
            notes.len() == 1 && notes[0].contains("--read NAME"),
            "{notes:?}"
        );
        let (second, notes) = drawn_noting("reads.slow5 --read r2", &held);
        assert!(second.unwrap().contains(">r2</text>"));
        assert!(notes.is_empty(), "{notes:?}");
        let (missing, _) = drawn_noting("reads.slow5 --read r9", &held);
        let error = missing.unwrap_err().to_string();
        assert!(error.contains("r1, r2"), "{error}");
    }

    /// A skyline in decimal years is a continuous time: its ruler writes the
    /// years it spans, its tooltips the times as written, and a place written
    /// in years is read in years.
    #[test]
    fn times_with_fractions_are_a_continuous_time() {
        let skyline = "year\tmedian\tlower\tupper\n2010.25\t100\t50\t200\n\
                       2012.5\t400\t300\t600\n2015.75\t900\t700\t1200\n";
        let held = [("sky.tsv", skyline)];
        let (svg, _) = drawn_noting("--phylodynamics sky.tsv", &held);
        let svg = svg.unwrap();
        assert!(svg.contains("<title>time 2010.25 |"), "{svg}");
        assert!(
            svg.contains(">2012</text>") || svg.contains(">2013</text>"),
            "{svg}"
        );
        assert!(!svg.contains("2,01"), "a year is not grouped: {svg}");
        // Counts in whole weeks beside it are read to the same thousandths.
        let (both, _) = drawn_noting(
            "--frequencies f.tsv --phylodynamics sky.tsv",
            &[
                ("f.tsv", "year\tlineage\tcount\ttotal\n2011\tA\t3\t9\n"),
                ("sky.tsv", skyline),
            ],
        );
        assert!(both.unwrap().contains("A | time 2011 |"));
        // A place is written in years.
        let (placed, _) = drawn_noting("year:2012-2013 --phylodynamics sky.tsv", &held);
        let placed = placed.unwrap();
        assert!(
            placed.contains("time 2012.5") && !placed.contains("time 2010.25"),
            "{placed}"
        );
    }

    /// A table with fractions on standard input is a continuous time as it is
    /// from a file: the figure reads every table's times before it draws any,
    /// and keeps what the pipe gave.
    #[test]
    fn fractions_on_standard_input_are_a_continuous_time() {
        let skyline = "year\tmedian\tlower\tupper\n2010.25\t100\t50\t200\n\
                       2012.5\t400\t300\t600\n2015.75\t900\t700\t1200\n";
        let piped = drawn_piped("--phylodynamics -", skyline, &[]).unwrap();
        assert!(piped.contains("<title>time 2010.25 |"), "{piped}");
    }

    /// A genetic map is drawn as a line of its rates, named on its own by the
    /// name the panels give their maps.
    #[test]
    fn a_genetic_map_is_a_line_of_rates() {
        let map = "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n\
                   1\t101\t2.5\t0\n1\t501\t40\t0.001\n1\t901\t0.5\t0.017\n";
        let held = [("genetic_map_chr1.txt", map)];
        let (svg, _) = drawn_noting("1:1-1000 genetic_map_chr1.txt", &held);
        let svg = svg.unwrap();
        assert!(
            svg.contains(">cM/Mb</text>"),
            "the axis says what it measures: {svg}"
        );
        assert!(svg.contains("<polyline") || svg.contains("<path"), "{svg}");
        let (flag, _) = drawn_noting("1:1-1000 --recombination genetic_map_chr1.txt", &held);
        assert_eq!(svg, flag.unwrap());
    }

    /// A table of windows is a heatmap of its samples, and read against each
    /// sample's own median it shows what changed rather than who was
    /// sequenced deeper.
    #[test]
    fn a_heatmap_draws_every_sample_in_its_windows_and_keys_its_ramp() {
        let depths = "chrom\tstart\tend\tdeep\tshallow\n\
                      c1\t0\t100\t100\t20\nc1\t100\t200\t100\t0\nc1\t200\t300\t200\t20\n";
        let held = [("d.tsv", depths)];
        let (svg, _) = drawn_noting("c1:1-300 --heatmap d.tsv", &held);
        let svg = svg.unwrap();
        assert!(svg.contains("deep, 3 of 3 windows called"), "{svg}");
        assert!(
            svg.contains(">200</text>"),
            "the key tops out at the largest value"
        );
        let (relative, _) = drawn_noting("c1:1-300 --heatmap d.tsv --relative", &held);
        let relative = relative.unwrap();
        // deep: 100, 100, 200 over a median of 100; shallow: 20, 0, 20 over 20.
        // Read either side of its usual depth: the loss in one hue, the gain
        // in the other, and the key writes the middle between them.
        assert!(
            relative.contains(">2×</text>")
                && relative.contains(">1×</text>")
                && relative.contains(">0×</text>"),
            "{relative}"
        );
        let theme = Theme::light();
        for hue in [theme.color(0), theme.color(1)] {
            assert!(
                relative.contains(&format!("fill=\"{hue}\"")),
                "no cell in {hue}"
            );
        }
        // A centre of its own, as for a log ratio, and a long table.
        let ratios = "c1\t0\t100\tS1\t-1\nc1\t100\t200\tS1\t0\nc1\t200\t300\tS1\t2\n";
        let (centred, _) =
            drawn_noting("c1:1-300 --heatmap r.tsv --center 0", &[("r.tsv", ratios)]);
        let centred = centred.unwrap();
        assert!(
            centred.contains(">-1</text>") && centred.contains(">2</text>"),
            "{centred}"
        );
    }

    /// Linkage is a triangle and a few pairs are arcs, unless told; an r² is
    /// read against one.
    #[test]
    fn pairs_are_a_triangle_where_most_were_measured_and_arcs_where_few_were() {
        let ld = " CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n\
                  c1 101 a c1 201 b 0.9\nc1 101 a c1 301 c 0.2\nc1 201 b c1 301 c 0.3\n";
        let few = "pos1\tpos2\tscore\n101\t901\t5\n301\t601\t2\n";
        let held = [("v.ld", ld), ("few.tsv", few)];
        let (triangle, _) = drawn_noting("c1:1-1000 v.ld", &held);
        let triangle = triangle.unwrap();
        assert_eq!(triangle.matches("<polygon").count(), 3, "{triangle}");
        assert!(
            triangle.contains(">1</text>"),
            "an r² is keyed to one: {triangle}"
        );
        let (arcs, _) = drawn_noting("c1:1-1000 --pairs few.tsv", &held);
        let arcs = arcs.unwrap();
        assert_eq!(arcs.matches("<polygon").count(), 0);
        assert!(
            arcs.contains("101 and 901: 5") && arcs.contains(">5</text>"),
            "{arcs}"
        );
        let (told, _) = drawn_noting(
            "c1:1-1000 --pairs v.ld --style arcs --threshold 0.25",
            &held,
        );
        let told = told.unwrap();
        assert!(
            !told.contains("<polygon") && !told.contains("101 and 301"),
            "{told}"
        );
    }

    /// A scan with the linkage of its lead is coloured by it, as LocusZoom
    /// draws one, and says so when the lead is not a variant it tested.
    #[test]
    fn a_scan_is_coloured_by_linkage_with_its_lead() {
        let scan = "CHR\tBP\tP\n1\t100\t1e-9\n1\t200\t1e-5\n1\t300\t0.2\n";
        // As PLINK's --ld-snp writes it: the lead in every row, by name.
        let ld = " CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n\
                  1 100 rs100 1 200 rs200 0.8\n1 100 rs100 1 300 rs300 0.05\n";
        let held = [("g.assoc", scan), ("lead.ld", ld)];
        let (svg, notes) = drawn_noting("1:1-400 g.assoc --ld lead.ld", &held);
        let svg = svg.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        // The scan names no variant, so the lead takes the name the table of
        // linkage gives it.
        assert!(
            svg.contains("lead variant rs100 at 100") && svg.contains("r² with rs100"),
            "{svg}"
        );
        // A scan that names its variants names the lead first.
        let named = "CHR\tSNP\tBP\tP\n1\tvar_a\t100\t1e-9\n1\tvar_b\t200\t1e-5\n";
        let (svg, _) = drawn_noting(
            "1:1-400 g.assoc --ld lead.ld",
            &[("g.assoc", named), ("lead.ld", ld)],
        );
        assert!(svg.unwrap().contains(">var_a</text>"));
        // A genetic map laid over the scan, read off a scale on the right.
        let map = "position rate\n50 2\n250 30\n390 1\n";
        let (svg, _) = drawn_noting(
            "1:1-400 g.assoc --ld lead.ld --with-recombination map.txt",
            &[("g.assoc", scan), ("lead.ld", ld), ("map.txt", map)],
        );
        let svg = svg.unwrap();
        assert!(svg.contains(" cM/Mb</text>"), "{svg}");
        assert!(svg.contains("recombination, cM/Mb"), "keyed: {svg}");
        // A lead the scan did not test has no diamond, and the figure says so.
        let held = [
            ("g.assoc", scan),
            (
                "lead.ld",
                "CHR_A BP_A CHR_B BP_B R2\n1 150 1 200 0.8\n1 150 1 300 0.1\n",
            ),
        ];
        let (svg, notes) = drawn_noting("1:1-400 g.assoc --ld lead.ld", &held);
        assert!(
            !svg.unwrap().contains("lead variant"),
            "no diamond, and none keyed"
        );
        assert!(notes.iter().any(|note| note.contains("150")), "{notes:?}");
    }

    /// The library ordered an alignment's rows by a tree and drew the tree
    /// beside them, and the command line could not ask for it.
    #[test]
    fn a_tree_named_with_the_rows_orders_them_and_is_drawn_beside_them() {
        let held = [("aln.fa", ALIGNMENT), ("t.nwk", ROWS_TREE)];
        // A panel of variable sites keeps its reference, the first row, on a
        // row of its own above the tree, and orders the rest.
        for (track, order) in [
            ("--msa", ["four", "two", "three", "one"]),
            ("--snps", ["one", "four", "two", "three"]),
        ] {
            let (plain, _) = drawn_noting(&format!("{track} aln.fa"), &held);
            let (ordered, notes) =
                drawn_noting(&format!("{track} aln.fa --with-tree t.nwk"), &held);
            let (plain, ordered) = (plain.unwrap(), ordered.unwrap());
            assert!(notes.is_empty(), "{notes:?}");
            assert_eq!(
                rows_drawn(&plain),
                ["one", "two", "three", "four"],
                "{track}"
            );
            assert_eq!(rows_drawn(&ordered), order, "{track}: {ordered}");
            // Every tip has a row, the reference's included.
            assert!(!ordered.contains("has no row"), "{track}: {ordered}");
            assert!(
                ordered.matches("<line").count() > plain.matches("<line").count(),
                "{track}: no branches beside the rows"
            );
        }
        // A tree that names none of the rows orders nothing, and is refused.
        let held = [("aln.fa", ALIGNMENT), ("t.nwk", "(X:1,Y:1);")];
        let (drawn, _) = drawn_noting("--msa aln.fa --with-tree t.nwk", &held);
        let error = drawn.unwrap_err().to_string();
        assert!(
            error.contains("no tip in this file names anything in the rows"),
            "{error}"
        );
        assert!(error.contains("X, Y"), "{error}");
    }

    /// A PAF writes each query's length, and a figure placed on the query by
    /// its name was refused as a place no file named.
    #[test]
    fn a_comparison_is_placed_on_its_query_by_name() {
        let paf = "asm1\t6000\t0\t2000\t+\tasm2\t6200\t0\t2000\t1990\t2000\t60\n\
                   asm1\t6000\t2000\t4000\t-\tasm2\t6200\t2100\t4100\t1990\t2000\t60\n";
        let held = [("a.paf", paf)];
        for line in ["asm1 a.paf", "asm1 --dotplot a.paf"] {
            let (svg, notes) = drawn_noting(line, &held);
            let svg = svg.unwrap();
            assert!(notes.is_empty(), "{line}: {notes:?}");
            assert!(svg.contains("asm1:1-6000"), "{line}: {svg}");
            // Which way round a block runs is keyed, as each is drawn.
            assert!(
                svg.contains(">same strand</text>") && svg.contains(">reversed</text>"),
                "{line}"
            );
        }
    }

    /// An alignment too wide for its letters paints them as colours, which a
    /// key now names, and it is not told to be drawn thousands of pixels wide.
    #[test]
    fn an_alignment_keys_its_colours_and_is_not_told_to_widen() {
        let rows: String = (0..4)
            .map(|row| format!(">r{row}\n{}\n", "ACGT".repeat(150)))
            .collect();
        let held = [("aln.fa", rows.as_str())];
        let (svg, notes) = drawn_noting("--msa aln.fa --style all", &held);
        let svg = svg.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        for base in ["A", "C", "G", "T"] {
            assert!(
                svg.contains(&format!(">{base}</text>")),
                "no key for {base}: {svg}"
            );
        }
    }

    /// PLINK writes 1 for the chromosome a FASTA calls NC_1, and every file
    /// had to call it one name before the two could be drawn together.
    #[test]
    fn a_file_that_calls_the_sequence_otherwise_is_read_by_that_name() {
        let scan = "CHR SNP BP A1 P\n1 rs1 150 A 0.5\n1 rs2 4800 A 1e-9\n";
        let fasta = format!(">NC_1\n{}\n", "ACGT".repeat(1_500));
        let held = [("gwas.assoc", scan), ("ref.fa", fasta.as_str())];
        // Named as the FASTA names it, the scan's rows are found under its own
        // name, and the figure runs the FASTA's length.
        let (svg, notes) = drawn_noting("NC_1 ref.fa gwas.assoc --rename 1=NC_1", &held);
        let svg = svg.unwrap();
        assert_eq!(locus_of(&svg), "NC_1:1-6000");
        assert_eq!(svg.matches("<circle").count(), 2, "{svg}");
        assert!(notes.is_empty(), "{notes:?}");
        // With nothing to give the length, the figure ends where the rows do,
        // and says so.
        let (svg, notes) = drawn_noting("NC_1 gwas.assoc --rename 1=NC_1", &held);
        assert_eq!(locus_of(&svg.unwrap()), "NC_1:1-4800");
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert!(
            notes[0].starts_with("NC_1 is drawn to 4,800, as far as gwas.assoc reaches"),
            "{notes:?}"
        );
        // The rename is given, so the FASTA is all the note asks for.
        assert!(!notes[0].contains("--rename"), "{notes:?}");
        // Placed on the table's own name, the note says how a FASTA that
        // calls it otherwise would draw it too.
        let (_, notes) = drawn_noting("1 gwas.assoc", &held);
        assert!(notes[0].ends_with("add --rename 1=THAT_NAME"), "{notes:?}");
        // Without --rename each refusal says which one to write.
        let error = drawn_from("NC_1 gwas.assoc", &held)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("If 1 is NC_1, add --rename 1=NC_1."),
            "{error}"
        );
        let error = drawn_from("NC_1:1-5000 gwas.assoc", &held)
            .unwrap_err()
            .to_string();
        assert!(
            error.ends_with("; if 1 is NC_1, add --rename 1=NC_1"),
            "{error}"
        );
    }

    /// The command line prints what the figure noted, so the files it reads
    /// through have to keep it.
    #[test]
    fn the_disk_keeps_what_the_figure_notes() {
        let mut disk = Disk::default();
        disk.note("chr1 is drawn to 4,800");
        assert_eq!(disk.notes, ["chr1 is drawn to 4,800"]);
        // Each thing once, in the order it was first said, on disk and in
        // memory alike.
        let mut held = Held::new();
        for files in [&mut disk as &mut dyn Files, &mut held] {
            for message in ["chr1 is drawn to 4,800", "a", "chr1 is drawn to 4,800", "a"] {
                files.note(message);
            }
        }
        assert_eq!(disk.notes, ["chr1 is drawn to 4,800", "a"]);
        assert_eq!(held.notes, ["chr1 is drawn to 4,800", "a"]);
    }

    /// A name that is the figure's own with a chr more or less is plainly
    /// that sequence; among several other names nothing is guessed.
    #[test]
    fn a_rename_is_offered_only_where_one_plainly_fits() {
        assert_eq!(
            rename_for(&[("7".into(), 3), ("8".into(), 2)], "chr7", false),
            Some(("7".to_string(), "chr7".to_string()))
        );
        assert_eq!(
            rename_for(&[("7".into(), 3), ("8".into(), 2)], "NC_1", true),
            None
        );
        assert_eq!(
            rename_for(&[("1".into(), 3)], "NC_1", true),
            Some(("1".to_string(), "NC_1".to_string()))
        );
        assert_eq!(rename_for(&[("1".into(), 3)], "NC_1", false), None);
        assert_eq!(rename_for(&[("NC_1".into(), 3)], "NC_1", true), None);
        // A gene named nearly right is a typo, and no rename is offered for it.
        let error = drawn_from("rpoC genes.gff3", &[("genes.gff3", GENES)]).unwrap_err();
        assert!(!error.to_string().contains("--rename"), "{error}");
    }

    const GENES: &str = "\
##gff-version 3
##sequence-region chr1 1 50000
chr1\t.\tregion\t1\t50000\t.\t+\t.\tID=chr1:1..50000;Name=ANONYMOUS
chr1\t.\tgene\t10001\t12000\t.\t+\t.\tID=gene-A;Name=rpoB;locus_tag=SYN_1
chr1\t.\tCDS\t10001\t11997\t.\t+\t0\tID=cds-A;Parent=gene-A;gene=rpoB
chr1\t.\tgene\t20001\t21000\t.\t-\t.\tID=gene-B;Name=katG
";

    /// The figure `line` builds from `held`, in `theme` and over `window`.
    fn built_from(line: &str, held: &[(&str, &str)], theme: Theme, window: Option<&str>) -> Built {
        let mut files = Held::new();
        for (name, text) in held {
            files.insert(*name, *text);
        }
        let window = window.map(|text| Region::parse(text).unwrap());
        build_figure(
            &invocation(line),
            &mut files,
            |_, _| None,
            theme,
            window.as_ref(),
        )
        .unwrap_or_else(|error| panic!("{line}: {error}"))
    }

    /// A command drawn over another window is the command with that window
    /// written in place of its own, and one placed by a gene keeps the gene
    /// as its title wherever it is moved to.
    #[test]
    fn a_command_is_drawn_over_the_window_a_page_moves_it_to() {
        let held = [("genes.gff3", GENES), ("depth.bg", "chr1\t0\t50000\t12\n")];
        let moved = built_from(
            "rpoB depth.bg genes.gff3",
            &held,
            Theme::light(),
            Some("chr1:11,001-13,000"),
        );
        let written =
            drawn_from("chr1:11,001-13,000 depth.bg genes.gff3 --title rpoB", &held).unwrap();
        assert_eq!(moved.figure.to_svg(), written);
        assert_eq!(moved.along, Region::parse("chr1:11,001-13,000").ok());
        let moved = built_from(
            "chr1:1-100 depth.bg",
            &held,
            Theme::light(),
            Some("chr1:201-300"),
        );
        assert_eq!(
            moved.figure.to_svg(),
            drawn_from("chr1:201-300 depth.bg", &held).unwrap()
        );
        // Where a sequence ends when no file says is a note about the figure
        // as the command draws it, and not about a window moved along it.
        let notes = |window: Option<&str>| {
            let mut files = Held::new();
            files.insert("depth.bg", "chr1\t0\t500\t12\n");
            let window = window.map(|text| Region::parse(text).unwrap());
            build_figure(
                &invocation("chr1 depth.bg"),
                &mut files,
                |_, _| None,
                Theme::light(),
                window.as_ref(),
            )
            .unwrap();
            files.notes
        };
        assert_eq!(notes(None).len(), 1, "{:?}", notes(None));
        assert!(notes(Some("chr1:101-200")).is_empty());
    }

    /// A figure along a genome says where, so a page can move it along; a
    /// figure that is its own place, or has none, says nothing.
    #[test]
    fn a_figure_says_whether_it_runs_along_a_genome() {
        let held = [
            ("genes.gff3", GENES),
            ("depth.bg", "chr1\t0\t50000\t12\n"),
            ("aln.fa", ALIGNMENT),
            ("t.nwk", ROWS_TREE),
            ("f.tsv", COUNTS),
        ];
        let along = |line: &str| {
            built_from(line, &held, Theme::light(), None)
                .along
                .map(|region| region.to_string())
        };
        assert_eq!(
            along("rpoB depth.bg genes.gff3").as_deref(),
            Some("chr1:9801-12200")
        );
        assert_eq!(along("chr1:1-100 depth.bg").as_deref(), Some("chr1:1-100"));
        for own in [
            "--msa aln.fa",
            "--tree t.nwk",
            "--frequencies f.tsv",
            "week:1-3 --frequencies f.tsv",
        ] {
            assert_eq!(along(own), None, "{own}");
        }
    }

    /// A page's own colours: the theme a command is built in is the theme it
    /// is drawn in, down to the ground under it.
    #[test]
    fn a_command_is_built_in_the_theme_it_is_given() {
        let held = [("depth.bg", "chr1\t0\t100\t12\n")];
        let dark = built_from("chr1:1-100 depth.bg", &held, Theme::dark(), None);
        assert_eq!(
            dark.figure.to_svg(),
            drawn_from("chr1:1-100 depth.bg --theme dark", &held).unwrap()
        );
        let light = built_from("chr1:1-100 depth.bg", &held, Theme::light(), None);
        assert_eq!(
            light.figure.to_svg(),
            drawn_from("chr1:1-100 depth.bg", &held).unwrap()
        );
        let mut page = Theme::light();
        page.background = "#fbfaff".to_string();
        let svg = built_from("chr1:1-100 depth.bg", &held, page, None)
            .figure
            .to_svg();
        // The first thing drawn after the definitions, under everything.
        let ground = |svg: &str| {
            svg.split("</defs>")
                .nth(1)
                .and_then(|drawn| drawn.split("<rect").nth(1))
                .and_then(|rect| rect.split("fill=\"").nth(1))
                .and_then(|fill| fill.split('"').next())
                .unwrap_or_default()
                .to_string()
        };
        assert_eq!(ground(&svg), "#fbfaff");
        assert_eq!(ground(&light.figure.to_svg()), Theme::light().background);
    }

    /// Chromosomes in the order a reader counts them: by number, with or
    /// without `chr`, then X, Y and the mitochondrion, then the rest.
    #[test]
    fn chromosomes_are_counted_by_number_then_x_y_and_the_mitochondrion() {
        let mut names = vec![
            "contig_12",
            "chr10",
            "MT",
            "chrX",
            "chr2",
            "contig_7",
            "Y",
            "chr1",
            "chrM",
        ];
        names.sort_by(|a, b| chromosome_order(a, b));
        assert_eq!(
            names,
            [
                "chr1",
                "chr2",
                "chr10",
                "chrX",
                "Y",
                "chrM",
                "MT",
                "contig_7",
                "contig_12"
            ]
        );
    }

    /// A scan read on its own is drawn across every sequence the tables
    /// name, end to end, each as long as the furthest position tested on it,
    /// each named under the scan in place of a ruler of positions nothing
    /// else uses. It was refused for having no region.
    #[test]
    fn a_scan_alone_is_drawn_across_the_whole_genome() {
        let first = "CHR\tBP\tP\nX\t150000000\t0.2\n2\t240000000\t1e-9\n\
                     1\t100000000\t0.3\n10\t130000000\t0.5\nX\t20000000\t0.3\n";
        let second = "CHR\tBP\tP\n3\t198000000\t0.01\n1\t248000000\t0.4\n";
        let held = [("a.assoc", first), ("b.assoc", second)];
        let svg = drawn_from(
            "--manhattan a.assoc --threshold genome-wide --manhattan b.assoc",
            &held,
        )
        .unwrap();
        // Each sequence under the scan says its name and length.
        let at = |name: &str| {
            svg.find(&format!("<title>{name}, "))
                .unwrap_or_else(|| panic!("{name} is not named: {svg}"))
        };
        assert!(at("1") < at("2") && at("2") < at("3") && at("3") < at("10") && at("10") < at("X"));
        // Every sequence of both tables, each as long as its furthest test.
        assert!(svg.contains("genome:1-966000000"), "{svg}");
        assert_eq!(svg.matches(">-log10 p</text>").count(), 2);
        // Every other chromosome a shade lighter, which is what tells a
        // reader where one ends and the next begins.
        let theme = Theme::light();
        let lighter = crate::theme::mix(&theme.muted, theme.surface(), 0.42);
        assert!(
            svg.contains(&format!("fill=\"{lighter}\"")),
            "no chromosome is a shade lighter: {svg}"
        );
        assert!(svg.contains("p = 5e-8"), "{svg}");
        for unit in [" kb</text>", " Mb</text>"] {
            assert!(!svg.contains(unit), "a ruler of positions: {svg}");
        }
        // A table of positions on no sequence has no place on a genome.
        let error =
            drawn_from("--manhattan a.assoc", &[("a.assoc", "BP\tP\n100\t0.5\n")]).unwrap_err();
        assert!(error.to_string().contains("names no sequence"), "{error}");
    }

    /// Where each sequence of a genome-wide figure is named under it, in the
    /// order they are drawn, and how long it is said to be.
    fn sequences_under(svg: &str) -> Vec<(String, String)> {
        svg.split("<title>")
            .skip(1)
            .filter_map(|title| title.split("</title>").next())
            .filter_map(|title| {
                let (name, length) = title.split_once(", ")?;
                let length = length.strip_suffix(" bp")?;
                Some((name.to_string(), length.to_string()))
            })
            .collect()
    }

    /// A bedGraph with no place is drawn across every sequence it names,
    /// end to end in the order chromosomes are counted, each as long as its
    /// furthest row, named under it in place of a ruler. It was refused for
    /// having no region.
    #[test]
    fn coverage_across_a_genome_lays_sequences_end_to_end_in_chromosome_order() {
        let depth = "chr10\t0\t500\t3\nchr2\t0\t1000\t4\nchr2\t1000\t1500\t9\n\
                     chr1\t200\t900\t5\nchrX\t0\t100\t1\n";
        let held = [("d.bg", depth)];
        let svg = drawn_from("d.bg", &held).unwrap();
        assert_eq!(
            sequences_under(&svg),
            [
                ("chr1".to_string(), "900".to_string()),
                ("chr2".to_string(), "1,500".to_string()),
                ("chr10".to_string(), "500".to_string()),
                ("chrX".to_string(), "100".to_string()),
            ]
        );
        assert!(svg.contains("genome:1-3000"), "{svg}");
        assert!(svg.contains(">d</text>"), "the track keeps its file's name");
        for unit in [" kb</text>", " bp</text>"] {
            assert!(!svg.contains(unit), "a ruler of positions: {svg}");
        }
        // Painted where each sequence is laid: chr2's last row from 900 on
        // the shared axis, after chr1's 900 bases.
        let built = built_from("d.bg", &held, Theme::light(), None);
        assert_eq!(built.along, None, "a genome is no window to move along");
        // samtools depth, and windows, are laid out the same way.
        let positions = "chr2\t5\t7\nchr1\t10\t4\n";
        let windows = "chr2\t0\t2000\t0.5\nchr1\t0\t3000\t-0.25\n";
        let held = [("d.depth", positions), ("w.bg", windows)];
        let svg = drawn_from("--coverage d.depth --format depth --windows w.bg", &held).unwrap();
        assert_eq!(
            sequences_under(&svg),
            [
                ("chr1".to_string(), "3,000".to_string()),
                ("chr2".to_string(), "2,000".to_string()),
            ]
        );
    }

    /// A sequence one file names and another does not is drawn as missing in
    /// the second, a gap in its line, rather than as a depth of nought, which
    /// would say every sequence a sample's file left out was lost. A gap
    /// between rows on a sequence the file does name stays nought.
    #[test]
    fn a_sequence_a_file_does_not_name_is_missing_not_zero() {
        // The windows name chr2, which the depth leaves out between chr1 and
        // chr3, and the windows draw blocks, so every line is the depth's.
        let depth = "chr1\t0\t1000\t20\nchr3\t0\t1000\t30\n";
        let windows = "chr2\t0\t1000\t1\n";
        let held = [("d.bg", depth), ("w.bg", windows)];
        let svg = drawn_from("d.bg --style line --windows w.bg", &held).unwrap();
        let lines = crate::track::polylines(&svg);
        assert_eq!(lines.len(), 2, "one line either side of chr2: {svg}");
        // A gap inside chr1 is nought, and the line runs down to it.
        let held = [
            ("d.bg", "chr1\t0\t100\t20\nchr1\t900\t1000\t20\n"),
            ("w.bg", windows),
        ];
        let svg = drawn_from("d.bg --style line --windows w.bg", &held).unwrap();
        assert_eq!(crate::track::polylines(&svg).len(), 1, "{svg}");
    }

    /// Sequences whose lengths sum past what a u64 holds lay every track
    /// with a saturating sum, as the genome's lengths are laid: a scan,
    /// windows, a segment table and a bigWig's windows each panicked on a
    /// plain sum, where a depth over the same rows was drawn.
    #[test]
    fn sequences_longer_together_than_a_u64_are_drawn_not_a_panic() {
        let far = u64::MAX - 615;
        let windows = format!("chr1\t0\t{far}\t1\nchr2\t0\t1000\t2\n");
        let segments =
            format!("chromosome\tstart\tend\tlog2\nchr1\t0\t{far}\t0\nchr2\t0\t1000\t1\n");
        let scan = format!("CHR\tBP\tP\nchr1\t{far}\t0.5\nchr2\t1000\t0.01\n");
        let held = [
            ("w.bg", windows.as_str()),
            ("t.cns", segments.as_str()),
            ("gwas.assoc", scan.as_str()),
        ];
        for line in [
            "w.bg",
            "--windows w.bg",
            "t.cns --ploidy 2",
            "gwas.assoc",
            "w.bg --windows w.bg t.cns --ploidy 2",
        ] {
            let svg = drawn_from(line, &held).unwrap_or_else(|error| panic!("{line}: {error}"));
            // chr2 is a sliver at the end of the axis, too thin to be named.
            let under = sequences_under(&svg);
            assert_eq!(
                under.first().map(|(name, _)| name.as_str()),
                Some("chr1"),
                "{line}"
            );
        }
        // A bigWig's windows after a sequence that long: chr3 starts where
        // the axis has run out, and its windows are laid there.
        let mut held = binaries();
        held.insert("far.bg", format!("chr1\t0\t{}\t1\n", u64::MAX - 100));
        held_figure(&mut held, "--windows signal.bw --windows far.bg").unwrap();
    }

    /// A bigWig with no place is drawn across every sequence its index
    /// names, each as long as the index says rather than as far as its
    /// values reach, painted from the zoom level the whole figure wants and
    /// read under the names `--rename` gives. A sequence another file names
    /// and the bigWig does not is a gap in its line.
    #[test]
    fn a_bigwig_alone_is_drawn_across_the_sequences_its_index_names() {
        let mut held = binaries();
        // The bigWig's own values, reaching as far as its index says each
        // sequence runs, so the text lays out the genome the bigWig does.
        let sizes = "chr1\t999\t1000\t0\nchr2\t499\t500\t0\nchr3\t59\t60\t0\n";
        held.insert("sizes.bg", sizes);
        held.insert("w.bg", "chr4\t0\t1000\t1\n");
        for (binary, text) in [
            (
                "signal.bw --windows sizes.bg",
                "signal.bedgraph --windows sizes.bg",
            ),
            (
                "--windows signal.bw --windows sizes.bg",
                "--windows signal.bedgraph --windows sizes.bg",
            ),
        ] {
            let drawn = held_figure(&mut held, binary).unwrap();
            assert_eq!(drawn, held_figure(&mut held, text).unwrap(), "{binary}");
            assert_eq!(
                sequences_under(&drawn),
                [
                    ("chr1".to_string(), "1,000".to_string()),
                    ("chr2".to_string(), "500".to_string()),
                    ("chr3".to_string(), "60".to_string()),
                ],
                "{binary}"
            );
        }
        // Alone, each sequence is as long as the index says, which is past
        // the furthest value on chr1 and chr3.
        let alone = held_figure(&mut held, "signal.bw --style line").unwrap();
        assert_eq!(
            sequences_under(&alone),
            [
                ("chr1".to_string(), "1,000".to_string()),
                ("chr2".to_string(), "500".to_string()),
                ("chr3".to_string(), "60".to_string()),
            ]
        );
        assert_eq!(crate::track::polylines(&alone).len(), 1, "{alone}");
        // chr3, renamed past chr4, which only the windows name: the line
        // breaks over chr4 and goes on over chrX, the bigWig's chr3.
        let svg = held_figure(
            &mut held,
            "signal.bw --style line --windows w.bg --rename chr3=chrX",
        )
        .unwrap();
        assert_eq!(
            sequences_under(&svg),
            [
                ("chr1".to_string(), "1,000".to_string()),
                ("chr2".to_string(), "500".to_string()),
                ("chr4".to_string(), "1,000".to_string()),
                ("chrX".to_string(), "60".to_string()),
            ]
        );
        assert_eq!(
            crate::track::polylines(&svg).len(),
            2,
            "one line either side of chr4: {svg}"
        );
        // Across 301,560 bases a pixel holds two of the finest level's bins
        // of 119, so that level is read, as it is over a window as long. One
        // of its bins made 99 where no base under it is past 7, and the
        // scale goes past 7.5 when the level is read and only then.
        let mut bent = STORED_BW.to_vec();
        let most = 920 + 20;
        assert_eq!(bent[most..most + 4], 7.0f32.to_le_bytes());
        bent[most..most + 4].copy_from_slice(&99.0f32.to_le_bytes());
        held.insert("most.bw", bent);
        held.insert("far.bg", "chr4\t0\t300000\t1\n");
        let seven = |svg: &str| svg.contains(">7.5</text>");
        let few = held_figure(&mut held, "most.bw --windows sizes.bg").unwrap();
        assert!(seven(&few), "a figure of few bases read the level");
        let many = held_figure(&mut held, "most.bw --windows far.bg").unwrap();
        assert!(
            !seven(&many),
            "a figure of many bases did not read the level"
        );
        // Read from a level of means, the scale still reaches the most any
        // base holds, as the bedGraph's does: a bin's mean is less.
        held.insert("steps.bw", STEPS_BW);
        held.insert(
            "steps.bedgraph",
            include_str!("../read/fixtures/steps.bedgraph"),
        );
        held.insert("farther.bg", "chr9\t0\t6000000\t1\n");
        let scale = |svg: &str| -> Vec<String> {
            svg.split("<text")
                .filter_map(|text| text.split_once('>'))
                .filter_map(|(_, rest)| rest.split_once("</text>"))
                .map(|(label, _)| label.to_string())
                .filter(|label| label.parse::<f64>().is_ok())
                .collect()
        };
        let line = |file: &str| format!("--coverage {file} --aggregate mean --windows farther.bg");
        let zoomed = held_figure(&mut held, &line("steps.bw")).unwrap();
        let written = held_figure(&mut held, &line("steps.bedgraph")).unwrap();
        assert_eq!(scale(&zoomed), scale(&written));
        assert!(scale(&zoomed).contains(&"30".to_string()), "{zoomed}");
    }

    /// Windows and segments across an assembly of many scaffolds are laid
    /// out looking each scaffold up by its name, as a depth is: a walk along
    /// every name before it for each of them took 14 seconds for windows and
    /// 20 for segments at 100,000 scaffolds, where the depth took a tenth of
    /// one. Timed against the depth, which is linear, so a machine that is
    /// slow, or busy, is slow at both.
    #[test]
    fn many_scaffolds_are_laid_out_in_linear_time() {
        let count = 50_000;
        let mut rows = String::new();
        let mut segments = String::from("chromosome\tstart\tend\tcn\n");
        for n in 0..count {
            rows.push_str(&format!("s{n}\t0\t1000\t{}\n", n % 5));
            segments.push_str(&format!("s{n}\t0\t1000\t{}\n", n % 5));
        }
        let held = [("d.bg", rows.as_str()), ("t.cns", segments.as_str())];
        let timed = |line: &str| {
            let started = std::time::Instant::now();
            drawn_from(line, &held).unwrap();
            started.elapsed()
        };
        let depth = timed("d.bg");
        for line in ["--windows d.bg", "t.cns --ploidy 2"] {
            let taken = timed(line);
            assert!(
                taken < depth * 10 + std::time::Duration::from_millis(100),
                "{line} in {taken:?}, the depth in {depth:?}"
            );
        }
    }

    /// A segment table with no place is drawn across every sequence it
    /// calls, at the ploidy it was given, for the one sample asked for, and
    /// no riser joins the level one sequence ends on to the next one's.
    #[test]
    fn copy_number_across_a_genome_needs_its_ploidy_and_one_sample() {
        let seg = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n\
                   T1\t2\t1\t3000\t9\t1.0\nT1\t1\t1\t2000\t9\t-1.0\n\
                   T2\t1\t1\t2000\t9\t0.0\nT2\t2\t1\t3000\t9\t0.0\n";
        let held = [("t.seg", seg)];
        let said = drawn_from("t.seg --ploidy 2", &held)
            .unwrap_err()
            .to_string();
        assert!(said.contains("holds T1, T2, and --sample says"), "{said}");
        let svg = drawn_from("t.seg --ploidy 2 --sample T1", &held).unwrap();
        assert_eq!(
            sequences_under(&svg),
            [
                ("1".to_string(), "2,000".to_string()),
                ("2".to_string(), "3,000".to_string()),
            ]
        );
        // Each segment where it is on its own sequence, and none of them a
        // riser from one copy at the end of 1 to four at the start of 2.
        assert!(svg.contains("<title>1:1 to 2,000, 1 copies"), "{svg}");
        assert!(svg.contains("<title>2:1 to 3,000, 4 copies"), "{svg}");
        let vertical = svg
            .split("<line ")
            .skip(1)
            .filter(|line| {
                let at = |name: &str| {
                    line.split(&format!("{name}=\""))
                        .nth(1)
                        .and_then(|rest| rest.split('"').next())
                };
                at("x1") == at("x2") && at("y1") != at("y2")
            })
            .count();
        assert_eq!(vertical, 0, "a riser across the join: {svg}");
        // A table of no calls at all has nothing to draw.
        let none = "chromosome\tstart\tend\tcn\n1\t0\t100\tNA\n";
        let said = drawn_from("--copy-number n.cns --ploidy 2", &[("n.cns", none)])
            .unwrap_err()
            .to_string();
        assert!(said.contains("called segments"), "{said}");
    }

    /// mosdepth's depth in windows is named `.regions.bed.gz`, a bedGraph by
    /// another name, and alone it is drawn across the genome; a `.bed` of
    /// genes alone is still refused for want of a place, as it always was,
    /// and beside a scan it is the file named.
    #[test]
    fn a_mosdepth_regions_file_alone_is_a_genome_wide_depth() {
        let regions = "chr1\t0\t500\t31.20\nchr1\t500\t1000\t29.80\nchr2\t0\t400\t30.05\n";
        let held = [("s.regions.bed", regions)];
        let svg = drawn_from("s.regions.bed", &held).unwrap();
        assert_eq!(sequences_under(&svg).len(), 2, "{svg}");
        assert!(svg.contains("genome:1-1400"), "{svg}");
    }

    #[test]
    fn a_bed_of_genes_alone_is_still_refused_for_want_of_a_place() {
        let genes = "chr1\t100\t900\tgeneA\t0\t+\nchr2\t50\t300\tgeneB\t0\t-\n";
        let alone = drawn_from("genes.bed", &[("genes.bed", genes)]).unwrap_err();
        assert!(matches!(
            alone,
            BuildError::Placeless(crate::cli::args::ArgError::NoRegion)
        ));
        assert_eq!(
            alone.to_string(),
            crate::cli::args::ArgError::NoRegion.to_string()
        );
        let table = "CHR\tBP\tP\n1\t100\t0.5\n";
        let beside = drawn_from(
            "gwas.assoc genes.bed",
            &[("genes.bed", genes), ("gwas.assoc", table)],
        )
        .unwrap_err()
        .to_string();
        assert!(
            beside.starts_with("--features genes.bed is drawn over a place"),
            "{beside}"
        );
    }

    /// Files that name one sequence are that sequence drawn whole, as long as
    /// they reach, with a ruler in its own bases: the figure the place
    /// written out draws, which a page can move along. A scan of one
    /// chromosome is drawn so too.
    #[test]
    fn one_sequence_alone_is_that_sequence_drawn_whole() {
        let depth = "NC_1\t0\t5000\t12\nNC_1\t5000\t9000\t30\n";
        let held = [("a.bg", depth), ("b.bg", "NC_1\t0\t9500\t20\n")];
        let alone = drawn_from("a.bg b.bg --same-scale", &held).unwrap();
        let written = drawn_from("NC_1:1-9,500 a.bg b.bg --same-scale", &held).unwrap();
        assert_eq!(alone, written);
        assert_eq!(
            built_from("a.bg b.bg", &held, Theme::light(), None)
                .along
                .map(|region| region.to_string())
                .as_deref(),
            Some("NC_1:1-9500")
        );
        let table = "CHR\tBP\tP\n3\t100\t0.5\n3\t4000\t1e-9\n";
        let held = [("gwas.assoc", table)];
        assert_eq!(
            drawn_from("gwas.assoc", &held).unwrap(),
            drawn_from("3:1-4,000 gwas.assoc", &held).unwrap()
        );
    }

    /// A scan and a depth whose files name a chromosome differently are laid
    /// on one sequence where `--rename` says they are one, named as the
    /// figure names it, and a shade there falls where both are drawn.
    #[test]
    fn a_scan_and_a_depth_share_one_genome() {
        let table = "CHR\tBP\tP\n1\t1000\t0.01\n1\t900000\t0.2\n2\t5000\t1e-9\n";
        let depth = "chr1\t0\t1000000\t30\nchr2\t0\t600000\t25\n";
        let held = [("gwas.assoc", table), ("d.bg", depth)];
        let apart = drawn_from("gwas.assoc d.bg", &held).unwrap();
        assert_eq!(sequences_under(&apart).len(), 4, "{apart}");
        let built = built_from(
            "gwas.assoc d.bg --rename 1=chr1 --rename 2=chr2 --shade chr2:1-100,000=start",
            &held,
            Theme::light(),
            None,
        );
        let svg = built.figure.to_svg();
        assert_eq!(
            sequences_under(&svg),
            [
                ("chr1".to_string(), "1,000,000".to_string()),
                ("chr2".to_string(), "600,000".to_string()),
            ]
        );
        let shades = built.figure.shades();
        assert_eq!((shades[0].start(), shades[0].end()), (1_000_000, 1_100_000));
        assert!(
            svg.contains("<title>start, chr2:1-100,000</title>"),
            "{svg}"
        );
        // Named by the table's own name for it, it is found too.
        let built = built_from(
            "gwas.assoc d.bg --rename 1=chr1 --rename 2=chr2 --shade 2:1-100,000",
            &held,
            Theme::light(),
            None,
        );
        assert_eq!(built.figure.shades()[0].start(), 1_000_000);
        // And the depth's names are read under the figure's as the table's are.
        let svg = drawn_from("gwas.assoc d.bg --rename chr1=1 --rename chr2=2", &held).unwrap();
        assert_eq!(
            sequences_under(&svg),
            [
                ("1".to_string(), "1,000,000".to_string()),
                ("2".to_string(), "600,000".to_string()),
            ]
        );
    }

    /// A shade on a depth across the genome is laid through the offsets the
    /// sequences are laid at, and held to the sequence it is on.
    #[test]
    fn a_shade_on_a_genome_wide_depth_falls_on_its_own_sequence() {
        let depth = "1\t0\t2000\t30\n2\t0\t3000\t25\n3\t0\t1000\t20\n";
        let held = [("d.bg", depth)];
        let built = built_from(
            "d.bg --shade 2:1,001-5,000=loss",
            &held,
            Theme::light(),
            None,
        );
        let shades = built.figure.shades();
        assert_eq!((shades[0].start(), shades[0].end()), (3_000, 5_000));
        assert!(built
            .figure
            .to_svg()
            .contains("<title>loss, 2:1,001-5,000</title>"));
        let said = drawn_from("d.bg --shade 4:1-10", &held)
            .unwrap_err()
            .to_string();
        assert!(said.contains("is on 4"), "{said}");
    }

    /// Several places are one panel each, the same tracks over each, the
    /// title over them all and the key once under them. A track with nothing
    /// in one place is a band there that says so, where it refused the whole
    /// figure; alone, a place with nothing on a track is refused as before.
    #[test]
    fn several_places_are_panels_of_the_same_tracks() {
        let genome = format!(">chr1\n{}\n", "ACGT".repeat(12_500));
        let calls = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chr1\t10500\t.\tA\tG\t.\t.\t.\n";
        let held = [
            ("genes.gff3", GENES),
            ("calls.vcf", calls),
            ("ref.fa", genome.as_str()),
        ];
        let svg = drawn_from(
            "rpoB katG genes.gff3 calls.vcf --sequence ref.fa --title Both",
            &held,
        )
        .unwrap();
        assert_eq!(
            svg.matches("<svg").count(),
            4,
            "the sheet, two panels and the key: {svg}"
        );
        assert_eq!(svg.matches(">Both</text>").count(), 1);
        assert!(svg.contains(">rpoB</text>") && svg.contains(">katG</text>"));
        assert_eq!(svg.matches(">no variants here</text>").count(), 1);
        assert_eq!(svg.matches(">T</text>").count(), 1, "the key once: {svg}");
        let alone = drawn_from("katG genes.gff3 calls.vcf", &held).unwrap_err();
        assert!(alone.to_string().contains("no variants"), "{alone}");
    }

    /// What is piped in is read once, for every panel of a sheet.
    /// A gene of two transcripts, as Ensembl writes one.
    const ISOFORMS: &str = "##gff-version 3
7\t.\tgene\t1001\t9000\t.\t+\t.\tID=g1;Name=GENE1
7\t.\tmRNA\t1001\t9000\t.\t+\t.\tID=t1;Parent=g1;Name=GENE1-201
7\t.\texon\t1001\t1500\t.\t+\t.\tParent=t1
7\t.\texon\t4001\t4500\t.\t+\t.\tParent=t1
7\t.\texon\t8001\t9000\t.\t+\t.\tParent=t1
7\t.\tCDS\t1201\t8600\t.\t+\t0\tParent=t1
7\t.\tmRNA\t1001\t9000\t.\t+\t.\tID=t2;Parent=g1;Name=GENE1-202
7\t.\texon\t1001\t1500\t.\t+\t.\tParent=t2
7\t.\texon\t8001\t9000\t.\t+\t.\tParent=t2
";

    #[test]
    fn a_gene_is_drawn_once_with_its_exons_or_each_transcript_on_request() {
        let held = [("genes.gff3", ISOFORMS)];
        let once = drawn_from("7:1-10,000 genes.gff3", &held).unwrap();
        assert!(
            once.contains(
                "<title>GENE1, 1,001 to 9,000, forward, 3 exons from 2 transcripts</title>"
            ),
            "{once}"
        );
        assert!(!once.contains("GENE1-20"), "{once}");
        let each = drawn_from("7:1-10,000 genes.gff3 --isoforms", &held).unwrap();
        assert!(each.contains("<title>GENE1-201 (GENE1), 1,001 to 9,000, forward, 3 exons</title>"));
        assert!(each.contains("<title>GENE1-202 (GENE1), 1,001 to 9,000, forward, 2 exons</title>"));
        // Placed by the gene's name, the transcripts are what is drawn there.
        let placed = drawn_from("GENE1 genes.gff3 --isoforms", &held).unwrap();
        assert_eq!(placed.matches("<title>GENE1-20").count(), 2, "{placed}");
        // Anywhere but after an annotation it is refused, by the track it was
        // written after.
        let args: Vec<String> = "7:1-10,000 --coverage genes.gff3 --isoforms"
            .split_whitespace()
            .map(String::from)
            .collect();
        assert!(parse(&args).is_err());
    }

    /// The depth of two samples, one reaching 97 and one 48, over the first
    /// kilobase, and the second sample alone over the next.
    const DEEP: &str = "c1\t0\t500\t97\nc1\t500\t1000\t60\nc1\t1000\t2000\t30\n";
    const SHALLOW: &str = "c1\t0\t1000\t48\nc1\t1000\t2000\t40\n";

    #[test]
    fn several_depths_share_one_scale_when_asked() {
        let held = [("deep.bedgraph", DEEP), ("shallow.bedgraph", SHALLOW)];
        let line = "c1:1-1,000 deep.bedgraph shallow.bedgraph";
        let own = drawn_from(line, &held).unwrap();
        let same = drawn_from(&format!("{line} --same-scale"), &held).unwrap();
        // As if each had been pinned by hand to the ceiling the deeper one
        // rounds to.
        let pinned = drawn_from(
            "c1:1-1,000 deep.bedgraph --max 100 shallow.bedgraph --max 100",
            &held,
        )
        .unwrap();
        assert_eq!(same, pinned);
        assert_ne!(same, own);
        // A pinned track keeps its pin, and the other has its own scale.
        let one_pinned = "c1:1-1,000 deep.bedgraph --max 200 shallow.bedgraph";
        assert_eq!(
            drawn_from(&format!("{one_pinned} --same-scale"), &held).unwrap(),
            drawn_from(one_pinned, &held).unwrap()
        );
    }

    #[test]
    fn the_panels_of_several_places_share_one_scale_when_asked() {
        let held = [("deep.bedgraph", DEEP)];
        // Over the first place the depth reaches 97, over the second 30: on
        // one scale the second panel reads off the first one's ceiling.
        let line = "c1:1-1,000 c1:1,001-2,000 deep.bedgraph";
        let same = drawn_from(&format!("{line} --same-scale"), &held).unwrap();
        let pinned = drawn_from(&format!("{line} --max 100"), &held).unwrap();
        assert_eq!(same, pinned);
        assert_ne!(same, drawn_from(line, &held).unwrap());
    }

    #[test]
    fn a_maximum_is_a_number_above_nought_for_a_track_with_a_ceiling() {
        let refused = |line: &str| {
            let args: Vec<String> = line.split_whitespace().map(String::from).collect();
            parse(&args).err().map(|error| error.to_string())
        };
        assert!(refused("c1:1-10 --features g.gff3 --max 5")
            .unwrap()
            .contains("--max means nothing to a features track"));
        for bad in ["0", "-4", "NaN", "ten"] {
            let said = refused(&format!("c1:1-10 --coverage d.bedgraph --max {bad}")).unwrap();
            assert!(said.contains("a number above nought"), "{said}");
        }
        assert_eq!(refused("c1:1-10 --manhattan g.assoc --max 12"), None);
        assert_eq!(refused("c1:1-10 --recombination m.txt --max 50"), None);
        assert_eq!(refused("c1:1-10 d.bedgraph --same-scale"), None);
        for track in [
            "--windows w.bg",
            "--matrix m.tsv",
            "--heatmap h.tsv",
            "--pairs p.ld",
        ] {
            assert_eq!(
                refused(&format!("c1:1-10 {track} --max 2")),
                None,
                "{track}"
            );
        }
    }

    /// Windows either side of a line of nought, with the furthest at 1.4.
    const SCORES: &str = "c1\t0\t2000\t1.4\nc1\t2000\t4000\t-0.6\nc1\t4000\t6000\t0.3\n\
                          c1\t6000\t8000\t-1.1\nc1\t8000\t10000\t0.9\n";

    /// The tick labels a document writes, in the order it writes them.
    fn texts(svg: &str) -> Vec<String> {
        svg.split("<text")
            .skip(1)
            .filter_map(|piece| piece.split_once('>'))
            .filter_map(|(_, rest)| rest.split_once("</text>"))
            .map(|(text, _)| text.to_string())
            .collect()
    }

    /// The top is the one number a band symmetric about its line has free,
    /// so `--max 2` reads -2 to 2, and it is what the library draws when its
    /// reach is pinned by hand.
    #[test]
    fn windows_take_a_maximum_either_side_of_their_line() {
        let held = [("w.bg", SCORES)];
        let own = drawn_from("c1:1-10,000 --windows w.bg", &held).unwrap();
        let pinned = drawn_from("c1:1-10,000 --windows w.bg --max 2", &held).unwrap();
        let ticks = texts(&pinned);
        for tick in ["2", "0", "-2"] {
            assert!(ticks.iter().any(|text| text == tick), "{tick}: {ticks:?}");
        }
        let own_ticks = texts(&own);
        for tick in ["1.5", "0", "-1.5"] {
            assert!(
                own_ticks.iter().any(|text| text == tick),
                "{tick}: {own_ticks:?}"
            );
        }
        assert!(!own_ticks.iter().any(|text| text == "-2"), "{own}");
        let windows =
            read::signal::windows(SCORES, &Region::new("c1", 0, 10_000).unwrap()).unwrap();
        let by_hand = Plot::over(Region::new("c1", 0, 10_000).unwrap())
            .add_track(
                WindowTrack::new(windows)
                    .style(WindowStyle::Steps)
                    .extent(2.0)
                    .label("w"),
            )
            .to_svg();
        assert_eq!(pinned, by_hand);
    }

    /// A heatmap read against each sample's usual depth runs from a full loss
    /// at nought to a gain at the pin, and a plain one from nought to it.
    #[test]
    fn a_heatmap_pinned_reads_its_colours_off_the_maximum() {
        let table = "chrom\tstart\tend\tS1\tS2\nc1\t0\t5000\t40\t60\nc1\t5000\t10000\t0\t131.4\n";
        let held = [("h.tsv", table)];
        let relative = drawn_from("c1:1-10,000 --heatmap h.tsv --relative --max 3", &held).unwrap();
        let ticks = texts(&relative);
        for end in ["0×", "1×", "3×"] {
            assert!(ticks.iter().any(|text| text == end), "{end}: {ticks:?}");
        }
        let own = texts(&drawn_from("c1:1-10,000 --heatmap h.tsv --relative", &held).unwrap());
        assert!(!own.iter().any(|text| text == "3×"), "{own:?}");
        let plain = texts(&drawn_from("c1:1-10,000 --heatmap h.tsv --max 200", &held).unwrap());
        assert!(plain.iter().any(|text| text == "200"), "{plain:?}");
        assert!(!plain.iter().any(|text| text == "131.4"), "{plain:?}");
        // A matrix of sites takes it the same way.
        let sites = "sample\t100\t900\nS1\t0.2\t0.4\nS2\t0.1\t0.3\n";
        let matrix =
            texts(&drawn_from("c1:1-1,000 --matrix m.tsv --max 1", &[("m.tsv", sites)]).unwrap());
        assert!(matrix.iter().any(|text| text == "1"), "{matrix:?}");
        assert!(!matrix.iter().any(|text| text == "0.4"), "{matrix:?}");
    }

    /// Linkage is read against an r² of one until told otherwise, and a
    /// pin is told otherwise: weak linkage read on a scale of its own.
    #[test]
    fn pairs_take_their_full_colour_from_max_even_for_linkage() {
        let ld = " CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n\
                   c1 100 a c1 400 b 0.31\n\
                   c1 100 a c1 700 c 0.12\n\
                   c1 400 b c1 700 c 0.44\n";
        let held = [("linkage.ld", ld)];
        let own = texts(&drawn_from("c1:1-1,000 linkage.ld", &held).unwrap());
        assert!(own.iter().any(|text| text == "1"), "{own:?}");
        let pinned = texts(&drawn_from("c1:1-1,000 linkage.ld --max 0.5", &held).unwrap());
        assert!(pinned.iter().any(|text| text == "0.5"), "{pinned:?}");
        assert!(!pinned.iter().any(|text| text == "1"), "{pinned:?}");
    }

    /// Two places in one region of chr1, and a gene on another sequence.
    const SHADED_GENES: &str = "##gff-version 3\n\
        c1\t.\tgene\t2001\t3000\t.\t+\t.\tID=g1;Name=GENE1\n\
        c1\t.\tCDS\t2001\t2997\t.\t+\t0\tID=cds1;Parent=g1;gene=GENE1\n\
        c2\t.\tgene\t501\t900\t.\t-\t.\tID=g2;Name=GENE2\n";

    /// A flat depth over two sequences.
    const SHADED_DEPTH: &str = "c1\t0\t10000\t30\nc2\t0\t10000\t20\n";

    /// The left and right of every wash rectangle under `said`.
    fn wash_of(svg: &str, said: &str) -> Vec<(f64, f64)> {
        let title = format!("<g><title>{said}</title>");
        let Some(at) = svg.find(&title) else {
            return Vec::new();
        };
        let group = &svg[at + title.len()..];
        let group = &group[..group.find("</g>").unwrap()];
        group
            .split("<rect")
            .skip(1)
            .map(|rect| {
                let number = |name: &str| -> f64 {
                    let key = format!(" {name}=\"");
                    let at = rect.find(&key).unwrap() + key.len();
                    rect[at..].split('"').next().unwrap().parse().unwrap()
                };
                (number("x"), number("x") + number("width"))
            })
            .collect()
    }

    #[test]
    fn a_shade_on_another_sequence_is_refused_and_names_the_places() {
        let held = [("d.bg", SHADED_DEPTH)];
        let said = drawn_from("c1:1-1,000 d.bg --shade c2:100-200", &held)
            .unwrap_err()
            .to_string();
        assert_eq!(
            said,
            "--shade c2:100-200 is on c2, and the figure is drawn over c1:1-1,000"
        );
    }

    #[test]
    fn a_shade_outside_the_window_is_a_note_and_draws_the_plain_figure() {
        let held = [("d.bg", SHADED_DEPTH)];
        let (svg, notes) = drawn_noting("c1:1-4,000 d.bg --shade c1:5,001-6,000=later", &held);
        let (plain, _) = drawn_noting("c1:1-4,000 d.bg", &held);
        assert_eq!(svg.unwrap(), plain.unwrap());
        assert_eq!(
            notes,
            ["--shade c1:5,001-6,000=later is outside c1:1-4,000, so it is not drawn"]
        );
        // In view, it is drawn and nothing is said of it.
        let (svg, notes) = drawn_noting("c1:1-6,000 d.bg --shade c1:5,001-6,000=later", &held);
        assert!(notes.is_empty(), "{notes:?}");
        assert_eq!(wash_of(&svg.unwrap(), "later, 5,001 to 6,000").len(), 1);
    }

    /// A gene is shaded from its own start to its own end: the figure placed
    /// on it has a margin either side, and the shade does not.
    #[test]
    fn a_gene_shade_covers_the_gene_and_not_its_margin() {
        let held = [("genes.gff3", SHADED_GENES), ("d.bg", SHADED_DEPTH)];
        let built = built_from(
            "GENE1 genes.gff3 d.bg --shade GENE1",
            &held,
            Theme::light(),
            None,
        );
        let region = built.figure.region().clone();
        assert!(region.start() < 2_000 && region.end() > 3_000, "{region}");
        let shades = built.figure.shades();
        assert_eq!(shades.len(), 1);
        // The gene row and its CDS are one place, at the gene's own ends.
        assert_eq!((shades[0].start(), shades[0].end()), (2_000, 3_000));
        let svg = built.figure.to_svg();
        let columns = wash_of(&svg, "2,001 to 3,000");
        assert_eq!(columns.len(), 1, "{svg}");
    }

    #[test]
    fn a_gene_to_shade_needs_an_annotation_and_a_name_it_has() {
        let said = drawn_from("c1:1-10,000 d.bg --shade GENE1", &[("d.bg", SHADED_DEPTH)])
            .unwrap_err()
            .to_string();
        assert!(
            said.contains("it has none; add the GFF3, GTF or BED"),
            "{said}"
        );
        let held = [("genes.gff3", SHADED_GENES), ("d.bg", SHADED_DEPTH)];
        let said = drawn_from("c1:1-10,000 genes.gff3 d.bg --shade GENE3", &held)
            .unwrap_err()
            .to_string();
        assert!(said.contains("did you mean GENE1 or GENE2?"), "{said}");
        // A gene the annotation puts on another sequence is elsewhere.
        let said = drawn_from("c1:1-10,000 genes.gff3 d.bg --shade GENE2", &held)
            .unwrap_err()
            .to_string();
        assert_eq!(
            said,
            "--shade GENE2 is on c2, and the figure is drawn over c1:1-10,000"
        );
    }

    /// An alignment is its own place, named after its file, and a span with
    /// no sequence is on its columns whatever the file is called.
    #[test]
    fn a_bare_span_shades_an_alignment_s_columns() {
        let held = [("aln.fa", ALIGNMENT)];
        let built = built_from("--msa aln.fa --shade 3-5", &held, Theme::light(), None);
        let shades = built.figure.shades();
        assert_eq!((shades[0].start(), shades[0].end()), (2, 5));
        assert_eq!(wash_of(&built.figure.to_svg(), "3 to 5").len(), 1);
    }

    /// A place in decimal years is read in years, and so is a shade on it.
    #[test]
    fn a_shade_on_a_continuous_time_is_rescaled_with_the_place() {
        let skyline = "year\tmedian\tlower\tupper\n2010.25\t100\t50\t200\n\
                       2012.5\t400\t300\t600\n2015.75\t900\t700\t1200\n";
        let held = [("sky.tsv", skyline)];
        for line in [
            "--phylodynamics sky.tsv --shade 2012-2013",
            "year:2010-2016 --phylodynamics sky.tsv --shade year:2012-2013",
        ] {
            let built = built_from(line, &held, Theme::light(), None);
            let shades = built.figure.shades();
            assert_eq!(
                (shades[0].start(), shades[0].end()),
                (2_012_000, 2_013_001),
                "{line}"
            );
            assert_eq!(
                wash_of(&built.figure.to_svg(), "2012 to 2013").len(),
                1,
                "{line}"
            );
        }
    }

    #[test]
    fn a_sheet_shades_each_panel_on_its_own_sequence() {
        let held = [("d.bg", SHADED_DEPTH)];
        let (svg, notes) = drawn_noting("c1:1-1,000 c2:1-1,000 d.bg --shade c2:100-200", &held);
        let svg = svg.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        assert_eq!(svg.matches("<title>100 to 200</title>").count(), 1, "{svg}");
        // The wash is in the second panel, after the first one's locus.
        let first = svg.find("c1:1-1000").unwrap();
        assert!(svg.find("<title>100 to 200</title>").unwrap() > first);
        let said = drawn_from("c1:1-1,000 c2:1-1,000 d.bg --shade c3:1-10", &held)
            .unwrap_err()
            .to_string();
        assert_eq!(
            said,
            "--shade c3:1-10 is on c3, and the figure is drawn over c1:1-1,000 and c2:1-1,000"
        );
        // On the right sequences and in neither window is a note.
        let (svg, notes) = drawn_noting("c1:1-1,000 c2:1-1,000 d.bg --shade c2:5,000", &held);
        assert!(svg.is_ok());
        assert_eq!(
            notes,
            ["--shade c2:5,000 is outside c2:1-1,000, so it is not drawn"]
        );
    }

    #[test]
    fn a_genome_wide_scan_shades_through_the_offsets() {
        let table = "CHR\tSNP\tBP\tP\n1\ta\t1000\t0.01\n1\tb\t900000\t0.2\n\
                     2\tc\t5000\t1e-9\n2\td\t1000000\t0.5\n";
        let held = [("gwas.assoc", table)];
        let built = built_from(
            "gwas.assoc --shade 2:1-1,000,000=peak",
            &held,
            Theme::light(),
            None,
        );
        let shades = built.figure.shades();
        // Sequence 2 starts where 1 ends, at its furthest marker, and each is
        // as long as its furthest marker reaches.
        let genome = crate::Genome::new([("1", 900_000u64), ("2", 1_000_000)]);
        let start = genome.offset("2").unwrap();
        assert_eq!(shades[0].start(), start);
        assert_eq!(shades[0].end(), start + 1_000_000);
        let svg = built.figure.to_svg();
        assert!(svg.contains("<title>peak, 2:1-1,000,000</title>"), "{svg}");
        // Past the furthest marker is a note, and a bare span is refused.
        let (_, notes) = drawn_noting("gwas.assoc --shade 2:2,000,001-2,000,100", &held);
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert!(
            notes[0].contains("past the furthest any file reaches on 2"),
            "{notes:?}"
        );
        let said = drawn_from("gwas.assoc --shade 1-1,000", &held)
            .unwrap_err()
            .to_string();
        assert!(said.contains("a figure across the whole genome"), "{said}");
        // A stretch past where a sequence ends in the figure stops there,
        // at its furthest marker, rather than running on into the next.
        let built = built_from(
            "gwas.assoc --shade 1:800,001-2,000,000",
            &held,
            Theme::light(),
            None,
        );
        let shades = built.figure.shades();
        let first = genome.offset("1").unwrap();
        assert_eq!(
            (shades[0].start(), shades[0].end()),
            (first + 800_000, first + 900_000)
        );
        // A sequence --rename names otherwise is found by that name too, and
        // shaded where the table's own name for it is laid.
        let built = built_from(
            "gwas.assoc --rename 2=chrB --shade chrB:1-500,000=p",
            &held,
            Theme::light(),
            None,
        );
        let shades = built.figure.shades();
        assert_eq!(
            (shades[0].start(), shades[0].end()),
            (start, start + 500_000)
        );
        assert!(built
            .figure
            .to_svg()
            .contains("<title>p, chrB:1-500,000</title>"));
    }

    /// A scan across the whole genome reads no annotation, and one beside it
    /// needs a place, so a gene is answered with the one form that shades
    /// it: its span on its sequence.
    #[test]
    fn a_gene_shade_on_a_genome_wide_scan_is_answered_with_its_span() {
        let table = "CHR\tSNP\tBP\tP\n1\ta\t1000\t0.01\n2\tc\t5000\t1e-9\n";
        let said = drawn_from("gwas.assoc --shade GENE1", &[("gwas.assoc", table)])
            .unwrap_err()
            .to_string();
        assert_eq!(
            said,
            "--shade GENE1: a figure across the whole genome reads no annotation to look a \
             gene up in; write the gene's span on its sequence, as 7:1,001-2,000"
        );
    }

    #[test]
    fn a_renamed_sequence_is_shaded_by_either_name() {
        let held = [("d.bg", "1\t0\t10000\t30\n")];
        for shade in ["NC_1:101-200", "1:101-200"] {
            let built = built_from(
                &format!("NC_1:1-1,000 d.bg --rename 1=NC_1 --shade {shade}"),
                &held,
                Theme::light(),
                None,
            );
            let shades = built.figure.shades();
            assert_eq!((shades[0].start(), shades[0].end()), (100, 200), "{shade}");
        }
    }

    #[test]
    fn a_figure_with_nothing_on_the_coordinates_refuses_a_shade() {
        let held = [("tree.nwk", "((a:1,b:1):1,c:2);")];
        let said = drawn_from("tree.nwk --shade c1:1-10", &held)
            .unwrap_err()
            .to_string();
        assert!(
            said.starts_with("--shade c1:1-10: nothing in this figure is laid on the coordinates"),
            "{said}"
        );
    }

    /// Two genomes' neighbourhoods, column one naming the genome, and what
    /// joins them.
    const SHADED_LOCI: &str = "H37Rv\t100\t900\tesxA\nH37Rv\t1000\t1800\tesxB\n\
                               BCG\t150\t950\tesxA_b\nBCG\t1100\t1900\tesxB_b\n";
    const SHADED_HITS: &str = "esxA\tesxA_b\t98.5\nesxB\tesxB_b\t97.0\n";

    /// A row of --loci is drawn on the figure's axis whatever genome its
    /// first column names, and its gene is shaded where it is drawn: held to
    /// the sequence of that name, it was refused as on H37Rv in the figure
    /// that drew it.
    #[test]
    fn a_gene_of_the_loci_is_shaded_where_its_row_draws_it() {
        let held = [("loci.bed", SHADED_LOCI), ("hits.tsv", SHADED_HITS)];
        for (line, span) in [
            (
                "locus:1-2,500 --loci loci.bed --links hits.tsv --shade esxA=A",
                (100, 900),
            ),
            (
                "H37Rv:1-2,500 --loci loci.bed --links hits.tsv --shade esxA_b",
                (150, 950),
            ),
        ] {
            let built = built_from(line, &held, Theme::light(), None);
            let shades = built.figure.shades();
            assert_eq!(shades.len(), 1, "{line}");
            assert_eq!((shades[0].start(), shades[0].end()), span, "{line}");
        }
        // Outside the window is a note, as for any gene.
        let (svg, notes) = drawn_noting(
            "locus:1,001-2,500 --loci loci.bed --links hits.tsv --shade esxA",
            &held,
        );
        assert!(svg.is_ok());
        assert_eq!(
            notes,
            ["--shade esxA is outside locus:1,001-2,500, so it is not drawn"]
        );
    }

    /// Every gene to shade is looked up in one read of each annotation, for
    /// every panel of a sheet: each shade read it again for each panel.
    #[test]
    fn gene_shades_read_the_annotation_once_for_every_panel() {
        let reads = |line: &str| {
            let args: Vec<String> = line.split_whitespace().map(String::from).collect();
            let Request::Draw(invocation) = parse(&args).unwrap() else {
                unreachable!("a figure")
            };
            let mut read = 0;
            build(&invocation, |source: &Source| {
                let Source::Path(path) = source else {
                    unreachable!("every source is a file")
                };
                match path.to_string_lossy().as_ref() {
                    "genes.gff3" => {
                        read += 1;
                        Ok(SHADED_GENES.to_string())
                    }
                    "d.bg" => Ok(SHADED_DEPTH.to_string()),
                    other => Err(io::Error::new(io::ErrorKind::NotFound, other.to_string())),
                }
            })
            .unwrap_or_else(|error| panic!("{line}: {error}"));
            read
        };
        let places = "c1:1-4,000 c1:4,001-8,000 c2:1-1,000 genes.gff3 d.bg";
        let plain = reads(places);
        assert_eq!(
            reads(&format!("{places} --shade GENE1 --shade GENE2")),
            plain + 1
        );
        // A figure of one place reads it once more, however many genes.
        let plain = reads("c1:1-4,000 genes.gff3 d.bg");
        assert_eq!(
            reads("c1:1-4,000 genes.gff3 d.bg --shade GENE1 --shade gene1=again"),
            plain + 1
        );
    }

    /// A skyline over whole years and one over fractions of them.
    const SHADED_YEARS: &str = "year\tmedian\tlower\tupper\n2010\t100\t50\t200\n\
                                2012\t400\t300\t600\n2015\t900\t700\t1200\n";
    const SHADED_FRACTIONS: &str = "year\tmedian\tlower\tupper\n2010.25\t100\t50\t200\n\
                                    2012.5\t400\t300\t600\n2015.75\t900\t700\t1200\n";

    /// A shade on a time is said as its ruler and its tooltips say a time,
    /// in whole years as in fractions: whole years were grouped as bases are,
    /// `2,012 to 2,013`, under a ruler reading 2012.
    #[test]
    fn a_shade_on_a_time_is_said_in_years_whole_or_not() {
        for table in [SHADED_YEARS, SHADED_FRACTIONS] {
            let svg = drawn_from(
                "--phylodynamics sky.tsv --shade 2012-2013=a",
                &[("sky.tsv", table)],
            )
            .unwrap();
            assert!(svg.contains("<title>a, 2012 to 2013</title>"), "{svg}");
            assert!(svg.contains(" Shaded: a, 2012 to 2013.</desc>"), "{svg}");
        }
        // Beside a depth the ruler counts bases and groups them, and the
        // shade is said as that ruler says a span.
        let svg = drawn_from(
            "c1:1-3,000 d.bg --phylodynamics sky.tsv --shade 2012-2013=a",
            &[("sky.tsv", SHADED_YEARS), ("d.bg", SHADED_DEPTH)],
        )
        .unwrap();
        assert!(svg.contains(">2,000</text>"), "{svg}");
        assert!(svg.contains("<title>a, 2,012 to 2,013</title>"), "{svg}");
    }

    /// The window a shade on a time is outside of, or refused against, is
    /// said in the table's own units: it was said in the thousandths a table
    /// with fractions is drawn at, as `year:2,010,251-2,015,751`.
    #[test]
    fn a_window_of_time_is_named_in_its_own_units() {
        for (line, table, window) in [
            (
                "--phylodynamics sky.tsv",
                SHADED_FRACTIONS,
                "year:2010.25-2015.75",
            ),
            (
                "year:2010-2016 --phylodynamics sky.tsv",
                SHADED_FRACTIONS,
                "year:2010-2016",
            ),
            ("--phylodynamics sky.tsv", SHADED_YEARS, "year:2010-2015"),
            (
                "year:2010-2016 --phylodynamics sky.tsv",
                SHADED_YEARS,
                "year:2010-2016",
            ),
        ] {
            let held = [("sky.tsv", table)];
            let (svg, notes) = drawn_noting(&format!("{line} --shade 2030-2040=a"), &held);
            assert!(svg.is_ok(), "{line}");
            assert_eq!(
                notes,
                [format!(
                    "--shade 2030-2040=a is outside {window}, so it is not drawn"
                )],
                "{line}"
            );
            let said = drawn_from(&format!("{line} --shade foo:1-2=a"), &held)
                .unwrap_err()
                .to_string();
            assert_eq!(
                said,
                format!("--shade foo:1-2=a is on foo, and the figure is drawn over {window}"),
                "{line}"
            );
        }
    }

    #[test]
    fn a_sheet_reads_what_is_piped_in_once() {
        let svg = drawn_piped("rpoB katG --features -", GENES, &[]).unwrap();
        assert_eq!(svg.matches("<svg").count(), 3, "{svg}");
        assert!(svg.contains(">rpoB</text>") && svg.contains(">katG</text>"));
    }

    fn locus_of(svg: &str) -> String {
        svg.split("<title id=\"karyon-title\">")
            .nth(1)
            .and_then(|rest| rest.split('<').next())
            .unwrap_or_default()
            .to_string()
    }

    /// A gene is where its annotation says, with a tenth of its length either
    /// side, and the figure is titled with its name. Its gene row and its CDS
    /// are one place.
    #[test]
    fn a_figure_placed_by_a_gene_is_the_gene_and_a_margin_around_it() {
        let svg = drawn_from("rpoB genes.gff3", &[("genes.gff3", GENES)]).unwrap();
        assert_eq!(locus_of(&svg), "rpoB, chr1:9801-12200");
        // Asked for in another case, it is titled the way the file spells it.
        let svg = drawn_from("RPOB genes.gff3", &[("genes.gff3", GENES)]).unwrap();
        assert!(locus_of(&svg).starts_with("rpoB,"), "{}", locus_of(&svg));
        // By its locus tag, too.
        let svg = drawn_from("SYN_1 genes.gff3", &[("genes.gff3", GENES)]).unwrap();
        assert!(
            locus_of(&svg).ends_with("chr1:9801-12200"),
            "{}",
            locus_of(&svg)
        );
        // And from a GTF or a BED.
        let gtf = "chr1\t.\tgene\t10001\t12000\t.\t+\t.\tgene_id \"G1\"; gene_name \"rpoB\";\n";
        let svg = drawn_from("rpoB genes.gtf", &[("genes.gtf", gtf)]).unwrap();
        assert_eq!(locus_of(&svg), "rpoB, chr1:9801-12200");
        let bed = "chr1\t10000\t12000\trpoB\t0\t+\n";
        let svg = drawn_from("rpoB genes.bed", &[("genes.bed", bed)]).unwrap();
        assert_eq!(locus_of(&svg), "rpoB, chr1:9801-12200");
    }

    #[test]
    fn a_gene_at_two_places_or_at_none_is_refused_saying_so() {
        let twice = format!("{GENES}chr2\t.\tgene\t5\t500\t.\t+\t.\tID=gene-C;Name=rpoB\n");
        let error = drawn_from("rpoB genes.gff3", &[("genes.gff3", &twice)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "rpoB is named at 2 places, chr1:10001-12000 and chr2:5-500; give the one you mean as a region"
        );
        let error = drawn_from("rpoC genes.gff3", &[("genes.gff3", GENES)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "no gene and no sequence is called rpoC; did you mean rpoB? genes.gff3 names the sequence chr1."
        );
    }

    /// A name no file has is most often one file's name for what another
    /// calls otherwise, as PLINK writes 1 for a chromosome a FASTA names in
    /// full, so the refusal says what each file does call its sequences. It
    /// said to add an annotation, and named only sequences a header gave.
    #[test]
    fn a_name_no_file_has_is_answered_with_what_each_file_calls_its_sequences() {
        let scan = "CHR SNP BP A1 P\n1 rs1 150 A 0.5\n1 rs2 4800 A 1e-9\n";
        let error = drawn_from("NC_1 gwas.assoc", &[("gwas.assoc", scan)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "no gene and no sequence is called NC_1. gwas.assoc names the sequence 1. If 1 is \
             NC_1, add --rename 1=NC_1. To find a gene by name, add the annotation that names \
             it, as --features genes.gff3."
        );
        // Files that name the same sequences are named together.
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   chr1\t50\t.\tA\tG\t.\t.\t.\n";
        let error = drawn_from(
            "chrZ genes.gff3 calls.vcf gwas.assoc",
            &[
                ("genes.gff3", GENES),
                ("calls.vcf", vcf),
                ("gwas.assoc", scan),
            ],
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "no gene and no sequence is called chrZ. genes.gff3 and calls.vcf name the \
             sequence chr1. gwas.assoc names the sequence 1."
        );
        // BOLT-LMM writes the chromosome second, behind the variant's name,
        // which the empty window's message took for the sequence.
        let bolt = "SNP\tCHR\tBP\tA1FREQ\tP_BOLT_LMM\nrs1\t1\t100\t0.1\t1e-3\n\
                    rs2\t2\t300\t0.1\t0.5\n";
        let error = drawn_from(
            "chr9:1-1000 --manhattan bolt.stats",
            &[("bolt.stats", bolt)],
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .ends_with("though the file holds 2 on 1 and 2"),
            "{error}"
        );
    }

    /// A sequence's name is the sequence, whole, as long as a file says it is,
    /// or as far as any file reaches on it.
    #[test]
    fn a_figure_placed_by_a_sequence_is_the_whole_of_it() {
        let svg = drawn_from("chr1 genes.gff3", &[("genes.gff3", GENES)]).unwrap();
        assert_eq!(locus_of(&svg), "chr1:1-50000");
        let fasta = format!(">ctg7 a contig\n{}\n", "ACGT".repeat(50));
        let svg = drawn_from("ctg7 ref.fa", &[("ref.fa", &fasta)]).unwrap();
        assert_eq!(locus_of(&svg), "ctg7:1-200");
        let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chrX,length=9000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrX\t50\t.\tA\tG\t.\t.\t.\n";
        let svg = drawn_from("chrX calls.vcf", &[("calls.vcf", vcf)]).unwrap();
        assert_eq!(locus_of(&svg), "chrX:1-9000");
        // Nothing says how long it is, so the scan says how far it reaches.
        let scan = "CHR SNP BP A1 P\n1 rs1 150 A 0.5\n1 rs2 4800 A 1e-9\n";
        let svg = drawn_from("1 gwas.assoc", &[("gwas.assoc", scan)]).unwrap();
        assert_eq!(locus_of(&svg), "1:1-4800");
    }

    /// A stem's height is the AF the VCF gives, and the axis says so; with
    /// no AF there is no axis, and nothing to title.
    #[test]
    fn the_calls_axis_is_titled_by_what_it_measures() {
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let with_af = format!("{header}chr1\t50\t.\tA\tG\t.\t.\tAF=0.4\n");
        let svg = drawn_from("chr1:1-100 calls.vcf", &[("calls.vcf", &with_af)]).unwrap();
        assert!(svg.contains(">AF</text>"), "{svg}");
        let without = format!("{header}chr1\t50\t.\tA\tG\t.\t.\t.\n");
        let svg = drawn_from("chr1:1-100 calls.vcf", &[("calls.vcf", &without)]).unwrap();
        assert!(!svg.contains(">AF</text>"), "{svg}");
    }

    /// The names in a figure's gutter, top to bottom.
    fn labels_of(svg: &str) -> Vec<String> {
        svg.split("font-weight=\"600\"")
            .skip(1)
            .filter_map(|piece| {
                piece
                    .split('>')
                    .nth(1)?
                    .split('<')
                    .next()
                    .map(str::to_string)
            })
            .collect()
    }

    /// Every track given no label is called after its file.
    #[test]
    fn a_track_is_called_after_its_file_unless_it_is_labelled() {
        let svg = drawn_from(
            "chr1:9,000-13,000 data/genes.gff3 --features genes.gff3 --label annotation",
            &[("data/genes.gff3", GENES), ("genes.gff3", GENES)],
        )
        .unwrap();
        let text = labels_of(&svg);
        assert_eq!(text, ["genes", "annotation"]);
        let unnamed = TrackSpec::new(Kind::Coverage, Some(Source::Path("/dev/fd/63".into())));
        assert_eq!(default_label(&unnamed), None, "a pipe the shell named");
        let piped = TrackSpec::new(Kind::Coverage, Some(Source::Stdin));
        assert_eq!(default_label(&piped), None);
        let compressed = TrackSpec::new(Kind::Variants, Some(Source::Path("calls.vcf.gz".into())));
        assert_eq!(default_label(&compressed).as_deref(), Some("calls"));
        // A BAM drawn as a line is its depth, which a name saying reads hid.
        let depth = TrackSpec::new(Kind::Coverage, Some(Source::Path("reads.bam".into())));
        assert_eq!(default_label(&depth).as_deref(), Some("reads depth"));
        let reads = TrackSpec::new(Kind::Pileup, Some(Source::Path("reads.bam".into())));
        assert_eq!(default_label(&reads).as_deref(), Some("reads"));
        // A tree is plain to see, and its file's name was a stray word beside it.
        let tree = TrackSpec::new(Kind::Tree, Some(Source::Path("tree.nwk".into())));
        assert_eq!(default_label(&tree), None);
    }

    /// modkit writes bedMethyl as `.bed`, and a `.bed` of four columns whose
    /// last is a number is a bedGraph.
    #[test]
    fn a_file_whose_name_hides_its_format_is_drawn_as_what_it_holds() {
        let methyl =
            "chr1\t100\t101\tm\t30\t+\t100\t101\t255,0,0\t30\t80.00\t24\t6\t0\t0\t0\t0\t0\n";
        let svg = drawn_from("chr1:1-500 calls.bed", &[("calls.bed", methyl)]).unwrap();
        assert!(svg.contains("methylation site"), "not drawn as methylation");
        let signal = "chr1\t0\t100\t5\nchr1\t100\t200\t9\n";
        assert_eq!(refine(Kind::Features, signal), Some(Kind::Coverage));
        assert_eq!(refine(Kind::Features, "chr1\t0\t100\tgeneA\n"), None);
    }

    /// A pileup given no reference of its own reads against the one the
    /// figure draws. It drew every read as agreeing, with the reference right
    /// above it, and the variant nowhere.
    #[test]
    fn a_pileup_reads_against_the_reference_the_figure_draws() {
        let reference = format!(">chr1\n{}\n", "ACGT".repeat(25));
        let sam = "r1\t0\tchr1\t11\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\n";
        let held = [("ref.fa", reference.as_str()), ("reads.sam", sam)];
        let implied = drawn_from("chr1:1-40 --sequence ref.fa --pileup reads.sam", &held).unwrap();
        let named = drawn_from(
            "chr1:1-40 --sequence ref.fa --pileup reads.sam --with-sequence ref.fa",
            &held,
        )
        .unwrap();
        let without = drawn_from("chr1:1-40 --pileup reads.sam", &held).unwrap();
        assert_eq!(implied, named);
        assert_ne!(
            implied.replace("ref</text>", ""),
            without,
            "no mismatch was painted"
        );
    }

    /// The key to a tree's colours goes under the ruler, and is left out when
    /// asked.
    #[test]
    fn the_key_to_the_colours_is_drawn_under_the_figure() {
        let tree = "((a:1,b:1):1,(c:1,d:1):1);";
        let sheet = "sample\tlineage\na\tL4\nb\tL4\nc\tL2\nd\tL1\n";
        let held = [("t.nwk", tree), ("s.tsv", sheet)];
        let svg = drawn_from("--tree t.nwk --traits s.tsv", &held).unwrap();
        for level in ["L4", "L2", "L1"] {
            assert!(
                svg.contains(&format!("lineage: {level}")),
                "no key for {level}"
            );
        }
        let bare = drawn_from("--tree t.nwk --traits s.tsv --no-legend", &held).unwrap();
        assert!(!bare.contains("lineage: L4"));
    }

    /// The colour of the swatch beside a key's label.
    fn key_colour(svg: &str, label: &str) -> String {
        let at = svg
            .find(&format!(">{label}</text>"))
            .unwrap_or_else(|| panic!("no key for {label}"));
        let before = &svg[..svg[..at].rfind("<text").expect("the label's element")];
        let fill = before.rfind("fill=\"#").expect("a swatch before the label");
        before[fill + 6..fill + 13].to_string()
    }

    /// Branches coloured by lineage beside a strip of countries painted a
    /// lineage and a country one colour; each column keeps its own stretch
    /// of the palette, whichever of them is drawn.
    #[test]
    fn a_tree_s_branches_and_its_strips_keep_to_their_own_colours() {
        let tree = "((a:1,b:1):1,(c:1,d:1):1);";
        let sheet =
            "sample\tlineage\tcountry\na\tL4\tKenya\nb\tL4\tSpain\nc\tL2\tKenya\nd\tL1\tVietnam\n";
        let held = [("t.nwk", tree), ("s.tsv", sheet)];
        let both = drawn_from("--tree t.nwk --traits s.tsv --color-by lineage", &held).unwrap();
        assert_ne!(
            key_colour(&both, "lineage: L4"),
            key_colour(&both, "country: Kenya")
        );
        let alone = drawn_from(
            "--tree t.nwk --traits s.tsv --columns country --color-by lineage",
            &held,
        )
        .unwrap();
        assert_eq!(
            key_colour(&alone, "country: Kenya"),
            key_colour(&both, "country: Kenya"),
            "a country's colour hung on which columns were drawn"
        );
        assert_ne!(
            key_colour(&alone, "lineage: L4"),
            key_colour(&alone, "country: Kenya")
        );
    }

    /// One sheet, one row name per row of every track that draws strips from
    /// a sheet, and the files each of those tracks reads.
    const STRIPPED: &[(&str, &str)] = &[
        ("s.tsv", "sample\tlineage\nA\tL1\nB\tL2\nC\tL1\nD\tL3\n"),
        (
            "cohort.vcf",
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n\
             chr\t50\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t0/1\n",
        ),
        (
            "m.tsv",
            "sample\t100\t200\nA\t1\t0\nB\t0\t1\nC\t1\t1\nD\t0\t0\n",
        ),
        (
            "h.tsv",
            "chrom\tstart\tend\tA\tB\tC\tD\n\
             chr\t0\t100\t1\t2\t3\t4\nchr\t100\t200\t2\t2\t3\t1\nchr\t200\t300\t1\t1\t1\t1\n",
        ),
        (
            "aln.fa",
            ">A\nACGTACGT\n>B\nACGTTCGT\n>C\nACGAACGT\n>D\nTCGTACGA\n",
        ),
        ("t.nwk", "((A:1,B:1):1,(C:1,D:1):1);"),
        (
            "c.gff",
            "SEQUENCE\tGUBBINS\tCDS\t101\t200\t0.000\t.\t0\tnode=\"N1\";taxa=\"A B\";\n",
        ),
        (
            "d.tsv",
            "A\tmd5\t300\tPfam\tPF1\tKinase\t10\t100\t1e-5\tT\t01-01-2026\n\
             B\tmd5\t300\tPfam\tPF1\tKinase\t20\t110\t1e-5\tT\t01-01-2026\n\
             C\tmd5\t300\tPfam\tPF2\tBinding\t150\t250\t1e-5\tT\t01-01-2026\n\
             D\tmd5\t300\tPfam\tPF1\tKinase\t30\t120\t1e-5\tT\t01-01-2026\n",
        ),
        (
            "l.bed",
            "A\t100\t900\tgA\t0\t+\nB\t150\t950\tgB\t0\t+\n\
             C\t120\t920\tgC\t0\t+\nD\t110\t910\tgD\t0\t+\n",
        ),
        ("links.tsv", "gA\tgB\t98\ngB\tgC\t97\ngC\tgD\t96\n"),
    ];

    /// Every track that draws strips from a sheet paints a level `--colors`
    /// chose in that colour, in its cells and in its key, and a level it did
    /// not choose in the palette's colour, the same in all of them.
    ///
    /// The colours reach a track through the sheet it is given, so a track
    /// added to the ones that take `--traits` and wired past `strip` would
    /// draw the palette here, and one left out of this list fails the count.
    #[test]
    fn colors_paint_a_level_alike_in_every_track_that_draws_strips() {
        // In the order `Kind::ALL` lists them.
        let lines = [
            (Kind::Genotypes, "chr:1-300 --genotypes cohort.vcf"),
            (Kind::Tree, "--tree t.nwk"),
            (Kind::Msa, "--msa aln.fa"),
            (Kind::Snps, "--snps aln.fa"),
            (Kind::Matrix, "chr:1-300 --matrix m.tsv"),
            (Kind::Heatmap, "chr:1-300 --heatmap h.tsv"),
            (
                Kind::Clades,
                "SEQUENCE:1-300 --clades c.gff --with-tree t.nwk",
            ),
            (Kind::Loci, "locus:1-1000 --loci l.bed --links links.tsv"),
            (Kind::Domains, "protein:1-300 --domains d.tsv"),
        ];
        let striped: Vec<Kind> = Kind::ALL
            .iter()
            .copied()
            .filter(|kind| kind.takes_traits())
            .collect();
        assert_eq!(
            striped,
            lines.map(|(kind, _)| kind),
            "a track that takes --traits is not drawn here"
        );
        let mut unchosen: Vec<String> = Vec::new();
        for (kind, line) in lines {
            let line = format!("{line} --traits s.tsv --colors lineage=L1:#aa0000,L2:#00aa00");
            let svg = drawn_from(&line, STRIPPED).unwrap_or_else(|error| panic!("{line}: {error}"));
            for row in ["A", "C"] {
                assert_eq!(
                    painted(&svg, &format!("{row}; lineage L1")),
                    ["#aa0000"],
                    "{kind:?} paints L1 its own colour"
                );
            }
            assert_eq!(painted(&svg, "B; lineage L2"), ["#00aa00"], "{kind:?}");
            assert_eq!(key_colour(&svg, "lineage: L1"), "#aa0000", "{kind:?}");
            unchosen.extend(painted(&svg, "D; lineage L3"));
        }
        // L3 was given no colour and takes the palette's, one in every track.
        assert_eq!(unchosen.len(), lines.len(), "{unchosen:?}");
        assert!(
            unchosen.iter().all(|colour| *colour == unchosen[0]),
            "{unchosen:?}"
        );
        assert_eq!(unchosen[0], crate::Theme::light().color(2));
    }

    /// The colours are the figure's: a tree and a matrix with sheets of their
    /// own paint a country alike, and the branches coloured by a column drawn
    /// as no strip take its colours too.
    #[test]
    fn colors_reach_every_sheet_of_a_figure_and_branches_with_no_strip() {
        let tree = "((a:1,b:1):1,(c:1,d:1):1);";
        let sheet =
            "sample\tlineage\tcountry\na\tL4\tKenya\nb\tL4\tSpain\nc\tL2\tKenya\nd\tL1\tPeru\n";
        // The matrix's sheet has no lineage, and says nothing of Peru.
        let other = "sample\tcountry\na\tKenya\nb\tSpain\nc\tKenya\nd\tSpain\n";
        let matrix = "sample\t100\t200\na\t1\t0\nb\t0\t1\nc\t1\t1\nd\t0\t0\n";
        let held = [
            ("t.nwk", tree),
            ("s.tsv", sheet),
            ("o.tsv", other),
            ("m.tsv", matrix),
        ];
        let colors = "--colors country=Kenya:#aa0000,Peru:#0000aa --colors lineage=L4:#00aa00";
        let svg = drawn_from(
            &format!(
                "chr:1-300 --tree t.nwk --traits s.tsv --matrix m.tsv --traits o.tsv {colors}"
            ),
            &held,
        )
        .unwrap();
        assert_eq!(painted(&svg, "a; country Kenya"), ["#aa0000", "#aa0000"]);
        assert_eq!(painted(&svg, "d; country Peru"), ["#0000aa"]);
        assert_eq!(painted(&svg, "a; lineage L4"), ["#00aa00"]);
        // Branches by a column the tree draws no strip of.
        let svg = drawn_from(
            &format!("--tree t.nwk --traits s.tsv --columns lineage --color-by country {colors}"),
            &held,
        )
        .unwrap();
        assert!(svg.contains("stroke=\"#aa0000\""), "no branch is Kenya's");
        assert!(svg.contains("stroke=\"#0000aa\""), "no branch is Peru's");
        assert_eq!(key_colour(&svg, "country: Kenya"), "#aa0000");
    }

    /// Seven countries in six colours are drawn as shapes, and the tree says
    /// so with the way out; given colours of their own they stay a strip, and
    /// the tree has nothing to say.
    #[test]
    fn colors_keep_a_column_of_seven_values_a_strip() {
        let names = ["a", "b", "c", "d", "e", "f", "g"];
        let tree = format!("({});", names.map(|name| format!("{name}:1")).join(","));
        let rows: String = names
            .iter()
            .enumerate()
            .map(|(i, name)| format!("{name}\tC{i}\n"))
            .collect();
        let sheet = format!("sample\tcountry\n{rows}");
        let held = [("t.nwk", tree.as_str()), ("s.tsv", sheet.as_str())];
        let (svg, notes) = drawn_noting("--tree t.nwk --traits s.tsv", &held);
        let said = "country: 7 values for 6 colours, so each is a shape as well; \
                    --colors gives them colours of their own";
        assert_eq!(notes, [format!("--tree t.nwk: {said}")]);
        assert!(svg.unwrap().contains(said), "the figure does not say it");
        let chosen: Vec<String> = (0..7).map(|i| format!("C{i}:#{i}{i}0000")).collect();
        let line = format!(
            "--tree t.nwk --traits s.tsv --colors country={}",
            chosen.join(",")
        );
        let (svg, notes) = drawn_noting(&line, &held);
        assert!(notes.is_empty(), "{notes:?}");
        let svg = svg.unwrap();
        for (i, name) in names.iter().enumerate() {
            assert_eq!(
                painted(&svg, &format!("{name}; country C{i}")),
                [format!("#{i}{i}0000")]
            );
        }
        assert!(!svg.contains("<polygon"), "a country is drawn as a shape");
    }

    /// A `--colors` that would paint nothing is refused, and says what the
    /// sheets hold instead.
    #[test]
    fn colors_that_would_paint_nothing_are_refused() {
        let tree = "((a:1,b:1):1,(c:1,d:1):1);";
        let sheet = "sample\tlineage\tyear\na\tL4\t2019\nb\tL4\t2020\nc\tL2\t2021\nd\tL1\t2020\n";
        let held = [("t.nwk", tree), ("s.tsv", sheet)];
        let refused = |colors: &str| {
            drawn_from(&format!("--tree t.nwk --traits s.tsv {colors}"), &held)
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            refused("--colors linage=L4:#aa0000"),
            "--colors names a column called linage, and s.tsv has none; it has lineage, year"
        );
        assert_eq!(
            refused("--colors lineage=L3:#aa0000"),
            "--colors names L3 in lineage, and no row of s.tsv holds it; lineage holds L1, L2, L4"
        );
        assert_eq!(
            refused("--colors year=2020:#aa0000"),
            "--colors names year, a column of numbers, which is drawn as a ramp; --colors \
             paints a column of words"
        );
        let undrawn = drawn_from(
            "--tree t.nwk --traits s.tsv --columns year --colors lineage=L4:#aa0000",
            &held,
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            undrawn,
            "--colors paints lineage, and no track draws it: name it in --columns, or colour \
             a tree's branches by it with --color-by lineage"
        );
        // Drawn as the branches' colour, it is drawn.
        assert!(drawn_from(
            "--tree t.nwk --traits s.tsv --columns year --color-by lineage \
             --colors lineage=L4:#aa0000",
            &held
        )
        .is_ok());
        // A column of nothing but gaps, and a sheet of a header alone, hold
        // no value to list, and say so rather than end at "holds".
        for (path, text) in [
            ("na.tsv", "sample\tlineage\na\tNA\nb\t.\nc\tNA\nd\tNA\n"),
            ("e.tsv", "sample\tlineage\n"),
        ] {
            let error = drawn_from(
                &format!("--tree t.nwk --traits {path} --colors lineage=L4:#aa0000"),
                &[("t.nwk", tree), (path, text)],
            )
            .unwrap_err()
            .to_string();
            assert_eq!(
                error,
                format!(
                    "--colors names L4 in lineage, and no row of {path} holds a value in \
                     lineage"
                )
            );
        }
    }

    /// The colours are checked over every sheet of the figure at once, and
    /// each rule over all of them: a column is drawn when a track whose own
    /// sheet has it draws it, and is a ramp when it is numbers in every sheet
    /// that has it.
    #[test]
    fn colors_are_checked_over_every_sheet_of_the_figure() {
        let tree = "((a:1,b:1):1,(c:1,d:1):1);";
        let matrix = "sample\t100\t200\na\t1\t0\nb\t0\t1\nc\t1\t1\nd\t0\t0\n";
        let sheet = "sample\tlineage\tgroup\na\tL4\tx\nb\tL4\ty\nc\tL2\tx\nd\tL1\ty\n";
        // The matrix's sheet has no lineage, so it draws every column it has
        // and no lineage among them.
        let other = "sample\tgroup\na\tx\nb\ty\nc\tx\nd\ty\n";
        // Lineage as numbers, a ramp beside the matrix.
        let numbers = "sample\tlineage\na\t4\nb\t4\nc\t2\nd\t1\n";
        let held = [
            ("t.nwk", tree),
            ("m.tsv", matrix),
            ("s.tsv", sheet),
            ("o.tsv", other),
            ("n.tsv", numbers),
        ];
        let undrawn = drawn_from(
            "chr:1-300 --tree t.nwk --traits s.tsv --columns group --matrix m.tsv \
             --traits o.tsv --colors lineage=L4:#aa0000",
            &held,
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            undrawn,
            "--colors paints lineage, and no track draws it: name it in --columns, or colour \
             a tree's branches by it with --color-by lineage"
        );
        // Words beside the tree and numbers beside the matrix: painted where
        // it is words.
        let svg = drawn_from(
            "chr:1-300 --tree t.nwk --traits s.tsv --matrix m.tsv --traits n.tsv \
             --colors lineage=L4:#aa0000",
            &held,
        )
        .unwrap();
        assert_eq!(painted(&svg, "a; lineage L4"), ["#aa0000"]);
    }

    /// Values a colour chose alike are named in the line under the tree,
    /// with their colour, and `--colors` is not offered as the way out: it is
    /// what joined them. The tree said "4 values for 6 colours".
    #[test]
    fn colors_chosen_alike_are_named_under_the_tree() {
        let held = [
            ("t.nwk", "((a:1,b:1):1,(c:1,d:1):1);"),
            ("s.tsv", "sample\tlineage\na\tL4\nb\tL2\nc\tL1\nd\tL3\n"),
        ];
        let (svg, notes) = drawn_noting(
            "--tree t.nwk --traits s.tsv --colors lineage=L1:#aa0000,L2:#aa0000",
            &held,
        );
        let said = "lineage: L1 and L2 are both #aa0000, so each is a shape as well";
        assert_eq!(notes, [format!("--tree t.nwk: {said}")]);
        assert!(svg.unwrap().contains(said), "the figure does not say it");
    }

    /// A value that holds a comma takes its colour, as a place written
    /// `Korea, Rep.` does; split at every comma, it was refused as if the
    /// flag were written wrong.
    #[test]
    fn colors_reach_a_value_that_holds_a_comma() {
        let mut files = Held::new();
        files.insert("t.nwk", "((a:1,b:1):1,(c:1,d:1):1);");
        files.insert(
            "s.tsv",
            "sample\tcountry\na\tKorea, Rep.\nb\tPeru\nc\tKorea, Rep.\nd\tPeru\n",
        );
        let line = [
            "--tree",
            "t.nwk",
            "--traits",
            "s.tsv",
            "--colors",
            "country=Korea, Rep.:#aa0000,Peru:#0000aa",
        ]
        .map(String::from);
        let Request::Draw(invocation) = parse(&line).unwrap() else {
            unreachable!("a figure")
        };
        let svg = build_files(&invocation, &mut files, |_, _| None).unwrap();
        assert_eq!(painted(&svg, "a; country Korea, Rep."), ["#aa0000"]);
        assert_eq!(painted(&svg, "b; country Peru"), ["#0000aa"]);
    }

    /// The values `--colors` names are looked up in a set of the column's
    /// values, not in a list of them: a column of a value per row took eight
    /// seconds at 120,000 rows, where drawing the figure took half a second.
    /// Timed against reading the sheet, which is linear, so a machine that
    /// is slow, or busy, is slow at both.
    #[test]
    fn colors_are_checked_against_a_column_of_many_values_in_linear_time() {
        let rows = 60_000;
        let mut sheet = String::from("sample\tbatch\n");
        for row in 0..rows {
            sheet.push_str(&format!("s{row}\tB{row}\n"));
        }
        let mut files = Held::new();
        files.insert("t.nwk", "((s0:1,s1:1):1,(s2:1,s3:1):1);");
        files.insert("s.tsv", sheet.as_str());
        let invocation = sheeted("--tree t.nwk --traits s.tsv --colors batch=B0:#aa0000");
        let started = std::time::Instant::now();
        let read = read::sheet::sheet(&sheet).unwrap();
        let reading = started.elapsed();
        assert_eq!(read.levels("batch").len(), rows);
        let started = std::time::Instant::now();
        colored_as_asked(&invocation, &mut files).unwrap();
        let checking = started.elapsed();
        assert!(
            checking < reading * 10 + std::time::Duration::from_millis(100),
            "checked in {checking:?}, read in {reading:?}"
        );
    }

    /// A phylogram draws its scale bar by default, and --no-scale-bar is how
    /// a figure goes without it.
    #[test]
    fn a_tree_has_its_scale_bar_unless_told_otherwise() {
        let held = [("t.nwk", "((a:0.1,b:0.1):0.05,(c:0.1,d:0.1):0.05);")];
        let bar = "branch length scale";
        assert!(drawn_from("--tree t.nwk", &held).unwrap().contains(bar));
        assert!(!drawn_from("--tree t.nwk --no-scale-bar", &held)
            .unwrap()
            .contains(bar));
    }

    /// Bases too narrow for their letters are blocks of colour, which the key
    /// names; with letters, or with nothing drawn, there is nothing to name.
    /// No simulated user could say which colour was which base.
    #[test]
    fn the_key_names_the_bases_while_they_are_blocks() {
        let fasta = format!(">chr1\n{}\n", "ACGT".repeat(500));
        let held = [("ref.fa", fasta.as_str())];
        let blocks = drawn_from("chr1:1-400 ref.fa", &held).unwrap();
        for base in ["A", "C", "G", "T"] {
            assert!(
                blocks.contains(&format!(">{base}</text>")),
                "no key for {base}"
            );
        }
        let quiet = drawn_from("chr1:1-400 ref.fa --no-legend", &held).unwrap();
        assert!(
            !quiet.contains(">A</text>"),
            "the key was drawn with --no-legend"
        );
        let whole = drawn_from("chr1:1-2000 ref.fa", &held).unwrap();
        assert!(whole.contains("zoom in to see bases"));
        assert!(
            !whole.contains(">A</text>"),
            "a key for bases nobody can see"
        );
    }

    /// A folded row is an internal node and carries nobody's metadata, so the
    /// strip beside it can only say what its tips agree on. Both directions,
    /// because saying nothing and picking one of two are the two ways to be
    /// wrong here.
    #[test]
    fn a_folded_clade_says_what_its_tips_agree_on_and_nothing_when_they_differ() {
        // (a, b) are both L4 and agree; (c, d) are L2 and L1 and do not.
        const TREE: &str = "((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);";
        const SHEET: &str = "name\tlineage\na\tL4\nb\tL4\nc\tL2\nd\tL1\n";
        let open = |source: &Source| -> io::Result<String> {
            Ok(match source {
                Source::Path(path) if path.to_string_lossy().ends_with(".nwk") => TREE.to_string(),
                _ => SHEET.to_string(),
            })
        };

        let invocation = sheeted("tree:1-1 --tree t.nwk --traits s.tsv --max-rows 2");
        let svg = build(&invocation, open).unwrap();

        assert!(
            svg.contains("a +1 more; lineage L4"),
            "a clade of two L4 tips is an L4 clade: {svg}"
        );
        assert!(
            svg.contains("c +1 more; lineage missing"),
            "a clade holding L2 and L1 is neither of them: {svg}"
        );
    }
    /// The whole road: an annotated Newick in, a clade marked out.
    ///
    /// What this pins is the road and not the clade rule. Marking only the
    /// branch a change happened on gives the same figure, because a colour is
    /// inherited down the tree and a tip with nothing of its own takes its
    /// ancestor's. The rule itself is pinned where it lives, in
    /// `a_change_is_carried_by_everything_below_where_it_happened`, and it is
    /// what anything asking the question directly gets back.
    #[test]
    fn asking_who_carries_a_change_marks_the_clade_that_does() {
        const TREE: &str = concat!(
            "((a[&muts=\"A1T\"]:0.1,b:0.1)[&muts=\"S:D614G\"]:0.1,",
            "(c[&muts=\"S:D614G\"]:0.1,d:0.1):0.1);"
        );
        let open = |_: &Source| -> io::Result<String> { Ok(TREE.to_string()) };

        let svg = build(
            &sheeted("tree:1-1 --tree t.nwk --mutations muts --carrying S:D614G"),
            open,
        )
        .unwrap();
        // A branch's tooltip names the change and not the tip, so each tip is
        // found by its label's own row and its branch read off at that height.
        let row_of = |tip: &str| -> String {
            let at = svg
                .find(&format!(">{tip}</text>"))
                .unwrap_or_else(|| panic!("{tip} is drawn: {svg}"));
            let y = svg[..at]
                .rmatch_indices("y=\"")
                .next()
                .map(|(start, _)| svg[start + 3..].split('"').next().unwrap_or("").to_string())
                .expect("the label has a y");
            // The label sits a third of a body below its own row, the body
            // being a point under the theme's.
            let row = y.parse::<f64>().unwrap() - (Theme::light().font_size - 1.0) * 0.35;
            let mark = format!("y1=\"{row}\"");
            let line = svg
                .split(&mark)
                .nth(1)
                .unwrap_or_else(|| panic!("no branch at {tip}'s row {row}: {svg}"));
            line.split("stroke=\"")
                .nth(1)
                .and_then(|piece| piece.split('"').next())
                .unwrap_or("")
                .to_string()
        };

        // a and b are under the branch it happened on, and c has it directly.
        // d is under neither.
        let marked = row_of("a");
        assert_ne!(
            marked,
            Theme::light().foreground,
            "a's branch is marked, not plain: {svg}"
        );
        assert_eq!(row_of("b"), marked, "b is under the same branch: {svg}");
        assert_eq!(row_of("c"), marked, "c has it of its own: {svg}");
        assert_eq!(
            row_of("d"),
            Theme::light().foreground,
            "d carries nothing and stays plain: {svg}"
        );

        // The tree the file holds is read with its annotations. It was read
        // without them for a long while, which meant --color-by and this could
        // never see anything a file said.
        let plain = build(&sheeted("tree:1-1 --tree t.nwk"), open).unwrap();
        assert!(
            !plain.contains("carries"),
            "nothing is marked unasked: {plain}"
        );

        // A change the tree has not got is refused against the ones it has,
        // rather than drawing a tree with nothing marked on it.
        let refused = build(
            &sheeted("tree:1-1 --tree t.nwk --mutations muts --carrying Z9Z"),
            open,
        )
        .unwrap_err()
        .to_string();
        assert!(refused.contains("no change called Z9Z"), "{refused}");
        assert!(
            refused.contains("S:D614G"),
            "and says what it has: {refused}"
        );

        // And a key with no changes under it is refused too, since a figure
        // with nothing marked looks like a figure of a tree that carries
        // nothing.
        let empty = build(
            &sheeted("tree:1-1 --tree t.nwk --mutations nothing --carrying S:D614G"),
            open,
        )
        .unwrap_err()
        .to_string();
        assert!(empty.contains("mutations under that key"), "{empty}");
    }

    /// `highlight_named` does nothing when it cannot find the name, which is a
    /// choice for a library call and a lie for a command line: the figure comes
    /// out with one clade fewer than was asked for and nothing says so.
    #[test]
    fn a_clade_that_is_not_there_is_refused_by_name() {
        const TREE: &str = "((a:1,b:1):1,(c:1,d:1):1);";
        let open = |_: &crate::cli::args::Source| Ok(TREE.to_string());
        let drawn = build(&sheeted("tree:1-1 --tree t.nwk --highlight a,c"), open);
        assert!(drawn.is_ok(), "two names that are there: {drawn:?}");

        let refused = build(&sheeted("tree:1-1 --tree t.nwk --highlight a,zzz"), open)
            .unwrap_err()
            .to_string();
        assert!(refused.contains("no clade called zzz"), "{refused}");
        assert!(refused.contains('a'), "and says what it has: {refused}");
    }

    /// A key no node carries colours no branch, and the tree came out exactly
    /// as it would have without the flag. Refused with the keys the tree does
    /// carry, the way a change it does not carry is, since a misspelt key is
    /// the usual way to get here.
    #[test]
    fn a_colour_key_the_tree_does_not_carry_is_refused_with_the_ones_it_does() {
        const TREE: &str = concat!(
            "((a[&lineage=\"L4\"]:0.1,b[&lineage=\"L4\"]:0.1)[&muts=\"C241T\"]:0.1,",
            "(c[&lineage=\"L2\"]:0.1,d:0.1):0.1);"
        );
        const SHEET: &str = "name\thost\na\tcattle\nb\thuman\nc\thuman\nd\tcattle\n";
        let open = |source: &Source| -> io::Result<String> {
            Ok(match source {
                Source::Path(path) if path.to_string_lossy().ends_with(".tsv") => SHEET,
                _ => TREE,
            }
            .to_string())
        };

        let refused = build(&sheeted("tree:1-1 --tree t.nwk --color-by linage"), open)
            .unwrap_err()
            .to_string();
        assert_eq!(
            refused,
            "--tree t.nwk has no annotation called linage; it has lineage, muts"
        );

        // A key only an internal node carries still colours, since a branch
        // takes its nearest annotated ancestor's value; and a column of the
        // sheet is on the tips by the time the key is looked for.
        for line in [
            "tree:1-1 --tree t.nwk --color-by lineage",
            "tree:1-1 --tree t.nwk --color-by muts",
            "tree:1-1 --tree t.nwk --traits s.tsv --color-by host",
        ] {
            assert!(build(&sheeted(line), open).is_ok(), "{line}");
        }

        // A tree that carries nothing at all says so rather than trailing off.
        let bare = |_: &Source| -> io::Result<String> { Ok("((a,b),(c,d));".to_string()) };
        let refused = build(&sheeted("tree:1-1 --tree t.nwk --color-by host"), bare)
            .unwrap_err()
            .to_string();
        assert_eq!(
            refused,
            "--tree t.nwk has no annotation called host; it has none"
        );
    }

    /// The parser has taken `--height` on a copy number track since the track
    /// arrived, and this arm never passed it on, so the band came out at its
    /// own height whatever was asked for, byte for byte the figure it would
    /// have been without the flag.
    #[test]
    fn a_copy_number_track_is_as_tall_as_it_is_asked_to_be() {
        let table = "chromosome\tstart\tend\tcn\nchr8\t0\t500\t3\nchr8\t500\t1000\t1\n";
        let tall = |px: u32| -> f64 {
            let line = format!("chr8:1-1000 --copy-number s.cns --ploidy 2 --height {px}");
            let svg = build(&sheeted(&line), |_| Ok(table.to_string())).unwrap();
            svg.split_once(" height=\"")
                .and_then(|(_, rest)| rest.split('"').next())
                .and_then(|number| number.parse().ok())
                .expect("the figure says how tall it is")
        };
        assert_eq!(tall(200) - tall(100), 100.0, "the band ignored --height");
    }

    /// A position the pileup could not call is skipped by the reader rather
    /// than drawn at nought per cent, and counted. The count stopped there: the
    /// band printed how many calls its floor held back and said nothing about
    /// the positions that never became calls at all.
    #[test]
    fn positions_nobody_could_call_are_counted_on_the_band() {
        let pileup = concat!(
            "chr1\t10\t11\tm\t40\t+\t10\t11\t0,0,0\t40\t95.00\t38\t2\t0\t0\t0\t0\t0\n",
            "chr1\t20\t21\tm\t0\t+\t20\t21\t0,0,0\t0\t0.00\t0\t0\t0\t0\t12\t0\t4\n",
            "chr1\t30\t31\tm\t0\t-\t30\t31\t0,0,0\t0\t0.00\t0\t0\t0\t0\t9\t0\t2\n",
            "chr1\t40\t41\tm\t3\t-\t40\t41\t0,0,0\t3\t66.67\t2\t1\t0\t0\t0\t0\t0\n",
        );
        let svg = build(&over("chr1:1-100", "--methylation", "calls.bed"), |_| {
            Ok(pileup.to_string())
        })
        .unwrap();
        // Beside the floor's own count, in the corner that already holds it.
        assert!(
            svg.contains(">1 under 5x, 2 with no coverage</text>"),
            "{svg}"
        );
    }

    const SIGNAL_BW: &[u8] = include_bytes!("../read/fixtures/signal.bw");
    const STORED_BW: &[u8] = include_bytes!("../read/fixtures/signal.unc.bw");
    const GENES_BB: &[u8] = include_bytes!("../read/fixtures/genes.bb");
    const SCORES_BB: &[u8] = include_bytes!("../read/fixtures/scores.bb");
    const REF_2BIT: &[u8] = include_bytes!("../read/fixtures/ref.2bit");
    const STEPS_BW: &[u8] = include_bytes!("../read/fixtures/steps.bw");

    /// Two reads on `chr1` of `ref.fa`, one over its run of N.
    const REF_READS: &str = "@SQ\tSN:chr1\tLN:120\n\
        r1\t0\tchr1\t3\t60\t10M\t*\t0\t0\tGTAAAATGAT\t*\n\
        r2\t16\tchr1\t31\t60\t12M\t*\t0\t0\tTCTTTCTCACGG\t*\n";

    /// The bigWig, bigBed and 2bit fixtures, and the text each was written
    /// from under the same name with the text's ending.
    fn binaries() -> Held {
        let mut held = Held::new();
        held.insert("signal.bw", SIGNAL_BW);
        held.insert("genes.bb", GENES_BB);
        held.insert("scores.bb", SCORES_BB);
        held.insert("ref.2bit", REF_2BIT);
        held.insert(
            "signal.bedgraph",
            include_str!("../read/fixtures/signal.bedgraph"),
        );
        held.insert("genes.bed", include_str!("../read/fixtures/genes.bed"));
        held.insert("scores.bed", include_str!("../read/fixtures/scores.bed"));
        held.insert("ref.fa", include_str!("../read/fixtures/ref.fa"));
        held.insert("reads.sam", REF_READS);
        held
    }

    fn held_figure(held: &mut Held, line: &str) -> Result<String, BuildError> {
        build_files(&invocation(line), held, |_, _| None)
    }

    /// A bigWig, a bigBed and a 2bit draw the figure the text they were
    /// written from draws, byte for byte, in every track that reads them: a
    /// signal, windows, scores under bases, genes and their exons, open
    /// reading frames, and a reference reads are read against.
    #[test]
    fn a_bigwig_a_bigbed_and_a_2bit_draw_what_their_text_draws() {
        let mut held = binaries();
        for (binary, text) in [
            (
                "chr1:1-1000 signal.bw genes.bb ref.2bit",
                "chr1:1-1000 signal.bedgraph genes.bed ref.fa",
            ),
            (
                "chr1:1-1000 genes.bb --isoforms --windows signal.bw",
                "chr1:1-1000 genes.bed --isoforms --windows signal.bedgraph",
            ),
            (
                "chr1:1-120 --orfs ref.2bit --dynseq signal.bw --with-sequence ref.2bit",
                "chr1:1-120 --orfs ref.fa --dynseq signal.bedgraph --with-sequence ref.fa",
            ),
            (
                "chr1:1-60 ref.2bit --pileup reads.sam",
                "chr1:1-60 ref.fa --pileup reads.sam",
            ),
            (
                "chr1:1-60 --pileup reads.sam --with-sequence ref.2bit",
                "chr1:1-60 --pileup reads.sam --with-sequence ref.fa",
            ),
            (
                "chr1:1-100 chr2:1-100 signal.bw genes.bb",
                "chr1:1-100 chr2:1-100 signal.bedgraph genes.bed",
            ),
            ("chr2 signal.bw", "chr2:1-500 signal.bedgraph"),
            (
                "geneA genes.bb ref.2bit",
                "chr1:101-1,000 genes.bed ref.fa --title geneA",
            ),
        ] {
            let drawn = held_figure(&mut held, binary).unwrap();
            assert_eq!(drawn, held_figure(&mut held, text).unwrap(), "{binary}");
        }
        // The pileup did read against the 2bit: a mismatch is drawn.
        let with = held_figure(&mut held, "chr1:1-60 ref.2bit --pileup reads.sam").unwrap();
        let without = held_figure(&mut held, "chr1:1-60 --pileup reads.sam").unwrap();
        assert_ne!(with, without);
    }

    const CONTACTS_HIC: &[u8] = include_bytes!("../read/fixtures/hic/contacts.hic");

    /// The cells `hictk dump --join` printed for one window of the `.hic`
    /// fixture, as BEDPE, by the line `make.sh` heads them with.
    fn dumped(window: &str) -> String {
        let dump = include_str!("../read/fixtures/hic/contacts.dump");
        let mut out = String::new();
        let mut inside = false;
        for line in dump.lines() {
            if let Some(heading) = line.strip_prefix("# ") {
                inside = heading == window;
            } else if inside {
                out.push_str(line);
                out.push('\n');
            }
        }
        assert!(!out.is_empty(), "{window}");
        out
    }

    /// The `.hic` fixture, and BEDPE files of what hictk dumped of it.
    fn contact_maps() -> Held {
        let mut held = Held::new();
        held.insert("contacts.hic", CONTACTS_HIC);
        held.insert("near.bedpe", dumped("1000 chr1 300000 420000"));
        held.insert("unaligned.bedpe", dumped("5000 chr1 612345 987654"));
        held.insert("dense.bedpe", dumped("250000 chr2 0 345678"));
        held
    }

    /// A `.hic` draws, byte for byte, the figure the BEDPE hictk dumps for
    /// the same window and resolution draws: a window inside one block at
    /// the finest resolution, one whose edges fall inside bins, and a dense
    /// block. And it is a triangle, where BEDPE as sparse as these cells is
    /// drawn as arcs unless told.
    #[test]
    fn a_hic_draws_what_the_bedpe_hictk_dumps_for_it_draws() {
        let mut held = contact_maps();
        for (hic, text) in [
            (
                "chr1:300,001-420,000 contacts.hic",
                "chr1:300,001-420,000 --pairs near.bedpe --style triangle --label contacts",
            ),
            (
                "chr1:612,346-987,654 --pairs contacts.hic --resolution 5000",
                "chr1:612,346-987,654 --pairs unaligned.bedpe --style triangle --label contacts",
            ),
            (
                "chr2:1-345,678 --pairs contacts.hic --resolution 250000 --log",
                "chr2:1-345,678 --pairs dense.bedpe --style triangle --label contacts --log",
            ),
        ] {
            let drawn = held_figure(&mut held, hic).unwrap();
            assert_eq!(drawn, held_figure(&mut held, text).unwrap(), "{hic}");
            assert!(drawn.contains("<polygon"), "{hic}");
        }
        let arcs = held_figure(&mut held, "chr1:300,001-420,000 --pairs near.bedpe").unwrap();
        assert!(!arcs.contains("<polygon"));
        assert!(held.notes.is_empty(), "{:?}", held.notes);
    }

    /// Without `--resolution`, a `.hic` is drawn at the finest of its
    /// resolutions that keeps the window to 250 bins, and a note says so
    /// where it holds finer; where none does, at its coarsest, which a note
    /// says too. Placed on a sequence's name, it is as long as the header
    /// says.
    #[test]
    fn a_hic_is_drawn_at_the_finest_resolution_that_keeps_the_window_to_250_bins() {
        let mut held = contact_maps();
        let whole = held_figure(&mut held, "chr1 contacts.hic").unwrap();
        assert_eq!(
            held.notes,
            [
                "contacts.hic is drawn at 5,000-base bins, the finest of its 6 resolutions that \
                 keeps the window to 250 bins; --resolution 2000 draws finer"
            ]
        );
        let region = Region::new("chr1", 0, 1_234_567).unwrap();
        let cells =
            read::hic::contacts(std::io::Cursor::new(CONTACTS_HIC), &region, 5_000).unwrap();
        assert_eq!(whole.matches("<polygon").count(), cells.len());
        assert!(whole.contains("1,230,001-1,234,567"), "the last bin, cut");
        let asked = held_figure(&mut held, "chr1 --pairs contacts.hic --resolution 5000").unwrap();
        assert_eq!(asked, whole);
        held.notes.clear();
        held_figure(&mut held, "chr1:1-100,000,000 contacts.hic").unwrap();
        assert_eq!(
            held.notes,
            [
                "contacts.hic is drawn at 250,000-base bins, its coarsest, which cut the window \
              into 400"
            ]
        );
        // A window the finest resolution keeps to 250 bins says nothing.
        held.notes.clear();
        held_figure(&mut held, "chr1:300,001-420,000 contacts.hic").unwrap();
        assert!(held.notes.is_empty(), "{:?}", held.notes);
    }

    #[test]
    fn a_resolution_a_hic_has_not_got_or_one_after_another_file_is_refused() {
        let mut held = contact_maps();
        let error = held_figure(&mut held, "chr1 --pairs contacts.hic --resolution 7000")
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "--pairs contacts.hic has no 7,000-base resolution; it holds 1,000, 2,000, 5,000, \
             10,000, 50,000 and 250,000"
        );
        let error = held_figure(
            &mut held,
            "chr1:1-1000 --pairs near.bedpe --resolution 1000",
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            error,
            "--pairs near.bedpe: --resolution picks one of the resolutions a .hic holds, and \
             near.bedpe is drawn at the bins it was written in"
        );
    }

    /// A `.hic` handed to a track that does not draw contacts names the one
    /// that does, one piped in asks for its name, and a window on a sequence
    /// it does not have names the ones it has, and among several places is
    /// a panel with nothing in it.
    #[test]
    fn a_hic_elsewhere_than_pairs_or_on_another_sequence_says_what_reads_it() {
        let mut held = contact_maps();
        let error = held_figure(&mut held, "chr1:1-1000 --coverage contacts.hic")
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "--coverage contacts.hic: the file is a contact map in Juicer's .hic, contacts \
             between the bins of a sequence, which --pairs draws"
        );
        let error = held_figure(&mut held, "chr3:1-1000 contacts.hic")
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "--pairs contacts.hic: the .hic has no sequence called chr3; it has chr1, chr2"
        );
        let both = held_figure(&mut held, "chr1:300,001-420,000 chr3:1-1000 contacts.hic").unwrap();
        assert!(both.contains("no pairs here"), "{both}");
        let piped = build(&over("chr1:1-1000", "--pairs", "-"), |_: &Source| {
            decoded(CONTACTS_HIC.to_vec(), None)
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            piped,
            "--pairs standard input: the file is a contact map in Juicer's .hic, which is read \
             through the index it holds, and a pipe cannot be read that way; name the file \
             instead"
        );
        held.insert("contacts.hic.gz", gzip_of(CONTACTS_HIC));
        let error = held_figure(&mut held, "chr1:1-1000 contacts.hic.gz")
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "--pairs contacts.hic.gz: the file is a contact map in Juicer's .hic compressed \
             with gzip, and a .hic is read through the index it holds, which the compression \
             hides; gunzip -k contacts.hic.gz writes contacts.hic, which karyon reads as it is"
        );
    }

    /// A `.hic` named for a file a track reads as text, the linkage `--ld`
    /// gives a scan, the links between loci or a sheet of traits, names the
    /// track that draws it, and is not asked for the name it was given. A
    /// file is a pipe by its own source and not by its track's: a bigWig
    /// named after `--ld`, with the scan piped in, is answered with the
    /// command that writes it as text, and a `.hic` piped in after `--ld` is
    /// asked for its name.
    #[test]
    fn a_hic_named_for_a_file_read_as_text_names_what_draws_it_and_not_a_pipe() {
        let scan = "CHR\tBP\tSNP\tP\nchr1\t5000\trs1\t0.001\nchr1\t50000\trs2\t0.2\n";
        let mut held = contact_maps();
        held.insert("scan.tsv", scan);
        held.insert("loci.bed", "chr1\t100\t2000\tg1\nchr1\t40000\t42000\tg2\n");
        held.insert("matrix.tsv", "gene\ts1\ts2\nA\t1\t2\nB\t3\t4\n");
        for (line, track) in [
            (
                "chr1:1-100,000 --manhattan scan.tsv --ld contacts.hic",
                "manhattan",
            ),
            (
                "chr1:1-100,000 --loci loci.bed --links contacts.hic",
                "loci",
            ),
            (
                "chr1:1-100,000 --matrix matrix.tsv --traits contacts.hic",
                "matrix",
            ),
        ] {
            let error = held_figure(&mut held, line).unwrap_err().to_string();
            assert_eq!(
                error,
                format!(
                    "--{track} contacts.hic: the file is a contact map in Juicer's .hic, \
                     contacts between the bins of a sequence, which --pairs draws"
                ),
                "{line}"
            );
        }
        let error = build(
            &invocation("chr1:1-100,000 --manhattan - --ld signal.bw"),
            |source: &Source| match source {
                Source::Stdin => Ok(scan.to_string()),
                Source::Path(path) => decoded(SIGNAL_BW.to_vec(), Some(path)),
            },
        )
        .unwrap_err()
        .to_string();
        let command = "bigWigToBedGraph -chrom=chr1 -start=0 -end=100000 signal.bw /dev/stdout";
        assert!(error.contains(&in_place(HOST, command, None)), "{error}");
        let error = build(
            &invocation("chr1:1-100,000 --manhattan scan.tsv --ld -"),
            |source: &Source| match source {
                Source::Stdin => decoded(CONTACTS_HIC.to_vec(), None),
                Source::Path(_) => Ok(scan.to_string()),
            },
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            error,
            "--manhattan standard input: the file is a contact map in Juicer's .hic, which is \
             read through the index it holds, and a pipe cannot be read that way; name the \
             file instead"
        );
    }

    /// The same files on disk draw what they draw held in memory.
    #[test]
    fn a_binary_on_disk_draws_what_it_draws_held_in_memory() {
        let dir = Scratch::new("binaries");
        let signal = dir.write("signal.bw", SIGNAL_BW);
        let genes = dir.write("genes.bb", GENES_BB);
        let reference = dir.write("ref.2bit", REF_2BIT);
        let mut held = Held::new();
        for (name, bytes) in [
            (&signal, SIGNAL_BW),
            (&genes, GENES_BB),
            (&reference, REF_2BIT),
        ] {
            held.insert(name.as_str(), bytes);
        }
        for line in [
            format!("chr1:1-1000 {signal} {genes} {reference}"),
            format!("chr2 {signal}"),
            format!("geneC {genes}"),
        ] {
            assert_eq!(
                held_figure(&mut held, &line).unwrap(),
                drawn_from_disk(&line).unwrap(),
                "{line}"
            );
        }
    }

    /// A sequence a bigWig, a bigBed or a 2bit names is as long as its index
    /// says, and a gene a bigBed names is placed with its margin.
    #[test]
    fn a_place_is_found_in_a_binary_file_s_index_and_rows() {
        let mut held = binaries();
        let along = |held: &mut Held, line: &str| {
            build_figure(&invocation(line), held, |_, _| None, Theme::light(), None)
                .unwrap()
                .along
                .map(|region| region.to_string())
        };
        assert_eq!(
            along(&mut held, "chr2 signal.bw").as_deref(),
            Some("chr2:1-500")
        );
        assert_eq!(
            along(&mut held, "chr3 ref.2bit").as_deref(),
            Some("chr3:1-12")
        );
        assert_eq!(
            along(&mut held, "chr2 genes.bb").as_deref(),
            Some("chr2:1-500")
        );
        // geneC runs from 10 to 400 on chr2, and a hundred bases each side
        // are cut at the sequence's start.
        assert_eq!(
            along(&mut held, "geneC genes.bb").as_deref(),
            Some("chr2:1-500")
        );
        assert_eq!(
            along(&mut held, "geneA genes.bb").as_deref(),
            Some("chr1:101-1000")
        );
        // A shade by a gene's name finds it in the bigBed too.
        assert!(held_figure(&mut held, "chr1:1-1000 genes.bb --shade geneB").is_ok());
        let error = held_figure(&mut held, "geneZ genes.bb")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("genes.bb names the sequences chr1 and chr2"),
            "{error}"
        );
        assert!(!error.contains("add the annotation"), "{error}");
    }

    /// A bigBed is read whole to find a gene by its name only up to a size,
    /// past which the figure says so rather than reading a gigabyte of rows.
    #[test]
    fn a_bigbed_too_large_to_read_whole_is_not_searched_for_a_name() {
        let mut held = binaries();
        let source = Source::Path("genes.bb".into());
        let rows = annotation_within(&mut held, &source, GENES_BB.len() as u64).unwrap();
        assert_eq!(rows, include_str!("../read/fixtures/genes.bb.bed"));
        assert!(held.notes.is_empty());
        assert_eq!(annotation_within(&mut held, &source, 1_000), None);
        assert_eq!(
            held.notes,
            [
                "genes.bb is a bigBed of 1 MB, which is not read whole to find a gene by its \
              name; write the gene's place, as chr1:1,001-2,000"
            ]
        );
        // And text is its own text, whatever the bound.
        let text = Source::Path("genes.bed".into());
        assert_eq!(
            annotation_within(&mut held, &text, 0).as_deref(),
            Some(include_str!("../read/fixtures/genes.bed"))
        );
    }

    /// A binary file handed to a track that does not draw what it holds is
    /// refused naming the tracks that do, and `--format` is refused, since
    /// the file says what it holds. On standard input it is refused asking
    /// for its name.
    #[test]
    fn a_binary_file_on_the_wrong_track_says_which_track_draws_it() {
        let mut held = binaries();
        for (line, said) in [
            (
                "chr1:1-100 --pileup signal.bw",
                "--pileup signal.bw: the file is bigWig, a signal along the sequence, which \
                 --coverage, --windows and --dynseq draw",
            ),
            (
                "chr1:1-100 --variants genes.bb",
                "--variants genes.bb: the file is bigBed, intervals along the sequence, which \
                 --features draws",
            ),
            (
                "chr1:1-100 --features ref.2bit",
                "--features ref.2bit: the file is 2bit, the bases of a reference, which \
                 --sequence and --orfs draw and --with-sequence reads",
            ),
            (
                "chr1:1-100 --dynseq signal.bedgraph --with-sequence signal.bw",
                "--dynseq signal.bw: the file is bigWig, a signal along the sequence, which \
                 --coverage, --windows and --dynseq draw",
            ),
            (
                "chr1:1-100 signal.bw --format bedgraph",
                "--coverage signal.bw: --format says what the columns of a text file are, \
                 and the file is bigWig, which says what it holds itself",
            ),
            (
                "chrX:1-100 signal.bw",
                "--coverage signal.bw: the bigWig has no sequence called chrX; it has chr1, \
                 chr2, chr3",
            ),
            (
                "chr2:101-200 ref.2bit",
                "--sequence ref.2bit: chr2 holds bases 1 to 37, and none of them is in \
                 chr2:101-200",
            ),
        ] {
            let error = held_figure(&mut held, line).unwrap_err().to_string();
            assert_eq!(error, said, "{line}");
        }
        let piped = build(&over("chr1:1-100", "--coverage", "-"), |_: &Source| {
            decoded(SIGNAL_BW.to_vec(), None)
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            piped,
            "--coverage standard input: the file is bigWig, which is read through the index \
             it holds, and a pipe cannot be read that way; name the file instead"
        );
    }

    /// A bigWig that calls the figure's sequence `1` is answered with the
    /// `--rename` that draws it, which then does.
    #[test]
    fn rename_reads_a_binary_file_under_its_own_names() {
        use crate::read::bigwig::fixture::{written, Items};
        let mut held = Held::new();
        let numbered = written(false, "1", 1000, &[Items::BedGraph(vec![(10, 20, 3.0)])]);
        held.insert("numbered.bw", numbered);
        let error = held_figure(&mut held, "chr1:1-100 numbered.bw").unwrap_err();
        assert_eq!(
            error.to_string(),
            "--coverage numbered.bw: the bigWig has no sequence called chr1; it has 1; if 1 \
             is chr1, add --rename 1=chr1"
        );
        let renamed = held_figure(&mut held, "chr1:1-100 numbered.bw --rename 1=chr1").unwrap();
        held.insert("numbered.bedgraph", "chr1\t10\t20\t3\n");
        assert_eq!(
            renamed,
            held_figure(&mut held, "chr1:1-100 numbered.bedgraph").unwrap()
        );
    }

    /// kent indexes only the sequences that hold data, so a bigWig or a
    /// bigBed written against every sequence's length can name one of them
    /// alone. Over several places, the place on a sequence it does not name
    /// holds none of its rows and says so there, as the text it was written
    /// from does, where the whole figure was refused; over that place alone,
    /// it is refused naming the sequences the file has.
    #[test]
    fn a_place_a_binary_file_does_not_name_holds_none_of_its_rows() {
        let mut held = binaries();
        let fixtures = [
            (
                "sparse.bw",
                &include_bytes!("../read/fixtures/sparse.bw")[..],
            ),
            ("sparse.bb", include_bytes!("../read/fixtures/sparse.bb")),
            (
                "sparse.bedgraph",
                include_bytes!("../read/fixtures/sparse.bedgraph"),
            ),
            ("sparse.bed", include_bytes!("../read/fixtures/sparse.bed")),
        ];
        for (name, bytes) in fixtures {
            held.insert(name, bytes);
        }
        let named = |bytes: &[u8]| -> Vec<String> {
            let names = match Binary::of(bytes, None) {
                Some(Binary::BigWig) => read::bigwig::sequences(io::Cursor::new(bytes)),
                _ => read::bigbed::sequences(io::Cursor::new(bytes)),
            };
            names.unwrap().into_iter().map(|(name, _)| name).collect()
        };
        assert_eq!(named(fixtures[0].1), ["chr1"]);
        assert_eq!(named(fixtures[1].1), ["chr1"]);
        for (binary, text) in [
            (
                "chr1:1-1000 chr2:1-500 sparse.bw sparse.bb --windows sparse.bw",
                "chr1:1-1000 chr2:1-500 sparse.bedgraph sparse.bed --windows sparse.bedgraph",
            ),
            (
                "geneA geneC genes.bed sparse.bw",
                "geneA geneC genes.bed sparse.bedgraph",
            ),
        ] {
            let drawn = held_figure(&mut held, binary).unwrap();
            assert_eq!(drawn, held_figure(&mut held, text).unwrap(), "{binary}");
        }
        let drawn = held_figure(
            &mut held,
            "chr1:1-1000 chr2:1-500 sparse.bw sparse.bb --windows sparse.bw",
        )
        .unwrap();
        for said in ["no values here", "no features here", "no windows here"] {
            assert_eq!(
                drawn.matches(&format!(">{said}</text>")).count(),
                1,
                "{said}"
            );
        }
        for (line, said) in [
            (
                "chr2:1-500 sparse.bw",
                "--coverage sparse.bw: the bigWig has no sequence called chr2; it has chr1; if \
                 chr1 is chr2, add --rename chr1=chr2",
            ),
            (
                "chr2:1-500 sparse.bb",
                "--features sparse.bb: the bigBed has no sequence called chr2; it has chr1; if \
                 chr1 is chr2, add --rename chr1=chr2",
            ),
        ] {
            let error = held_figure(&mut held, line).unwrap_err();
            assert!(matches!(error, BuildError::Absent { .. }), "{error:?}");
            assert_eq!(error.to_string(), said, "{line}");
        }
        // A file that names the sequence and holds nothing in the window is
        // refused as it was.
        assert!(matches!(
            held_figure(&mut held, "chr1:31-99 sparse.bw").unwrap_err(),
            BuildError::Empty { .. }
        ));
    }

    /// Drawn from a zoom level by its bins' means or their least, a bigWig
    /// is scaled to the most of its values under the window, as its values
    /// as written scale it. Scaled to what the bins were painted with, the
    /// top of this file's scale came down from 30 to 5 by their means, and
    /// by their least, every one of which has a gap and so is nought or
    /// less, there was no scale and no profile at all.
    #[test]
    fn a_bigwig_drawn_from_its_zoom_level_is_scaled_as_its_values_are() {
        let mut held = Held::new();
        held.insert("steps.bw", STEPS_BW);
        held.insert(
            "steps.bedgraph",
            include_str!("../read/fixtures/steps.bedgraph"),
        );
        // 6,000,000 bases over 900 pixels is 6,667 to a pixel, which holds
        // two of the finest level's bins of 3,184.
        let region = Region::parse("chr1:1-6,000,000").unwrap();
        let signal =
            read::bigwig::window(io::Cursor::new(STEPS_BW), &region, 6_667.0, Aggregate::Mean)
                .unwrap();
        assert_eq!(signal.zoom, Some(3_184));
        let scale = |svg: &str| -> Vec<String> {
            svg.split("<text")
                .filter_map(|text| text.split_once('>'))
                .filter_map(|(_, rest)| rest.split_once("</text>"))
                .map(|(label, _)| label.to_string())
                .filter(|label| label.parse::<f64>().is_ok())
                .collect()
        };
        for aggregate in ["max", "mean", "min"] {
            let line = |file: &str| {
                format!("chr1:1-6,000,000 --coverage {file} --aggregate {aggregate} --label steps")
            };
            let zoomed = held_figure(&mut held, &line("steps.bw")).unwrap();
            let written = held_figure(&mut held, &line("steps.bedgraph")).unwrap();
            assert_eq!(scale(&zoomed), scale(&written), "{aggregate}");
            assert!(scale(&zoomed).contains(&"30".to_string()), "{aggregate}");
            assert!(zoomed.contains("<path"), "{aggregate}: no profile");
        }
    }

    /// A value a bigWig holds as not a number, or as infinite, is read back
    /// as one, so it is drawn as the bedGraph that holds it draws it: as no
    /// value, a gap in the profile. Dropped, the bases under it were nought.
    #[test]
    fn a_bigwig_value_that_is_not_a_number_is_drawn_as_its_bedgraph_draws_it() {
        use crate::read::bigwig::fixture::{written, Items};
        let spans = vec![
            (0, 10, 1.0),
            (10, 20, f32::NAN),
            (20, 30, f32::INFINITY),
            (30, 40, 2.0),
        ];
        let bytes = written(false, "chr1", 100, &[Items::BedGraph(spans)]);
        let region = Region::parse("chr1:1-100").unwrap();
        let signal =
            read::bigwig::window(io::Cursor::new(&bytes), &region, 1.0, Aggregate::Max).unwrap();
        assert_eq!(
            read::bigwig::bedgraph("chr1", &signal.spans),
            "chr1\t0\t10\t1\nchr1\t10\t20\tNaN\nchr1\t20\t30\tinf\nchr1\t30\t40\t2\n"
        );
        let mut held = Held::new();
        held.insert("n.bw", bytes);
        held.insert(
            "n.bedgraph",
            "chr1\t0\t10\t1\nchr1\t10\t20\tnan\nchr1\t20\t30\tinf\nchr1\t30\t40\t2\n",
        );
        held.insert("dropped.bedgraph", "chr1\t0\t10\t1\nchr1\t30\t40\t2\n");
        let drawn = held_figure(&mut held, "chr1:1-100 n.bw").unwrap();
        assert_eq!(
            drawn,
            held_figure(&mut held, "chr1:1-100 n.bedgraph").unwrap()
        );
        assert_ne!(
            drawn,
            held_figure(&mut held, "chr1:1-100 dropped.bedgraph --label n").unwrap()
        );
    }

    /// A bigWig, a bigBed or a 2bit compressed with gzip is read as it is
    /// once it is out of the wrapper, and says so, with the `gunzip` that
    /// takes it out. It was told it was not text and given the UCSC command
    /// that writes its text, which refuses a compressed file too.
    #[test]
    fn a_gzipped_bigwig_bigbed_or_2bit_is_answered_with_the_gunzip_that_reads_it() {
        let mut held = binaries();
        for (name, bytes, line, said) in [
            (
                "signal.bw",
                SIGNAL_BW,
                "chr1:1-1000 signal.bw.gz",
                "--coverage signal.bw.gz: the file is bigWig compressed with gzip, and a bigWig \
                 is read through the index it holds, which the compression hides; gunzip -k \
                 signal.bw.gz writes signal.bw, which karyon reads as it is",
            ),
            (
                "genes.bb",
                GENES_BB,
                "chr1:1-1000 --features genes.bb.gz",
                "--features genes.bb.gz: the file is bigBed compressed with gzip, and a bigBed \
                 is read through the index it holds, which the compression hides; gunzip -k \
                 genes.bb.gz writes genes.bb, which karyon reads as it is",
            ),
            (
                "ref.2bit",
                REF_2BIT,
                "chr1:1-60 ref.2bit.gz",
                "--sequence ref.2bit.gz: the file is 2bit compressed with gzip, and a 2bit is \
                 read through the index it holds, which the compression hides; gunzip -k \
                 ref.2bit.gz writes ref.2bit, which karyon reads as it is",
            ),
        ] {
            held.insert(format!("{name}.gz"), gzip_of(bytes));
            let error = held_figure(&mut held, line).unwrap_err().to_string();
            assert_eq!(error, said, "{line}");
            // Done as it says, the file it names draws.
            assert!(held_figure(&mut held, &line.replace(".gz", "")).is_ok());
        }
        let piped = build(&over("chr1:1-100", "--coverage", "-"), |_: &Source| {
            decoded(gzip_of(SIGNAL_BW), None)
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            piped,
            "--coverage standard input: the file is bigWig compressed with gzip, and a bigWig \
             is read through the index it holds, which the compression hides; decompress it \
             into a file and name that file, which karyon reads as it is"
        );
    }

    /// Drawn many bases to a pixel, a bigWig's coverage is read from its zoom
    /// level, and drawn few, from its values as written. The first bin of
    /// this file's one level is made to say 99 where its values say 7, so
    /// which of the two a figure was drawn from shows in its scale.
    #[test]
    fn a_bigwig_drawn_many_bases_to_a_pixel_reads_its_zoom_level() {
        let mut bent = STORED_BW.to_vec();
        // The most of the bin from 129 to 248 made 99, and the least of the
        // one from 486 to 605 made -50, where kent wrote 7 for both. That one
        // has no gaps and sits in a column every base of which has a value,
        // so its least is the least the column is drawn with.
        let (most, least) = (920 + 20, 1016 + 16);
        assert_eq!(&bent[920..932], &[0, 0, 0, 0, 129, 0, 0, 0, 248, 0, 0, 0]);
        assert_eq!(
            &bent[1016..1032],
            &[0, 0, 0, 0, 230, 1, 0, 0, 93, 2, 0, 0, 119, 0, 0, 0]
        );
        for at in [most, least] {
            assert_eq!(bent[at..at + 4], 7.0f32.to_le_bytes());
        }
        bent[most..most + 4].copy_from_slice(&99.0f32.to_le_bytes());
        let mut held = binaries();
        held.insert("most.bw", bent.clone());
        bent[least..least + 4].copy_from_slice(&(-50.0f32).to_le_bytes());
        held.insert("signal.bw", bent);
        // The level's bins are 119 bases, so a pixel holds two of them from
        // 238 bases on: 300,000 bases over 900 pixels is 333.
        let zoomed = held_figure(&mut held, "chr1:1-300,000 signal.bw").unwrap();
        let written = held_figure(&mut held, "chr1:1-1,000 signal.bw").unwrap();
        // Values up to 7 put a tick at 7.5 on the scale, and up to 99 do not.
        let seven = |svg: &str| svg.contains(">7.5</text>");
        assert!(!seven(&zoomed), "the zoomed figure did not read the level");
        assert!(seven(&written), "the figure of few bases read the level");
        // Drawn by its least, the bin paints its least, and the figure is
        // not the one the file without that change draws. Painted with its
        // most instead, that column's least would be 7 either way.
        let least = |held: &mut Held, file: &str| {
            let line = format!("chr1:1-300,000 --coverage {file} --aggregate min --label signal");
            held_figure(held, &line).unwrap()
        };
        assert_ne!(
            least(&mut held, "signal.bw"),
            least(&mut held, "most.bw"),
            "the least of the bin was not painted"
        );
        // Windows are drawn whole and read as written at any scale.
        for locus in ["chr1:1-300,000", "chr1:1-1,000"] {
            assert_eq!(
                held_figure(&mut held, &format!("{locus} --windows signal.bw")).unwrap(),
                held_figure(
                    &mut held,
                    &format!("{locus} --windows signal.bedgraph --label signal")
                )
                .unwrap(),
                "{locus}"
            );
        }
    }

    /// Rows written from a binary file are what they were written as: a
    /// bigBed of numbered rows is not taken for a signal, and a bigWig over
    /// a few hundred bases is not told it could be drawn as reads.
    #[test]
    fn rows_written_from_a_binary_file_are_not_guessed_at() {
        let mut held = binaries();
        let bigbed = held_figure(&mut held, "chr1:1-100 scores.bb").unwrap();
        assert_eq!(
            bigbed,
            held_figure(&mut held, "chr1:1-100 --features scores.bed").unwrap()
        );
        assert_ne!(
            bigbed,
            held_figure(&mut held, "chr1:1-100 scores.bed").unwrap()
        );
        held_figure(&mut held, "chr1:1-300 signal.bw").unwrap();
        assert!(held.notes.is_empty(), "{:?}", held.notes);
    }

    /// Rows written from a binary file are read as what they were written
    /// as, not guessed at: a bigWig's spans that overlap are drawn as a
    /// bedGraph given `--format bedgraph` draws them, where a text file's are
    /// refused as two samples' depth, and a bigBed row whose seventh column is
    /// a dot is BED, where in a text file the dot says GFF3.
    #[test]
    fn rows_written_from_a_binary_file_are_read_as_what_they_were_written_as() {
        use crate::read::bigbed::fixture::written as bigbed;
        use crate::read::bigwig::fixture::{written as bigwig, Items};
        let mut held = Held::new();
        let overlapping = vec![(10, 20, 1.0), (15, 25, 2.0)];
        held.insert(
            "over.bw",
            bigwig(false, "chr1", 100, &[Items::BedGraph(overlapping)]),
        );
        held.insert("over.bedgraph", "chr1\t10\t20\t1\nchr1\t15\t25\t2\n");
        let drawn = held_figure(&mut held, "chr1:1-100 over.bw").unwrap();
        assert_eq!(
            drawn,
            held_figure(
                &mut held,
                "chr1:1-100 --coverage over.bedgraph --format bedgraph --label over"
            )
            .unwrap()
        );
        assert!(held_figure(&mut held, "chr1:1-100 --coverage over.bedgraph").is_err());
        held.insert(
            "dots.bb",
            bigbed("chr1", 100, &[(10, 40, "geneD\t0\t+\t.")], 7),
        );
        held.insert("dots.bed", "chr1\t10\t40\tgeneD\t0\t+\t.\n");
        let drawn = held_figure(&mut held, "chr1:1-100 dots.bb").unwrap();
        assert!(drawn.contains("geneD"), "{drawn}");
        assert_eq!(
            drawn,
            held_figure(
                &mut held,
                "chr1:1-100 --features dots.bed --format bed --label dots"
            )
            .unwrap()
        );
        assert!(held_figure(&mut held, "chr1:1-100 --features dots.bed").is_err());
    }

    /// Four samples called at three sites of chr1, named as `ROWS_TREE`
    /// names its tips.
    const COHORT_VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\tthree\tfour
chr1\t100\trs1\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t./.
chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t0/0\t1/1
chr1\t700\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/0\t0/1\t0/0
";

    /// Named on its own, a VCF is its calls, and a cohort's says the samples
    /// can be drawn a row each. One sample's would be a row repeating its
    /// sites, a track asked for by name needs no telling, and a first
    /// attempt read again under a --rename says it once.
    #[test]
    fn a_vcf_of_several_samples_named_on_its_own_says_its_genotypes_can_be_drawn() {
        let held = [("cohort.vcf", COHORT_VCF)];
        let (svg, notes) = drawn_noting("chr1:1-1000 cohort.vcf", &held);
        assert!(svg.unwrap().contains("<circle"), "drawn as calls");
        assert_eq!(
            notes,
            ["cohort.vcf is drawn as its calls; --genotypes cohort.vcf draws its 4 samples, a row each"]
        );
        let (_, notes) = drawn_noting("chr1:1-1000 --variants cohort.vcf", &held);
        assert!(notes.is_empty(), "{notes:?}");
        // Drawn as both, the figure already shows what the note would offer.
        let (svg, notes) = drawn_noting("chr1:1-1000 cohort.vcf --genotypes cohort.vcf", &held);
        assert!(svg.is_ok());
        assert!(notes.is_empty(), "{notes:?}");

        let single = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                      chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t1\n";
        let (_, notes) = drawn_noting("chr1:1-1000 one.vcf", &[("one.vcf", single)]);
        assert!(notes.is_empty(), "{notes:?}");

        let numbered = COHORT_VCF.replace("\nchr1\t", "\n1\t");
        let (svg, notes) = drawn_noting(
            "chr1:1-1000 cohort.vcf --rename 1=chr1",
            &[("cohort.vcf", numbered.as_str())],
        );
        assert!(svg.is_ok());
        assert_eq!(notes.len(), 1, "{notes:?}");
    }

    /// A name the VCF has not got is refused with the names it has, the
    /// first five of them: a cohort's header can name thousands.
    #[test]
    fn sample_names_a_sample_the_vcf_has_not_got() {
        let names: Vec<String> = (1..=40).map(|n| format!("S{n:02}")).collect();
        let vcf = format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n\
             chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t{}\n",
            names.join("\t"),
            vec!["0"; 40].join("\t")
        );
        let error = drawn_from(
            "chr1:1-1000 --genotypes c.vcf --sample S03,S99",
            &[("c.vcf", vcf.as_str())],
        )
        .unwrap_err();
        assert!(
            matches!(&error, BuildError::Unnamed { what: "sample", wanted, .. } if wanted == "S99"),
            "{error:?}"
        );
        assert_eq!(
            error.to_string(),
            "--genotypes c.vcf has no sample called S99; it has S01, S02, S03, S04, S05 and \
             35 more samples"
        );
    }

    /// `--matrix` read a VCF's first record as its header of sites and
    /// refused it for a word where a position goes.
    #[test]
    fn a_vcf_given_to_matrix_names_genotypes() {
        let error = drawn_from(
            "chr1:1-1000 --matrix cohort.vcf",
            &[("cohort.vcf", COHORT_VCF)],
        )
        .unwrap_err();
        let said = error.to_string();
        assert!(said.contains("this is a VCF"), "{said}");
        assert!(said.contains("--genotypes"), "{said}");
    }

    /// Rows on another sequence are rows elsewhere, and the refusal says
    /// where, with the --rename that would draw them.
    #[test]
    fn genotypes_on_another_sequence_say_where_the_rows_are() {
        let error = drawn_from(
            "chr2:1-1000 --genotypes cohort.vcf",
            &[("cohort.vcf", COHORT_VCF)],
        )
        .unwrap_err();
        assert!(
            matches!(
                &error,
                BuildError::Elsewhere {
                    wanted: "genotypes",
                    held: 3,
                    ..
                }
            ),
            "{error:?}"
        );
        assert!(
            error.to_string().contains("add --rename chr1=chr2"),
            "{error}"
        );
    }

    /// One place of several with no record on it is a band saying so,
    /// rather than the whole figure refused.
    #[test]
    fn a_place_of_several_with_no_genotypes_says_so_there() {
        let svg = drawn_from(
            "chr1:1-1000 chr1:5000-6000 --genotypes cohort.vcf",
            &[("cohort.vcf", COHORT_VCF)],
        )
        .unwrap();
        assert_eq!(svg.matches(">no genotypes here</text>").count(), 1);
    }

    /// A bgzipped VCF is the text inside, read the same way.
    #[test]
    fn genotypes_from_a_bgzipped_vcf_match_the_plain_one() {
        let dir = Scratch::new("genotypes-gz");
        let plain = dir.write("cohort.vcf", COHORT_VCF.as_bytes());
        let gz = dir.write("cohort.vcf.gz", &gzip_of(COHORT_VCF.as_bytes()));
        let svg = drawn_from_disk(&format!("chr1:1-1000 --genotypes {plain} --label c")).unwrap();
        let unwrapped =
            drawn_from_disk(&format!("chr1:1-1000 --genotypes {gz} --label c")).unwrap();
        assert!(svg.contains("one, 3 of 3 sites called"), "{svg}");
        assert_eq!(svg, unwrapped);
    }

    /// Drawn under the calls of the same file, two bands both called `calls`
    /// would be left for the reader to tell apart.
    #[test]
    fn the_label_of_a_genotypes_track_says_genotypes() {
        let held = [("calls.vcf.gz", COHORT_VCF)];
        let svg = drawn_from("chr1:1-1000 calls.vcf.gz --genotypes calls.vcf.gz", &held).unwrap();
        assert!(svg.contains(">calls genotypes</text>"), "{svg}");
        assert!(svg.contains(">calls</text>"), "{svg}");
    }

    /// The tree puts the rows in the order of its tips, and --sample picks
    /// which rows there are, in its own order when there is no tree.
    #[test]
    fn with_tree_orders_genotype_rows() {
        let held = [("c.vcf", COHORT_VCF), ("t.nwk", ROWS_TREE)];
        let plain = drawn_from("chr1:1-1000 --genotypes c.vcf", &held).unwrap();
        assert_eq!(rows_drawn(&plain), ["one", "two", "three", "four"]);
        let ordered = drawn_from("chr1:1-1000 --genotypes c.vcf --with-tree t.nwk", &held).unwrap();
        assert_eq!(rows_drawn(&ordered), ["four", "two", "three", "one"]);
        let picked = drawn_from("chr1:1-1000 --genotypes c.vcf --sample three,one", &held).unwrap();
        assert_eq!(rows_drawn(&picked), ["three", "one"]);
        let both = drawn_from(
            "chr1:1-1000 --genotypes c.vcf --sample three,one,two --with-tree t.nwk",
            &held,
        )
        .unwrap();
        assert_eq!(rows_drawn(&both), ["two", "three", "one"]);
        assert!(both.contains("1 tip of the tree has no row"), "{both}");
    }

    /// The options of a track of rows reach the genotype rows: two rows of
    /// four, the other two counted, no names, and rows six pixels tall.
    #[test]
    fn row_options_after_genotypes_reach_the_rows() {
        let held = [("c.vcf", COHORT_VCF)];
        let svg = drawn_from(
            "chr1:1-1000 --genotypes c.vcf --max-rows 2 --row-height 6",
            &held,
        )
        .unwrap();
        assert_eq!(rows_drawn(&svg), ["one", "two"]);
        assert!(svg.contains(">+2 more</text>"), "{svg}");
        assert!(svg.contains(" height=\"6\""), "rows six pixels tall: {svg}");
        assert!(
            !svg.contains(" height=\"11\""),
            "no row at the default: {svg}"
        );
        let named = drawn_from("chr1:1-1000 --genotypes c.vcf", &held).unwrap();
        assert!(named.contains(">three</text>"), "{named}");
        let unnamed = drawn_from("chr1:1-1000 --genotypes c.vcf --no-names", &held).unwrap();
        assert!(!unnamed.contains(">three</text>"), "{unnamed}");
    }

    /// A sheet beside the rows is joined through the same door as a matrix's,
    /// so a column it has not got is refused with the ones it has, and its
    /// levels are keyed under the figure.
    #[test]
    fn traits_beside_genotype_rows_are_a_sheet_joined_by_name() {
        let sheet = "sample\tlineage\none\tL1\ntwo\tL2\nthree\tL1\nfour\tL2\n";
        let held = [("c.vcf", COHORT_VCF), ("s.tsv", sheet)];
        let svg = drawn_from(
            "chr1:1-1000 --genotypes c.vcf --traits s.tsv --columns lineage",
            &held,
        )
        .unwrap();
        assert!(svg.contains("one; lineage L1"), "{svg}");
        assert!(
            svg.contains(">lineage: L2</text>"),
            "the key names the levels: {svg}"
        );
        let error = drawn_from(
            "chr1:1-1000 --genotypes c.vcf --traits s.tsv --columns linage",
            &held,
        )
        .unwrap_err();
        assert!(
            matches!(error, BuildError::Unnamed { what: "column", .. }),
            "{error:?}"
        );
    }

    /// Files on disk that say which of them were read whole, as the command
    /// line reads them, for telling a window read through its index from the
    /// same figure drawn from the whole file.
    struct Watched {
        disk: Disk,
        whole: Vec<String>,
    }

    impl Files for Watched {
        fn text(&mut self, source: &Source) -> io::Result<String> {
            self.whole.push(called(source));
            self.disk.text(source)
        }

        fn seekable(&mut self, source: &Source) -> io::Result<Option<Box<dyn Seekable>>> {
            self.disk.seekable(source)
        }

        fn beside(&mut self, source: &Source, ending: &str) -> io::Result<Option<Beside>> {
            self.disk.beside(source, ending)
        }

        fn note(&mut self, message: &str) {
            self.disk.note(message);
        }
    }

    /// What a command line draws from disk, or the refusal as it is printed,
    /// with the files it read whole and what it noted.
    fn watched(line: &str) -> (String, Vec<String>, Vec<String>) {
        let mut files = Watched {
            disk: Disk::default(),
            whole: Vec::new(),
        };
        let drawn = build_files(&invocation(line), &mut files, |_, _| None)
            .unwrap_or_else(|error| format!("refused: {error}"));
        (drawn, files.whole, files.disk.notes)
    }

    /// Copies the bgzipped files `src/read/fixtures/indexed/make.py` wrote,
    /// each with the indexes `with` names after it, into `dir`, the file
    /// first, so its index is no older than it.
    fn indexed_into(dir: &Scratch, names: &[&str], with: &[&str]) -> Vec<String> {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/read/fixtures/indexed");
        names
            .iter()
            .map(|name| {
                let path = dir.write(name, &fs::read(fixtures.join(name)).unwrap());
                for ending in with {
                    let index = format!("{name}{ending}");
                    if fixtures.join(&index).is_file() {
                        dir.write(&index, &fs::read(fixtures.join(&index)).unwrap());
                    }
                }
                path
            })
            .collect()
    }

    /// Every kind read through an index draws what it draws from the whole
    /// file, over windows wide and narrow, empty, on a sequence the file does
    /// not name and on one it calls otherwise, and refuses what it refuses in
    /// the same words; and it is read through the index, which a figure
    /// cannot show, so the files read whole are counted too. A GTF, a GFF3
    /// whose exons name transcripts it has no row for, a bedMethyl with no
    /// code named, structural calls, a long table of windows and a BED its
    /// first row calls GFF3 are read whole with their index beside them, and
    /// so is a scan with no header over a window whose values all look like
    /// p-values, and each draws the same.
    #[test]
    fn an_indexed_file_draws_what_the_whole_file_draws() {
        let dir = Scratch::new("tabix-parity");
        let names = [
            "cohort.vcf.gz",
            "depth.bedgraph.gz",
            "genes.gff3.gz",
            "genes.bed.gz",
            "exons.gtf.gz",
            "exons.gff3.gz",
            "methyl.bed.gz",
            "sj.tab.gz",
            "scan1.tsv.gz",
            "scan2.tsv.gz",
            "windows.tsv.gz",
            "sv.vcf.gz",
            "scan0.tsv.gz",
            "long.tsv.gz",
            "odd.bed.gz",
        ];
        let paths = indexed_into(&dir, &names, &[".tbi"]);
        let at = |name: &str| paths[names.iter().position(|n| *n == name).unwrap()].clone();
        // chr2 longer than the bedGraph's rows reach, and chr3, on which it
        // has none.
        let fasta = format!(
            ">chr2\n{}\n>chr3\n{}\n",
            "ACGT".repeat(175_000),
            "ACGT".repeat(2_500)
        );
        let reference = dir.write("ref.fa", fasta.as_bytes());
        // Each command, whether it is read through the index, and the places
        // it is drawn over.
        let chromosomes = [
            "chr1:1-2,000",
            "chr1:15,001-17,000",
            "chr1:400,001-460,000",
            "chr1:1,000,001-1,003,000",
            "chr1:2,400,001-2,500,000",
            "chr2:1-600,000",
            "chr2:100,001-100,500",
            "chr3:1-1,000",
            "chr9:1-100",
            "1:400,001-460,000",
        ];
        let numbers = [
            "1:1-2,000",
            "1:400,001-460,000",
            "2:1-800,000",
            "chr1:1-10,000",
        ];
        let cases: Vec<(String, bool, &[&str])> = vec![
            (at("cohort.vcf.gz"), true, &chromosomes),
            (
                format!("--genotypes {}", at("cohort.vcf.gz")),
                true,
                &chromosomes,
            ),
            (at("depth.bedgraph.gz"), true, &chromosomes),
            (
                format!("--windows {}", at("depth.bedgraph.gz")),
                true,
                &chromosomes,
            ),
            (
                format!(
                    "--dynseq {} --with-sequence {reference}",
                    at("depth.bedgraph.gz")
                ),
                true,
                &[
                    "chr2:100,001-100,600",
                    "chr2:1-3,000",
                    "chr2:650,001-651,000",
                    "chr3:1-1,000",
                ],
            ),
            (at("genes.gff3.gz"), true, &chromosomes),
            (
                format!("{} --isoforms", at("genes.gff3.gz")),
                true,
                &chromosomes,
            ),
            (at("genes.bed.gz"), true, &chromosomes),
            (
                format!("--methylation {} --modification m", at("methyl.bed.gz")),
                true,
                &chromosomes,
            ),
            (
                format!("--junctions {}", at("sj.tab.gz")),
                true,
                &chromosomes,
            ),
            (
                format!("--manhattan {}", at("scan1.tsv.gz")),
                true,
                &numbers,
            ),
            (
                format!("--manhattan {}", at("scan2.tsv.gz")),
                true,
                &numbers,
            ),
            (
                format!("--heatmap {}", at("windows.tsv.gz")),
                true,
                &chromosomes,
            ),
            (at("exons.gtf.gz"), false, &chromosomes),
            // The same exons as GFF3, naming transcripts it has no row for:
            // through the index, a window inside an intron drew nothing and
            // one over an exon drew it alone.
            (at("exons.gff3.gz"), false, &chromosomes),
            (
                format!("--methylation {}", at("methyl.bed.gz")),
                false,
                &chromosomes,
            ),
            (
                format!("--structural {}", at("sv.vcf.gz")),
                false,
                &["chr1:400,001-400,100", "chr1:1-1,000,000"],
            ),
            // Every value over these places lies between nought and one, and
            // past 800 kb on 1 they do not: the file is no table of p-values,
            // and a window that looks like one is refused and read whole.
            (
                format!("--manhattan {}", at("scan0.tsv.gz")),
                false,
                &["1:1-200,000", "1:400,001-700,000", "2:1-300,000"],
            ),
            (
                format!("--manhattan {}", at("scan0.tsv.gz")),
                true,
                &["1:700,001-1,000,000"],
            ),
            // S3 has no row before 200 kb, and is a row of the heatmap there.
            (
                format!("--heatmap {}", at("long.tsv.gz")),
                false,
                &["chr1:1-100,000", "chr1:150,001-300,000"],
            ),
            // Its first row says GFF3 in column seven, and its others BED.
            (at("odd.bed.gz"), false, &["chr1:100,001-200,000"]),
        ];
        for (track, through, places) in &cases {
            let mut drawn = 0;
            let file = track
                .split_whitespace()
                .find(|word| word.ends_with(".gz"))
                .unwrap();
            for place in places.iter() {
                let line = format!("{place} {track}");
                let (indexed, whole, notes) = watched(&line);
                assert!(
                    !notes.iter().any(|note| note.contains(".tbi")),
                    "{line}: {notes:?}"
                );
                assert_eq!(
                    !whole.iter().any(|read| read == file),
                    *through,
                    "{line}: read whole {whole:?}"
                );
                let index = format!("{file}.tbi");
                fs::rename(&index, format!("{index}.aside")).unwrap();
                let (from_whole, _, _) = watched(&line);
                fs::rename(format!("{index}.aside"), &index).unwrap();
                assert!(indexed == from_whole, "{line}:\n{indexed}\n{from_whole}");
                drawn += usize::from(!indexed.starts_with("refused"));
            }
            // Each draws over some of its places, and is refused alike over
            // the rest, but for a bedMethyl of two codes and none named, and
            // a BED taken whole for GFF3, refused over all of them.
            let refused = (track.starts_with("--methylation") && !track.contains("--modification"))
                || track.ends_with("odd.bed.gz");
            assert_eq!(drawn == 0, refused, "{track} drew over {drawn} places");
        }
    }

    /// A track that says how many rows a file held elsewhere is read whole
    /// from an index that does not count them, which the format allows: the
    /// rows over a window are not the file's, and an empty window of
    /// junctions said the file held none.
    #[test]
    fn an_index_that_counts_no_rows_leaves_the_count_to_the_whole_file() {
        let text = "chr1\t101\t200\t1\t1\t1\t5\t0\t10\nchr1\t301\t400\t1\t1\t1\t3\t0\t10\n";
        let data = crate::read::bgzf::fixture::blocks(text.as_bytes(), &[]);
        let end = ((data.len() as u64) - 31) << 16;
        let index = crate::read::index::fixture::rooted(&["chr1"], [0, 1, 2, 3], (0, end));
        let line = "chr9:1-100 --junctions sj.tab.gz";
        let drawn = |with: bool| {
            let mut held = Held::new();
            held.insert("sj.tab.gz", data.clone());
            if with {
                held.insert("sj.tab.gz.tbi", index.clone());
            }
            let drawn = build_files(&invocation(line), &mut held, |_, _| None);
            drawn.map_err(|error| error.to_string())
        };
        let whole = drawn(false).unwrap_err();
        assert_eq!(
            whole,
            "--junctions sj.tab.gz: no junctions in chr9:1-100, though the file holds 2"
        );
        assert_eq!(drawn(true).unwrap_err(), whole);
        // And over the window the index finds the rows of, it draws.
        let mut held = Held::new();
        held.insert("sj.tab.gz", data.clone());
        held.insert("sj.tab.gz.tbi", index.clone());
        let line = "chr1:1-1,000 --junctions sj.tab.gz";
        assert!(build_files(&invocation(line), &mut held, |_, _| None).is_ok());
    }

    /// A `.bed.gz` named on its own is told by the file's first row, as the
    /// whole file tells it, though the window holds none of its rows: a
    /// bedGraph so named is a signal, refused over an empty window as a
    /// signal is, and read through its index as one.
    #[test]
    fn a_bed_named_on_its_own_is_told_by_the_file_s_first_row() {
        let dir = Scratch::new("tabix-told");
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/read/fixtures/indexed");
        let bed = dir.write(
            "signal.bed.gz",
            &fs::read(fixtures.join("depth.bedgraph.gz")).unwrap(),
        );
        dir.write(
            "signal.bed.gz.tbi",
            &fs::read(fixtures.join("depth.bedgraph.gz.tbi")).unwrap(),
        );
        for place in ["chr3:1-1,000", "chr1:400,001-460,000"] {
            let line = format!("{place} {bed}");
            let (through, whole, _) = watched(&line);
            assert!(whole.is_empty(), "{line}: {whole:?}");
            fs::rename(format!("{bed}.tbi"), format!("{bed}.aside")).unwrap();
            assert_eq!(watched(&line).0, through, "{line}");
            fs::rename(format!("{bed}.aside"), format!("{bed}.tbi")).unwrap();
        }
        let (refused, _, _) = watched(&format!("chr3:1-1,000 {bed}"));
        assert!(refused.contains("--coverage"), "{refused}");
    }

    /// A `.csi`, as `bcftools index` and `tabix -C` write one, is read as a
    /// `.tbi` is, and in place of one beside it, as htslib reads them.
    #[test]
    fn a_csi_is_read_as_a_tbi_is_and_before_one() {
        let dir = Scratch::new("tabix-csi");
        let paths = indexed_into(&dir, &["cohort.vcf.gz", "genes.gff3.gz"], &[".csi"]);
        let line = format!("chr1:400,001-460,000 --variants {} {}", paths[0], paths[1]);
        let (through_csi, whole, notes) = watched(&line);
        assert!(whole.is_empty() && notes.is_empty(), "{whole:?} {notes:?}");
        // A damaged `.tbi` beside them is not read, as tabix does not read
        // it while there is a `.csi`.
        for path in &paths {
            dir.write(&format!("{}.tbi", shortened(path)), b"not an index");
        }
        assert_eq!(
            watched(&line),
            (through_csi.clone(), Vec::new(), Vec::new())
        );
        for path in &paths {
            fs::remove_file(format!("{path}.csi")).unwrap();
            fs::remove_file(format!("{path}.tbi")).unwrap();
        }
        assert_eq!(watched(&line).0, through_csi);
    }

    /// A window of a file with an index reads the rows over the window and
    /// no others: a VCF whose row on chr2 has seven columns, which tabix
    /// indexes and the reader of calls refuses, draws chr1 through the index
    /// and is refused read whole. It is drawn through `build_files`, so a
    /// `Files` between the command line and the disk that kept the file's
    /// bytes or its index to itself would turn it off, and fail here.
    #[test]
    fn an_indexed_file_is_read_only_over_the_window() {
        let dir = Scratch::new("tabix-seven");
        let seven = &indexed_into(&dir, &["seven.vcf.gz"], &[".tbi"])[0];
        let line = format!("chr1:1-1,000 {seven}");
        assert!(drawn_from_disk(&line).is_ok());
        fs::remove_file(format!("{seven}.tbi")).unwrap();
        let error = drawn_from_disk(&line).unwrap_err().to_string();
        assert!(
            error.contains("line 5: a VCF line has at least 8 columns"),
            "{error}"
        );
    }

    /// A circle written from base 1 is checked against the lengths the files
    /// state in their headers and indexes, and reads no file whole for how
    /// far its rows reach: a VCF with an index and no `##contig` was read
    /// whole, every row on every other sequence, for a reach then thrown away.
    /// An annotation is checked by its header, and refused there.
    #[test]
    fn a_circle_written_from_base_1_reads_an_indexed_file_over_itself_alone() {
        let dir = Scratch::new("tabix-circle");
        let paths = indexed_into(&dir, &["seven.vcf.gz", "genes.gff3.gz"], &[".tbi"]);
        let (drawn, whole, _) = watched(&format!("chr1:1-1,000 --circular {}", paths[0]));
        assert!(drawn.contains("chr1, 1,000 bases"), "{drawn}");
        assert!(whole.is_empty(), "read whole: {whole:?}");
        // As a band over the same span reads it.
        assert!(watched(&format!("chr1:1-1,000 {}", paths[0])).1.is_empty());
        // The GFF3 says chr1 is 2,000,000 bases in its header.
        let (refused, whole, _) = watched(&format!("chr1:1-1,000 --circular {}", paths[1]));
        assert!(
            refused.contains("closes the circle at 1,000, and the files say chr1 is 2,000,000"),
            "{refused}"
        );
        assert!(whole.is_empty(), "read whole: {whole:?}");
    }

    /// A refusal of a row over the window names its line in the file, and not
    /// its line in the window's text: the track is read again whole for it.
    #[test]
    fn a_bad_row_over_the_window_is_refused_on_its_line_in_the_file() {
        let dir = Scratch::new("tabix-line");
        let seven = &indexed_into(&dir, &["seven.vcf.gz"], &[".tbi"])[0];
        let line = format!("chr2:1-1,000 {seven}");
        let (through, whole, _) = watched(&line);
        assert_eq!(
            whole,
            std::slice::from_ref(seven),
            "read whole once the window was refused"
        );
        fs::remove_file(format!("{seven}.tbi")).unwrap();
        assert_eq!(through, watched(&line).0);
        assert!(through.contains("line 5:"), "{through}");
    }

    /// The commonest empty window, a file that calls the sequence otherwise,
    /// is answered with the counts the index keeps, word for word what the
    /// whole file says, and without reading it whole.
    #[test]
    fn the_empty_window_says_what_the_index_counts() {
        let dir = Scratch::new("tabix-elsewhere");
        let paths = indexed_into(&dir, &["cohort.vcf.gz", "sj.tab.gz"], &[".tbi"]);
        for line in [
            format!("1:1-1,000 {}", paths[0]),
            format!("chr9:1-1,000 --junctions {}", paths[1]),
            format!("chr1:1-100 --junctions {}", paths[1]),
        ] {
            let (through, whole, _) = watched(&line);
            assert!(whole.is_empty(), "{line}: {whole:?}");
            assert!(through.starts_with("refused"), "{through}");
            for path in &paths {
                fs::rename(format!("{path}.tbi"), format!("{path}.aside")).unwrap();
            }
            let (read_whole, _, _) = watched(&line);
            for path in &paths {
                fs::rename(format!("{path}.aside"), format!("{path}.tbi")).unwrap();
            }
            assert_eq!(through, read_whole, "{line}");
        }
        let (said, _, _) = watched(&format!("1:1-1,000 {}", paths[0]));
        assert!(
            said.ends_with("holds 450 on chr1, chr2 and chr3; if chr1 is 1, add --rename chr1=1"),
            "{said}"
        );
    }

    /// A figure placed by a gene's name reads the calls beside the annotation
    /// for their header alone, which says how long each sequence is; their
    /// rows are read for how far they reach only where nothing else places
    /// the figure, and once. A name nowhere is answered with the sequences
    /// the index names.
    #[test]
    fn a_gene_named_figure_reads_an_indexed_vcf_by_its_header() {
        let dir = Scratch::new("tabix-place");
        let paths = indexed_into(&dir, &["genes.gff3.gz", "cohort.vcf.gz"], &[".tbi"]);
        let named = |name: &str| format!("{name} {} {}", paths[0], paths[1]);
        let (gene, whole, _) = watched(&named("chr1g3"));
        assert!(gene.starts_with("<svg"), "{gene}");
        assert_eq!(
            whole,
            [paths[0].clone()],
            "the annotation alone is read whole"
        );
        let (sequence, whole, _) = watched(&format!("chr3 --variants {}", paths[1]));
        assert_eq!(
            locus_of(&sequence),
            "chr3:1-200000",
            "the length the header says"
        );
        assert!(whole.is_empty(), "{whole:?}");
        // A name that is neither reads the calls whole once, for how far
        // they reach, and the refusal names their sequences from the index.
        let (nowhere, whole, _) = watched(&named("nosuch"));
        assert_eq!(whole, paths, "{nowhere}");
        assert!(
            nowhere.ends_with(&format!(
                "{} names the sequences chr1, chr2 and chr3.",
                paths[1]
            )),
            "{nowhere}"
        );
        // A bedGraph says no length, so its rows are read for how far they
        // reach, once nothing else has placed the figure.
        let depth = &indexed_into(&dir, &["depth.bedgraph.gz"], &[".tbi"])[0];
        let (reached, whole, _) = watched(&format!("chr2 {depth}"));
        assert_eq!(whole, std::slice::from_ref(depth));
        assert!(locus_of(&reached).starts_with("chr2:1-"), "{reached}");
        for path in paths.iter().chain([depth]) {
            fs::rename(format!("{path}.tbi"), format!("{path}.aside")).unwrap();
        }
        assert_eq!(watched(&named("chr1g3")).0, gene);
        assert_eq!(
            watched(&format!("chr3 --variants {}", paths[1])).0,
            sequence
        );
        assert_eq!(watched(&named("nosuch")).0, nowhere);
        assert_eq!(watched(&format!("chr2 {depth}")).0, reached);
    }

    /// An index older than its file may be for an earlier version of it, so
    /// it is not trusted: the file is read whole, which draws the same, and
    /// a note says why and how to write the index again, once however many
    /// panels read the file.
    #[cfg(unix)]
    #[test]
    fn an_index_older_than_its_file_is_not_trusted() {
        let dir = Scratch::new("tabix-older");
        let paths = indexed_into(&dir, &["seven.vcf.gz", "cohort.vcf.gz"], &[".tbi"]);
        made_older(
            &paths
                .iter()
                .map(|path| format!("{path}.tbi"))
                .collect::<Vec<_>>(),
        );
        // Read whole, the row of seven columns is refused.
        assert!(drawn_from_disk(&format!("chr1:1-1,000 {}", paths[0])).is_err());
        let line = format!("chr1:1-2,000 chr1:400,001-460,000 --variants {}", paths[1]);
        let (drawn, whole, notes) = watched(&line);
        assert_eq!(whole.len(), 2, "read whole for each panel");
        assert_eq!(
            notes,
            [format!(
                "{0}.tbi is older than {0}, so it was not trusted and the file was read \
                 whole; tabix -f -p vcf {0} writes it again",
                paths[1]
            )]
        );
        fs::remove_file(format!("{}.tbi", paths[1])).unwrap();
        assert_eq!(watched(&line).0, drawn);
    }

    /// Dates each index in `paths` to the year 2000, older than the file
    /// beside it.
    #[cfg(unix)]
    fn made_older(paths: &[String]) {
        for path in paths {
            let touched = std::process::Command::new("touch")
                .args(["-t", "200001010000", path])
                .status()
                .unwrap();
            assert!(touched.success());
        }
    }

    /// The note for an index older than its file gives the command that
    /// writes it as it was written: `-C` for a `.csi`, which is looked for
    /// before a `.tbi`, so one written again without it left the older `.csi`
    /// found first and the same note said again; tabix's word for a GFF3; and
    /// the columns of a table tabix has no word for.
    #[cfg(unix)]
    #[test]
    fn an_older_index_is_asked_for_as_it_was_written() {
        let dir = Scratch::new("tabix-again");
        let csi = indexed_into(&dir, &["cohort.vcf.gz"], &[".csi"]);
        let tbi = indexed_into(
            &dir,
            &["genes.gff3.gz", "windows.tsv.gz", "sj.tab.gz"],
            &[".tbi"],
        );
        made_older(&[
            format!("{}.csi", csi[0]),
            format!("{}.tbi", tbi[0]),
            format!("{}.tbi", tbi[1]),
            format!("{}.tbi", tbi[2]),
        ]);
        for (line, index, options) in [
            (
                format!("chr1:400,001-460,000 --variants {}", csi[0]),
                "csi",
                "-C -p vcf",
            ),
            (format!("chr1:400,001-460,000 {}", tbi[0]), "tbi", "-p gff"),
            (
                format!("chr1:400,001-460,000 --heatmap {}", tbi[1]),
                "tbi",
                "-s1 -b2 -e3 -0 -S1",
            ),
            (
                format!("chr1:1-20,000 --junctions {}", tbi[2]),
                "tbi",
                "-s1 -b2 -e3",
            ),
        ] {
            let file = line.split_whitespace().last().unwrap().to_string();
            let (_, whole, notes) = watched(&line);
            assert_eq!(whole, std::slice::from_ref(&file), "{line}");
            assert_eq!(
                notes,
                [format!(
                    "{file}.{index} is older than {file}, so it was not trusted and the file \
                     was read whole; tabix -f {options} {file} writes it again"
                )],
                "{line}"
            );
        }
    }

    /// An index older than a file read whole whatever its index, a GTF, a
    /// GFF3 whose exons name transcripts it has no row for, or a bedMethyl
    /// drawn with no code named, is not asked for again: written again, it
    /// left the file read whole as before, and the note came back on the
    /// next figure of a GFF3 that has the same rows. The same index beside a
    /// bedMethyl drawn with a code named is asked for.
    #[cfg(unix)]
    #[test]
    fn an_older_index_beside_a_file_read_whole_anyway_is_not_asked_for() {
        let dir = Scratch::new("tabix-older-whole");
        let paths = indexed_into(
            &dir,
            &["exons.gtf.gz", "exons.gff3.gz", "methyl.bed.gz"],
            &[".tbi"],
        );
        made_older(
            &paths
                .iter()
                .map(|path| format!("{path}.tbi"))
                .collect::<Vec<_>>(),
        );
        let place = "chr1:400,001-460,000";
        for track in [
            paths[0].clone(),
            paths[1].clone(),
            format!("--methylation {}", paths[2]),
        ] {
            let (_, whole, notes) = watched(&format!("{place} {track}"));
            assert_eq!(whole.len(), 1, "{track}: {whole:?}");
            assert!(notes.is_empty(), "{track}: {notes:?}");
        }
        let line = format!("{place} --methylation {} --modification m", paths[2]);
        let (_, _, notes) = watched(&line);
        assert_eq!(
            notes,
            [format!(
                "{0}.tbi is older than {0}, so it was not trusted and the file was read \
                 whole; tabix -f -p bed {0} writes it again",
                paths[2]
            )]
        );
    }

    /// A GFF3 whose exons name a transcript it has no row for is read whole
    /// where the rows over the window show one, though its first rows have
    /// every transcript's: through the index the transcript was drawn as the
    /// one exon over the window, where the whole file draws it from its
    /// first exon to its last with the line of its intron. A window whose
    /// rows have their transcripts is read through the index still.
    #[test]
    fn exons_naming_no_row_over_the_window_send_the_file_whole() {
        let text = "##gff-version 3\n\
            chr1\t.\tgene\t1001\t9000\t.\t+\t.\tID=g1\n\
            chr1\t.\tmRNA\t1001\t9000\t.\t+\t.\tID=t1;Parent=g1\n\
            chr1\t.\texon\t1001\t2000\t.\t+\t.\tParent=t1\n\
            chr1\t.\texon\t8001\t9000\t.\t+\t.\tParent=t1\n\
            chr1\t.\texon\t300001\t301000\t.\t+\t.\tParent=t9\n\
            chr1\t.\texon\t330001\t331000\t.\t+\t.\tParent=t9\n";
        let data = crate::read::bgzf::fixture::blocks(text.as_bytes(), &[]);
        let end = ((data.len() as u64) - 31) << 16;
        // Generic columns counted from one, as `tabix -p gff` writes them.
        let index = crate::read::index::fixture::rooted(&["chr1"], [0, 1, 4, 5], (0, end));
        let dir = Scratch::new("tabix-parentless");
        let gff = dir.write("mixed.gff3.gz", &data);
        dir.write("mixed.gff3.gz.tbi", &index);
        for (place, through) in [("chr1:300,501-300,600", false), ("chr1:1-9,000", true)] {
            let line = format!("{place} {gff}");
            let (indexed, whole, notes) = watched(&line);
            fs::rename(format!("{gff}.tbi"), format!("{gff}.aside")).unwrap();
            let (from_whole, _, _) = watched(&line);
            fs::rename(format!("{gff}.aside"), format!("{gff}.tbi")).unwrap();
            assert!(indexed == from_whole, "{line}:\n{indexed}\n{from_whole}");
            assert!(notes.is_empty(), "{line}: {notes:?}");
            assert_eq!(whole.is_empty(), through, "{line}: read whole {whole:?}");
        }
    }

    /// An index whose counts add up to more than a count holds, which only a
    /// damaged one says, is read as one that counts nothing: added up, they
    /// stopped the program with an overflow, over a window and over an empty
    /// one alike. Each figure, and each refusal, is the whole file's.
    #[test]
    fn an_index_whose_counts_overflow_is_read_as_one_that_counts_none() {
        let dir = Scratch::new("tabix-overflow");
        let cohort = &indexed_into(&dir, &["cohort.vcf.gz"], &[])[0];
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/read/fixtures/indexed");
        let packed = fs::read(fixtures.join("cohort.vcf.gz.tbi")).unwrap();
        let mut tbi = crate::read::gzip::decompress(&packed).unwrap();
        // The pseudo-bin of each sequence, bin 37450 with two pairs, the
        // second of which opens with how many rows the sequence has.
        let pseudo = [0x4a, 0x92, 0, 0, 2, 0, 0, 0];
        let found: Vec<usize> = (0..tbi.len().saturating_sub(8))
            .filter(|at| tbi[*at..*at + 8] == pseudo)
            .collect();
        assert_eq!(found.len(), 3, "one for each sequence");
        for at in found {
            tbi[at + 24..at + 32].copy_from_slice(&(1u64 << 63).to_le_bytes());
        }
        dir.write("cohort.vcf.gz.tbi", &tbi);
        for place in ["chr1:400,001-460,000", "chr9:1-100"] {
            let line = format!("{place} --variants {cohort}");
            let (through, _, notes) = watched(&line);
            assert!(notes.is_empty(), "{line}: {notes:?}");
            fs::rename(format!("{cohort}.tbi"), format!("{cohort}.aside")).unwrap();
            let (from_whole, _, _) = watched(&line);
            fs::rename(format!("{cohort}.aside"), format!("{cohort}.tbi")).unwrap();
            assert_eq!(through, from_whole, "{line}");
        }
    }

    /// An index for another file, one that is not an index, and one beside a
    /// file compressed with gzip rather than bgzip are each read past: the
    /// file is read whole, and a note says why, once over two panels.
    #[test]
    fn an_index_that_does_not_fit_its_file_is_read_past_and_said() {
        let dir = Scratch::new("tabix-other");
        let paths = indexed_into(&dir, &["cohort.vcf.gz", "depth.bedgraph.gz"], &[]);
        let (cohort, depth) = (&paths[0], &paths[1]);
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/read/fixtures/indexed");
        let line = format!("chr1:1-2,000 chr1:400,001-460,000 --variants {cohort}");
        let (plain, _, _) = watched(&line);
        let index = |name: &str| fs::read(fixtures.join(name)).unwrap();
        dir.write("cohort.vcf.gz.tbi", &index("depth.bedgraph.gz.tbi"));
        let (drawn, whole, notes) = watched(&line);
        assert_eq!((drawn.as_str(), whole.len()), (plain.as_str(), 2));
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert!(
            notes[0].starts_with(&format!(
                "{cohort}.tbi does not describe {cohort}: the index"
            )),
            "{notes:?}"
        );
        assert!(
            notes[0].ends_with(&format!(
                "so the file was read whole; tabix -f -p bed {cohort} writes it again"
            )),
            "{notes:?}"
        );
        dir.write("cohort.vcf.gz.tbi", b"TBI\x01 and nothing an index holds");
        let (drawn, _, notes) = watched(&line);
        assert_eq!(drawn, plain);
        assert!(
            notes.len() == 1 && notes[0].contains("cannot be read as an index"),
            "{notes:?}"
        );
        // The rows of the bedGraph, compressed with gzip in one member, with
        // the index bgzip's blocks would have had beside it.
        let text = crate::read::gzip::decompress(&fs::read(depth).unwrap()).unwrap();
        let first: Vec<u8> = text
            .split_inclusive(|byte| *byte == b'\n')
            .take(400)
            .flatten()
            .copied()
            .collect();
        dir.write("depth.bedgraph.gz", &gzip_of(&first));
        dir.write("depth.bedgraph.gz.tbi", &index("depth.bedgraph.gz.tbi"));
        let (_, whole, notes) = watched(&format!("chr1:1-20,000 {depth}"));
        assert_eq!(whole, std::slice::from_ref(depth));
        assert_eq!(
            notes,
            [format!(
                "{depth} is compressed with gzip rather than bgzip, so {depth}.tbi beside it \
                 has no blocks to point to, and the file was read whole; gunzip it, bgzip it \
                 and index it again to read it a window at a time"
            )]
        );
    }

    /// Each BCF of the fixtures, with its CSI, beside the text `bcftools view`
    /// prints for it under the same name as a VCF.
    fn calls_held(with_index: bool) -> Held {
        use crate::read::bcf::fixture as bcf;
        let mut held = Held::new();
        for (name, bytes, index, view) in [
            ("cohort", bcf::COHORT, bcf::COHORT_CSI, bcf::COHORT_VIEW),
            ("tiny", bcf::TINY, bcf::TINY_CSI, bcf::TINY_VIEW),
            ("sv", bcf::SV, bcf::SV_CSI, bcf::SV_VIEW),
        ] {
            held.insert(format!("{name}.bcf"), bytes);
            held.insert(format!("{name}.vcf"), view);
            if with_index {
                held.insert(format!("{name}.bcf.csi"), index);
            }
        }
        held
    }

    /// What `line` draws from BCF and from the same calls as VCF, `.bcf`
    /// written as `.vcf` for the second, each with the notes it made, those
    /// of the VCF said of the BCF.
    fn both(held: &mut Held, line: &str) -> [(String, Vec<String>); 2] {
        [line.to_string(), line.replace(".bcf", ".vcf")].map(|line| {
            held.notes.clear();
            let drawn =
                held_figure(held, &line).unwrap_or_else(|error| format!("refused: {error}"));
            let notes = held
                .notes
                .iter()
                .map(|note| note.replace(".vcf", ".bcf"))
                .collect();
            (drawn.replace(".vcf", ".bcf"), notes)
        })
    }

    /// A BCF draws the figure the text `bcftools view` prints for it draws,
    /// byte for byte, through the CSI beside it and without one: its calls,
    /// its samples' genotypes, both of one file in one figure, structural
    /// calls whose mates lie outside the window, windows over every sequence
    /// with a record across two blocks, and a value of every type and every
    /// way of being missing. A cohort's calls named on their own say how many
    /// samples `--genotypes` draws, from the header, which the sites read for
    /// them do not name.
    #[test]
    fn a_bcf_draws_the_figure_its_vcf_draws() {
        for with_index in [true, false] {
            let mut held = calls_held(with_index);
            let mut drawn = 0;
            for window in [
                "chr1:1-200,000",
                "chr1:6,647-6,647",
                "chr1:1,000,000-1,400,000",
                "chr2:1-1,000,000",
                "chr3:50,000-60,000",
            ] {
                for tracks in [
                    "cohort.bcf",
                    "--genotypes cohort.bcf",
                    "--variants cohort.bcf --genotypes cohort.bcf --sample S3,S1",
                ] {
                    let [bcf, vcf] = both(&mut held, &format!("{window} {tracks}"));
                    assert_eq!(bcf, vcf, "{window} {tracks}");
                    drawn += usize::from(bcf.0.starts_with("<svg"));
                }
            }
            assert!(drawn >= 12, "{drawn}");
            // chr1:60,000 has `AF=.`, chr2:80,000 a GT left out of a sample
            // under `DP:GT`, which bcftools prints `4:.`, and chr2:400,000 an
            // END and an SVLEN left missing, each drawn as from the VCF.
            for line in [
                "chr1:1-100,000 tiny.bcf",
                "chr1:1-100,000 --genotypes tiny.bcf",
                "chr2:1-90,000 --genotypes tiny.bcf",
                "chr1:1-1,000,000 --structural sv.bcf",
                "chr2:1-500,000 --structural sv.bcf",
                "chr1:390,000-410,000 --structural sv.bcf",
            ] {
                let [bcf, vcf] = both(&mut held, line);
                assert!(bcf.0.starts_with("<svg"), "{line}: {}", bcf.0);
                assert_eq!(bcf, vcf, "{line}");
            }
            let [(_, notes), _] = both(&mut held, "chr1:1-200,000 cohort.bcf");
            assert_eq!(
                notes,
                ["cohort.bcf is drawn as its calls; --genotypes cohort.bcf draws its 3 samples, \
                  a row each"]
            );
        }
    }

    /// A BCF's empty window says what the file holds, from the counts its CSI
    /// keeps or from its records, as a VCF's does from its rows, and from the
    /// CSI alone where the header names a sequence no record is on; a place
    /// is a sequence its header names, as long as the header says, under the
    /// name `--rename` gives it, or as far as its records reach; and a name
    /// no file has is answered with the sequences the header names.
    #[test]
    fn a_bcf_is_placed_and_its_empty_windows_explained_as_its_vcf_is() {
        for with_index in [true, false] {
            let mut held = calls_held(with_index);
            for line in [
                "chr9:1-100 --variants cohort.bcf",
                "chr1:1-5 --genotypes cohort.bcf",
                "chr1:2,999,000-3,000,000 cohort.bcf",
                "chr2 cohort.bcf",
                "chrM tiny.bcf",
                "chr9 cohort.bcf",
                "chr1:20,000,000-20,001,000 --structural sv.bcf",
            ] {
                let [bcf, vcf] = both(&mut held, line);
                assert_eq!(bcf, vcf, "{line}");
            }
            let [(refused, _), _] = both(&mut held, "chr9:1-100 --variants cohort.bcf");
            assert_eq!(
                refused,
                "refused: --variants cohort.bcf: no variants in chr9:1-100, though the file \
                 holds 450 on chr1, chr2 and chr3"
            );
            let [(_, notes), _] = both(&mut held, "chrM tiny.bcf");
            assert!(
                notes[0].starts_with("chrM is drawn to 100, as far as tiny.bcf reaches"),
                "{notes:?}"
            );
            // A sequence placed by the name --rename gives it is as long as
            // the header says under its own name.
            let [bcf, vcf] = both(&mut held, "2 tiny.bcf --rename chr2=2");
            assert!(bcf.0.starts_with("<svg"), "{}", bcf.0);
            assert_eq!(bcf, vcf);
            assert!(
                !bcf.1.iter().any(|note| note.contains("as far as")),
                "{:?}",
                bcf.1
            );
        }
        // sv.bcf's header names scaffold_9, which no record is on and the CSI
        // gives no bins: the counts are the CSI's all the same, and a block
        // of records damaged, which every record read whole would stop at, is
        // never read for them.
        use crate::read::bcf::fixture::{SV, SV_CSI, SV_VIEW};
        let mut starts = Vec::new();
        let mut at = 0;
        while at + 18 <= SV.len() {
            starts.push(at);
            at += usize::from(u16::from_le_bytes([SV[at + 16], SV[at + 17]])) + 1;
        }
        let mut bent = SV.to_vec();
        bent[starts[starts.len() - 2] + 30] ^= 0xff;
        let mut held = Held::new();
        held.insert("sv.bcf", bent.clone());
        held.insert("sv.bcf.csi", SV_CSI);
        held.insert("sv.vcf", SV_VIEW);
        let [(bcf, _), (vcf, _)] = both(&mut held, "chr1:900,000-950,000 --variants sv.bcf");
        assert_eq!(
            bcf,
            "refused: --variants sv.bcf: no variants in chr1:900000-950000, though the file \
             holds 9 on chr1 and chr2"
        );
        assert_eq!(bcf, vcf);
        let mut whole = Held::new();
        whole.insert("sv.bcf", bent);
        assert!(held_figure(&mut whole, "chr1:1-1,000,000 --variants sv.bcf").is_err());
    }

    /// A BCF is read over the window through its CSI: a block far from the
    /// window, damaged, is never read through the index, and refuses the
    /// figure read without one.
    #[test]
    fn a_bcf_window_reads_only_the_blocks_its_csi_points_to() {
        use crate::read::bcf::fixture::{COHORT, COHORT_CSI};
        let mut starts = Vec::new();
        let mut at = 0;
        while at + 18 <= COHORT.len() {
            starts.push(at);
            at += usize::from(u16::from_le_bytes([COHORT[at + 16], COHORT[at + 17]])) + 1;
        }
        let mut bent = COHORT.to_vec();
        bent[starts[starts.len() - 2] + 30] ^= 0xff;
        let mut held = Held::new();
        held.insert("cohort.bcf", bent);
        held.insert("cohort.vcf", crate::read::bcf::fixture::COHORT_VIEW);
        held.insert("cohort.bcf.csi", COHORT_CSI);
        let [bcf, vcf] = both(&mut held, "chr1:1-200,000 cohort.bcf");
        assert_eq!(bcf, vcf);
        let mut whole = Held::new();
        whole.insert("cohort.bcf", held.files["cohort.bcf"].as_ref().to_vec());
        assert!(held_figure(&mut whole, "chr1:1-200,000 cohort.bcf").is_err());
    }

    /// An index beside a BCF that does not fit it is read past, the file read
    /// whole, which draws the same figure, and the note says why and the
    /// command that writes it again: the index of another file, one that
    /// does not read as an index, and a tabix index named as a CSI. An empty
    /// window counts the file's own records, and not the index's.
    #[test]
    fn a_csi_that_does_not_fit_its_bcf_is_read_past_and_said() {
        use crate::read::bcf::fixture::{COHORT, COHORT_VIEW, TINY_CSI};
        for (index, said) in [
            (
                TINY_CSI,
                "cohort.bcf.csi does not describe cohort.bcf: the index puts the first record",
            ),
            (
                &b"not an index"[..],
                "cohort.bcf.csi cannot be read as an index (not an index: it starts with none \
                 of BAI's, TBI's and CSI's magic), so cohort.bcf was read whole; bcftools \
                 index -f cohort.bcf writes it again",
            ),
            (
                &crate::read::index::fixture::ROWS_TBI[..],
                "cohort.bcf.csi is not the CSI bcftools index writes for a BCF, so cohort.bcf \
                 was read whole; bcftools index -f cohort.bcf writes it again",
            ),
        ] {
            let mut held = Held::new();
            held.insert("cohort.bcf", COHORT);
            held.insert("cohort.bcf.csi", index);
            held.insert("cohort.vcf", COHORT_VIEW);
            let [bcf, vcf] = both(&mut held, "chr1:1-200,000 --variants cohort.bcf");
            assert_eq!(bcf.0, vcf.0);
            assert_eq!(bcf.1.len(), 1, "{:?}", bcf.1);
            assert!(bcf.1[0].starts_with(said), "{:?}", bcf.1);
            assert!(bcf.1[0].ends_with("bcftools index -f cohort.bcf writes it again"));
            // An empty window counts the file's own records, and not the
            // records the index counts, which are another file's or none.
            let [(bcf, _), (vcf, _)] =
                both(&mut held, "chr1:2,999,000-3,000,000 --variants cohort.bcf");
            assert!(
                bcf.ends_with("though the file holds 450 on chr1, chr2 and chr3"),
                "{bcf}"
            );
            assert_eq!(bcf, vcf);
        }
    }

    /// A CSI older than its BCF is not trusted, as a tabix index older than
    /// its file is not, and the figure is the one the file draws read whole.
    #[cfg(unix)]
    #[test]
    fn a_csi_older_than_its_bcf_is_not_trusted() {
        use crate::read::bcf::fixture::{COHORT, COHORT_CSI};
        let dir = Scratch::new("bcf-older");
        let bcf = dir.write("cohort.bcf", COHORT);
        let csi = dir.write("cohort.bcf.csi", COHORT_CSI);
        let line = format!("chr1:1-200,000 --variants {bcf}");
        let (sound, _, notes) = watched(&line);
        assert!(sound.starts_with("<svg") && notes.is_empty(), "{notes:?}");
        made_older(std::slice::from_ref(&csi));
        let (drawn, _, notes) = watched(&line);
        assert_eq!(drawn, sound);
        assert_eq!(
            notes,
            [format!(
                "{csi} is older than {bcf}, so it was not trusted and the file was read \
                 whole; bcftools index -f {bcf} writes it again"
            )]
        );
    }

    /// A BCF handed to a track that draws no calls is refused naming the
    /// tracks that do, as is one given a `--format`; one on standard input
    /// is refused asking for its name; and one whose record a reader refuses
    /// names the record by its place, since it is no line of the file.
    #[test]
    fn a_bcf_where_its_calls_cannot_be_drawn_says_why() {
        let mut held = calls_held(true);
        for (line, said) in [
            (
                "chr1:1-100 --coverage cohort.bcf",
                "--coverage cohort.bcf: the file is BCF, calls along the sequence, which \
                 --variants, --genotypes and --structural draw",
            ),
            (
                "chr1:1-100 --coverage cohort.bcf --format bedgraph",
                "--coverage cohort.bcf: --format says what the columns of a text file are, and \
                 the file is BCF, which says what it holds itself",
            ),
        ] {
            assert_eq!(held_figure(&mut held, line).unwrap_err().to_string(), said);
        }
        let piped = build(&over("chr1:1-100", "--variants", "-"), |_: &Source| {
            decoded(crate::read::bcf::fixture::COHORT.to_vec(), None)
        })
        .unwrap_err()
        .to_string();
        assert_eq!(
            piped,
            "--variants standard input: the file is BCF, which karyon reads from a file named \
             on the command line and not from a pipe; name the file instead, or pipe in the \
             text bcftools view writes"
        );
        // Sample A's 0/1 at chr1:10 made 0/3, on a site of one alternate.
        let mut odd = crate::read::bcf::fixture::TINY_RAW.to_vec();
        let at = odd
            .windows(5)
            .position(|bytes| bytes == [0x11, 0x0d, 0x21, 0x02, 0x04])
            .unwrap();
        odd[at + 4] = 0x08;
        held.insert("odd.bcf", odd);
        assert_eq!(
            held_figure(&mut held, "chr1:1-100 --genotypes odd.bcf")
                .unwrap_err()
                .to_string(),
            "--genotypes odd.bcf: the record at chr1:10: sample A's GT 0/3 names allele 3, \
             and the line has 1 alternate allele"
        );
    }

    /// A BAM is read through the CSI `samtools index -c` writes as through
    /// its BAI, and the CSI is looked for first, as htslib looks: a BAI
    /// beside it that is not one is never read, by either name the CSI goes
    /// by, and a CSI that is not one is refused.
    #[test]
    fn a_bam_is_read_through_the_csi_beside_it_before_a_bai() {
        use crate::read::bam::fixture::{BAI, BAM, CSI};
        let line = "chr1:10-35 --coverage reads.bam --pileup reads.bam";
        let mut through_bai = Held::new();
        through_bai.insert("reads.bam", BAM);
        through_bai.insert("reads.bam.bai", BAI);
        let drawn = held_figure(&mut through_bai, line).unwrap();
        for name in ["reads.bam.csi", "reads.csi"] {
            let mut held = Held::new();
            held.insert("reads.bam", BAM);
            held.insert(name, CSI);
            held.insert("reads.bam.bai", "not an index");
            assert_eq!(held_figure(&mut held, line).unwrap(), drawn, "{name}");
        }
        let mut held = Held::new();
        held.insert("reads.bam", BAM);
        held.insert("reads.bam.csi", BAI);
        held.insert("reads.bam.bai", BAI);
        let error = held_figure(&mut held, line).unwrap_err().to_string();
        assert!(error.contains("reads.bam.csi is not the CSI"), "{error}");
        held.insert("reads.bam.csi", &CSI[..40]);
        assert!(held_figure(&mut held, line).is_err());
    }

    /// The figure the documentation's Start here page shows, drawn from
    /// `docs/data`, where `calls.vcf.gz` has its `.tbi` beside it: read
    /// through the index, it is the committed figure byte for byte.
    #[test]
    fn the_docs_reads_figure_is_the_same_through_calls_tbi() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/data");
        let at = |name: &str| data.join(name).display().to_string();
        let line = format!(
            "rpoB {} {} {} --width 720 --background #fbfaff",
            at("reads.bam"),
            at("genes.gff3"),
            at("calls.vcf.gz")
        );
        let (drawn, whole, notes) = watched(&line);
        assert!(!whole.contains(&at("calls.vcf.gz")), "{whole:?}");
        assert!(notes.is_empty(), "{notes:?}");
        let committed = fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/assets/start/reads.svg"),
        )
        .unwrap();
        assert!(drawn == committed, "the Start here figure changed");
    }

    /// Two genes, one either way round, each written as a gene over its CDS
    /// as NCBI writes a gene with no introns.
    const CODING: &str = "##gff-version 3
chr1\t.\tgene\t1001\t1300\t.\t+\t.\tID=gF;Name=geneF
chr1\t.\tCDS\t1001\t1300\t.\t+\t0\tParent=gF
chr1\t.\tgene\t2001\t2300\t.\t-\t.\tID=gR;Name=geneR
chr1\t.\tCDS\t2001\t2300\t.\t-\t0\tParent=gR
";

    /// Three kilobases of C with the first two codons of each gene written
    /// in: ATG and TGA, read forwards from base 1,001 and backwards from base
    /// 2,300, so the second codon is a stop or a tryptophan by the table.
    fn coding_reference() -> String {
        let mut bases = vec![b'C'; 3_000];
        bases[1_000..1_006].copy_from_slice(b"ATGTGA");
        bases[2_294..2_300].copy_from_slice(b"TCACAT");
        format!(">chr1\n{}\n", String::from_utf8(bases).unwrap())
    }

    /// The figure `line` draws from the two genes and their reference, and
    /// what it noted, less what the reference says of its own letters.
    fn coded(line: &str, genes: &str) -> (Result<String, BuildError>, Vec<String>) {
        let reference = coding_reference();
        let (drawn, notes) = drawn_noting(line, &[("genes.gff3", genes), ("ref.fa", &reference)]);
        let notes = notes
            .into_iter()
            .filter(|note| !note.starts_with("the bases are blocks of colour"))
            .collect();
        (drawn, notes)
    }

    #[test]
    fn codons_on_a_gene_placed_figure_count_from_its_start_codon() {
        let (drawn, notes) = coded("geneF genes.gff3 ref.fa --codons", CODING);
        let svg = drawn.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        assert!(
            svg.contains("<title>codon 1, 1,001 to 1,003, M</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>codon 2, 1,004 to 1,006, *</title>"));
        assert!(svg.contains("<title>codon 100, 1,298 to 1,300, P</title>"));
        assert!(!svg.contains("<title>codon 101,"), "a codon past the CDS");
        // Named for the gene, in the gutter: once more than the figure
        // without the ruler names it, in its title and on the gene.
        let without = coded("geneF genes.gff3 ref.fa", CODING).0.unwrap();
        assert_eq!(
            svg.matches(">geneF</text>").count(),
            without.matches(">geneF</text>").count() + 1,
            "{svg}"
        );
    }

    #[test]
    fn codons_on_a_reverse_gene_count_right_to_left() {
        let svg = coded("geneR genes.gff3 ref.fa --codons", CODING).0.unwrap();
        // Codon 1 is at the highest coordinate, read off the other strand.
        assert!(
            svg.contains("<title>codon 1, 2,298 to 2,300, M</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>codon 2, 2,295 to 2,297, *</title>"));
        assert!(svg.contains("<title>codon 100, 2,001 to 2,003, G</title>"));
    }

    #[test]
    fn codons_in_a_written_place_take_the_one_gene_coding_there() {
        let svg = coded("chr1:1101-1200 genes.gff3 ref.fa --codons", CODING)
            .0
            .unwrap();
        assert!(
            svg.contains("<title>codon 35, 1,103 to 1,105, P</title>"),
            "{svg}"
        );
        assert!(svg.contains(">geneF</text>"));
        // A codon half in the window has its letter, from the bases read
        // past either edge of it: base 1,100 on the left, 1,201 on the right.
        assert!(
            svg.contains("<title>codon 34, 1,100 to 1,102, P</title>"),
            "{svg}"
        );
        assert!(
            svg.contains("<title>codon 67, 1,199 to 1,201, P</title>"),
            "{svg}"
        );
    }

    /// A page that moves the window runs the command it was given over the
    /// window it moved to, and the ruler keeps counting the gene of the place
    /// written, though the window is now past the end of it.
    #[test]
    fn codons_keep_the_gene_of_the_place_written_wherever_the_window_goes() {
        let line: Vec<String> = "chr1:1101-1200 genes.gff3 --codons"
            .split_whitespace()
            .map(String::from)
            .collect();
        let Request::Draw(invocation) = parse(&line).unwrap() else {
            unreachable!("a figure")
        };
        let mut files = Held::new();
        files.insert("genes.gff3", CODING);
        let mut moved = |to: &str| {
            let to = Region::parse(to).unwrap();
            build_figure(
                &invocation,
                &mut files,
                |_, _| None,
                Theme::light(),
                Some(&to),
            )
            .unwrap()
            .figure
            .to_svg()
        };
        let svg = moved("chr1:1251-1350");
        assert!(
            svg.contains("<title>codon 100, 1,298 to 1,300</title>"),
            "{svg}"
        );
        // Moved onto the other gene, the ruler is still geneF's, and none of
        // it is in view.
        let svg = moved("chr1:2001-2100");
        assert!(!svg.contains("<title>codon"), "{svg}");
    }

    #[test]
    fn codons_in_a_place_two_genes_code_in_are_refused_naming_both() {
        let error = coded("chr1:1201-2100 genes.gff3 --codons", CODING)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::Several { genes, .. })
                if genes == &["geneF", "geneR"]),
            "{error:?}"
        );
        let said = error.to_string();
        assert!(said.starts_with("geneF and geneR code in chr1:"), "{said}");
        assert!(said.contains("as karyon geneF in place of"), "{said}");
        // The ruler first, as the annotation would refuse a place it has
        // nothing in before the ruler was asked.
        let none = coded("chr1:1501-1800 --codons genes.gff3", CODING)
            .0
            .unwrap_err();
        assert!(
            matches!(none, BuildError::Uncounted(CodonRefusal::NoneCodes { .. })),
            "{none:?}"
        );
        // A place on a gene's untranslated end is on a gene that has a CDS,
        // just not there.
        let leader = CODING.replacen("gene\t1001\t1300", "gene\t901\t1300", 1);
        let error = coded("chr1:911-990 --codons genes.gff3", &leader)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::NoneCodes { place, files })
                if place == "chr1:911-990" && files == &["genes.gff3"]),
            "{error:?}"
        );
    }

    #[test]
    fn codons_without_an_annotation_or_a_gene_are_refused() {
        let error = coded("chr1:1101-1200 ref.fa --codons", CODING)
            .0
            .unwrap_err();
        assert!(
            matches!(error, BuildError::Uncounted(CodonRefusal::NoAnnotation)),
            "{error:?}"
        );
        // A sequence drawn whole is a place and no gene.
        let error = coded("chr1 genes.gff3 ref.fa --codons", CODING)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::WholeSequence { sequence })
                if sequence == "chr1"),
            "{error:?}"
        );
    }

    /// A gene row alone says where a gene is and not where it starts coding,
    /// whether the figure is placed on it by name or over a stretch of it.
    #[test]
    fn codons_of_a_gene_with_no_cds_are_refused() {
        let genes = "##gff-version 3\nchr1\t.\tgene\t1001\t1300\t.\t+\t.\tID=gF;Name=geneF\n";
        for line in [
            "geneF genes.gff3 --codons",
            "chr1:1101-1200 genes.gff3 --codons",
        ] {
            let error = coded(line, genes).0.unwrap_err();
            assert!(
                matches!(&error, BuildError::Uncounted(CodonRefusal::NoCds { gene, files })
                    if gene == "geneF" && files == &["genes.gff3"]),
                "{line}: {error:?}"
            );
            assert!(error
                .to_string()
                .starts_with("geneF has no CDS in genes.gff3"));
        }
    }

    #[test]
    fn a_cds_split_by_an_intron_is_refused_as_spliced() {
        let genes = "##gff-version 3
chr1\t.\tgene\t1001\t2000\t.\t+\t.\tID=g1;Name=split
chr1\t.\tmRNA\t1001\t2000\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\texon\t1001\t1300\t.\t+\t.\tParent=t1
chr1\t.\texon\t1701\t2000\t.\t+\t.\tParent=t1
chr1\t.\tCDS\t1001\t1300\t.\t+\t0\tParent=t1
chr1\t.\tCDS\t1701\t1800\t.\t+\t0\tParent=t1
";
        let error = coded("split genes.gff3 --codons", genes).0.unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::Spliced { gene, pieces: 2 })
                if gene == "split"),
            "{error:?}"
        );
        assert!(error.to_string().contains("introns"), "{error}");
        // A BED12 says the same with its thick span across its blocks.
        let bed = "chr1\t1000\t2000\tsplit\t0\t+\t1000\t1800\t0\t2\t300,300,\t0,700,\n";
        let error = drawn_from("split genes.bed --codons", &[("genes.bed", bed)]).unwrap_err();
        assert!(
            matches!(
                error,
                BuildError::Uncounted(CodonRefusal::Spliced { pieces: 2, .. })
            ),
            "{error:?}"
        );
    }

    #[test]
    fn transcripts_that_code_different_stretches_are_refused_naming_them() {
        let genes = "##gff-version 3
chr1\t.\tgene\t1001\t1300\t.\t+\t.\tID=g1;Name=geneF
chr1\t.\tmRNA\t1001\t1300\t.\t+\t.\tID=t1;Name=geneF-201;Parent=g1
chr1\t.\tCDS\t1001\t1300\t.\t+\t0\tParent=t1
chr1\t.\tmRNA\t1001\t1300\t.\t+\t.\tID=t2;Name=geneF-202;Parent=g1
chr1\t.\tCDS\t1031\t1300\t.\t+\t0\tParent=t2
";
        let error = coded("geneF genes.gff3 --codons", genes).0.unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::Isoforms { gene, transcripts })
                if gene == "geneF" && transcripts == &["geneF-201", "geneF-202"]),
            "{error:?}"
        );
        // Placed on one of them by its name, it is counted.
        let svg = coded("geneF-202 genes.gff3 ref.fa --codons", genes)
            .0
            .unwrap();
        assert!(
            svg.contains("<title>codon 1, 1,031 to 1,033, P</title>"),
            "{svg}"
        );
    }

    #[test]
    fn codons_translate_with_the_table_the_cds_names_unless_told_otherwise() {
        let four = CODING.replacen("Parent=gF", "Parent=gF;transl_table=4", 1);
        let (drawn, notes) = coded("geneF genes.gff3 ref.fa --codons", &four);
        let svg = drawn.unwrap();
        assert!(notes.is_empty(), "{notes:?}");
        assert!(
            svg.contains("<title>codon 2, 1,004 to 1,006, W</title>"),
            "{svg}"
        );
        // The flag wins, and says what it overruled.
        let (drawn, notes) = coded("geneF genes.gff3 ref.fa --codons --genetic-code 1", &four);
        assert!(drawn
            .unwrap()
            .contains("<title>codon 2, 1,004 to 1,006, *</title>"));
        assert_eq!(
            notes,
            ["--genetic-code 1 reads geneF, whose CDS in genes.gff3 names table 4"]
        );
        // And reaches a gene whose CDS names none.
        let svg = coded("geneR genes.gff3 ref.fa --codons --genetic-code 2", CODING)
            .0
            .unwrap();
        assert!(
            svg.contains("<title>codon 2, 2,295 to 2,297, W</title>"),
            "{svg}"
        );
        // A table NCBI retired is refused rather than read as the standard one.
        let seven = CODING.replacen("Parent=gF", "Parent=gF;transl_table=7", 1);
        let error = coded("geneF genes.gff3 ref.fa --codons", &seven)
            .0
            .unwrap_err();
        assert!(
            matches!(
                error,
                BuildError::Uncounted(CodonRefusal::UnknownTable { table: 7, .. })
            ),
            "{error:?}"
        );
        assert!(
            coded("geneF genes.gff3 ref.fa --codons --genetic-code 4", &seven)
                .0
                .is_ok()
        );
    }

    #[test]
    fn a_cds_that_is_not_whole_codons_says_what_is_left_off() {
        let short = CODING.replacen("1300\t.\t+\t0\tParent=gF", "1299\t.\t+\t0\tParent=gF", 1);
        let (drawn, notes) = coded("geneF genes.gff3 --codons", &short);
        let svg = drawn.unwrap();
        assert_eq!(
            notes,
            [
                "the last 2 bases of the CDS of geneF are not a whole codon, and are left off \
              the ruler"
            ]
        );
        assert!(
            svg.contains("<title>codon 99, 1,295 to 1,297</title>"),
            "{svg}"
        );
        assert!(!svg.contains("<title>codon 100,"), "{svg}");
    }

    #[test]
    fn a_gene_on_no_strand_is_refused() {
        let unstranded = CODING.replace("\t+\t", "\t.\t");
        let error = coded("geneF genes.gff3 --codons", &unstranded)
            .0
            .unwrap_err();
        assert!(
            matches!(error, BuildError::Uncounted(CodonRefusal::NoStrand { .. })),
            "{error:?}"
        );
        assert!(
            error
                .to_string()
                .ends_with(": write + or - in the strand column of its rows"),
            "{error}"
        );
    }

    /// The table, the colour and the label each reach the ruler.
    #[test]
    fn a_codon_ruler_takes_its_label_and_its_colour() {
        let svg = coded(
            "geneF genes.gff3 --codons --label residues --color #d55e00",
            CODING,
        )
        .0
        .unwrap();
        assert!(svg.contains(">residues</text>"), "{svg}");
        assert!(svg.contains("#d55e00"), "{svg}");
    }

    /// An annotation that will not open, or will not read, is said as its
    /// own track says it, though the ruler is built first: the ruler's want
    /// of a gene is not what is wrong.
    #[test]
    fn a_codon_ruler_first_says_what_is_wrong_with_the_annotation() {
        let error = coded("chr1:1101-1200 --codons nowhere.gff3", CODING)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Open { track: "features", path, .. }
                if path == "nowhere.gff3"),
            "{error:?}"
        );
        let broken = "##gff-version 3\nchr1\t.\tgene\t1001\tlots\t.\t+\t.\tID=g\n";
        let error = coded("chr1:1101-1200 --codons genes.gff3", broken)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Parse { track: "features", path, .. }
                if path == "genes.gff3"),
            "{error:?}"
        );
    }

    /// A CDS that does not begin on a codon, by its phase or by NCBI's word
    /// that it goes on past its 5' end, has no start codon in the annotation
    /// to count from, and is refused on either strand rather than numbered
    /// from its first base: phase 1 here puts the first codon, ATG, at base
    /// 1,001, and on the reverse strand at 2,298, one base in from the end.
    #[test]
    fn a_cds_that_does_not_begin_on_its_start_codon_is_refused() {
        let phased = "##gff-version 3
chr1\t.\tgene\t1000\t1300\t.\t+\t.\tID=gF;Name=geneF
chr1\t.\tCDS\t1000\t1300\t.\t+\t1\tParent=gF;partial=true
chr1\t.\tgene\t2001\t2301\t.\t-\t.\tID=gR;Name=geneR
chr1\t.\tCDS\t2001\t2301\t.\t-\t1\tParent=gR;partial=true
";
        for line in [
            "geneF genes.gff3 ref.fa --codons",
            "chr1:1000-1020 genes.gff3 ref.fa --codons",
            "geneR genes.gff3 ref.fa --codons",
            "chr1:2280-2301 genes.gff3 ref.fa --codons",
        ] {
            let error = coded(line, phased).0.unwrap_err();
            assert!(
                matches!(
                    &error,
                    BuildError::Uncounted(CodonRefusal::Partial { phase: Some(1), .. })
                ),
                "{line}: {error:?}"
            );
            assert!(
                error
                    .to_string()
                    .contains("has phase 1, so its first base ends a codon that begins before it"),
                "{error}"
            );
        }
        // Phase 0 and a 5' end NCBI marks as going on past the contig.
        let ranged = CODING.replacen("Parent=gF", "Parent=gF;partial=true;start_range=.,1001", 1);
        let error = coded("geneF genes.gff3 ref.fa --codons", &ranged)
            .0
            .unwrap_err();
        assert!(
            matches!(
                &error,
                BuildError::Uncounted(CodonRefusal::Partial { phase: None, .. })
            ),
            "{error:?}"
        );
        let ranged = CODING.replacen("Parent=gR", "Parent=gR;partial=true;end_range=2300,.", 1);
        assert!(coded("geneR genes.gff3 ref.fa --codons", &ranged)
            .0
            .is_err());
        // A 3' end cut short still starts on its start codon.
        let open = CODING.replacen("Parent=gF", "Parent=gF;partial=true;end_range=1300,.", 1);
        let svg = coded("geneF genes.gff3 ref.fa --codons", &open).0.unwrap();
        assert!(
            svg.contains("<title>codon 1, 1,001 to 1,003, M</title>"),
            "{svg}"
        );
    }

    /// NCBI writes a ribosomal slippage as CDS rows that share a base.
    /// Joined into one stretch, it was counted straight through, and every
    /// codon past the slip was in the frame the ribosome left, with no word
    /// of it.
    #[test]
    fn cds_rows_that_shift_the_frame_are_refused_where_they_shift_it() {
        let slipped = "##gff-version 3
chr1\t.\tgene\t1001\t1600\t.\t+\t.\tID=gene-1;Name=gag-pol
chr1\t.\tCDS\t1001\t1150\t.\t+\t0\tID=cds-1;Parent=gene-1;exception=ribosomal slippage
chr1\t.\tCDS\t1150\t1600\t.\t+\t0\tID=cds-1;Parent=gene-1
chr1\t.\tgene\t1001\t1250\t.\t+\t.\tID=gene-2;Name=gag
chr1\t.\tCDS\t1001\t1250\t.\t+\t0\tParent=gene-2
";
        let error = coded("gag-pol genes.gff3 ref.fa --codons", slipped)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::Frameshift { gene, at: 1_149 })
                if gene == "gag-pol"),
            "{error:?}"
        );
        assert!(
            error.to_string().starts_with(
                "the CDS of gag-pol changes frame at 1,150, where its rows overlap or meet out \
                 of frame"
            ),
            "{error}"
        );
        // The gene beside it, which codes in one frame, is counted.
        let svg = coded("gag genes.gff3 ref.fa --codons", slipped).0.unwrap();
        assert!(svg.contains("<title>codon 1, 1,001 to 1,003, M</title>"));
        // On the reverse strand the rows are read from the right, and the
        // frame changes at the left one's 5' end.
        let reverse = "##gff-version 3
chr1\t.\tgene\t2001\t2300\t.\t-\t.\tID=g;Name=geneR
chr1\t.\tCDS\t2001\t2150\t.\t-\t0\tParent=g
chr1\t.\tCDS\t2150\t2300\t.\t-\t0\tParent=g
";
        let error = coded("geneR genes.gff3 ref.fa --codons", reverse)
            .0
            .unwrap_err();
        assert!(
            matches!(
                &error,
                BuildError::Uncounted(CodonRefusal::Frameshift { at: 2_149, .. })
            ),
            "{error:?}"
        );
        // Rows that meet are one frame where the phase carries it on, as a
        // GTF writes its stop codon after the CDS, and a shift where it does
        // not: 151 bases leave a codon two bases short.
        let met = |first: &str, phase: &str| {
            format!(
                "##gff-version 3
chr1\t.\tgene\t1001\t1300\t.\t+\t.\tID=gF;Name=geneF
chr1\t.\tCDS\t1001\t{first}\t.\t+\t0\tParent=gF
chr1\t.\tCDS\t{}\t1300\t.\t+\t{phase}\tParent=gF
",
                first.parse::<u64>().unwrap() + 1
            )
        };
        let svg = coded("geneF genes.gff3 ref.fa --codons", &met("1150", "0"))
            .0
            .unwrap();
        assert!(svg.contains("<title>codon 100, 1,298 to 1,300, P</title>"));
        assert!(coded("geneF genes.gff3 ref.fa --codons", &met("1151", "2"))
            .0
            .is_ok());
        let error = coded("geneF genes.gff3 ref.fa --codons", &met("1151", "0"))
            .0
            .unwrap_err();
        assert!(
            matches!(
                &error,
                BuildError::Uncounted(CodonRefusal::Frameshift { at: 1_151, .. })
            ),
            "{error:?}"
        );
        let gtf = "chr1\t.\tCDS\t1001\t1297\t.\t+\t0\tgene_id \"geneF\"; transcript_id \"t1\";
chr1\t.\tstart_codon\t1001\t1003\t.\t+\t0\tgene_id \"geneF\"; transcript_id \"t1\";
chr1\t.\tstop_codon\t1298\t1300\t.\t+\t0\tgene_id \"geneF\"; transcript_id \"t1\";
";
        let held = [("genes.gtf", gtf), ("ref.fa", &coding_reference())];
        let svg = drawn_from("geneF genes.gtf ref.fa --codons", &held).unwrap();
        assert!(svg.contains("<title>codon 100, 1,298 to 1,300, P</title>"));
    }

    /// GFF3 writes one feature in pieces as rows under one `ID=`, and a CDS
    /// under no gene written so is spliced, where the last row alone was
    /// counted from its own first base as codon 1.
    #[test]
    fn a_cds_written_in_rows_under_one_id_is_every_row_of_it() {
        let pieces = "##gff-version 3
chr1\t.\tCDS\t1001\t1100\t.\t+\t0\tID=c1;Name=abcA
chr1\t.\tCDS\t1201\t1300\t.\t+\t2\tID=c1;Name=abcA
";
        let error = coded("chr1:1050-1250 genes.gff3 ref.fa --codons", pieces)
            .0
            .unwrap_err();
        assert!(
            matches!(&error, BuildError::Uncounted(CodonRefusal::Spliced { gene, pieces: 2 })
                if gene == "abcA"),
            "{error:?}"
        );
    }

    /// A refusal is one line however wide the place: five names and how
    /// many more, as the files of the figure are named.
    #[test]
    fn a_refusal_names_five_genes_or_transcripts_and_counts_the_rest() {
        let mut genes = String::from("##gff-version 3\n");
        for at in 0..7 {
            let start = 1_001 + at * 100;
            genes.push_str(&format!(
                "chr1\t.\tCDS\t{start}\t{}\t.\t+\t0\tID=c{at};Name=gene{at}\n",
                start + 29
            ));
        }
        let error = coded("chr1:1001-2000 genes.gff3 --codons", &genes)
            .0
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with(
                "gene0, gene1, gene2, gene3, gene4 and 2 more genes code in chr1:1001-2000"
            ),
            "{error}"
        );
        assert!(!error.contains("gene6"), "{error}");
        let mut isoforms =
            String::from("##gff-version 3\nchr1\t.\tgene\t1001\t1300\t.\t+\t.\tID=g;Name=geneF\n");
        for at in 0..7 {
            isoforms.push_str(&format!(
                "chr1\t.\tmRNA\t1001\t1300\t.\t+\t.\tID=t{at};Parent=g\n\
                 chr1\t.\tCDS\t{}\t1300\t.\t+\t0\tParent=t{at}\n",
                1_001 + at * 3
            ));
        }
        let error = coded("geneF genes.gff3 --codons", &isoforms)
            .0
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with(
                "the transcripts of geneF, t0, t1, t2, t3, t4 and 2 more transcripts, code"
            ),
            "{error}"
        );
    }

    /// A place written over an annotation with a tabix index beside it is
    /// read through the index, as the annotation's own track reads it, and
    /// not whole, and is refused or drawn in the same words as the whole
    /// file draws it.
    #[test]
    fn codons_read_an_indexed_annotation_over_the_place_alone() {
        let dir = Scratch::new("tabix-codons");
        let genes = &indexed_into(&dir, &["genes.gff3.gz"], &[".tbi"])[0];
        let line = format!("chr1:9,001-9,100 {genes} --codons");
        let (through, whole, _) = watched(&line);
        assert!(whole.is_empty(), "read whole: {whole:?}");
        assert!(through.contains("chr1g1"), "{through}");
        fs::remove_file(format!("{genes}.tbi")).unwrap();
        let (read, whole, _) = watched(&line);
        assert!(whole.contains(genes), "{whole:?}");
        assert_eq!(through, read);
        // A CDS under no gene, in two rows under one ID 49 kb apart: the
        // window holds one of them, and the file is read whole for the other
        // rather than counting the one as the whole CDS. The ruler first,
        // since the features track draws the last row of an ID alone.
        let pieces = &indexed_into(&dir, &["pieces.gff3.gz"], &[".tbi"])[0];
        let line = format!("chr1:1,001-1,100 --codons {pieces}");
        let (through, _, _) = watched(&line);
        assert!(
            through.starts_with("refused: abcA codes in 2 pieces with introns between them"),
            "{through}"
        );
        fs::remove_file(format!("{pieces}.tbi")).unwrap();
        assert_eq!(through, watched(&line).0);
    }

    /// The documentation's own annotation writes rpoB's CDS, so its files
    /// draw the ruler the reads page offers: base 761,154 under codon 450,
    /// as the library's own example of the gene places it, and the table the
    /// CDS names, 11, read as the standard one.
    #[test]
    fn the_docs_genes_number_rpob_with_s450_at_base_761_154() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/data");
        let at = |name: &str| data.join(name).display().to_string();
        let svg = drawn_from_disk(&format!(
            "NC_000962.3:761,081-761,200 {} {} --codons",
            at("genes.gff3"),
            at("ref.fa")
        ))
        .unwrap();
        let library = crate::CodonTrack::new(759_806, 763_325, crate::Strand::Forward);
        assert_eq!(library.codon_of(761_153), Some(450));
        assert!(
            svg.contains("<title>codon 450, 761,154 to 761,156, "),
            "{svg}"
        );
        assert!(svg.contains(">rpoB</text>"), "{svg}");
        let text = fs::read_to_string(at("genes.gff3")).unwrap();
        assert_eq!(
            read::interval::translation_table(&text, "NC_000962.3", 761_153, 761_156),
            Some(11)
        );
    }

    /// A chromosome of a hundred kilobases with a gene on each strand, a
    /// GFF3 that says how long it is, calls whose header says so too, and a
    /// depth in windows of a kilobase that loses a stretch and doubles one.
    fn circle_files() -> Vec<(&'static str, String)> {
        let genes = "##gff-version 3\n##sequence-region chrC 1 100000\n\
                     chrC\t.\tgene\t1001\t3000\t.\t+\t.\tName=gA\n\
                     chrC\t.\tgene\t50001\t52000\t.\t-\t.\tName=gB\n"
            .to_string();
        let calls = "##fileformat=VCFv4.2\n##contig=<ID=chrC,length=100000>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chrC\t1500\t.\tC\tT\t60\tPASS\tANN=T|missense_variant|MODERATE|gA\n\
                     chrC\t2500\t.\tC\tT\t60\tPASS\tANN=T|synonymous_variant|LOW|gA\n"
            .to_string();
        let depth = |low: f64, high: f64| -> String {
            (0..100u64)
                .map(|window| {
                    let value = match window {
                        20..=29 => low,
                        60..=69 => high,
                        _ => 30.0,
                    };
                    format!(
                        "chrC\t{}\t{}\t{value}\n",
                        window * 1000,
                        (window + 1) * 1000
                    )
                })
                .collect()
        };
        vec![
            ("genes.gff3", genes),
            ("calls.vcf", calls),
            ("depth.bg", depth(0.0, 60.0)),
            ("deeper.bg", depth(0.0, 300.0)),
        ]
    }

    /// A command line drawn as a circle over [`circle_files`] and `more`.
    fn circle(line: &str, more: &[(&str, &str)]) -> Result<String, BuildError> {
        let mut held = Held::new();
        for (name, text) in circle_files() {
            held.insert(name, text);
        }
        for (name, text) in more {
            held.insert(*name, *text);
        }
        build_files(&invocation(line), &mut held, |_, _| None)
    }

    /// The rings a command line's circle is made of, outside in.
    fn rings_of(line: &str) -> Vec<Ringed> {
        let mut held = Held::new();
        for (name, text) in circle_files() {
            held.insert(name, text);
        }
        let invocation = invocation(line);
        let region = whole_sequence(&invocation, &mut held).unwrap();
        circle_rings(&invocation, &region, &Theme::light(), &mut held)
            .unwrap_or_else(|error| panic!("{line}: {error}"))
    }

    /// Where each of `needles` first appears in `svg`, in the order given.
    fn first_at(svg: &str, needles: &[&str]) -> Vec<usize> {
        needles
            .iter()
            .map(|needle| {
                svg.find(needle)
                    .unwrap_or_else(|| panic!("no {needle} in {svg}"))
            })
            .collect()
    }

    /// Each track is a ring, in the order written, the first outermost
    /// inside the ruler, and each ring is one tooltip under its name.
    #[test]
    fn a_circle_draws_one_ring_per_track_outermost_first() {
        let svg = circle("chrC --circular genes.gff3 calls.vcf depth.bg", &[]).unwrap();
        let at = first_at(
            &svg,
            &[
                "<g><title>genes</title>",
                "<g><title>calls</title>",
                "<g><title>depth</title>",
            ],
        );
        assert!(at[0] < at[1] && at[1] < at[2], "{at:?}");
        // Drawn in that order from the outside in, which the rings say for
        // themselves: the ruler's circle is the widest, then the annotation.
        let rings = rings_of("chrC --circular genes.gff3 calls.vcf depth.bg");
        let kinds: Vec<Kind> = rings.iter().map(|ringed| ringed.kind).collect();
        assert_eq!(
            kinds,
            [Kind::Axis, Kind::Features, Kind::Variants, Kind::Coverage]
        );
        assert!(svg.contains("<title id=\"karyon-title\">chrC, 100,000 bases</title>"));
        // A ruler written among the tracks is drawn there, and none outside.
        let rings = rings_of("chrC --circular genes.gff3 --axis calls.vcf");
        let kinds: Vec<Kind> = rings.iter().map(|ringed| ringed.kind).collect();
        assert_eq!(kinds, [Kind::Features, Kind::Axis, Kind::Variants]);
        let rings = rings_of("chrC --circular genes.gff3 --no-axis");
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].kind, Kind::Features);
    }

    /// The key under the circle names each ring outside in, a line each, with
    /// what its colours mean, and `--no-legend` leaves it out.
    #[test]
    fn the_key_under_a_circle_names_each_ring_outside_in() {
        let svg = circle("chrC --circular genes.gff3 calls.vcf depth.bg", &[]).unwrap();
        let line = |name: &str| -> f64 {
            let tail = format!(">{name}</text>");
            let end = svg
                .find(&tail)
                .unwrap_or_else(|| panic!("no key line {name}"));
            let start = svg[..end].rfind("<text").unwrap();
            let y = svg[start..end].split(" y=\"").nth(1).unwrap();
            y.split('"').next().unwrap().parse().unwrap()
        };
        let (genes, calls, depth) = (line("genes"), line("calls"), line("depth"));
        assert!(genes < calls && calls < depth, "{genes} {calls} {depth}");
        // Under the square the circle is drawn in, which the image grows to hold.
        assert!(genes > 668.0, "{genes}");
        for said in [
            "forward strand, outer half",
            "reverse strand, inner half",
            "missense_variant",
            "synonymous_variant",
            "above its median, 30",
            "below it",
        ] {
            assert!(svg.contains(&format!(">{said}</text>")), "{said}");
        }
        let bare = circle(
            "chrC --circular genes.gff3 calls.vcf depth.bg --no-legend",
            &[],
        )
        .unwrap();
        assert!(!bare.contains(">forward strand, outer half</text>"));
        assert!(bare.contains("height=\"668\""), "{}", &bare[..200]);
    }

    /// A depth is cut into arcs, neighbouring arcs of one value drawn as one,
    /// and read either side of its median: the stretch lost lies inside the
    /// line and the one carried twice outside it.
    #[test]
    fn a_ring_of_coverage_is_binned_and_read_either_side_of_its_median() {
        let rings = rings_of("chrC:1-100,000 --circular depth.bg");
        let RingOf::Signal(ring) = &rings[1].ring else {
            panic!("a depth is a signal ring");
        };
        assert_eq!(ring.baseline_value(), 30.0);
        let spans: Vec<(u64, u64, f64)> = ring
            .windows()
            .iter()
            .map(|window| (window.start, window.end, window.value))
            .collect();
        assert_eq!(
            spans,
            [
                (0, 20_000, 30.0),
                (20_000, 30_000, 0.0),
                (30_000, 60_000, 30.0),
                (60_000, 70_000, 60.0),
                (70_000, 100_000, 30.0),
            ]
        );
        // A depth file of a line a base is still at most a thousand arcs.
        let per_base: String = (0..100_000u64)
            .map(|pos| format!("chrC\t{}\t{}\n", pos + 1, 30 + pos % 7))
            .collect();
        let mut held = Held::new();
        held.insert("calls.vcf", circle_files()[1].1.clone());
        held.insert("base.depth", per_base);
        let invocation =
            invocation("chrC --circular calls.vcf --coverage base.depth --aggregate min");
        let region = whole_sequence(&invocation, &mut held).unwrap();
        let rings = circle_rings(&invocation, &region, &Theme::light(), &mut held).unwrap();
        let RingOf::Signal(ring) = &rings[2].ring else {
            panic!("a depth is a signal ring");
        };
        assert!(ring.windows().len() <= ARCS, "{}", ring.windows().len());
        assert!(ring.windows().iter().all(|window| window.value == 30.0));
    }

    /// `--same-scale` gives every ring of depth one reach, as it puts every
    /// band of depth on one ceiling; without it each reaches its own.
    #[test]
    fn same_scale_gives_every_coverage_ring_one_reach() {
        let reaches = |line: &str| -> Vec<f64> {
            rings_of(line)
                .iter()
                .filter_map(|ringed| match &ringed.ring {
                    RingOf::Signal(ring) => Some(ring.reach()),
                    RingOf::Other(_) => None,
                })
                .collect()
        };
        let apart = reaches("chrC:1-100,000 --circular depth.bg deeper.bg");
        assert!(apart[0] < apart[1], "{apart:?}");
        let together = reaches("chrC:1-100,000 --circular depth.bg deeper.bg --same-scale");
        assert_eq!(together[0], together[1]);
        assert_eq!(together[1], apart[1]);
    }

    /// A circle closes where its sequence ends, so a sequence no file gives
    /// the length of is refused rather than closed where a bedGraph's last
    /// row happens to stop; written from base 1, the length is the one
    /// written.
    #[test]
    fn a_circle_without_a_stated_length_is_refused_rather_than_closed_early() {
        let refused = circle("chrC --circular depth.bg", &[]).unwrap_err();
        assert!(
            matches!(
                &refused,
                BuildError::Uncircled(CircleRefusal::NoLength { sequence, reach: 100_000, .. })
                    if sequence == "chrC"
            ),
            "{refused:?}"
        );
        let said = refused.to_string();
        assert!(
            said.contains("##contig") && said.contains("chrC:1-LENGTH"),
            "{said}"
        );
        let drawn = circle("chrC:1-100,000 --circular depth.bg", &[]).unwrap();
        assert!(drawn.contains("chrC, 100,000 bases"));
    }

    /// A span written from base 1 is the whole of a sequence that long, and
    /// refused where a file says the sequence is another length.
    #[test]
    fn a_written_span_the_files_disagree_with_is_refused() {
        let refused = circle("chrC:1-50,000 --circular genes.gff3", &[]).unwrap_err();
        assert!(
            matches!(
                refused,
                BuildError::Uncircled(CircleRefusal::Length {
                    written: 50_000,
                    stated: 100_000,
                    ..
                })
            ),
            "{refused:?}"
        );
        assert!(circle("chrC:1-100,000 --circular genes.gff3", &[]).is_ok());
    }

    /// Two files that give a sequence different lengths are refused, each
    /// named with its length, whichever is written first and whether the
    /// sequence is named or written from base 1. The circle was closed at the
    /// first file's length, so the same files drew it at 50 kb one way round
    /// and 100 kb the other, and a span was checked against the first alone.
    #[test]
    fn a_circle_the_files_give_two_lengths_is_refused_in_either_order() {
        let short = "##fileformat=VCFv4.2\n##contig=<ID=chrC,length=50000>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chrC\t1500\t.\tC\tT\t60\tPASS\t.\n";
        let short_first = [("short.vcf", 50_000), ("genes.gff3", 100_000)];
        let genes_first = [("genes.gff3", 100_000), ("short.vcf", 50_000)];
        for (line, order) in [
            ("chrC --circular short.vcf genes.gff3 depth.bg", short_first),
            ("chrC --circular genes.gff3 short.vcf depth.bg", genes_first),
            ("chrC:1-50,000 --circular short.vcf genes.gff3", short_first),
            ("chrC:1-50,000 --circular genes.gff3 short.vcf", genes_first),
            (
                "chrC:1-100,000 --circular short.vcf genes.gff3",
                short_first,
            ),
        ] {
            let refused = circle(line, &[("short.vcf", short)]).unwrap_err();
            let said: Vec<(String, u64)> = order
                .iter()
                .map(|(file, length)| (file.to_string(), *length))
                .collect();
            assert!(
                matches!(
                    &refused,
                    BuildError::Uncircled(CircleRefusal::Lengths { sequence, said: got })
                        if sequence == "chrC" && *got == said
                ),
                "{line}: {refused:?}"
            );
        }
        let refused = circle(
            "chrC --circular short.vcf genes.gff3",
            &[("short.vcf", short)],
        )
        .unwrap_err();
        assert_eq!(
            refused.to_string(),
            "a circle closes where chrC ends, and the files disagree on where that is: \
             short.vcf says 50,000 bases and genes.gff3 says 100,000 bases; draw it from \
             files that agree on how long chrC is"
        );
        // Files that agree are one circle, however many of them say so.
        let agreeing = short.replace("50000", "100000");
        let svg = circle(
            "chrC --circular agree.vcf genes.gff3 calls.vcf",
            &[("agree.vcf", &agreeing)],
        )
        .unwrap();
        assert!(svg.contains("chrC, 100,000 bases"), "{svg}");
    }

    /// A consequence is the colour its rank gives it, not the colour of its
    /// place in the window. Dealt by first appearance, missense was the
    /// second colour across the whole gene, whose first call is synonymous,
    /// and the first over a window of missense calls alone.
    #[test]
    fn a_consequence_keeps_its_colour_in_a_zoom_into_its_gene() {
        let calls = "##fileformat=VCFv4.2\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chr1\t100\t.\tC\tT\t.\t.\tANN=T|synonymous_variant|LOW|g\n\
                     chr1\t500\t.\tC\tT\t.\t.\tANN=T|missense_variant|MODERATE|g\n\
                     chr1\t900\t.\tC\tT\t.\t.\tANN=T|synonymous_variant|LOW|g\n";
        let colour = |line: &str, category: &str| -> String {
            let mut held = Held::new();
            held.insert("calls.vcf", calls);
            let svg = held_figure(&mut held, line).unwrap();
            let rest = svg
                .split(&format!("<title>{category}, colour "))
                .nth(1)
                .unwrap_or_else(|| panic!("no {category} in the key of {line}: {svg}"));
            rest[..rest.find("</title>").unwrap()].to_string()
        };
        let theme = Theme::default();
        let whole = colour("chr1:1-1000 calls.vcf", "missense_variant");
        assert_eq!(whole, theme.color(0));
        for zoom in ["chr1:450-550", "chr1:50-600", "chr1:450-1000"] {
            assert_eq!(
                colour(&format!("{zoom} calls.vcf"), "missense_variant"),
                whole
            );
        }
        assert_eq!(
            colour("chr1:1-1000 calls.vcf", "synonymous_variant"),
            colour("chr1:50-600 calls.vcf", "synonymous_variant")
        );
    }

    /// Each line of the key is the colour its ring paints what it names: a
    /// gene on each strand, and a call of each consequence, the first named
    /// again after the second. The key's colours and the ring's are worked
    /// out apart, so a key that named the reverse strand in the forward
    /// strand's colour drew a figure every other test passed.
    #[test]
    fn the_key_under_a_circle_is_the_colour_of_what_it_names() {
        let calls = "##fileformat=VCFv4.2\n##contig=<ID=chrC,length=100000>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chrC\t1500\t.\tC\tT\t60\tPASS\tANN=T|missense_variant|MODERATE|gA\n\
                     chrC\t2500\t.\tC\tT\t60\tPASS\tANN=T|synonymous_variant|LOW|gA\n\
                     chrC\t3500\t.\tC\tT\t60\tPASS\tANN=T|missense_variant|MODERATE|gA\n";
        for theme in ["light", "dark"] {
            let line = format!("chrC --circular genes.gff3 three.vcf --theme {theme}");
            let svg = circle(&line, &[("three.vcf", calls)]).unwrap();
            // The colour a line of the key says it is.
            let key = |said: &str| -> String {
                let rest = svg
                    .split(&format!("<title>{said}, colour "))
                    .nth(1)
                    .unwrap_or_else(|| panic!("{theme}: no key line {said}: {svg}"));
                rest[..rest.find("</title>").unwrap()].to_string()
            };
            // The colour a feature is painted in.
            let painted = |title: &str| -> String {
                let rest = svg
                    .split(&format!("<title>{title}"))
                    .nth(1)
                    .unwrap_or_else(|| panic!("{theme}: no {title}: {svg}"));
                let rest = rest.split(" fill=\"").nth(1).unwrap();
                rest[..rest.find('"').unwrap()].to_string()
            };
            assert_eq!(
                painted("gA, "),
                key("forward strand, outer half"),
                "{theme}"
            );
            assert_eq!(
                painted("gB, "),
                key("reverse strand, inner half"),
                "{theme}"
            );
            let ring = svg.split("<g><title>three</title>").nth(1).unwrap();
            let ring = &ring[..ring.find("</g>").unwrap()];
            let ticks: Vec<String> = ring
                .split(" stroke=\"")
                .skip(1)
                .map(|rest| rest[..rest.find('"').unwrap()].to_string())
                .collect();
            assert_eq!(
                ticks,
                [
                    key("missense_variant"),
                    key("synonymous_variant"),
                    key("missense_variant")
                ],
                "{theme}"
            );
            assert_ne!(ticks[0], ticks[1], "{theme}");
            assert_ne!(painted("gA, "), painted("gB, "), "{theme}");
        }
    }

    /// The middle of the circle names the sequence and says how long it is,
    /// as a figure's locus does at its top right; a title takes the name's
    /// place and keeps both under it, and `--no-region-label` leaves them out.
    #[test]
    fn the_middle_names_the_sequence_and_its_length_or_the_title() {
        let middle = |line: &str| -> Vec<String> {
            let svg = circle(line, &[]).unwrap();
            svg.split("text-anchor=\"middle\"")
                .skip(1)
                .filter_map(|rest| {
                    let text = rest.split('>').nth(1)?.split('<').next()?;
                    (!rest.starts_with(" aria-hidden")).then(|| text.to_string())
                })
                .filter(|text| text.contains("chrC") || text.contains("bases") || text == "Plasmid")
                .collect()
        };
        assert_eq!(
            middle("chrC --circular genes.gff3"),
            ["chrC", "100,000 bases"]
        );
        assert_eq!(
            middle("chrC --circular genes.gff3 --title Plasmid"),
            ["Plasmid", "chrC, 100,000 bases"]
        );
        assert!(middle("chrC --circular genes.gff3 --no-region-label").is_empty());
    }

    /// A gene's name places a figure on part of a sequence, and a circle is
    /// all of one.
    #[test]
    fn a_gene_is_not_a_circle() {
        let refused = circle("gA --circular genes.gff3", &[]).unwrap_err();
        assert_eq!(
            refused.to_string(),
            "gA is a gene, and a circle is a whole sequence: name the sequence it is on, as \
             karyon chrC --circular"
        );
    }

    /// A FASTA is drawn as its GC skew, in windows of a thousandth of the
    /// sequence, and says how long the sequence is.
    #[test]
    fn a_fasta_on_a_circle_is_its_gc_skew_and_its_length() {
        let bases: String = (0..100_000)
            .map(|at| if at < 50_000 { 'G' } else { 'C' })
            .collect();
        let fasta = format!(">chrC\n{bases}\n");
        let mut held = Held::new();
        held.insert("ref.fa", fasta);
        let invocation = invocation("chrC --circular ref.fa");
        let region = whole_sequence(&invocation, &mut held).unwrap();
        assert_eq!(region.end(), 100_000);
        let rings = circle_rings(&invocation, &region, &Theme::light(), &mut held).unwrap();
        assert_eq!(rings[1].name.as_deref(), Some("ref GC skew"));
        let RingOf::Signal(ring) = &rings[1].ring else {
            panic!("a reference is a signal ring");
        };
        assert_eq!(ring.windows().len(), 1_000);
        assert_eq!(ring.windows()[0].end - ring.windows()[0].start, 100);
        assert!(ring.windows()[..500]
            .iter()
            .all(|window| window.value == 1.0));
        assert!(ring.windows()[500..]
            .iter()
            .all(|window| window.value == -1.0));
    }

    /// A deletion is an arc on its ring, and a breakend join on the sequence
    /// a chord across the middle; a join to another sequence has nowhere to
    /// go on a circle of one.
    #[test]
    fn a_breakend_join_is_a_chord_and_a_deletion_is_an_arc() {
        let calls = "##fileformat=VCFv4.2\n##contig=<ID=chrC,length=100000>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chrC\t20000\tdel1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=25000\n\
                     chrC\t40000\tbnd1\tN\tN[chrC:70000[\t60\tPASS\tSVTYPE=BND\n\
                     chrC\t70000\tbnd2\tN\t]chrC:40000]N\t60\tPASS\tSVTYPE=BND\n\
                     chrC\t80000\tbnd3\tN\tN[chrD:100[\t60\tPASS\tSVTYPE=BND\n";
        let svg = circle("chrC --circular --structural sv.vcf", &[("sv.vcf", calls)]).unwrap();
        assert_eq!(svg.matches("<title>link, source").count(), 1, "{svg}");
        // The mate's position, 70,000, is the last base of the join, as the
        // file writes it, not the base after it.
        assert!(
            svg.contains("<title>link, source 40,001, target 70,000</title>"),
            "{svg}"
        );
        // And the chord names the two bases the band of the same file does.
        let band = circle("chrC:1-100,000 --structural sv.vcf", &[("sv.vcf", calls)]).unwrap();
        let said = band
            .split("<title>translocation, ")
            .nth(1)
            .unwrap_or_else(|| panic!("no join on the band: {band}"));
        let said = &said[..said.find("</title>").unwrap()];
        let (first, last) = said.split_once(" to ").unwrap();
        assert!(
            svg.contains(&format!(
                "<title>link, source {first}, target {last}</title>"
            )),
            "the band says {said}: {svg}"
        );
        assert!(
            svg.contains("<title>deletion del1, 20,001 to 25,000"),
            "{svg}"
        );
        assert!(
            svg.contains(">translocation, across the middle</text>"),
            "{svg}"
        );
        assert!(svg.contains(">deletion</text>"), "{svg}");
        // A join is a ribbon a pixel and a half wide, drawn nearly opaque.
        let chord = svg.split("<title>link, source").nth(1).unwrap();
        let chord = &chord[..chord.find("</g>").unwrap()];
        assert!(chord.contains("fill-opacity=\"0.8\""), "{chord}");
    }

    /// A ring of windows is each arc's mean, read either side of 0 as a band
    /// of windows is, and a stretch no window covers is no arc rather than
    /// an arc of nought.
    #[test]
    fn a_ring_of_windows_is_their_mean_and_a_gap_is_no_arc() {
        let mut held = Held::new();
        held.insert("calls.vcf", circle_files()[1].1.clone());
        // Windows of 50 bases, alternating 1 and -1, over the first half.
        let windows: String = (0..1_000u64)
            .map(|at| {
                let value = if at % 2 == 0 { 1.0 } else { -1.0 };
                format!("chrC\t{}\t{}\t{value}\n", at * 50, at * 50 + 50)
            })
            .collect();
        held.insert("skew.bg", windows);
        let invocation = invocation("chrC --circular calls.vcf --windows skew.bg");
        let region = whole_sequence(&invocation, &mut held).unwrap();
        let rings = circle_rings(&invocation, &region, &Theme::light(), &mut held).unwrap();
        let RingOf::Signal(ring) = &rings[2].ring else {
            panic!("windows are a signal ring");
        };
        assert_eq!(ring.baseline_value(), 0.0);
        // Two windows to an arc of a hundred bases, one of each sign: the
        // mean, not the larger, and nothing past the last window.
        assert_eq!(ring.windows(), [crate::Window::new(0, 50_000, 0.0)]);
        let said: Vec<&str> = rings[2]
            .legend
            .items()
            .iter()
            .map(|item| match item {
                crate::LegendItem::Key { label, .. } => label.as_str(),
                crate::LegendItem::Ramp { label, .. } => label.as_str(),
            })
            .collect();
        assert_eq!(said, ["above 0", "below 0"]);
    }

    /// Calls are coloured by consequence from the most damaging down, as a
    /// band of calls deals them, a call no annotator named by its shape, as
    /// the reader names it, after every consequence.
    #[test]
    fn calls_on_a_ring_are_coloured_by_consequence_from_the_most_damaging_down() {
        let theme = Theme::light();
        let keyed = |line: &str, more: &[(&str, &str)]| -> Vec<(String, String)> {
            let mut held = Held::new();
            for (name, text) in circle_files() {
                held.insert(name, text);
            }
            for (name, text) in more {
                held.insert(*name, *text);
            }
            let invocation = invocation(line);
            let region = whole_sequence(&invocation, &mut held).unwrap();
            let rings = circle_rings(&invocation, &region, &theme, &mut held).unwrap();
            rings
                .last()
                .unwrap()
                .legend
                .items()
                .iter()
                .map(|item| match item {
                    crate::LegendItem::Key { label, color, .. } => (label.clone(), color.clone()),
                    crate::LegendItem::Ramp { label, .. } => (label.clone(), String::new()),
                })
                .collect()
        };
        assert_eq!(
            keyed("chrC --circular calls.vcf", &[]),
            [
                ("missense_variant".to_string(), theme.color(0).to_string()),
                ("synonymous_variant".to_string(), theme.color(1).to_string()),
            ]
        );
        let mixed = "##contig=<ID=chrC,length=100000>\n\
                     chrC\t100\t.\tC\tT\t.\t.\t.\n\
                     chrC\t200\t.\tC\tT\t.\t.\tANN=T|stop_gained|HIGH|gA\n";
        assert_eq!(
            keyed(
                "chrC --circular --variants mixed.vcf",
                &[("mixed.vcf", mixed)]
            ),
            [
                ("stop_gained".to_string(), theme.color(0).to_string()),
                ("substitution".to_string(), theme.color(1).to_string()),
            ]
        );
    }

    /// Names are written on a ring of a few named features and left off one
    /// of many, which would be a wheel of unreadable text.
    #[test]
    fn a_ring_of_many_named_features_leaves_their_names_off() {
        let genes = |count: u64| -> String {
            let mut text = "##sequence-region chrC 1 100000\n".to_string();
            for at in 0..count {
                text.push_str(&format!(
                    "chrC\t.\tgene\t{}\t{}\t.\t+\t.\tName=g{at}\n",
                    at * 4_000 + 1,
                    at * 4_000 + 2_000
                ));
            }
            text
        };
        let few = genes(NAMED_ON_A_RING as u64);
        let many = genes(NAMED_ON_A_RING as u64 + 1);
        let svg = circle("chrC --circular few.gff3", &[("few.gff3", &few)]).unwrap();
        assert!(svg.contains(">g0</text>"), "{svg}");
        let svg = circle("chrC --circular many.gff3", &[("many.gff3", &many)]).unwrap();
        assert!(!svg.contains(">g0</text>"), "{svg}");
        assert!(svg.contains("<title>g0, "), "each gene keeps its tooltip");
    }

    /// A file that calls the sequence by another name is read under the one
    /// `--rename` gives it, as a figure along the sequence reads it.
    #[test]
    fn a_ring_reads_a_file_under_the_name_rename_gives_it() {
        let depth = "1\t0\t50000\t30\n1\t50000\t100000\t60\n";
        let svg = circle(
            "chrC --circular calls.vcf --coverage one.bg --rename 1=chrC",
            &[("one.bg", depth)],
        )
        .unwrap();
        assert!(svg.contains("<g><title>one</title>"), "{svg}");
        assert!(circle(
            "chrC --circular calls.vcf --coverage one.bg",
            &[("one.bg", depth)]
        )
        .is_err());
    }

    /// A `.bed` named on its own is told by what it holds, and one that holds
    /// methylation has no ring; one that holds depth in windows does.
    #[test]
    fn a_file_that_turns_out_to_have_no_ring_is_refused_by_what_it_holds() {
        let methylation =
            "chrC\t100\t101\tm\t20\t+\t100\t101\t255,0,0\t20\t85.0\t17\t3\t0\t0\t0\t0\t0\n";
        let refused = circle(
            "chrC --circular calls.vcf calls.bed",
            &[("calls.bed", methylation)],
        )
        .unwrap_err();
        assert!(
            matches!(
                &refused,
                BuildError::Uncircled(CircleRefusal::NoRing {
                    kind: Kind::Methylation,
                    ..
                })
            ),
            "{refused:?}"
        );
        assert!(refused.to_string().contains("--methylation"), "{refused}");
        let windows = "chrC\t0\t50000\t30\nchrC\t50000\t100000\t60\n";
        assert!(circle(
            "chrC --circular calls.vcf sample.regions.bed",
            &[("sample.regions.bed", windows)]
        )
        .is_ok());
    }

    /// A circle is not a figure along a sequence, and the builders of those
    /// say so rather than drawing its tracks as bands.
    #[test]
    fn a_circle_is_refused_by_the_builders_of_figures() {
        let mut held = Held::new();
        for (name, text) in circle_files() {
            held.insert(name, text);
        }
        let circular = invocation("chrC --circular genes.gff3");
        let refused = build_figure(&circular, &mut held, |_, _| None, Theme::light(), None);
        assert!(matches!(
            refused,
            Err(BuildError::Uncircled(CircleRefusal::NotAFigure))
        ));
        let refused = build_sheet(&circular, &mut held, |_, _| None, Theme::light());
        assert!(matches!(
            refused,
            Err(BuildError::Uncircled(CircleRefusal::NotAFigure))
        ));
    }

    /// A short feature on a long sequence is still an arc a pointer can find,
    /// and a feature ring is painted a strand a colour, as a band is.
    #[test]
    fn the_annotation_ring_names_its_genes_and_paints_them_by_strand() {
        let svg = circle("chrC --circular genes.gff3", &[]).unwrap();
        let theme = Theme::light();
        let painted = |gene: &str| -> &str {
            let at = svg.find(&format!("<title>{gene}, ")).unwrap();
            let fill = svg[at..].split("fill=\"").nth(1).unwrap();
            fill.split('"').next().unwrap()
        };
        assert_eq!(painted("gA"), theme.color(0));
        assert_eq!(painted("gB"), theme.color(1));
        assert!(
            svg.contains(">gA</text>") && svg.contains(">gB</text>"),
            "{svg}"
        );
        let unnamed = circle("chrC --circular genes.gff3 --no-names", &[]).unwrap();
        assert!(!unnamed.contains(">gA</text>"));
        let one = circle("chrC --circular genes.gff3 --color #123456", &[]).unwrap();
        assert!(one.matches("fill=\"#123456\"").count() >= 2);
    }
}
