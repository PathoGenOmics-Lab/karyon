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
use std::io::{self, Read as _};
use std::path::Path;

use crate::{
    Aggregate, BisulfiteTrack, CladeTrack, CopyNumberTrack, CoverageTrack, DomainTrack,
    DotplotTrack, DynseqTrack, FeatureTrack, Figure, IdeogramTrack, JunctionTrack, LocusTrack,
    LogoTrack, ManhattanTrack, MatrixTrack, MethylationTrack, MsaSequence, MsaTrack, OrfTrack,
    PairStyle, PairTrack, PhylodynamicScale, PhylodynamicTrack, PileupTrack, Plot, Region,
    SelectionEvidence, SelectionTrack, SequenceTrack, SnpTrack, SplitReadTrack, SquiggleTrack,
    StructuralTrack, SurveillanceTrack, SyntenyTrack, TanglegramTrack, Theme, Track, Tree,
    TreeTrack, VariantTrack, WindowStyle, WindowTrack,
};

use crate::cli::args::{
    Invocation, Kind, Palette, Place, Source, Style, Threshold, TrackSpec, TreeSupport,
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
            } => match instead {
                Some(command) => write!(
                    f,
                    "--{track} {path}: {binary}; write <({command}) where its name is, \
                     or turn it into text first"
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
            BuildError::Tree { flag, path, cause } => write!(f, "{flag} {path}: {cause}"),
        }
    }
}

impl std::error::Error for BuildError {}

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
    let mut names = Vec::with_capacity(places.len());
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
        let built = build_one(&one, &mut kept, &mut parsed, theme.clone(), None, true)?;
        gather(&mut legend, &built.legend);
        sheet = sheet.push_bare(&built.figure);
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
    build_one(invocation, files, parsed, theme, window, false)
}

/// The same, where `tolerant` draws a track with nothing in the place as a
/// band that says so rather than refusing the figure, for a panel of a sheet
/// of several places.
fn build_one(
    invocation: &Invocation,
    files: &mut dyn Files,
    mut parsed: impl FnMut(&str, &str) -> Option<Tree>,
    theme: Theme,
    window: Option<&Region>,
    tolerant: bool,
) -> Result<Built, BuildError> {
    let mut kept = KeptStdin { files, stdin: None };
    let files: &mut dyn Files = &mut kept;
    if invocation.genome_wide() && window.is_none() {
        return build_genome(invocation, files, theme);
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
        let context = Context {
            region,
            theme: &theme,
            reference: reference.as_ref(),
            decimals,
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
                    (None, BuildError::Empty { wanted, .. }) if tolerant => Box::new(Nothing {
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
    if let Some(counting) = counting.as_ref().filter(|_| invocation.axis) {
        figure = figure.push_ruler(
            crate::AxisTrack::new()
                .counting()
                .decimals(counting.decimals)
                .label(&counting.unit),
        );
    }
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
    if sequence && bases.items().iter().all(|item| key.items().contains(item)) {
        // And how wide the figure would have to be for the letters, which a
        // reader asked to show a sequence came for.
        let px = figure.px_per_bp();
        let span = region.len() as f64;
        let now = figure.dimensions().0;
        let wanted =
            ((now + (crate::track::sequence::LETTER_PX - px) * span) / 100.0).ceil() * 100.0;
        if wanted <= 100_000.0 {
            files.note(&format!(
                "the bases are blocks of colour at this width, too narrow for their \
                 letters; --width {} draws the letters",
                wanted as u64
            ));
        }
    }
    gather(&mut legend, &key);
    if invocation.legend && !legend.is_empty() {
        figure = figure.push(crate::track::legend::LegendTrack::new(legend.clone()));
    }
    Ok(Built {
        figure,
        along,
        legend,
    })
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
fn thresholded(
    track: ManhattanTrack,
    spec: &TrackSpec,
    p_values: bool,
    path: &str,
) -> Result<ManhattanTrack, BuildError> {
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

/// A scan across a whole genome: every sequence the `--manhattan` tables
/// name, end to end, in the order a reader counts chromosomes, with the
/// bands a genome-wide plot is read by and each sequence named under the
/// scan in place of a ruler of positions no file uses.
///
/// Each sequence is as long as the furthest position a table tests on it,
/// which is how a scan is drawn: an association table says where its
/// markers are and not how long the chromosomes they are on run.
fn build_genome(
    invocation: &Invocation,
    files: &mut dyn Files,
    theme: Theme,
) -> Result<Built, BuildError> {
    let mut scans = Vec::with_capacity(invocation.tracks.len());
    for spec in &invocation.tracks {
        let name = spec.kind.flag();
        let Some(source) = spec.source.as_ref() else {
            continue;
        };
        let (text, path) = fetch(name, source, files)?;
        let read = wrap(name, &path, read::point::genome_associations(&text))?;
        if read.sequences.iter().all(|(_, points)| points.is_empty()) {
            return Err(BuildError::Empty {
                track: name,
                path,
                wanted: "association statistics",
            });
        }
        scans.push((spec, path, read));
    }
    let mut lengths: Vec<(String, u64)> = Vec::new();
    for (_, _, read) in &scans {
        for (sequence, points) in &read.sequences {
            let end = points
                .iter()
                .map(|point| point.pos.saturating_add(1))
                .max()
                .unwrap_or(1);
            match lengths.iter_mut().find(|(named, _)| named == sequence) {
                Some((_, length)) => *length = (*length).max(end),
                None => lengths.push((sequence.clone(), end)),
            }
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
    for (spec, path, read) in scans {
        let points: Vec<crate::Association> = read
            .sequences
            .iter()
            .flat_map(|(sequence, points)| {
                let offset = genome.offset(sequence).unwrap_or(0);
                points
                    .iter()
                    .map(move |point| crate::Association::new(offset + point.pos, point.value))
            })
            .collect();
        let mut track = ManhattanTrack::new(points).bands(genome.boundaries());
        if read.p_values {
            track = track.axis_title("-log10 p");
        }
        let mut track = thresholded(track, spec, read.p_values, &path)?;
        if let Some(height) = spec.height {
            track = track.height(height);
        }
        let label = spec.label.clone().or_else(|| default_label(spec));
        plot = plot.add_track(named(track, label, ManhattanTrack::label));
    }
    let mut figure = plot.add_genome(genome).into_figure();
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

/// Files that name the same sequences, and those sequences, each with how
/// many rows the files give it.
type Naming = (Vec<String>, Vec<(String, usize)>);

/// Where a place named by a word is, and whether it was a gene.
struct Placed {
    region: Region,
    /// The gene's name as its annotation spells it, where the place is a gene.
    gene: Option<String>,
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
    let mut lengths: Vec<(String, u64)> = Vec::new();
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut furthest: Option<(u64, String)> = None;
    let mut annotated = false;
    let mut spelled: Option<String> = None;
    let aliases = called_by(invocation, name);
    let open_error = |spec: &TrackSpec, source: &Source, cause: io::Error| BuildError::Open {
        track: spec.kind.flag(),
        path: called(source),
        cause,
    };
    for spec in &invocation.tracks {
        let sources = [spec.source.as_ref(), spec.second.as_ref()];
        for source in sources.into_iter().flatten() {
            if let Some(held) = files
                .sequences(source)
                .map_err(|cause| open_error(spec, source, cause))?
            {
                lengths.extend(held.into_iter().map(|(n, l)| (renamed(invocation, n), l)));
                continue;
            }
            // A file that will not open is its track's to report.
            let Ok(text) = files.text(source) else {
                continue;
            };
            lengths.extend(
                sequence_lengths(&text)
                    .into_iter()
                    .map(|(n, l)| (renamed(invocation, n), l)),
            );
            // A PAF writes the length of every query it aligns, which is the
            // sequence a synteny figure or a dot plot is drawn along. It was
            // not asked, so a figure placed on its own query was refused.
            if matches!(spec.kind, Kind::Synteny | Kind::Dotplot)
                && spec.source.as_ref() == Some(source)
            {
                lengths.extend(
                    paf_query_lengths(&text)
                        .into_iter()
                        .map(|(n, l)| (renamed(invocation, n), l)),
                );
            }
            if matches!(spec.kind, Kind::Features | Kind::Loci)
                && spec.source.as_ref() == Some(source)
            {
                annotated = true;
                let found = read::interval::named(&text, name);
                spelled = spelled.or(found.spelled);
                spans.extend(
                    found
                        .spans
                        .into_iter()
                        .map(|(sequence, start, end)| (renamed(invocation, sequence), start, end)),
                );
                names.extend(found.names);
            }
            // As far as its rows reach on the sequence, under any name the
            // figure gives it.
            let reach = aliases
                .iter()
                .filter_map(|alias| {
                    if spec.kind == Kind::Manhattan {
                        read::point::association_table(
                            &text,
                            &Region::new(*alias, 0, 1 << 28)
                                .unwrap_or_else(|_| Region::new("x", 0, 1).expect("a window")),
                        )
                        .ok()
                        .and_then(|table| table.points.iter().map(|point| point.pos + 1).max())
                    } else {
                        sequence_column(spec.kind)
                            .and_then(|column| reach_on(&text, &column, alias))
                    }
                })
                .max();
            if let Some(reach) = reach {
                let further = match &furthest {
                    None => true,
                    Some((known, _)) => reach > *known,
                };
                if further {
                    furthest = Some((reach, called(source)));
                }
            }
        }
    }

    if let Some((_, length)) = lengths.iter().find(|(sequence, _)| sequence == name) {
        let region = Region::new(name, 0, (*length).max(1))
            .map_err(|_| nowhere(name, invocation, files, &names, annotated))?;
        return Ok(Placed {
            region,
            gene: None,
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
                .find(|(named, _)| named == sequence)
                .map(|(_, length)| *length);
            let stop = (end + margin).min(length.unwrap_or(u64::MAX));
            let region = Region::new(sequence, start.saturating_sub(margin), stop.max(start + 1))
                .map_err(|_| nowhere(name, invocation, files, &names, annotated))?;
            return Ok(Placed {
                region,
                reached: None,
                gene: Some(spelled.unwrap_or_else(|| name.to_string())),
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

    if let Some((end, file)) = furthest {
        if let Ok(region) = Region::new(name, 0, end.max(1)) {
            return Ok(Placed {
                region,
                gene: None,
                reached: Some(file),
            });
        }
    }
    Err(nowhere(name, invocation, files, &names, annotated))
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
        let mut sequences: Vec<(String, usize)> = match files.sequences(source) {
            Ok(Some(lengths)) => lengths.into_iter().map(|(name, _)| (name, 0)).collect(),
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
    // Names within two edits, or the same letters in another case.
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
        near: near
            .into_iter()
            .take(3)
            .map(|(_, candidate)| candidate.clone())
            .collect(),
        annotated,
    }
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
fn strip(
    spec: &TrackSpec,
    sheet: Option<&(read::sheet::Sheet, String)>,
    rows: &[String],
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
    Ok(Some(Traits::from_sheet(held).strips(wanted)))
}

/// A track's own file, what it was called, and whether it arrived converted
/// from a binary alignment file rather than read as it is.
///
/// A track that draws reads, or the depth they add up to, asks the files for
/// that first: a BAM on disk answers with the reads over the window, or their
/// depth as bedGraph, and anything that cannot answer is read as text.
fn slurp(
    spec: &TrackSpec,
    region: &Region,
    files: &mut dyn Files,
) -> Result<(String, String, bool), BuildError> {
    let Some(source) = spec.source.as_ref() else {
        return Ok((String::new(), String::new(), false));
    };
    let converted = match spec.kind {
        Kind::Coverage => files.depth(source, region),
        Kind::Pileup | Kind::SplitReads => files.reads(source, region),
        _ => Ok(None),
    }
    .map_err(|cause| BuildError::Open {
        track: spec.kind.flag(),
        path: called(source),
        cause,
    })?;
    match converted {
        Some(text) => Ok((text, called(source), true)),
        None => {
            let (text, name) = fetch(spec.kind.flag(), source, files)?;
            Ok((text, name, false))
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
fn shortened(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
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
        inside
    } else {
        bytes
    };
    String::from_utf8(bytes).map_err(|error| {
        // What a genomics file is when it is not text is nearly always one of
        // a few formats, and naming it is what turns "stream did not contain
        // valid UTF-8" into the command that reads it.
        match Binary::of(error.as_bytes(), path) {
            Some(binary) => io::Error::new(io::ErrorKind::InvalidData, binary),
            None => io::Error::new(io::ErrorKind::InvalidData, error),
        }
    })
}

/// Where a figure's files come from.
///
/// Text is the one question every source answers: a shell reads a path, a
/// page looks a name up among its buffers, a test hands over a literal, and
/// any closure from a source to its text is a `Files`. The other questions
/// are for sources that are not text, and a reader that cannot answer them
/// says `None`, which sends the track to [`Files::text`] instead. [`Disk`]
/// answers them for a BAM, reading only the reads over the window through
/// its index.
pub trait Files {
    /// The text a source holds.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn text(&mut self, source: &Source) -> io::Result<String>;

    /// The depth of the reads a binary alignment file holds over `region`, as
    /// bedGraph, or `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn depth(&mut self, _source: &Source, _region: &Region) -> io::Result<Option<String>> {
        Ok(None)
    }

    /// The reads a binary alignment file holds over `region`, as SAM text, or
    /// `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn reads(&mut self, _source: &Source, _region: &Region) -> io::Result<Option<String>> {
        Ok(None)
    }

    /// The sequences a binary file names, each with its length, or `None`
    /// for a source that is not one this can read.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn sequences(&mut self, _source: &Source) -> io::Result<Option<Vec<(String, u64)>>> {
        Ok(None)
    }

    /// The records of one read in a binary alignment file, by its name, as
    /// SAM text, or `None` for a source this cannot read that way.
    ///
    /// # Errors
    ///
    /// Whatever stopped it being read.
    fn named_read(&mut self, _source: &Source, _name: &str) -> io::Result<Option<String>> {
        Ok(None)
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

/// The files a command line names, read from disk.
///
/// Text is read whole, and compressed text is taken out of its wrapper, by
/// [`open_from_disk`]. A BAM is read a window at a time, through the `.bai`
/// beside it where there is one, so a figure of one gene reads the blocks that
/// gene is in.
///
/// A pipe the shell named, as `<(zcat genes.gff3.gz)`, cannot be read twice,
/// so it is kept once read, since placing a figure by a gene's name reads the
/// annotation before the track does. Standard input is read once and kept by
/// [`build_files`], which every figure is drawn through. A file on disk is
/// read again instead: kept, every file was held twice while its track was
/// built, once here and once in the text handed out, and a depth file of
/// 176 MB took 435 MB to draw where it now takes 259.
#[derive(Debug, Default)]
pub struct Disk {
    kept: std::collections::HashMap<Source, String>,
    /// What [`Files::note`] was told, for the command line to print.
    pub notes: Vec<String>,
}

impl Disk {
    /// The header and the reads over `region`, for a source that is a BAM.
    fn bam(
        &mut self,
        source: &Source,
        region: &Region,
    ) -> io::Result<Option<(read::bam::Header, Vec<read::bam::Record>)>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        if !is_bam(path)? {
            return Ok(None);
        }
        let index = bam_index(path)?;
        let file = io::BufReader::new(fs::File::open(path)?);
        read::bam::window(file, index.as_ref(), region)
            .map(Some)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
    }
}

impl Files for Disk {
    fn note(&mut self, message: &str) {
        self.notes.push(message.to_string());
    }

    fn text(&mut self, source: &Source) -> io::Result<String> {
        if let Some(text) = self.kept.get(source) {
            return Ok(text.clone());
        }
        let text = open_from_disk(source)?;
        if let Source::Path(path) = source {
            if !fs::metadata(path).is_ok_and(|meta| meta.is_file()) {
                self.kept.insert(source.clone(), text.clone());
            }
        }
        Ok(text)
    }

    fn depth(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(self.bam(source, region)?.map(|(_, records)| {
            read::bam::bedgraph(region.seq(), region, &read::bam::depth(&records, region))
        }))
    }

    fn reads(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(self
            .bam(source, region)?
            .map(|(header, records)| read::bam::sam(&header, &records)))
    }

    fn named_read(&mut self, source: &Source, name: &str) -> io::Result<Option<String>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        if !is_bam(path)? {
            return Ok(None);
        }
        let file = io::BufReader::new(fs::File::open(path)?);
        read::bam::named(file, name)
            .map(|(header, records)| Some(read::bam::sam(&header, &records)))
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
    }

    fn sequences(&mut self, source: &Source) -> io::Result<Option<Vec<(String, u64)>>> {
        let Source::Path(path) = source else {
            return Ok(None);
        };
        if !is_bam(path)? {
            return Ok(None);
        }
        let file = io::BufReader::new(fs::File::open(path)?);
        read::bam::header_of(file)
            .map(|header| Some(header.references))
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
    }
}

/// Whether a file is a BAM: bgzip on the outside, and BAM's magic inside.
fn is_bam(path: &Path) -> io::Result<bool> {
    let mut first = [0u8; 3];
    let mut file = fs::File::open(path)?;
    if file.read(&mut first)? < first.len() || !read::gzip::is_gzip(&first) {
        return Ok(false);
    }
    Ok(read::bam::header_of(io::BufReader::new(fs::File::open(path)?)).is_ok())
}

/// The `.bai` beside a BAM, as `reads.bam.bai` or `reads.bai`, if there is one.
fn bam_index(path: &Path) -> io::Result<Option<read::bam::Index>> {
    let mut beside = path.as_os_str().to_owned();
    beside.push(".bai");
    for candidate in [std::path::PathBuf::from(beside), path.with_extension("bai")] {
        if candidate.is_file() {
            let bytes = fs::read(&candidate)?;
            return read::bam::index(&bytes)
                .map(Some)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()));
        }
    }
    Ok(None)
}

/// The files a command line names, held in memory by name.
///
/// What [`Disk`] reads from a path, read from bytes a caller already has: a
/// page that fetched them, a service that was sent them, a test. Compressed
/// text is taken out of its wrapper, a BAM is read a window at a time through
/// the `.bai` held beside it, and one read is found by its name, so a figure
/// drawn from these files is the figure a shell draws from the same files on
/// disk. A file is found by its name as the command line writes it.
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
    files: std::collections::BTreeMap<String, Vec<u8>>,
    /// What [`Files::note`] was told.
    pub notes: Vec<String>,
}

impl Held {
    /// Holds no files yet.
    pub fn new() -> Self {
        Held::default()
    }

    /// Holds `bytes` as the file called `name`, in place of any held before.
    pub fn insert(&mut self, name: impl Into<String>, bytes: impl Into<Vec<u8>>) {
        self.files.insert(name.into(), bytes.into());
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
            Some(bytes) => Ok((path.as_path(), bytes.as_slice())),
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

    /// A source's bytes where it is a BAM: bgzip on the outside, and BAM's
    /// header inside. `None` for any other file, and for one not held, which
    /// [`Files::text`] then says.
    fn bam<'s>(&self, source: &'s Source) -> Option<(&'s Path, &[u8])> {
        let (path, bytes) = self.held(source).ok()?;
        let bam =
            read::gzip::is_gzip(bytes) && read::bam::header_of(io::Cursor::new(bytes)).is_ok();
        bam.then_some((path, bytes))
    }

    /// The `.bai` held beside a BAM, as `reads.bam.bai` or `reads.bai`, if
    /// there is one.
    fn bam_index(&self, path: &Path) -> io::Result<Option<read::bam::Index>> {
        let mut beside = path.as_os_str().to_owned();
        beside.push(".bai");
        for candidate in [std::path::PathBuf::from(beside), path.with_extension("bai")] {
            if let Some(bytes) = self.files.get(&candidate.display().to_string()) {
                return read::bam::index(bytes).map(Some).map_err(|error| {
                    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
                });
            }
        }
        Ok(None)
    }

    /// The header and the reads over `region`, for a source that is a BAM.
    fn window(
        &self,
        source: &Source,
        region: &Region,
    ) -> io::Result<Option<(read::bam::Header, Vec<read::bam::Record>)>> {
        let Some((path, bytes)) = self.bam(source) else {
            return Ok(None);
        };
        let index = self.bam_index(path)?;
        read::bam::window(io::Cursor::new(bytes), index.as_ref(), region)
            .map(Some)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
    }
}

impl Files for Held {
    fn note(&mut self, message: &str) {
        self.notes.push(message.to_string());
    }

    fn text(&mut self, source: &Source) -> io::Result<String> {
        let (path, bytes) = self.held(source)?;
        decoded(bytes.to_vec(), Some(path))
    }

    fn depth(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(self.window(source, region)?.map(|(_, records)| {
            read::bam::bedgraph(region.seq(), region, &read::bam::depth(&records, region))
        }))
    }

    fn reads(&mut self, source: &Source, region: &Region) -> io::Result<Option<String>> {
        Ok(self
            .window(source, region)?
            .map(|(header, records)| read::bam::sam(&header, &records)))
    }

    fn named_read(&mut self, source: &Source, name: &str) -> io::Result<Option<String>> {
        let Some((_, bytes)) = self.bam(source) else {
            return Ok(None);
        };
        read::bam::named(io::Cursor::new(bytes), name)
            .map(|(header, records)| Some(read::bam::sam(&header, &records)))
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
    }

    fn sequences(&mut self, source: &Source) -> io::Result<Option<Vec<(String, u64)>>> {
        let Some((_, bytes)) = self.bam(source) else {
            return Ok(None);
        };
        read::bam::header_of(io::Cursor::new(bytes))
            .map(|header| Some(header.references))
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
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
    /// BCF, which is bgzip by its bytes and BCF by its name.
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
            Binary::Hic => {
                "hic2cool convert FILE.hic FILE.cool -r N writes one resolution as a .cool, \
                 and cooler dump --join -r REGION FILE.cool writes that as the BEDPE --pairs \
                 reads"
            }
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
                    // A pipe has no name to write a command in place of.
                    let instead = matches!(spec.source, Some(Source::Path(_)))
                        .then(|| binary.reader(spec.kind, &path, region))
                        .flatten();
                    BuildError::NotText {
                        track,
                        path,
                        binary,
                        instead,
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
            // drew has no need to know.
            let held = match (region, spec.source.as_ref()) {
                (Some(_), Some(source @ Source::Path(_))) => files
                    .text(source)
                    .map(|text| rows_on(spec.kind, &text))
                    .unwrap_or_default(),
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
fn track(
    spec: &TrackSpec,
    context: &Context<'_>,
    files: &mut dyn Files,
    parsed: &mut dyn FnMut(&str, &str) -> Option<Tree>,
    legend: &mut crate::track::legend::Legend,
) -> Result<Box<dyn Track>, BuildError> {
    let region = context.region;
    let theme = context.theme;
    let (text, path, converted) = slurp(spec, region, files)?;
    // A file named on its own was placed by its name, and a few names hide
    // another format: modkit writes its bedMethyl as `.bed`.
    let refined;
    let spec = match spec.guessed.then(|| refine(spec.kind, &text)).flatten() {
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

    let empty = |wanted: &'static str| BuildError::Empty {
        track: name,
        path: path.clone(),
        wanted,
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
                // Depth worked out from a BAM arrives as bedGraph, whatever
                // `--format` said about a text file.
                read::signal::fold_spans(
                    &text,
                    region,
                    if converted {
                        Some(crate::Format::BedGraph)
                    } else {
                        spec.format
                    },
                    |start, end, value| painted.paint(start, end, value),
                ),
            )?;
            drop(text);
            // Named on its own, a BAM is its depth; over a window a few reads
            // wide the reads are what a reader came for, and nothing said
            // they could be drawn.
            if converted && spec.guessed && region.len() <= READS_WINDOW {
                files.note(&format!(
                    "{path} is drawn as its depth; --pileup {path} draws its reads"
                ));
            }
            if spans == 0 {
                return Err(empty("values"));
            }
            let mut track = painted.aggregate(spec.aggregate.unwrap_or(Aggregate::Max));
            if let Some(style) = spec.style.and_then(|style| style.coverage()) {
                track = track.style(style);
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
            Box::new(named(track, label, CoverageTrack::label))
        }
        Kind::Dynseq => {
            let Some(source) = spec.second.as_ref() else {
                return Err(BuildError::MissingSecond { track: name });
            };
            let reference = second_sequence(name, source, region, files)?;

            let found = wrap(name, &path, read::dynseq::scores(&text, region))?;
            if found.records == 0 {
                return Err(empty("scores"));
            }
            if found.spans.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "scores",
                    held: found.records,
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
            if found.records == 0 {
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
                    held: found.records,
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
            let (from, bases) = sequence(name, &path, &text, region)?.clip(region)?;
            let mut track = SequenceTrack::new(from, bases);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, SequenceTrack::label))
        }
        Kind::Features => {
            let features = wrap(
                name,
                &path,
                read::interval::features(&text, region, spec.format),
            )?;
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
            let variants = wrap(name, &path, read::point::variants(&text, region))?;
            if variants.is_empty() {
                return Err(empty("variants"));
            }
            // A stem is as tall as the allele fraction the VCF's AF gives, and
            // the axis says so; a file with no AF draws no axis to title.
            let mut track = VariantTrack::new(variants).axis_title("AF");
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
            Box::new(named(track, label, VariantTrack::label))
        }
        Kind::Windows => {
            let windows = wrap(name, &path, read::signal::windows(&text, region))?;
            if windows.is_empty() {
                return Err(empty("windows"));
            }
            let mut track = WindowTrack::new(windows).style(
                spec.style
                    .and_then(|s| s.window())
                    .unwrap_or(WindowStyle::Steps),
            );
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, WindowTrack::label))
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
            let held = strip(spec, sheet.as_ref(), &leaves)?;

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
                track = track.traits(held);
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
            if found.records == 0 {
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
                    held: found.records,
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
            let found = wrap(name, &path, read::structural::variants(&text, region))?;
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            let (from, bases) = sequence(name, &path, &text, region)?.clip(region)?;
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            let held = wrap(name, &path, read::segments::samples(&text))?;
            if spec.sample.is_some() && held.is_empty() {
                // A flag accepted and then ignored gives a figure that is not
                // the one asked for and does not look wrong: this table names
                // no samples, so the whole of it would be drawn under a name
                // the command asked to pick out of it.
                return Err(BuildError::Ambiguous {
                    track: name,
                    path: path.clone(),
                    flag: "--sample",
                    choices: vec!["no sample column".to_string()],
                });
            }
            if held.len() > 1 && spec.sample.is_none() {
                return Err(BuildError::Ambiguous {
                    track: name,
                    path: path.clone(),
                    flag: "--sample",
                    choices: held,
                });
            }
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
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_row_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            if let Some(color) = &spec.color {
                track = track.color(color);
            }
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, CoverageTrack::label))
        }
        Kind::Pairs => {
            let (pairs, measured) = wrap(name, &path, read::pairs::pairs(&text, region))?;
            if pairs.is_empty() {
                return Err(empty("pairs"));
            }
            // A triangle where most of the pairs the places could make were
            // measured, as linkage and contacts are, and arcs where a few were.
            // Linkage is a triangle however PLINK's window filtered it.
            let correlation = measured.as_deref().is_some_and(read::pairs::is_correlation);
            let style = spec.style.and_then(Style::pairs).unwrap_or_else(|| {
                if correlation {
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
            if let Some(px) = spec.row_height {
                track = track.row_height(px);
            }
            if spec.no_names {
                track = track.show_row_names(false);
            }
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
        let device = Source::Path("/dev/null".into());
        assert_eq!(disk.text(&device).unwrap(), "");
        assert!(disk.kept.contains_key(&device), "a device was not kept");
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
    #[test]
    fn a_file_that_is_not_text_is_answered_with_the_command_that_reads_it() {
        let line = |text: &str| -> Invocation {
            let args: Vec<String> = text.split_whitespace().map(String::from).collect();
            match parse(&args).unwrap() {
                Request::Draw(invocation) => *invocation,
                other => panic!("expected a figure, got {other:?}"),
            }
        };
        for (command, binary, wanted) in [
            (
                "chr1:1-5000 --variants calls.vcf.gz",
                Binary::Gzip,
                "<(gzip -dc calls.vcf.gz)",
            ),
            (
                "chr1:1-5000 --pileup reads.bam",
                Binary::Bam,
                "<(samtools view -h reads.bam chr1:1-5000)",
            ),
            (
                "chr1:1-5000 --coverage reads.bam",
                Binary::Bam,
                "<(samtools depth -a -r chr1:1-5000 reads.bam)",
            ),
            (
                "chr1:1-5000 --variants calls.bcf",
                Binary::Bcf,
                "<(bcftools view calls.bcf)",
            ),
            (
                "chr1:1-5000 --coverage depth.bw",
                Binary::BigWig,
                "<(bigWigToBedGraph -chrom=chr1 -start=0 -end=5000 depth.bw /dev/stdout)",
            ),
            (
                "chr1:1-5000 contacts.cool",
                Binary::Cool,
                "<(cooler dump --join -r chr1:1-5000 contacts.cool)",
            ),
            // No one command: the resolution has to be picked first.
            (
                "chr1:1-5000 contacts.mcool",
                Binary::Mcool,
                "cooler ls lists its resolutions",
            ),
            (
                "chr1:1-5000 --pairs contacts.hic",
                Binary::Hic,
                "hic2cool convert FILE.hic FILE.cool -r N",
            ),
        ] {
            let error = build(&line(command), |_| {
                Err(io::Error::new(io::ErrorKind::InvalidData, binary))
            })
            .unwrap_err()
            .to_string();
            assert!(error.contains(wanted), "{command}: {error}");
            assert!(error.contains("karyon reads text"), "{error}");
        }
        // A pipe has no name to write a command in place of.
        let error = build(&line("chr1:1-5000 --pileup -"), |_| {
            Err(io::Error::new(io::ErrorKind::InvalidData, Binary::Bam))
        })
        .unwrap_err()
        .to_string();
        assert!(error.contains("pipe it through the tool"), "{error}");
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
        assert_eq!(
            from_memory(&line, &mut held),
            drawn_from_disk(&line).unwrap()
        );
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
}
