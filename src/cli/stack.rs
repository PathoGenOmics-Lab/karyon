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

use crate::{
    Aggregate, BisulfiteTrack, CladeTrack, CopyNumberTrack, CoverageTrack, DomainTrack,
    DotplotTrack, FeatureTrack, IdeogramTrack, LocusTrack, LogoTrack, ManhattanTrack, MatrixTrack,
    MethylationTrack, MsaSequence, MsaTrack, OrfTrack, PileupTrack, Plot, Region, SequenceTrack,
    SnpTrack, SplitReadTrack, StructuralTrack, SyntenyTrack, TanglegramTrack, Theme, Track, Tree,
    TreeTrack, VariantTrack, WindowStyle, WindowTrack,
};

use crate::cli::args::{Invocation, Kind, Palette, Source, TrackSpec};
use crate::read;
use crate::track::traits::Traits;

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
    /// A column was asked for by name and the sheet has no such column.
    ///
    /// Not [`BuildError::Ambiguous`], which is a file holding several things
    /// and a command naming none of them. Here the command named one and the
    /// file has not got it, which is nearly always a spelling, so the columns
    /// it does have are worth printing beside it.
    Unnamed {
        /// Which track wanted it.
        track: &'static str,
        /// What the sheet was called.
        path: String,
        /// The column that is not in it.
        wanted: String,
        /// The columns that are.
        held: Vec<String>,
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
            } => {
                write!(f, "--{track} {path}: no {wanted} in {region}")?;
                write!(f, ", though the file holds {held}")?;
                if !named.is_empty() {
                    write!(f, " on {named}")?;
                }
                Ok(())
            }
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
                wanted,
                held,
            } => write!(
                f,
                "--{track} {path} has no column called {wanted}; it has {}",
                held.join(", ")
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
    mut open: impl FnMut(&Source) -> io::Result<String>,
) -> Result<String, BuildError> {
    let region = &invocation.region;
    let mut plot = Plot::over(region.clone());
    if let Some(title) = &invocation.title {
        plot = plot.title(title);
    }
    if let Some(width) = invocation.width {
        plot = plot.width(width);
    }
    if invocation.theme == Palette::Dark {
        plot = plot.theme(Theme::dark());
    }
    if !invocation.region_label {
        plot = plot.remove_region_label();
    }
    if !invocation.axis {
        plot = plot.remove_axis();
    }

    for spec in &invocation.tracks {
        // The ruler is the one track that reads nothing, and the one the plot
        // has to be told about so it does not append a second.
        if spec.kind == Kind::Axis {
            let mut axis = plot.add_axis();
            if let Some(label) = &spec.label {
                axis = axis.label(label);
            }
            if let Some(height) = spec.height {
                axis = axis.adjust(|track| track.height(height));
            }
            plot = axis.done();
            continue;
        }
        plot = plot.add_boxed(track(spec, region, &mut open)?);
    }
    Ok(plot.to_svg())
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
    open: &mut dyn FnMut(&Source) -> io::Result<String>,
) -> Result<Option<(read::sheet::Sheet, String)>, BuildError> {
    let Some(source) = spec.traits.as_ref() else {
        return Ok(None);
    };
    let name = spec.kind.flag();
    let (text, path) = fetch(name, source, open)?;
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

    Ok(Some(Traits::new(held.rows.clone()).spread(wanted)))
}

fn slurp(
    spec: &TrackSpec,
    open: &mut dyn FnMut(&Source) -> io::Result<String>,
) -> Result<(String, String), BuildError> {
    let Some(source) = spec.source.as_ref() else {
        return Ok((String::new(), String::new()));
    };
    fetch(spec.kind.flag(), source, open)
}

/// Reads one source, and says what it was called.
///
/// Split out from [`slurp`] because a track drawn from two files opens the
/// second the same way it opened the first, and both belong to the same track
/// as far as an error message is concerned.
fn fetch(
    track: &'static str,
    source: &Source,
    open: &mut dyn FnMut(&Source) -> io::Result<String>,
) -> Result<(String, String), BuildError> {
    let name = match source {
        Source::Path(path) => path.display().to_string(),
        Source::Stdin => "standard input".to_string(),
    };
    let text = open(source).map_err(|cause| BuildError::Open {
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
    match source {
        Source::Path(path) => fs::read_to_string(path),
        Source::Stdin => {
            let mut text = String::new();
            io::stdin().read_to_string(&mut text)?;
            Ok(text)
        }
    }
}

/// Builds one track from its flags and its file.
fn track(
    spec: &TrackSpec,
    region: &Region,
    open: &mut dyn FnMut(&Source) -> io::Result<String>,
) -> Result<Box<dyn Track>, BuildError> {
    let (text, path) = slurp(spec, open)?;
    let name = spec.kind.flag();
    let label = spec.label.clone();
    let height = spec.height;

    let empty = |wanted: &'static str| BuildError::Empty {
        track: name,
        path: path.clone(),
        wanted,
    };
    // Read before the match rather than inside the arms, because two arms
    // shadow `text` with a second file of their own and a sheet fetched after
    // that would be read out of the wrong one.
    let sheet = sheet(spec, open)?;

    let built: Box<dyn Track> = match spec.kind {
        Kind::Coverage => {
            let pairs = wrap(name, &path, read::signal::dense(&text, region, spec.format))?;
            if pairs.is_empty() {
                return Err(empty("values"));
            }
            let mut track = CoverageTrack::from_pairs(region, pairs)
                .aggregate(spec.aggregate.unwrap_or(Aggregate::Max));
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
        Kind::Sequence => {
            let records = wrap(name, &path, read::seq::fasta(&text))?;
            let (_, bases) = records
                .into_iter()
                .next()
                .ok_or_else(|| empty("sequence"))?;
            let mut track = SequenceTrack::new(region.start(), clip(&bases, region));
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
            let mut track = VariantTrack::new(variants);
            if let Some(height) = height {
                track = track.height(height);
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
            let points = wrap(name, &path, read::point::associations(&text, region))?;
            if points.is_empty() {
                return Err(empty("association statistics"));
            }
            let mut track = ManhattanTrack::new(points);
            if let Some(height) = height {
                track = track.height(height);
            }
            Box::new(named(track, label, ManhattanTrack::label))
        }
        Kind::Tree => {
            let tree = Tree::parse_newick(text.trim()).map_err(|cause| BuildError::Tree {
                flag: "--tree",
                path: path.clone(),
                cause,
            })?;
            Box::new(named(TreeTrack::new(tree), label, TreeTrack::label))
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
            let (other, right_path) = fetch(name, source, open)?;
            let parse = |text: &str, flag, path: &str| {
                Tree::parse_newick(text.trim()).map_err(|cause| BuildError::Tree {
                    flag,
                    path: path.to_string(),
                    cause,
                })
            };
            let left = parse(&text, "--tanglegram", &path)?;
            let right = parse(&other, "--against", &right_path)?;
            // Named, because two phylogenies side by side with nothing over
            // them do not say which is which, and which is which is the whole
            // of what a tanglegram is read for.
            let track =
                TanglegramTrack::new(left, right).names(shortened(&path), shortened(&right_path));
            Box::new(named(track, label, TanglegramTrack::label))
        }
        // One bedMethyl is one track only when it counted one modification. A
        // dual-mode run writes m and h at the same cytosine, and stacked on one
        // axis those are two marks at one position with nothing naming either.
        Kind::Methylation => {
            let held = wrap(name, &path, read::methyl::codes(&text))?;
            if held.is_empty() {
                return Err(empty("modified bases"));
            }
            let code = chosen(name, &path, "--modification", &held, &spec.selects)?;

            let found = wrap(name, &path, read::methyl::sites(&text, region, &code))?;
            if found.records == 0 {
                return Err(empty("modified bases"));
            }
            if found.sites.is_empty() {
                return Err(BuildError::Elsewhere {
                    track: name,
                    path: path.clone(),
                    wanted: "modified bases",
                    held: found.records,
                    named: code.clone(),
                    region: region.to_string(),
                });
            }

            let mut track = MethylationTrack::new(found.sites);
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
                });
            }
            let mut track = StructuralTrack::new(found.variants);
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
                });
            }
            let track = SplitReadTrack::new(found.reads);
            Box::new(named(track, label, SplitReadTrack::label))
        }
        // One row per molecule and one column per site. The reader builds the
        // grid by position, since a call written into the wrong column is a
        // methylation pattern that never existed, drawn as cleanly as one that
        // did, and nothing downstream could tell.
        Kind::Bisulfite => {
            let held = wrap(name, &path, read::bisulfite::contexts(&text))?;
            if held.is_empty() {
                return Err(empty("methylation calls"));
            }
            let context = chosen(name, &path, "--context", &held, &spec.selects)?;

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
                });
            }

            let track = BisulfiteTrack::new(found.sites, found.molecules);
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
                track = track.traits(traits);
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
            let (newick, tree_path) = fetch(name, source, open)?;
            let tree = Tree::parse_newick(newick.trim()).map_err(|cause| BuildError::Tree {
                flag: "--with-tree",
                path: tree_path,
                cause,
            })?;

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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
                });
            }

            let (text, link_path) = fetch(name, source, open)?;
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
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
            let records = wrap(name, &path, read::seq::fasta(&text))?;
            let (_, bases) = records.into_iter().next().ok_or_else(|| empty("orfs"))?;
            let track = OrfTrack::new(region.start(), clip(&bases, region));
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
                track = track.traits(traits);
            }
            Box::new(named(track, label, MsaTrack::label))
        }
        Kind::Snps => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            let names: Vec<String> = sequences.iter().map(|row| row.name.clone()).collect();
            let mut track = SnpTrack::from_alignment(0, &sequences);
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
                track = track.traits(traits);
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
            let ploidy = spec.ploidy.unwrap_or(2.0);
            let held = wrap(name, &path, read::segments::samples(&text))?;
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
                });
            }

            let track = CopyNumberTrack::at_ploidy(found.segments, ploidy);
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
            if let Some(traits) = strip(spec, sheet.as_ref(), &names)? {
                track = track.traits(traits);
            }
            Box::new(named(track, label, MatrixTrack::label))
        }
        Kind::Pileup => {
            let reads = wrap(name, &path, read::align::sam(&text, region))?;
            if reads.is_empty() {
                return Err(empty("reads"));
            }
            Box::new(named(PileupTrack::new(reads), label, PileupTrack::label))
        }
        Kind::Axis => unreachable!("the ruler is added by build"),
    };
    Ok(built)
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

/// Cuts a whole reference down to the region on display.
///
/// A FASTA record is the whole sequence and the figure is a window on it, so
/// the bases the region names are the ones the track gets.
fn clip(bases: &[u8], region: &Region) -> Vec<u8> {
    let start = region.start() as usize;
    if start >= bases.len() {
        return Vec::new();
    }
    let end = (region.end() as usize).min(bases.len());
    bases[start..end].to_vec()
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
        let from_library = Figure::new(region.clone())
            .push(OrfTrack::new(100, bases[100..400].to_vec()))
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
        // A column of an alignment is its own coordinate, so it anchors at nought.
        let right = Figure::new(Region::parse("aln:1-10").unwrap())
            .push(LogoTrack::from_sequences(0, &rows))
            .push(crate::AxisTrack::new())
            .to_svg();
        let wrong = Figure::new(Region::parse("aln:1-10").unwrap())
            .push(LogoTrack::from_sequences(1, &rows))
            .push(crate::AxisTrack::new())
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
    fn over(locus: &str, flag: &str, path: &str) -> Invocation {
        let args = vec![locus.to_string(), flag.to_string(), path.to_string()];
        match parse(&args).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    #[test]
    fn a_reference_is_cut_down_to_the_region() {
        let region = Region::parse("chr1:11-20").unwrap();
        let bases: Vec<u8> = (b'A'..=b'Z').collect();
        assert_eq!(clip(&bases, &region), b"KLMNOPQRST".to_vec());
    }

    #[test]
    fn a_reference_shorter_than_the_region_is_not_an_index_panic() {
        let region = Region::parse("chr1:1-1000").unwrap();
        assert_eq!(clip(b"ACGT", &region), b"ACGT".to_vec());
        let past = Region::parse("chr1:100-200").unwrap();
        assert!(clip(b"ACGT", &past).is_empty());
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
}
