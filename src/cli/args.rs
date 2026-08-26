//! The command line, which is the stack.
//!
//! A figure is a list of tracks in the order they are drawn, and `argv` is a
//! list in the order it was written, so the grammar is the obvious one: each
//! `--<track>` flag starts a track, and the flags after it describe that track
//! until the next one starts. It is [`crate::Plot`] with spaces instead of
//! dots, and the correspondence is exact:
//!
//! ```text
//! --coverage depth.bg --label depth --aggregate min
//!
//!     .add_coverage(..).label("depth").adjust(|t| t.aggregate(Aggregate::Min))
//! ```
//!
//! Figure flags such as `--title` are not attached to anything and may sit
//! anywhere.
//!
//! # What position alone cannot say
//!
//! Where a word sits is the whole of what binds it to a track, and two things
//! fall outside that.
//!
//! One is a track whose data is not one file. A `--<track>` flag takes one
//! path, and a tanglegram is two trees, so the second arrives by a name of its
//! own: `--tanglegram left.nwk --against right.nwk`. It is spelled by what the
//! file means rather than by where it sits, because a second phylogeny and a
//! second table are not interchangeable, and a track that takes one is refused
//! without it. That refusal is the point: a tanglegram of one tree against
//! itself has no crossings at all, which is what a perfect answer looks like.
//! Eighteen of the crate's thirty-three track types are what the command line
//! reaches.
//!
//! The other is a modifier the track before it has no use for, which the order
//! of the words does nothing to prevent. Every modifier therefore carries the
//! tracks it says anything to and is refused by name anywhere else: a flag
//! written and then ignored gives a figure that is not the one asked for and
//! does not look wrong.
//!
//! # Where `--format` overrules a guess
//!
//! A signal file is bedGraph, `samtools depth` or a bare column of values
//! according to how many columns it has, and an interval file is BED or GFF3
//! according to a pragma and to column seven. Either guess can be wrong without
//! failing, and a wrong one moves every coordinate in the figure by a base,
//! which is what `--format` settles. It is the one modifier not checked here
//! against the track it follows, its words spelling formats for two readers at
//! once: coverage and features ask for a format, and anywhere else it parses
//! and goes nowhere.

use std::fmt;
use std::path::PathBuf;

use crate::read::locus::Identity;
use crate::{Aggregate, CoverageStyle, Format, Region, WindowStyle};

/// What went wrong before anything was read.
#[derive(Debug)]
pub enum ArgError {
    /// No arguments at all.
    NoArguments,
    /// A flag that is not in the grammar.
    UnknownFlag(String),
    /// A flag whose value is missing.
    MissingValue(&'static str),
    /// A value that is not one of the words the flag takes.
    BadValue {
        /// The flag.
        flag: &'static str,
        /// What was written.
        given: String,
        /// What it takes.
        expected: &'static str,
    },
    /// A track modifier written before any track.
    NoTrackYet(&'static str),
    /// A modifier the track it follows has no use for.
    WrongTrack {
        /// The modifier.
        flag: &'static str,
        /// The track it landed on.
        track: &'static str,
    },
    /// The first argument was not a locus string.
    BadRegion(crate::Error),
    /// A locus was given twice, or a positional argument came after one.
    ExtraRegion(String),
    /// No locus at all.
    NoRegion,
    /// A locus whose span is larger than a figure is drawn over.
    HugeRegion {
        /// The locus as it was written.
        given: String,
        /// How many bases it spans.
        span: u64,
    },
    /// More than one track wanted standard input.
    StdinTwice,
    /// A track drawn from two files was given one.
    MissingSecond {
        /// The flag that names the other file.
        flag: &'static str,
        /// The track that is short of it.
        track: &'static str,
    },
}

/// The widest figure that is drawn, in pixels.
///
/// A width becomes a column of pixels per band, and no screen or page goes
/// past this, so a larger number is a typo rather than an intention.
const MAX_WIDTH: f64 = 100_000.0;

/// The longest region that is drawn, in bases.
///
/// A per-base track keeps one value for every base of the window, so the span
/// is what a figure costs in memory. This is well above the longest sequence
/// anyone assembles and below the point where the buffer cannot be allocated.
const MAX_SPAN: u64 = 1 << 32;

impl fmt::Display for ArgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ArgError::NoArguments => write!(f, "no arguments: try karyon --help"),
            ArgError::UnknownFlag(flag) => write!(f, "unknown flag {flag}"),
            ArgError::MissingValue(flag) => write!(f, "{flag} needs a value"),
            ArgError::BadValue {
                flag,
                given,
                expected,
            } => write!(f, "{flag} does not take {given:?}, only {expected}"),
            ArgError::NoTrackYet(flag) => write!(
                f,
                "{flag} describes the track before it, and no track has been given yet"
            ),
            ArgError::WrongTrack { flag, track } => {
                write!(f, "{flag} means nothing to a {track} track")
            }
            ArgError::BadRegion(error) => write!(f, "{error}"),
            ArgError::ExtraRegion(extra) => {
                write!(f, "one region per figure, and {extra:?} is a second one")
            }
            ArgError::NoRegion => write!(
                f,
                "the first argument is the region, as in NC_000962.3:761,000-763,000"
            ),
            ArgError::HugeRegion { given, span } => write!(
                f,
                "{given:?} spans {span} bases, and a figure is drawn over at most {MAX_SPAN}"
            ),
            ArgError::StdinTwice => write!(f, "only one track can read from standard input"),
            ArgError::MissingSecond { flag, track } => write!(
                f,
                "a {track} track is drawn from two files, and {flag} names the second"
            ),
        }
    }
}

impl std::error::Error for ArgError {}

/// Where a track's data comes from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Source {
    /// A file on disk.
    Path(PathBuf),
    /// Standard input, written as `-`.
    Stdin,
}

/// Which track a `--<track>` flag asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Per-base signal from bedGraph or `samtools depth`.
    Coverage,
    /// Reference bases from FASTA.
    Sequence,
    /// Intervals from BED or GFF3.
    Features,
    /// Point calls from VCF.
    Variants,
    /// A statistic in windows, from bedGraph.
    Windows,
    /// Association statistics from a table.
    Manhattan,
    /// A phylogeny from Newick.
    Tree,
    /// An alignment from aligned FASTA.
    Msa,
    /// Variable sites of an alignment, from aligned FASTA.
    Snps,
    /// Cytogenetic bands from a cytoBand table.
    Ideogram,
    /// A value per sample per site, from a table.
    Matrix,
    /// Aligned reads from SAM text.
    Pileup,
    /// Alignment ribbons between two sequences, from PAF.
    Synteny,
    /// The same alignments as a dot plot, from the same PAF.
    Dotplot,
    /// Open reading frames in six frames, computed from FASTA bases.
    Orfs,
    /// A sequence logo, counted from aligned FASTA.
    Logo,
    /// Two phylogenies face to face, from two Newick files.
    Tanglegram,
    /// Spans carried by named taxa, painted onto a phylogeny.
    Clades,
    /// Gene neighbourhoods from several genomes, with their homologies drawn
    /// between them.
    Loci,
    /// The coordinate ruler, which reads nothing.
    Axis,
}

impl Kind {
    /// Every track the command line can draw, in the order the help text
    /// lists them.
    ///
    /// A caller offering a choice of tracks, a browser front end among them,
    /// wants the list rather than a copy of it that goes stale. The help text
    /// is checked against this, so a track added without a line in it is a
    /// failing test rather than a flag nobody can find.
    pub const ALL: [Kind; 20] = [
        Kind::Coverage,
        Kind::Sequence,
        Kind::Features,
        Kind::Variants,
        Kind::Windows,
        Kind::Manhattan,
        Kind::Tree,
        Kind::Msa,
        Kind::Snps,
        Kind::Ideogram,
        Kind::Matrix,
        Kind::Pileup,
        Kind::Synteny,
        Kind::Dotplot,
        Kind::Orfs,
        Kind::Logo,
        Kind::Tanglegram,
        Kind::Clades,
        Kind::Loci,
        Kind::Axis,
    ];

    /// The flag that asks for this track, without the dashes.
    pub fn flag(self) -> &'static str {
        match self {
            Kind::Coverage => "coverage",
            Kind::Sequence => "sequence",
            Kind::Features => "features",
            Kind::Variants => "variants",
            Kind::Windows => "windows",
            Kind::Manhattan => "manhattan",
            Kind::Tree => "tree",
            Kind::Msa => "msa",
            Kind::Snps => "snps",
            Kind::Ideogram => "ideogram",
            Kind::Matrix => "matrix",
            Kind::Pileup => "pileup",
            Kind::Synteny => "synteny",
            Kind::Dotplot => "dotplot",
            Kind::Orfs => "orfs",
            Kind::Logo => "logo",
            Kind::Tanglegram => "tanglegram",
            Kind::Clades => "clades",
            Kind::Loci => "loci",
            Kind::Axis => "axis",
        }
    }

    /// The flag as it is written, for the errors that name one.
    ///
    /// Exhaustive on purpose. [`Kind::flag`] spells the track and this spells
    /// the word someone typed, and when the two were kept in step by a
    /// fallback arm instead of by the compiler, four tracks added later fell
    /// through it and `--synteny` with no path reported that `--axis` needed a
    /// value.
    pub fn dashed(self) -> &'static str {
        match self {
            Kind::Coverage => "--coverage",
            Kind::Sequence => "--sequence",
            Kind::Features => "--features",
            Kind::Variants => "--variants",
            Kind::Windows => "--windows",
            Kind::Manhattan => "--manhattan",
            Kind::Tree => "--tree",
            Kind::Msa => "--msa",
            Kind::Snps => "--snps",
            Kind::Ideogram => "--ideogram",
            Kind::Matrix => "--matrix",
            Kind::Pileup => "--pileup",
            Kind::Synteny => "--synteny",
            Kind::Dotplot => "--dotplot",
            Kind::Orfs => "--orfs",
            Kind::Logo => "--logo",
            Kind::Tanglegram => "--tanglegram",
            Kind::Clades => "--clades",
            Kind::Loci => "--loci",
            Kind::Axis => "--axis",
        }
    }

    /// Whether `--aggregate` and `--log` mean anything here.
    fn takes_aggregate(self) -> bool {
        matches!(self, Kind::Coverage)
    }

    /// Whether the track has a height of its own, rather than one that follows
    /// from how many rows its data needs.
    fn takes_height(self) -> bool {
        matches!(
            self,
            Kind::Coverage
                | Kind::Sequence
                | Kind::Variants
                | Kind::Windows
                | Kind::Manhattan
                | Kind::Ideogram
                | Kind::Synteny
                | Kind::Dotplot
                | Kind::Axis
        )
    }

    /// The flag that names this track's second file, where one file is not all
    /// of its data.
    ///
    /// Spelled by what the second file means and not by where it sits. A
    /// second phylogeny and a second table are different files asking
    /// different questions, so one spelling for both would be a flag that is
    /// accepted everywhere and correct in one place, which is the thing the
    /// module refuses to do with every other modifier.
    ///
    /// A track that takes a second file needs it. Drawing one from the single
    /// file it was given is not a smaller version of the figure asked for, it
    /// is a different figure: a tanglegram of one tree against itself has no
    /// crossings and looks like a perfect result.
    ///
    /// Public because it gates a public field: a caller building a
    /// [`TrackSpec`] by hand, which is what a front end that is not a shell
    /// does, has to be able to ask which tracks want a second file and what to
    /// call the control that asks for it.
    pub fn second_flag(self) -> Option<&'static str> {
        // Exhaustive for the reason `dashed` is. A fallback arm here does not
        // give a wrong error message, it gives a track drawn without a file it
        // cannot do without, and each of the three has a figure that looks
        // finished when that happens.
        match self {
            Kind::Tanglegram => Some("--against"),
            Kind::Clades => Some("--with-tree"),
            Kind::Loci => Some("--links"),
            Kind::Coverage
            | Kind::Sequence
            | Kind::Features
            | Kind::Variants
            | Kind::Windows
            | Kind::Manhattan
            | Kind::Tree
            | Kind::Msa
            | Kind::Snps
            | Kind::Ideogram
            | Kind::Matrix
            | Kind::Pileup
            | Kind::Synteny
            | Kind::Dotplot
            | Kind::Orfs
            | Kind::Logo
            | Kind::Axis => None,
        }
    }
}

/// Either kind of style flag, since `--style` spells both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Style {
    /// A filled area under the values.
    Area,
    /// A line through them.
    Line,
    /// One bar per point.
    Bars,
    /// A step function, which only a window track draws.
    Steps,
}

impl Style {
    /// The coverage spelling, or `None` for a style coverage does not have.
    pub fn coverage(self) -> Option<CoverageStyle> {
        Some(match self {
            Style::Area => CoverageStyle::Area,
            Style::Line => CoverageStyle::Line,
            Style::Bars => CoverageStyle::Bars,
            Style::Steps => return None,
        })
    }

    /// The window spelling, or `None` for a style a window track has not got.
    pub fn window(self) -> Option<WindowStyle> {
        Some(match self {
            Style::Steps => WindowStyle::Steps,
            Style::Line => WindowStyle::Line,
            Style::Area | Style::Bars => return None,
        })
    }
}

/// One track, and everything written about it before the next one started.
#[derive(Debug, Clone)]
pub struct TrackSpec {
    /// Which track.
    pub kind: Kind,
    /// Where its data is, or `None` for the axis.
    pub source: Option<Source>,
    /// The other file, for a track whose data is not one file.
    ///
    /// Which flag fills it is [`Kind::second_flag`], and a track that has one
    /// is refused without it, so this is `None` only where it means nothing.
    pub second: Option<Source>,
    /// `--label`.
    pub label: Option<String>,
    /// `--height`.
    pub height: Option<f64>,
    /// `--aggregate`.
    pub aggregate: Option<Aggregate>,
    /// `--style`.
    pub style: Option<Style>,
    /// `--log`.
    pub log: bool,
    /// `--color`.
    pub color: Option<String>,
    /// `--format`, when the file cannot be told by looking at it.
    pub format: Option<Format>,
    /// `--identity`, when a homology file's third column could be either unit.
    pub identity: Option<Identity>,
}

impl TrackSpec {
    /// A track with nothing said about it yet.
    fn new(kind: Kind, source: Option<Source>) -> Self {
        TrackSpec {
            kind,
            source,
            second: None,
            label: None,
            height: None,
            aggregate: None,
            style: None,
            log: false,
            color: None,
            format: None,
            identity: None,
        }
    }
}

/// Which theme was asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Palette {
    /// The default.
    Light,
    /// A selected set of colours, not an inversion of the light one.
    Dark,
}

/// A whole command line, parsed and not yet acted on.
#[derive(Debug)]
pub struct Invocation {
    /// The region every track is drawn over.
    pub region: Region,
    /// The tracks, in the order they were written and will be drawn.
    pub tracks: Vec<TrackSpec>,
    /// `--title`.
    pub title: Option<String>,
    /// `--width`.
    pub width: Option<f64>,
    /// `--theme`.
    pub theme: Palette,
    /// Cleared by `--no-axis`.
    pub axis: bool,
    /// Cleared by `--no-region-label`.
    pub region_label: bool,
    /// `-o`, or standard output when absent.
    pub output: Option<PathBuf>,
}

/// What the command line asked for, which is not always a figure.
#[derive(Debug)]
pub enum Request {
    /// Draw this.
    Draw(Box<Invocation>),
    /// Print the help text.
    Help,
    /// Print the version.
    Version,
}

/// Reads `args`, which must not include the program name.
///
/// # Errors
///
/// Returns the first thing that does not fit the grammar. Nothing is read from
/// disk here, so an error means the command line was wrong rather than the data.
pub fn parse(args: &[String]) -> Result<Request, ArgError> {
    if args.is_empty() {
        return Err(ArgError::NoArguments);
    }
    if args.iter().any(|a| a == "-h" || a == "--help") {
        return Ok(Request::Help);
    }
    if args.iter().any(|a| a == "-V" || a == "--version") {
        return Ok(Request::Version);
    }

    let mut region: Option<Region> = None;
    let mut tracks: Vec<TrackSpec> = Vec::new();
    let mut title = None;
    let mut width = None;
    let mut theme = Palette::Light;
    let mut axis = true;
    let mut region_label = true;
    let mut output = None;

    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let mut value = |flag: &'static str| rest.next().ok_or(ArgError::MissingValue(flag));

        // A track flag starts a track. `--axis` is the one that reads nothing.
        let track = match arg.as_str() {
            "--coverage" => Some((Kind::Coverage, true)),
            "--sequence" => Some((Kind::Sequence, true)),
            "--features" => Some((Kind::Features, true)),
            "--variants" => Some((Kind::Variants, true)),
            "--windows" => Some((Kind::Windows, true)),
            "--manhattan" => Some((Kind::Manhattan, true)),
            "--tree" => Some((Kind::Tree, true)),
            "--msa" => Some((Kind::Msa, true)),
            "--snps" => Some((Kind::Snps, true)),
            "--ideogram" => Some((Kind::Ideogram, true)),
            "--matrix" => Some((Kind::Matrix, true)),
            "--pileup" => Some((Kind::Pileup, true)),
            "--synteny" => Some((Kind::Synteny, true)),
            "--dotplot" => Some((Kind::Dotplot, true)),
            "--orfs" => Some((Kind::Orfs, true)),
            "--logo" => Some((Kind::Logo, true)),
            "--tanglegram" => Some((Kind::Tanglegram, true)),
            "--clades" => Some((Kind::Clades, true)),
            "--loci" => Some((Kind::Loci, true)),
            "--axis" => Some((Kind::Axis, false)),
            _ => None,
        };
        if let Some((kind, reads)) = track {
            let source = if reads {
                let word = value(kind.dashed())?;
                Some(if word == "-" {
                    Source::Stdin
                } else {
                    Source::Path(PathBuf::from(word))
                })
            } else {
                None
            };
            if matches!(source, Some(Source::Stdin)) && stdin_taken(&tracks) {
                return Err(ArgError::StdinTwice);
            }
            if kind == Kind::Axis {
                axis = false;
            }
            tracks.push(TrackSpec::new(kind, source));
            continue;
        }

        match arg.as_str() {
            "--label" => {
                let text = value("--label")?.clone();
                last(&mut tracks, "--label")?.label = Some(text);
            }
            "--height" => {
                let text = value("--height")?;
                let px = text.parse::<f64>().map_err(|_| ArgError::BadValue {
                    flag: "--height",
                    given: text.clone(),
                    expected: "a number of pixels",
                })?;
                let track = last(&mut tracks, "--height")?;
                if !track.kind.takes_height() {
                    return Err(ArgError::WrongTrack {
                        flag: "--height",
                        track: track.kind.flag(),
                    });
                }
                track.height = Some(px);
            }
            "--aggregate" => {
                let text = value("--aggregate")?;
                let how = match text.as_str() {
                    "max" => Aggregate::Max,
                    "mean" => Aggregate::Mean,
                    "min" => Aggregate::Min,
                    _ => {
                        return Err(ArgError::BadValue {
                            flag: "--aggregate",
                            given: text.clone(),
                            expected: "max, mean or min",
                        })
                    }
                };
                let track = last(&mut tracks, "--aggregate")?;
                if !track.kind.takes_aggregate() {
                    return Err(ArgError::WrongTrack {
                        flag: "--aggregate",
                        track: track.kind.flag(),
                    });
                }
                track.aggregate = Some(how);
            }
            "--style" => {
                let text = value("--style")?;
                let style = match text.as_str() {
                    "area" => Style::Area,
                    "line" => Style::Line,
                    "bars" => Style::Bars,
                    "steps" => Style::Steps,
                    _ => {
                        return Err(ArgError::BadValue {
                            flag: "--style",
                            given: text.clone(),
                            expected: "area, line, bars or steps",
                        })
                    }
                };
                let track = last(&mut tracks, "--style")?;
                // The two tracks that take a style do not take the same words.
                let fits = match track.kind {
                    Kind::Coverage => style.coverage().is_some(),
                    Kind::Windows => style.window().is_some(),
                    _ => false,
                };
                if !fits {
                    return Err(ArgError::BadValue {
                        flag: "--style",
                        given: text.clone(),
                        expected: match track.kind {
                            Kind::Coverage => "area, line or bars for a coverage track",
                            Kind::Windows => "steps or line for a window track",
                            _ => "nothing: this track has no style",
                        },
                    });
                }
                track.style = Some(style);
            }
            "--log" => {
                let track = last(&mut tracks, "--log")?;
                if !track.kind.takes_aggregate() {
                    return Err(ArgError::WrongTrack {
                        flag: "--log",
                        track: track.kind.flag(),
                    });
                }
                track.log = true;
            }
            "--color" => {
                let text = value("--color")?.clone();
                // The colour goes into an SVG attribute as it was written, so a
                // value carrying a quote or an angle bracket would end the
                // attribute early and leave a document that will not parse. No
                // spelling of a paint value needs one of these characters.
                if text.contains(['"', '\'', '<', '>', '&']) {
                    return Err(ArgError::BadValue {
                        flag: "--color",
                        given: text,
                        expected: "a colour, as in '#d55e00'",
                    });
                }
                let track = last(&mut tracks, "--color")?;
                if !matches!(track.kind, Kind::Coverage | Kind::Features) {
                    return Err(ArgError::WrongTrack {
                        flag: "--color",
                        track: track.kind.flag(),
                    });
                }
                track.color = Some(text);
            }
            flag @ ("--against" | "--with-tree" | "--links") => {
                // One arm for every second path, because the mechanism is one
                // mechanism; only the spelling changes, and the spelling is
                // what says which file it is.
                let flag: &'static str = match flag {
                    "--with-tree" => "--with-tree",
                    "--links" => "--links",
                    _ => "--against",
                };
                let word = value(flag)?;
                // Checked before the track is borrowed, and against both
                // fields: a pipe can be read once, and it is now possible for
                // one track to ask for it twice by itself.
                let stdin = word == "-";
                if stdin && stdin_taken(&tracks) {
                    return Err(ArgError::StdinTwice);
                }
                let source = if stdin {
                    Source::Stdin
                } else {
                    Source::Path(PathBuf::from(word))
                };
                let track = last(&mut tracks, flag)?;
                if track.kind.second_flag() != Some(flag) {
                    return Err(ArgError::WrongTrack {
                        flag,
                        track: track.kind.flag(),
                    });
                }
                track.second = Some(source);
            }
            "--identity" => {
                let text = value("--identity")?;
                let unit = Identity::parse(text).ok_or_else(|| ArgError::BadValue {
                    flag: "--identity",
                    given: text.clone(),
                    expected: "percent or fraction",
                })?;
                let track = last(&mut tracks, "--identity")?;
                if track.kind != Kind::Loci {
                    return Err(ArgError::WrongTrack {
                        flag: "--identity",
                        track: track.kind.flag(),
                    });
                }
                track.identity = Some(unit);
            }
            "--format" => {
                let text = value("--format")?;
                let format = Format::parse(text).ok_or_else(|| ArgError::BadValue {
                    flag: "--format",
                    given: text.clone(),
                    expected: "bedgraph, depth, values, bed or gff3",
                })?;
                last(&mut tracks, "--format")?.format = Some(format);
            }
            "--title" => title = Some(value("--title")?.clone()),
            "--width" => {
                let text = value("--width")?;
                let px = text.parse::<f64>().map_err(|_| ArgError::BadValue {
                    flag: "--width",
                    given: text.clone(),
                    expected: "a number of pixels",
                })?;
                // A width becomes one column of pixels per band, so a figure
                // wider than any screen or page asks for an allocation that
                // fails rather than a figure. A small width is raised to a
                // drawable one further down, so only the top end is refused.
                if !px.is_finite() || px > MAX_WIDTH {
                    return Err(ArgError::BadValue {
                        flag: "--width",
                        given: text.clone(),
                        expected: "a number of pixels, at most 100000",
                    });
                }
                width = Some(px);
            }
            "--theme" => {
                let text = value("--theme")?;
                theme = match text.as_str() {
                    "light" => Palette::Light,
                    "dark" => Palette::Dark,
                    _ => {
                        return Err(ArgError::BadValue {
                            flag: "--theme",
                            given: text.clone(),
                            expected: "light or dark",
                        })
                    }
                };
            }
            "--no-axis" => axis = false,
            "--no-region-label" => region_label = false,
            "-o" | "--output" => output = Some(PathBuf::from(value("-o")?)),
            flag if flag.starts_with('-') && flag != "-" => {
                return Err(ArgError::UnknownFlag(flag.to_string()))
            }
            locus => {
                if region.is_some() {
                    return Err(ArgError::ExtraRegion(locus.to_string()));
                }
                let parsed = Region::parse(locus).map_err(ArgError::BadRegion)?;
                // A track that keeps one value per base of the window sizes its
                // buffer from the span, so a span past every sequence anyone
                // has is an allocation that fails rather than a figure.
                if parsed.len() > MAX_SPAN {
                    return Err(ArgError::HugeRegion {
                        given: locus.to_string(),
                        span: parsed.len(),
                    });
                }
                region = Some(parsed);
            }
        }
    }

    // Late, because the flag that fills it comes after the track flag and may
    // be anywhere before the next one.
    for spec in &tracks {
        if let Some(flag) = spec.kind.second_flag() {
            if spec.second.is_none() {
                return Err(ArgError::MissingSecond {
                    flag,
                    track: spec.kind.flag(),
                });
            }
        }
    }

    let region = region.ok_or(ArgError::NoRegion)?;
    Ok(Request::Draw(Box::new(Invocation {
        region,
        tracks,
        title,
        width,
        theme,
        axis,
        region_label,
        output,
    })))
}

/// Whether standard input has already been spoken for.
///
/// Both fields, since a track drawn from two files could otherwise ask for the
/// pipe twice and read the same text as each of them.
fn stdin_taken(tracks: &[TrackSpec]) -> bool {
    tracks
        .iter()
        .flat_map(|t| [t.source.as_ref(), t.second.as_ref()])
        .any(|source| matches!(source, Some(Source::Stdin)))
}

/// The track a modifier belongs to, which is the one before it.
fn last<'a>(
    tracks: &'a mut [TrackSpec],
    flag: &'static str,
) -> Result<&'a mut TrackSpec, ArgError> {
    tracks.last_mut().ok_or(ArgError::NoTrackYet(flag))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(line: &str) -> Vec<String> {
        line.split_whitespace().map(String::from).collect()
    }

    fn draw(line: &str) -> Invocation {
        match parse(&args(line)).unwrap() {
            Request::Draw(invocation) => *invocation,
            other => panic!("expected a figure, got {other:?}"),
        }
    }

    #[test]
    fn the_order_of_the_flags_is_the_order_of_the_stack() {
        let it = draw("chr1:1-1000 --features g.bed --coverage d.bg --variants v.vcf");
        let kinds: Vec<Kind> = it.tracks.iter().map(|t| t.kind).collect();
        assert_eq!(kinds, vec![Kind::Features, Kind::Coverage, Kind::Variants]);
    }

    #[test]
    fn a_modifier_lands_on_the_track_before_it() {
        let it = draw("chr1:1-1000 --coverage d.bg --label depth --features g.bed --label genes");
        assert_eq!(it.tracks[0].label.as_deref(), Some("depth"));
        assert_eq!(it.tracks[1].label.as_deref(), Some("genes"));
    }

    #[test]
    fn a_modifier_before_any_track_is_an_error() {
        let err = parse(&args("chr1:1-1000 --label depth")).unwrap_err();
        assert!(matches!(err, ArgError::NoTrackYet("--label")));
    }

    #[test]
    fn a_modifier_the_track_has_no_use_for_is_an_error() {
        let err = parse(&args("chr1:1-1000 --features g.bed --aggregate min")).unwrap_err();
        assert!(matches!(
            err,
            ArgError::WrongTrack {
                flag: "--aggregate",
                track: "features"
            }
        ));
    }

    #[test]
    fn the_region_can_sit_anywhere_but_only_once() {
        let it = draw("--coverage d.bg chr1:1-1000");
        assert_eq!(it.region.seq(), "chr1");
        let err = parse(&args("chr1:1-1000 chr2:1-1000")).unwrap_err();
        assert!(matches!(err, ArgError::ExtraRegion(_)));
    }

    /// The grammar gives one path per flag, and a tanglegram is two trees. The
    /// second arrives by a name rather than by a position, because a second
    /// phylogeny and a second table are different files, and a modifier that
    /// means one thing here and another there is the thing this module refuses
    /// to have.
    /// `ALL` is a list, and a list is the kind of thing that goes stale. The
    /// compiler checks its length and nothing else, so this checks that every
    /// entry is a different track and that each one is a flag the parser
    /// actually answers to.
    #[test]
    fn the_list_of_tracks_is_every_track_exactly_once() {
        let mut spellings: Vec<&str> = Kind::ALL.iter().map(|k| k.dashed()).collect();
        spellings.sort_unstable();
        spellings.dedup();
        assert_eq!(spellings.len(), Kind::ALL.len(), "a track is listed twice");

        for kind in Kind::ALL {
            let line = format!("chr1:1-1000 {} f.txt --against t.nwk", kind.dashed());
            let parsed = parse(&args(&line));
            // Every one of them starts a track, so none is an unknown flag.
            assert!(
                !matches!(parsed, Err(ArgError::UnknownFlag(_))),
                "{} is listed and the parser does not know it",
                kind.dashed()
            );
        }
    }

    #[test]
    fn a_second_path_arrives_by_name_and_lands_on_the_track_before_it() {
        let it = draw("chr1:1-1000 --tanglegram before.nwk --against after.nwk --label topology");
        assert_eq!(it.tracks.len(), 1, "the second path started a second track");
        assert_eq!(
            it.tracks[0].source,
            Some(Source::Path("before.nwk".into())),
            "the first path moved"
        );
        assert_eq!(
            it.tracks[0].second,
            Some(Source::Path("after.nwk".into())),
            "the second path did not arrive"
        );
        // And an ordinary modifier still lands after it, so the second path
        // did not end the track.
        assert_eq!(it.tracks[0].label.as_deref(), Some("topology"));
    }

    /// A tanglegram of one tree against itself has no crossings, which is what
    /// a perfect result looks like, so the missing half has to be refused
    /// rather than filled in. The check is late because `--against` may sit
    /// anywhere before the next track flag.
    #[test]
    fn a_track_drawn_from_two_files_is_refused_with_one() {
        let error = parse(&args("chr1:1-1000 --tanglegram before.nwk")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::MissingSecond {
                    flag: "--against",
                    track: "tanglegram"
                }
            ),
            "{error}"
        );

        // Written after other flags, which is the whole reason for a late pass.
        let it = draw("chr1:1-1000 --tanglegram a.nwk --label x --against b.nwk --title t");
        assert_eq!(it.tracks[0].second, Some(Source::Path("b.nwk".into())));
    }

    /// Every other modifier is refused by name where it says nothing, and this
    /// one is no different: a coverage track has no second file, so `--against`
    /// on one is a flag that would be read and then ignored.
    #[test]
    fn the_second_path_flag_is_refused_where_a_track_has_no_second_file() {
        let error = parse(&args("chr1:1-1000 --coverage d.bg --against t.nwk")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--against",
                    track: "coverage"
                }
            ),
            "{error}"
        );

        let error = parse(&args("chr1:1-1000 --against t.nwk")).unwrap_err();
        assert!(
            matches!(error, ArgError::NoTrackYet("--against")),
            "{error}"
        );
    }

    /// A pipe can be read once, and a track drawn from two files is the first
    /// thing that can ask for it twice on its own.
    #[test]
    fn one_track_cannot_take_standard_input_for_both_of_its_files() {
        let error = parse(&args("chr1:1-1000 --tanglegram - --against -")).unwrap_err();
        assert!(matches!(error, ArgError::StdinTwice), "{error}");

        // And the second path is counted against the other tracks too, in
        // both directions.
        let error = parse(&args(
            "chr1:1-1000 --tanglegram a.nwk --against - --coverage -",
        ));
        assert!(matches!(error, Err(ArgError::StdinTwice)), "{error:?}");
        let error = parse(&args(
            "chr1:1-1000 --coverage - --tanglegram a.nwk --against -",
        ));
        assert!(matches!(error, Err(ArgError::StdinTwice)), "{error:?}");

        // One of the two, on the other hand, is the ordinary case.
        let it = draw("chr1:1-1000 --tanglegram - --against after.nwk");
        assert_eq!(it.tracks[0].source, Some(Source::Stdin));
    }

    /// [`Kind::flag`] spells the track and [`Kind::dashed`] spells the word,
    /// and a fallback arm once let four tracks fall through so that
    /// `--synteny` with no path reported that `--axis` needed a value.
    #[test]
    fn a_track_flag_missing_its_path_names_itself() {
        for flag in [
            "--synteny",
            "--dotplot",
            "--orfs",
            "--logo",
            "--tanglegram",
            "--clades",
            "--loci",
        ] {
            let error = parse(&args(&format!("chr1:1-1000 {flag}"))).unwrap_err();
            let ArgError::MissingValue(named) = error else {
                panic!("{flag} without a path gave {error}");
            };
            assert_eq!(named, flag, "{flag} reported itself as {named}");
        }
    }

    /// The mechanism is one mechanism and the spellings are three, because a
    /// second phylogeny, a second table of names and the right-hand tree of a
    /// tanglegram are three different files. Each is refused by name anywhere
    /// it says nothing, which is what the spelling is for.
    #[test]
    fn every_track_drawn_from_two_files_asks_for_the_second_by_its_own_name() {
        for (track, flag, other) in [
            ("--tanglegram", "--against", "--links"),
            ("--clades", "--with-tree", "--against"),
            ("--loci", "--links", "--with-tree"),
        ] {
            let it = draw(&format!("chr1:1-1000 {track} one.txt {flag} two.txt"));
            assert_eq!(it.tracks.len(), 1, "{flag} started a second track");
            assert_eq!(
                it.tracks[0].second,
                Some(Source::Path("two.txt".into())),
                "{track} did not take {flag}"
            );

            // Without it, refused by name rather than drawn from one file.
            let error = parse(&args(&format!("chr1:1-1000 {track} one.txt"))).unwrap_err();
            let ArgError::MissingSecond { flag: named, .. } = error else {
                panic!("{track} with one file gave {error}");
            };
            assert_eq!(named, flag);

            // And another track's spelling means nothing here.
            let error = parse(&args(&format!(
                "chr1:1-1000 {track} one.txt {other} two.txt"
            )))
            .unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { flag, .. } if flag == other),
                "{track} accepted {other}: {error}"
            );

            // One pipe, whichever of the two files asks for it.
            let error = parse(&args(&format!("chr1:1-1000 {track} - {flag} -"))).unwrap_err();
            assert!(matches!(error, ArgError::StdinTwice), "{track}: {error}");
        }
    }

    /// `--identity` says whether a homology file's third column is a
    /// percentage or a fraction, which no other track has an opinion about.
    #[test]
    fn the_identity_unit_means_nothing_to_a_track_that_is_not_loci() {
        let it = draw("chr1:1-1000 --loci g.bed --links h.tsv --identity fraction");
        assert_eq!(it.tracks[0].identity, Some(Identity::Fraction));

        let error = parse(&args("chr1:1-1000 --coverage d.bg --identity percent")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--identity",
                    ..
                }
            ),
            "{error}"
        );

        let error = parse(&args(
            "chr1:1-1000 --loci g.bed --links h.tsv --identity half",
        ))
        .unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::BadValue {
                    flag: "--identity",
                    ..
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn a_dash_means_standard_input_and_only_one_track_gets_it() {
        let it = draw("chr1:1-1000 --coverage -");
        assert_eq!(it.tracks[0].source, Some(Source::Stdin));
        let err = parse(&args("chr1:1-1000 --coverage - --features -")).unwrap_err();
        assert!(matches!(err, ArgError::StdinTwice));
    }

    #[test]
    fn an_explicit_axis_cancels_the_automatic_one() {
        assert!(draw("chr1:1-1000 --coverage d.bg").axis);
        assert!(!draw("chr1:1-1000 --axis --coverage d.bg").axis);
        assert!(!draw("chr1:1-1000 --coverage d.bg --no-axis").axis);
    }

    #[test]
    fn figure_flags_are_not_attached_to_a_track() {
        let it = draw("chr1:1-1000 --title one --coverage d.bg --width 1200 --theme dark");
        assert_eq!(it.title.as_deref(), Some("one"));
        assert_eq!(it.width, Some(1200.0));
        assert_eq!(it.theme, Palette::Dark);
        assert!(it.tracks[0].label.is_none());
    }

    #[test]
    fn a_word_a_flag_does_not_take_says_what_it_does_take() {
        let err = parse(&args("chr1:1-1000 --coverage d.bg --aggregate median")).unwrap_err();
        assert_eq!(
            err.to_string(),
            "--aggregate does not take \"median\", only max, mean or min"
        );
    }

    #[test]
    fn a_flag_with_no_value_says_so() {
        let err = parse(&args("chr1:1-1000 --coverage")).unwrap_err();
        assert!(matches!(err, ArgError::MissingValue("--coverage")));
    }

    #[test]
    fn help_and_version_win_over_everything_else() {
        assert!(matches!(parse(&args("--help")).unwrap(), Request::Help));
        assert!(matches!(
            parse(&args("nonsense --version")).unwrap(),
            Request::Version
        ));
    }

    #[test]
    fn a_colour_carrying_an_xml_metacharacter_is_refused() {
        // The value is written into an SVG attribute as it stands, so this is
        // attribute injection and not merely a malformed file: `blue"` closed
        // the attribute early, and a value carrying `<` went further and put
        // elements of the caller's choosing into the document.
        for given in ["blue\"", "<script>", "red\" onload=\"x", "a&b", "it's"] {
            let args = vec![
                "chr1:1-1000".to_string(),
                "--coverage".to_string(),
                "d.bg".to_string(),
                "--color".to_string(),
                given.to_string(),
            ];
            let err = parse(&args).unwrap_err();
            assert!(
                matches!(
                    &err,
                    ArgError::BadValue {
                        flag: "--color",
                        ..
                    }
                ),
                "{given:?} was taken: {err}"
            );
        }
    }

    #[test]
    fn every_spelling_of_a_colour_still_goes_through() {
        for given in ["#d55e00", "blue", "rgb(0,0,255)", "url(#g)"] {
            let args = vec![
                "chr1:1-1000".to_string(),
                "--coverage".to_string(),
                "d.bg".to_string(),
                "--color".to_string(),
                given.to_string(),
            ];
            match parse(&args).unwrap() {
                Request::Draw(it) => assert_eq!(it.tracks[0].color.as_deref(), Some(given)),
                other => panic!("expected a figure, got {other:?}"),
            }
        }
    }

    #[test]
    fn a_width_no_renderer_could_draw_is_an_error_rather_than_a_panic() {
        // 1e308 used to reach an allocation of one column per pixel and abort
        // the process with `capacity overflow`.
        for given in ["1e308", "inf", "nan", "1e15", "100001"] {
            let err = parse(&args(&format!(
                "chr1:1-200 --coverage d.bg --width {given}"
            )))
            .unwrap_err();
            assert!(
                matches!(
                    &err,
                    ArgError::BadValue {
                        flag: "--width",
                        ..
                    }
                ),
                "{given:?} was taken: {err}"
            );
        }
        assert_eq!(draw("chr1:1-200 --width 100000").width, Some(100_000.0));
    }

    #[test]
    fn a_region_longer_than_a_figure_is_drawn_over_is_an_error() {
        // The whole u64 range used to reach `vec![0.0; region.len() as usize]`
        // and abort the process with `capacity overflow`.
        let err = parse(&args("chr1:1-18446744073709551615 --coverage d.bg")).unwrap_err();
        assert_eq!(
            err.to_string(),
            "\"chr1:1-18446744073709551615\" spans 18446744073709551615 bases, \
             and a figure is drawn over at most 4294967296"
        );
        // A whole large sequence is an ordinary figure and stays one.
        assert_eq!(draw("chr1:1-248956422").region.len(), 248_956_422);
    }

    #[test]
    fn a_missing_region_is_its_own_message() {
        let err = parse(&args("--coverage d.bg")).unwrap_err();
        assert!(err.to_string().contains("the first argument is the region"));
    }
}
