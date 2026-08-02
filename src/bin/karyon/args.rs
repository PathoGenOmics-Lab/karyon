//! The command line, which is the stack.
//!
//! A figure is a list of tracks in the order they are drawn, and `argv` is a
//! list in the order it was written, so the grammar is the obvious one: each
//! `--<track>` flag starts a track, and the flags after it describe that track
//! until the next one starts. It is [`karyon::Plot`] with spaces instead of
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
//! fall outside that. One is a track whose data is not one file: a `--<track>`
//! flag takes one path, so a tanglegram, which is two trees, has nowhere to put
//! the second, and thirteen of the crate's twenty-nine track types are all the
//! command line reaches.
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

use karyon::{Aggregate, CoverageStyle, Region, WindowStyle};

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
    BadRegion(karyon::Error),
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
    /// The coordinate ruler, which reads nothing.
    Axis,
}

impl Kind {
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
            Kind::Axis => "axis",
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
                | Kind::Axis
        )
    }
}

/// How a file should be read, when the shape of it is not enough to tell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// `chrom start end value`, 0-based half-open.
    BedGraph,
    /// `chrom pos depth`, the 1-based output of `samtools depth`.
    Depth,
    /// One value per line, starting at the left edge of the region.
    Values,
    /// `chrom start end [name] [score] [strand]`, 0-based half-open.
    Bed,
    /// Nine columns, 1-based inclusive, attributes in the ninth.
    Gff3,
}

impl Format {
    /// Parses the word given to `--format`.
    fn parse(word: &str) -> Option<Format> {
        Some(match word {
            "bedgraph" | "bg" => Format::BedGraph,
            "depth" => Format::Depth,
            "values" => Format::Values,
            "bed" => Format::Bed,
            "gff3" | "gff" | "gtf" => Format::Gff3,
            _ => return None,
        })
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
}

impl TrackSpec {
    /// A track with nothing said about it yet.
    fn new(kind: Kind, source: Option<Source>) -> Self {
        TrackSpec {
            kind,
            source,
            label: None,
            height: None,
            aggregate: None,
            style: None,
            log: false,
            color: None,
            format: None,
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
            "--axis" => Some((Kind::Axis, false)),
            _ => None,
        };
        if let Some((kind, reads)) = track {
            let source = if reads {
                let word = value(leak(kind.flag()))?;
                Some(if word == "-" {
                    Source::Stdin
                } else {
                    Source::Path(PathBuf::from(word))
                })
            } else {
                None
            };
            if matches!(source, Some(Source::Stdin))
                && tracks
                    .iter()
                    .any(|t| matches!(t.source, Some(Source::Stdin)))
            {
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

/// The track a modifier belongs to, which is the one before it.
fn last<'a>(
    tracks: &'a mut [TrackSpec],
    flag: &'static str,
) -> Result<&'a mut TrackSpec, ArgError> {
    tracks.last_mut().ok_or(ArgError::NoTrackYet(flag))
}

/// A track flag's name, for the one error message that needs it as a literal.
fn leak(flag: &str) -> &'static str {
    // The set is closed and small, so a match beats leaking a String.
    match flag {
        "coverage" => "--coverage",
        "sequence" => "--sequence",
        "features" => "--features",
        "variants" => "--variants",
        "windows" => "--windows",
        "manhattan" => "--manhattan",
        "tree" => "--tree",
        "msa" => "--msa",
        "snps" => "--snps",
        "ideogram" => "--ideogram",
        "matrix" => "--matrix",
        "pileup" => "--pileup",
        _ => "--axis",
    }
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
