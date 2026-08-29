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
//! Twenty-eight of the crate's thirty-six track types are what the command
//! line reaches.
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
use crate::track::tree::TreeProjection;
use crate::track::MsaDisplay;
use crate::{Aggregate, CoverageStyle, Format, Region, VariantStyle, WindowStyle};

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
    /// A list of columns on a track that was given no sheet to take them from.
    ///
    /// A copy number track with no ploidy to read its levels against.
    ///
    /// Required rather than defaulted, because where balanced sits is not in
    /// the file and this crate does not know what it is drawing. Two copies is
    /// right for a human autosome and wrong for most of what else is handed to
    /// it, and a rule in the wrong place does not merely mis-scale the ladder,
    /// it swaps every gain for a loss and says so confidently.
    MissingPloidy,
    /// Checked once the whole line is read rather than where the flag sits, so
    /// that `--columns` before `--traits` and after it are the same command.
    /// Two modifiers of one track are not in an order, and refusing one of the
    /// two spellings would be refusing a line that says exactly what it means.
    Unsourced {
        /// The track that has columns and no sheet.
        track: &'static str,
    },
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
/// is what a figure costs in memory: eight bytes a base, and
/// [`CoverageTrack::from_pairs`](crate::CoverageTrack::from_pairs) allocates
/// the whole of it at once.
///
/// The number is the largest such buffer a 32-bit target can hold. A `Vec` may
/// not exceed `isize::MAX` bytes, so `1 << 28` values of eight bytes each is
/// exactly one byte too many, and asking for it is a capacity overflow rather
/// than an allocation failure: a panic, and in a build that aborts on one, a
/// trap. That is what this limit is for, and the old one, `1 << 32`, was
/// thirty-two gigabytes and above the ceiling on every target there is.
///
/// It leaves every chromosome anyone sequences inside it. The longest human
/// one is 248,956,422 bases, which is eight per cent under this.
const MAX_SPAN: u64 = (1 << 28) - 1;

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
                write!(f, "{flag} means nothing to {} {track} track", article(track))
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
            ArgError::MissingPloidy => write!(
                f,
                "--copy-number needs --ploidy, since where balanced sits is not in the file"
            ),
            ArgError::Unsourced { track } => write!(
                f,
                "--columns picks out of the sheet --traits names, and this {track} track was given no sheet"
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
    /// Segmented copy number from a caller's segment table.
    CopyNumber,
    /// Per-base model attribution from a bedGraph, over a reference.
    Dynseq,
    /// Splice junctions from an aligner's SJ.out.tab.
    Junctions,
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
    /// Modified bases per strand, from a bedMethyl pileup.
    Methylation,
    /// Structural calls as arcs between their breakpoints, from a VCF.
    Structural,
    /// Molecules that aligned in pieces, from SAM and its SA tag.
    SplitReads,
    /// Methylation one molecule at a time, from a Bismark extractor file.
    Bisulfite,
    /// Protein domains, from an InterProScan table.
    Domains,
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
    pub const ALL: [Kind; 28] = [
        Kind::Coverage,
        Kind::CopyNumber,
        Kind::Dynseq,
        Kind::Junctions,
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
        Kind::Methylation,
        Kind::Structural,
        Kind::SplitReads,
        Kind::Bisulfite,
        Kind::Domains,
        Kind::Axis,
    ];

    /// The flag that asks for this track, without the dashes.
    pub fn flag(self) -> &'static str {
        match self {
            Kind::Coverage => "coverage",
            Kind::CopyNumber => "copy-number",
            Kind::Dynseq => "dynseq",
            Kind::Junctions => "junctions",
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
            Kind::Methylation => "methylation",
            Kind::Structural => "structural",
            Kind::SplitReads => "split-reads",
            Kind::Bisulfite => "bisulfite",
            Kind::Domains => "domains",
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
            Kind::CopyNumber => "--copy-number",
            Kind::Dynseq => "--dynseq",
            Kind::Junctions => "--junctions",
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
            Kind::Methylation => "--methylation",
            Kind::Structural => "--structural",
            Kind::SplitReads => "--split-reads",
            Kind::Bisulfite => "--bisulfite",
            Kind::Domains => "--domains",
            Kind::Axis => "--axis",
        }
    }

    /// Whether `--aggregate` and `--log` mean anything here.
    fn takes_aggregate(self) -> bool {
        matches!(self, Kind::Coverage)
    }

    /// Whether the track has a height of its own, rather than one that follows
    /// from how many rows its data needs.
    /// Whether a sheet of metadata means anything to this track.
    ///
    /// The tracks drawn as a row per named thing, which are the ones a strip
    /// can sit beside and line up with. A pileup has rows too, and they are
    /// reads rather than samples: nobody keeps a sheet keyed by read name, and
    /// a flag accepted there would be a flag that draws nothing.
    pub fn takes_traits(self) -> bool {
        matches!(
            self,
            Kind::Matrix | Kind::Msa | Kind::Snps | Kind::Clades | Kind::Domains | Kind::Loci
        )
    }

    /// Whether `--max-rows` means anything here.
    ///
    /// The tracks that stack one row per record and cap themselves. A track
    /// whose rows come from its data rather than from a cap, a feature track
    /// packing what fits, is not one of them: it has no cap to move.
    ///
    /// A tree is here and answers differently from the other four. They stop
    /// opening rows and count what they left out, which a tree cannot do: a
    /// tip is not interchangeable with the tip below it and cutting the list
    /// would cut a clade in half. So it collapses the smallest clades instead
    /// until it fits, and every tip stays on the figure inside a triangle that
    /// says how many it holds.
    fn takes_max_rows(self) -> bool {
        matches!(
            self,
            Kind::Pileup | Kind::Msa | Kind::Snps | Kind::Bisulfite | Kind::Tree
        )
    }

    /// Whether `--no-names` means anything here.
    ///
    /// The tracks that write a name for each thing they draw, which is not the
    /// same set as the tracks that have rows: an ideogram writes a name on each
    /// band and a feature track writes one on each feature. The builder behind
    /// it is called `show_names` in most of them, `show_row_names` in the
    /// matrix and `show_band_names` in the ideogram, and the flag is named for
    /// what the reader sees rather than for what the three of them are called.
    ///
    /// A coverage track has nothing to name, and a track whose only name is the
    /// one in the gutter is turned off with `--label ''` instead.
    ///
    /// An ideogram is not here, and it is the one case worth explaining.
    /// `IdeogramTrack::show_band_names` defaults to false, so there are no band
    /// names on the figure for this flag to take off, and accepting it would be
    /// a flag that draws the same picture. What the ideogram is missing is the
    /// opposite flag, one that turns the names on, and this grammar has no
    /// positive form: `--no-axis` and `--no-region-label` both take away
    /// something that is there by default. Adding one is a decision about the
    /// shape of the grammar rather than about this flag.
    fn takes_no_names(self) -> bool {
        matches!(
            self,
            Kind::Features
                | Kind::Msa
                | Kind::Snps
                | Kind::Matrix
                | Kind::SplitReads
                | Kind::Structural
                | Kind::Bisulfite
                | Kind::Domains
                | Kind::Loci
                | Kind::Clades
        )
    }

    /// Whether `--threshold` means anything here.
    ///
    /// One track, because one track draws a line across itself and colours
    /// what stands above it. The other numbers that decide what a track
    /// believes are floors on the data rather than a line on the figure, and
    /// giving them the same name would say they were the same thing.
    /// Whether `--projection` means anything here.
    ///
    /// One track. A tanglegram is two trees facing each other and a clade
    /// track paints spans onto one, and both of those are rectangular by
    /// construction: there is no second axis to bend into a circle.
    fn takes_projection(self) -> bool {
        matches!(self, Kind::Tree)
    }

    /// Whether `--focus` means anything here.
    fn takes_focus(self) -> bool {
        matches!(self, Kind::Tree)
    }

    /// Whether `--compare-to` means anything here.
    ///
    /// The two tracks that read every row against one of them. An alignment
    /// compares against the consensus until told otherwise, which its own
    /// documentation calls deliberate; a variable-site panel compares against
    /// whichever record the file happened to hold first, which is not a
    /// decision anyone made. A logo takes a nought in stack.rs too and is not
    /// here: that nought is a start coordinate.
    fn takes_compare_to(self) -> bool {
        matches!(self, Kind::Msa | Kind::Snps)
    }

    /// Whether `--no-counts` means anything here.
    ///
    /// The two tracks that print a number beside the thing it counts: the
    /// per-sample difference count down the right of a variable-site panel and
    /// the read count over each junction arc. Named for what the reader sees,
    /// the way `--no-names` is, and for the same reason: they are the same mark
    /// on the page even where they count different things.
    fn takes_no_counts(self) -> bool {
        matches!(self, Kind::Snps | Kind::Junctions)
    }

    /// Whether `--min-reads` means anything here.
    ///
    /// The two tracks with a read floor of their own, both of which say how
    /// many they dropped. The library calls them `min_coverage` and
    /// `min_reads`, and the difference is only the local noun: after
    /// `--methylation` the reads behind the site, after `--junctions` the reads
    /// across the intron. A structural track is not here, because support is a
    /// per-record setter there and the track has no floor, so the flag would
    /// filter nothing.
    fn takes_min_reads(self) -> bool {
        matches!(self, Kind::Methylation | Kind::Junctions)
    }

    /// Whether `--row-height` means anything here.
    ///
    /// The tracks whose height is a number of rows times the height of one,
    /// which is most of what `Kind::takes_height` refuses for that reason. The
    /// builder is `row_height` in ten of them, `read_height` on a pileup and
    /// `lane_height` on an ORF track, and the flag is named for what the
    /// reader sees, the way `--no-names` is over its three spellings.
    ///
    /// A locus track is not here and is the case worth explaining. Its row is
    /// two heights, `gene_height` and `link_height`, and one number cannot set
    /// both without this flag deciding the ratio between them, which is a
    /// decision about the figure rather than about its size.
    fn takes_row_height(self) -> bool {
        matches!(
            self,
            Kind::Features
                | Kind::Msa
                | Kind::Snps
                | Kind::Matrix
                | Kind::Pileup
                | Kind::Orfs
                | Kind::Tree
                | Kind::Tanglegram
                | Kind::Clades
                | Kind::SplitReads
                | Kind::Bisulfite
                | Kind::Domains
        )
    }

    fn takes_threshold(self) -> bool {
        matches!(self, Kind::Manhattan)
    }

    fn takes_height(self) -> bool {
        matches!(
            self,
            Kind::Coverage
                | Kind::CopyNumber
                | Kind::Dynseq
                | Kind::Sequence
                | Kind::Variants
                | Kind::Windows
                | Kind::Manhattan
                | Kind::Ideogram
                | Kind::Synteny
                | Kind::Dotplot
                | Kind::Methylation
                | Kind::Structural
                | Kind::Junctions
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
    /// A second file this track will take and can do without.
    ///
    /// A pileup colours a base that disagrees with the reference, and the
    /// reference is the only place it can learn what a base should be. Without
    /// one it draws every read agreeing, which is a figure worth having and is
    /// what every pileup this command line has ever drawn looks like, so the
    /// file is offered rather than demanded. [`Kind::second_flag`] is the other
    /// half of this: the tracks there are refused without their second file.
    pub fn optional_second(self) -> Option<&'static str> {
        match self {
            Kind::Pileup => Some("--with-sequence"),
            _ => None,
        }
    }

    /// The flag that fills [`TrackSpec::second`] for a track that cannot be
    /// drawn without it, and `None` for one that can.
    pub fn second_flag(self) -> Option<&'static str> {
        // Exhaustive for the reason `dashed` is. A fallback arm here does not
        // give a wrong error message, it gives a track drawn without a file it
        // cannot do without, and each of the three has a figure that looks
        // finished when that happens.
        match self {
            Kind::Dynseq => Some("--with-sequence"),
            Kind::Tanglegram => Some("--against"),
            Kind::Clades => Some("--with-tree"),
            Kind::Loci => Some("--links"),
            Kind::Coverage
            | Kind::CopyNumber
            | Kind::Junctions
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
            | Kind::Methylation
            | Kind::Structural
            | Kind::SplitReads
            | Kind::Bisulfite
            | Kind::Domains
            | Kind::Axis => None,
        }
    }

    /// The flag that says which of the several things a file holds to draw.
    ///
    /// Three formats hold more than one dataset in one file and one of them is
    /// a track: a pileup counting two modifications of the same cytosine, an
    /// extractor file carrying three sequence contexts, a domain table where a
    /// dozen member databases describe the same kinase in their own words.
    /// Drawn together they are several marks in one place naming none of them.
    ///
    /// Spelled by what is being chosen, as [`Kind::second_flag`] is, and
    /// exhaustive for the same reason: a track added without one silently
    /// stacks whatever it found.
    pub fn selector(self) -> Option<&'static str> {
        match self {
            Kind::Methylation => Some("--modification"),
            Kind::Bisulfite => Some("--context"),
            Kind::Domains => Some("--analysis"),
            Kind::Coverage
            | Kind::CopyNumber
            | Kind::Dynseq
            | Kind::Junctions
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
            | Kind::Tanglegram
            | Kind::Clades
            | Kind::Loci
            | Kind::Structural
            | Kind::SplitReads
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
    /// A plain vertical mark, which only a variant track draws. The right
    /// choice once the calls are dense enough that lollipop heads would
    /// overlap into a smear, and the only way a command line reader could
    /// reach the library's answer to that.
    Tick,
    /// A stem with a head, which only a variant track draws, and its default.
    Lollipop,
    /// Only the cells that disagree with the row an alignment is compared
    /// against, which is what it draws until it is told otherwise.
    ///
    /// Not to be confused with [`Style::ALL`], the associated constant one
    /// keystroke away, which is the list of every style. This one is a style.
    Differences,
    /// Every cell of an alignment, agreements included.
    All,
}

impl Style {
    /// Every style the parser accepts, with the word it is written as.
    ///
    /// Here rather than in the help text, so a style wired up and left out of
    /// the help is a failing test instead of a value nobody can find.
    pub const ALL: [(Style, &'static str); 8] = [
        (Style::Area, "area"),
        (Style::Line, "line"),
        (Style::Bars, "bars"),
        (Style::Steps, "steps"),
        (Style::Tick, "tick"),
        (Style::Lollipop, "lollipop"),
        (Style::Differences, "differences"),
        (Style::All, "all"),
    ];

    /// The coverage spelling, or `None` for a style coverage does not have.
    pub fn coverage(self) -> Option<CoverageStyle> {
        Some(match self {
            Style::Area => CoverageStyle::Area,
            Style::Line => CoverageStyle::Line,
            Style::Bars => CoverageStyle::Bars,
            Style::Steps | Style::Tick | Style::Lollipop | Style::Differences | Style::All => {
                return None
            }
        })
    }

    /// The window spelling, or `None` for a style a window track has not got.
    pub fn window(self) -> Option<WindowStyle> {
        Some(match self {
            Style::Steps => WindowStyle::Steps,
            Style::Line => WindowStyle::Line,
            Style::Area
            | Style::Bars
            | Style::Tick
            | Style::Lollipop
            | Style::Differences
            | Style::All => return None,
        })
    }

    /// The variant spelling, or `None` for a style a variant track has not got.
    pub fn variant(self) -> Option<VariantStyle> {
        Some(match self {
            Style::Tick => VariantStyle::Tick,
            Style::Lollipop => VariantStyle::Lollipop,
            Style::Area
            | Style::Line
            | Style::Bars
            | Style::Steps
            | Style::Differences
            | Style::All => return None,
        })
    }

    /// The alignment spelling, or `None` for a style an alignment has not got.
    ///
    /// The word is `all` rather than the library's `Bases`, because the track
    /// draws proteins as readily as nucleotides and the variant's own doc
    /// comment says "every residue", which is the library correcting its own
    /// name. A library name may be loose; a word a reader types may not.
    pub fn msa(self) -> Option<MsaDisplay> {
        Some(match self {
            Style::Differences => MsaDisplay::Differences,
            Style::All => MsaDisplay::Bases,
            Style::Area
            | Style::Line
            | Style::Bars
            | Style::Steps
            | Style::Tick
            | Style::Lollipop => return None,
        })
    }
}

/// The article a track name takes, which is not decided by its first letter.
///
/// Most of them are plain: a coverage track, a features track. Four begin with
/// a vowel and take `an`, and two more are read out letter by letter and begin
/// with a vowel sound while beginning with a consonant on the page: an msa
/// track is said "an em-ess-ay track", and an snps track "an ess-en-pees
/// track". Spelling this from the first character alone gets those two wrong,
/// which is why they are named here rather than derived.
fn article(track: &str) -> &'static str {
    let vowel = matches!(track.chars().next(), Some('a' | 'e' | 'i' | 'o' | 'u'));
    if vowel || matches!(track, "msa" | "snps") {
        "an"
    } else {
        "a"
    }
}

/// What `--max-rows` was given.
///
/// The tracks that stack rows cap themselves, at forty, and say how many they
/// left out. The cap is there because a pile a thousand deep draws off the
/// bottom of the page, and it is worth moving because forty is a guess about
/// the reader's screen rather than about their data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowCap {
    /// A number of rows, at least one.
    Rows(usize),
    /// `all`, which lifts the cap and draws every row there is.
    All,
}

impl RowCap {
    /// Parses what `--max-rows` takes.
    pub fn parse(word: &str) -> Option<RowCap> {
        if word == "all" {
            return Some(RowCap::All);
        }
        match word.parse::<usize>() {
            Ok(rows) if rows >= 1 => Some(RowCap::Rows(rows)),
            _ => None,
        }
    }

    /// What the track builders take, where `None` is no cap at all.
    pub fn rows(self) -> Option<usize> {
        match self {
            RowCap::Rows(rows) => Some(rows),
            RowCap::All => None,
        }
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
    /// Which of the several things a file holds to draw, named by the flag
    /// [`Kind::selector`] gives this track.
    pub selects: Option<String>,
    /// `--max-rows`, for the tracks that stack their data into rows.
    pub max_rows: Option<RowCap>,
    /// `--no-names`, which leaves out the names a track writes beside or on
    /// what it draws. Not the track's own name, which is `--label`.
    pub no_names: bool,
    /// `--threshold`, the line a scan is read against.
    pub threshold: Option<f64>,
    /// `--compare-to`, the row every other row is read against.
    pub compare_to: Option<String>,
    /// `--projection`, the shape a phylogeny is laid out in.
    pub projection: Option<TreeProjection>,
    /// `--focus`, the one clade of a phylogeny to draw.
    pub focus: Option<Vec<String>>,
    /// `--no-counts`, which leaves out the number printed beside the thing it
    /// counts.
    pub no_counts: bool,
    /// `--min-reads`, the fewest reads behind a record for it to be drawn.
    pub min_reads: Option<u32>,
    /// `--fade-by-mapq`, which draws a read at less than full strength the
    /// lower its mapping quality.
    pub fade_by_mapq: bool,
    /// `--row-height`, for the tracks whose height follows from their rows.
    ///
    /// The complement of [`TrackSpec::height`], and the two never both apply:
    /// `Kind::takes_height` and `Kind::takes_row_height` are disjoint.
    pub row_height: Option<f64>,
    /// `--ploidy`, where balanced sits on a copy number ladder.
    pub ploidy: Option<f64>,
    /// `--sample`, for a file holding more than one.
    pub sample: Option<String>,
    /// `--traits`, the sample sheet whose columns are drawn beside the rows.
    pub traits: Option<Source>,
    /// `--columns`, the columns of that sheet to draw and the order to draw
    /// them in. `None` draws every column the sheet has, in the order its
    /// header named them.
    pub columns: Option<Vec<String>>,
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
            selects: None,
            max_rows: None,
            no_names: false,
            threshold: None,
            row_height: None,
            compare_to: None,
            projection: None,
            focus: None,
            no_counts: false,
            min_reads: None,
            fade_by_mapq: false,
            ploidy: None,
            sample: None,
            traits: None,
            columns: None,
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
            "--copy-number" => Some((Kind::CopyNumber, true)),
            "--dynseq" => Some((Kind::Dynseq, true)),
            "--junctions" => Some((Kind::Junctions, true)),
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
            "--methylation" => Some((Kind::Methylation, true)),
            "--structural" => Some((Kind::Structural, true)),
            "--split-reads" => Some((Kind::SplitReads, true)),
            "--bisulfite" => Some((Kind::Bisulfite, true)),
            "--domains" => Some((Kind::Domains, true)),
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
            "--ploidy" => {
                let text = value("--ploidy")?;
                let copies = text
                    .parse::<f64>()
                    .ok()
                    // Greater than nought, not merely not negative. At nought
                    // a log ratio becomes nought copies everywhere, and every
                    // one of them lands on the rule that means unchanged, so a
                    // called amplification is drawn as a quiet arm.
                    .filter(|copies| copies.is_finite() && *copies > 0.0)
                    .ok_or_else(|| ArgError::BadValue {
                        flag: "--ploidy",
                        given: text.clone(),
                        expected: "a number of copies above nought, as in 2",
                    })?;
                let track = last(&mut tracks, "--ploidy")?;
                if track.kind != Kind::CopyNumber {
                    return Err(ArgError::WrongTrack {
                        flag: "--ploidy",
                        track: track.kind.flag(),
                    });
                }
                track.ploidy = Some(copies);
            }
            "--sample" => {
                let text = value("--sample")?.clone();
                let track = last(&mut tracks, "--sample")?;
                if track.kind != Kind::CopyNumber {
                    return Err(ArgError::WrongTrack {
                        flag: "--sample",
                        track: track.kind.flag(),
                    });
                }
                track.sample = Some(text);
            }
            "--traits" => {
                let word = value("--traits")?.clone();
                // Checked before the track is borrowed, and against every
                // source: a pipe can be read once, and a third file per track
                // is a third way for one command to ask for it twice.
                let stdin = word == "-";
                if stdin && stdin_taken(&tracks) {
                    return Err(ArgError::StdinTwice);
                }
                let source = if stdin {
                    Source::Stdin
                } else {
                    Source::Path(PathBuf::from(word))
                };
                let track = last(&mut tracks, "--traits")?;
                if !track.kind.takes_traits() {
                    return Err(ArgError::WrongTrack {
                        flag: "--traits",
                        track: track.kind.flag(),
                    });
                }
                track.traits = Some(source);
            }
            "--columns" => {
                let text = value("--columns")?.clone();
                let track = last(&mut tracks, "--columns")?;
                if !track.kind.takes_traits() {
                    return Err(ArgError::WrongTrack {
                        flag: "--columns",
                        track: track.kind.flag(),
                    });
                }
                let wanted: Vec<String> = text
                    .split(',')
                    .map(str::trim)
                    .filter(|name| !name.is_empty())
                    .map(str::to_string)
                    .collect();
                if wanted.is_empty() {
                    return Err(ArgError::BadValue {
                        flag: "--columns",
                        given: text,
                        expected: "one or more column names separated by commas",
                    });
                }
                track.columns = Some(wanted);
            }
            "--no-names" => {
                let track = last(&mut tracks, "--no-names")?;
                if !track.kind.takes_no_names() {
                    return Err(ArgError::WrongTrack {
                        flag: "--no-names",
                        track: track.kind.flag(),
                    });
                }
                track.no_names = true;
            }
            "--threshold" => {
                let text = value("--threshold")?;
                // The word is here because the number behind it is
                // -log10(5e-8), which nobody types from memory, and because
                // asking for it by name is the only way to ask for it: the
                // track has no default threshold, deliberately, since a
                // Bonferroni correction for a million tests is the wrong
                // number wherever a million tests were not run.
                let value = if text == "genome-wide" {
                    -(5e-8f64).log10()
                } else {
                    match text.parse::<f64>() {
                        Ok(number) if number.is_finite() => number,
                        _ => {
                            return Err(ArgError::BadValue {
                                flag: "--threshold",
                                given: text.clone(),
                                expected: "a number, or genome-wide",
                            })
                        }
                    }
                };
                let track = last(&mut tracks, "--threshold")?;
                if !track.kind.takes_threshold() {
                    return Err(ArgError::WrongTrack {
                        flag: "--threshold",
                        track: track.kind.flag(),
                    });
                }
                track.threshold = Some(value);
            }
            "--max-rows" => {
                let text = value("--max-rows")?;
                let cap = RowCap::parse(text).ok_or_else(|| ArgError::BadValue {
                    flag: "--max-rows",
                    given: text.clone(),
                    expected: "a number of rows, or all",
                })?;
                let track = last(&mut tracks, "--max-rows")?;
                if !track.kind.takes_max_rows() {
                    return Err(ArgError::WrongTrack {
                        flag: "--max-rows",
                        track: track.kind.flag(),
                    });
                }
                track.max_rows = Some(cap);
            }
            "--projection" => {
                let text = value("--projection")?;
                let projection = match text.as_str() {
                    "rectangular" => TreeProjection::Rectangular,
                    "circular" => TreeProjection::Circular,
                    "unrooted" => TreeProjection::Unrooted,
                    _ => {
                        return Err(ArgError::BadValue {
                            flag: "--projection",
                            given: text.clone(),
                            expected: "rectangular, circular or unrooted",
                        })
                    }
                };
                let track = last(&mut tracks, "--projection")?;
                if !track.kind.takes_projection() {
                    return Err(ArgError::WrongTrack {
                        flag: "--projection",
                        track: track.kind.flag(),
                    });
                }
                track.projection = Some(projection);
            }
            "--focus" => {
                let text = value("--focus")?;
                let names: Vec<String> = text
                    .split(',')
                    .map(str::trim)
                    .filter(|name| !name.is_empty())
                    .map(str::to_string)
                    .collect();
                if names.is_empty() || names.len() > 2 {
                    return Err(ArgError::BadValue {
                        flag: "--focus",
                        given: text.clone(),
                        expected: "a clade name, a tip name, or two tip names separated by a comma",
                    });
                }
                let track = last(&mut tracks, "--focus")?;
                if !track.kind.takes_focus() {
                    return Err(ArgError::WrongTrack {
                        flag: "--focus",
                        track: track.kind.flag(),
                    });
                }
                track.focus = Some(names);
            }
            "--compare-to" => {
                let name = value("--compare-to")?.clone();
                let track = last(&mut tracks, "--compare-to")?;
                if !track.kind.takes_compare_to() {
                    return Err(ArgError::WrongTrack {
                        flag: "--compare-to",
                        track: track.kind.flag(),
                    });
                }
                track.compare_to = Some(name);
            }
            "--no-counts" => {
                let track = last(&mut tracks, "--no-counts")?;
                if !track.kind.takes_no_counts() {
                    return Err(ArgError::WrongTrack {
                        flag: "--no-counts",
                        track: track.kind.flag(),
                    });
                }
                track.no_counts = true;
            }
            "--min-reads" => {
                let text = value("--min-reads")?;
                let reads = text.parse::<u32>().map_err(|_| ArgError::BadValue {
                    flag: "--min-reads",
                    given: text.clone(),
                    expected: "a number of reads",
                })?;
                let track = last(&mut tracks, "--min-reads")?;
                if !track.kind.takes_min_reads() {
                    return Err(ArgError::WrongTrack {
                        flag: "--min-reads",
                        track: track.kind.flag(),
                    });
                }
                track.min_reads = Some(reads);
            }
            "--fade-by-mapq" => {
                // One track, so no predicate: the flags several kinds share
                // have one, and the flags one kind has compare inline, the way
                // `--ploidy` and `--identity` do.
                let track = last(&mut tracks, "--fade-by-mapq")?;
                if track.kind != Kind::Pileup {
                    return Err(ArgError::WrongTrack {
                        flag: "--fade-by-mapq",
                        track: track.kind.flag(),
                    });
                }
                track.fade_by_mapq = true;
            }
            "--row-height" => {
                let text = value("--row-height")?;
                let px = text
                    .parse::<f64>()
                    .ok()
                    // Refused rather than clamped, and refused at both ends.
                    // Every builder behind this flag ends in `height.max(n)`,
                    // and `f64::max` lets an infinity through and answers with
                    // the floor for a NaN. An infinite row height draws a
                    // figure whose clip is one pixel tall and whose rows are
                    // empty groups, and a NaN is indistinguishable from asking
                    // for the minimum. Both exit nought, which is the worst
                    // way for either to end.
                    .filter(|px| px.is_finite() && *px > 0.0)
                    .ok_or_else(|| ArgError::BadValue {
                        flag: "--row-height",
                        given: text.clone(),
                        expected: "a number of pixels above nought, as in 20",
                    })?;
                let track = last(&mut tracks, "--row-height")?;
                if !track.kind.takes_row_height() {
                    return Err(ArgError::WrongTrack {
                        flag: "--row-height",
                        track: track.kind.flag(),
                    });
                }
                track.row_height = Some(px);
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
                let style = match Style::ALL.iter().find(|(_, word)| *word == text) {
                    Some((style, _)) => *style,
                    None => {
                        return Err(ArgError::BadValue {
                            flag: "--style",
                            given: text.clone(),
                            expected: "area, line, bars, steps, tick, lollipop, differences or all",
                        })
                    }
                };
                let track = last(&mut tracks, "--style")?;
                // A track with no styles at all is the wrong track for the
                // flag rather than a track given the wrong word, and saying so
                // needs a different error. It used to fall through to the one
                // below and print `--style does not take "bars", only nothing:
                // this track has no style`, which is not a sentence.
                if !matches!(
                    track.kind,
                    Kind::Coverage | Kind::Windows | Kind::Variants | Kind::Msa
                ) {
                    return Err(ArgError::WrongTrack {
                        flag: "--style",
                        track: track.kind.flag(),
                    });
                }
                // The tracks that take a style do not take the same words.
                let fits = match track.kind {
                    Kind::Coverage => style.coverage().is_some(),
                    Kind::Windows => style.window().is_some(),
                    Kind::Variants => style.variant().is_some(),
                    Kind::Msa => style.msa().is_some(),
                    _ => unreachable!("every styled track is named above"),
                };
                if !fits {
                    return Err(ArgError::BadValue {
                        flag: "--style",
                        given: text.clone(),
                        expected: match track.kind {
                            Kind::Coverage => "area, line or bars for a coverage track",
                            Kind::Windows => "steps or line for a window track",
                            Kind::Variants => "tick or lollipop for a variant track",
                            Kind::Msa => "differences or all for an alignment track",
                            _ => unreachable!("every styled track is named above"),
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
                if !matches!(
                    track.kind,
                    Kind::Coverage | Kind::Features | Kind::Junctions
                ) {
                    return Err(ArgError::WrongTrack {
                        flag: "--color",
                        track: track.kind.flag(),
                    });
                }
                track.color = Some(text);
            }
            flag @ ("--against" | "--with-tree" | "--links" | "--with-sequence") => {
                // One arm for every second path, because the mechanism is one
                // mechanism; only the spelling changes, and the spelling is
                // what says which file it is.
                let flag: &'static str = match flag {
                    "--with-tree" => "--with-tree",
                    "--links" => "--links",
                    "--with-sequence" => "--with-sequence",
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
                if track.kind.second_flag() != Some(flag)
                    && track.kind.optional_second() != Some(flag)
                {
                    return Err(ArgError::WrongTrack {
                        flag,
                        track: track.kind.flag(),
                    });
                }
                track.second = Some(source);
            }
            flag @ ("--modification" | "--context" | "--analysis") => {
                let flag: &'static str = match flag {
                    "--context" => "--context",
                    "--analysis" => "--analysis",
                    _ => "--modification",
                };
                let chosen = value(flag)?.clone();
                let track = last(&mut tracks, flag)?;
                if track.kind.selector() != Some(flag) {
                    return Err(ArgError::WrongTrack {
                        flag,
                        track: track.kind.flag(),
                    });
                }
                track.selects = Some(chosen);
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
        // Late for the same reason `--against` is checked late: the flag may
        // sit anywhere after the track it describes.
        if spec.kind == Kind::CopyNumber && spec.ploidy.is_none() {
            return Err(ArgError::MissingPloidy);
        }
        // Late for the same reason: `--columns` may be written before the
        // `--traits` it picks from, and both orders describe one track.
        if spec.columns.is_some() && spec.traits.is_none() {
            return Err(ArgError::Unsourced {
                track: spec.kind.flag(),
            });
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
        .flat_map(|t| [t.source.as_ref(), t.second.as_ref(), t.traits.as_ref()])
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
    fn the_flags_that_take_no_value_land_only_where_they_mean_something() {
        for (flag, good, bad) in [
            (
                "--no-counts",
                ["--snps a.fa", "--junctions j.tab"],
                ["--msa a.fa", "--coverage d.bg"],
            ),
            (
                "--fade-by-mapq",
                ["--pileup r.sam", "--pileup -"],
                ["--msa a.fa", "--split-reads s.sam"],
            ),
        ] {
            for track in good {
                let line = format!("chr1:1-1000 {track} {flag}");
                assert!(parse(&args(&line)).is_ok(), "{track} should take {flag}");
            }
            for track in bad {
                let line = format!("chr1:1-1000 {track} {flag}");
                let error = parse(&args(&line)).unwrap_err();
                assert!(
                    matches!(&error, ArgError::WrongTrack { flag: f, .. } if *f == flag),
                    "{track} should refuse {flag}, got {error:?}"
                );
            }
        }
    }

    #[test]
    fn a_read_floor_is_a_whole_number_of_reads() {
        let it = draw("chr1:1-1000 --junctions j.tab --min-reads 12");
        assert_eq!(it.tracks[0].min_reads, Some(12));
        assert_eq!(
            draw("chr1:1-1000 --junctions j.tab").tracks[0].min_reads,
            None
        );

        for word in ["-1", "2.5", "many"] {
            let line = format!("chr1:1-1000 --junctions j.tab --min-reads {word}");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::BadValue { flag, .. } if flag == "--min-reads"),
                "{word} should be refused, got {error:?}"
            );
        }
        let error = parse(&args("chr1:1-1000 --structural v.vcf --min-reads 5")).unwrap_err();
        assert!(
            matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--min-reads"),
            "a structural track has no floor to raise, {error:?}"
        );
    }

    #[test]
    fn an_alignment_takes_two_style_words_and_no_others() {
        for word in ["differences", "all"] {
            let line = format!("chr1:1-1000 --msa a.fa --style {word}");
            assert!(
                parse(&args(&line)).is_ok(),
                "an alignment should take {word}"
            );
        }
        for word in ["area", "bars", "tick", "steps"] {
            let line = format!("chr1:1-1000 --msa a.fa --style {word}");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::BadValue { flag, .. } if flag == "--style"),
                "{word} is not an alignment style, got {error:?}"
            );
        }
        // A track with no styles at all is the wrong track for the flag, not a
        // track handed the wrong word, and it used to be told the second thing
        // in a sentence that did not parse as English.
        let error = parse(&args("chr1:1-1000 --logo a.fa --style bars")).unwrap_err();
        assert_eq!(error.to_string(), "--style means nothing to a logo track");
    }

    #[test]
    fn only_the_tracks_that_compare_rows_take_a_row_to_compare_against() {
        let it = draw("chr1:1-1000 --msa a.fa --compare-to H37Rv");
        assert_eq!(it.tracks[0].compare_to.as_deref(), Some("H37Rv"));
        for bad in ["--matrix m.tsv", "--logo a.fa", "--pileup r.sam"] {
            let line = format!("chr1:1-1000 {bad} --compare-to X");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--compare-to"),
                "{bad} should refuse --compare-to, got {error:?}"
            );
        }
    }

    #[test]
    fn a_row_height_is_a_finite_number_of_pixels_above_nought() {
        let it = draw("chr1:1-1000 --msa a.fa --row-height 24");
        assert_eq!(it.tracks[0].row_height, Some(24.0));
        assert_eq!(draw("chr1:1-1000 --msa a.fa").tracks[0].row_height, None);

        // All four are refused rather than clamped. Every builder behind this
        // flag ends in `height.max(n)`, which lets an infinity through and
        // answers with the floor for a NaN, so a clamp would draw a figure
        // nobody asked for and say nothing.
        for word in ["inf", "-inf", "NaN", "0", "-5", "tall"] {
            let line = format!("chr1:1-1000 --msa a.fa --row-height {word}");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::BadValue { flag, .. } if flag == "--row-height"),
                "{word} should be refused, got {error:?}"
            );
        }
    }

    #[test]
    fn height_and_row_height_do_not_overlap() {
        // Two flags for one thing would be a grammar that cannot be read. A
        // track sizes itself by rows or it takes a height, and no track does
        // both, so every kind answers to exactly one of them or to neither.
        for kind in Kind::ALL {
            assert!(
                !(kind.takes_height() && kind.takes_row_height()),
                "{:?} takes both --height and --row-height",
                kind
            );
        }
        let error = parse(&args("chr1:1-1000 --coverage d.bg --row-height 20")).unwrap_err();
        assert!(
            matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--row-height"),
            "{error:?}"
        );
        let error = parse(&args("chr1:1-1000 --msa a.fa --height 200")).unwrap_err();
        assert!(
            matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--height"),
            "{error:?}"
        );
    }

    #[test]
    fn a_threshold_is_a_number_or_the_convention_by_name() {
        let it = draw("chr1:1-1000 --manhattan a.tsv --threshold 5.5");
        assert_eq!(it.tracks[0].threshold, Some(5.5));

        // The word is worth having only if it is the number the field means by
        // it, so this pins the number and not just that a word was accepted.
        let named = draw("chr1:1-1000 --manhattan a.tsv --threshold genome-wide");
        let wanted = -(5e-8f64).log10();
        assert!((named.tracks[0].threshold.unwrap() - wanted).abs() < 1e-12);
        assert!((wanted - 7.301_029_995_663_981).abs() < 1e-9);

        assert_eq!(
            draw("chr1:1-1000 --manhattan a.tsv").tracks[0].threshold,
            None
        );
    }

    #[test]
    fn a_threshold_that_is_not_a_number_is_refused() {
        for word in ["lots", "genome_wide", "genomewide", "inf", "NaN", ""] {
            let line = format!("chr1:1-1000 --manhattan a.tsv --threshold {word}");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(&error, ArgError::BadValue { flag, .. } if *flag == "--threshold")
                    || matches!(&error, ArgError::MissingValue { .. }),
                "{word:?} should be refused, got {error:?}"
            );
        }
        let error = parse(&args("chr1:1-1000 --coverage d.bg --threshold 5")).unwrap_err();
        assert!(
            matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--threshold"),
            "{error:?}"
        );
    }

    #[test]
    fn a_pileup_takes_a_reference_and_does_not_demand_one() {
        // Every pileup drawn from a command line before now had no reference,
        // so demanding one would refuse the commands that already exist.
        let plain = draw("chr1:1-1000 --pileup r.sam");
        assert_eq!(plain.tracks[0].second, None);

        let with = draw("chr1:1-1000 --pileup r.sam --with-sequence ref.fa");
        assert!(with.tracks[0].second.is_some());
        assert_eq!(Kind::Pileup.optional_second(), Some("--with-sequence"));
        assert_eq!(Kind::Pileup.second_flag(), None);

        // A dynseq track is the other half: the same flag, and refused without.
        assert_eq!(Kind::Dynseq.second_flag(), Some("--with-sequence"));
        let error = parse(&args("chr1:1-1000 --dynseq s.bg")).unwrap_err();
        assert!(matches!(error, ArgError::MissingSecond { .. }), "{error:?}");
    }

    #[test]
    fn a_second_file_still_lands_only_where_it_belongs() {
        let error = parse(&args("chr1:1-1000 --coverage d.bg --with-sequence ref.fa")).unwrap_err();
        assert!(
            matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--with-sequence"),
            "{error:?}"
        );
    }

    #[test]
    fn a_refusal_names_the_track_with_the_right_article() {
        let refusal = |line: &str| parse(&args(line)).unwrap_err().to_string();
        assert_eq!(
            refusal("chr1:1-1000 --msa a.fa --aggregate max"),
            "--aggregate means nothing to an msa track"
        );
        assert_eq!(
            refusal("chr1:1-1000 --snps a.fa --aggregate max"),
            "--aggregate means nothing to an snps track"
        );
        assert_eq!(
            refusal("chr1:1-1000 --ideogram b.txt --no-names"),
            "--no-names means nothing to an ideogram track"
        );
        assert_eq!(
            refusal("chr1:1-1000 --logo a.fa --log"),
            "--log means nothing to a logo track"
        );
    }

    #[test]
    fn only_the_tracks_that_write_names_take_no_names() {
        for good in [
            "--msa a.fa",
            "--snps a.fa",
            "--matrix m.tsv",
            "--features g.bed",
            "--split-reads s.sam",
            "--structural v.vcf",
            "--bisulfite c.txt",
            "--domains d.tsv",
        ] {
            let line = format!("chr1:1-1000 {good} --no-names");
            assert!(parse(&args(&line)).is_ok(), "{good} should take --no-names");
        }
        // A coverage track has nothing to name, and an ideogram draws no band
        // names to begin with, so taking them off would draw the same figure.
        for bad in [
            "--coverage d.bg",
            "--ideogram b.txt",
            "--pileup r.sam",
            "--logo a.fa",
        ] {
            let line = format!("chr1:1-1000 {bad} --no-names");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--no-names"),
                "{bad} should refuse --no-names, got {error:?}"
            );
        }
    }

    /// The playground's own list of the flags that open a track, held against
    /// this one.
    ///
    /// `docs/assets/karyon-wasm.js` needs to know where one track's words end
    /// and the next track's begin, or a control cannot tell one track's
    /// `--style` from another's. It cannot ask the parser: the controls are
    /// drawn before the program has finished loading, so the list is written
    /// out there and this is what stops it falling behind. Both directions,
    /// for the reason the help text is checked both ways: a track added to the
    /// parser and not to the page fails here, and so does a spelling on the
    /// page that the parser has never heard of.
    #[test]
    fn the_playground_knows_every_track_the_parser_does() {
        const SHIM: &str = include_str!("../../docs/assets/karyon-wasm.js");

        let opened = SHIM
            .split_once("var TRACKS = [")
            .expect("the playground has no TRACKS list")
            .1
            .split_once("];")
            .expect("the TRACKS list does not end")
            .0;
        let listed: Vec<&str> = opened
            .split('"')
            .filter(|piece| piece.starts_with("--"))
            .collect();

        for kind in Kind::ALL {
            assert!(
                listed.contains(&kind.dashed()),
                "the playground's TRACKS does not know {}",
                kind.dashed()
            );
        }
        for word in &listed {
            assert!(
                Kind::ALL.iter().any(|kind| kind.dashed() == *word),
                "the playground's TRACKS has {word}, which opens no track here"
            );
        }
        assert_eq!(
            listed.len(),
            Kind::ALL.len(),
            "the playground lists {} tracks and the parser has {}",
            listed.len(),
            Kind::ALL.len()
        );
    }

    #[test]
    fn a_row_cap_is_a_number_or_the_word_all() {
        let it = draw("chr1:1-1000 --pileup r.sam --max-rows 12");
        assert_eq!(it.tracks[0].max_rows, Some(RowCap::Rows(12)));
        assert_eq!(it.tracks[0].max_rows.unwrap().rows(), Some(12));

        let lifted = draw("chr1:1-1000 --pileup r.sam --max-rows all");
        assert_eq!(lifted.tracks[0].max_rows, Some(RowCap::All));
        // `None` is what the track builders take for no cap at all, so the
        // word has to arrive as that and not as a very large number.
        assert_eq!(lifted.tracks[0].max_rows.unwrap().rows(), None);

        assert_eq!(draw("chr1:1-1000 --pileup r.sam").tracks[0].max_rows, None);
    }

    #[test]
    fn a_row_cap_of_none_is_refused_rather_than_drawn_empty() {
        // Zero rows is a figure with the data left out, and the caller almost
        // certainly meant `all`. The builders clamp it to one, which would
        // draw one row and say nothing about the other thousand.
        for word in ["0", "-3", "lots", "40.5"] {
            let error = parse(&args(&format!(
                "chr1:1-1000 --pileup r.sam --max-rows {word}"
            )))
            .unwrap_err();
            assert!(
                matches!(error, ArgError::BadValue { flag, .. } if flag == "--max-rows"),
                "{word} should be refused, got {error:?}"
            );
        }
    }

    #[test]
    fn only_the_tracks_that_cap_themselves_take_a_row_cap() {
        for good in [
            "--pileup r.sam",
            "--msa a.fa",
            "--snps a.fa",
            "--bisulfite c.txt",
            // A tree takes it and answers differently: it collapses the
            // smallest clades until it fits rather than dropping the rest,
            // because a tip is not interchangeable with the tip below it and
            // cutting the list would cut a clade in half.
            "--tree t.nwk",
        ] {
            let line = format!("chr1:1-1000 {good} --max-rows 5");
            assert!(parse(&args(&line)).is_ok(), "{good} should take --max-rows");
        }
        // A feature track packs into as many rows as its data needs and has no
        // cap to move, so the flag would be accepted and do nothing.
        for bad in [
            "--features g.bed",
            "--coverage d.bg",
            "--logo a.fa",
            "--tanglegram t.nwk",
        ] {
            let line = format!("chr1:1-1000 {bad} --max-rows 5");
            let error = parse(&args(&line)).unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { flag, .. } if flag == "--max-rows"),
                "{bad} should refuse --max-rows, got {error:?}"
            );
        }
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

        // A sheet is a third file per track, so it is a third way to ask for
        // the pipe twice. The list `stdin_taken` walks is an array literal,
        // and a field left out of it is a silent second read of a drained
        // pipe, which comes back as an empty file rather than as an error.
        let error = parse(&args("chr1:1-1000 --matrix - --traits -"));
        assert!(matches!(error, Err(ArgError::StdinTwice)), "{error:?}");
        let error = parse(&args("chr1:1-1000 --matrix m.tsv --traits - --coverage -"));
        assert!(matches!(error, Err(ArgError::StdinTwice)), "{error:?}");
        let error = parse(&args("chr1:1-1000 --coverage - --matrix m.tsv --traits -"));
        assert!(matches!(error, Err(ArgError::StdinTwice)), "{error:?}");

        // One of the two, on the other hand, is the ordinary case.
        let it = draw("chr1:1-1000 --tanglegram - --against after.nwk");
        assert_eq!(it.tracks[0].source, Some(Source::Stdin));
        let it = draw("chr1:1-1000 --matrix - --traits sheet.tsv");
        assert_eq!(it.tracks[0].source, Some(Source::Stdin));
    }

    /// Where balanced sits is not in a segment table, and a rule in the wrong
    /// place does not mis-scale a copy number ladder, it inverts it.
    #[test]
    fn a_copy_number_track_is_refused_without_a_ploidy() {
        let error = parse(&args("chr8:1-1000 --copy-number seg.cns")).unwrap_err();
        assert!(matches!(error, ArgError::MissingPloidy), "{error}");

        // Late, so the flag may sit anywhere after the track it describes.
        let it = draw("chr8:1-1000 --copy-number seg.cns --label copies --ploidy 2");
        assert_eq!(it.tracks[0].ploidy, Some(2.0));

        let error = parse(&args("chr8:1-1000 --copy-number seg.cns --ploidy two")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::BadValue {
                    flag: "--ploidy",
                    ..
                }
            ),
            "{error}"
        );

        let error = parse(&args("chr8:1-1000 --coverage d.bg --ploidy 2")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--ploidy",
                    track: "coverage"
                }
            ),
            "{error}"
        );
    }

    #[test]
    fn a_sheet_of_metadata_is_refused_by_a_track_that_has_no_rows() {
        // The refusal is the point. A coverage track has nothing for a strip
        // to line up with, and a flag accepted and then ignored gives a figure
        // that is not the one asked for and does not look wrong.
        let error = parse(&args("chr1:1-1000 --coverage d.bg --traits sheet.tsv")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--traits",
                    track: "coverage"
                }
            ),
            "{error}"
        );

        for line in [
            "chr1:1-1000 --matrix m.tsv --traits s.tsv",
            "chr1:1-1000 --msa a.fa --traits s.tsv",
            "chr1:1-1000 --snps a.fa --traits s.tsv",
            "chr1:1-1000 --clades b.gff --with-tree t.nwk --traits s.tsv",
            "chr1:1-1000 --domains d.tsv --traits s.tsv",
            "chr1:1-1000 --loci l.gff --links h.tsv --traits s.tsv",
        ] {
            let it = draw(line);
            assert_eq!(
                it.tracks[0].traits,
                Some(Source::Path(PathBuf::from("s.tsv"))),
                "{line}"
            );
        }
    }

    #[test]
    fn a_ploidy_of_nought_is_refused_and_the_height_is_not() {
        // At nought copies a log ratio becomes nought everywhere, and every one
        // of them lands on the rule that means unchanged.
        let error = parse(&args("chr8:1-1000 --copy-number s.cns --ploidy 0")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::BadValue {
                    flag: "--ploidy",
                    ..
                }
            ),
            "{error}"
        );

        // Both new tracks size themselves by a field rather than by rows, so
        // both take a height, and a flag refused where it means something is
        // as wrong as one accepted where it does not.
        for line in [
            "chr8:1-1000 --copy-number s.cns --ploidy 2 --height 90",
            "chr1:1-1000 --dynseq d.bg --with-sequence r.fa --height 90",
        ] {
            let it = draw(line);
            assert_eq!(it.tracks[0].height, Some(90.0), "{line}");
        }
    }

    #[test]
    fn the_columns_of_a_sheet_may_be_named_before_the_sheet_is() {
        // Two modifiers of one track are not in an order, so the check that
        // one needs the other waits until the whole line has been read.
        let it = draw("chr1:1-1000 --matrix m.tsv --columns host,lineage --traits s.tsv");
        assert_eq!(
            it.tracks[0].columns,
            Some(vec!["host".to_string(), "lineage".to_string()])
        );

        let error = parse(&args("chr1:1-1000 --matrix m.tsv --columns host")).unwrap_err();
        assert!(
            matches!(error, ArgError::Unsourced { track: "matrix" }),
            "{error}"
        );

        let error = parse(&args("chr1:1-1000 --coverage d.bg --columns host")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--columns",
                    track: "coverage"
                }
            ),
            "{error}"
        );

        let error = parse(&args(
            "chr1:1-1000 --matrix m.tsv --traits s.tsv --columns ,,",
        ))
        .unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::BadValue {
                    flag: "--columns",
                    ..
                }
            ),
            "{error}"
        );
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
            "--methylation",
            "--structural",
            "--split-reads",
            "--bisulfite",
            "--domains",
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

    /// Three formats hold more than one dataset in one file and one of them is
    /// a track, and each names what it is choosing among. The mechanism is one
    /// mechanism, and each spelling is refused everywhere it says nothing.
    #[test]
    fn every_track_that_picks_one_dataset_names_what_it_is_picking_among() {
        for (track, flag, other) in [
            ("--methylation", "--modification", "--context"),
            ("--bisulfite", "--context", "--analysis"),
            ("--domains", "--analysis", "--modification"),
        ] {
            let it = draw(&format!("chr1:1-1000 {track} f.txt {flag} chosen"));
            assert_eq!(it.tracks[0].selects.as_deref(), Some("chosen"), "{track}");

            let error = parse(&args(&format!("chr1:1-1000 {track} f.txt {other} x"))).unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { flag, .. } if flag == other),
                "{track} accepted {other}: {error}"
            );

            // And a track with nothing to choose among refuses all three.
            let error = parse(&args(&format!("chr1:1-1000 --coverage d.bg {flag} x"))).unwrap_err();
            assert!(
                matches!(error, ArgError::WrongTrack { .. }),
                "{flag} landed on a coverage track: {error}"
            );
        }
    }

    /// Neither of the last two sizes itself from a height, so a height flag on
    /// one would parse and go nowhere, and the figure would come out at the
    /// default with nothing saying the number was dropped.
    #[test]
    fn the_last_two_tracks_take_no_second_path_and_no_height() {
        for flag in ["--bisulfite", "--domains"] {
            let it = draw(&format!("chr1:1-1000 {flag} f.txt"));
            assert_eq!(it.tracks[0].second, None);

            for refused in ["--against x", "--height 90"] {
                let error =
                    parse(&args(&format!("chr1:1-1000 {flag} f.txt {refused}"))).unwrap_err();
                assert!(
                    matches!(error, ArgError::WrongTrack { .. }),
                    "{flag} accepted {refused}: {error}"
                );
            }
        }
    }

    /// The three that came last take one file each, so none of them may be
    /// given a second, and only the two that do not size themselves from their
    /// data take a height.
    #[test]
    fn the_single_file_tracks_take_no_second_path_and_the_right_heights() {
        for flag in ["--methylation", "--structural", "--split-reads"] {
            let it = draw(&format!("chr1:1-1000 {flag} f.txt"));
            assert_eq!(it.tracks[0].second, None);
            let error = parse(&args(&format!("chr1:1-1000 {flag} f.txt --against x"))).unwrap_err();
            assert!(
                matches!(
                    error,
                    ArgError::WrongTrack {
                        flag: "--against",
                        ..
                    }
                ),
                "{flag} accepted a second path: {error}"
            );
        }

        // A height flag that parses and goes nowhere is a figure that is not
        // the one asked for and does not look wrong, so the two that have no
        // height of their own refuse it by name.
        for flag in ["--methylation", "--structural"] {
            let it = draw(&format!("chr1:1-1000 {flag} f.txt --height 90"));
            assert_eq!(it.tracks[0].height, Some(90.0), "{flag}");
        }
        let error = parse(&args("chr1:1-1000 --split-reads r.sam --height 90")).unwrap_err();
        assert!(
            matches!(
                error,
                ArgError::WrongTrack {
                    flag: "--height",
                    ..
                }
            ),
            "{error}"
        );
    }

    /// The span is what a figure costs in memory, eight bytes a base, and the
    /// buffer is allocated in one piece. On a 32-bit target a `Vec` may not
    /// exceed `isize::MAX` bytes, so one base past this is a capacity overflow,
    /// which is a panic, which in the build the documentation site runs is a
    /// trap with no `Result` to carry it. The old limit was thirty-two
    /// gigabytes and refused nothing that mattered.
    #[test]
    fn the_longest_region_is_one_a_thirty_two_bit_target_can_hold() {
        assert_eq!(MAX_SPAN, 268_435_455);
        // Eight bytes a base has to stay inside what a Vec can address.
        assert!(MAX_SPAN * 8 <= i32::MAX as u64);

        // Every chromosome anyone sequences is inside it. The longest human
        // one is chromosome 1.
        assert!(draw("chr1:1-248,956,422 --coverage d.bg").tracks.len() == 1);

        let error = parse(&args("chr1:1-268,435,456 --coverage d.bg")).unwrap_err();
        assert!(matches!(error, ArgError::HugeRegion { .. }), "{error}");
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
        // and abort the process with `capacity overflow`. The limit that was
        // put in front of it stopped that and was still four gigabases, which
        // is sixteen times what a 32-bit target can hold, so the same abort
        // came back the moment the crate was compiled for one.
        let err = parse(&args("chr1:1-18446744073709551615 --coverage d.bg")).unwrap_err();
        assert_eq!(
            err.to_string(),
            "\"chr1:1-18446744073709551615\" spans 18446744073709551615 bases, \
             and a figure is drawn over at most 268435455"
        );
        // A whole large sequence is an ordinary figure and stays one.
        assert_eq!(draw("chr1:1-248956422").region.len(), 248_956_422);
    }

    #[test]
    fn a_missing_region_is_its_own_message() {
        let err = parse(&args("--coverage d.bg")).unwrap_err();
        assert!(err.to_string().contains("the first argument is the region"));
    }
    #[test]
    fn focus_takes_a_clade_a_tip_or_a_pair_and_refuses_the_rest() {
        let focus = |value: &str| {
            parse(&args(&format!("tree:1-1 --tree t.nwk --focus {value}"))).map(|request| {
                match request {
                    Request::Draw(invocation) => invocation.tracks[0].focus.clone(),
                    other => panic!("expected a figure, got {other:?}"),
                }
            })
        };
        assert_eq!(
            focus("outbreak").unwrap(),
            Some(vec!["outbreak".to_string()])
        );
        assert_eq!(
            focus("A,B").unwrap(),
            Some(vec!["A".to_string(), "B".to_string()])
        );
        // Spaces round a comma survive, for the reader who quotes them. It
        // is spelled out as one word here because a shell would have done
        // that, and a helper that splits on spaces would not.
        let spaced: Vec<String> = ["tree:1-1", "--tree", "t.nwk", "--focus", "A, B"]
            .iter()
            .map(|word| word.to_string())
            .collect();
        match parse(&spaced).unwrap() {
            Request::Draw(invocation) => assert_eq!(
                invocation.tracks[0].focus,
                Some(vec!["A".to_string(), "B".to_string()])
            ),
            other => panic!("expected a figure, got {other:?}"),
        }
        // Three names name no clade: two pick one out and a third can only
        // agree or contradict.
        let refused = focus("A,B,C").unwrap_err().to_string();
        assert!(refused.contains("--focus"), "{refused}");
        assert!(
            refused.contains("two tip names"),
            "the refusal says what it wanted: {refused}"
        );
        assert!(focus(",").is_err(), "a comma alone names nothing");
    }

    #[test]
    fn focus_means_nothing_to_a_track_that_is_not_a_tree() {
        let refused = parse(&args(
            "NC_000962.3:1-1000 --coverage depth.bedgraph --focus A",
        ))
        .unwrap_err()
        .to_string();
        assert!(refused.contains("--focus"), "{refused}");
        assert!(refused.contains("coverage"), "{refused}");
    }
}
