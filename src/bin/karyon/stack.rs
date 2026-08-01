//! Turning a parsed command line into a figure.
//!
//! Every arm here is one track of the stack, built in the order the flags were
//! written, which is the whole reason the grammar in [`crate::args`] looks the
//! way it does.
//!
//! The tracks are built and then handed over with
//! [`Plot::add_track`](karyon::Plot::add_track) rather than through the `add_`
//! methods, because a command line always has the settings in hand before the
//! track exists: `--label` and `--height` have already been read by the time
//! this runs, so there is nothing left for `Plot::label` to do afterwards.

use std::fmt;
use std::fs;
use std::io::{self, Read as _};

use karyon::{
    Aggregate, CoverageTrack, FeatureTrack, IdeogramTrack, ManhattanTrack, MatrixTrack,
    MsaSequence, MsaTrack, PileupTrack, Plot, Region, SequenceTrack, SnpTrack, Theme, Track, Tree,
    TreeTrack, VariantTrack, WindowStyle, WindowTrack,
};

use crate::args::{Invocation, Kind, Palette, Source, TrackSpec};
use crate::read;

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
    /// The tree would not parse.
    Tree(karyon::Error),
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
            BuildError::Tree(cause) => write!(f, "--tree: {cause}"),
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
/// by [`crate::args::parse`], so everything here is about the data.
pub fn build(invocation: &Invocation) -> Result<String, BuildError> {
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
        plot = plot.add_boxed(track(spec, region)?);
    }
    Ok(plot.to_svg())
}

/// Reads a track's source, whether that is a file or standard input.
fn slurp(spec: &TrackSpec) -> Result<(String, String), BuildError> {
    let track = spec.kind.flag();
    match spec.source.as_ref() {
        Some(Source::Path(path)) => {
            let name = path.display().to_string();
            let text = fs::read_to_string(path).map_err(|cause| BuildError::Open {
                track,
                path: name.clone(),
                cause,
            })?;
            Ok((text, name))
        }
        Some(Source::Stdin) => {
            let mut text = String::new();
            io::stdin()
                .read_to_string(&mut text)
                .map_err(|cause| BuildError::Open {
                    track,
                    path: "standard input".to_string(),
                    cause,
                })?;
            Ok((text, "standard input".to_string()))
        }
        None => Ok((String::new(), String::new())),
    }
}

/// Builds one track from its flags and its file.
fn track(spec: &TrackSpec, region: &Region) -> Result<Box<dyn Track>, BuildError> {
    let (text, path) = slurp(spec)?;
    let name = spec.kind.flag();
    let label = spec.label.clone();
    let height = spec.height;

    let empty = |wanted: &'static str| BuildError::Empty {
        track: name,
        path: path.clone(),
        wanted,
    };

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
            let tree = Tree::parse_newick(text.trim()).map_err(BuildError::Tree)?;
            Box::new(named(TreeTrack::new(tree), label, TreeTrack::label))
        }
        Kind::Msa => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            Box::new(named(MsaTrack::new(sequences), label, MsaTrack::label))
        }
        Kind::Snps => {
            let sequences = msa(wrap(name, &path, read::seq::alignment(&text))?, &empty)?;
            let track = SnpTrack::from_alignment(0, &sequences);
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
        Kind::Matrix => {
            let (sites, rows) = wrap(name, &path, read::table::matrix(&text, region))?;
            if rows.is_empty() {
                return Err(empty("samples"));
            }
            Box::new(named(
                MatrixTrack::new(sites, rows),
                label,
                MatrixTrack::label,
            ))
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
    use crate::args::{parse, Request};

    fn invocation(line: &str) -> Invocation {
        let args: Vec<String> = line.split_whitespace().map(String::from).collect();
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
        let svg = build(&invocation("chr1:1-1000")).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.ends_with("</svg>"));
    }

    #[test]
    fn figure_flags_reach_the_figure() {
        let svg = build(&invocation("chr1:1-1000 --title one --width 1200")).unwrap();
        assert!(svg.contains("width=\"1200\""));
        assert!(svg.contains(">one<"));
    }

    #[test]
    fn a_file_that_is_not_there_says_which_flag_wanted_it() {
        let error = build(&invocation("chr1:1-1000 --features nowhere.bed")).unwrap_err();
        let message = error.to_string();
        assert!(message.starts_with("--features nowhere.bed:"), "{message}");
    }
}
