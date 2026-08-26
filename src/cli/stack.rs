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
    Aggregate, CoverageTrack, FeatureTrack, IdeogramTrack, LogoTrack, ManhattanTrack, MatrixTrack,
    MsaSequence, MsaTrack, OrfTrack, PileupTrack, Plot, Region, SequenceTrack, SnpTrack, Theme,
    Track, Tree, TreeTrack, VariantTrack, WindowStyle, WindowTrack,
};

use crate::cli::args::{Invocation, Kind, Palette, Source, TrackSpec};
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
    Tree {
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
            BuildError::Tree { path, cause } => write!(f, "--tree {path}: {cause}"),
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
fn slurp(
    spec: &TrackSpec,
    open: &mut dyn FnMut(&Source) -> io::Result<String>,
) -> Result<(String, String), BuildError> {
    let Some(source) = spec.source.as_ref() else {
        return Ok((String::new(), String::new()));
    };
    let name = match source {
        Source::Path(path) => path.display().to_string(),
        Source::Stdin => "standard input".to_string(),
    };
    let text = open(source).map_err(|cause| BuildError::Open {
        track: spec.kind.flag(),
        path: name.clone(),
        cause,
    })?;
    Ok((text, name))
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
                path: path.clone(),
                cause,
            })?;
            Box::new(named(TreeTrack::new(tree), label, TreeTrack::label))
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
            // No sample lines at all, and a header whose every site lies
            // outside the window, both leave nothing to draw: the second one
            // used to give a lane of names beside no cells.
            if rows.is_empty() || sites.is_empty() {
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
