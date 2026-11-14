//! The public types as a crate that depends on this one sees them.
//!
//! Most of what is checked here is checked by the compiler: a type that stops
//! being `Send` or `Sync` fails this file's build, and nothing has to run.
//! That is the promise a caller leans on when a figure is built on a worker
//! thread and drawn on another, or kept by a server, or when an error is sent
//! up a channel or boxed with `?`.

use std::error::Error as StdError;

use karyon::cli::args::ArgError;
use karyon::cli::stack::{BuildError, CircleRefusal, CodonRefusal, ShadeRefusal};
use karyon::read::ReadError;
use karyon::*;

/// Compiles only for a type that can be sent to another thread, shared
/// between threads, and holds no borrow.
fn send_sync<T: Send + Sync + 'static>() {}

/// A figure, a plot and a circle hold their tracks and rings as trait
/// objects, which are `Send + Sync` because [`Track`] and [`Ring`] ask for
/// it. Without that, these three are what stops compiling first.
#[test]
fn a_figure_a_plot_and_a_circle_can_cross_threads() {
    send_sync::<Figure>();
    send_sync::<Plot>();
    send_sync::<Plot<CoverageTrack>>();
    send_sync::<Rings>();
    send_sync::<Box<dyn Track>>();
    send_sync::<Box<dyn Ring>>();

    // And one is drawn on a thread it was not built on.
    let figure = plot("chr1:1-100")
        .unwrap()
        .add_coverage(vec![9.0; 100])
        .into_figure();
    let svg = std::thread::spawn(move || figure.to_svg()).join().unwrap();
    assert!(svg.starts_with("<svg"));
}

#[test]
fn every_track_the_crate_ships_can_cross_threads() {
    send_sync::<AxisTrack>();
    send_sync::<BisulfiteTrack>();
    send_sync::<CladeTrack>();
    send_sync::<CodonTrack>();
    send_sync::<CopyNumberTrack>();
    send_sync::<CoverageTrack>();
    send_sync::<DomainTrack>();
    send_sync::<DotplotTrack>();
    send_sync::<DynseqTrack>();
    send_sync::<FeatureTrack>();
    send_sync::<GenomeTrack>();
    send_sync::<GenotypeTrack>();
    send_sync::<IdeogramTrack>();
    send_sync::<JunctionTrack>();
    send_sync::<LegendTrack>();
    send_sync::<LocusTrack>();
    send_sync::<LogoTrack>();
    send_sync::<ManhattanTrack>();
    send_sync::<MatrixTrack>();
    send_sync::<MethylationTrack>();
    send_sync::<MsaTrack>();
    send_sync::<OrfTrack>();
    send_sync::<PairTrack>();
    send_sync::<PhylodynamicTrack>();
    send_sync::<PileupTrack>();
    send_sync::<SelectionTrack>();
    send_sync::<SequenceTrack>();
    send_sync::<SnpTrack>();
    send_sync::<SplitReadTrack>();
    send_sync::<SquiggleTrack>();
    send_sync::<StructuralTrack>();
    send_sync::<SurveillanceTrack>();
    send_sync::<SyntenyTrack>();
    send_sync::<TanglegramTrack>();
    send_sync::<TranscriptionUnitTrack>();
    send_sync::<TreeTrack>();
    send_sync::<VariantTrack>();
    send_sync::<WindowTrack>();
}

/// The list above is a copy of the one `plot.rs` implements `Slot` for, which
/// is every track the crate ships, and a copy is the kind of thing that falls
/// behind. This holds the two together: a track added there and not here
/// fails here, rather than going unchecked.
#[test]
fn the_list_of_tracks_above_is_the_plot_s_list() {
    const PLOT: &str = include_str!("../src/plot.rs");
    const HERE: &str = include_str!("api.rs");
    let list = PLOT
        .split_once("\ntracks![\n")
        .and_then(|(_, rest)| rest.split_once("\n];"))
        .expect("plot.rs lists its tracks in tracks![...]")
        .0;
    let tracks: Vec<&str> = list
        .split(',')
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .collect();
    assert!(!tracks.is_empty(), "no track found in tracks![...]");
    for track in tracks {
        assert!(
            HERE.contains(&format!("send_sync::<{track}>();")),
            "{track} is a track plot.rs ships and this file does not check"
        );
    }
}

#[test]
fn the_rings_the_other_drawings_and_the_values_they_are_built_from_can_cross_threads() {
    send_sync::<AxisRing>();
    send_sync::<FeatureRing>();
    send_sync::<MarkerRing>();
    send_sync::<SignalRing>();
    send_sync::<Panels>();
    send_sync::<Map>();
    send_sync::<PhyloMap>();
    send_sync::<Pdf>();
    send_sync::<Region>();
    send_sync::<Genome>();
    send_sync::<Theme>();
    send_sync::<Tree>();
    send_sync::<Sheet>();
    send_sync::<Shade>();
    send_sync::<Legend>();
    send_sync::<Feature>();
    send_sync::<Variant>();
    send_sync::<Read>();
}

#[test]
fn every_error_can_cross_threads_and_goes_into_a_boxed_error() {
    send_sync::<Error>();
    send_sync::<ReadError>();
    send_sync::<ArgError>();
    send_sync::<BuildError>();
    send_sync::<CircleRefusal>();
    send_sync::<CodonRefusal>();
    send_sync::<ShadeRefusal>();

    // What `?` does with each, into the box a program's `main` returns.
    fn boxed<E: StdError + Send + Sync + 'static>(error: E) -> Box<dyn StdError + Send + Sync> {
        error.into()
    }
    let words = |line: &str| {
        line.split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>()
    };
    let errors = [
        boxed(Region::parse("chr1").unwrap_err()),
        boxed(read::seq::fasta("ACGT\n").unwrap_err()),
        boxed(cli::args::parse(&words("chr1:1-100 --no-such-flag")).unwrap_err()),
        boxed(
            cli::stack::build(
                match cli::args::parse(&words("chr1:1-100 --coverage depth.txt")).unwrap() {
                    cli::args::Request::Draw(invocation) => invocation,
                    other => panic!("that line draws, and it was read as {other:?}"),
                }
                .as_ref(),
                |_| Ok("chr1\tten\t20\t7\n".to_string()),
            )
            .unwrap_err(),
        ),
    ];
    for error in errors {
        assert!(!error.to_string().is_empty(), "{error:?}");
    }

    // The reasons a build gives inside a BuildError are errors of their own.
    let reasons = [
        boxed(CircleRefusal::NotAFigure),
        boxed(CodonRefusal::NoAnnotation),
        boxed(ShadeRefusal::NoSuchGene { near: Vec::new() }),
    ];
    for reason in reasons {
        assert!(!reason.to_string().is_empty(), "{reason:?}");
    }
}
