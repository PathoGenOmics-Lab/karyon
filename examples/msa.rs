//! Renders the alignment figures used in the README.
//!
//! ```text
//! cargo run --example msa -- assets
//! ```
//!
//! The coordinates here are alignment columns, not genomic positions, so the
//! region spans the width of the alignment and the ruler counts columns.
//!
//! The conservation logo above the alignment is built from the same sequences
//! the alignment track draws, which is the point of both taking a list of rows.

use std::env;
use std::path::PathBuf;

use karyon::{plot, LogoScore, LogoTrack, MsaColoring, MsaDisplay, MsaSequence};

/// Width of the nucleotide alignment in columns.
const COLUMNS: usize = 120;

fn main() -> std::io::Result<()> {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let rows = nucleotide_alignment();
    let plain: Vec<String> = rows
        .iter()
        .map(|row| String::from_utf8_lossy(&row.residues).into_owned())
        .collect();

    let figure = plot(&format!("alignment:1-{COLUMNS}"))?
        .title("An alignment, and what disagrees in it")
        .width(940.0)
        .remove_region_label()
        // The same sequences, scored for conservation. Shrinkage is on
        // because twelve rows is not many, and a logo drawn from twelve
        // sequences should say so.
        .add_track(
            LogoTrack::from_sequences(0, &plain)
                .alphabet_size(4)
                .score(LogoScore::InformationContent)
                .stabilize()
                .label("conservation")
                .height(52.0),
        )
        .add_msa(rows.clone())
        .label("isolates")
        .adjust(|track| track.compare_to(0).row_height(13.0))
        .save(out.join("example-msa.svg"))?
        .into_figure();

    // A protein alignment, coloured by physicochemical class and zoomed far
    // enough in for the residues themselves.
    let protein = vec![
        MsaSequence::new("KatG_H37Rv", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_CDC1551", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_S315T", b"MPEQHPPITETTTGATSNGCPV".to_vec()),
        MsaSequence::new("KatG_Beijing", b"MPEQHPPVTETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_Erdman", b"MPEQ-PPITETTTGAASNGCPV".to_vec()),
    ];
    let residues = plot("alignment:1-22")?
        .title("The same idea on a protein, coloured by residue class")
        .width(720.0)
        .remove_region_label()
        .add_msa(protein)
        .label("KatG")
        .adjust(|track| {
            track
                .display(MsaDisplay::Bases)
                .coloring(MsaColoring::Residue)
                .compare_to(0)
                .row_height(15.0)
        })
        .add_axis()
        .adjust(|axis| axis.center_on_bases(true))
        .save(out.join("example-msa-protein.svg"))?
        .into_figure();

    let (w1, h1) = figure.dimensions();
    let (w2, h2) = residues.dimensions();
    println!("example-msa.svg         {w1:.0} x {h1:.0}");
    println!("example-msa-protein.svg {w2:.0} x {h2:.0}");
    Ok(())
}

/// Twelve isolates that agree about most of a locus.
fn nucleotide_alignment() -> Vec<MsaSequence> {
    let mut rng = Lcg::new(90_210);
    let backbone: Vec<u8> = (0..COLUMNS)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();

    (0..12)
        .map(|isolate| {
            let mut row = backbone.clone();
            // A clade of four shares two substitutions and a deletion.
            if isolate % 3 == 1 {
                row[31] = b'T';
                row[32] = b'G';
                for cell in row.iter_mut().take(76).skip(70) {
                    *cell = b'-';
                }
            }
            // Everyone picks up a little private variation.
            for _ in 0..2 {
                let column = (rng.next() as usize) % COLUMNS;
                row[column] = b"ACGT"[(rng.next() % 4) as usize];
            }
            MsaSequence::new(format!("isolate_{:02}", isolate + 1), row)
        })
        .collect()
}

/// A linear congruential generator, so the figure is reproducible without a
/// dependency.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
}
