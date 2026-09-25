//! The figures `cargo run --example msa` writes, as functions.
//!
//! The example writes them to files, and the documentation site's playground
//! includes this file too and draws them in the page, in the page's colours,
//! at the width of the column and over whatever columns the reader has moved
//! to. One copy of each figure serves both, so the two cannot drift apart.

use karyon::{
    Drawing, LogoScore, LogoTrack, MsaColoring, MsaDisplay, MsaSequence, Plot, Region, Theme,
};

/// Width of the nucleotide alignment in columns.
const COLUMNS: usize = 120;

/// `example-msa.svg`: an alignment, and what disagrees in it.
///
/// `theme`, `width` and `region` replace the light theme, the 940 pixels and
/// the whole alignment the committed figure is drawn with. The region counts
/// alignment columns rather than bases.
pub fn example_msa(theme: &Theme, width: Option<f64>, region: Option<&Region>) -> Box<dyn Drawing> {
    let rows = nucleotide_alignment();
    let plain: Vec<String> = rows
        .iter()
        .map(|row| String::from_utf8_lossy(&row.residues).into_owned())
        .collect();

    let own = Region::parse(&format!("alignment:1-{COLUMNS}"))
        .expect("the alignment in this example has columns");
    let figure = Plot::over(region.cloned().unwrap_or(own))
        .title("An alignment, and what disagrees in it")
        .theme(theme.clone())
        .width(width.unwrap_or(940.0))
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
        .add_msa(rows)
        .label("isolates")
        .adjust(|track| track.compare_to(0).row_height(13.0))
        .into_figure();
    Box::new(figure)
}

/// `example-msa-protein.svg`: the same idea on a protein, coloured by residue
/// class.
///
/// `theme`, `width` and `region` replace the light theme, the 720 pixels and
/// the whole alignment the committed figure is drawn with.
pub fn example_msa_protein(
    theme: &Theme,
    width: Option<f64>,
    region: Option<&Region>,
) -> Box<dyn Drawing> {
    // A protein alignment, coloured by physicochemical class and zoomed far
    // enough in for the residues themselves.
    let protein = vec![
        MsaSequence::new("KatG_H37Rv", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_CDC1551", b"MPEQHPPITETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_S315T", b"MPEQHPPITETTTGATSNGCPV".to_vec()),
        MsaSequence::new("KatG_Beijing", b"MPEQHPPVTETTTGAASNGCPV".to_vec()),
        MsaSequence::new("KatG_Erdman", b"MPEQ-PPITETTTGAASNGCPV".to_vec()),
    ];
    let own = Region::parse("alignment:1-22").expect("the alignment in this example has columns");
    let residues = Plot::over(region.cloned().unwrap_or(own))
        .title("The same idea on a protein, coloured by residue class")
        .theme(theme.clone())
        .width(width.unwrap_or(720.0))
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
        .into_figure();
    Box::new(residues)
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
