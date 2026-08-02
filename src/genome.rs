//! Several sequences on one axis.
//!
//! A [`Figure`](crate::Figure) is one [`Region`] on one sequence, which is
//! right for a locus and wrong for a genome: an assembly is two hundred
//! contigs, a eukaryote is a couple of dozen chromosomes, and a bacterium with
//! a plasmid is two replicons. A [`Genome`] lays the sequences end to end and
//! hands back the single region that covers them all, so a genome-wide picture
//! is an ordinary figure over an unusually long axis.
//!
//! # Every track works unchanged
//!
//! Positions cross onto the shared axis through [`Genome::at`] on the way in,
//! and nothing after that knows the difference: the coverage track, the feature
//! track, the variant track and the association track are all drawing on one
//! axis that happens to be several sequences long. The classic genome-wide
//! association figure is a [`ManhattanTrack`](crate::ManhattanTrack) over a
//! `Genome`, with [`Genome::boundaries`] handed to
//! [`ManhattanTrack::bands`](crate::ManhattanTrack::bands).
//!
//! # The name is the join
//!
//! Everything crosses onto the shared axis through a sequence name, and nothing
//! here normalises one, so a reference that calls it `chr1` and a variant file
//! that calls it `1` never meet. Stripping a prefix would be right for one pair
//! of files and wrong for the next, so [`Genome::at`] and [`Genome::map`] both
//! refuse to guess, and both say so.
//!
//! # What it costs
//!
//! The axis is a concatenation, so a distance measured across a boundary is not
//! a distance. Two points a pixel apart may be the last base of one contig and
//! the first of the next, which are not neighbours in any sense that matters.
//! [`GenomeTrack`](crate::GenomeTrack) draws where the joins are, and a ruler
//! of global coordinates would be a ruler of a coordinate system nothing else
//! uses, so it labels each sequence instead.

use crate::error::Error;
use crate::region::Region;

/// One sequence of a genome: a chromosome, a plasmid, a contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chromosome {
    /// Name as it appears in the FASTA header or the BAM.
    pub name: String,
    /// Length in bases.
    pub length: u64,
}

impl Chromosome {
    /// A sequence `length` bases long.
    pub fn new(name: impl Into<String>, length: u64) -> Self {
        Chromosome {
            name: name.into(),
            length: length.max(1),
        }
    }
}

/// Several sequences laid end to end on one coordinate axis.
///
/// ```
/// use karyon::{Genome, GenomeTrack, Figure, ManhattanTrack, Association};
///
/// let genome = Genome::new([("chr1", 1_000_000u64), ("chr2", 600_000), ("chr3", 400_000)]);
/// assert_eq!(genome.total(), 2_000_000);
///
/// // A hit on chr2 becomes a position on the shared axis.
/// let at = genome.at("chr2", 1_000).unwrap();
/// assert_eq!(at, 1_001_000);
/// assert_eq!(genome.locate(at), Some(("chr2", 1_000)));
///
/// let svg = Figure::new(genome.region())
///     .push(ManhattanTrack::new(vec![Association::new(at, 8.1)]).bands(genome.boundaries()))
///     .push(GenomeTrack::new(genome))
///     .to_svg();
/// assert!(svg.contains("chr2"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genome {
    sequences: Vec<Chromosome>,
    gap: u64,
}

impl Genome {
    /// A genome of sequences in the order given.
    ///
    /// The order is the drawing order and is not sorted: a reference is
    /// conventionally ordered by chromosome number and an assembly by
    /// descending contig length, and neither is something to guess at.
    pub fn new(sequences: impl IntoIterator<Item = (impl Into<String>, u64)>) -> Self {
        Genome {
            sequences: sequences
                .into_iter()
                .map(|(name, length)| Chromosome::new(name, length))
                .collect(),
            gap: 0,
        }
    }

    /// Puts `bases` of blank axis between one sequence and the next.
    ///
    /// Zero by default. A gap makes the joins obvious without any drawing, at
    /// the cost of coordinates that belong to no sequence:
    /// [`Genome::locate`] returns `None` inside one, which is the honest
    /// answer.
    pub fn gap(mut self, bases: u64) -> Self {
        self.gap = bases;
        self
    }

    /// The sequences.
    pub fn sequences(&self) -> &[Chromosome] {
        &self.sequences
    }

    /// How many sequences the genome has.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Whether the genome has no sequences.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Length of the whole axis, gaps included.
    pub fn total(&self) -> u64 {
        let bases: u64 = self.sequences.iter().map(|seq| seq.length).sum();
        bases + self.gap * self.sequences.len().saturating_sub(1) as u64
    }

    /// The region covering every sequence, for the figure to be built over.
    pub fn region(&self) -> Region {
        Region::new("genome", 0, self.total().max(1)).expect("the total is at least one")
    }

    /// Where a sequence starts on the shared axis.
    pub fn offset(&self, name: &str) -> Option<u64> {
        let mut at = 0u64;
        for sequence in &self.sequences {
            if sequence.name == name {
                return Some(at);
            }
            at += sequence.length + self.gap;
        }
        None
    }

    /// A position on one sequence, as a position on the shared axis.
    ///
    /// `None` when the genome has no such sequence, which is what a name
    /// mismatch between a VCF and a reference looks like, and `None` when the
    /// position is past that sequence's own end, which is what a file called
    /// against a different reference build, or 1-based positions handed in
    /// unconverted, look like. Better to notice either than to draw the point
    /// somewhere plausible, which here means inside the next sequence along.
    ///
    /// A caller that wants the shared-axis coordinate of the end of a
    /// half-open interval, where `position` may equal the length, wants
    /// [`Genome::offset`] instead.
    pub fn at(&self, name: &str, position: u64) -> Option<u64> {
        let mut offset = 0u64;
        for sequence in &self.sequences {
            if sequence.name == name {
                return (position < sequence.length).then_some(offset + position);
            }
            offset += sequence.length + self.gap;
        }
        None
    }

    /// The sequence and offset a shared-axis position falls on.
    ///
    /// The inverse of [`Genome::at`]. `None` inside a gap, or past the end.
    pub fn locate(&self, position: u64) -> Option<(&str, u64)> {
        let mut at = 0u64;
        for sequence in &self.sequences {
            if position < at + sequence.length {
                return Some((&sequence.name, position - at));
            }
            at += sequence.length;
            if position < at + self.gap {
                return None;
            }
            at += self.gap;
        }
        None
    }

    /// Where each sequence starts on the shared axis.
    ///
    /// Hand these to [`ManhattanTrack::bands`](crate::ManhattanTrack::bands) to
    /// get the alternating chromosome colouring a genome-wide plot is read by.
    pub fn boundaries(&self) -> Vec<u64> {
        let mut at = 0u64;
        let mut starts = Vec::with_capacity(self.sequences.len());
        for sequence in &self.sequences {
            starts.push(at);
            at += sequence.length + self.gap;
        }
        starts
    }

    /// Each sequence as its name and its span on the shared axis.
    pub fn spans(&self) -> Vec<(&str, u64, u64)> {
        let mut at = 0u64;
        let mut out = Vec::with_capacity(self.sequences.len());
        for sequence in &self.sequences {
            out.push((sequence.name.as_str(), at, at + sequence.length));
            at += sequence.length + self.gap;
        }
        out
    }

    /// Maps a list of per-sequence positions onto the shared axis.
    ///
    /// Anything naming a sequence the genome does not have is dropped, as is
    /// anything sitting past the end of the sequence it names, and the count of
    /// those comes back with the rest so a caller can say how many went missing
    /// rather than wondering why the figure looks thin.
    pub fn map<T>(
        &self,
        items: impl IntoIterator<Item = (String, u64, T)>,
    ) -> (Vec<(u64, T)>, usize) {
        let mut mapped = Vec::new();
        let mut dropped = 0usize;
        for (name, position, value) in items {
            match self.at(&name, position) {
                Some(at) => mapped.push((at, value)),
                None => dropped += 1,
            }
        }
        (mapped, dropped)
    }

    /// Builds a genome, failing when a name is repeated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::DuplicateSequence`] when two sequences share a name,
    /// since [`Genome::at`] would then silently place everything on the first
    /// of them.
    pub fn checked(
        sequences: impl IntoIterator<Item = (impl Into<String>, u64)>,
    ) -> Result<Self, Error> {
        let genome = Genome::new(sequences);
        for (index, sequence) in genome.sequences.iter().enumerate() {
            if genome.sequences[..index]
                .iter()
                .any(|earlier| earlier.name == sequence.name)
            {
                return Err(Error::DuplicateSequence {
                    name: sequence.name.clone(),
                });
            }
        }
        Ok(genome)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome() -> Genome {
        Genome::new([("chr1", 1_000u64), ("chr2", 600), ("chr3", 400)])
    }

    #[test]
    fn sequences_lie_end_to_end_on_one_axis() {
        let genome = genome();
        assert_eq!(genome.total(), 2_000);
        assert_eq!(genome.region().len(), 2_000);
        assert_eq!(genome.offset("chr1"), Some(0));
        assert_eq!(genome.offset("chr2"), Some(1_000));
        assert_eq!(genome.offset("chr3"), Some(1_600));
        assert_eq!(genome.offset("chrX"), None);
    }

    #[test]
    fn a_position_maps_there_and_back() {
        let genome = genome();
        let at = genome.at("chr2", 250).unwrap();
        assert_eq!(at, 1_250);
        assert_eq!(genome.locate(at), Some(("chr2", 250)));
        assert_eq!(genome.locate(0), Some(("chr1", 0)));
        assert_eq!(genome.locate(1_999), Some(("chr3", 399)));
        assert_eq!(genome.locate(2_000), None, "past the end");
    }

    #[test]
    fn a_name_the_genome_does_not_have_comes_back_empty() {
        // A name mismatch between a VCF and a reference is common and quiet.
        // Better to notice it than to draw the point somewhere plausible.
        assert_eq!(genome().at("chr9", 10), None);
    }

    #[test]
    fn a_position_past_the_end_of_its_sequence_is_refused_rather_than_landing_on_the_next() {
        // chr1 is 1,000 bases, so 1,500 is 500 bases into chr2. Placing it
        // there would draw a chr1 point inside chr2's band and count it as
        // mapped, which is what a file called against another reference build,
        // or 1-based positions handed in unconverted, look like.
        let genome = genome();
        assert_eq!(genome.spans()[1], ("chr2", 1_000, 1_600));
        assert_eq!(genome.at("chr1", 1_500), None);
        assert_eq!(genome.at("chr1", 1_000), None, "the length is past the end");
        assert_eq!(genome.at("chr1", 999), Some(999), "the last base is not");
        assert_eq!(
            genome.map([("chr1".to_string(), 1_500u64, 9.9)]),
            (Vec::new(), 1),
            "and it is counted, not quietly placed"
        );
    }

    #[test]
    fn every_position_at_accepts_comes_back_from_locate() {
        // locate is documented as the inverse of at, so the two have to agree
        // on which positions exist at all.
        for genome in [genome(), genome().gap(100)] {
            for sequence in genome.sequences() {
                for position in [0, 1, sequence.length / 2, sequence.length - 1] {
                    let at = genome
                        .at(&sequence.name, position)
                        .expect("in range on its own sequence");
                    assert_eq!(
                        genome.locate(at),
                        Some((sequence.name.as_str(), position)),
                        "{} at {position}",
                        sequence.name
                    );
                }
                assert_eq!(genome.at(&sequence.name, sequence.length), None);
            }
        }
    }

    #[test]
    fn a_gap_belongs_to_no_sequence() {
        let spaced = genome().gap(100);
        assert_eq!(spaced.total(), 2_000 + 200);
        assert_eq!(spaced.offset("chr2"), Some(1_100));
        assert_eq!(spaced.locate(1_050), None, "inside the gap");
        assert_eq!(spaced.locate(1_100), Some(("chr2", 0)));
    }

    #[test]
    fn boundaries_are_where_each_sequence_starts() {
        assert_eq!(genome().boundaries(), vec![0, 1_000, 1_600]);
        assert_eq!(genome().gap(50).boundaries(), vec![0, 1_050, 1_700]);
    }

    #[test]
    fn spans_carry_the_name_and_both_ends() {
        assert_eq!(
            genome().spans(),
            vec![
                ("chr1", 0, 1_000),
                ("chr2", 1_000, 1_600),
                ("chr3", 1_600, 2_000)
            ]
        );
    }

    #[test]
    fn mapping_a_list_reports_what_it_could_not_place() {
        let (mapped, dropped) = genome().map([
            ("chr1".to_string(), 10u64, 1.0),
            ("chr3".to_string(), 20, 2.0),
            ("chrZ".to_string(), 30, 3.0),
        ]);
        assert_eq!(mapped, vec![(10, 1.0), (1_620, 2.0)]);
        assert_eq!(dropped, 1, "a silent drop is worse than a thin figure");
    }

    #[test]
    fn a_repeated_name_is_refused_rather_than_silently_merged() {
        assert!(Genome::checked([("chr1", 10u64), ("chr1", 20)]).is_err());
        assert!(Genome::checked([("chr1", 10u64), ("chr2", 20)]).is_ok());
        // The unchecked constructor still builds one, and everything lands on
        // the first of the two.
        let sloppy = Genome::new([("chr1", 10u64), ("chr1", 20)]);
        assert_eq!(sloppy.at("chr1", 0), Some(0));
    }

    #[test]
    fn an_empty_genome_still_gives_a_region_to_draw_on() {
        let empty = Genome::new(Vec::<(String, u64)>::new());
        assert!(empty.is_empty());
        assert_eq!(empty.total(), 0);
        assert_eq!(empty.region().len(), 1);
        assert!(empty.boundaries().is_empty());
    }

    #[test]
    fn a_zero_length_sequence_is_widened_rather_than_vanishing() {
        let genome = Genome::new([("empty", 0u64), ("chr", 100)]);
        assert_eq!(genome.sequences()[0].length, 1);
        assert_eq!(genome.offset("chr"), Some(1));
    }
}
