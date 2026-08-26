//! Molecules that aligned in pieces, from SAM and its `SA` tag.
//!
//! One read that visits three places in the reference is one primary alignment
//! and two supplementary ones, and the `SA` tag on each of them lists the
//! others:
//!
//! ```text
//! SA:Z:rname,pos,strand,CIGAR,mapQ,NM;...
//! ```
//!
//! `pos` inside that tag is 1-based, like column four and unlike everything
//! else this reader touches, so both are moved back one and an audit fixture
//! pins the two on the same base. Reading one and not the other puts every
//! supplementary alignment exactly one base from where the aligner put it,
//! which at a whole-chromosome scale is invisible and at a breakpoint is the
//! thing being looked at.
//!
//! # The order of the pieces is the whole claim
//!
//! [`SplitRead`] takes its segments in read order and the figure draws them in
//! that order: the connectors say in what order and in which direction one
//! molecule visited those coordinates, and
//! [`SplitRead::goes_backwards`](crate::SplitRead::goes_backwards) is what
//! separates a read crossing an inversion from a read crossing a deletion.
//!
//! Handed the same segments sorted by reference position instead, the same read
//! draws a different rearrangement, and both figures look finished. So the order
//! is computed rather than assumed, from where each alignment sat on the read.
//!
//! That is not the leading clip of its CIGAR. A CIGAR is written along the
//! reference, so a reverse-strand alignment's clips are the far end of the
//! molecule, and its position on the read is measured from the other side:
//! `length - end` to `length - start`. Skipping that step reverses the order of
//! every read that crosses an inversion, which is exactly the read an inversion
//! figure is made of.
//!
//! # Each piece once
//!
//! A supplementary alignment is a line of its own *and* an entry in the primary's
//! `SA`, so reading every line and every tag counts each piece twice and draws a
//! hop from a place to itself. Only primary lines are read here, and their
//! segments come from the line and its own tag, which also recovers the pieces a
//! region-restricted `samtools view` never wrote out.
//!
//! The cost is that a read whose primary alignment is not in the text is not
//! read at all. `samtools view -f 0x800` alone therefore draws nothing, and that
//! is better than drawing each of those reads with the one piece that happened
//! to be in view.

use std::collections::BTreeMap;

use crate::{CigarOp, Region, SplitRead, SplitSegment, Strand};

use super::align::cigar;
use super::{columns, lines, number, ReadError};

/// The reads a SAM holds, and every read it declined to draw.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitReads {
    /// One read per molecule, its segments in read order.
    pub reads: Vec<SplitRead>,
    /// Primary alignments in the file, before any filter.
    pub records: usize,
    /// Reads that aligned in one piece, which is not a split read.
    pub not_split: usize,
    /// Reads with a piece on another sequence.
    ///
    /// Dropped whole rather than drawn with the pieces that are here. A read
    /// that went from this sequence to another and back, drawn without the
    /// middle, is two pieces side by side with a connector between them, which
    /// is the picture of a deletion junction and not of a translocation.
    pub elsewhere: usize,
    /// Reads none of whose pieces touch the window.
    pub off_region: usize,
    /// Reads whose alignments disagree about how long the molecule is.
    ///
    /// Every alignment of one read spends the same number of bases on clips,
    /// matches and insertions, so a disagreement means the order cannot be
    /// worked out, and an order that is guessed at is a rearrangement that is
    /// guessed at.
    pub ambiguous_order: usize,
}

/// Reads the split alignments out of SAM text.
///
/// Column four and the `SA` tag's own position are both 1-based, so both come
/// back one lower. Segments are put in read order, which for a reverse-strand
/// alignment is measured from the far end of the molecule.
///
/// Only primary alignments are read, since a supplementary one is already an
/// entry in its primary's tag.
///
/// # Errors
///
/// Returns the first row that is not SAM, whose CIGAR will not parse, whose
/// strand is neither `+` nor `-`, or whose `SA` entry covers no reference at
/// all. A file that parses and holds no split read in the window is not an
/// error here: [`SplitReads`] says why, and the caller decides.
pub fn reads(text: &str, region: &Region) -> Result<SplitReads, ReadError> {
    let mut found = SplitReads {
        reads: Vec::new(),
        records: 0,
        not_split: 0,
        elsewhere: 0,
        off_region: 0,
        ambiguous_order: 0,
    };

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 11 {
            return Err(ReadError::at(
                at,
                format!(
                    "a SAM row has at least 11 columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        let flag: u16 = number(cols[1], "FLAG", at)?;
        // 0x100 secondary, 0x800 supplementary. Both are already an entry in
        // some primary's tag, and read as well they double every piece.
        if flag & 0x900 != 0 {
            continue;
        }
        if cols[2] == "*" {
            continue;
        }
        found.records += 1;

        let mut pieces = vec![piece(
            cols[2],
            number::<u64>(cols[3], "POS", at)?,
            if flag & 0x10 != 0 {
                Strand::Reverse
            } else {
                Strand::Forward
            },
            cols[5],
            number::<u8>(cols[4], "MAPQ", at)?,
            at,
        )?];

        for entry in tag(&cols, "SA:Z:").unwrap_or_default().split(';') {
            if entry.trim().is_empty() {
                continue;
            }
            let field: Vec<&str> = entry.split(',').collect();
            if field.len() < 5 {
                return Err(ReadError::at(
                    at,
                    format!(
                        "an SA entry is rname,pos,strand,CIGAR,mapQ,NM, and this one is {entry:?}"
                    ),
                ));
            }
            let strand = match field[2] {
                "+" => Strand::Forward,
                "-" => Strand::Reverse,
                other => {
                    return Err(ReadError::at(
                        at,
                        format!("the strand of an SA entry is + or -, this one is {other:?}"),
                    ))
                }
            };
            pieces.push(piece(
                field[0],
                number::<u64>(field[1], "the position of an SA entry", at)?,
                strand,
                field[3],
                number::<u8>(field[4], "the mapping quality of an SA entry", at)?,
                at,
            )?);
        }

        if pieces.len() < 2 {
            found.not_split += 1;
            continue;
        }
        if pieces.iter().any(|piece| piece.sequence != region.seq()) {
            found.elsewhere += 1;
            continue;
        }
        // Every alignment of one molecule accounts for every one of its bases,
        // so a disagreement means the read positions below are measured against
        // different rulers and the order that comes out of them is invented.
        let length = pieces[0].length;
        if pieces.iter().any(|piece| piece.length != length) {
            found.ambiguous_order += 1;
            continue;
        }
        if !pieces
            .iter()
            .any(|piece| piece.end > region.start() && piece.start < region.end())
        {
            found.off_region += 1;
            continue;
        }

        pieces.sort_by_key(|piece| (piece.read_from, piece.start));
        let segments = pieces
            .iter()
            .map(|piece| {
                let mut segment = SplitSegment::new(piece.start, piece.end, piece.strand)
                    .read_span(piece.read_from, piece.read_to);
                // 255 is the SAM spelling of "the aligner did not say". Passed
                // through it is a segment at full strength, drawn exactly like
                // one the aligner was certain of.
                if piece.mapq != 255 {
                    segment = segment.mapq(piece.mapq);
                }
                segment
            })
            .collect();

        let mut read = SplitRead::new(segments);
        if cols[0] != "*" {
            read = read.name(cols[0]);
        }
        found.reads.push(read);
    }

    Ok(found)
}

/// The read names in the file and how many primary alignments each has.
///
/// Ordered, so the same file answers the same way twice. For a caller wanting
/// to say how many molecules a figure of forty rows was chosen out of.
pub fn names(text: &str) -> Result<BTreeMap<String, usize>, ReadError> {
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for (_, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 11 {
            continue;
        }
        let Ok(flag) = cols[1].parse::<u16>() else {
            continue;
        };
        if flag & 0x900 == 0 && cols[2] != "*" {
            *seen.entry(cols[0].to_string()).or_insert(0) += 1;
        }
    }
    Ok(seen)
}

/// One alignment of one molecule, placed on both the reference and the read.
#[derive(Debug)]
struct Piece {
    sequence: String,
    start: u64,
    end: u64,
    strand: Strand,
    mapq: u8,
    /// Where this alignment sat on the molecule, counted from its own 5' end.
    read_from: u64,
    read_to: u64,
    /// How long the whole molecule is, by this alignment's reckoning.
    length: u64,
}

/// Places one alignment, on the reference and on the read.
fn piece(
    sequence: &str,
    pos: u64,
    strand: Strand,
    text: &str,
    mapq: u8,
    at: usize,
) -> Result<Piece, ReadError> {
    if pos == 0 {
        return Err(ReadError::at(
            at,
            "SAM counts from 1, so 0 is not a position",
        ));
    }
    let operations = cigar(text, at)?;

    let span: u64 = operations.iter().map(|op| op.reference_len()).sum();
    if span == 0 {
        return Err(ReadError::at(
            at,
            format!("the alignment {text:?} covers no reference bases"),
        ));
    }

    // Along the read as it is stored, which is the reference orientation. A
    // hard clip counts here: those bases were cut out of the stored sequence
    // and are still part of the molecule.
    let clipped = |op: &CigarOp| matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_));
    let consumed = |op: &CigarOp| match op {
        CigarOp::HardClip(n) => *n as u64,
        other => other.query_len(),
    };
    let length: u64 = operations.iter().map(consumed).sum();
    let leading: u64 = operations
        .iter()
        .take_while(|op| clipped(op))
        .map(consumed)
        .sum();
    let aligned: u64 = operations
        .iter()
        .filter(|op| !clipped(op))
        .map(|op| op.query_len())
        .sum();

    // A CIGAR runs along the reference, so a reverse-strand alignment's clips
    // are at the far end of the molecule and its place on the read is measured
    // from the other side. Without this the pieces of every read that crosses
    // an inversion come out in the opposite order, and the figure draws a
    // rearrangement that is the mirror of the one that happened.
    let (read_from, read_to) = if strand == Strand::Reverse {
        (
            length.saturating_sub(leading + aligned),
            length.saturating_sub(leading),
        )
    } else {
        (leading, leading + aligned)
    };

    Ok(Piece {
        sequence: sequence.to_string(),
        start: pos - 1,
        end: pos - 1 + span,
        strand,
        mapq,
        read_from,
        read_to,
        length,
    })
}

/// The value of an optional tag, from column twelve onwards.
fn tag<'a>(cols: &[&'a str], prefix: &str) -> Option<&'a str> {
    cols.iter()
        .skip(11)
        .find_map(|field| field.strip_prefix(prefix))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The specification's own chimeric read, both lines of it. The molecule
    /// aligns at 29 on the reverse strand and then at 9 on the forward one.
    const SAM: &str = "\
r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:ref,29,-,6H5M,17,0;
r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t*\tSA:Z:ref,9,+,5S6M,30,0;
";

    fn window() -> Region {
        Region::new("ref", 0, 100).unwrap()
    }

    #[test]
    fn both_positions_are_moved_back_one_because_both_count_from_one() {
        // Column four and the tag's own position are the same convention, and
        // reading one and not the other is a whole supplementary alignment
        // sitting one base from where the aligner put it.
        let found = reads(SAM, &window()).unwrap();
        assert_eq!(found.reads.len(), 1);
        let places: Vec<(u64, u64)> = found.reads[0]
            .segments()
            .iter()
            .map(|s| (s.start(), s.end()))
            .collect();
        // 29 with 5M is 28..33, and 9 with 6M is 8..14.
        assert!(places.contains(&(28, 33)), "{places:?}");
        assert!(places.contains(&(8, 14)), "{places:?}");
    }

    #[test]
    fn a_reverse_alignment_is_placed_on_the_read_from_its_far_end() {
        // This is the whole of the ordering. The reverse piece carries 6H5M, so
        // along the reference its clip leads; along the molecule those six
        // bases are the tail, and the piece is the first five bases of the read.
        let found = reads(SAM, &window()).unwrap();
        let segments = found.reads[0].segments();
        assert_eq!(
            (segments[0].start(), segments[0].end()),
            (28, 33),
            "the reverse piece is not first in read order"
        );
        assert_eq!((segments[1].start(), segments[1].end()), (8, 14));

        // And the figure says what that means: the molecule runs backwards
        // along the reference. Sorted by reference position instead it would
        // say the opposite, and look just as finished.
        assert!(found.reads[0].goes_backwards());
    }

    #[test]
    fn a_supplementary_line_is_not_read_twice() {
        // Both lines of the read are in the file and each names the other, so
        // reading every line and every tag gives four pieces and a hop from a
        // place to itself.
        let found = reads(SAM, &window()).unwrap();
        assert_eq!(found.records, 1, "a supplementary line was read as a read");
        assert_eq!(found.reads[0].segments().len(), 2);

        // And the primary alone reads the same, since its tag holds the rest.
        let primary: String = SAM.lines().next().unwrap().to_string();
        let alone = reads(&format!("{primary}\n"), &window()).unwrap();
        assert_eq!(alone.reads, found.reads);
    }

    #[test]
    fn a_read_that_aligned_in_one_piece_is_not_a_split_read() {
        let text = "r1\t0\tref\t9\t30\t11M\t*\t0\t0\tGCCTAAGCTAA\t*\n";
        let found = reads(text, &window()).unwrap();
        assert!(found.reads.is_empty());
        assert_eq!(found.not_split, 1);
        assert_eq!(found.records, 1);
    }

    #[test]
    fn a_read_with_a_piece_elsewhere_is_dropped_whole_and_counted() {
        // Drawn with only its local pieces it becomes two bars and a connector,
        // which is the picture of a deletion junction. The molecule crossed to
        // another sequence and back, which is a different event entirely.
        let text = "\
r1\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:other,29,-,6H5M,17,0;
";
        let found = reads(text, &window()).unwrap();
        assert!(found.reads.is_empty());
        assert_eq!(found.elsewhere, 1);

        // And the case a two-piece read cannot show: three pieces with the
        // middle one elsewhere. Keeping the two that are here leaves them side
        // by side with a connector between them, which is the picture of a
        // deletion junction drawn from a molecule that left the sequence.
        let three = "\
r2\t0\tref\t9\t60\t5S6M\t*\t0\t0\t*\t*\tSA:Z:other,100,+,5H4M2S,60,0;ref,40,+,9H2M,60,0;
";
        let found = reads(three, &window()).unwrap();
        assert!(
            found.reads.is_empty(),
            "a read that left the sequence was drawn without the part that left"
        );
        assert_eq!(found.elsewhere, 1);
    }

    /// 255 is the SAM spelling of "the aligner did not say", and it is not a
    /// quality of 255.
    ///
    /// What this pins is that the figure is the same either way. The track
    /// fades a segment below a quality of 20 and nothing above it, so today an
    /// unavailable quality and a high one already draw alike, and no assertion
    /// on the output can tell the guard from its absence. It is here so that a
    /// later fade curve, or a getter on the segment, cannot start reading 255
    /// as a measurement without this failing.
    #[test]
    fn a_mapping_quality_the_aligner_did_not_give_draws_as_one_that_is_absent() {
        let draw = |q: &str| {
            let text = format!(
                "r1\t0\tref\t9\t{q}\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:ref,29,-,6H5M,{q},0;\n"
            );
            let found = reads(&text, &window()).unwrap();
            crate::Figure::new(window())
                .push(crate::SplitReadTrack::new(found.reads))
                .to_svg()
        };
        assert_eq!(draw("255"), draw("60"), "255 is drawn as its own number");
        assert_ne!(draw("255"), draw("3"), "a low quality stopped fading");
    }

    #[test]
    fn an_sa_entry_that_covers_no_reference_stops_the_read() {
        // 150H is a clip and nothing else, so the alignment has no span, and
        // start == end draws a visible arrowhead labelled as one base.
        let text = "r1\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCC\t*\tSA:Z:ref,50,+,150H,60,0;\n";
        let error = reads(text, &window()).unwrap_err();
        assert!(error.reason.contains("covers no reference"), "{error}");
    }

    #[test]
    fn a_strand_that_is_neither_stops_the_read_rather_than_becoming_forward() {
        let text = "r1\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCC\t*\tSA:Z:ref,29,?,6H5M,17,0;\n";
        let error = reads(text, &window()).unwrap_err();
        assert!(error.reason.contains("+ or -"), "{error}");
    }

    #[test]
    fn alignments_that_disagree_about_the_molecule_cannot_be_ordered() {
        // 5S6M is eleven bases and 6H9M is fifteen, so one of the two is not
        // this read, and the read positions are measured against two rulers.
        let text = "r1\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCC\t*\tSA:Z:ref,29,-,6H9M,17,0;\n";
        let found = reads(text, &window()).unwrap();
        assert!(found.reads.is_empty());
        assert_eq!(found.ambiguous_order, 1);
    }

    #[test]
    fn a_read_with_nothing_in_the_window_is_counted_rather_than_given_a_row() {
        let found = reads(SAM, &Region::new("ref", 500, 600).unwrap()).unwrap();
        assert!(found.reads.is_empty());
        assert_eq!(found.off_region, 1);
        assert_eq!(found.records, 1, "the file did hold a read");
    }

    #[test]
    fn a_file_says_which_molecules_it_holds_before_anything_is_drawn() {
        assert_eq!(names(SAM).unwrap().get("r003"), Some(&1));
    }

    #[test]
    fn a_line_that_is_not_sam_stops_the_read_and_says_which() {
        let error = reads("r1\t0\tref\t9\n", &window()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("at least 11 columns"), "{error}");
    }
}
