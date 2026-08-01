//! SAM text, as `samtools view` writes it.

use karyon::{CigarOp, Read, Region, Strand};

use super::{columns, lines, number, ReadError};

/// The flag bit that says the record was never placed on a reference.
const UNMAPPED: u16 = 4;

/// The flag bit that says the read came off the reverse strand.
const REVERSE: u16 = 16;

/// How many columns a record has before the optional tags start.
const MANDATORY: usize = 11;

/// The MAPQ an aligner writes when it has no mapping quality to give.
const NO_QUALITY: u32 = 255;

/// Reads aligned records.
///
/// Column four is `POS`, 1-based, so the read starts at `POS - 1`. Column six
/// is the CIGAR, whose operations map to `karyon::CigarOp`: `M`, `=` and `X`
/// are all matches, since the track compares the sequences itself; `I`, `D`,
/// `N`, `S` and `H` are themselves; `P` is padding and is dropped.
///
/// The strand is bit 16 of the flag in column two. An unmapped record, bit 4,
/// is skipped, as is a record with `*` for a CIGAR.
///
/// Column ten is the read sequence, kept when it is not `*` so the track can
/// paint mismatches against the reference.
///
/// Records on another sequence than `region.seq()` are skipped.
pub fn sam(text: &str, region: &Region) -> Result<Vec<Read>, ReadError> {
    let mut reads = Vec::new();

    for (line, row) in lines(text) {
        let fields = columns(row);
        if fields.len() < MANDATORY {
            return Err(ReadError::at(
                line,
                format!(
                    "a SAM record has {MANDATORY} columns, this one has {}",
                    fields.len()
                ),
            ));
        }

        let flag: u16 = number(fields[1], "the flag", line)?;
        // A record with no place on the reference and a record with no shape
        // are both left out. Neither is malformed: they are what a SAM file
        // carries for reads the aligner could not put anywhere.
        if flag & UNMAPPED != 0 {
            continue;
        }
        if fields[2] != region.seq() {
            continue;
        }
        let shape = fields[5];
        if shape == "*" {
            continue;
        }

        let position: u64 = number(fields[3], "POS", line)?;
        if position == 0 {
            return Err(ReadError::at(
                line,
                "POS is 1-based, and 0 is a record that was never placed",
            ));
        }

        let strand = if flag & REVERSE != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        // The 1-based POS becomes the 0-based start the whole crate works in.
        let mut read = Read::new(position - 1, cigar(shape, line)?).strand(strand);

        let quality: u32 = number(fields[4], "MAPQ", line)?;
        if quality > NO_QUALITY {
            return Err(ReadError::at(
                line,
                format!("MAPQ runs from 0 to {NO_QUALITY}, and this one is {quality}"),
            ));
        }
        if quality != NO_QUALITY {
            read = read.mapping_quality(quality as u8);
        }

        let sequence = fields[9];
        if sequence != "*" && !sequence.is_empty() {
            read = read.sequence(sequence.as_bytes().to_vec());
        }

        // A read outside the window is not an error, it is simply not in this
        // figure, the same as a read on another sequence.
        if read.overlaps(region) {
            reads.push(read);
        }
    }

    Ok(reads)
}

/// Parses a CIGAR string: a run of digits, then one letter, repeated.
///
/// `P` is padding between reads in a multiple alignment of a pileup. It moves
/// along neither the reference nor the read, so it is dropped here rather than
/// carried as an operation that draws nothing.
fn cigar(text: &str, line: usize) -> Result<Vec<CigarOp>, ReadError> {
    if text.is_empty() {
        return Err(ReadError::at(line, "the CIGAR column is empty"));
    }

    let mut operations = Vec::new();
    let mut length: Option<u32> = None;

    for symbol in text.chars() {
        if let Some(digit) = symbol.to_digit(10) {
            let so_far = length.unwrap_or(0);
            length = Some(
                so_far
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(digit))
                    .ok_or_else(|| {
                        ReadError::at(
                            line,
                            format!("a CIGAR operation is longer than a read can be: {text:?}"),
                        )
                    })?,
            );
            continue;
        }

        let Some(length) = length.take() else {
            return Err(ReadError::at(
                line,
                format!("a CIGAR operation needs a length in front of it: {text:?}"),
            ));
        };
        operations.push(match symbol {
            // The track finds mismatches by comparing the sequences, so the
            // three letters that all mean aligned bases arrive as one.
            'M' | '=' | 'X' => CigarOp::Match(length),
            'I' => CigarOp::Insertion(length),
            'D' => CigarOp::Deletion(length),
            'N' => CigarOp::Skip(length),
            'S' => CigarOp::SoftClip(length),
            'H' => CigarOp::HardClip(length),
            'P' => continue,
            _ => {
                return Err(ReadError::at(
                    line,
                    format!("{symbol:?} is not a CIGAR operation, in {text:?}"),
                ))
            }
        });
    }

    if length.is_some() {
        return Err(ReadError::at(
            line,
            format!("a CIGAR ends with a length and no operation: {text:?}"),
        ));
    }
    Ok(operations)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record, written out column by column so a test can change one of them.
    fn record(
        name: &str,
        flag: u16,
        seq: &str,
        pos: u64,
        mapq: u32,
        cigar: &str,
        bases: &str,
    ) -> String {
        format!("{name}\t{flag}\t{seq}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{bases}\t*")
    }

    fn region() -> Region {
        // 1-based inclusive on the way in, so this is 3899 to 4100 inside.
        Region::parse("NC_002516.2:3,900-4,100").unwrap()
    }

    #[test]
    fn pos_is_one_based_and_the_read_starts_one_before_it() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "10M", "ACGTACGTAC");
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].start, 4000, "POS 4001 is 0-based 4000");
        assert_eq!(reads[0].end(), 4010);
    }

    #[test]
    fn the_first_base_of_a_sequence_is_pos_one_and_offset_zero() {
        let region = Region::new("chrIV", 0, 100).unwrap();
        let text = record("read1", 0, "chrIV", 1, 60, "5M", "ACGTA");
        let reads = sam(&text, &region).unwrap();
        assert_eq!(reads[0].start, 0);
    }

    #[test]
    fn the_cigar_is_walked_operation_by_operation() {
        let text = record(
            "read1",
            0,
            "NC_002516.2",
            4001,
            60,
            "3S5M2I4M1D6M2H",
            "AAAGGGGGTTCCCCTTTTTT",
        );
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(
            reads[0].cigar,
            vec![
                CigarOp::SoftClip(3),
                CigarOp::Match(5),
                CigarOp::Insertion(2),
                CigarOp::Match(4),
                CigarOp::Deletion(1),
                CigarOp::Match(6),
                CigarOp::HardClip(2),
            ]
        );
        // Only what advances along the reference counts towards the span.
        assert_eq!(reads[0].reference_span(), 16);
    }

    #[test]
    fn the_three_letters_for_aligned_bases_all_arrive_as_matches() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "4=1X5M", "ACGTACGTAC");
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(
            reads[0].cigar,
            vec![CigarOp::Match(4), CigarOp::Match(1), CigarOp::Match(5)]
        );
    }

    #[test]
    fn padding_is_dropped_and_does_not_shift_what_follows() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "5M3P5M", "ACGTACGTAC");
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(
            reads[0].cigar,
            vec![CigarOp::Match(5), CigarOp::Match(5)],
            "padding moves along neither the read nor the reference"
        );
        assert_eq!(reads[0].end(), 4010);
    }

    #[test]
    fn a_multi_digit_length_is_one_operation() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "151M", "");
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads[0].cigar, vec![CigarOp::Match(151)]);
    }

    #[test]
    fn bit_sixteen_is_the_reverse_strand() {
        let text = format!(
            "{}\n{}\n",
            record("forward", 0, "NC_002516.2", 4001, 60, "10M", "ACGTACGTAC"),
            record("reverse", 16, "NC_002516.2", 4021, 60, "10M", "ACGTACGTAC"),
        );
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads[0].strand, Strand::Forward);
        assert_eq!(reads[1].strand, Strand::Reverse);
    }

    #[test]
    fn an_unmapped_record_is_skipped_whatever_else_it_carries() {
        // Bit 4 set, and the leftover POS and CIGAR that samtools still writes.
        let text = record("read1", 4, "NC_002516.2", 4001, 0, "10M", "ACGTACGTAC");
        assert!(sam(&text, &region()).unwrap().is_empty());
    }

    #[test]
    fn a_record_with_no_cigar_is_skipped() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "*", "ACGTACGTAC");
        assert!(sam(&text, &region()).unwrap().is_empty());
    }

    #[test]
    fn a_record_on_another_sequence_is_skipped() {
        let text = format!(
            "{}\n{}\n",
            record("elsewhere", 0, "NC_002517.1", 4001, 60, "10M", "ACGTACGTAC"),
            record("here", 0, "NC_002516.2", 4001, 60, "10M", "ACGTACGTAC"),
        );
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].start, 4000);
    }

    #[test]
    fn a_read_off_the_side_of_the_window_is_left_out_and_one_across_it_is_not() {
        let text = format!(
            "{}\n{}\n{}\n",
            // Ends at 3800, before the window starts at 3899.
            record("before", 0, "NC_002516.2", 3751, 60, "50M", "*"),
            // Starts at 3880 and runs to 3930, so it crosses the left edge.
            record("across", 0, "NC_002516.2", 3881, 60, "50M", "*"),
            // Starts at 5000, well past the right edge at 4100.
            record("after", 0, "NC_002516.2", 5001, 60, "50M", "*"),
        );
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].start, 3880);
    }

    #[test]
    fn the_sequence_is_kept_unless_the_record_has_none() {
        let with = record("read1", 0, "NC_002516.2", 4001, 60, "10M", "ACGTACGTAC");
        let without = record("read2", 0, "NC_002516.2", 4001, 60, "10M", "*");
        assert_eq!(sam(&with, &region()).unwrap()[0].sequence, b"ACGTACGTAC");
        assert!(sam(&without, &region()).unwrap()[0].sequence.is_empty());
    }

    #[test]
    fn the_sequence_lines_up_with_the_reference_through_the_cigar() {
        // Soft clipped bases stay in the record, so the first aligned base is
        // read offset 3 and reference 4000.
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "3S7M", "TTTACGTACGT");
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads[0].base_at(4000), Some(b'A'));
        // Seven aligned bases, so the last one is read offset 9 at 4006.
        assert_eq!(reads[0].base_at(4006), Some(b'G'));
        assert_eq!(reads[0].base_at(4007), None, "past the end of the read");
    }

    #[test]
    fn a_mapping_quality_is_kept_and_the_missing_one_is_not() {
        let given = record("read1", 0, "NC_002516.2", 4001, 60, "10M", "*");
        let missing = record("read2", 0, "NC_002516.2", 4001, 255, "10M", "*");
        assert_eq!(sam(&given, &region()).unwrap()[0].mapping_quality, Some(60));
        assert_eq!(sam(&missing, &region()).unwrap()[0].mapping_quality, None);
    }

    #[test]
    fn the_header_is_not_data() {
        let text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:NC_002516.2\tLN:6264404\n{}\n",
            record("read1", 0, "NC_002516.2", 4001, 60, "10M", "ACGTACGTAC"),
        );
        let reads = sam(&text, &region()).unwrap();
        assert_eq!(reads.len(), 1);
    }

    #[test]
    fn a_cigar_that_is_not_one_says_which_line_it_was_on() {
        let text = format!(
            "@HD\tVN:1.6\n{}\n{}\n",
            record("read1", 0, "NC_002516.2", 4001, 60, "10M", "*"),
            record("read2", 0, "NC_002516.2", 4021, 60, "10Z", "*"),
        );
        let error = sam(&text, &region()).unwrap_err();
        assert_eq!(error.line, 3);
        assert!(
            error.to_string().contains("not a CIGAR operation"),
            "{error}"
        );
    }

    #[test]
    fn a_cigar_operation_without_a_length_is_an_error() {
        let text = record("read1", 0, "NC_002516.2", 4001, 60, "M10", "*");
        let error = sam(&text, &region()).unwrap_err();
        assert!(error.to_string().contains("needs a length"), "{error}");

        let trailing = record("read1", 0, "NC_002516.2", 4001, 60, "10M5", "*");
        let error = sam(&trailing, &region()).unwrap_err();
        assert!(error.to_string().contains("no operation"), "{error}");
    }

    #[test]
    fn a_record_with_columns_missing_is_an_error_and_not_a_panic() {
        let text = "read1\t0\tNC_002516.2\t4001\t60\t10M\n";
        let error = sam(text, &region()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("11 columns"), "{error}");
    }

    #[test]
    fn a_flag_or_a_position_that_is_not_a_number_says_so() {
        let flag = "read1\tforward\tNC_002516.2\t4001\t60\t10M\t*\t0\t0\t*\t*";
        assert!(sam(flag, &region())
            .unwrap_err()
            .to_string()
            .contains("the flag is not a number"));

        let position = record("read1", 0, "NC_002516.2", 0, 60, "10M", "*");
        let error = sam(&position, &region()).unwrap_err();
        assert!(error.to_string().contains("POS is 1-based"), "{error}");
    }

    #[test]
    fn an_empty_file_is_no_reads_rather_than_an_error() {
        assert!(sam("", &region()).unwrap().is_empty());
        assert!(sam("@HD\tVN:1.6\n", &region()).unwrap().is_empty());
    }
}
