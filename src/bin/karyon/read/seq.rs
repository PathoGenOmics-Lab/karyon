//! FASTA, both the single sequence kind and the aligned kind.
//!
//! FASTA carries no coordinates of its own, so there is nothing to convert
//! here. A record starts at its own first base, which means byte `n` of what
//! comes back is the base at 0-based position `n`, and the 1-based position `p`
//! of a locus string is byte `p - 1`. The cut down to the region on display is
//! done later, in `stack`, by indexing what these functions return, so an
//! off-by-one introduced here would move every base of the figure.

use super::ReadError;

/// One record, with the line its header sat on.
///
/// The line number is not part of what the two public functions return, but an
/// error about a record has to say where that record is, so it is carried this
/// far and dropped at the end.
struct Record {
    /// The line the `>` was on, counting from one.
    line: usize,
    /// The header up to the first whitespace, without the `>`.
    name: String,
    /// The sequence lines, joined.
    residues: Vec<u8>,
}

/// Reads every record of a FASTA file, in the order they appear.
///
/// The name is everything on the header line up to the first whitespace.
/// Sequence lines are joined with no separator and their case is kept, since a
/// soft-masked reference says something by being lower case.
pub fn fasta(text: &str) -> Result<Vec<(String, Vec<u8>)>, ReadError> {
    Ok(records(text)?
        .into_iter()
        .map(|record| (record.name, record.residues))
        .collect())
}

/// Reads an alignment, checking that it is one.
///
/// Every record has to be the same length, which is what makes it an alignment
/// rather than a set of sequences, and a file where they are not says which
/// record broke it and by how much.
pub fn alignment(text: &str) -> Result<Vec<(String, Vec<u8>)>, ReadError> {
    let parsed = records(text)?;
    // The first record sets the width. Nothing to compare against when the file
    // held no records at all, and an empty stack is the caller's to complain
    // about, since it knows which flag asked for the file.
    let Some(first) = parsed.first() else {
        return Ok(Vec::new());
    };
    let width = first.residues.len();
    let reference = first.name.clone();
    for record in parsed.iter().skip(1) {
        let length = record.residues.len();
        if length == width {
            continue;
        }
        let (direction, by) = if length < width {
            ("shorter", width - length)
        } else {
            ("longer", length - width)
        };
        return Err(ReadError::at(
            record.line,
            format!(
                "an alignment has every record the same length, and {:?} is {by} {direction} than {reference:?}, which is {width} columns",
                record.name
            ),
        ));
    }
    Ok(parsed
        .into_iter()
        .map(|record| (record.name, record.residues))
        .collect())
}

/// The records of a FASTA file, with the line each header was on.
///
/// Blank lines, and the comment lines a generated file sometimes carries, are
/// dropped by [`super::lines`] before they get here, which is why a record can
/// have a blank line in the middle of it and still read as one sequence.
fn records(text: &str) -> Result<Vec<Record>, ReadError> {
    let mut parsed: Vec<Record> = Vec::new();
    for (number, line) in super::lines(text) {
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("");
            if name.is_empty() {
                return Err(ReadError::at(
                    number,
                    "a FASTA header needs a name after the >",
                ));
            }
            parsed.push(Record {
                line: number,
                name: name.to_string(),
                residues: Vec::new(),
            });
            continue;
        }
        // Whitespace is never a residue, so a line padded at either end is the
        // same sequence as one that is not.
        let residues = line.trim();
        match parsed.last_mut() {
            Some(record) => record.residues.extend_from_slice(residues.as_bytes()),
            None => {
                return Err(ReadError::at(
                    number,
                    format!(
                        "a FASTA file starts with a > header, and this line is {}",
                        preview(residues)
                    ),
                ))
            }
        }
    }
    // A header with nothing after it is a truncated file rather than a record
    // of no bases, and saying so is more use than a track that draws nothing.
    for record in &parsed {
        if record.residues.is_empty() {
            return Err(ReadError::at(
                record.line,
                format!("record {:?} has a header and no sequence", record.name),
            ));
        }
    }
    Ok(parsed)
}

/// The start of a line, quoted, for an error message.
///
/// A FASTA file can hold a whole chromosome on one line, and an error is easier
/// to read than a megabase of it.
fn preview(line: &str) -> String {
    const KEEP: usize = 24;
    if line.chars().count() <= KEEP {
        format!("{line:?}")
    } else {
        let head: String = line.chars().take(KEEP).collect();
        format!("{head:?}...")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Wrapped at 60 columns, which is what most references are written at.
    /// The second line starts "CG", so the base at 1-based position 62 is a G.
    const YEAST: &str = "\
>chrI Saccharomyces cerevisiae S288C
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
CGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
";

    #[test]
    fn the_lines_of_a_record_are_joined_with_nothing_between_them() {
        let records = fasta(">lambda\nACGT\nTTGA\n").unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].0, "lambda");
        assert_eq!(records[0].1, b"ACGTTTGA".to_vec());
    }

    #[test]
    fn a_soft_masked_reference_keeps_its_lower_case() {
        let records = fasta(">chr2L Drosophila melanogaster\nACGTacgtACGT\n").unwrap();
        assert_eq!(records[0].1, b"ACGTacgtACGT".to_vec());
    }

    #[test]
    fn the_name_stops_at_the_first_whitespace() {
        let records = fasta(">NC_000962.3 Mycobacterium tuberculosis H37Rv\nACGT\n").unwrap();
        assert_eq!(records[0].0, "NC_000962.3");
    }

    #[test]
    fn the_base_at_1_based_position_62_is_byte_61() {
        let records = fasta(YEAST).unwrap();
        let bases = &records[0].1;
        // Two lines of sixty, so the join has to be exact for this to hold.
        assert_eq!(bases.len(), 120);
        assert_eq!(bases[0], b'A', "1-based position 1 is byte 0");
        assert_eq!(bases[59], b'A', "1-based position 60 is byte 59");
        assert_eq!(bases[60], b'C', "1-based position 61 is byte 60");
        assert_eq!(bases[61], b'G', "1-based position 62 is byte 61");
    }

    #[test]
    fn several_records_come_back_in_the_order_they_were_written() {
        let text = "\
>plasmid_a
AAAA
>plasmid_b
CCCC
>plasmid_c
GGGG
";
        let names: Vec<String> = fasta(text).unwrap().into_iter().map(|(n, _)| n).collect();
        assert_eq!(names, vec!["plasmid_a", "plasmid_b", "plasmid_c"]);
    }

    #[test]
    fn a_file_that_does_not_start_with_a_header_says_which_line() {
        let error = fasta("ACGTACGT\n>chr1\nACGT\n").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(
            error.to_string().contains("starts with a > header"),
            "{error}"
        );
    }

    #[test]
    fn a_header_with_no_sequence_after_it_names_the_record() {
        let text = "\
>Pf3D7_01_v3
ACGT
>Pf3D7_02_v3
>Pf3D7_03_v3
ACGT
";
        let error = fasta(text).unwrap_err();
        // The header of the record that is empty, not the one that follows it.
        assert_eq!(error.line, 3);
        assert!(error.to_string().contains("Pf3D7_02_v3"), "{error}");
    }

    #[test]
    fn a_header_at_the_end_of_the_file_is_the_same_error() {
        let error = fasta(">chrM\nACGT\n>chrX\n").unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.to_string().contains("no sequence"), "{error}");
    }

    #[test]
    fn a_header_with_no_name_is_an_error_rather_than_a_nameless_record() {
        let error = fasta(">\nACGT\n").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("needs a name"), "{error}");
    }

    #[test]
    fn carriage_returns_do_not_become_bases() {
        let records = fasta(">Bacillus_subtilis\r\nACGT\r\nTTGA\r\n").unwrap();
        assert_eq!(records[0].0, "Bacillus_subtilis");
        assert_eq!(records[0].1, b"ACGTTTGA".to_vec());
    }

    #[test]
    fn a_blank_line_inside_a_record_does_not_end_it() {
        let text = "\
>Arabidopsis_thaliana_chr1
ACGT

TTGA

CCGG
";
        let records = fasta(text).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].1, b"ACGTTTGACCGG".to_vec());
    }

    #[test]
    fn an_empty_file_is_no_records_rather_than_an_error() {
        assert!(fasta("").unwrap().is_empty());
        assert!(fasta("\n\n").unwrap().is_empty());
        assert!(alignment("").unwrap().is_empty());
    }

    #[test]
    fn an_alignment_of_equal_records_reads_like_any_other_fasta() {
        let text = "\
>sample_01
ACGT-ACGT
>sample_02
ACGTTACGT
>sample_03
ACGT-ACGA
";
        let records = alignment(text).unwrap();
        assert_eq!(records.len(), 3);
        assert!(records.iter().all(|(_, row)| row.len() == 9));
        // Column 5 of the alignment, which is byte 4 of every row.
        assert_eq!(records[1].1[4], b'T');
    }

    #[test]
    fn a_short_record_is_named_and_the_difference_given() {
        let text = "\
>gyrA_ref
ACGTACGTAC
>gyrA_alt
ACGTACG
";
        let error = alignment(text).unwrap_err();
        assert_eq!(error.line, 3);
        let message = error.to_string();
        assert!(message.contains("gyrA_alt"), "{message}");
        assert!(message.contains("3 shorter"), "{message}");
        assert!(message.contains("gyrA_ref"), "{message}");
    }

    #[test]
    fn a_long_record_is_named_too() {
        let text = ">a\nACGT\n>b\nACGT\n>c\nACGTAC\n";
        let error = alignment(text).unwrap_err();
        assert_eq!(error.line, 5);
        let message = error.to_string();
        assert!(message.contains("\"c\" is 2 longer"), "{message}");
    }

    #[test]
    fn an_alignment_reads_its_records_the_same_way_fasta_does() {
        let text = ">one\nac\ngt\n>two\nAC\nGT\n";
        assert_eq!(alignment(text).unwrap(), fasta(text).unwrap());
    }

    #[test]
    fn a_comment_line_is_not_a_base() {
        let records = fasta("# written by some pipeline\n>chr1\nACGT\n").unwrap();
        assert_eq!(records[0].1, b"ACGT".to_vec());
    }

    #[test]
    fn a_long_line_is_not_printed_whole_in_an_error() {
        let line = "A".repeat(4000);
        let error = fasta(&format!("{line}\n")).unwrap_err();
        assert!(error.reason.len() < 100, "{}", error.reason);
        assert!(error.reason.contains("..."), "{}", error.reason);
    }
}
