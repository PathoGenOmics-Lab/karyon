//! Point events: VCF calls and association statistics.
//!
//! Both formats count from one, so both lose a base on the way in: a VCF `POS`
//! and the position column of an association table land at `POS - 1`. Zero is
//! not a position in a file that starts at one, and taking one off it would
//! wrap rather than fail, so a zero is a malformed line and is refused as one.
//!
//! # Where a call sits and how far it reaches
//!
//! A call is drawn at one base and can be about several. A deletion is written
//! one base to the left of the bases it removes, so the window is tested
//! against what `REF` spells and not against the anchor alone; filtering on the
//! anchor would take a deletion that eats the first bases of the window out of
//! the figure entirely. A row with several alternates becomes one call per
//! alternate at the same position, each with its own fraction when `AF` gives
//! one per allele and the shared one when it gives a single number, since that
//! is what a caller writing one `AF` for the row means. A row with no `AF` at
//! all is drawn full height: a call is a call whether or not anyone measured
//! what proportion of the reads carried it.
//!
//! What a call is called is the annotator's word when there is one, from `ANN`
//! or from `BCSQ`, and otherwise the shape of the call read off `REF` against
//! `ALT`, which says whether it is a substitution, an insertion or a deletion
//! without needing anything annotated at all.
//!
//! # What is dropped and what is refused
//!
//! Rows on another sequence, rows outside the window, and rows with no
//! alternate allele go past without a word. The last of those is most of a
//! gVCF, and a reference block is a statement that nothing happened rather than
//! a call. A row that does not parse stops the read on its line, because a VCF
//! that cannot be read is not a VCF and a figure short of the calls it should
//! have shows nothing wrong on its face.
//!
//! An association table of two columns names no sequence, so every row in one
//! is taken to be about the sequence on display; three columns puts the name in
//! front and it is matched. A word where a position belongs is a header, but
//! only on the first line worth looking at, so the same word further down the
//! file is an error rather than a row that disappears.

use std::cmp::Ordering;

use karyon::{Association, Region, Variant};

use super::{columns, lines, number, ReadError};

/// Reads calls from VCF text.
///
/// VCF `POS` is 1-based, so the variant lands at `POS - 1`. The category is
/// taken from the `ANN` or `BCSQ` consequence when one is there, and otherwise
/// from the shape of the call: a substitution, an insertion or a deletion,
/// which is `REF` against `ALT` and needs no annotation.
///
/// The value is the allele fraction, from `AF` in `INFO` when present, and 1.0
/// when not, since a call with no fraction is a call.
///
/// Rows on another sequence than `region.seq()` are skipped.
pub fn variants(text: &str, region: &Region) -> Result<Vec<Variant>, ReadError> {
    let mut calls = Vec::new();
    for (line, row) in lines(text) {
        let fields = columns(row);
        if fields.len() < 8 {
            return Err(ReadError::at(
                line,
                format!(
                    "a VCF line has at least 8 columns, this one has {}",
                    fields.len()
                ),
            ));
        }
        if fields[0] != region.seq() {
            continue;
        }
        let pos = position(fields[1], "POS", line)?;
        // A row with no alternate allele is a reference block, which is most of
        // what a gVCF holds, and it is not a call.
        if fields[4] == "." {
            continue;
        }

        let (reference, info) = (fields[3], fields[7]);
        // A whole genome VCF is a normal thing to hand over when only a window
        // is being drawn, so the rest of it stops here rather than being
        // carried into a track that would not draw it. What is kept is what
        // REF spells, not the anchor base alone: a deletion is written one base
        // to the left of the bases it removes, so a call anchored just outside
        // the window can still be a call about the window.
        let spelled = pos + reference.len().max(1) as u64;
        if pos >= region.end() || spelled <= region.start() {
            continue;
        }
        let alternates: Vec<&str> = fields[4].split(',').collect();
        let fractions = allele_fractions(info, alternates.len(), line)?;
        for (index, alt) in alternates.iter().enumerate() {
            let category = consequence(info, alt).unwrap_or_else(|| shape(reference, alt));
            calls.push(Variant::new(pos).value(fractions[index]).category(category));
        }
    }
    Ok(calls)
}

/// Reads association statistics from a table.
///
/// Two or three columns: a position, a value, and optionally a sequence name
/// in front of them. A header line naming the columns is allowed and skipped.
/// Positions are 1-based, as every association tool writes them.
///
/// The value is taken as a p-value and handed over as it is; the track is what
/// decides to draw it on a log scale.
pub fn associations(text: &str, region: &Region) -> Result<Vec<Association>, ReadError> {
    let mut points = Vec::new();
    let mut first = true;
    for (line, row) in lines(text) {
        let fields = columns(row);
        let (sequence, at, value) = match fields.as_slice() {
            [at, value] => (None, *at, *value),
            [sequence, at, value] => (Some(*sequence), *at, *value),
            other => {
                return Err(ReadError::at(
                    line,
                    format!(
                        "an association table is a position and a value, or a sequence \
                         name and those two, and this line has {} columns",
                        other.len()
                    ),
                ))
            }
        };

        // A header is told from data by its position column: a word there is a
        // column name and a number is a position. Only the first line worth
        // looking at gets to be a header, so a word further down the file is
        // still an error rather than a line that disappears.
        if first {
            first = false;
            if at.parse::<u64>().is_err() {
                continue;
            }
        }

        if let Some(sequence) = sequence {
            if sequence != region.seq() {
                continue;
            }
        }
        let pos = position(at, "position", line)?;
        if !region.contains(pos) {
            continue;
        }
        points.push(Association::new(pos, number::<f64>(value, "value", line)?));
    }
    Ok(points)
}

/// A 1-based coordinate as the 0-based one the rest of the crate counts in.
///
/// This is the subtraction the file exists to get right, so it is written once
/// rather than at each call site. Zero is not a position in a 1-based file and
/// taking one off it would wrap, so it is a malformed line and says so.
fn position(field: &str, what: &str, line: usize) -> Result<u64, ReadError> {
    let one_based: u64 = number(field, what, line)?;
    one_based
        .checked_sub(1)
        .ok_or_else(|| ReadError::at(line, format!("{what} is 1-based and starts at 1, not 0")))
}

/// The value of one `INFO` key, or `None` when the key is not there.
///
/// The keys have to be matched whole. `AF` is the end of `MLEAF` and of
/// `AF_ESP`, so a search for the substring finds the wrong number in a file
/// written by GATK or annotated against a population database.
fn info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';')
        .filter_map(|entry| entry.split_once('='))
        .find(|(name, _)| *name == key)
        .map(|(_, value)| value)
}

/// The fraction each alternate allele gets.
///
/// `AF` carries one number per alternate allele, so a multi-allelic row hands
/// each alternate its own. A single number is spread over all of them, since
/// that is what a caller writing one `AF` for the row means. Any other count
/// is a row that does not add up, and saying so beats guessing which allele
/// the numbers belong to.
fn allele_fractions(info: &str, alternates: usize, line: usize) -> Result<Vec<f64>, ReadError> {
    let Some(field) = info_field(info, "AF") else {
        // A call with no fraction is a call, and it gets a full height stem.
        return Ok(vec![1.0; alternates]);
    };
    let mut values = Vec::with_capacity(alternates);
    for value in field.split(',') {
        values.push(number::<f64>(value, "AF", line)?);
    }
    match values.len() {
        1 => Ok(vec![values[0]; alternates]),
        given if given == alternates => Ok(values),
        given => Err(ReadError::at(
            line,
            format!("AF has {given} values for {alternates} alternate alleles"),
        )),
    }
}

/// The consequence an annotator wrote for this alternate allele, if any.
///
/// `ANN`, which snpEff and VEP both write, is one entry per allele per feature,
/// entries separated by commas and fields by pipes, with the allele in field
/// one and the consequence in field two. The entry naming this allele is the
/// one to read, and the first entry is the fallback for an annotator that
/// trimmed the allele down to something the VCF row does not spell.
///
/// `BCSQ`, from `bcftools csq`, puts the consequence first instead, marks a
/// consequence it is unsure of with a leading `*`, and writes pointer entries
/// such as `@761154` that carry no fields of their own.
fn consequence(info: &str, alt: &str) -> Option<String> {
    let word = |text: &str| {
        let text = text.trim().trim_start_matches('*');
        (!text.is_empty()).then(|| text.to_string())
    };

    if let Some(value) = info_field(info, "ANN") {
        let mut annotated = value.split(',').filter(|entry| entry.contains('|'));
        let first = annotated.next();
        let chosen = value
            .split(',')
            .find(|entry| entry.contains('|') && entry.split('|').next() == Some(alt))
            .or(first);
        if let Some(found) = chosen
            .and_then(|entry| entry.split('|').nth(1))
            .and_then(word)
        {
            return Some(found);
        }
    }

    if let Some(value) = info_field(info, "BCSQ") {
        return value
            .split(',')
            .filter(|entry| !entry.starts_with('@'))
            .filter_map(|entry| entry.split('|').next())
            .find_map(word);
    }

    None
}

/// The category a call gets when nothing annotated it, which is what `REF` and
/// `ALT` say between them.
///
/// A symbolic allele such as `<DEL>` carries no sequence to measure, so the tag
/// inside the brackets answers instead of the lengths, which would call a five
/// character `<DEL>` against a one base `REF` an insertion. A breakend names
/// the other side of a join in square brackets and is neither, and the star
/// allele is the one an overlapping deletion took away.
fn shape(reference: &str, alt: &str) -> String {
    if alt.contains('[') || alt.contains(']') {
        return "breakend".to_string();
    }
    if let Some(tag) = alt.strip_prefix('<').and_then(|alt| alt.strip_suffix('>')) {
        let tag = tag.split(':').next().unwrap_or(tag);
        return match tag {
            "DEL" => "deletion".to_string(),
            "INS" => "insertion".to_string(),
            other => other.to_lowercase(),
        };
    }
    if alt == "*" {
        return "deletion".to_string();
    }
    match alt.len().cmp(&reference.len()) {
        Ordering::Equal => "substitution".to_string(),
        Ordering::Greater => "insertion".to_string(),
        Ordering::Less => "deletion".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One annotated call in rpoB, the position the whole conversion turns on.
    const RPOB: &str = "\
##fileformat=VCFv4.2
##contig=<ID=NC_000962.3,length=4411532>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_000962.3\t761155\t.\tC\tT\t900\tPASS\tDP=54;AF=0.98;ANN=T|missense_variant|MODERATE|rpoB
";

    fn rpob_region() -> Region {
        Region::parse("NC_000962.3:761,000-763,000").unwrap()
    }

    #[test]
    fn a_vcf_position_is_one_base_to_the_left_of_what_the_file_says() {
        let calls = variants(RPOB, &rpob_region()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(
            calls[0].pos, 761_154,
            "POS 761155 is 761154 counting from zero"
        );
    }

    #[test]
    fn the_annotated_consequence_and_the_allele_fraction_come_through() {
        let calls = variants(RPOB, &rpob_region()).unwrap();
        assert_eq!(calls[0].category.as_deref(), Some("missense_variant"));
        assert_eq!(calls[0].value, Some(0.98));
    }

    #[test]
    fn an_allele_that_no_annotation_describes_is_read_off_ref_and_alt() {
        let text = "\
NC_045512.2\t21563\t.\tA\tG\t.\t.\t.
NC_045512.2\t21570\t.\tA\tAGGT\t.\t.\tDP=8
NC_045512.2\t21580\t.\tACGT\tA\t.\t.\tDP=9
";
        let region = Region::parse("NC_045512.2:21,000-22,000").unwrap();
        let calls = variants(text, &region).unwrap();
        let categories: Vec<&str> = calls
            .iter()
            .map(|call| call.category.as_deref().unwrap())
            .collect();
        assert_eq!(categories, vec!["substitution", "insertion", "deletion"]);
        assert_eq!(calls[1].pos, 21_569);
    }

    #[test]
    fn a_symbolic_allele_is_not_measured_by_its_spelling() {
        assert_eq!(shape("C", "<DEL>"), "deletion");
        assert_eq!(shape("C", "<INS:ME:ALU>"), "insertion");
        assert_eq!(shape("C", "<INV>"), "inv");
        assert_eq!(shape("C", "C[chrIV:200["), "breakend");
        assert_eq!(shape("CTTT", "*"), "deletion");
    }

    #[test]
    fn every_alternate_allele_is_its_own_variant_at_the_same_place() {
        let text = "chrIV\t900\t.\tC\tT,G\t.\t.\tAF=0.7,0.2\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants(text, &region).unwrap();
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].pos, 899);
        assert_eq!(calls[1].pos, 899);
        assert_eq!(calls[0].value, Some(0.7));
        assert_eq!(calls[1].value, Some(0.2));
    }

    #[test]
    fn one_fraction_written_for_a_multi_allelic_row_is_shared_by_its_alleles() {
        let text = "chrIV\t900\t.\tC\tT,G\t.\t.\tAF=0.5\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants(text, &region).unwrap();
        assert_eq!(calls[0].value, Some(0.5));
        assert_eq!(calls[1].value, Some(0.5));
    }

    #[test]
    fn a_row_whose_fractions_do_not_match_its_alleles_says_so() {
        let text = "chrIV\t900\t.\tC\tT,G,A\t.\t.\tAF=0.5,0.2\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let error = variants(text, &region).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("AF has 2 values"), "{error}");
    }

    #[test]
    fn a_call_with_no_fraction_is_still_a_call() {
        let text = "chrIV\t900\t.\tC\tT\t.\t.\tDP=30\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        assert_eq!(variants(text, &region).unwrap()[0].value, Some(1.0));
    }

    #[test]
    fn a_key_ending_in_af_is_not_af() {
        assert_eq!(info_field("DP=8;MLEAF=0.3", "AF"), None);
        assert_eq!(info_field("MLEAF=0.3;AF=0.9", "AF"), Some("0.9"));
        assert_eq!(info_field("AF=0.9;AF_ESP=0.01", "AF"), Some("0.9"));
    }

    #[test]
    fn a_fraction_that_is_not_a_number_is_an_error_and_not_a_default() {
        let text = "chrIV\t900\t.\tC\tT\t.\t.\tAF=high\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let error = variants(text, &region).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("AF"), "{error}");
    }

    #[test]
    fn the_annotation_naming_this_allele_is_the_one_read() {
        let info = "ANN=T|missense_variant|MODERATE|katG,G|stop_gained|HIGH|katG";
        assert_eq!(consequence(info, "G").as_deref(), Some("stop_gained"));
        assert_eq!(consequence(info, "T").as_deref(), Some("missense_variant"));
        // An annotator that trimmed the allele leaves the first entry to use.
        assert_eq!(consequence(info, "A").as_deref(), Some("missense_variant"));
    }

    #[test]
    fn bcftools_csq_puts_the_consequence_first_instead() {
        let text = "chrIV\t900\t.\tC\tT\t.\t.\tBCSQ=@899,*missense|YDL0|YDL0W|protein_coding\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants(text, &region).unwrap();
        assert_eq!(calls[0].category.as_deref(), Some("missense"));
    }

    #[test]
    fn rows_on_another_sequence_or_outside_the_window_are_not_data_for_this_figure() {
        let text = "\
chrI\t900\t.\tC\tT\t.\t.\t.
chrIV\t50\t.\tC\tT\t.\t.\t.
chrIV\t900\t.\tC\tT\t.\t.\t.
chrIV\t5000\t.\tC\tT\t.\t.\t.
";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants(text, &region).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].pos, 899);
    }

    #[test]
    fn a_reference_block_is_not_a_variant() {
        let text = "chrIV\t900\t.\tC\t.\t.\t.\tEND=950\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        assert!(variants(text, &region).unwrap().is_empty());
    }

    #[test]
    fn a_line_too_short_to_be_a_vcf_line_says_which_line_it_was() {
        let text = "chrIV\t900\t.\tC\tT\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let error = variants(text, &region).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("at least 8 columns"), "{error}");
    }

    #[test]
    fn a_position_of_zero_is_not_a_position_in_a_one_based_file() {
        let text = "chrIV\t0\t.\tC\tT\t.\t.\t.\n";
        let region = Region::parse("chrIV:1-1000").unwrap();
        let error = variants(text, &region).unwrap_err();
        assert!(error.to_string().contains("1-based"), "{error}");
    }

    #[test]
    fn an_association_position_is_one_base_to_the_left_of_what_the_file_says() {
        let text = "\
pos\tp
4100\t3.2e-9
";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let points = associations(text, &region).unwrap();
        assert_eq!(points.len(), 1);
        assert_eq!(
            points[0].pos, 4_099,
            "position 4100 is 4099 counting from zero"
        );
        assert_eq!(points[0].value, 3.2e-9);
    }

    #[test]
    fn a_table_with_no_header_reads_from_its_first_line() {
        let text = "4100\t3.2e-9\n4150\t0.4\n";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let points = associations(text, &region).unwrap();
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].pos, 4_099);
        assert_eq!(points[1].pos, 4_149);
    }

    #[test]
    fn a_third_column_in_front_is_the_sequence_and_is_matched_against_the_region() {
        let text = "\
chrom\tpos\tpvalue
Pf3D7_07_v3\t4100\t3.2e-9
Pf3D7_08_v3\t4110\t1.0e-12
Pf3D7_07_v3\t4150\t0.4
";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let points = associations(text, &region).unwrap();
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].pos, 4_099);
        assert_eq!(points[1].pos, 4_149);
    }

    #[test]
    fn points_outside_the_window_are_left_out() {
        let text = "10\t0.5\n4100\t3.2e-9\n900000\t0.2\n";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let points = associations(text, &region).unwrap();
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].pos, 4_099);
    }

    #[test]
    fn a_word_where_a_position_belongs_is_a_header_only_on_the_first_line() {
        let text = "\
pos\tp
4100\t3.2e-9
locus\t0.4
";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let error = associations(text, &region).unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.to_string().contains("position"), "{error}");
    }

    #[test]
    fn a_value_that_is_not_a_number_is_an_error_and_not_a_missing_point() {
        let text = "pos\tp\n4100\tNA\n";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let error = associations(text, &region).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("value"), "{error}");
    }

    #[test]
    fn a_table_that_is_not_two_or_three_columns_says_how_wide_it_was() {
        let text = "Pf3D7_07_v3\t4100\t3.2e-9\tA\tT\n";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let error = associations(text, &region).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("5 columns"), "{error}");
    }
}
