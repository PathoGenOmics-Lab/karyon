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

use crate::{Association, Region, Variant};

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

/// What an association table held inside the window.
#[derive(Debug, Clone, PartialEq)]
pub struct Associations {
    /// The tested positions in the window, each with the statistic drawn: the
    /// value as written, or `-log10` of it where the file wrote p-values.
    pub points: Vec<Association>,
    /// Whether the value column held p-values, converted on the way in.
    ///
    /// A threshold is given in the units the file is in, so a threshold for a
    /// file of p-values is a p-value too and wants the same conversion.
    pub p_values: bool,
}

/// Reads association statistics from a table.
///
/// The points of [`association_table`], which says what the value column held
/// as well.
pub fn associations(text: &str, region: &Region) -> Result<Vec<Association>, ReadError> {
    association_table(text, region).map(|table| table.points)
}

/// Reads association statistics from a table, and says whether they were
/// p-values.
///
/// Two or three columns: a position, a value, and optionally a sequence name
/// in front of them. A header line naming the columns is allowed and skipped.
/// Positions are 1-based, as every association tool writes them.
///
/// A scan is drawn higher meaning stronger, so what comes back is a statistic
/// that grows with the evidence. The header says what the value column holds:
/// a column named as p-values are, `P`, `p_wald`, `P_BOLT_LMM`, `p.value`,
/// `P-value`, or a q-value or FDR column, is converted to `-log10` with
/// [`Association::from_p_value`], and anything else is drawn as written, a
/// `-log10(p)` or `LOG10P` column included. Every tool that writes a p-value
/// column names it one of these ways.
///
/// # Errors
///
/// The line that does not read: a position that is not one, a value that is
/// not a number, a table that is not two or three columns wide, or, in a
/// column of p-values, a value outside nought to one. And the whole file when
/// it has no header and every value in it lies between nought and one: that is
/// how a column of p-values looks, and drawn as written it put the strongest
/// hit at the bottom of the figure and exited nought, while nothing in the
/// file can say whether it was one.
pub fn association_table(text: &str, region: &Region) -> Result<Associations, ReadError> {
    // A table of more columns, as association tools write them, is read by
    // the names its header gives the columns.
    if let Some((line, wide)) = hashed_header(text) {
        return wide.read(text, region, line);
    }
    if let Some((line, head)) = lines(text).next() {
        let names = columns(head);
        if names.len() > 3 {
            return match Wide::of(&names) {
                Some(wide) => wide.read(text, region, line),
                None => Err(ReadError::at(
                    line,
                    format!(
                        "an association table of {} columns is read by its header, and none \
                         of these names a position (BP, POS) and a p-value or its logarithm \
                         (P, LOG10P): {}",
                        names.len(),
                        names
                            .iter()
                            .map(|name| name.trim())
                            .collect::<Vec<_>>()
                            .join(" ")
                    ),
                )),
            };
        }
    }
    let mut points: Vec<(u64, f64, usize)> = Vec::new();
    let mut first = true;
    let mut named: Option<&str> = None;
    // Whether every value in the file, read or not, lies between nought and
    // one, which is what decides a table with no header.
    let mut stated = false;
    let mut probabilities = true;
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
                named = Some(value.trim());
                continue;
            }
        }

        // Every row counts here, on any sequence and in any window, and read
        // leniently: a row this figure does not draw is not otherwise parsed,
        // and a scan whose peaks are all on another chromosome is still a scan.
        if let Ok(number) = value.trim().parse::<f64>() {
            if number.is_finite() {
                stated = true;
                probabilities &= (0.0..=1.0).contains(&number);
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
        points.push((pos, number::<f64>(value, "value", line)?, line));
    }

    let p_values = match named {
        Some(name) => names_p_values(name),
        None if stated && probabilities => {
            return Err(ReadError::whole(
                "every value lies between 0 and 1, as p-values do, and the table has no \
                 header to say whether they are; a scan is drawn as -log10(p), so name the \
                 value column on a first line: P to have the values drawn as -log10(P), or \
                 what they are, such as mlog10p, to draw them as written",
            ))
        }
        None => false,
    };

    let mut converted = Vec::with_capacity(points.len());
    for (pos, value, line) in points {
        if !p_values {
            converted.push(Association::new(pos, value));
            continue;
        }
        if !(0.0..=1.0).contains(&value) {
            return Err(ReadError::at(
                line,
                format!(
                    "the value column is named as p-values are, and a p-value lies \
                     between 0 and 1, not {value}; a column of -log10(p) is drawn as \
                     written when its name says so, as mlog10p does"
                ),
            ));
        }
        converted.push(Association::from_p_value(pos, value));
    }
    Ok(Associations {
        points: converted,
        p_values,
    })
}

/// A header written behind a `#`, as PLINK 2 writes `#CHROM POS ID ... P`,
/// and the line it is on.
///
/// [`lines`] passes over every line starting with `#` as a comment, which is
/// what one is above a table of two or three columns. Here it is the header
/// only when it names what a scan is drawn from, so a comment stays one.
fn hashed_header(text: &str) -> Option<(usize, Wide)> {
    let (line, head) = text
        .strip_prefix('\u{feff}')
        .unwrap_or(text)
        .lines()
        .enumerate()
        .map(|(index, line)| (index + 1, line.trim_end_matches('\r')))
        .find(|(_, line)| !line.trim().is_empty() && !line.starts_with("##"))?;
    let names = columns(head);
    if !head.starts_with('#') || names.len() <= 3 {
        return None;
    }
    Wide::of(&names).map(|wide| (line, wide))
}

/// The sequences an association table's rows are on, each once with how many
/// rows it has, in the order the file first names them; nothing for a table
/// that names none.
///
/// The column is the one the header names, which for BOLT-LMM is the second,
/// behind the variant's own name, and the first in a table of three.
pub fn association_sequences(text: &str) -> Vec<(String, usize)> {
    let (column, header) = match hashed_header(text) {
        Some((line, wide)) => match wide.sequence {
            Some(at) => (at, Some(line)),
            None => return Vec::new(),
        },
        None => {
            let Some((line, head)) = lines(text).next() else {
                return Vec::new();
            };
            let names = columns(head);
            match names.len() {
                3 => (0, names[1].trim().parse::<u64>().is_err().then_some(line)),
                wide if wide > 3 => match Wide::of(&names).and_then(|wide| wide.sequence) {
                    Some(at) => (at, Some(line)),
                    None => return Vec::new(),
                },
                _ => return Vec::new(),
            }
        }
    };
    let mut held: Vec<(String, usize)> = Vec::new();
    for (line, row) in lines(text) {
        if Some(line) == header {
            continue;
        }
        let fields = columns(row);
        let Some(name) = fields.get(column).map(|name| name.trim()) else {
            continue;
        };
        match held.iter_mut().find(|(seen, _)| seen == name) {
            Some((_, count)) => *count += 1,
            None => held.push((name.to_string(), 1)),
        }
    }
    held
}

/// Where a wide association table keeps what a scan draws, by the names its
/// header gives the columns.
///
/// PLINK writes `CHR SNP BP A1 ... P`, PLINK 2 `#CHROM POS ID ... P`, REGENIE
/// `CHROM GENPOS ... LOG10P`, BOLT `SNP CHR BP ... P_BOLT_LMM`, GEMMA
/// `chr rs ps ... p_wald`, SAIGE `CHR POS ... p.value`, and the GWAS Catalog
/// `chromosome base_pair_location ... p_value`: a position and a p-value, or
/// its logarithm, are in every one of them, and are found by their names.
struct Wide {
    sequence: Option<usize>,
    position: usize,
    value: usize,
    p_values: bool,
}

impl Wide {
    fn of(names: &[&str]) -> Option<Wide> {
        let lower: Vec<String> = names
            .iter()
            .map(|name| name.trim().to_ascii_lowercase())
            .collect();
        let find = |set: &[&str]| lower.iter().position(|name| set.contains(&name.as_str()));
        let position = find(&[
            "bp",
            "pos",
            "position",
            "base_pair_location",
            "genpos",
            "ps",
            "bp_hg19",
            "bp_hg38",
        ])?;
        let sequence = find(&[
            "chr",
            "#chr",
            "chrom",
            "#chrom",
            "chromosome",
            "seqname",
            "contig",
        ]);
        // A p-value first, since a table holding both was written to be read
        // by it, and the logarithm of one after.
        let (value, p_values) = names
            .iter()
            .position(|name| names_p_values(name.trim()))
            .map(|at| (at, true))
            .or_else(|| {
                lower
                    .iter()
                    .position(|name| name.contains("log") && name.contains('p'))
                    .map(|at| (at, false))
            })?;
        Some(Wide {
            sequence,
            position,
            value,
            p_values,
        })
    }

    /// The rows after the header, which is on line `header`.
    fn read(&self, text: &str, region: &Region, header: usize) -> Result<Associations, ReadError> {
        let widest = self
            .position
            .max(self.value)
            .max(self.sequence.unwrap_or(0));
        let mut points = Vec::new();
        for (line, row) in lines(text).filter(|(line, _)| *line != header) {
            let fields = columns(row);
            if fields.len() <= widest {
                return Err(ReadError::at(
                    line,
                    format!(
                        "this row has {} columns, and the header puts what is drawn in column {}",
                        fields.len(),
                        widest + 1
                    ),
                ));
            }
            if let Some(at) = self.sequence {
                if fields[at].trim() != region.seq() {
                    continue;
                }
            }
            let pos = position(fields[self.position].trim(), "position", line)?;
            if !region.contains(pos) {
                continue;
            }
            // A test the tool could not run is written as NA, and has nothing
            // to draw.
            let value = fields[self.value].trim();
            if matches!(value, "" | "." | "-" | "NA" | "na" | "nan" | "NaN") {
                continue;
            }
            let value: f64 = number(value, "value", line)?;
            if !self.p_values {
                points.push(Association::new(pos, value));
                continue;
            }
            if !(0.0..=1.0).contains(&value) {
                return Err(ReadError::at(
                    line,
                    format!("a p-value lies between 0 and 1, not {value}"),
                ));
            }
            points.push(Association::from_p_value(pos, value));
        }
        Ok(Associations {
            points,
            p_values: self.p_values,
        })
    }
}

/// Whether a column name says it holds p-values, or q-values, which a scan
/// draws as `-log10` of themselves.
///
/// A name that mentions a logarithm has had it taken already, so `LOG10P`,
/// `-log10(p)` and `mlog10p` are drawn as written.
fn names_p_values(name: &str) -> bool {
    let lower = name.to_ascii_lowercase();
    if lower.contains("log") {
        return false;
    }
    let squashed: String = lower.chars().filter(char::is_ascii_alphanumeric).collect();
    matches!(
        squashed.as_str(),
        "p" | "pval" | "pvals" | "pvalue" | "pvalues" | "q" | "qval" | "qvalue" | "fdr" | "padj"
    ) || lower.starts_with("p_")
        || lower.starts_with("p.")
        || lower.starts_with("p-")
        || squashed.ends_with("pvalue")
        || squashed.ends_with("pval")
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
        // A column called p holds p-values, drawn as -log10 of themselves.
        assert!((points[0].value - -(3.2e-9f64).log10()).abs() < 1e-12);
    }

    #[test]
    fn a_value_comes_through_as_written_since_the_track_draws_it_as_given() {
        // The table the format guide shows, spaces and all. Nothing between
        // the file and the track takes a logarithm, so a scan arrives as
        // -log10(p) already, in the units --threshold is given in.
        let text = "\
chrom         pos   neglog10p
Pf3D7_07_v3   4100  8.49
Pf3D7_08_v3   4110  12.0
Pf3D7_07_v3   4150  0.40
";
        let region = Region::parse("Pf3D7_07_v3:4,000-4,200").unwrap();
        let values: Vec<f64> = associations(text, &region)
            .unwrap()
            .iter()
            .map(|point| point.value)
            .collect();
        assert_eq!(values, vec![8.49, 0.40]);
    }

    #[test]
    fn a_table_with_no_header_reads_from_its_first_line() {
        let text = "4100\t8.5\n4150\t0.4\n";
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
        let text = "10\t0.5\n4100\t8.5\n900000\t0.2\n";
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

    /// The column's name says what it holds, and every tool that writes
    /// p-values names the column one of these ways. Drawn as written, a
    /// column of them put the strongest hit on the floor of the scan and the
    /// weakest at the top, and the command exited nought.
    #[test]
    fn a_column_named_as_p_values_are_is_drawn_as_minus_log10() {
        let region = Region::parse("chr1:1-1000").unwrap();
        for name in [
            "P",
            "p",
            "p_wald",
            "P_BOLT_LMM",
            "p.value",
            "P-value",
            "PVAL",
            "frequentist_add_pvalue",
            "FDR",
            "q",
            "padj",
        ] {
            let text = format!("chr\tpos\t{name}\nchr1\t100\t1e-8\nchr1\t200\t0.5\n");
            let table = association_table(&text, &region).unwrap();
            assert!(table.p_values, "{name} is a p-value column");
            let values: Vec<f64> = table.points.iter().map(|point| point.value).collect();
            assert!((values[0] - 8.0).abs() < 1e-12, "{name}: {values:?}");
            assert!(values[0] > values[1], "{name}: the strong hit is lower");
        }
    }

    #[test]
    fn a_column_named_as_anything_else_is_drawn_as_written() {
        let region = Region::parse("chr1:1-1000").unwrap();
        // The last three are named after a p-value and hold its logarithm,
        // which is why a name mentioning one is read as written first.
        for name in [
            "LOG10P",
            "-log10(p)",
            "mlog10p",
            "neglog10p",
            "score",
            "PIP",
            "pos_prob",
            "log10_pvalue",
            "-log10(pval)",
            "P_log10",
        ] {
            let text = format!("chr\tpos\t{name}\nchr1\t100\t0.9\nchr1\t200\t0.2\n");
            let table = association_table(&text, &region).unwrap();
            assert!(!table.p_values, "{name} is not a p-value column");
            let values: Vec<f64> = table.points.iter().map(|point| point.value).collect();
            assert_eq!(values, [0.9, 0.2], "{name}");
        }
    }

    #[test]
    fn a_p_value_outside_nought_to_one_is_refused_on_its_line() {
        // A column named P that holds -log10 values: converting 8.49 would
        // draw a point below the floor.
        let text = "pos\tP\n100\t0.3\n200\t8.49\n";
        let region = Region::parse("chr1:1-1000").unwrap();
        let error = association_table(text, &region).unwrap_err();
        assert_eq!(error.line, 3);
        assert!(error.to_string().contains("between 0 and 1"), "{error}");
    }

    #[test]
    fn a_table_with_no_header_and_only_values_under_one_is_refused() {
        // Nothing in the file says whether these are p-values, and they look
        // like them. Drawn as written, the strongest hit was the lowest.
        let region = Region::parse("chr1:1-1000").unwrap();
        let error = association_table("chr1\t100\t1e-8\nchr1\t200\t0.5\n", &region).unwrap_err();
        assert_eq!(error.line, 0, "the whole file, not one line");
        assert!(
            error.to_string().contains("P to have the values drawn"),
            "{error}"
        );
        // One value above one anywhere in the file, even on a sequence this
        // figure does not draw, and it is a statistic to draw as written.
        let text = "chr1\t100\t0.9\nchr1\t200\t0.5\nchr2\t300\t12.5\n";
        let table = association_table(text, &region).unwrap();
        assert!(!table.p_values);
        assert_eq!(table.points.len(), 2);
    }

    /// PLINK's own output, spaced as PLINK spaces it, with a test it could
    /// not run: the table is read by its header, and the NA is left out.
    #[test]
    fn an_association_tool_s_own_table_is_read_by_its_header() {
        let region = Region::parse("1:1-1000").unwrap();
        let plink = " CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR \n   1     rs1        100    A   0.3357   0.4500    G        3.776       0.1        1.625 \n   1     rs2        200    A   0.3357   0.4500    G           NA          NA           NA \n   1     rs3        300    A   0.3357   0.4500    G        3.776       1e-9       1.625 \n   2     rs4        100    A   0.3357   0.4500    G        3.776       0.5        1.625 \n";
        let table = association_table(plink, &region).unwrap();
        assert!(table.p_values);
        let points: Vec<(u64, f64)> = table.points.iter().map(|p| (p.pos, p.value)).collect();
        assert_eq!(points.len(), 2, "{points:?}");
        assert_eq!(points[0].0, 99);
        assert!((points[1].1 - 9.0).abs() < 1e-9);

        // REGENIE writes the logarithm, and it is drawn as written.
        let regenie = "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA\n1 150 rs1 A G 0.1 900 ADD 0.1 0.01 3.2 7.5 NA\n";
        let table = association_table(regenie, &region).unwrap();
        assert!(!table.p_values);
        assert_eq!(table.points[0].value, 7.5);

        // With no sequence column every row is on the figure's sequence, and
        // the header is still not a row.
        let unplaced = "SNP BP A1 P\nrs1 100 A 0.5\nrs2 250 A 1e-4\n";
        let table = association_table(unplaced, &region).unwrap();
        assert_eq!(table.points.len(), 2);

        // A wide header naming nothing to draw says what it did name.
        let error = association_table("A B C D\n1 2 3 4\n", &region).unwrap_err();
        assert!(error.to_string().contains("A B C D"), "{error}");
    }

    /// PLINK 2 writes its header behind a `#`, where every other reader here
    /// sees a comment. It read the first row of numbers as the header and
    /// refused the file.
    #[test]
    fn a_plink_2_header_behind_a_hash_is_read_as_the_header() {
        let region = Region::parse("1:1-1000").unwrap();
        let plink2 = "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n\
                      1\t100\trs1\tA\tG\tG\tADD\t500\t0.1\t0.02\t5.0\t1e-3\n\
                      1\t200\trs2\tA\tG\tG\tADD\t500\t0.3\t0.02\t9.0\t1e-9\n\
                      1\t300\trs3\tA\tG\tG\tADD\t500\t0.0\t0.02\t0.1\tNA\n\
                      2\t150\trs4\tA\tG\tG\tADD\t500\t0.3\t0.02\t9.0\t0.5\n";
        let table = association_table(plink2, &region).unwrap();
        assert!(table.p_values);
        let positions: Vec<u64> = table.points.iter().map(|point| point.pos).collect();
        assert_eq!(positions, [99, 199]);
        // Above a table of three columns, a line behind a `#` is a comment.
        let commented = "# a scan\nchrom pos P\n1 100 0.001\n";
        let table = association_table(commented, &region).unwrap();
        assert_eq!(table.points.len(), 1);
        assert!(table.p_values);
    }

    #[test]
    fn an_association_table_says_which_sequences_its_rows_are_on() {
        let owned = |pairs: &[(&str, usize)]| -> Vec<(String, usize)> {
            pairs
                .iter()
                .map(|(name, count)| (name.to_string(), *count))
                .collect()
        };
        let plink = "CHR SNP BP A1 P\n1 rs1 150 A 0.5\n1 rs2 480 A 1e-9\n2 rs3 5 A 0.1\n";
        assert_eq!(association_sequences(plink), owned(&[("1", 2), ("2", 1)]));
        let plink2 = "#CHROM POS ID P\n7 100 rs1 0.5\n";
        assert_eq!(association_sequences(plink2), owned(&[("7", 1)]));
        // BOLT-LMM keeps the chromosome second, behind the variant's name.
        let bolt = "SNP\tCHR\tBP\tP_BOLT_LMM\nrs1\t3\t100\t0.5\n";
        assert_eq!(association_sequences(bolt), owned(&[("3", 1)]));
        let named = "chrom pos P\nchr1 100 0.5\n";
        assert_eq!(association_sequences(named), owned(&[("chr1", 1)]));
        let bare = "chr1 100 0.5\nchr2 200 0.1\n";
        assert_eq!(
            association_sequences(bare),
            owned(&[("chr1", 1), ("chr2", 1)])
        );
        assert!(association_sequences("100 0.5\n200 0.1\n").is_empty());
        assert!(association_sequences("SNP BP A1 P\nrs1 100 A 0.5\n").is_empty());
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
