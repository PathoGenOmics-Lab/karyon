//! Protein domains, from the table `InterProScan` writes.
//!
//! Eleven columns and up to four more, tab separated, no header:
//!
//! ```text
//! 1  protein accession       8  stop, 1-based inclusive
//! 2  sequence MD5            9  score, often a dot
//! 3  sequence length        10  status
//! 4  analysis               11  date
//! 5  signature accession    12  InterPro accession   (may be absent)
//! 6  signature description  13  InterPro description (may be absent)
//! 7  start, 1-based inclusive
//! ```
//!
//! Start and stop count from one and include their end, so a start comes back
//! one lower and a stop is unchanged, the same conversion GFF3 gets.
//!
//! Split on tabs and never on whitespace, because column six is a sentence:
//! `BRCA2, oligonucleotide/oligosaccharide-binding, domain 1` is one field and
//! seven words. A row with no tab in it is refused rather than taken apart on
//! its spaces.
//!
//! # The axis is residues
//!
//! A domain is at a place in a protein, not at a place in a genome, so the
//! window a figure of these is drawn over is a residue range: `Region::parse`
//! takes `P00533:1-1210` as readily as it takes a chromosome, and the ruler
//! underneath counts amino acids.
//!
//! # Column one names the row rather than selecting it
//!
//! As in [`locus`](super::locus), and for the same reason: the figure is the
//! comparison. Every protein in the file is a row and they share one residue
//! axis, so the architectures can be read against each other, which is what
//! makes a domain gained or lost visible at all.
//!
//! # A length that is absent is not a length of nought
//!
//! Column three is how long the protein is, and it is what the row's backbone
//! is drawn from. The furthest domain is not that number: a protein whose last
//! annotated domain ends at residue 300 may run to 800, and drawing the
//! backbone to 300 says the domain reaches the C terminus, which is a claim
//! about the protein rather than about the annotation. A length of nought
//! removes the backbone and every domain on the row, leaving the name standing
//! over an empty line, so both are refused here rather than drawn.
//!
//! # One protein, many analyses
//!
//! `InterProScan` runs a dozen member databases at once and they annotate the
//! same residues: `Pfam`, `PANTHER`, `SUPERFAMILY`, `Gene3D` and the rest all
//! describe one kinase domain in their own words. Drawn together they are five
//! boxes over one region of the protein, each named differently, and the
//! architecture is unreadable. So a read names its analysis: [`analyses`] says
//! what a file holds and [`architectures`] takes the one to draw.
//!
//! # What is not read here
//!
//! `hmmscan --domtblout`, the other route to Pfam. Which of its columns holds
//! the protein depends on whether `hmmscan` or `hmmsearch` wrote it, and the
//! only thing that says so is a comment line at the foot of the file, which is
//! the first thing to go when output is piped or pasted. Read with the wrong
//! orientation the row is built against the length of the model rather than of
//! the protein, and every domain past the model's end is clipped off a backbone
//! that is hundreds of residues too short. `Pfam` is reachable through this
//! reader as one of the analyses instead.

use std::collections::BTreeMap;

use crate::{DomainArchitecture, DomainFeature, Region};

use super::{lines, number, ReadError};

/// The architectures a file holds for one analysis, and what it held besides.
#[derive(Debug, Clone, PartialEq)]
pub struct Domains {
    /// One row per protein, in the order the proteins were first seen.
    pub rows: Vec<DomainArchitecture>,
    /// The proteins the file named, in first-seen order.
    pub proteins: Vec<String>,
    /// Rows in the file, before any filter.
    pub records: usize,
    /// Rows from an analysis this read did not ask for.
    pub other_analysis: usize,
    /// Rows whose domain touches no residue of the window.
    pub off_region: usize,
}

/// Every analysis in the file, with how many rows each has.
///
/// Ordered, so the same file answers the same way twice. The question that
/// comes first, since a dozen member databases annotating one protein are a
/// dozen boxes over one region of it and no architecture at all.
pub fn analyses(text: &str) -> Result<BTreeMap<String, usize>, ReadError> {
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for (at, line) in lines(text) {
        let cols = fields(line, at)?;
        *seen.entry(cols[3].to_string()).or_insert(0) += 1;
    }
    Ok(seen)
}

/// The architectures for one analysis, inside one window.
///
/// Start and stop count from one and include their end, so a start comes back
/// one lower. Column one names the protein whose row a domain belongs to rather
/// than selecting which rows are drawn, so every protein in the file is a row.
///
/// # Errors
///
/// Returns the first row that is not an `InterProScan` row, whose start is
/// nought, whose stop is before its start, whose length will not parse or is
/// nought, or that gives a length for a protein that differs from the one an
/// earlier row gave. A file that parses and holds nothing for this analysis in
/// this window is not an error here: [`Domains`] says which of the reasons it
/// was.
pub fn architectures(text: &str, region: &Region, wanted: &str) -> Result<Domains, ReadError> {
    let mut found = Domains {
        rows: Vec::new(),
        proteins: Vec::new(),
        records: 0,
        other_analysis: 0,
        off_region: 0,
    };
    let mut order: Vec<String> = Vec::new();
    let mut lengths: BTreeMap<String, u64> = BTreeMap::new();
    let mut features: BTreeMap<String, Vec<DomainFeature>> = BTreeMap::new();

    for (at, line) in lines(text) {
        let cols = fields(line, at)?;
        found.records += 1;

        let protein = cols[0];
        let length: u64 = number(cols[2], "sequence length", at)?;
        if length == 0 {
            return Err(ReadError::at(
                at,
                "a protein of no length has no row to draw a domain on",
            ));
        }
        // Recorded for every row, whatever analysis it came from, so the file
        // is checked against itself rather than against the subset that is
        // being drawn.
        match lengths.get(protein) {
            None => {
                lengths.insert(protein.to_string(), length);
                order.push(protein.to_string());
                features.entry(protein.to_string()).or_default();
            }
            Some(before) if *before != length => {
                return Err(ReadError::at(
                    at,
                    format!("{protein} is {length} long here and {before} long earlier"),
                ))
            }
            Some(_) => {}
        }

        if cols[3] != wanted {
            found.other_analysis += 1;
            continue;
        }

        let stated: u64 = number(cols[6], "start", at)?;
        let stop: u64 = number(cols[7], "stop", at)?;
        let Some(start) = stated.checked_sub(1) else {
            return Err(ReadError::at(
                at,
                "InterProScan counts from 1, so 0 is not a start",
            ));
        };
        // An inverted span is drawn as nothing at all here, which is the same
        // as a protein nobody annotated, and swapping the two would invent a
        // domain out of a row that names none.
        if stop < stated {
            return Err(ReadError::at(at, "stop is before start"));
        }
        if start >= region.end() || stop <= region.start() {
            found.off_region += 1;
            continue;
        }

        let mut feature = DomainFeature::new(start, stop);
        // The signature description, then its accession. A dot or an empty
        // field is what the table writes where there is nothing, and left in it
        // becomes the name of a domain: the track draws the string it is given
        // and calls an unnamed one `domain`, so both spellings of nothing have
        // to go before it gets there.
        if let Some(label) = [cols[5], cols[4]]
            .iter()
            .map(|field| field.trim())
            .find(|field| !field.is_empty() && *field != "." && *field != "-")
        {
            feature = feature.label(label);
        }
        features
            .entry(protein.to_string())
            .or_default()
            .push(feature);
    }

    found.proteins = order.clone();
    found.rows = order
        .into_iter()
        .map(|name| {
            let length = lengths[&name];
            let mut row = DomainArchitecture::new(name.clone(), length);
            for feature in features.remove(&name).unwrap_or_default() {
                row = row.feature(feature);
            }
            row
        })
        .collect();
    Ok(found)
}

/// The columns of one row, split on tabs alone.
///
/// The shared splitter falls back to whitespace when a line holds no tab, and
/// column six is a sentence, so that fallback turns one description into seven
/// fields and every column after it into the wrong one.
fn fields(line: &str, at: usize) -> Result<Vec<&str>, ReadError> {
    if !line.contains('\t') {
        return Err(ReadError::at(
            at,
            "an InterProScan row is tab separated, and this line holds no tab",
        ));
    }
    let cols: Vec<&str> = line.split('\t').collect();
    // Eleven as standard, and up to fifteen with the InterPro entry, the GO
    // terms and the pathways. A row that stops short of eleven is missing the
    // date or the status, and one of the columns before them is then not the
    // column this reader thinks it is.
    if cols.len() < 11 {
        return Err(ReadError::at(
            at,
            format!(
                "an InterProScan row is at least 11 columns, this one has {}",
                cols.len()
            ),
        ));
    }
    Ok(cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two proteins, two analyses, and a signature whose description is a
    /// sentence with commas and a slash in it.
    const TSV: &str = "\
P1\tmd5\t500\tPfam\tPF00069\tProtein kinase domain\t11\t275\t1.2e-40\tT\t01-01-2026\tIPR000719\tProtein kinase
P1\tmd5\t500\tPfam\tPF03793\tPASTA domain\t341\t400\t3.0e-10\tT\t01-01-2026\tIPR005543\tPASTA
P1\tmd5\t500\tPANTHER\tPTHR43289\tSERINE/THREONINE KINASE, PUTATIVE\t5\t480\t0.0\tT\t01-01-2026
P2\tmd5\t420\tPfam\tPF00069\tProtein kinase domain\t20\t280\t8.0e-38\tT\t01-01-2026\tIPR000719\tProtein kinase
";

    fn window() -> Region {
        Region::new("protein", 0, 500).unwrap()
    }

    #[test]
    fn a_start_moves_back_one_and_a_stop_does_not() {
        let found = architectures(TSV, &window(), "Pfam").unwrap();
        let kinase = &found.rows[0].features[0];
        assert_eq!((kinase.start, kinase.end), (10, 275));
        // 11 to 275 written 1-based inclusive is 265 residues, one more than
        // the two numbers suggest.
        assert_eq!(kinase.end - kinase.start, 265);
    }

    #[test]
    fn column_one_names_the_row_rather_than_selecting_it() {
        // The figure is the comparison, so every protein is a row on one
        // residue axis, whatever the region happens to be called.
        let found = architectures(TSV, &window(), "Pfam").unwrap();
        assert_eq!(found.rows.len(), 2);
        assert_eq!(found.rows[0].name, "P1");
        assert_eq!(found.rows[1].name, "P2");
        assert_eq!(found.proteins, vec!["P1".to_string(), "P2".to_string()]);
    }

    #[test]
    fn the_length_is_the_column_and_never_the_furthest_domain() {
        // A protein whose last domain ends at 400 may run to 500, and a
        // backbone drawn to 400 says the domain reaches the C terminus.
        let found = architectures(TSV, &window(), "Pfam").unwrap();
        assert_eq!(found.rows[0].length, 500);
        assert_eq!(found.rows[1].length, 420);
        let furthest = found.rows[0].features.iter().map(|f| f.end).max().unwrap();
        assert!(furthest < found.rows[0].length);
    }

    #[test]
    fn a_length_of_nought_is_refused_rather_than_drawn_as_an_empty_row() {
        // It removes the backbone and every domain on it, leaving the protein's
        // name standing over nothing, which is what a protein nobody annotated
        // looks like.
        let text = "P\tmd5\t0\tPfam\tPF1\tkinase\t11\t275\t.\tT\t01-01-2026\n";
        let error = architectures(text, &window(), "Pfam").unwrap_err();
        assert!(error.reason.contains("no length"), "{error}");
    }

    #[test]
    fn a_file_that_disagrees_with_itself_about_a_protein_stops_the_read() {
        let text = "\
P\tmd5\t500\tPfam\tPF1\tkinase\t11\t275\t.\tT\t01-01-2026
P\tmd5\t900\tPfam\tPF2\tPASTA\t341\t400\t.\tT\t01-01-2026
";
        let error = architectures(text, &window(), "Pfam").unwrap_err();
        assert!(error.reason.contains("900 long here"), "{error}");
    }

    #[test]
    fn a_description_with_spaces_in_it_stays_one_field() {
        // The shared splitter falls back to whitespace, and column six is a
        // sentence, so that fallback would make every column after it the
        // wrong one and the start of a domain would be read out of a word.
        let found = architectures(TSV, &window(), "PANTHER").unwrap();
        assert_eq!(
            found.rows[0].features[0].label.as_deref(),
            Some("SERINE/THREONINE KINASE, PUTATIVE")
        );
        assert_eq!(
            (
                found.rows[0].features[0].start,
                found.rows[0].features[0].end
            ),
            (4, 480)
        );
    }

    #[test]
    fn a_name_that_is_a_dot_is_not_a_name() {
        // The track draws whatever string it is given, and calls an unnamed
        // domain `domain`, so a dot left in becomes a domain called `.`.
        let text = "P\tmd5\t500\tPfam\tPF00069\t.\t11\t275\t.\tT\t01-01-2026\n";
        let found = architectures(text, &window(), "Pfam").unwrap();
        assert_eq!(
            found.rows[0].features[0].label.as_deref(),
            Some("PF00069"),
            "the accession is the name when the description is a dot"
        );

        let neither = "P\tmd5\t500\tPfam\t-\t-\t11\t275\t.\tT\t01-01-2026\n";
        let found = architectures(neither, &window(), "Pfam").unwrap();
        assert_eq!(found.rows[0].features[0].label, None);
    }

    #[test]
    fn one_file_is_one_analysis_and_says_which_others_it_holds() {
        let held = analyses(TSV).unwrap();
        assert_eq!(held.get("Pfam"), Some(&3));
        assert_eq!(held.get("PANTHER"), Some(&1));

        let found = architectures(TSV, &window(), "Pfam").unwrap();
        assert_eq!(found.other_analysis, 1);
        assert_eq!(
            found.rows[0].features.len(),
            2,
            "PANTHER reached a Pfam row"
        );
    }

    #[test]
    fn a_protein_named_only_by_another_analysis_is_still_a_row() {
        // Its length was stated and its backbone is a fact about the protein,
        // so it is a row with nothing on it rather than a row nobody drew.
        let text = "P9\tmd5\t300\tPANTHER\tPTHR1\twhatever\t5\t200\t.\tT\t01-01-2026\n";
        let found = architectures(text, &window(), "Pfam").unwrap();
        assert_eq!(found.rows.len(), 1);
        assert_eq!(found.rows[0].length, 300);
        assert!(found.rows[0].features.is_empty());
        assert_eq!(found.other_analysis, 1);
    }

    #[test]
    fn a_domain_outside_the_window_is_counted_rather_than_clipped_in_silence() {
        // Kept, it is drawn to the edge of the axis, which is the idiom for a
        // domain that runs to the C terminus.
        let narrow = architectures(TSV, &Region::new("protein", 0, 300).unwrap(), "Pfam").unwrap();
        assert_eq!(narrow.off_region, 1, "the PASTA domain is past 300");
        assert_eq!(narrow.rows[0].features.len(), 1);
    }

    #[test]
    fn an_inverted_span_is_refused_the_way_an_interval_reader_refuses_one() {
        let text = "P\tmd5\t500\tPfam\tPF1\tkinase\t275\t11\t.\tT\t01-01-2026\n";
        let error = architectures(text, &window(), "Pfam").unwrap_err();
        assert!(error.reason.contains("stop is before start"), "{error}");

        let zero = "P\tmd5\t500\tPfam\tPF1\tkinase\t0\t11\t.\tT\t01-01-2026\n";
        let error = architectures(zero, &window(), "Pfam").unwrap_err();
        assert!(error.reason.contains("counts from 1"), "{error}");
    }

    #[test]
    fn a_line_that_is_not_an_interproscan_row_stops_the_read_and_says_which() {
        let error = architectures("P md5 500 Pfam\n", &window(), "Pfam").unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("holds no tab"), "{error}");

        let short = "P\tmd5\t500\tPfam\tPF1\tkinase\t11\t275\n";
        let error = architectures(short, &window(), "Pfam").unwrap_err();
        assert!(error.reason.contains("at least 11 columns"), "{error}");
    }
}
