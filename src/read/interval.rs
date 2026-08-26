//! Intervals: BED, GFF3 and cytoBand.
//!
//! Gene models, the other annotated spans that go in a feature track, and the
//! cytogenetic bands an ideogram is made of. BED and cytoBand are 0-based and
//! half-open, the convention the rest of the crate counts in, so their
//! coordinates pass straight through. GFF3 is 1-based and inclusive, so its
//! start comes back one and its end stays where it is, having already been one
//! past the last base once the count began at zero. A GFF3 span is one base
//! longer than the BED span written with the same two numbers, and that is the
//! whole of the conversion.
//!
//! # Telling the two apart without being told
//!
//! A `##gff-version` pragma settles it, and most GFF3 files have one. Failing
//! that the answer is in column seven, which GFF3 spends on the strand and a
//! BED of nine or more columns spends on `thickStart`, a number; a row with
//! fewer than seven columns is BED, since a GFF3 row always has nine. The guess
//! is worth making because reading one format as the other moves every feature
//! by a base without failing, and `--format` is there because a guess that can
//! be wrong needs a way to be overruled.
//!
//! # What refuses, and what passes by
//!
//! The sequence name is checked before the rest of a row is understood, so a
//! whole genome annotation costs one comparison per row it does not need, and
//! the FASTA section some GFF3 files carry at the end goes past rather than
//! failing as a run of broken features. A feature that touches no base of the
//! window is dropped here too: it would not be drawn, and it would still take a
//! row in the track's layout, which is what decides how tall the track is.
//!
//! A row that cannot be read stops the file on its line. So does a span whose
//! end is before its start, which is not malformed anywhere else but is here:
//! [`Feature::new`] widens an inverted span into one base at the start and
//! [`Band::new`] clamps it to nothing, so a gene would be drawn a whole
//! interval from where either coordinate put it, or a band would vanish while
//! its end still set the length of the chromosome.

use crate::{Band, Feature, Region, Stain, Strand};

use super::Format;
use super::{columns, lines, number, ReadError};

/// Which of the two interval formats a file turned out to be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Flavour {
    /// `chrom start end [name] [score] [strand]`, 0-based half-open.
    Bed,
    /// Nine columns, 1-based inclusive, attributes in the ninth.
    Gff3,
}

/// Reads gene models and other intervals.
///
/// BED is 0-based half-open and passes straight through. GFF3 is 1-based
/// inclusive, so its start becomes `start - 1` and its end is unchanged.
///
/// The two are told apart by column seven, which is a strand in GFF3 and a
/// number in a BED with that many columns, and by the `##gff-version` line.
/// `format` overrides the guess.
///
/// A name comes from BED column four, or from `Name=`, failing that `gene=`,
/// failing that `ID=` in the GFF3 attributes. A strand comes from BED column
/// six or GFF3 column seven.
///
/// Rows on another sequence than `region.seq()` are skipped.
pub fn features(
    text: &str,
    region: &Region,
    format: Option<Format>,
) -> Result<Vec<Feature>, ReadError> {
    let flavour = flavour(text, format);
    let mut features = Vec::new();

    for (at, line) in lines(text) {
        let cols = columns(line);
        // The sequence is checked before the rest of the line is understood, so
        // that a whole genome annotation costs one comparison per row it does
        // not need. It also means the FASTA section some GFF3 files carry at the
        // end goes past without being read as a broken feature.
        if cols.first().copied().unwrap_or_default() != region.seq() {
            continue;
        }
        let feature = match flavour {
            Flavour::Bed => bed(&cols, at)?,
            Flavour::Gff3 => gff3(&cols, at)?,
        };
        // A feature the window does not touch is not drawn, and it would still
        // take a row in the track's layout, which is what decides how tall the
        // track is. So it is dropped here rather than carried.
        if feature.end <= region.start() || feature.start >= region.end() {
            continue;
        }
        features.push(feature);
    }

    Ok(features)
}

/// Reads a cytoBand table: `chrom start end name stain`, 0-based half-open.
///
/// Returns the length of the sequence, taken as the highest end seen, and the
/// bands. The stain words are the UCSC ones: `gneg`, `gpos25`, `gpos50`,
/// `gpos75`, `gpos100`, `acen`, `gvar` and `stalk`.
pub fn cytoband(text: &str, sequence: &str) -> Result<(u64, Vec<Band>), ReadError> {
    let mut bands = Vec::new();
    let mut length = 0;

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.first().copied().unwrap_or_default() != sequence {
            continue;
        }
        if cols.len() < 3 {
            return Err(ReadError::at(
                at,
                "expected at least 3 columns: chrom start end [name] [stain]",
            ));
        }
        // cytoBand is 0-based and half-open, so both ends pass through.
        let start: u64 = number(cols[1].trim(), "start", at)?;
        let end: u64 = number(cols[2].trim(), "end", at)?;
        // `Band::new` clamps an inverted span to zero width, so the band would
        // vanish from the ideogram while its end still set the chromosome
        // length. Both of those are silent, so the line is refused instead.
        if end < start {
            return Err(ReadError::at(at, "end is before start"));
        }
        // An unknown or missing stain is the palest band rather than a guess,
        // which is what Stain::from_name already does.
        let mut band = Band::new(
            start,
            end,
            Stain::from_name(cols.get(4).copied().unwrap_or("")),
        );
        if let Some(name) = label(cols.get(3).copied()) {
            band = band.name(name);
        }
        length = length.max(band.end);
        bands.push(band);
    }

    Ok((length, bands))
}

/// Decides whether a file is BED or GFF3.
///
/// In order: an explicit `--format` wins, then a `##gff-version` line, then
/// column seven of the first data row, which GFF3 spends on the strand and a
/// BED of nine or more columns spends on `thickStart`, a number. A row with
/// fewer than seven columns is BED, since a GFF3 row always has nine.
pub(crate) fn flavour(text: &str, format: Option<Format>) -> Flavour {
    match format {
        Some(Format::Gff3) => return Flavour::Gff3,
        Some(Format::Bed) => return Flavour::Bed,
        // The other words name signal formats, so they say nothing about which
        // of these two an interval file is and the guess still runs.
        _ => {}
    }
    if text
        .lines()
        .any(|line| line.trim_start().starts_with("##gff-version"))
    {
        return Flavour::Gff3;
    }
    let Some((_, first)) = lines(text).next() else {
        return Flavour::Bed;
    };
    let seventh = columns(first)
        .get(6)
        .map(|column| column.trim().to_string());
    match seventh.as_deref() {
        Some("+" | "-" | "." | "?") => Flavour::Gff3,
        _ => Flavour::Bed,
    }
}

/// One BED row, whose coordinates are already the ones the crate uses.
pub(crate) fn bed(cols: &[&str], at: usize) -> Result<Feature, ReadError> {
    if cols.len() < 3 {
        return Err(ReadError::at(
            at,
            "expected at least 3 columns: chrom start end",
        ));
    }
    let start: u64 = number(cols[1].trim(), "start", at)?;
    let end: u64 = number(cols[2].trim(), "end", at)?;
    // `Feature::new` widens an inverted span into one base at the start, which
    // would draw a gene a whole interval away from where the file put it.
    if end < start {
        return Err(ReadError::at(at, "end is before start"));
    }

    let mut feature = Feature::new(start, end);
    if let Some(name) = label(cols.get(3).copied()) {
        feature = feature.name(name);
    }
    if let Some(strand) = cols.get(5) {
        feature = feature.strand(strand_of(strand));
    }
    Ok(feature)
}

/// One GFF3 row, whose start counts from one and whose end is included.
pub(crate) fn gff3(cols: &[&str], at: usize) -> Result<Feature, ReadError> {
    // The last four columns are not needed to place a feature, so a file that
    // stops short of nine still reads. The first five are.
    if cols.len() < 5 {
        return Err(ReadError::at(
            at,
            "expected 9 columns: seqid source type start end score strand phase attributes",
        ));
    }
    let start: u64 = number(cols[3].trim(), "start", at)?;
    let end: u64 = number(cols[4].trim(), "end", at)?;
    if start == 0 {
        return Err(ReadError::at(
            at,
            "GFF3 counts from 1, so 0 is not a start position",
        ));
    }
    // 1-based inclusive to 0-based half-open: the start moves back one, the end
    // stays where it is because it was already one past the last base once the
    // count started at zero.
    let mut feature = Feature::new(start - 1, end);
    if let Some(name) = cols.get(8).and_then(|attributes| gff3_name(attributes)) {
        feature = feature.name(name);
    }
    if let Some(strand) = cols.get(6) {
        feature = feature.strand(strand_of(strand));
    }
    Ok(feature)
}

/// The name a GFF3 record goes by.
///
/// `Name=` is what a browser shows, `gene=` is what the gene is called when the
/// record has no display name, and `ID=` is there in every record and so is the
/// last resort rather than the first choice.
fn gff3_name(attributes: &str) -> Option<String> {
    ["Name", "gene", "ID"]
        .into_iter()
        .find_map(|key| attribute(attributes, key))
}

/// One `key=value` out of the ninth column, exactly as it was written.
///
/// Undecoded on purpose, for a value that is a list. The column spends `,` on
/// its own list syntax, so a name holding a comma arrives as `%2C`, and
/// decoding the whole value before splitting it turns one name into two
/// without anything failing. A caller reading a list splits first and calls
/// [`percent_decode`] on each piece.
pub(crate) fn raw_attribute<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    attributes.split(';').find_map(|pair| {
        let (found, value) = pair.trim().split_once('=')?;
        (found.trim() == key).then(|| value.trim())
    })
}

/// One `key=value` out of the ninth column, decoded.
fn attribute(attributes: &str, key: &str) -> Option<String> {
    label(Some(&percent_decode(raw_attribute(attributes, key)?)))
}

/// Decodes the `%XX` escapes of a GFF3 attribute value.
///
/// The column spends `;`, `=` and `,` on its own syntax, so a value holding one
/// of them arrives escaped, and `%20` for a space is just as common. Anything
/// that is not a complete escape is left as it was written, since a stray `%`
/// in a name is more likely than a truncated one.
pub(crate) fn percent_decode(value: &str) -> String {
    if !value.contains('%') {
        return value.to_string();
    }
    let bytes = value.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 3 <= bytes.len() {
            let decoded = std::str::from_utf8(&bytes[i + 1..i + 3])
                .ok()
                .and_then(|hex| u8::from_str_radix(hex, 16).ok());
            if let Some(byte) = decoded {
                out.push(byte);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// A name column, with the `.` every one of these formats uses for "none" read
/// as none rather than as a feature called `.`.
pub(crate) fn label(field: Option<&str>) -> Option<String> {
    let text = field?.trim();
    if text.is_empty() || text == "." {
        None
    } else {
        Some(text.to_string())
    }
}

/// The strand column of either format.
fn strand_of(field: &str) -> Strand {
    field
        .trim()
        .chars()
        .next()
        .map_or(Strand::Unknown, Strand::from_symbol)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Arabidopsis genes, as a BED with a UCSC track line on top.
    const BED: &str = "\
track name=genes description=\"TAIR10\"
Chr1\t3630\t5899\tAT1G01010\t0\t+
Chr1\t6787\t9130\tAT1G01020\t0\t-
Chr2\t3000\t4000\tAT2G01010\t0\t+
";

    /// A GFF3 of the rifampicin resistance locus, pragma and all.
    const GFF3: &str = "\
##gff-version 3
#!genome-build H37Rv
NC_000962.3\tRefSeq\tgene\t759807\t763325\t.\t+\t.\tID=gene-Rv0667;Name=rpoB
NC_000962.3\tRefSeq\tgene\t763370\t767320\t.\t+\t.\tID=gene-Rv0668;Name=rpoC
";

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    fn read(text: &str, locus: &str, format: Option<Format>) -> Vec<Feature> {
        features(text, &region(locus), format).unwrap()
    }

    #[test]
    fn bed_coordinates_pass_straight_through() {
        let genes = read(BED, "Chr1:1-10000", None);
        assert_eq!(genes.len(), 2);
        // The BED row says 3630 5899 and the feature says the same, because both
        // count from zero and leave the end out.
        assert_eq!(genes[0].start, 3630);
        assert_eq!(genes[0].end, 5899);
        assert_eq!(genes[0].name.as_deref(), Some("AT1G01010"));
        assert_eq!(genes[0].strand, Strand::Forward);
        assert_eq!(genes[1].strand, Strand::Reverse);
    }

    #[test]
    fn a_gff3_start_moves_back_one_and_its_end_does_not() {
        // The row is 759807..763325, 1-based and inclusive. Both ends name the
        // same two bases afterwards, counted from zero with the end left out.
        let genes = read(GFF3, "NC_000962.3:759000-764000", None);
        assert_eq!(genes[0].start, 759_806);
        assert_eq!(genes[0].end, 763_325);
        assert_eq!(genes[0].name.as_deref(), Some("rpoB"));
        assert_eq!(genes[0].strand, Strand::Forward);
    }

    #[test]
    fn the_pinned_gff3_conversion_is_start_minus_one_and_end_unchanged() {
        let text = "##gff-version 3\nchrX\tsource\tgene\t100\t200\t.\t+\t.\tID=x\n";
        let feature = &read(text, "chrX:1-1000", None)[0];
        assert_eq!((feature.start, feature.end), (99, 200));
        // A hundred and one bases, which is what 100..200 inclusive holds.
        assert_eq!(feature.len(), 101);
    }

    #[test]
    fn column_seven_tells_the_two_apart_without_a_pragma() {
        // A GFF3 with no pragma: column seven is a strand.
        let gff = "Chr1\tphytozome\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010\n";
        assert_eq!(flavour(gff, None), Flavour::Gff3);
        assert_eq!(read(gff, "Chr1:1-10000", None)[0].start, 3630);

        // A BED9: column seven is thickStart, a number.
        let bed = "chr7\t127471196\t127472363\tPos1\t0\t+\t127471196\t127472363\t255,0,0\n";
        assert_eq!(flavour(bed, None), Flavour::Bed);
        assert_eq!(
            read(bed, "chr7:127471000-127473000", None)[0].start,
            127_471_196
        );
    }

    #[test]
    fn a_row_with_fewer_than_seven_columns_is_bed() {
        let bed = "NC_045512.2\t265\t13468\n";
        assert_eq!(flavour(bed, None), Flavour::Bed);
        let feature = &read(bed, "NC_045512.2:1-30000", None)[0];
        assert_eq!((feature.start, feature.end), (265, 13_468));
        assert_eq!(feature.name, None);
        assert_eq!(feature.strand, Strand::Unknown);
    }

    #[test]
    fn the_format_flag_wins_over_the_guess() {
        // Six columns, so the guess says BED and the second column is not a
        // number. Told it is GFF3, the same row reads as one.
        let text = "NC_045512.2\tfeature\tgene\t266\t13468\t.\n";
        assert_eq!(flavour(text, None), Flavour::Bed);
        assert!(features(text, &region("NC_045512.2:1-30000"), None).is_err());

        let feature = &read(text, "NC_045512.2:1-30000", Some(Format::Gff3))[0];
        assert_eq!((feature.start, feature.end), (265, 13_468));
    }

    #[test]
    fn a_pragma_is_enough_on_its_own() {
        assert_eq!(flavour(GFF3, None), Flavour::Gff3);
        assert_eq!(flavour(BED, None), Flavour::Bed);
        assert_eq!(flavour("", None), Flavour::Bed);
    }

    #[test]
    fn the_name_is_name_then_gene_then_id() {
        let all = "chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene0;gene=katG;Name=catalase\n";
        assert_eq!(
            read(all, "chr1:1-100", None)[0].name.as_deref(),
            Some("catalase")
        );

        let no_display = "chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene0;gene=katG\n";
        assert_eq!(
            read(no_display, "chr1:1-100", None)[0].name.as_deref(),
            Some("katG")
        );

        let bare = "chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene0\n";
        assert_eq!(
            read(bare, "chr1:1-100", None)[0].name.as_deref(),
            Some("gene0")
        );

        let none = "chr1\t.\tgene\t1\t9\t.\t+\t.\tParent=x\n";
        assert_eq!(read(none, "chr1:1-100", None)[0].name, None);
    }

    #[test]
    fn attribute_values_are_percent_decoded() {
        let text = "\
##gff-version 3
NC_002516.2\t.\tgene\t1\t9\t.\t+\t.\tName=chromosomal%20replication%2C%20initiator
";
        assert_eq!(
            read(text, "NC_002516.2:1-100", None)[0].name.as_deref(),
            Some("chromosomal replication, initiator")
        );
        // A per cent sign that starts nothing is a per cent sign.
        assert_eq!(percent_decode("100%"), "100%");
        assert_eq!(percent_decode("50%2Fhour"), "50/hour");
        assert_eq!(percent_decode("%zz"), "%zz");
    }

    #[test]
    fn rows_on_another_sequence_are_skipped() {
        // The Chr2 gene sits inside the coordinates on display and is still not
        // in the figure, because it is not on the sequence being drawn.
        let genes = read(BED, "Chr1:1-10000", None);
        assert!(genes.iter().all(|g| g.name.as_deref() != Some("AT2G01010")));
        let other = read(BED, "Chr2:1-10000", None);
        assert_eq!(other.len(), 1);
        assert_eq!(other[0].name.as_deref(), Some("AT2G01010"));
    }

    #[test]
    fn rows_outside_the_region_are_skipped() {
        let genes = read(BED, "Chr1:5000-6000", None);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name.as_deref(), Some("AT1G01010"));

        // A feature that ends exactly where the region starts touches no base
        // of it. Chr1:3631-4000 is 3630..4000, and 0..3630 stops one short.
        let touching = "Chr1\t0\t3630\tupstream\n";
        assert!(read(touching, "Chr1:3631-4000", None).is_empty());
    }

    #[test]
    fn a_malformed_row_carries_its_line_number() {
        let text = "\
##gff-version 3
NC_000962.3\tRefSeq\tgene\t759807\t763325\t.\t+\t.\tID=a
NC_000962.3\tRefSeq\tgene\tninety\t763325\t.\t+\t.\tID=b
";
        let error = features(text, &region("NC_000962.3:1-800000"), None).unwrap_err();
        assert_eq!(error.line, 3);
        assert_eq!(
            error.to_string(),
            "line 3: start is not a number: \"ninety\""
        );
    }

    #[test]
    fn a_row_too_short_to_be_a_feature_says_what_was_expected() {
        let error = features("Chr1\t3630\n", &region("Chr1:1-10000"), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("chrom start end"), "{error}");

        let text = "##gff-version 3\nChr1\tRefSeq\tgene\t100\n";
        let error = features(text, &region("Chr1:1-10000"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("seqid source type"), "{error}");
    }

    #[test]
    fn a_gff3_start_of_zero_is_an_error_and_not_an_underflow() {
        let text = "##gff-version 3\nchr1\t.\tgene\t0\t200\t.\t+\t.\tID=x\n";
        let error = features(text, &region("chr1:1-1000"), None).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("counts from 1"), "{error}");
    }

    #[test]
    fn a_file_with_nothing_in_it_is_not_an_error() {
        assert!(read("", "chr1:1-1000", None).is_empty());
        assert!(read("# only a comment\n\n", "chr1:1-1000", None).is_empty());
    }

    #[test]
    fn a_gff3_carrying_its_sequence_at_the_end_still_reads() {
        // Assemblers write the FASTA after the annotation. Those lines are on no
        // sequence the region names, so they go past rather than failing.
        let text = "\
##gff-version 3
NC_045512.2\tprokka\tgene\t266\t13468\t.\t+\t.\tID=ORF1ab;Name=ORF1ab
##FASTA
>NC_045512.2
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
";
        let genes = read(text, "NC_045512.2:1-30000", None);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].start, 265);
    }

    /// Human chromosome 21, with a neighbour in the file to be ignored.
    const CYTOBAND: &str = "\
chr21\t0\t2800000\tp13\tgvar
chr21\t2800000\t6970000\tp12\tstalk
chr21\t10900000\t12000000\tp11.1\tacen
chr21\t12000000\t46709983\tq22.3\tgneg
chr20\t0\t64444167\tp13\tgneg
";

    #[test]
    fn cytoband_coordinates_pass_straight_through() {
        let (length, bands) = cytoband(CYTOBAND, "chr21").unwrap();
        assert_eq!(bands.len(), 4);
        // The row says 0 2800000 and so does the band: cytoBand is BED.
        assert_eq!((bands[0].start, bands[0].end), (0, 2_800_000));
        assert_eq!(bands[0].name.as_deref(), Some("p13"));
        // The length is the far end of the last band on this sequence, and not
        // the far end of the longer chromosome further down the file.
        assert_eq!(length, 46_709_983);
    }

    #[test]
    fn the_stain_words_are_the_ucsc_ones() {
        let (_, bands) = cytoband(CYTOBAND, "chr21").unwrap();
        assert_eq!(bands[0].stain, Stain::Gvar);
        assert_eq!(bands[1].stain, Stain::Stalk);
        assert_eq!(bands[2].stain, Stain::Acen);
        assert_eq!(bands[3].stain, Stain::Gneg);
    }

    #[test]
    fn a_sequence_the_table_does_not_hold_gives_no_bands() {
        let (length, bands) = cytoband(CYTOBAND, "chr22").unwrap();
        assert_eq!(length, 0);
        assert!(bands.is_empty());
    }

    #[test]
    fn a_cytoband_row_may_stop_after_the_coordinates() {
        // Mouse, from a table with no stain column at all.
        let (length, bands) = cytoband("chr19\t0\t61420004\n", "chr19").unwrap();
        assert_eq!(length, 61_420_004);
        assert_eq!(bands[0].stain, Stain::Gneg);
        assert_eq!(bands[0].name, None);
    }

    #[test]
    fn a_malformed_cytoband_row_carries_its_line_number() {
        let text = "chr21\t0\t2800000\tp13\tgvar\nchr21\t2800000\tsix\tp12\tstalk\n";
        let error = cytoband(text, "chr21").unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.to_string().contains("end is not a number"), "{error}");

        let short = cytoband("chr21\t0\n", "chr21").unwrap_err();
        assert_eq!(short.line, 1);
        assert!(short.to_string().contains("chrom start end"), "{short}");
    }
}
