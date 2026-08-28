//! Segmented copy number, read by column name off a required header.
//!
//! ```text
//! chromosome  start      end        gene  log2   cn  cn1  cn2
//! chr8        127200000  127740000  MYC   1.86   7   5    2
//! chr8        127740000  129100000  -     -0.02  2   1    1
//! chr17       7565000    7590000    TP53  -1.04  1   1    0
//! ```
//!
//! # The header is not optional, and it is not a courtesy
//!
//! Everywhere else in this module a format is told by its shape: a BED row has
//! a number in column two, a GFF3 row has nine columns. Every segment table has
//! a chromosome, two coordinates and some numbers, and the writers disagree
//! about the order, the count and the names. Sniffing by column count would
//! read `nMajor` where `nMinor` was written, which is a figure claiming loss of
//! heterozygosity in the arms that kept it, drawn confidently.
//!
//! So the columns are found by name, and a file whose header names none of the
//! shapes below is refused rather than guessed at.
//!
//! # Which files count from one
//!
//! CNVkit's `.cns` is BED-like: 0-based and half-open, passed straight through.
//! ASCAT's segment table and the `.seg` file IGV and GISTIC2 read are 1-based
//! and inclusive, so one comes off the start and the end is unchanged. The
//! header is what decides which, since the two spell the same two coordinates
//! by different names.
//!
//! # A log ratio is not a copy number until somebody says what balanced is
//!
//! A `.seg` file carries `seg.mean`, a log2 ratio against a reference, and a
//! `.cns` carries `log2` for the same reason. Turning either into copies is
//! `ploidy * 2^log2`, and the ploidy is not in the file. It is the caller's,
//! which is why it is an argument here rather than a default: two is right for
//! a human autosome and wrong for everything else this crate is used on.
//!
//! A file that carries a called integer copy number is read from that instead,
//! since it is what the caller concluded rather than what this arithmetic
//! would infer.

use crate::{CopyNumber, CopyNumberSegment, Region};

use super::{columns, lines, ReadError};

/// The segments a table holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Segmentation {
    /// The segments, in the order the file listed them.
    pub segments: Vec<CopyNumberSegment>,
    /// Segment rows in the file, before any filter.
    pub records: usize,
    /// Rows whose copy number is a stated no call.
    ///
    /// Separate from a row that would not parse: a caller writing `NA` has
    /// said something, and a figure short of those segments is short of them
    /// for a reason worth reporting.
    pub no_call: usize,
    /// Rows naming another sequence, where the file named more than one.
    pub other_sequence: usize,
    /// Rows on this sequence that touch no base of the window.
    pub off_region: usize,
    /// The samples the file named, in the order they were first seen.
    pub samples: Vec<String>,
}

/// Which table this is, and therefore where its numbers are.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shape {
    /// CNVkit `.cns`, 0-based half-open.
    Cns,
    /// ASCAT segments, 1-based inclusive.
    Ascat,
    /// The `.seg` file, 1-based inclusive.
    Seg,
}

/// Where each thing this reader needs lives in a row.
#[derive(Debug, Clone, Copy)]
struct Layout {
    shape: Shape,
    sequence: usize,
    start: usize,
    end: usize,
    total: Option<usize>,
    major: Option<usize>,
    minor: Option<usize>,
    ratio: Option<usize>,
    sample: Option<usize>,
}

/// Reads a segment table, converting log ratios at `ploidy` copies.
///
/// `sample` picks one out of a table holding several. A table naming more than
/// one sample and no choice made is not read: the segments of two samples in
/// one band are two step functions drawn over each other, which looks like one
/// sample with a great many breakpoints.
///
/// # Errors
///
/// Returns the first row that is not the shape its header promised, a
/// coordinate that is not a number, a 1-based coordinate of zero, or a span
/// whose end is before its start. A header naming none of the known shapes is
/// refused before any row is read.
pub fn copy_numbers(
    text: &str,
    region: &Region,
    ploidy: f64,
    sample: Option<&str>,
) -> Result<Segmentation, ReadError> {
    let mut rows = lines(text);
    let (at, head) = rows.next().ok_or_else(|| {
        ReadError::whole("a segment table begins with a header and this file is empty")
    })?;
    let layout = layout(&columns(head), at)?;

    let mut found = Segmentation::default();
    let mut wanted = 0usize;

    for (at, line) in rows {
        let cols = columns(line);
        let width = [
            layout.sequence,
            layout.start,
            layout.end,
            layout.total.unwrap_or(0),
            layout.major.unwrap_or(0),
            layout.minor.unwrap_or(0),
            layout.ratio.unwrap_or(0),
            layout.sample.unwrap_or(0),
        ]
        .into_iter()
        .max()
        .unwrap_or(0);
        if cols.len() <= width {
            return Err(ReadError::at(
                at,
                format!(
                    "the header names a column at position {} and this row has {}",
                    width + 1,
                    cols.len()
                ),
            ));
        }
        found.records += 1;

        if let Some(index) = layout.sample {
            let named = cols[index].to_string();
            if !found.samples.contains(&named) {
                found.samples.push(named.clone());
            }
            if let Some(wanted_sample) = sample {
                if named != wanted_sample {
                    continue;
                }
            }
            wanted += 1;
        }

        if cols[layout.sequence] != region.seq() {
            found.other_sequence += 1;
            continue;
        }

        let raw_start: u64 = super::number(cols[layout.start], "start", at)?;
        let raw_end: u64 = super::number(cols[layout.end], "end", at)?;
        // Compared as the file wrote them, before either is converted. Checked
        // afterwards, a 1-based row whose end is one below its start survives
        // the subtraction and becomes a segment covering no bases.
        if raw_end < raw_start {
            return Err(ReadError::at(at, "end is before start"));
        }
        let end = raw_end;
        let start = if layout.shape == Shape::Cns {
            raw_start
        } else {
            if raw_start == 0 {
                return Err(ReadError::at(
                    at,
                    "this table counts from 1, so 0 is not a start",
                ));
            }
            raw_start - 1
        };
        if end <= region.start() || start >= region.end() {
            found.off_region += 1;
            continue;
        }

        let Some(copy) = call(&cols, &layout, ploidy) else {
            found.no_call += 1;
            continue;
        };
        found.segments.push(CopyNumberSegment { start, end, copy });
    }

    // Two samples in one band are two step functions drawn over each other,
    // which reads as one sample with a great many breakpoints.
    if sample.is_none() && found.samples.len() > 1 {
        return Err(ReadError::whole(format!(
            "this table holds {} samples and --sample says which to draw",
            found.samples.len()
        )));
    }
    let _ = wanted;

    Ok(found)
}

/// The samples a table names, for a caller that has to choose between them.
pub fn samples(text: &str) -> Result<Vec<String>, ReadError> {
    let mut rows = lines(text);
    let (at, head) = rows.next().ok_or_else(|| {
        ReadError::whole("a segment table begins with a header and this file is empty")
    })?;
    let layout = layout(&columns(head), at)?;
    let Some(index) = layout.sample else {
        return Ok(Vec::new());
    };

    let mut named: Vec<String> = Vec::new();
    for (_, line) in rows {
        let cols = columns(line);
        if let Some(sample) = cols.get(index) {
            let sample = sample.to_string();
            if !named.contains(&sample) {
                named.push(sample);
            }
        }
    }
    Ok(named)
}

/// Finds each column this reader needs, or says the header is none of the
/// shapes it knows.
fn layout(head: &[&str], at: usize) -> Result<Layout, ReadError> {
    let lower: Vec<String> = head.iter().map(|name| name.trim().to_lowercase()).collect();
    let find = |wanted: &[&str]| {
        lower
            .iter()
            .position(|name| wanted.contains(&name.as_str()))
    };

    let sequence = find(&["chromosome", "chrom", "chr", "seqnames"]);
    let sample = find(&["sample", "id", "sampleid", "sample_id", "name"]);

    if let (Some(sequence), Some(start), Some(end)) = (
        sequence,
        find(&["start", "startpos", "start_pos"]),
        find(&["end", "endpos", "end_pos"]),
    ) {
        let major = find(&["cn1", "nmajor", "nmaj", "major"]);
        let minor = find(&["cn2", "nminor", "nmin", "minor"]);
        let total = find(&["cn", "total_cn", "copy_number", "copies"]);
        let ratio = find(&["log2", "log2ratio", "logr"]);
        // ASCAT counts from one and CNVkit does not, and the two are told
        // apart by the names they use for the same two coordinates. The list
        // has to be the same list `find` accepted, or a table spelled one of
        // the aliases is read in the other convention and every segment lands
        // a base to the right with the guard against a start of nought
        // stepped over on the way.
        let one_based = ["startpos", "start_pos", "endpos", "end_pos"];
        let shape = if lower.iter().any(|name| one_based.contains(&name.as_str())) {
            Shape::Ascat
        } else {
            Shape::Cns
        };
        if total.is_some() || ratio.is_some() || (major.is_some() && minor.is_some()) {
            return Ok(Layout {
                shape,
                sequence,
                start,
                end,
                total,
                major,
                minor,
                ratio,
                sample,
            });
        }
    }

    if let (Some(sequence), Some(start), Some(end), Some(ratio)) = (
        sequence,
        find(&["loc.start", "loc_start"]),
        find(&["loc.end", "loc_end"]),
        find(&["seg.mean", "seg_mean", "segmean"]),
    ) {
        return Ok(Layout {
            shape: Shape::Seg,
            sequence,
            start,
            end,
            total: None,
            major: None,
            minor: None,
            ratio: Some(ratio),
            sample,
        });
    }

    Err(ReadError::at(
        at,
        format!(
            "this header names none of the segment tables this reads: {}",
            head.join(", ")
        ),
    ))
}

/// What a row calls, or nothing where it calls nothing.
fn call(cols: &[&str], layout: &Layout, ploidy: f64) -> Option<CopyNumber> {
    // The allele split first, because it is the stronger statement and the one
    // a total alone cannot be turned back into.
    if let (Some(major), Some(minor)) = (layout.major, layout.minor) {
        if let (Some(a), Some(b)) = (number(cols.get(major)), number(cols.get(minor))) {
            return Some(CopyNumber::allelic(a, b));
        }
    }
    if let Some(total) = layout.total {
        if let Some(copies) = number(cols.get(total)) {
            return Some(CopyNumber::Total(copies));
        }
    }
    // Last, because it is inferred rather than concluded: a caller that wrote a
    // copy number wrote what it decided, and this is arithmetic on a ratio.
    if let Some(ratio) = layout.ratio {
        if let Some(log2) = number(cols.get(ratio)) {
            return Some(CopyNumber::Total(ploidy * 2f64.powf(log2)));
        }
    }
    None
}

/// One number, or nothing where the field says there is none.
///
/// The three spellings the rest of this module reads as nothing, plus a value
/// that parses and is not a number: an `NaN` copy number is a copy number
/// nobody called, and drawn it would be a segment at whatever the clamp
/// decided.
fn number(field: Option<&&str>) -> Option<f64> {
    let field = field?.trim();
    if field.is_empty() || field == "." || field == "NA" || field == "-" {
        return None;
    }
    let value = field.parse::<f64>().ok()?;
    value.is_finite().then_some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region() -> Region {
        Region::new("chr8", 0, 10_000_000).unwrap()
    }

    const CNS: &str = "\
chromosome\tstart\tend\tgene\tlog2\tcn\tcn1\tcn2
chr8\t0\t2000000\t-\t-0.02\t2\t1\t1
chr8\t2000000\t3500000\t-\t-1.04\t1\t1\t0
chr8\t7000000\t7400000\t-\tNA\tNA\tNA\tNA
chr17\t0\t100\t-\t0.5\t3\t2\t1
";

    #[test]
    fn a_cns_is_read_where_it_lies_and_counted_where_it_does_not() {
        let found = copy_numbers(CNS, &region(), 2.0, None).expect("segments");
        assert_eq!(found.records, 4);
        assert_eq!(found.segments.len(), 2);
        assert_eq!(found.no_call, 1);
        assert_eq!(found.other_sequence, 1);
        // BED-like, so the start passes straight through.
        assert_eq!(found.segments[0].start, 0);
        assert_eq!(found.segments[0].end, 2_000_000);
    }

    #[test]
    fn the_allele_split_is_read_before_the_total_and_the_total_before_the_ratio() {
        // A caller that wrote both wrote the split on purpose, and a total is
        // what it concluded rather than what this arithmetic would infer.
        let found = copy_numbers(CNS, &region(), 2.0, None).expect("segments");
        assert_eq!(found.segments[1].copy.minor(), Some(0.0));
        assert_eq!(found.segments[1].copy.total(), 1.0);
        assert_eq!(found.segments[1].loh(), Some(true));
    }

    #[test]
    fn a_seg_file_counts_from_one_and_its_ratio_becomes_copies_at_the_ploidy() {
        let text = "\
ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean
S1\tchr8\t1\t2000000\t1043\t0.0
S1\tchr8\t2000001\t3000000\t998\t1.0
";
        let found = copy_numbers(text, &region(), 2.0, None).expect("segments");
        assert_eq!(found.segments[0].start, 0, "1-based start not taken down");
        assert_eq!(found.segments[0].end, 2_000_000);
        assert_eq!(found.segments[0].copy.total(), 2.0);
        // One doubling above the reference, at two copies, is four.
        assert_eq!(found.segments[1].copy.total(), 4.0);
        // And the ploidy is what says so: the same file at one copy is two.
        let haploid = copy_numbers(text, &region(), 1.0, None).expect("segments");
        assert_eq!(haploid.segments[1].copy.total(), 2.0);
    }

    #[test]
    fn a_table_of_several_samples_is_not_drawn_until_one_is_chosen() {
        // Two step functions in one band read as one sample with a great many
        // breakpoints, which is a figure nobody would question.
        let text = "\
ID\tchrom\tloc.start\tloc.end\tseg.mean
S1\tchr8\t1\t2000000\t0.0
S2\tchr8\t1\t2000000\t1.0
";
        let error = copy_numbers(text, &region(), 2.0, None).unwrap_err();
        assert!(error.reason.contains("2 samples"), "{error}");

        let found = copy_numbers(text, &region(), 2.0, Some("S2")).expect("segments");
        assert_eq!(found.segments.len(), 1);
        assert_eq!(found.segments[0].copy.total(), 4.0);
        assert_eq!(samples(text).expect("samples"), ["S1", "S2"]);
    }

    #[test]
    fn a_header_naming_none_of_the_known_tables_is_refused() {
        // Rather than guessed at. Reading nMajor where nMinor was written is a
        // figure claiming lost heterozygosity in the arms that kept it.
        let error = copy_numbers("a\tb\tc\n1\t2\t3\n", &region(), 2.0, None).unwrap_err();
        assert!(
            error.reason.contains("names none of the segment tables"),
            "{error}"
        );
    }

    #[test]
    fn ascat_counts_from_one_and_cnvkit_does_not() {
        let ascat = "\
sample\tchr\tstartpos\tendpos\tnMajor\tnMinor
T1\tchr8\t1\t2000000\t2\t0
";
        let found = copy_numbers(ascat, &region(), 2.0, None).expect("segments");
        assert_eq!(found.segments[0].start, 0);
        assert_eq!(found.segments[0].copy.minor(), Some(0.0));

        let zero = "\
sample\tchr\tstartpos\tendpos\tnMajor\tnMinor
T1\tchr8\t0\t2000000\t2\t0
";
        let error = copy_numbers(zero, &region(), 2.0, None).unwrap_err();
        assert!(error.reason.contains("counts from 1"), "{error}");
    }

    #[test]
    fn a_stated_no_call_is_counted_and_never_becomes_a_number() {
        for spelling in ["NA", ".", "", "-", "NaN"] {
            let text = format!("chromosome\tstart\tend\tcn\nchr8\t0\t100\t{spelling}\n");
            let found = copy_numbers(&text, &region(), 2.0, None).expect("segments");
            assert_eq!(found.segments.len(), 0, "{spelling:?} became a copy number");
            assert_eq!(found.no_call, 1, "{spelling:?} was not counted");
        }
    }

    #[test]
    fn a_row_shorter_than_its_header_is_refused_by_line() {
        let error = copy_numbers(
            "chromosome\tstart\tend\tcn\nchr8\t0\n",
            &region(),
            2.0,
            None,
        )
        .unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.reason.contains("this row has 2"), "{error}");
    }

    #[test]
    fn an_end_before_its_start_is_refused() {
        let error = copy_numbers(
            "chromosome\tstart\tend\tcn\nchr8\t500\t100\t2\n",
            &region(),
            2.0,
            None,
        )
        .unwrap_err();
        assert!(error.reason.contains("end is before start"), "{error}");
    }
}
