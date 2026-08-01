//! Coordinate audit: one known base per format, start and end checked apart.
//!
//! Every fixture below covers exactly one base, or a span whose length is known
//! from the specification, so that a shift of one shows up as a failed
//! assertion rather than as a figure that looks fine.
//!
//! The base everything aims at is 1-based position 100 on `chr1`, which is
//! 0-based 99.

#![cfg(test)]

use karyon::{Region, Strand};

use super::{align, interval, point, signal, table};

/// 0-based 99, which every fixture in this file is written to land on.
const TARGET: u64 = 99;

fn region(locus: &str) -> Region {
    Region::parse(locus).unwrap()
}

// ---------------------------------------------------------------------------
// One base, per format.
// ---------------------------------------------------------------------------

#[test]
fn audit_bedgraph_one_base() {
    // bedGraph is 0-based half-open: 99..100 is the single base 1-based 100.
    let pairs = signal::dense("chr1\t99\t100\t7\n", &region("chr1:1-200"), None).unwrap();
    assert_eq!(pairs, vec![(TARGET, 7.0)]);
}

#[test]
fn audit_depth_one_base() {
    // samtools depth is 1-based: position 100 is 0-based 99.
    let pairs = signal::dense("chr1\t100\t7\n", &region("chr1:1-200"), None).unwrap();
    assert_eq!(pairs, vec![(TARGET, 7.0)]);
}

#[test]
fn audit_values_one_base() {
    // A bare column starts at the left edge of the region. chr1:100-200 is
    // 0-based 99..200, so the first value is the base at 99.
    let pairs = signal::dense("7\n", &region("chr1:100-200"), None).unwrap();
    assert_eq!(pairs, vec![(TARGET, 7.0)]);
}

#[test]
fn audit_bed_one_base() {
    // BED is 0-based half-open: 99..100 is one base.
    let genes = interval::features("chr1\t99\t100\tone\n", &region("chr1:1-200"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (TARGET, TARGET + 1));
    assert_eq!(genes[0].len(), 1);
}

#[test]
fn audit_gff3_one_base() {
    // GFF3 is 1-based inclusive: 100..100 is one base, 0-based 99..100.
    let text = "##gff-version 3\nchr1\t.\tgene\t100\t100\t.\t+\t.\tID=one\n";
    let genes = interval::features(text, &region("chr1:1-200"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (TARGET, TARGET + 1));
    assert_eq!(genes[0].len(), 1);
    assert_eq!(genes[0].strand, Strand::Forward);
}

#[test]
fn audit_cytoband_one_base() {
    // cytoBand is BED: 0-based half-open.
    let (length, bands) = interval::cytoband("chr1\t99\t100\tp1\tgneg\n", "chr1").unwrap();
    assert_eq!((bands[0].start, bands[0].end), (TARGET, TARGET + 1));
    assert_eq!(bands[0].len(), 1);
    assert_eq!(length, 100);
}

#[test]
fn audit_vcf_one_base() {
    let text = "chr1\t100\t.\tC\tT\t.\t.\t.\n";
    let calls = point::variants(text, &region("chr1:1-200")).unwrap();
    assert_eq!(calls[0].pos, TARGET);
}

#[test]
fn audit_association_one_base() {
    let points = point::associations("100\t1e-9\n", &region("chr1:1-200")).unwrap();
    assert_eq!(points[0].pos, TARGET);
}

#[test]
fn audit_sam_one_base() {
    // POS 100, one aligned base, so the read is exactly 0-based 99..100.
    let text = "r1\t0\tchr1\t100\t60\t1M\t*\t0\t0\tA\t*\n";
    let reads = align::sam(text, &region("chr1:1-200")).unwrap();
    assert_eq!(reads[0].start, TARGET);
    assert_eq!(reads[0].end(), TARGET + 1);
    assert_eq!(reads[0].base_at(TARGET), Some(b'A'));
}

#[test]
fn audit_matrix_one_site() {
    let (sites, rows) = table::matrix("sample\t100\nS1\t1\n", &region("chr1:1-200")).unwrap();
    assert_eq!(sites, vec![TARGET]);
    assert_eq!(rows[0].value(0), Some(1.0));
}

// ---------------------------------------------------------------------------
// The end coordinate, checked apart from the start.
// ---------------------------------------------------------------------------

#[test]
fn audit_bedgraph_end_is_exclusive_and_every_base_is_filled() {
    // 100..200 is 100 bases: 100 through 199, and 200 belongs to no one.
    let pairs = signal::dense("chr1\t100\t200\t5\n", &region("chr1:1-1000"), None).unwrap();
    assert_eq!(pairs.len(), 100, "a bedGraph span fills every base in it");
    assert_eq!(pairs.first().copied(), Some((100, 5.0)));
    assert_eq!(pairs.last().copied(), Some((199, 5.0)));
    assert!(pairs.iter().all(|(_, value)| *value == 5.0));
}

#[test]
fn audit_bed_end_is_exclusive() {
    let genes = interval::features("chr1\t100\t200\tx\n", &region("chr1:1-1000"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (100, 200));
    assert_eq!(genes[0].len(), 100);
}

#[test]
fn audit_gff3_end_is_inclusive_so_the_span_is_one_longer() {
    let text = "##gff-version 3\nchr1\t.\tgene\t100\t200\t.\t+\t.\tID=x\n";
    let genes = interval::features(text, &region("chr1:1-1000"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (99, 200));
    assert_eq!(genes[0].len(), 101);
}

#[test]
fn audit_cytoband_end_is_exclusive() {
    let (length, bands) = interval::cytoband("chr1\t100\t200\tp1\tgneg\n", "chr1").unwrap();
    assert_eq!((bands[0].start, bands[0].end), (100, 200));
    assert_eq!(bands[0].len(), 100);
    assert_eq!(length, 200);
}

#[test]
fn audit_sam_end_walks_the_cigar() {
    // POS 100 is 0-based 99, and 101M reaches 0-based 200 exclusive.
    let text = "r1\t0\tchr1\t100\t60\t101M\t*\t0\t0\t*\t*\n";
    let reads = align::sam(text, &region("chr1:1-1000")).unwrap();
    assert_eq!(reads[0].start, 99);
    assert_eq!(reads[0].end(), 200);
    assert_eq!(reads[0].reference_span(), 101);
}

#[test]
fn audit_window_end_is_exclusive() {
    let found = signal::windows("chr1\t100\t200\t5\n", &region("chr1:1-1000")).unwrap();
    assert_eq!((found[0].start, found[0].end), (100, 200));
}

// ---------------------------------------------------------------------------
// The region edges: an overlap is kept, only what is entirely outside is not.
// ---------------------------------------------------------------------------

/// The window every edge test below is measured against: 0-based 100..200.
fn window() -> Region {
    region("chr1:101-200")
}

#[test]
fn audit_bed_keeps_a_feature_that_overlaps_either_edge() {
    let read = |text: &str| interval::features(text, &window(), None).unwrap();

    // Ends at 101, so its last base is 100, the first base of the window.
    assert_eq!(
        read("chr1\t0\t101\tleft\n").len(),
        1,
        "overlaps the left edge"
    );
    // Ends at 100, so its last base is 99, one short of the window.
    assert!(read("chr1\t0\t100\toutside\n").is_empty(), "entirely left");
    // Starts at 199, the last base of the window.
    assert_eq!(
        read("chr1\t199\t500\tright\n").len(),
        1,
        "overlaps the right edge"
    );
    // Starts at 200, one past the last base of the window.
    assert!(
        read("chr1\t200\t500\toutside\n").is_empty(),
        "entirely right"
    );
}

#[test]
fn audit_gff3_keeps_a_feature_that_overlaps_either_edge() {
    let read = |text: &str| {
        interval::features(&format!("##gff-version 3\n{text}"), &window(), None).unwrap()
    };

    // 1-based 1..101 is 0-based 0..101, whose last base is 100.
    assert_eq!(read("chr1\t.\tgene\t1\t101\t.\t+\t.\tID=a\n").len(), 1);
    // 1-based 1..100 is 0-based 0..100, one short of the window.
    assert!(read("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=b\n").is_empty());
    // 1-based 200 is 0-based 199, the last base of the window.
    assert_eq!(read("chr1\t.\tgene\t200\t500\t.\t+\t.\tID=c\n").len(), 1);
    // 1-based 201 is 0-based 200, one past it.
    assert!(read("chr1\t.\tgene\t201\t500\t.\t+\t.\tID=d\n").is_empty());
}

#[test]
fn audit_bedgraph_fills_only_the_part_of_a_span_inside_the_window() {
    // 0..101 overlaps the window by one base, the one at 100.
    let pairs = signal::dense("chr1\t0\t101\t5\n", &window(), None).unwrap();
    assert_eq!(pairs, vec![(100, 5.0)]);
    // 0..100 stops one base short.
    let pairs = signal::dense("chr1\t0\t100\t5\n", &window(), None).unwrap();
    assert!(pairs.is_empty());
    // 199..500 overlaps by the last base of the window.
    let pairs = signal::dense("chr1\t199\t500\t5\n", &window(), None).unwrap();
    assert_eq!(pairs, vec![(199, 5.0)]);
    // 200..500 starts one past it.
    let pairs = signal::dense("chr1\t200\t500\t5\n", &window(), None).unwrap();
    assert!(pairs.is_empty());
}

#[test]
fn audit_windows_keeps_a_window_that_overlaps_either_edge() {
    let read = |text: &str| signal::windows(text, &window()).unwrap();
    assert_eq!(read("chr1\t0\t101\t5\n").len(), 1);
    assert!(read("chr1\t0\t100\t5\n").is_empty());
    assert_eq!(read("chr1\t199\t500\t5\n").len(), 1);
    assert!(read("chr1\t200\t500\t5\n").is_empty());
}

#[test]
fn audit_sam_keeps_a_read_that_overlaps_either_edge() {
    let read = |pos: u64, cigar: &str| {
        let text = format!("r\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\n");
        align::sam(&text, &window()).unwrap()
    };
    // POS 1 is 0-based 0; 101M reaches 0-based 100, the first base of the window.
    assert_eq!(read(1, "101M").len(), 1, "overlaps the left edge");
    assert!(read(1, "100M").is_empty(), "stops one base short");
    // POS 200 is 0-based 199, the last base of the window.
    assert_eq!(read(200, "50M").len(), 1, "overlaps the right edge");
    assert!(read(201, "50M").is_empty(), "starts one past it");
}

#[test]
fn audit_point_readers_use_the_half_open_end_of_the_window() {
    let vcf = |pos: u64| {
        let text = format!("chr1\t{pos}\t.\tC\tT\t.\t.\t.\n");
        point::variants(&text, &window()).unwrap().len()
    };
    // 1-based 101 is 0-based 100, the first base of the window.
    assert_eq!(vcf(101), 1);
    assert_eq!(vcf(100), 0);
    // 1-based 200 is 0-based 199, the last base.
    assert_eq!(vcf(200), 1);
    assert_eq!(vcf(201), 0);

    let assoc = |pos: u64| {
        let text = format!("{pos}\t0.5\n");
        point::associations(&text, &window()).unwrap().len()
    };
    assert_eq!(assoc(101), 1);
    assert_eq!(assoc(100), 0);
    assert_eq!(assoc(200), 1);
    assert_eq!(assoc(201), 0);

    let site = |pos: u64| {
        let text = format!("sample\t{pos}\nS1\t1\n");
        table::matrix(&text, &window()).unwrap().0.len()
    };
    assert_eq!(site(101), 1);
    assert_eq!(site(100), 0);
    assert_eq!(site(200), 1);
    assert_eq!(site(201), 0);
}

// ---------------------------------------------------------------------------
// Rows the readers do get right, pinned so a later change cannot move them.
// ---------------------------------------------------------------------------

#[test]
fn audit_bedgraph_row_before_the_window_does_not_underflow() {
    let pairs = signal::dense("chr1\t0\t10\t5\n", &window(), None).unwrap();
    assert!(pairs.is_empty());
}

#[test]
fn audit_values_run_out_before_the_window_does() {
    // Four values into a five base window, so the last base stays unset.
    let pairs = signal::dense("1\n2\n3\n4\n", &region("chr1:100-104"), None).unwrap();
    assert_eq!(pairs, vec![(99, 1.0), (100, 2.0), (101, 3.0), (102, 4.0)]);
}

#[test]
fn audit_gff3_without_a_pragma_whose_seventh_column_is_a_dot() {
    // An unstranded GFF3 record is still GFF3, so its start still moves back one.
    let text = "chr1\tRefSeq\tregion\t100\t200\t.\t.\t.\tID=x\n";
    let genes = interval::features(text, &region("chr1:1-1000"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (99, 200));
    assert_eq!(genes[0].strand, Strand::Unknown);
}

#[test]
fn audit_a_bed12_row_read_without_a_pragma() {
    // BED12: column seven is thickStart, so the guess should say BED and the
    // coordinates should pass straight through.
    let text = "chr1\t99\t100\tx\t0\t+\t99\t100\t0,0,0\t1\t1\t0\n";
    let genes = interval::features(text, &region("chr1:1-1000"), None).unwrap();
    assert_eq!((genes[0].start, genes[0].end), (99, 100));
}

#[test]
fn audit_a_gff3_pragma_carrying_a_length_is_not_a_feature() {
    let text =
        "##gff-version 3\n##sequence-region chr1 1 1000\nchr1\t.\tgene\t100\t200\t.\t+\t.\tID=x\n";
    let genes = interval::features(text, &region("chr1:1-1000"), None).unwrap();
    assert_eq!(genes.len(), 1);
    assert_eq!((genes[0].start, genes[0].end), (99, 200));
}

// ---------------------------------------------------------------------------
// Format identification, which is where the readers went wrong rather than in
// the conversions. Every one of these was a figure that came out wrong with no
// error, which is the only kind of bug worth this much fixture.
// ---------------------------------------------------------------------------

#[test]
fn samtools_depth_over_two_bams_is_not_guessed_at_as_a_bedgraph() {
    // `samtools depth a.bam b.bam` writes `chrom pos depth_a depth_b`, four
    // columns, and `Shape::sniff` calls four columns bedGraph. The position
    // becomes a start, the first depth becomes an end and the second depth
    // becomes the value, so one line turns into a run of bases.
    //
    // A deep amplicon is where it never errors: the depth is above the position
    // on every line, so `end < start` never trips.
    let amplicon: String = (1..=5)
        .map(|pos| format!("amplicon\t{pos}\t3000\t2900\n"))
        .collect();
    // Read as a bedGraph the intervals overlap, which no bedGraph does, so the
    // guess is refused and the message names the flag that reads it right.
    let error = signal::dense(&amplicon, &region("amplicon:1-1500"), None).unwrap_err();
    assert!(error.to_string().contains("--format depth"), "{error}");

    // Five lines of a depth file are five positions, 0-based 0 to 4.
    let pairs = signal::dense(
        &amplicon,
        &region("amplicon:1-1500"),
        Some(crate::args::Format::Depth),
    )
    .unwrap();
    assert_eq!(
        pairs,
        vec![
            (0, 3000.0),
            (1, 3000.0),
            (2, 3000.0),
            (3, 3000.0),
            (4, 3000.0)
        ]
    );
}

#[test]
fn a_three_column_file_is_samtools_depth_and_nothing_else() {
    // Three columns is `samtools depth`, and a BED3 handed to a coverage track
    // is read as one: interval 100..200 becomes a single position at 0-based 99
    // carrying a depth of 200. Nothing distinguishes the two, since a BED3 has
    // no value column to notice the absence of, so this is written down rather
    // than guessed at. The help text says a three column file is depth.
    let pairs = signal::dense("chr1\t100\t200\n", &region("chr1:1-1000"), None).unwrap();
    assert_eq!(pairs, vec![(99, 200.0)]);
}

#[test]
fn an_inverted_bed_span_is_refused_the_way_signal_refuses_it() {
    // `signal::dense` and `signal::windows` both call `end < start` an error on
    // the same four columns. `interval::features` does not check, so
    // `Feature::new` widens 200..100 into a one base feature at 200: a gene
    // drawn 100 bases from where either coordinate says it is.
    let read = interval::features("chr1\t200\t100\tbad\n", &region("chr1:1-1000"), None);
    assert!(
        read.is_err(),
        "expected the same error signal.rs gives, got {:?}",
        read.map(|f| (f[0].start, f[0].end))
    );
}

#[test]
fn an_inverted_cytoband_span_is_refused_rather_than_collapsed() {
    let read = interval::cytoband("chr1\t200\t100\tp1\tgneg\n", "chr1");
    assert!(
        read.is_err(),
        "expected an error, got {:?}",
        read.map(|(length, bands)| (length, bands[0].start, bands[0].end))
    );
}

#[test]
fn a_sequence_named_track_keeps_its_rows_whatever_the_separator() {
    // `lines` drops anything starting with `track ` or `browser `. A tab
    // separated file survives, because the prefix is `track\t`. A space
    // separated one does not, and the row disappears with no error.
    let spaced = interval::features("track 99 100 x\n", &region("track:1-200"), None).unwrap();
    let tabbed = interval::features("track\t99\t100\tx\n", &region("track:1-200"), None).unwrap();
    assert_eq!(tabbed.len(), 1, "the tab separated row survives");
    assert_eq!(
        spaced.len(),
        1,
        "the same row separated by spaces was dropped as a UCSC track line"
    );
}

#[test]
fn a_deletion_reaching_into_the_window_is_a_call_about_the_window() {
    // VCF POS 100 is 0-based 99, one base left of a window starting at 100, and
    // REF spells 99..104. The deletion takes four bases out of the window and
    // the track is filtered on the anchor alone, so nothing is drawn.
    let text = "chr1\t100\t.\tCAAAA\tC\t.\t.\t.\n";
    let calls = point::variants(text, &window()).unwrap();
    assert_eq!(calls.len(), 1, "the deletion covers 100..104 of the window");
}
