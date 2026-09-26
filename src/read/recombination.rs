//! Recombination rates along a chromosome: a genetic map, as HapMap and the
//! imputation panels write one, or a bedGraph of rates.
//!
//! A genetic map gives a rate at each of its positions, in centimorgans per
//! megabase, and the rate holds from that position to the next one. HapMap
//! writes a chromosome, a position, the rate and the map distance, under a
//! header naming them:
//!
//! ```text
//! Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
//! chr22       16051347      8.096992     0.000000
//! ```
//!
//! The maps that come with IMPUTE2 and SHAPEIT are one chromosome to a file and
//! leave the chromosome out: `position COMBINED_rate(cM/Mb) Genetic_Map(cM)`.
//! Both are read by the names in their header, and their positions count from
//! one. A file with no header is a bedGraph of rates, 0-based and half-open.

use crate::read::signal::spans;
use crate::region::Region;
use crate::Format;

use super::{fields, lines, number, ReadError};

/// The names a column of positions goes by.
const POSITION: &[&str] = &["position(bp)", "position", "pos", "bp", "physical_position"];
/// The names a column of rates goes by.
const RATE: &[&str] = &[
    "rate(cm/mb)",
    "combined_rate(cm/mb)",
    "rate",
    "cm/mb",
    "recombination_rate",
];
/// The names a column of chromosomes goes by.
const CHROMOSOME: &[&str] = &["chromosome", "chrom", "chr"];

/// The column whose header is one of `names`, compared without case.
fn column(header: &[&str], names: &[&str]) -> Option<usize> {
    header
        .iter()
        .position(|field| names.iter().any(|name| field.eq_ignore_ascii_case(name)))
}

/// Reads the rates over the region, as 0-based half-open `(start, end, rate)`
/// spans in centimorgans per megabase.
///
/// A genetic map is sorted by position before its spans are made, so a map
/// whose rows are out of order still gives each rate the stretch up to the
/// next position. Its last position ends the map, so the rate written there
/// covers that one base. A rate that is empty or `NA` is a stretch the map says
/// nothing about, and is left out rather than drawn at nought.
pub fn rates(text: &str, region: &Region) -> Result<Vec<(u64, u64, f64)>, ReadError> {
    let Some((_, head)) = lines(text).next() else {
        return Ok(Vec::new());
    };
    let header = fields(head);
    let (Some(at), Some(rate)) = (column(&header, POSITION), column(&header, RATE)) else {
        // No header naming a position and a rate: a bedGraph of rates.
        return spans(text, region, Some(Format::BedGraph));
    };
    let chromosome = column(&header, CHROMOSOME);
    let mut points: Vec<(u64, Option<f64>)> = Vec::new();
    for (line, row) in lines(text).skip(1) {
        let row = fields(row);
        let field = |index: usize| row.get(index).copied().unwrap_or_default();
        if chromosome.is_some_and(|c| field(c) != region.seq()) {
            continue;
        }
        let position: u64 = number(field(at), "the position", line)?;
        if position == 0 {
            return Err(ReadError::at(
                line,
                "a position of 0: a genetic map counts positions from 1",
            ));
        }
        let value = field(rate);
        let value = if value.is_empty() || value.eq_ignore_ascii_case("na") {
            None
        } else {
            let value: f64 = number(value, "the rate", line)?;
            if !value.is_finite() || value < 0.0 {
                return Err(ReadError::at(
                    line,
                    format!("a rate of {value}: a rate is a number of centimorgans per megabase, nought or more"),
                ));
            }
            Some(value)
        };
        points.push((position - 1, value));
    }
    points.sort_by_key(|(start, _)| *start);
    let mut out = Vec::new();
    for (index, (start, value)) in points.iter().enumerate() {
        let end = points
            .get(index + 1)
            .map_or(start + 1, |(next, _)| (*next).max(start + 1));
        let Some(value) = value else {
            continue;
        };
        if end > region.start() && *start < region.end() {
            out.push((*start, end, *value));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(locus: &str) -> Region {
        Region::parse(locus).unwrap()
    }

    #[test]
    fn a_hapmap_map_holds_each_rate_to_the_next_position() {
        let map = "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n\
                   chr22\t101\t2.5\t0.0\n\
                   chr22\t201\t40.0\t0.00025\n\
                   chr22\t301\t0.0\t0.0042\n\
                   chr21\t150\t9.0\t0.0\n";
        let read = rates(map, &region("chr22:1-1,000")).unwrap();
        assert_eq!(
            read,
            vec![(100, 200, 2.5), (200, 300, 40.0), (300, 301, 0.0)]
        );
    }

    #[test]
    fn an_impute_map_has_no_chromosome_and_a_bedgraph_no_header() {
        let map = "position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n301 1.5 0.0\n101 0.5 0.0\n";
        let read = rates(map, &region("chr1:1-1,000")).unwrap();
        assert_eq!(read, vec![(100, 300, 0.5), (300, 301, 1.5)], "sorted first");
        let bedgraph = "chr1\t0\t100\t3.5\nchr1\t100\t200\t0.2\n";
        let read = rates(bedgraph, &region("chr1:1-1,000")).unwrap();
        assert_eq!(read, vec![(0, 100, 3.5), (100, 200, 0.2)]);
    }

    #[test]
    fn a_missing_rate_is_left_out_and_a_negative_one_refused() {
        let map = "position rate\n101 NA\n201 3\n301 1\n";
        let read = rates(map, &region("chr1:1-1,000")).unwrap();
        assert_eq!(read, vec![(200, 300, 3.0), (300, 301, 1.0)]);
        let error = rates("position rate\n101 -2\n", &region("chr1:1-1,000")).unwrap_err();
        assert!(error.reason.contains("nought or more"), "{}", error.reason);
    }
}
