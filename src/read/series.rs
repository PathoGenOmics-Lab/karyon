//! Tables whose rows are not places on a sequence: counts over time, an
//! estimate over time, a test at each site of a gene, and the current of a
//! nanopore read.
//!
//! Each is its own axis. A week, a codon or a sample is a whole unit counted
//! from one, and a figure of one of these is drawn over the units the file
//! holds rather than over a region of a genome. Every reader here gives each
//! row the coordinate one less than the number it is written as, so the ruler,
//! which numbers from one, prints the number the file wrote.
//!
//! # Read by their headers
//!
//! The tools that write these tables do not agree on a column order, and most
//! of them agree on names, so each column is found by what its header calls
//! it, without regard to case: `week` or `time`, `lineage` or `mutation`,
//! `mean` or `estimate`, `alpha` and `beta` as HyPhy writes the rates. Tabs,
//! commas and runs of spaces all separate columns, since R writes commas and
//! most command lines tabs.

use crate::read::{fields, lines, ReadError};
use crate::track::phylodynamics::PhylodynamicPoint;
use crate::track::selection::SelectionSite;
use crate::track::surveillance::SurveillanceObservation;

/// The names a column of times goes by.
const TIME: &[&str] = &[
    "time",
    "week",
    "day",
    "month",
    "year",
    "date",
    "t",
    "t_end",
    "passage",
    "generation",
];

/// The column whose header is one of `names`, compared without case.
fn column(header: &[&str], names: &[&str]) -> Option<usize> {
    header
        .iter()
        .position(|field| names.iter().any(|name| field.eq_ignore_ascii_case(name)))
}

/// A column the table has to have, or the error naming what it could be
/// called and what the header does call its columns.
fn needed(line: usize, header: &[&str], what: &str, names: &[&str]) -> Result<usize, ReadError> {
    column(header, names).ok_or_else(|| {
        ReadError::at(
            line,
            format!(
                "no column of {what}: it is found by its header, as {}, and this header is {}",
                names
                    .iter()
                    .take(4)
                    .map(|name| format!("`{name}`"))
                    .collect::<Vec<_>>()
                    .join(", "),
                header.join(" ")
            ),
        )
    })
}

/// A whole unit of time, counted from one or from any later number.
///
/// A week, a day, a month or a year: a fraction of one is refused with the
/// way round it, since the axis counts whole units and a skyline in decimal
/// years would be drawn rounded without a word.
fn time(value: &str, line: usize) -> Result<u64, ReadError> {
    if let Ok(whole) = value.parse::<u64>() {
        if whole == 0 {
            return Err(ReadError::at(
                line,
                "a time of 0: times are counted from 1, as week 1, day 1 or a year",
            ));
        }
        return Ok(whole - 1);
    }
    match value.parse::<f64>() {
        Ok(number) if number.is_finite() => Err(ReadError::at(
            line,
            format!(
                "a time of {value}: times are whole units, a week, a day, a month or a \
                 year; count in a smaller unit, as days, to keep the fraction"
            ),
        )),
        // A date is the usual time that is not a number, and says so by
        // its separators.
        _ if value.contains(['-', '/']) => Err(ReadError::at(
            line,
            format!(
                "a time of {value:?}: times are whole numbers, as week 12 or day 340; \
                 count dates from a start, as days or weeks since the first sample"
            ),
        )),
        _ => Err(ReadError::at(
            line,
            format!("a time of {value:?}, which is not a number of whole units"),
        )),
    }
}

/// A number in a column that has to hold one.
fn value(field: &str, what: &str, line: usize) -> Result<f64, ReadError> {
    field
        .parse::<f64>()
        .ok()
        .filter(|number| number.is_finite())
        .ok_or_else(|| ReadError::at(line, format!("{what} is not a number: {field:?}")))
}

/// A number in a column that may be empty or `NA`.
fn optional(field: Option<&str>, what: &str, line: usize) -> Result<Option<f64>, ReadError> {
    match field {
        None => Ok(None),
        Some(field) if field.is_empty() || field.eq_ignore_ascii_case("na") => Ok(None),
        Some(field) => value(field, what, line).map(Some),
    }
}

/// How many of each group were seen at each time, out of how many: the
/// lineages of a surveillance programme, or the reads carrying each mutation
/// in each sample of an evolving population.
///
/// The columns are a time (`week`, `day`, `time`...), a group (`lineage`,
/// `mutation`, `variant`...), a count (`count`, `n`, `alt`...) and the total
/// it is out of (`total`, `depth`...). Returns the rows and the name the time
/// column goes by, which is what the ruler under them counts.
pub fn counts(text: &str) -> Result<(Vec<SurveillanceObservation>, String), ReadError> {
    let mut rows = lines(text);
    let (line, head) = rows
        .next()
        .ok_or_else(|| ReadError::whole("an empty table: a header and a row per count"))?;
    let header = fields(head);
    let at = needed(line, &header, "times", TIME)?;
    let group = needed(
        line,
        &header,
        "groups",
        &[
            "lineage", "mutation", "variant", "group", "clade", "genotype", "allele", "name",
        ],
    )?;
    let count = needed(
        line,
        &header,
        "counts",
        &[
            "count", "n", "alt", "observed", "reads", "cases", "carriers",
        ],
    )?;
    let total = needed(
        line,
        &header,
        "totals",
        &[
            "total",
            "depth",
            "denominator",
            "sequenced",
            "samples",
            "all",
        ],
    )?;
    let unit = header[at].to_ascii_lowercase();
    let mut out = Vec::new();
    for (line, row) in rows {
        let row = fields(row);
        let field = |index: usize| row.get(index).copied().unwrap_or_default();
        let seen = value(field(count), "the count", line)?;
        let of = value(field(total), "the total", line)?;
        if seen < 0.0 || of < 0.0 || seen.fract() != 0.0 || of.fract() != 0.0 {
            return Err(ReadError::at(
                line,
                format!("a count and a total are whole numbers of nought or more: {seen} of {of}"),
            ));
        }
        if seen > of {
            return Err(ReadError::at(
                line,
                format!("a count of {seen} out of a total of {of}"),
            ));
        }
        out.push(SurveillanceObservation::new(
            time(field(at), line)?,
            field(group),
            seen as u64,
            of as u64,
        ));
    }
    if out.is_empty() {
        return Err(ReadError::whole("a header and no row under it"));
    }
    Ok((out, unit))
}

/// An estimate at each time with the interval around it, as a skyline, an
/// effective population size or a reproductive number is written.
///
/// The columns are a time, an estimate (`estimate`, `mean`, `median`...)
/// and, where the file has them, the two ends of its interval (`lower` and
/// `upper`, or the quantiles EpiEstim writes). Returns the points and the name
/// of the time column.
pub fn estimates(text: &str) -> Result<(Vec<PhylodynamicPoint>, String), ReadError> {
    let mut rows = lines(text);
    let (line, head) = rows
        .next()
        .ok_or_else(|| ReadError::whole("an empty table: a header and a row per time"))?;
    let header = fields(head);
    let at = needed(line, &header, "times", TIME)?;
    let estimate = needed(
        line,
        &header,
        "estimates",
        &[
            "estimate",
            "mean",
            "median",
            "value",
            "ne",
            "r",
            "mean(r)",
            "median(r)",
        ],
    )?;
    let lower = column(
        &header,
        &[
            "lower",
            "lower95",
            "hpd_lower",
            "lower_hpd",
            "low",
            "quantile.0.025(r)",
            "quantile 0.025(r)",
            "q2.5",
        ],
    );
    let upper = column(
        &header,
        &[
            "upper",
            "upper95",
            "hpd_upper",
            "upper_hpd",
            "high",
            "quantile.0.975(r)",
            "quantile 0.975(r)",
            "q97.5",
        ],
    );
    let unit = header[at].to_ascii_lowercase();
    let mut out = Vec::new();
    for (line, row) in rows {
        let row = fields(row);
        let field = |index: usize| row.get(index).copied().unwrap_or_default();
        let mut point = PhylodynamicPoint::new(
            time(field(at), line)?,
            value(field(estimate), "the estimate", line)?,
        );
        let low = optional(lower.map(field), "the lower end", line)?;
        let high = optional(upper.map(field), "the upper end", line)?;
        if let (Some(low), Some(high)) = (low, high) {
            point = point.interval(low, high);
        }
        out.push(point);
    }
    if out.is_empty() {
        return Err(ReadError::whole("a header and no row under it"));
    }
    Ok((out, unit))
}

/// A test of selection at each site of a gene, as HyPhy's FEL and MEME, a
/// PAML site model or a table of your own write it.
///
/// A site column (`site`, `codon`...) where there is one, else the rows in
/// order from site 1. The rates as `alpha` and `beta` (or `dS` and `dN`), or
/// their ratio as `omega`; the evidence as a `p-value` or a `posterior`; and
/// MEME's `beta-`, `beta+` and `p+` where they are there.
pub fn selection(text: &str) -> Result<Vec<SelectionSite>, ReadError> {
    let mut rows = lines(text);
    let (line, head) = rows
        .next()
        .ok_or_else(|| ReadError::whole("an empty table: a header and a row per site"))?;
    let header = fields(head);
    let site = column(&header, &["site", "codon", "pos", "position"]);
    let synonymous = column(&header, &["alpha", "ds", "syn", "synonymous"]);
    let nonsynonymous = column(&header, &["beta", "dn", "nonsyn", "nonsynonymous"]);
    let omega = column(&header, &["omega", "dn/ds", "w"]);
    let p_value = column(&header, &["p-value", "pvalue", "p", "p value"]);
    let posterior = column(
        &header,
        &[
            "posterior",
            "prob",
            "pr(w>1)",
            "posterior probability",
            "beb",
        ],
    );
    let (minus, plus, weight) = (
        column(&header, &["beta-"]),
        column(&header, &["beta+"]),
        column(&header, &["p+"]),
    );
    if (synonymous.is_none() || nonsynonymous.is_none()) && omega.is_none() {
        return Err(ReadError::at(
            line,
            format!(
                "no rates: a site is read as `alpha` and `beta`, `dS` and `dN`, or their ratio \
                 `omega`, and this header is {}",
                header.join(" ")
            ),
        ));
    }
    let mut out = Vec::new();
    for (index, (line, row)) in rows.enumerate() {
        let row = fields(row);
        let field = |at: usize| row.get(at).copied();
        let number = match site {
            Some(at) => field(at)
                .unwrap_or_default()
                .parse::<u64>()
                .ok()
                .filter(|site| *site > 0)
                .ok_or_else(|| {
                    ReadError::at(
                        line,
                        "a site is a whole number from 1, as HyPhy counts them",
                    )
                })?,
            None => index as u64 + 1,
        };
        let mut tested = SelectionSite::new(number - 1);
        match (synonymous, nonsynonymous) {
            (Some(s), Some(n)) => {
                let s = optional(field(s), "the synonymous rate", line)?;
                let n = optional(field(n), "the nonsynonymous rate", line)?;
                if let (Some(s), Some(n)) = (s, n) {
                    tested = tested.rates(s, n);
                }
            }
            _ => {
                if let Some(omega) = optional(omega.and_then(field), "omega", line)? {
                    tested = tested.rates(1.0, omega);
                }
            }
        }
        if let Some(p) = optional(p_value.and_then(field), "the p-value", line)? {
            tested = tested.p_value(p);
        }
        if let Some(p) = optional(posterior.and_then(field), "the posterior", line)? {
            tested = tested.posterior(p);
        }
        if let (Some(minus), Some(plus), Some(weight)) = (
            optional(minus.and_then(field), "beta-", line)?,
            optional(plus.and_then(field), "beta+", line)?,
            optional(weight.and_then(field), "p+", line)?,
        ) {
            tested = tested.episodic_rates(minus, plus, weight);
        }
        out.push(tested);
    }
    if out.is_empty() {
        return Err(ReadError::whole("a header and no row under it"));
    }
    Ok(out)
}

/// The current of one nanopore read, as [`squiggle`] reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct Signal {
    /// The samples, in picoamperes, in the order they were measured.
    pub samples: Vec<f64>,
    /// The read they belong to, or `signal` for a column of numbers.
    pub read: String,
    /// How many reads the file holds, this one among them.
    pub reads: usize,
}

/// The current of one nanopore read, in picoamperes, and what it is called.
///
/// A SLOW5 file as `slow5tools view` writes it, its raw signal put into
/// picoamperes with the read's own digitisation, offset and range; the read
/// named `read`, or the first. Or plain numbers, one sample after another,
/// which are taken as picoamperes already.
pub fn squiggle(text: &str, read: Option<&str>) -> Result<Signal, ReadError> {
    let slow5 = text
        .lines()
        .any(|line| line.starts_with("#read_id") || line.starts_with("#slow5_version"));
    if !slow5 {
        let mut signal = Vec::new();
        for (line, row) in lines(text) {
            for field in fields(row) {
                if !field.is_empty() {
                    signal.push(value(field, "a sample of the signal", line)?);
                }
            }
        }
        if signal.is_empty() {
            return Err(ReadError::whole(
                "no samples: a number per sample of the signal",
            ));
        }
        return Ok(Signal {
            samples: signal,
            read: "signal".to_string(),
            reads: 1,
        });
    }
    let header_line = text
        .lines()
        .find(|line| line.starts_with("#read_id"))
        .ok_or_else(|| ReadError::whole("a SLOW5 file with no #read_id header"))?;
    let header: Vec<&str> = header_line.trim_start_matches('#').split('\t').collect();
    let at = |name: &str| header.iter().position(|field| *field == name);
    let (Some(id), Some(raw)) = (at("read_id"), at("raw_signal")) else {
        return Err(ReadError::whole(
            "a SLOW5 header with no read_id or no raw_signal column",
        ));
    };
    let (digitisation, offset, range) = (at("digitisation"), at("offset"), at("range"));
    let rows = || {
        text.lines().enumerate().filter(|(_, line)| {
            !line.starts_with('#') && !line.starts_with('@') && !line.trim().is_empty()
        })
    };
    let name_of = |line: &str| line.split('\t').nth(id).unwrap_or_default().to_string();
    let reads = rows().count();
    for (index, line) in rows() {
        let row: Vec<&str> = line.split('\t').collect();
        let name = row.get(id).copied().unwrap_or_default();
        if read.is_some_and(|wanted| wanted != name) {
            continue;
        }
        let number = |at: Option<usize>, fallback: f64| {
            at.and_then(|at| row.get(at))
                .and_then(|field| field.parse::<f64>().ok())
                .unwrap_or(fallback)
        };
        let (scale, shift) = match (digitisation, range) {
            (Some(_), Some(_)) => {
                let steps = number(digitisation, 1.0);
                (number(range, steps) / steps, number(offset, 0.0))
            }
            _ => (1.0, 0.0),
        };
        let mut signal = Vec::new();
        for sample in row.get(raw).copied().unwrap_or_default().split(',') {
            let raw: f64 = sample.trim().parse().map_err(|_| {
                ReadError::at(
                    index + 1,
                    format!("a raw sample is not a number: {sample:?}"),
                )
            })?;
            signal.push((raw + shift) * scale);
        }
        return Ok(Signal {
            samples: signal,
            read: name.to_string(),
            reads,
        });
    }
    Err(ReadError::whole(match read {
        Some(wanted) => {
            let names: Vec<String> = rows().take(3).map(|(_, line)| name_of(line)).collect();
            let more = if reads > names.len() {
                format!(" and {} more", reads - names.len())
            } else {
                String::new()
            };
            format!(
                "no read called {wanted}; the file holds {}{more}",
                names.join(", ")
            )
        }
        None => "a SLOW5 file with no read in it".to_string(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_are_found_by_their_headers_and_numbered_as_written() {
        let text = "week,lineage,count,total\n1,BA.2,30,100\n2,BA.2,55,110\n2,BA.5,40,110\n";
        let (rows, unit) = counts(text).unwrap();
        assert_eq!(unit, "week");
        assert_eq!(rows.len(), 3);
        // Week 1 is coordinate 0, which the ruler prints as 1.
        assert_eq!(rows[0].time, 0);
        assert_eq!(rows[2].lineage, "BA.5");
        let error = counts("week\tlineage\tcount\n1\tA\t3\n").unwrap_err();
        assert!(error.reason.contains("totals"), "{}", error.reason);
        let error = counts("week\tlineage\tcount\ttotal\n1.5\tA\t3\t9\n").unwrap_err();
        assert!(error.reason.contains("whole units"), "{}", error.reason);
        let error = counts("week\tlineage\tcount\ttotal\n1\tA\t12\t9\n").unwrap_err();
        assert!(error.reason.contains("out of a total"), "{}", error.reason);
    }

    #[test]
    fn epiestim_is_read_as_r_writes_it_and_a_date_is_told_what_to_do() {
        // EpiEstim's estimate_R table through write.csv, quotes and all.
        let text = "\"t_start\",\"t_end\",\"Mean(R)\",\"Std(R)\",\"Quantile.0.025(R)\",\
                    \"Median(R)\",\"Quantile.0.975(R)\"\n2,8,1.4,0.2,1.1,1.39,1.8\n";
        let (points, unit) = estimates(text).unwrap();
        assert_eq!(unit, "t_end");
        assert_eq!(points[0].time, 7);
        assert_eq!(points[0].bounds(), Some((1.1, 1.8)));
        let error = counts("date,lineage,count,total\n2024-03-01,A,3,9\n").unwrap_err();
        assert!(
            error.reason.contains("days or weeks since"),
            "{}",
            error.reason
        );
    }

    #[test]
    fn estimates_take_the_interval_where_the_file_has_one() {
        let text = "year\tmean\tlower\tupper\n2015\t120\t80\t190\n2016\t150\tNA\tNA\n";
        let (points, unit) = estimates(text).unwrap();
        assert_eq!(unit, "year");
        assert_eq!(points[0].time, 2014);
        assert_eq!(points.len(), 2);
    }

    #[test]
    fn a_site_is_read_as_hyphy_writes_it() {
        let fel = "alpha,beta,alpha=beta,LRT,p-value\n1,0.2,0.6,3.1,0.08\n0.5,2.5,1.5,6.2,0.01\n";
        let sites = selection(fel).unwrap();
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[1].pos, 1);
        let error = selection("site\tp-value\n1\t0.2\n").unwrap_err();
        assert!(error.reason.contains("no rates"), "{}", error.reason);
    }

    #[test]
    fn a_squiggle_is_plain_numbers_or_a_slow5_read_in_picoamperes() {
        let plain = squiggle("80.5\n81\n79.2 90\n", None).unwrap();
        assert_eq!((plain.samples.len(), plain.read.as_str()), (4, "signal"));
        let slow5 = "#slow5_version\t0.2.0\n#num_read_groups\t1\n\
                     #read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\tlen_raw_signal\traw_signal\n\
                     r1\t0\t2048\t10\t1024\t4000\t3\t100,110,120\n\
                     r2\t0\t2048\t10\t1024\t4000\t2\t200,210\n";
        let second = squiggle(slow5, Some("r2")).unwrap();
        assert_eq!((second.read.as_str(), second.reads), ("r2", 2));
        // (raw + offset) * range / digitisation: (200 + 10) * 1024 / 2048.
        assert_eq!(second.samples, [105.0, 110.0]);
        assert_eq!(squiggle(slow5, None).unwrap().read, "r1");
        let error = squiggle(slow5, Some("r9")).unwrap_err();
        assert!(error.reason.contains("r1, r2"), "{}", error.reason);
    }
}
