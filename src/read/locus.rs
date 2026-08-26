//! Gene neighbourhoods from several genomes, and what joins one to the next.
//!
//! Two files, because they are two different things. The genes are intervals
//! and come from the same BED or GFF3 anything else does; the homologies are a
//! search result and come from whatever ran the search. [`loci`] reads the
//! first and [`links`] the second, and [`links`] needs the output of [`loci`]
//! because of what a [`Homology`] is.
//!
//! # Column one names the row rather than selecting it
//!
//! This is the one place an interval file is read differently from everywhere
//! else in this directory. [`features`](super::interval::features) takes column
//! one as the sequence and drops every row that is not the one being drawn.
//! Here every row is drawn and column one says which genome it belongs to, so
//! the file is the concatenation an ordinary shell produces:
//!
//! ```text
//! cat H37Rv.bed CDC1551.bed Erdman.bed > loci.bed
//! ```
//!
//! Rows keep their file order and so do the genomes, because the track draws
//! the rows in the order it is given them and a homology names its row by
//! number. Sorting anything after this reader has run silently re-points every
//! ribbon on the figure.
//!
//! Coordinates are whatever the interval format says: BED 0-based half-open and
//! passed straight through, GFF3 1-based inclusive and moved back one at the
//! start. Loci from different genomes share the figure's one axis, so their
//! coordinates are compared as they arrive; [`Locus::offset`] is how a row is
//! lined up with its neighbour, and that is a decision about the figure rather
//! than about the file, so it is not made here.
//!
//! # A homology names genes and the track counts them
//!
//! [`Homology`] refers to genes by their position in a [`Locus`], and every
//! file that produces homologies names them instead. Somebody has to do that
//! join, and doing it wrongly is close to undetectable: a name that resolved to
//! the wrong gene draws a well formed ribbon with a confident label, and a name
//! that resolved to nothing at all draws no ribbon and *removes the unmatched
//! outline from the gene at the other end*, so the figure claims a match it
//! cannot show.
//!
//! So the join is here, it is exact, and every way it can fail is counted
//! rather than resolved. A name that is not a gene is listed in
//! [`Links::unjoined`] and its row is dropped without a [`Homology`] being
//! built at all. A name that is more than one gene stops the read, because
//! there is no answer to which one was meant.
//!
//! # Nothing joined is the loudest wrong answer this crate can draw
//!
//! [`LocusTrack::mark_unmatched`](crate::LocusTrack::mark_unmatched) is on by
//! default and outlines every gene no homology reaches, which is right: the
//! finding in a comparison of gene neighbourhoods is what is missing. A links
//! file whose names did not join therefore produces a figure in which *every
//! gene in every genome* carries the heaviest mark the track owns, and that
//! reads as a discovery rather than as a mistake.
//!
//! The names in a search result come from the FASTA it was run against, and the
//! names in an annotation come from its ninth column, and those are routinely
//! not the same strings: `lcl|NC_000962.3_cds_NP_215181.1_667` against `rpoB`.
//! The join rate is the one fact the figure has no way to show, which is why
//! [`Links::unjoined`] exists and why a caller is expected to refuse on it.

use std::collections::{BTreeMap, BTreeSet};

use crate::{Feature, Homology, Locus, Region};

use super::interval::{bed, flavour, gff3, Flavour};
use super::{columns, lines, number, Format, ReadError};

/// The loci a file holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq)]
pub struct Loci {
    /// One locus per distinct name in column one, in first-seen order.
    pub loci: Vec<Locus>,
    /// Gene records in the file, before any filter.
    pub records: usize,
    /// Records that touch no base of the window.
    pub off_region: usize,
}

/// The homologies a file holds, and every row that did not become one.
#[derive(Debug, Clone, PartialEq)]
pub struct Links {
    /// The links, both ends resolved, in the order the pairs were first seen.
    pub links: Vec<Homology>,
    /// Rows in the file, comments and headers already gone.
    pub records: usize,
    /// Names no locus had a gene by, de-duplicated, in first-seen order.
    ///
    /// The count that decides whether the figure means anything. A links file
    /// against the wrong names joins nothing, and nothing joined draws as every
    /// gene in every genome being unique to it.
    pub unjoined: Vec<String>,
    /// Rows whose two genes were in the same locus.
    pub self_hits: usize,
    /// Rows joining two loci that are not next to each other in the stack.
    pub not_adjacent: usize,
    /// Rows folded onto a pair already seen.
    pub merged: usize,
}

/// Whether an identity column is a percentage or a fraction.
///
/// The tools disagree and the file does not say. BLAST and DIAMOND write
/// column three as a percentage from 0 to 100; some others write the same
/// column as a fraction from 0 to 1. Guessed wrongly in one direction every
/// ribbon becomes a perfect match, and in the other every ribbon becomes the
/// palest shade on the ramp, and neither fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identity {
    /// 0 to 100, as BLAST and DIAMOND write it.
    Percent,
    /// 0 to 1.
    Fraction,
}

impl Identity {
    /// Parses a unit name, as `--identity` spells it.
    pub fn parse(word: &str) -> Option<Identity> {
        Some(match word {
            "percent" | "percentage" => Identity::Percent,
            "fraction" | "proportion" => Identity::Fraction,
            _ => return None,
        })
    }

    /// The largest value this unit allows.
    fn ceiling(self) -> f64 {
        match self {
            Identity::Percent => 100.0,
            Identity::Fraction => 1.0,
        }
    }
}

/// Reads gene neighbourhoods, one row per distinct name in column one.
///
/// BED is 0-based half-open and passes straight through; GFF3 is 1-based
/// inclusive and its start comes back one lower. The two are told apart the
/// same way [`features`](super::interval::features) tells them apart, and
/// `format` overrules the guess.
///
/// Column one names the genome rather than selecting it, which is the one place
/// this reader differs from every other interval reader here.
///
/// # Errors
///
/// Returns the first row that will not parse, on its line. A file that parses
/// and holds nothing inside the window is not an error: [`Loci::records`] and
/// [`Loci::off_region`] say which of the two happened.
pub fn loci(text: &str, region: &Region, format: Option<Format>) -> Result<Loci, ReadError> {
    let flavour = flavour(text, format);
    let mut found = Loci {
        loci: Vec::new(),
        records: 0,
        off_region: 0,
    };
    // The order genomes were first seen in, which is the order they stack.
    let mut order: Vec<String> = Vec::new();
    let mut genes: Vec<Vec<Feature>> = Vec::new();

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.is_empty() {
            continue;
        }
        found.records += 1;
        let feature = match flavour {
            Flavour::Bed => bed(&cols, at)?,
            Flavour::Gff3 => gff3(&cols, at)?,
        };
        if feature.end <= region.start() || feature.start >= region.end() {
            found.off_region += 1;
            continue;
        }
        let row = match order.iter().position(|name| name == cols[0]) {
            Some(row) => row,
            None => {
                order.push(cols[0].to_string());
                genes.push(Vec::new());
                order.len() - 1
            }
        };
        genes[row].push(feature);
    }

    found.loci = order
        .into_iter()
        .zip(genes)
        .map(|(name, genes)| Locus::new(name, genes))
        .collect();
    Ok(found)
}

/// Reads homologies and joins their names to the genes `loci` holds.
///
/// Takes BLAST tabular output (`-outfmt 6`, and `-outfmt 7`, whose extra lines
/// are comments), which DIAMOND and others write too: query name, subject name,
/// identity, then nine more columns this reader does not need. A three column
/// `query subject identity` file and a two column `query subject` file are read
/// as well, the second with every identity unstated.
///
/// No coordinate is read. A [`Homology`] holds none: it names two genes, and
/// where those genes are is already in `loci`.
///
/// `identity` says whether column three is a percentage or a fraction. Left as
/// `None` it is worked out from the values, and refused where they cannot say.
///
/// # Errors
///
/// Returns a row that will not parse, an identity outside its own unit, an
/// identity column that could be either unit, or a name that more than one gene
/// answers to. A row that names a gene nobody has is not an error: it is
/// counted in [`Links::unjoined`], and the caller decides what a low join rate
/// means.
pub fn links(text: &str, loci: &[Locus], identity: Option<Identity>) -> Result<Links, ReadError> {
    let index = index(loci)?;
    let rows = rows(text)?;

    let unit = match identity {
        Some(unit) => Some(unit),
        None => infer(&rows)?,
    };

    let mut found = Links {
        links: Vec::new(),
        records: rows.len(),
        unjoined: Vec::new(),
        self_hits: 0,
        not_adjacent: 0,
        merged: 0,
    };
    let mut missing: BTreeSet<String> = BTreeSet::new();
    let mut made: BTreeSet<(usize, usize, usize)> = BTreeSet::new();

    for row in &rows {
        let (Some(&query), Some(&subject)) = (index.get(&row.query), index.get(&row.subject))
        else {
            // Listed rather than resolved. A name pointed at gene nought draws
            // a ribbon into the leftmost gene and strips its outline, and there
            // is nothing on the figure to say the ribbon was invented.
            for name in [&row.query, &row.subject] {
                if !index.contains_key(name) && missing.insert(name.clone()) {
                    found.unjoined.push(name.clone());
                }
            }
            continue;
        };

        if query.0 == subject.0 {
            // A search run against everything reports every gene against
            // itself, and a paralogue pair inside one genome looks the same.
            // Neither is a link between two rows, and the track has no way to
            // draw one that stays inside a row.
            found.self_hits += 1;
            continue;
        }
        if query.0.abs_diff(subject.0) != 1 {
            // A ribbon is between neighbouring rows, so a hit two genomes apart
            // has nowhere to go. Drawn at `row = min` it would name the wrong
            // genome in its own tooltip.
            found.not_adjacent += 1;
            continue;
        }

        // Which of the two is the upper row is the file's business and not the
        // figure's, so a reciprocal pair folds onto one key here.
        let (upper, lower) = if query.0 < subject.0 {
            (query, subject)
        } else {
            (subject, query)
        };
        let key = (upper.0, upper.1, lower.1);
        if !made.insert(key) {
            // Several rows for one pair: the fragments of one alignment, or the
            // same pair reported from both ends. Left in, they paint the same
            // ribbon several times and the last one drawn is the one seen, so
            // the identity on the page is an artefact of the file's order.
            found.merged += 1;
            continue;
        }

        found.links.push(match row.identity {
            Some(value) => Homology::new(key.0, key.1, key.2, scale(value, unit, row.line)?),
            None => Homology::unstated(key.0, key.1, key.2),
        });
    }

    Ok(found)
}

/// One row of a homology file, before anything is joined.
#[derive(Debug)]
struct Row {
    query: String,
    subject: String,
    identity: Option<f64>,
    line: usize,
}

/// Every gene name, and the one gene it is.
///
/// A name that is two genes is refused rather than resolved: the figure draws
/// two genes with the same name identically, so nothing on it could show which
/// one a ribbon meant.
fn index(loci: &[Locus]) -> Result<BTreeMap<String, (usize, usize)>, ReadError> {
    let mut index: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for (row, locus) in loci.iter().enumerate() {
        for (at, gene) in locus.genes.iter().enumerate() {
            // A gene with no name cannot be named by a homology file, and an
            // empty name is not a name: left in the index it would answer to a
            // blank field and match whatever came first.
            let Some(name) = gene
                .name
                .as_deref()
                .map(str::trim)
                .filter(|n| !n.is_empty())
            else {
                continue;
            };
            if index.insert(name.to_string(), (row, at)).is_some() {
                return Err(ReadError::whole(format!(
                    "{name:?} names more than one gene, so a homology naming it does not say which"
                )));
            }
        }
    }
    Ok(index)
}

/// The rows of a homology file, with the identity column still in its own unit.
fn rows(text: &str) -> Result<Vec<Row>, ReadError> {
    let mut rows = Vec::new();
    for (at, line) in lines(text) {
        let cols = columns(line);
        // Twelve is BLAST tabular. Three is a pair and a number, two is a pair.
        // Anything between them is a file that is not one of those, and reading
        // column three of it as an identity is how a length becomes 100%.
        let identity = match cols.len() {
            2 => None,
            3 => Some(cols[2]),
            n if n >= 12 => Some(cols[2]),
            n => {
                return Err(ReadError::at(
                    at,
                    format!(
                        "a homology row is 2 or 3 columns, or the 12 of BLAST tabular, and this one has {n}"
                    ),
                ))
            }
        };
        let identity = match identity {
            // A search that had nothing to say writes it several ways, and all
            // of them mean the same thing: no number was reported.
            Some(word) if matches!(word.trim(), "" | "." | "NA" | "na" | "N/A" | "*") => None,
            Some(word) => {
                let value: f64 = number(word, "identity", at)?;
                if !value.is_finite() {
                    return Err(ReadError::at(
                        at,
                        format!("identity {word:?} is not a number a search reported"),
                    ));
                }
                Some(value)
            }
            None => None,
        };
        rows.push(Row {
            query: cols[0].to_string(),
            subject: cols[1].to_string(),
            identity,
            line: at,
        });
    }
    Ok(rows)
}

/// Works out whether an identity column is a percentage, or refuses to.
///
/// A value above one can only be a percentage. A file whose values are all at
/// or below one could be either, since 0.95 is a fraction and 0.95% is a
/// number a search can report, so it is refused by name rather than guessed at.
fn infer(rows: &[Row]) -> Result<Option<Identity>, ReadError> {
    let stated: Vec<f64> = rows.iter().filter_map(|row| row.identity).collect();
    if stated.is_empty() {
        return Ok(None);
    }
    if stated.iter().any(|value| *value > 1.0) {
        return Ok(Some(Identity::Percent));
    }
    Err(ReadError::whole(
        "every identity in this file is at or below 1, which is a fraction or a very small \
         percentage, so --identity has to say which",
    ))
}

/// Turns an identity into the fraction a [`Homology`] takes.
fn scale(value: f64, unit: Option<Identity>, line: usize) -> Result<f64, ReadError> {
    // `unit` is `None` only where no row stated an identity, and this is
    // reached from a row that did.
    let unit = unit.unwrap_or(Identity::Fraction);
    let ceiling = unit.ceiling();
    if value < 0.0 || value > ceiling {
        return Err(ReadError::at(
            line,
            format!("identity {value} is outside 0 to {ceiling}"),
        ));
    }
    Ok(value / ceiling)
}

#[cfg(test)]
mod tests {
    use super::*;

    const BED: &str = "\
H37Rv\t0\t1200\tespA\t0\t+
H37Rv\t1300\t2100\tespC\t0\t+
H37Rv\t2200\t3000\tespD\t0\t-
CDC1551\t0\t1200\tespA2\t0\t+
CDC1551\t2200\t3000\tespD2\t0\t-
";

    fn window() -> Region {
        Region::new("ESX-1", 0, 3_000).unwrap()
    }

    fn read() -> Vec<Locus> {
        loci(BED, &window(), None).unwrap().loci
    }

    #[test]
    fn column_one_names_the_row_rather_than_selecting_it() {
        // Every other interval reader here drops a row whose first column is
        // not the sequence being drawn. This one draws them all and takes that
        // column as the name of the genome.
        let found = loci(BED, &window(), None).unwrap();
        assert_eq!(found.loci.len(), 2);
        assert_eq!(found.loci[0].name, "H37Rv");
        assert_eq!(found.loci[1].name, "CDC1551");
        assert_eq!(found.loci[0].genes.len(), 3);
        assert_eq!(found.loci[1].genes.len(), 2);
        assert_eq!(found.records, 5);
    }

    #[test]
    fn rows_and_genomes_keep_their_file_order_because_a_link_names_a_row_by_number() {
        let reversed = "\
Zed\t0\t100\tg1\t0\t+
Alpha\t0\t100\tg2\t0\t+
";
        let found = loci(reversed, &window(), None).unwrap();
        assert_eq!(found.loci[0].name, "Zed", "the rows were sorted");
        assert_eq!(found.loci[1].name, "Alpha");
    }

    #[test]
    fn a_bed_passes_its_coordinates_through_and_a_gff3_moves_its_start() {
        let bed = loci("g\t100\t200\tx\n", &window(), None).unwrap();
        assert_eq!(
            (bed.loci[0].genes[0].start, bed.loci[0].genes[0].end),
            (100, 200)
        );

        let gff = "##gff-version 3\ng\t.\tgene\t100\t200\t.\t+\t.\tName=x\n";
        let found = loci(gff, &window(), None).unwrap();
        assert_eq!(
            (found.loci[0].genes[0].start, found.loci[0].genes[0].end),
            (99, 200),
            "the same two numbers in the other convention are one base apart"
        );
    }

    #[test]
    fn genes_outside_the_window_are_counted_rather_than_dropped_in_silence() {
        let found = loci(BED, &Region::new("ESX-1", 2_500, 3_000).unwrap(), None).unwrap();
        assert_eq!(found.records, 5);
        assert_eq!(found.off_region, 3);
        assert_eq!(found.loci.len(), 2);
    }

    #[test]
    fn a_file_with_nothing_in_it_is_no_loci_rather_than_an_error() {
        let found = loci("# nothing here\n", &window(), None).unwrap();
        assert!(found.loci.is_empty());
        assert_eq!(found.records, 0);
    }

    // -----------------------------------------------------------------------
    // The join
    // -----------------------------------------------------------------------

    #[test]
    fn a_percent_identity_is_divided_by_a_hundred() {
        // Left as it was written it goes through Homology::new, which clamps,
        // so 78.4 and 99.9 both become a perfect match and nothing fails.
        let text = "espA\tespA2\t78.4\t100\t0\t0\t1\t100\t1\t100\t1e-50\t150\n";
        let found = links(text, &read(), None).unwrap();
        assert_eq!(found.links.len(), 1);
        assert_eq!(found.links[0].identity, Some(0.784));
    }

    #[test]
    fn an_identity_column_that_could_be_either_unit_stops_the_read_and_names_the_flag() {
        let text = "espA\tespA2\t0.98\n";
        let error = links(text, &read(), None).unwrap_err();
        assert!(error.reason.contains("--identity"), "{error}");

        // Told which it is, the same file reads.
        let found = links(text, &read(), Some(Identity::Fraction)).unwrap();
        assert_eq!(found.links[0].identity, Some(0.98));
        // And read as a percentage it is a different number, not an error.
        let found = links(text, &read(), Some(Identity::Percent)).unwrap();
        assert_eq!(found.links[0].identity, Some(0.0098));
    }

    #[test]
    fn an_identity_outside_its_unit_is_refused_and_not_clamped() {
        let text = "espA\tespA2\t140.0\n";
        let error = links(text, &read(), Some(Identity::Percent)).unwrap_err();
        assert!(error.reason.contains("outside 0 to 100"), "{error}");

        let text = "espA\tespA2\t1.4\n";
        let error = links(text, &read(), Some(Identity::Fraction)).unwrap_err();
        assert!(error.reason.contains("outside 0 to 1"), "{error}");
    }

    #[test]
    fn an_identity_nobody_stated_is_unstated_and_not_a_zero() {
        // A zero is the palest ribbon on the ramp, which is a weak match, which
        // is a claim about two genes nobody compared.
        for word in [".", "NA", "*", ""] {
            let text = format!("espA\tespA2\t{word}\n");
            let found = links(&text, &read(), None).unwrap();
            assert_eq!(found.links[0].identity, None, "{word:?}");
        }
        // And a bare pair file is every link unstated.
        let found = links("espA\tespA2\n", &read(), None).unwrap();
        assert_eq!(found.links[0].identity, None);
    }

    #[test]
    fn a_name_no_locus_has_is_listed_rather_than_resolved_to_a_gene() {
        // The failure this reader exists to stop. Resolved to index nought it
        // draws a ribbon into espA with a confident label and takes espA's
        // unmatched outline off, and nothing on the figure says so.
        let text = "lcl|NC_000962.3_cds_667\tespA2\t98.0\n";
        let found = links(text, &read(), None).unwrap();
        assert!(
            found.links.is_empty(),
            "a name that is not a gene was joined"
        );
        assert_eq!(found.unjoined, vec!["lcl|NC_000962.3_cds_667".to_string()]);
        assert_eq!(found.records, 1);
    }

    #[test]
    fn the_join_is_exact_so_a_name_that_only_looks_right_is_listed() {
        for name in ["ESPA", "espA ", " espA", "espa"] {
            let text = format!("{name}\tespA2\t98.0\n");
            let found = links(&text, &read(), None).unwrap();
            assert!(
                found.links.is_empty() && !found.unjoined.is_empty(),
                "{name:?} joined when it is not the gene's name"
            );
        }
    }

    #[test]
    fn a_name_that_is_more_than_one_gene_stops_the_read_and_says_which() {
        let twice = "\
A\t0\t100\tsame\t0\t+
B\t0\t100\tsame\t0\t+
";
        let loci = loci(twice, &window(), None).unwrap().loci;
        let error = links("same\tsame\t98.0\n", &loci, None).unwrap_err();
        assert!(error.reason.contains("more than one gene"), "{error}");
    }

    #[test]
    fn a_gene_with_no_name_is_kept_and_drawn_and_cannot_be_joined() {
        let unnamed = "A\t0\t100\nB\t0\t100\tg\t0\t+\n";
        let found = loci(unnamed, &window(), None).unwrap();
        assert_eq!(found.loci[0].genes.len(), 1);
        assert_eq!(found.loci[0].genes[0].name, None);
        // Two unnamed genes are not one name that is two genes, so the read
        // does not stop; they simply cannot be named by a homology file.
        let both = "A\t0\t100\nB\t0\t100\n";
        let loci = loci(both, &window(), None).unwrap().loci;
        assert!(links("x\ty\n", &loci, None).is_ok());
    }

    #[test]
    fn a_hit_between_two_genes_of_one_locus_is_counted_and_not_drawn() {
        // An all-against-all search reports every gene against itself, and a
        // paralogue pair inside one genome looks exactly the same.
        let text = "espA\tespA\t100.0\nespA\tespC\t72.0\n";
        let found = links(text, &read(), None).unwrap();
        assert!(found.links.is_empty());
        assert_eq!(found.self_hits, 2);
        // And counted as what they are. A hit inside one genome and a hit
        // across two rows that are not neighbours are different mistakes with
        // the same outcome, and only one of them means the rows are in the
        // wrong order, so a reader told the wrong one looks in the wrong place.
        assert_eq!(
            found.not_adjacent, 0,
            "a hit inside one row was blamed on the order of the rows"
        );
    }

    #[test]
    fn a_hit_between_loci_that_are_not_neighbours_is_counted_and_not_drawn() {
        let three = "\
A\t0\t100\ta1\t0\t+
B\t0\t100\tb1\t0\t+
C\t0\t100\tc1\t0\t+
";
        let loci = loci(three, &window(), None).unwrap().loci;
        let found = links("a1\tc1\t98.0\n", &loci, None).unwrap();
        assert!(found.links.is_empty(), "a ribbon skipped a row");
        assert_eq!(found.not_adjacent, 1);
    }

    #[test]
    fn a_reciprocal_pair_becomes_one_link_whichever_way_round_it_is_written() {
        // Two rows painting one ribbon means the last drawn is the one seen, so
        // the identity on the page would follow the file's order.
        let forward = links("espA\tespA2\t95.0\nespA2\tespA\t75.0\n", &read(), None).unwrap();
        let backward = links("espA2\tespA\t95.0\nespA\tespA2\t75.0\n", &read(), None).unwrap();
        assert_eq!(forward.links.len(), 1);
        assert_eq!(backward.links.len(), 1);
        assert_eq!(forward.merged, 1);
        // The first row of a pair is the one kept, and both files put the
        // same pair first, so both say 95 and neither says 75.
        assert_eq!(forward.links[0].identity, Some(0.95));
        assert_eq!(backward.links[0].identity, Some(0.95));
        assert_eq!(forward.links[0].row, backward.links[0].row);
        assert_eq!(forward.links[0].from, backward.links[0].from);
        assert_eq!(forward.links[0].to, backward.links[0].to);
    }

    #[test]
    fn several_hsps_of_one_pair_become_one_link() {
        let text = "\
espA\tespA2\t98.0\t400\t0\t0\t1\t400\t1\t400\t1e-99\t700
espA\tespA2\t91.0\t120\t0\t0\t500\t620\t500\t620\t1e-30\t200
espA\tespA2\t88.0\t90\t0\t0\t700\t790\t700\t790\t1e-20\t150
";
        let found = links(text, &read(), None).unwrap();
        assert_eq!(found.links.len(), 1);
        assert_eq!(found.merged, 2);
        assert_eq!(found.records, 3);
        assert_eq!(found.links[0].identity, Some(0.98));
    }

    #[test]
    fn a_links_file_carries_no_coordinate_and_this_reader_reads_none() {
        // A Homology names two genes and holds no span, so the eight coordinate
        // columns of BLAST tabular must not reach the figure by any route.
        let one = "espA\tespA2\t98.0\t400\t0\t0\t1\t400\t1\t400\t1e-99\t700\n";
        let other = "espA\tespA2\t98.0\t400\t0\t0\t9\t99\t2\t7\t1e-99\t700\n";
        assert_eq!(
            links(one, &read(), None).unwrap(),
            links(other, &read(), None).unwrap()
        );
    }

    #[test]
    fn the_comment_lines_of_outfmt_seven_read_the_same_as_outfmt_six() {
        let seven = "\
# BLASTN 2.16.0+
# Query: espA
# 1 hits found
espA\tespA2\t98.0\t400\t0\t0\t1\t400\t1\t400\t1e-99\t700
";
        let six = "espA\tespA2\t98.0\t400\t0\t0\t1\t400\t1\t400\t1e-99\t700\n";
        assert_eq!(
            links(seven, &read(), None).unwrap(),
            links(six, &read(), None).unwrap()
        );
    }

    #[test]
    fn a_row_that_is_not_a_homology_row_stops_the_read_and_says_which() {
        let five = "espA\tespA2\t98.0\t400\t0\n";
        let error = links(five, &read(), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("2 or 3 columns"), "{error}");
    }

    #[test]
    fn a_link_resolved_at_both_ends_names_the_genes_it_meant() {
        let text = "espD\tespD2\t99.1\n";
        let found = links(text, &read(), Some(Identity::Percent)).unwrap();
        // espD is gene 2 of row 0, espD2 is gene 1 of row 1.
        assert_eq!(found.links[0].row, 0);
        assert_eq!(found.links[0].from, 2);
        assert_eq!(found.links[0].to, 1);
        assert!(found.unjoined.is_empty());
        assert_eq!(found.self_hits, 0);
        assert_eq!(found.not_adjacent, 0);

        // And the track agrees: the genes it joined are the ones no longer
        // outlined as having nothing to match.
        let track = crate::LocusTrack::new(read()).links(found.links);
        assert_eq!(track.unmatched(0), vec![0, 1]);
        assert_eq!(track.unmatched(1), vec![0]);
    }
}
