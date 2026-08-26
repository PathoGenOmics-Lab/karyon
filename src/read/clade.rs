//! Stretches of reference carried by a named set of taxa.
//!
//! A [`CladeBlock`] is a span plus the taxa that carry it, and one file shape
//! in ordinary use holds both in one record: the GFF3 that recombination
//! detection writes, `Gubbins` above all, whose ninth column carries the taxa
//! as an attribute. Every other route to the same claim is two files that have
//! to be joined by hand.
//!
//! ```text
//! SEQUENCE  GUBBINS  CDS  1531  1670  0.000  .  0  node="N7";taxa="s1 s2 s3";
//! ```
//!
//! # 1-based inclusive in, 0-based half-open out
//!
//! GFF3 counts from one and includes its end, so a start comes back one lower
//! and an end stays where it is, having already been one past the last base
//! once the count began at zero. That is the whole of the conversion, and it is
//! the same one [`interval`](super::interval) does for a gene.
//!
//! # A list is split before it is decoded
//!
//! The ninth column spends `;`, `=` and `,` on its own syntax, so a name
//! holding one of them arrives percent encoded. Decoding the whole value and
//! then splitting it turns `taxa=A%2CB,C` into three names where the file said
//! two, and neither of the wrong two is in any tree, so the block reports as
//! carried by nobody. The value is therefore split first and each name decoded
//! after. A taxon whose name holds a space or a comma cannot be written in this
//! format at all, which is worth knowing and is not worth guessing about.
//!
//! # The sequence a clade file names is usually not one
//!
//! `Gubbins` writes the literal `SEQUENCE` in column one whatever the reference
//! was called. Filtering on the region's own sequence name the way a feature
//! reader does therefore keeps nothing, and keeping nothing is not an error
//! anywhere: it draws the phylogeny with no blocks on it, which reads as a
//! confident statement that there was no recombination here.
//!
//! So the filter is conditional. A file naming exactly one sequence is about
//! that sequence whatever it is called, and [`Blocks::sequences`] reports the
//! name so a caller can say which. A file naming several is a whole genome, and
//! then the region's own name selects among them and the rest are counted.
//!
//! # What is refused, and what is counted
//!
//! A row that will not parse stops the read on its line. So does a block with
//! no taxa: an empty carrier list is not a block carried by nobody, it is a
//! file that does not hold what the flag claimed, and the difference is the
//! whole of what the figure would be saying.
//!
//! Everything else that does not reach the blocks is counted rather than
//! dropped in silence, because a caller drawing four blocks out of four hundred
//! should be able to say so, and because an empty answer and an empty file look
//! identical once they are drawn.

use std::collections::BTreeSet;

use crate::{CladeBlock, Region};

use super::interval::{percent_decode, raw_attribute};
use super::{columns, lines, number, ReadError};

/// The blocks a clade file holds, and what it held that is not in them.
#[derive(Debug, Clone, PartialEq)]
pub struct Blocks {
    /// The blocks, in the order the file listed them.
    pub blocks: Vec<CladeBlock>,
    /// Block records in the file, before any filter.
    ///
    /// The field that separates a file saying nothing from a file that was not
    /// read. Both give no blocks, and only one of them is worth drawing.
    pub records: usize,
    /// The sequences the file named, in the order they were first seen.
    pub sequences: Vec<String>,
    /// Records on this sequence that touch no base of the window.
    pub off_region: usize,
    /// Records naming another sequence, where the file named more than one.
    pub passed_over: usize,
    /// Taxon names dropped as repeats of one already in their own block.
    ///
    /// Not a detail. A block listing a taxon twice claims more carriers than
    /// the tree has rows, and the count of rows a clade covers is then taken
    /// from a subtraction that has no answer.
    pub duplicate_taxa: usize,
}

/// Reads clade blocks out of a GFF3 carrying a `taxa` attribute.
///
/// GFF3 is 1-based and inclusive, so a start comes back one lower and an end
/// is unchanged. The taxa are the `taxa` attribute, split on whitespace and on
/// commas, each name percent decoded afterwards. A `node` attribute becomes the
/// block's name, and `Name` and `ID` are accepted after it.
///
/// Rows naming another sequence than `region.seq()` are skipped, but only where
/// the file names more than one sequence: a file about a single sequence is
/// about it whatever column one calls it.
///
/// # Errors
///
/// Returns the first row that is not a GFF3 row, whose span is inverted or
/// empty, or that carries no taxa. A file that parses and holds no block on
/// this sequence is not an error here: [`Blocks::records`] says which of the
/// two happened, and the caller decides.
pub fn blocks(text: &str, region: &Region) -> Result<Blocks, ReadError> {
    let named = sequences(text)?;
    // One sequence in the file is the sequence the file is about. More than one
    // is a whole genome, and then the region picks.
    let only = (named.len() == 1).then(|| named[0].0.clone());

    let mut found = Blocks {
        blocks: Vec::new(),
        records: 0,
        sequences: named.iter().map(|(name, _)| name.clone()).collect(),
        off_region: 0,
        passed_over: 0,
        duplicate_taxa: 0,
    };

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 9 {
            return Err(ReadError::at(
                at,
                format!("a GFF3 row has nine columns, this one has {}", cols.len()),
            ));
        }
        found.records += 1;
        if only.is_none() && cols[0] != region.seq() {
            found.passed_over += 1;
            continue;
        }

        let start: u64 = number(cols[3], "start", at)?;
        let end: u64 = number(cols[4], "end", at)?;
        if start == 0 {
            return Err(ReadError::at(at, "GFF3 counts from 1, so 0 is not a start"));
        }
        if end < start {
            return Err(ReadError::at(at, "end is before start"));
        }
        // A GFF3 span is inclusive at both ends, so the shortest one anybody
        // can write is a single base and `end == start - 1` never happens. The
        // conversion below turns `start..end` into `start - 1..end`, so an
        // empty block cannot arrive and is not guarded against here.
        let (start, end) = (start - 1, end);

        let taxa = taxa(cols[8], at, &mut found.duplicate_taxa)?;
        if start >= region.end() || end <= region.start() {
            found.off_region += 1;
            continue;
        }

        let mut block = CladeBlock::new(start, end, taxa);
        if let Some(name) = ["node", "Name", "ID"]
            .iter()
            .find_map(|key| raw_attribute(cols[8], key).map(unquote))
            .map(|value| percent_decode(&value))
            .filter(|value| !value.is_empty() && value != ".")
        {
            block = block.name(name);
        }
        found.blocks.push(block);
    }

    Ok(found)
}

/// The sequences a clade file names, most records first, name to settle a tie.
///
/// The question that comes before drawing anything, because the name in column
/// one of this format is often not a sequence name at all. Ordered, so the same
/// file answers the same way twice.
pub fn sequences(text: &str) -> Result<Vec<(String, usize)>, ReadError> {
    let mut seen: Vec<(String, usize)> = Vec::new();
    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 9 {
            return Err(ReadError::at(
                at,
                format!("a GFF3 row has nine columns, this one has {}", cols.len()),
            ));
        }
        match seen.iter_mut().find(|(name, _)| name == cols[0]) {
            Some((_, count)) => *count += 1,
            None => seen.push((cols[0].to_string(), 1)),
        }
    }
    seen.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    Ok(seen)
}

/// The carriers named in one attribute column.
///
/// Split first and decoded after, so a name holding an escaped separator stays
/// one name. Repeats are dropped and counted: a block claiming one taxon twice
/// claims more carriers than the tree has leaves, and what a clade covers is
/// then worked out by taking one number from a smaller one.
fn taxa(attributes: &str, at: usize, duplicates: &mut usize) -> Result<Vec<String>, ReadError> {
    let Some(raw) = raw_attribute(attributes, "taxa") else {
        return Err(ReadError::at(
            at,
            "a clade block names the taxa carrying it, and this row has no taxa attribute",
        ));
    };
    let raw = unquote(raw);

    let mut seen = BTreeSet::new();
    let mut taxa = Vec::new();
    for word in raw.split([' ', '\t', ',']) {
        let name = percent_decode(word);
        let name = name.trim();
        // A blank name is not a taxon. Left in, it matches a leaf with no name
        // of its own, and an unlabelled tip would then be drawn as a carrier of
        // a clade nobody said it was in.
        if name.is_empty() {
            continue;
        }
        if seen.insert(name.to_string()) {
            taxa.push(name.to_string());
        } else {
            *duplicates += 1;
        }
    }

    if taxa.is_empty() {
        return Err(ReadError::at(
            at,
            "a clade block names the taxa carrying it, and this row names none",
        ));
    }
    Ok(taxa)
}

/// Strips the quotes a GTF-style attribute value is written inside.
///
/// GFF3 leaves values bare and the recombination tools quote them, so both
/// spellings reach here and only one of them can be split into names.
fn unquote(value: &str) -> String {
    value.trim().trim_matches('"').trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two blocks on one sequence, as `Gubbins` writes them.
    const GFF: &str = "\
##gff-version 3
SEQUENCE\tGUBBINS\tCDS\t100\t200\t0.000\t.\t0\tnode=\"N7\";taxa=\"s1 s2 s3\";snp_count=\"8\";
SEQUENCE\tGUBBINS\tCDS\t1531\t1670\t0.000\t.\t0\tnode=\"N2\";taxa=\"s4 s5\";snp_count=\"3\";
";

    fn window() -> Region {
        Region::new("NC_011900.1", 0, 5_000).unwrap()
    }

    #[test]
    fn a_start_moves_back_one_and_an_end_does_not() {
        let found = blocks(GFF, &window()).unwrap();
        assert_eq!(found.blocks.len(), 2);
        // 100..200 written 1-based inclusive is 99..200 counted from zero, and
        // that is 101 bases, one more than the two numbers suggest.
        assert_eq!(found.blocks[0].start(), 99);
        assert_eq!(found.blocks[0].end(), 200);
        assert_eq!(found.blocks[0].end() - found.blocks[0].start(), 101);
    }

    #[test]
    fn the_taxa_attribute_is_read_and_the_node_names_the_block() {
        let found = blocks(GFF, &window()).unwrap();
        assert_eq!(found.blocks[1].taxa(), ["s4", "s5"]);
        // The whole block, so the name the node attribute carries is pinned
        // too. CladeBlock keeps no getter for it and does keep an equality.
        assert_eq!(
            found.blocks[0],
            CladeBlock::new(99, 200, ["s1", "s2", "s3"]).name("N7")
        );
    }

    #[test]
    fn a_list_is_split_before_it_is_decoded() {
        // The whole reason this reader does not use the decoding accessor. The
        // file names two taxa, one of which holds a comma. Decoding first gives
        // three names, none of which anybody has.
        let text = "chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"A%2CB,C\";\n";
        let found = blocks(text, &Region::new("chr", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks[0].taxa(), ["A,B", "C"]);

        // And a space that arrived escaped stays inside its name too.
        let spaced = "chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"sample%20one sample%20two\";\n";
        let found = blocks(spaced, &Region::new("chr", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks[0].taxa(), ["sample one", "sample two"]);
    }

    #[test]
    fn runs_of_whitespace_do_not_become_taxa_that_are_not_there() {
        let text = "chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"  s3   s4  s2 \";\n";
        let found = blocks(text, &Region::new("chr", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks[0].taxa(), ["s3", "s4", "s2"]);
    }

    #[test]
    fn a_repeated_taxon_is_dropped_and_counted() {
        // Left in, the block claims four carriers where the tree has three
        // rows, and what the clade covers is then a subtraction with no answer.
        let text = "chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"s1 s2 s1 s3\";\n";
        let found = blocks(text, &Region::new("chr", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks[0].taxa(), ["s1", "s2", "s3"]);
        assert_eq!(found.duplicate_taxa, 1);
    }

    #[test]
    fn a_block_naming_no_taxa_stops_the_read_rather_than_carrying_none() {
        for attributes in [
            "node=\"N1\";",
            "taxa=\"\";",
            "taxa=\"   \";",
            "taxa=\",,\";",
        ] {
            let text = format!("chr\t.\tCDS\t1\t9\t.\t.\t0\t{attributes}\n");
            let error = blocks(&text, &Region::new("chr", 0, 100).unwrap()).unwrap_err();
            assert_eq!(error.line, 1, "{attributes}");
            assert!(error.reason.contains("taxa"), "{attributes}: {error}");
        }
    }

    #[test]
    fn a_blank_taxon_is_never_matched_against_a_leaf_with_no_name() {
        // `taxa="s1  s2"` holds two names and one gap. A gap kept as a name
        // matches an unlabelled tip, and that tip is then drawn as a carrier.
        let text = "chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"s1 , s2\";\n";
        let found = blocks(text, &Region::new("chr", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks[0].taxa(), ["s1", "s2"]);
        assert!(found.blocks[0].taxa().iter().all(|t| !t.is_empty()));
    }

    #[test]
    fn an_inverted_or_zero_span_stops_the_read() {
        let region = Region::new("chr", 0, 1_000).unwrap();
        // CladeBlock::new swaps an inverted span silently, so a 3 kb block
        // would be conjured out of two numbers written the wrong way round.
        let inverted = "chr\t.\tCDS\t400\t100\t.\t.\t0\ttaxa=\"s1\";\n";
        let error = blocks(inverted, &region).unwrap_err();
        assert!(error.reason.contains("end is before start"), "{error}");

        let zero = "chr\t.\tCDS\t0\t100\t.\t.\t0\ttaxa=\"s1\";\n";
        let error = blocks(zero, &region).unwrap_err();
        assert!(error.reason.contains("counts from 1"), "{error}");
    }

    #[test]
    fn a_file_holding_nothing_is_counted_rather_than_answered_with_a_bare_vector() {
        // A tree drawn with no blocks on it is a figure that says there was no
        // recombination here. It is not the figure a file with nothing in it
        // should produce, and the count is what tells the two apart.
        let empty = blocks("##gff-version 3\n", &window()).unwrap();
        assert_eq!(empty.records, 0);
        assert!(empty.blocks.is_empty());

        let culled = blocks(GFF, &Region::new("SEQUENCE", 4_000, 5_000).unwrap()).unwrap();
        assert_eq!(culled.records, 2, "the file did hold blocks");
        assert_eq!(culled.off_region, 2);
        assert!(culled.blocks.is_empty());
    }

    #[test]
    fn a_file_naming_one_sequence_is_read_whatever_that_sequence_is_called() {
        // Gubbins writes the literal SEQUENCE, and a region is named after the
        // accession. Filtering on the name keeps nothing and says nothing.
        let found = blocks(GFF, &window()).unwrap();
        assert_eq!(found.blocks.len(), 2);
        assert_eq!(found.sequences, vec!["SEQUENCE".to_string()]);
        assert_eq!(found.passed_over, 0);
    }

    #[test]
    fn a_file_naming_several_sequences_filters_on_the_region_and_counts_the_rest() {
        let text = "\
chrA\t.\tCDS\t10\t20\t.\t.\t0\ttaxa=\"s1\";
chrB\t.\tCDS\t10\t20\t.\t.\t0\ttaxa=\"s2\";
chrA\t.\tCDS\t30\t40\t.\t.\t0\ttaxa=\"s3\";
";
        let found = blocks(text, &Region::new("chrA", 0, 100).unwrap()).unwrap();
        assert_eq!(found.blocks.len(), 2);
        assert_eq!(found.passed_over, 1);
        assert_eq!(found.records, 3);

        assert_eq!(
            sequences(text).unwrap(),
            vec![("chrA".to_string(), 2), ("chrB".to_string(), 1)]
        );
    }

    #[test]
    fn a_line_that_is_not_a_gff3_row_stops_the_read_and_says_which() {
        let short = "chr\t.\tCDS\t1\t9\n";
        let error = blocks(short, &window()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("nine columns"), "{error}");

        let text = "\
chr\t.\tCDS\t1\t9\t.\t.\t0\ttaxa=\"s1\";
chr\t.\tCDS\tx\t9\t.\t.\t0\ttaxa=\"s1\";
";
        let error = blocks(text, &Region::new("chr", 0, 100).unwrap()).unwrap_err();
        assert_eq!(error.line, 2);
        assert!(error.reason.contains("start"), "{error}");
    }
}
