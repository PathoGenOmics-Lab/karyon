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

use std::collections::BTreeMap;

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
/// A gene comes back once, with the exons of all its transcripts and the
/// stretches any of them codes, and says how many transcripts it merged; a
/// BED12 row comes back with its blocks and its thick span. [`transcripts`]
/// reads each transcript as a feature of its own instead.
///
/// Rows on another sequence than `region.seq()` are skipped.
pub fn features(
    text: &str,
    region: &Region,
    format: Option<Format>,
) -> Result<Vec<Feature>, ReadError> {
    models(text, region, format, Level::Gene).map(drawn)
}

/// Reads gene models one transcript at a time: each isoform a feature of
/// its own, named for itself and naming the gene it belongs to.
///
/// Everything else is read as [`features`] reads it. A gene with no
/// transcripts under it, as a bacterial gene written over its CDS, is one
/// feature either way, and so is every BED row.
pub fn transcripts(
    text: &str,
    region: &Region,
    format: Option<Format>,
) -> Result<Vec<Feature>, ReadError> {
    models(text, region, format, Level::Transcript).map(drawn)
}

/// Reads each transcript as it codes, for a ruler that counts its codons: its
/// exons, and its CDS exactly as the file writes it, row by row.
///
/// [`transcripts`] reads the same models for drawing, and a drawing has no
/// use for a coding stretch from one end of a feature to the other: it looks
/// the same as none, so it is left out there. Here the two are different
/// answers, a gene that codes from end to end and one that codes nowhere, so
/// it is kept. A CDS row under nothing, as some annotations write each one
/// beside its gene, codes over itself, and over every other row under its
/// `ID=`, as GFF3 writes one feature in pieces. A BED row codes over its
/// thick span. Rows on another sequence than `region.seq()` are skipped.
///
/// A drawing joins the rows of a CDS that touch, and so does
/// [`Feature::coding`], which is right for a picture and wrong for a count:
/// two rows sharing a base, as NCBI writes a ribosomal slippage, are one
/// stretch on the page and two frames to a ruler. So each transcript comes
/// back with [`Cds::rows`] beside it, the rows apart, each with its phase.
///
/// ```
/// use karyon::read::interval::coding;
/// use karyon::Region;
///
/// let gff = "##gff-version 3\n\
///            chr1\t.\tgene\t101\t400\t.\t-\t.\tID=g1;Name=abcD\n\
///            chr1\t.\tCDS\t101\t400\t.\t-\t0\tParent=g1\n";
/// let found = coding(gff, &Region::parse("chr1:1-1000")?, None)?;
/// assert_eq!(found[0].feature.coding, [(100, 400)]);
/// assert_eq!(found[0].rows[0].phase, Some(0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn coding(text: &str, region: &Region, format: Option<Format>) -> Result<Vec<Cds>, ReadError> {
    Ok(models(text, region, format, Level::Coding)?
        .into_iter()
        .map(|(feature, rows)| Cds { feature, rows })
        .collect())
}

/// One transcript as [`coding`] reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct Cds {
    /// The transcript: its name, its strand, the gene it belongs to, its
    /// exons, and its CDS joined into the stretches a drawing shows.
    pub feature: Feature,
    /// The rows its CDS is written in, apart and in the order of their
    /// starts: a GFF3's CDS rows, a GTF's with its `stop_codon` rows where
    /// they lie outside them, and a BED's thick span. A GTF's `start_codon`
    /// rows lie inside the CDS and say nothing more, so they are left out,
    /// and so is a `stop_codon` row inside a CDS row, as AUGUSTUS writes one.
    /// Empty for a transcript that codes nowhere.
    pub rows: Vec<CdsRow>,
}

/// One row of a coding sequence as the file writes it, 0-based and
/// half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdsRow {
    /// First base, 0-based.
    pub start: u64,
    /// One past the last base.
    pub end: u64,
    /// Its phase, the eighth column of a GFF3 or a GTF: how many bases at
    /// its 5' end, its start on the forward strand and its end on the
    /// reverse, finish a codon begun before it, 0, 1 or 2. `None` where the
    /// row writes `.`, and for a BED, which has no such column.
    ///
    /// A CDS whose 5'-most row has a phase of 1 or 2 does not begin on a
    /// codon, so the annotation does not hold its start codon: NCBI writes
    /// one so at the edge of a contig.
    pub phase: Option<u8>,
    /// Whether the row says the CDS goes on past its 5' end, as NCBI writes
    /// `start_range=.,N` on the forward strand and `end_range=N,.` on the
    /// reverse, so its first codon is not the start codon whatever its
    /// phase.
    pub partial: bool,
}

impl CdsRow {
    /// The row as `cols` write it, with the span `feature` read from them.
    fn of(cols: &[&str], feature: &Feature) -> Self {
        let phase = match cols.get(7).map(|phase| phase.trim()) {
            Some("0") => Some(0),
            Some("1") => Some(1),
            Some("2") => Some(2),
            _ => None,
        };
        let attributes = cols.get(8).copied().unwrap_or_default();
        let range = |key: &str| raw_attribute(attributes, key).is_some();
        let partial = match feature.strand {
            Strand::Forward => range("start_range"),
            Strand::Reverse => range("end_range"),
            Strand::Unknown => range("start_range") || range("end_range"),
        };
        CdsRow {
            start: feature.start,
            end: feature.end,
            phase,
            partial,
        }
    }
}

/// The NCBI translation table the CDS rows over `start..end` of `sequence`
/// name with `transl_table=`, as NCBI writes each one, 0-based and half-open:
/// the first such row's, or `None` where none of them names one.
///
/// A BED names no table, and neither does a GTF, so both give `None`.
///
/// ```
/// use karyon::read::interval::translation_table;
///
/// let gff = "##gff-version 3\n\
///            chrM\t.\tCDS\t3307\t4262\t.\t+\t0\tID=cds-ND1;transl_table=2\n";
/// assert_eq!(translation_table(gff, "chrM", 3_400, 3_500), Some(2));
/// assert_eq!(translation_table(gff, "chrM", 5_000, 5_100), None);
/// ```
pub fn translation_table(text: &str, sequence: &str, start: u64, end: u64) -> Option<u8> {
    if flavour(text, None) == Flavour::Bed {
        return None;
    }
    lines(text).find_map(|(_, line)| {
        let cols = columns(line);
        if cols.first().copied().unwrap_or_default() != sequence
            || part(cols.get(2)?.trim()) != Some(Part::Coding)
        {
            return None;
        }
        let from = cols.get(3)?.trim().parse::<u64>().ok()?.checked_sub(1)?;
        let to = cols.get(4)?.trim().parse::<u64>().ok()?;
        if to <= start || from >= end {
            return None;
        }
        attribute(cols.get(8)?, "transl_table")?.parse().ok()
    })
}

/// Whether a CDS row of `text` on `sequence` is under nothing and has an
/// `ID=` of its own, which [`coding`] reads as one piece of every row under
/// that ID.
///
/// The rows of a gene lie under the gene's own row, so a window widened by
/// [`reach`] holds every one of them; the pieces of a CDS under nothing lie
/// under no row, and a window through an index holds the pieces over it
/// alone, which counted on their own are a CDS cut short. A GTF row always
/// names its transcript, so it is never one of these.
pub(crate) fn pieced(text: &str, sequence: &str) -> bool {
    lines(text).any(|(_, line)| {
        let cols = columns(line);
        cols.first() == Some(&sequence)
            && cols
                .get(2)
                .is_some_and(|kind| part(kind.trim()) == Some(Part::Coding))
            && parents_of(&cols).is_empty()
            && declared_as(&cols).is_some()
    })
}

/// Whether a gene is drawn once or each of its transcripts is, or each
/// transcript is read for what it codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Level {
    Gene,
    Transcript,
    Coding,
}

/// The models over `region` at `level`, each with the rows its CDS is
/// written in, which only [`Level::Coding`] reads and the others leave empty.
fn models(
    text: &str,
    region: &Region,
    format: Option<Format>,
    level: Level,
) -> Result<Vec<(Feature, Vec<CdsRow>)>, ReadError> {
    let flavour = flavour(text, format);
    // The sequence is checked before the rest of the line is understood, so
    // that a whole genome annotation costs one comparison per row it does not
    // need. It also means the FASTA section some GFF3 files carry at the end
    // goes past without being read as a broken feature.
    let rows = lines(text)
        .map(|(at, line)| (at, columns(line)))
        .filter(|(_, cols)| cols.first().copied().unwrap_or_default() == region.seq());

    let mut features = Vec::new();
    if flavour == Flavour::Bed {
        for (at, cols) in rows {
            if level == Level::Coding {
                let feature = bed_as_written(&cols, at)?;
                let rows = feature
                    .coding
                    .iter()
                    .map(|&(start, end)| CdsRow {
                        start,
                        end,
                        phase: None,
                        partial: false,
                    })
                    .collect();
                keep((feature, rows), region, &mut features);
            } else {
                keep((bed(&cols, at)?, Vec::new()), region, &mut features);
            }
        }
        return Ok(features);
    }

    // A GFF3 or GTF annotation writes a gene once for every level of it, the
    // gene, each transcript, each exon and each CDS, and a feature track drew
    // every one as a feature of its own: five rows for one gene, the gene
    // beside its own CDS. So the levels are put back together: the exons and
    // the CDS of a transcript become its structure, and the transcripts of a
    // gene become the gene, or each a feature of its own.
    let rows: Vec<(usize, Vec<&str>)> = rows.collect();
    for model in assemble(&rows, level)? {
        keep(model, region, &mut features);
    }
    Ok(features)
}

/// The features of [`models`], without the rows of their CDS.
fn drawn(models: Vec<(Feature, Vec<CdsRow>)>) -> Vec<Feature> {
    models.into_iter().map(|(feature, _)| feature).collect()
}

/// What a row is to the transcript it belongs to, for the rows that are
/// pieces of one rather than things of their own.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Part {
    /// An exon, whatever its kind: `exon`, `noncoding_exon`, `pseudogenic_exon`.
    Exon,
    /// What codes: a CDS, and the start and stop codons GTF writes apart
    /// from it, the stop outside it.
    Coding,
    /// An untranslated region, which is part of an exon.
    Untranslated,
}

/// The part a row of this type is, if it is one.
fn part(kind: &str) -> Option<Part> {
    let kind = kind.to_ascii_lowercase();
    match kind.as_str() {
        "exon" => Some(Part::Exon),
        "cds" | "start_codon" | "stop_codon" => Some(Part::Coding),
        "five_prime_utr" | "three_prime_utr" | "utr" | "5utr" | "3utr" | "utr5" | "utr3" => {
            Some(Part::Untranslated)
        }
        _ if kind.ends_with("_exon") => Some(Part::Exon),
        _ => None,
    }
}

/// Whether a row of this type is a gene, and so not a transcript of anything:
/// a gene under an operon makes the operon a cluster of genes, which has no
/// exons of its own to draw.
fn gene_level(kind: &str) -> bool {
    let kind = kind.to_ascii_lowercase();
    kind == "gene" || kind == "pseudogene" || kind.ends_with("_gene")
}

/// What a node is known by: the name a row declares or a part names, or,
/// for a row that declares none, where it is in the file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Key<'a> {
    Named(&'static str, &'a str),
    Row(usize),
}

/// A gene, a transcript, or something else a row declares, with what the
/// file says is under it.
#[derive(Debug)]
struct Node<'a> {
    /// Its own row, read, or none for a parent only its parts name.
    row: Option<Feature>,
    /// Its row's type, empty for a parent with no row.
    kind: &'a str,
    /// The first row that mentions it, which is the order genes come back in.
    first: usize,
    /// The pieces under it: exons, coding stretches, untranslated ones.
    parts: Vec<(Part, u64, u64)>,
    /// Its coding pieces as their rows write them, for [`Level::Coding`],
    /// each with whether it is a `start_codon` or a `stop_codon` row, which a
    /// CDS row may hold already.
    rows: Vec<(bool, CdsRow)>,
    /// What is under it that is not a piece: the transcripts of a gene.
    children: Vec<Key<'a>>,
    /// What it is under.
    parents: Vec<Key<'a>>,
    /// For a parent with no row: where its members reach, their strand, and
    /// the names they give it.
    reach: Option<(u64, u64)>,
    strand: Strand,
    names: Vec<String>,
}

impl<'a> Node<'a> {
    fn new(first: usize) -> Self {
        Node {
            row: None,
            kind: "",
            first,
            parts: Vec::new(),
            rows: Vec::new(),
            children: Vec::new(),
            parents: Vec::new(),
            reach: None,
            strand: Strand::Unknown,
            names: Vec::new(),
        }
    }

    /// Takes in a member's span and strand, for a node with no row of its own.
    fn reaches(&mut self, start: u64, end: u64, strand: Strand) {
        self.reach = Some(match self.reach {
            Some((from, to)) => (from.min(start), to.max(end)),
            None => (start, end),
        });
        if self.strand == Strand::Unknown {
            self.strand = strand;
        }
    }

    /// The node as a feature with nothing under it yet: its own row, or, for
    /// a parent only its members name, their reach under the name they agree
    /// on, or the one they call it by.
    fn feature(&self, key: Key<'_>) -> Feature {
        if let Some(row) = &self.row {
            return row.clone();
        }
        let (start, end) = self.reach.unwrap_or((0, 1));
        let mut feature = Feature::new(start, end).strand(self.strand);
        let agreed = self
            .names
            .first()
            .filter(|first| self.names.iter().all(|name| name == *first));
        match (agreed, key) {
            (Some(name), _) => feature = feature.name(name.clone()),
            (None, Key::Named(_, id)) => feature = feature.name(percent_decode(id)),
            (None, Key::Row(_)) => {}
        }
        feature
    }
}

/// The exons and the coding stretches of one transcript, from its parts.
///
/// Exons come from its exon rows. A transcript written without any, as a
/// gene over its CDS alone, is made of its coding and untranslated pieces,
/// reaching out to its own ends, which are untranslated where nothing says
/// otherwise rather than introns: a gene row longer than its CDS is a gene
/// with ends that do not code.
/// The exons of a transcript and the stretches of it that code.
struct Shape {
    exons: Vec<(u64, u64)>,
    coding: Vec<(u64, u64)>,
    /// The rows of its CDS, as [`Cds::rows`] holds them.
    rows: Vec<CdsRow>,
}

fn structure(node: &Node<'_>, span: (u64, u64)) -> Shape {
    let of = |which: &[Part]| -> Vec<(u64, u64)> {
        node.parts
            .iter()
            .filter(|(part, _, _)| which.contains(part))
            .map(|&(_, start, end)| (start, end))
            .collect()
    };
    let coding = of(&[Part::Coding]);
    let mut exons = of(&[Part::Exon]);
    if exons.is_empty() {
        let mut pieces = of(&[Part::Coding, Part::Untranslated]);
        pieces.sort_unstable();
        for (start, end) in pieces {
            match exons.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => exons.push((start, end)),
            }
        }
        if let Some(first) = exons.first_mut() {
            first.0 = first.0.min(span.0);
        }
        if let Some(last) = exons.last_mut() {
            last.1 = last.1.max(span.1);
        }
    }
    // A start or stop codon written inside a CDS row adds nothing to it, and
    // kept, it would be a second row over the same bases, which is how a
    // ribosomal slippage is written.
    let mut rows: Vec<CdsRow> = node
        .rows
        .iter()
        .filter(|(codon, row)| {
            !codon
                || !node
                    .rows
                    .iter()
                    .any(|(other, cds)| !other && cds.start <= row.start && row.end <= cds.end)
        })
        .map(|&(_, row)| row)
        .collect();
    rows.sort_unstable_by_key(|row| (row.start, row.end));
    rows.dedup();
    Shape {
        exons,
        coding,
        rows,
    }
}

/// A feature with nothing to show in its structure is left as one piece: one
/// exon from end to end, and a coding stretch from end to end or none.
fn one_piece(mut feature: Feature) -> Feature {
    if feature.exons[..] == [(feature.start, feature.end)] {
        feature.exons.clear();
    }
    if feature.exons.is_empty() && feature.coding[..] == [(feature.start, feature.end)] {
        feature.coding.clear();
    }
    feature
}

/// Puts the levels of a GFF3 or GTF annotation back together into features.
///
/// A row that is a part (an exon, a CDS, an untranslated region) belongs to
/// the transcript it names, and a transcript to the gene it names, whether or
/// not the file has a row for either: a GTF of exons alone, as a table browser
/// writes it, names its transcript and its gene in every row, and those two
/// are what it draws. A GFF3 names parents with `Parent=` and `Derives_from=`,
/// a GTF with `transcript_id` and `gene_id`.
///
/// A gene is its own row with the exons of every transcript under it, or at
/// [`Level::Transcript`] each transcript is. Anything deeper than a transcript
/// that is not a part, a polypeptide or an intron row, is under it and not
/// drawn, as it was not before. A row under nothing is drawn as it stands.
fn assemble(
    rows: &[(usize, Vec<&str>)],
    level: Level,
) -> Result<Vec<(Feature, Vec<CdsRow>)>, ReadError> {
    let mut nodes: BTreeMap<Key<'_>, Node<'_>> = BTreeMap::new();
    for (index, (at, cols)) in rows.iter().enumerate() {
        // Read before anything else, so a broken row stops the file on its
        // line whether or not it would have been drawn.
        let feature = gff3(cols, *at)?;
        if describes_sequence(cols) {
            continue;
        }
        let kind = cols.get(2).map_or("", |kind| kind.trim());
        let parents: Vec<Key<'_>> = parents_of(cols)
            .into_iter()
            .map(|(space, name)| Key::Named(space, name))
            .collect();
        let attributes = cols.get(8).copied().unwrap_or_default();

        if let (Some(part), false) = (part(kind), parents.is_empty()) {
            // A part belongs to every transcript a GFF3 row names, and to the
            // transcript of a GTF row, or its gene where it names none.
            let owners: Vec<Key<'_>> =
                if parents.iter().any(|key| matches!(key, Key::Named("id", _))) {
                    parents.clone()
                } else {
                    parents.iter().take(1).copied().collect()
                };
            for owner in owners {
                let node = nodes.entry(owner).or_insert_with(|| Node::new(index));
                node.parts.push((part, feature.start, feature.end));
                if part == Part::Coding && level == Level::Coding {
                    let codon = !kind.eq_ignore_ascii_case("cds");
                    node.rows.push((codon, CdsRow::of(cols, &feature)));
                }
                if node.row.is_none() {
                    node.reaches(feature.start, feature.end, feature.strand);
                    // What the part calls its transcript, where it says: a
                    // GTF row always does, a GFF3 exon mostly does not.
                    let name = match owner {
                        Key::Named("transcript", id) => Some(
                            gtf_attribute(attributes, "transcript_name")
                                .unwrap_or(id)
                                .to_string(),
                        ),
                        Key::Named("gene", id) => Some(
                            gtf_attribute(attributes, "gene_name")
                                .unwrap_or(id)
                                .to_string(),
                        ),
                        _ => feature.name.clone(),
                    };
                    node.names.extend(name);
                    // A GTF transcript with no row still belongs to its gene.
                    if let Key::Named("transcript", _) = owner {
                        for parent in parents.iter().skip(1) {
                            if !node.parents.contains(parent) {
                                node.parents.push(*parent);
                            }
                        }
                    }
                }
            }
            continue;
        }

        let key =
            declared_as(cols).map_or(Key::Row(index), |(space, name)| Key::Named(space, name));
        let node = nodes.entry(key).or_insert_with(|| Node::new(index));
        node.first = node.first.min(index);
        node.kind = kind;
        match (level, part(kind), &mut node.row) {
            // A CDS row under nothing is a coding sequence of its own, which
            // only a reader of what codes has to say, and GFF3 writes one in
            // pieces as rows under one `ID=`: each a piece of it, where the
            // last row alone was a CDS cut short and counted from the wrong
            // base.
            (Level::Coding, Some(Part::Coding), row) => {
                node.parts.push((Part::Coding, feature.start, feature.end));
                let codon = !kind.eq_ignore_ascii_case("cds");
                node.rows.push((codon, CdsRow::of(cols, &feature)));
                match row {
                    Some(row) => {
                        row.start = row.start.min(feature.start);
                        row.end = row.end.max(feature.end);
                    }
                    None => *row = Some(feature),
                }
            }
            (_, _, row) => *row = Some(feature),
        }
        node.names.clear();
        node.reach = None;
        node.parents = parents;
    }

    // Each node under the ones it names. A GTF gene the file never writes is
    // made from its transcripts, since every GTF row names its gene and a
    // file without gene rows still has genes. A GFF3 parent the file never
    // writes is only a word, and the row naming it stands on its own, as it
    // always did: a gene under an operon that is not in the file is a gene.
    let keys: Vec<Key<'_>> = nodes.keys().copied().collect();
    for key in keys {
        let parents: Vec<Key<'_>> = nodes[&key]
            .parents
            .iter()
            .copied()
            .filter(|parent| nodes.contains_key(parent) || matches!(parent, Key::Named("gene", _)))
            .collect();
        if let Some(node) = nodes.get_mut(&key) {
            node.parents = parents.clone();
        }
        let first = nodes[&key].first;
        let (start, end, strand, name) = {
            let feature = nodes[&key].feature(key);
            (feature.start, feature.end, feature.strand, feature.name)
        };
        for parent in parents {
            let node = nodes.entry(parent).or_insert_with(|| Node::new(first));
            node.first = node.first.min(first);
            node.children.push(key);
            if node.row.is_none() {
                node.reaches(start, end, strand);
                let named = match (parent, rows.get(first)) {
                    (Key::Named("gene", id), Some((_, cols))) => Some(
                        gtf_attribute(cols.get(8).copied().unwrap_or_default(), "gene_name")
                            .unwrap_or(id)
                            .to_string(),
                    ),
                    _ => name.clone(),
                };
                node.names.extend(named);
            }
        }
    }

    let mut tops: Vec<(usize, Key<'_>)> = nodes
        .iter()
        .filter(|(_, node)| node.parents.is_empty())
        .map(|(key, node)| (node.first, *key))
        .collect();
    tops.sort_unstable();

    let mut out = Vec::new();
    for (_, key) in tops {
        let top = &nodes[&key];
        let feature = top.feature(key);
        // A cluster of genes, as an operon over the genes it holds, has no
        // exons of its own, and neither has anything with no parts under it.
        let children: Vec<(usize, Key<'_>)> = {
            let mut children: Vec<(usize, Key<'_>)> = top
                .children
                .iter()
                .map(|child| (nodes[child].first, *child))
                .collect();
            children.sort_unstable();
            children.dedup();
            children
        };
        if children
            .iter()
            .any(|(_, child)| gene_level(nodes[child].kind))
        {
            out.push((feature, Vec::new()));
            continue;
        }
        let mut transcripts: Vec<(Feature, Shape)> = Vec::new();
        if !top.parts.is_empty() {
            let shape = structure(top, (feature.start, feature.end));
            transcripts.push((feature.clone(), shape));
        }
        for (_, child) in &children {
            let node = &nodes[child];
            let own = node.feature(*child);
            let shape = structure(node, (own.start, own.end));
            transcripts.push((own, shape));
        }
        if transcripts.is_empty() {
            out.push((feature, Vec::new()));
            continue;
        }
        match level {
            Level::Gene => {
                let count = transcripts.len();
                let mut exons = Vec::new();
                let mut coding = Vec::new();
                for (own, mut shape) in transcripts {
                    // A transcript with nothing under it is all exon.
                    if shape.exons.is_empty() {
                        shape.exons.push((own.start, own.end));
                    }
                    exons.extend(shape.exons);
                    coding.extend(shape.coding);
                }
                let mut gene = feature.exons(exons).coding(coding);
                if count > 1 {
                    gene = gene.transcripts(count);
                }
                out.push((one_piece(gene), Vec::new()));
            }
            Level::Transcript | Level::Coding => {
                let gene = feature.name.clone();
                for (own, shape) in transcripts {
                    let mut transcript = own.exons(shape.exons).coding(shape.coding);
                    if let Some(gene) = gene
                        .clone()
                        .filter(|gene| Some(gene) != transcript.name.as_ref())
                    {
                        transcript = transcript.gene(gene);
                    }
                    out.push(if level == Level::Coding {
                        (transcript, shape.rows)
                    } else {
                        (one_piece(transcript), Vec::new())
                    });
                }
            }
        }
    }
    Ok(out)
}

/// How far the gene models over a window reach either side of it: the window
/// widened to the least start and the greatest end of every GFF3 or GTF row
/// that touches it, 0-based and half-open.
///
/// A gene model is put back together from all its rows, and the rows of one
/// over the window need not be over it themselves: a window inside an intron
/// touches the gene's row and its transcripts' and none of their exons. The
/// rows an index finds over a window are the ones over it, so a window read
/// through one is read again over this reach, which holds every row of every
/// model the window touches, since the gene's own row spans them all. A model
/// with no row of its own, exons naming a transcript the file never writes,
/// is spanned by nothing, and [`parentless`] says when the rows read hold
/// one. A row
/// that describes the sequence it is on, as NCBI's `region` does from its
/// first base to its last, is left out: it is no model, and it would widen
/// every window to its whole chromosome. A BED row is a model of its own, so
/// a BED reaches no further than the window.
///
/// Rows that do not read as a span are passed over, since the reader that
/// draws the window says what is wrong with them.
///
/// ```
/// use karyon::read::interval::reach;
/// use karyon::Region;
///
/// let gff = "##gff-version 3\n\
///            chr1\t.\tgene\t101\t900\t.\t+\t.\tID=g1\n\
///            chr1\t.\texon\t101\t200\t.\t+\t.\tParent=g1\n\
///            chr1\t.\texon\t801\t900\t.\t+\t.\tParent=g1\n";
/// // A window in the intron reaches both exons.
/// let window = Region::parse("chr1:401-500")?;
/// assert_eq!(reach(gff, &window, None), (100, 900));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn reach(text: &str, region: &Region, format: Option<Format>) -> (u64, u64) {
    let mut reach = (region.start(), region.end());
    if flavour(text, format) == Flavour::Bed {
        return reach;
    }
    for (_, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 5 || cols[0] != region.seq() || describes_sequence(&cols) {
            continue;
        }
        let (Ok(start), Ok(end)) = (cols[3].trim().parse::<u64>(), cols[4].trim().parse::<u64>())
        else {
            continue;
        };
        // 1-based and inclusive, so the start comes down by one.
        let start = start.saturating_sub(1);
        if end <= region.start() || start >= region.end() {
            continue;
        }
        reach = (reach.0.min(start), reach.1.max(end));
    }
    reach
}

/// Whether a GFF3 or GTF row over `region`, or anywhere in `text` for `None`,
/// belongs to a model with no row of its own in `text`: a part, an exon, a
/// coding or an untranslated stretch, naming a transcript no row on its
/// sequence is, or any row naming a GTF gene no row on its sequence is.
///
/// The readers put such a model together from all its parts on the sequence
/// and draw it from the first of them to the last, so how far it reaches is
/// in no row of its own, and the rows over a window, widened by [`reach`], do
/// not hold the parts outside it. A parent that has a row spans the rows
/// under it, so it lies over any window they lie over and is among the rows
/// an index finds there: only a parent with no row can be missing. A row
/// naming a parent the file does not write that is not a part stands on its
/// own and is drawn as it stands, so it is no such model.
///
/// ```
/// use karyon::read::interval::parentless;
/// use karyon::Region;
///
/// let window = Region::parse("chr1:101-200")?;
/// let exons = "##gff-version 3\n\
///              chr1\t.\texon\t101\t200\t.\t+\t.\tParent=t1\n\
///              chr1\t.\texon\t801\t900\t.\t+\t.\tParent=t1\n";
/// assert!(parentless(exons, Some(&window)));
/// let rooted = format!("{exons}chr1\t.\tmRNA\t101\t900\t.\t+\t.\tID=t1\n");
/// assert!(!parentless(&rooted, Some(&window)));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn parentless(text: &str, region: Option<&Region>) -> bool {
    let mut declared = std::collections::BTreeSet::new();
    let mut named = Vec::new();
    for (_, line) in lines(text) {
        let cols = columns(line);
        // A model is put together from the rows on its own sequence, so the
        // rows on another declare nothing it could be under.
        let elsewhere = region.is_some_and(|region| cols.first() != Some(&region.seq()));
        if cols.len() < 5 || elsewhere || describes_sequence(&cols) {
            continue;
        }
        let sequence = cols[0];
        if let Some((space, name)) = declared_as(&cols) {
            declared.insert((sequence, space, name));
        }
        if let Some(region) = region {
            let span = (cols[3].trim().parse::<u64>(), cols[4].trim().parse::<u64>());
            let (Ok(start), Ok(end)) = span else {
                continue;
            };
            // 1-based and inclusive, so the start comes down by one.
            if end <= region.start() || start.saturating_sub(1) >= region.end() {
                continue;
            }
        }
        let a_part = part(cols[2].trim()).is_some();
        named.extend(
            parents_of(&cols)
                .into_iter()
                .filter(|(space, _)| a_part || *space == "gene")
                .map(|(space, name)| (sequence, space, name)),
        );
    }
    named.iter().any(|key| !declared.contains(key))
}

/// Keeps a feature that touches the window.
///
/// A feature the window does not touch is not drawn, and it would still take a
/// row in the track's layout, which is what decides how tall the track is. So
/// it is dropped here rather than carried.
fn keep(model: (Feature, Vec<CdsRow>), region: &Region, into: &mut Vec<(Feature, Vec<CdsRow>)>) {
    if model.0.end <= region.start() || model.0.start >= region.end() {
        return;
    }
    into.push(model);
}

/// Whether a GFF3 row describes the sequence it is on rather than something
/// on it.
///
/// NCBI opens each sequence with a `region` row from its first base to its
/// last, named `ANONYMOUS` as often as not, Ensembl with a `chromosome` or a
/// `scaffold` row, and a converted GenBank file with a `source` or a
/// `databank_entry`. Drawn, it was a feature as long as the chromosome under
/// the window, named as if it were a gene. A row of one of these types that
/// does not start at base 1 marks part of a sequence, and is drawn.
fn describes_sequence(cols: &[&str]) -> bool {
    const SEQUENCES: &[&str] = &[
        "region",
        "chromosome",
        "scaffold",
        "supercontig",
        "databank_entry",
        "source",
    ];
    matches!(cols.get(3).map(|start| start.trim()), Some("1"))
        && cols
            .get(2)
            .is_some_and(|kind| SEQUENCES.contains(&kind.trim()))
}

/// The name another row could call this one by, and the kind of name it is.
///
/// GFF3 names every record it means to be a parent with `ID=`. GTF has no
/// such key: a gene is known by its `gene_id` and a transcript by its
/// `transcript_id`, and the two are kept apart because a file may give a gene
/// and its one transcript the same word.
fn declared_as<'a>(cols: &[&'a str]) -> Option<(&'static str, &'a str)> {
    let attributes = cols.get(8)?;
    if let Some(id) = raw_attribute(attributes, "ID") {
        return Some(("id", id));
    }
    match cols.get(2)?.trim() {
        "gene" => gtf_attribute(attributes, "gene_id").map(|id| ("gene", id)),
        "transcript" | "mRNA" => {
            gtf_attribute(attributes, "transcript_id").map(|id| ("transcript", id))
        }
        _ => None,
    }
}

/// The rows this one is a part of, by the names [`declared_as`] gives them.
///
/// GFF3 says so with `Parent=`, a list, and `Derives_from=`, which a
/// polypeptide uses to point at its transcript. A GTF row is part of its
/// transcript and of its gene, a transcript of its gene, and a gene of
/// nothing.
fn parents_of<'a>(cols: &[&'a str]) -> Vec<(&'static str, &'a str)> {
    let Some(attributes) = cols.get(8) else {
        return Vec::new();
    };
    let mut parents: Vec<(&'static str, &'a str)> = ["Parent", "Derives_from"]
        .into_iter()
        .filter_map(|key| raw_attribute(attributes, key))
        .flat_map(|list| list.split(','))
        .map(|parent| ("id", parent.trim()))
        .collect();
    if !parents.is_empty() {
        return parents;
    }
    let kind = cols.get(2).map_or("", |kind| kind.trim());
    if kind != "gene" {
        if kind != "transcript" && kind != "mRNA" {
            if let Some(transcript) = gtf_attribute(attributes, "transcript_id") {
                parents.push(("transcript", transcript));
            }
        }
        if let Some(gene) = gtf_attribute(attributes, "gene_id") {
            parents.push(("gene", gene));
        }
    }
    parents
}

/// Where an annotation names a gene, for placing a figure by the gene's name.
#[derive(Debug, Default)]
pub struct Named {
    /// Each row that goes by the name: its sequence, and its span, 0-based
    /// and half-open.
    pub spans: Vec<(String, u64, u64)>,
    /// Every name the annotation gives, each once, for suggesting the nearest
    /// when none is the name asked for.
    pub names: Vec<String>,
    /// The name as the annotation spells it, where a row goes by it: `rpoB`
    /// for a figure asked for as `rpob`.
    pub spelled: Option<String>,
}

/// The rows of an annotation, on any sequence, that go by `name`.
///
/// A GFF3 row goes by its `Name`, `gene`, `locus_tag` and `ID`, a GTF row by
/// its `gene_name`, `gene_id`, `transcript_name` and `transcript_id`, and a
/// BED row by its fourth column. The name is matched exactly, and failing
/// that in any case, so `RPOB` finds `rpoB` where nothing is called `RPOB`.
pub fn named(text: &str, name: &str) -> Named {
    let mut each = named_each(text, &[name]);
    Named {
        spans: each.spans.pop().unwrap_or_default(),
        names: each.names,
        spelled: each.spelled.pop().flatten(),
    }
}

/// Where an annotation names each of several genes, found in one pass.
#[derive(Debug, Default)]
pub(crate) struct NamedEach {
    /// For each name asked for, in the order asked, the rows that go by it,
    /// as [`Named::spans`] gives them for one.
    pub(crate) spans: Vec<Vec<(String, u64, u64)>>,
    /// Every name the annotation gives, each once.
    pub(crate) names: Vec<String>,
    /// For each name asked for, the name as the annotation spells it.
    pub(crate) spelled: Vec<Option<String>>,
}

/// The rows that go by each of `wanted`, as [`named`] finds them for one,
/// in one pass over the annotation.
///
/// One pass for all of them because the pass is the cost: a figure that
/// shaded eight genes of a 63 MB GFF3 read it eight times, a second each.
pub(crate) fn named_each(text: &str, wanted: &[&str]) -> NamedEach {
    let flavour = flavour(text, None);
    let mut names = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    let mut exact: Vec<Vec<(String, u64, u64)>> = vec![Vec::new(); wanted.len()];
    let mut loose: Vec<Vec<(String, u64, u64, String)>> = vec![Vec::new(); wanted.len()];
    for (_, line) in lines(text) {
        let cols = columns(line);
        let (start, end, calls): (Option<u64>, Option<u64>, Vec<String>) = match flavour {
            Flavour::Bed => (
                cols.get(1).and_then(|v| v.trim().parse().ok()),
                cols.get(2).and_then(|v| v.trim().parse().ok()),
                cols.get(3)
                    .map(|v| vec![v.trim().to_string()])
                    .unwrap_or_default(),
            ),
            Flavour::Gff3 => {
                let attributes = cols.get(8).copied().unwrap_or_default();
                let mut calls: Vec<String> = ["Name", "gene", "locus_tag", "ID"]
                    .into_iter()
                    .filter_map(|key| attribute(attributes, key))
                    .collect();
                calls.extend(
                    ["gene_name", "gene_id", "transcript_name", "transcript_id"]
                        .into_iter()
                        .filter_map(|key| gtf_attribute(attributes, key).map(str::to_string)),
                );
                (
                    cols.get(3)
                        .and_then(|v| v.trim().parse::<u64>().ok())
                        .map(|v| v.saturating_sub(1)),
                    cols.get(4).and_then(|v| v.trim().parse().ok()),
                    calls,
                )
            }
        };
        let (Some(start), Some(end), Some(sequence)) = (start, end, cols.first()) else {
            continue;
        };
        if describes_sequence(&cols) {
            continue;
        }
        for call in &calls {
            if seen.insert(call.clone()) {
                names.push(call.clone());
            }
        }
        let end = end.max(start + 1);
        for (index, name) in wanted.iter().enumerate() {
            if calls.iter().any(|call| call == name) {
                exact[index].push((sequence.to_string(), start, end));
            } else if let Some(call) = calls.iter().find(|call| call.eq_ignore_ascii_case(name)) {
                loose[index].push((sequence.to_string(), start, end, call.clone()));
            }
        }
    }
    let mut found = NamedEach {
        names,
        ..NamedEach::default()
    };
    for ((name, exact), loose) in wanted.iter().zip(exact).zip(loose) {
        if exact.is_empty() {
            found
                .spelled
                .push(loose.first().map(|(_, _, _, call)| call.clone()));
            found.spans.push(
                loose
                    .into_iter()
                    .map(|(sequence, start, end, _)| (sequence, start, end))
                    .collect(),
            );
        } else {
            found.spelled.push(Some((*name).to_string()));
            found.spans.push(exact);
        }
    }
    found
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
    bed_as_written(cols, at).map(one_piece)
}

/// One BED row with its blocks and its thick span as the row writes them, a
/// thick span from end to end included, which [`bed`] reads as one piece.
fn bed_as_written(cols: &[&str], at: usize) -> Result<Feature, ReadError> {
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
    Ok(feature
        .exons(blocks(cols, start, end))
        .coding(thick(cols, start, end)))
}

/// The coding stretch of a BED row of eight columns or more: its thickStart
/// and thickEnd, where they are a stretch of the row.
///
/// Read only where they make sense, because narrowPeak and the other formats
/// grown out of BED spend the same two columns on a signal and a p-value, and a
/// peak coloured as if it coded is a claim nobody made. A thickStart equal to
/// the thickEnd is how BED says nothing codes, so it is no stretch at all.
fn thick(cols: &[&str], start: u64, end: u64) -> Option<(u64, u64)> {
    let at = |index: usize| cols.get(index)?.trim().parse::<u64>().ok();
    let (from, to) = (at(6)?, at(7)?);
    (start <= from && from < to && to <= end).then_some((from, to))
}

/// The exons of a BED12 row: blockCount, blockSizes and blockStarts, the
/// starts counted from the row's own start.
///
/// All or nothing, and only where every block lies inside the row: a list
/// one short, or a block past the end, is a file this reader cannot vouch for,
/// and the row is drawn as the one interval its first three columns say.
fn blocks(cols: &[&str], start: u64, end: u64) -> Vec<(u64, u64)> {
    let list = |index: usize| -> Option<Vec<u64>> {
        cols.get(index)?
            .trim()
            .split(',')
            .filter(|value| !value.trim().is_empty())
            .map(|value| value.trim().parse::<u64>().ok())
            .collect()
    };
    let read = || -> Option<Vec<(u64, u64)>> {
        let count = cols.get(9)?.trim().parse::<usize>().ok()?;
        let (sizes, starts) = (list(10)?, list(11)?);
        if count == 0 || sizes.len() != count || starts.len() != count {
            return None;
        }
        starts
            .iter()
            .zip(&sizes)
            .map(|(&offset, &size)| {
                let from = start.checked_add(offset)?;
                let to = from.checked_add(size)?;
                (to <= end).then_some((from, to))
            })
            .collect()
    };
    read().unwrap_or_default()
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
    // The same refusal `bed` makes, and it was missing here. `Feature::new`
    // widens an inverted span into one base at the start, so 400 to 100 came
    // back as a one-base gene at 400: a real gene, drawn confidently, three
    // hundred bases from where either number put it.
    if end < start {
        return Err(ReadError::at(at, "end is before start"));
    }
    // 1-based inclusive to 0-based half-open: the start moves back one, the end
    // stays where it is because it was already one past the last base once the
    // count started at zero.
    let mut feature = Feature::new(start - 1, end);
    let kind = cols.get(2).map_or("", |kind| kind.trim());
    if let Some(name) = cols
        .get(8)
        .and_then(|attributes| gff3_name(attributes).or_else(|| gtf_name(kind, attributes)))
    {
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

/// The name a GTF record goes by.
///
/// GTF writes its ninth column as `key "value";` pairs, which no GFF3 key
/// reads, so every gene of a GTF file was drawn with no name at all. A gene is
/// called by its `gene_name` and failing that its `gene_id`; a transcript by
/// its own name first, since the isoforms of one gene share the gene's; and
/// anything else by the gene it belongs to.
fn gtf_name(kind: &str, attributes: &str) -> Option<String> {
    let keys: &[&str] = match kind {
        "gene" => &["gene_name", "gene_id"],
        "transcript" | "mRNA" => &["transcript_name", "transcript_id", "gene_name", "gene_id"],
        _ => &["gene_name", "gene_id", "transcript_name", "transcript_id"],
    };
    keys.iter()
        .find_map(|key| gtf_attribute(attributes, key))
        .and_then(|value| label(Some(value)))
}

/// One `key "value"` out of a GTF ninth column, without its quotes.
fn gtf_attribute<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    attributes.split(';').find_map(|pair| {
        let (found, value) = pair.trim().split_once(char::is_whitespace)?;
        (found == key).then(|| value.trim().trim_matches('"'))
    })
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

    /// Two genes of two transcripts and several exons each, as NCBI writes
    /// them, with the `region` row that opens the sequence.
    const MODELS: &str = "\
##gff-version 3
chr1\tRefSeq\tregion\t1\t50000\t.\t+\t.\tID=chr1:1..50000
chr1\t.\tgene\t1001\t9000\t.\t+\t.\tID=g1;Name=alpha
chr1\t.\tmRNA\t1001\t9000\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\texon\t1001\t1500\t.\t+\t.\tParent=t1
chr1\t.\texon\t4001\t4500\t.\t+\t.\tParent=t1
chr1\t.\texon\t8001\t9000\t.\t+\t.\tParent=t1
chr1\t.\tmRNA\t1201\t8500\t.\t+\t.\tID=t2;Parent=g1
chr1\t.\texon\t1201\t1500\t.\t+\t.\tParent=t2
chr1\t.\texon\t8001\t8500\t.\t+\t.\tParent=t2
chr1\t.\tgene\t20001\t30000\t.\t-\t.\tID=g2;Name=beta
chr1\t.\tmRNA\t20001\t30000\t.\t-\t.\tID=t3;Parent=g2
chr1\t.\texon\t20001\t21000\t.\t-\t.\tParent=t3
chr1\t.\texon\t29001\t30000\t.\t-\t.\tParent=t3
";

    /// The rows of `MODELS` that touch `[start, end)`, as an index finds
    /// them, behind the header.
    fn over(start: u64, end: u64) -> String {
        let mut text = String::from("##gff-version 3\n");
        for line in MODELS.lines().skip(1) {
            let cols: Vec<&str> = line.split('\t').collect();
            let (from, to) = (
                cols[3].parse::<u64>().unwrap() - 1,
                cols[4].parse::<u64>().unwrap(),
            );
            if from < end && to > start {
                text.push_str(line);
                text.push('\n');
            }
        }
        text
    }

    /// A window over any part of a gene, an intron included, reaches every
    /// row of it, and the genes read from the rows over that reach are the
    /// ones read from the whole file, with every exon. Read from the rows
    /// over the window alone, a window in the intron of alpha drew it with
    /// none of its exons.
    #[test]
    fn reach_covers_every_exon_of_a_gene_over_the_window() {
        for (start, end) in [
            (2_000u64, 2_100u64),
            (1_499, 1_501),
            (8_999, 9_000),
            (12_000, 13_000),
            (25_000, 25_100),
            (5_000, 21_000),
            (40_000, 45_000),
        ] {
            let region = Region::new("chr1", start, end).unwrap();
            let (from, to) = reach(&over(start, end), &region, None);
            assert!(from <= start && to >= end, "{region}");
            for level in [features, transcripts] {
                let whole = level(MODELS, &region, None).unwrap();
                assert_eq!(
                    level(&over(from, to), &region, None).unwrap(),
                    whole,
                    "{region}"
                );
            }
        }
        let intron = Region::new("chr1", 2_000, 2_100).unwrap();
        assert_eq!(reach(&over(2_000, 2_100), &intron, None), (1_000, 9_000));
        assert!(features(&over(2_000, 2_100), &intron, None).unwrap()[0]
            .exons
            .is_empty());
    }

    /// NCBI's `region` row spans the whole sequence and is no gene model,
    /// so it widens nothing; counted, every window would reach the whole
    /// chromosome and be read whole.
    #[test]
    fn reach_leaves_out_a_row_that_describes_the_sequence() {
        let gap = Region::new("chr1", 12_000, 13_000).unwrap();
        assert_eq!(reach(MODELS, &gap, None), (12_000, 13_000));
        // A BED row is its own model, and a BED reaches no further.
        let bed = "chr1\t0\t50000\tall\nchr1\t1000\t9000\talpha\n";
        assert_eq!(reach(bed, &gap, None), (12_000, 13_000));
        assert_eq!(reach(MODELS, &gap, Some(Format::Bed)), (12_000, 13_000));
    }

    /// Exons naming a transcript no row is make a model the reach of a window
    /// does not hold, and `parentless` finds them and only them: the
    /// transcript is drawn from its first exon to its last, so the rows over
    /// a window by one exon drew that exon alone. A gene naming an operon
    /// the file does not write stands on its own, and a row describing the
    /// sequence is no parent a part can belong to.
    #[test]
    fn parentless_finds_the_models_a_reach_does_not_hold() {
        let exons = "##gff-version 3\n\
            chr1\t.\texon\t1001\t2000\t.\t+\t.\tParent=t9\n\
            chr1\t.\texon\t50001\t51000\t.\t+\t.\tParent=t9\n";
        let window = Region::new("chr1", 1_500, 1_600).unwrap();
        let (from, to) = reach(exons, &window, None);
        let reached = Region::new("chr1", from, to).unwrap();
        let rows: String = exons
            .lines()
            .take(2)
            .map(|line| line.to_string() + "\n")
            .collect();
        assert_ne!(
            features(&rows, &window, None).unwrap(),
            features(exons, &window, None).unwrap()
        );
        assert!(parentless(&rows, Some(&reached)));
        assert!(parentless(exons, None));
        let between = Region::new("chr1", 10_000, 20_000).unwrap();
        assert!(!parentless(&rows, Some(&between)), "no row over it");
        for (start, end) in [(2_000, 2_100), (12_000, 13_000), (40_000, 45_000)] {
            let region = Region::new("chr1", start, end).unwrap();
            let (from, to) = reach(&over(start, end), &region, None);
            let reached = Region::new("chr1", from, to).unwrap();
            assert!(!parentless(&over(from, to), Some(&reached)), "{region}");
        }
        assert!(!parentless(MODELS, None));
        let operon = "chr1\t.\tgene\t1001\t2000\t.\t+\t.\tID=g1;Parent=op1\n";
        assert!(!parentless(operon, None));
        let described = "chr1\tRefSeq\tregion\t1\t90000\t.\t+\t.\tID=chr1\n\
                         chr1\t.\texon\t1001\t1500\t.\t+\t.\tParent=chr1\n";
        assert!(parentless(described, None));
    }

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

    /// The module doc has always said an inverted span stops the read, and
    /// only two of the three readers did it. GFF3 handed the pair to
    /// `Feature::new`, which widens it into a single base at the start, so the
    /// gene was drawn a whole interval from where either coordinate put it.
    #[test]
    fn an_inverted_span_stops_the_read_in_both_formats() {
        let gff = "##gff-version 3\nchr1\t.\tgene\t400\t100\t.\t+\t.\tID=x\n";
        let error = features(gff, &Region::parse("chr1:1-1000").unwrap(), None).unwrap_err();
        // Line two: the pragma is line one, and the numbering counts it.
        assert_eq!(error.line, 2);
        assert!(error.reason.contains("end is before start"), "{error}");

        let bed = "chr1\t400\t100\tx\n";
        let error = features(bed, &Region::parse("chr1:1-1000").unwrap(), None).unwrap_err();
        assert!(error.reason.contains("end is before start"), "{error}");

        // A single-base GFF3 gene, where start equals end, is still a gene.
        let one = "##gff-version 3\nchr1\t.\tgene\t100\t100\t.\t+\t.\tID=x\n";
        let genes = features(one, &Region::parse("chr1:1-1000").unwrap(), None).unwrap();
        assert_eq!(genes[0].len(), 1);
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

    fn names(features: &[Feature]) -> Vec<&str> {
        features
            .iter()
            .map(|feature| feature.name.as_deref().unwrap_or("?"))
            .collect()
    }

    /// NCBI writes a gene five times: the sequence it is on, the gene, its
    /// transcript, each exon and the CDS. Every one was drawn, the first as a
    /// feature as long as the chromosome named ANONYMOUS.
    #[test]
    fn a_gene_written_at_every_level_is_drawn_once() {
        let text = "\
##gff-version 3
chr1\tRefSeq\tregion\t1\t20000\t.\t+\t.\tID=chr1:1..20000;Name=ANONYMOUS;genome=chromosome
chr1\tRefSeq\tgene\t1500\t2900\t.\t+\t.\tID=gene-A;Name=geneA
chr1\tRefSeq\tmRNA\t1500\t2900\t.\t+\t.\tID=rna-A;Parent=gene-A;gene=geneA
chr1\tRefSeq\texon\t1500\t1900\t.\t+\t.\tID=exon-A-1;Parent=rna-A;gene=geneA
chr1\tRefSeq\texon\t2300\t2900\t.\t+\t.\tID=exon-A-2;Parent=rna-A;gene=geneA
chr1\tRefSeq\tCDS\t1550\t1900\t.\t+\t0\tID=cds-A;Parent=rna-A;Name=NP_1
chr1\tRefSeq\tCDS\t2300\t2800\t.\t+\t2\tID=cds-A;Parent=rna-A;Name=NP_1
chr1\tRefSeq\tpseudogene\t5000\t6000\t.\t-\t.\tID=gene-B;Name=geneB
";
        let features = read(text, "chr1:1-20000", None);
        assert_eq!(names(&features), ["geneA", "geneB"]);
        assert_eq!((features[0].start, features[0].end), (1_499, 2_900));
        // Drawn once, and with what the other levels said about it: two
        // exons, and the CDS inside them.
        assert_eq!(features[0].exons, [(1_499, 1_900), (2_299, 2_900)]);
        assert_eq!(features[0].coding, [(1_549, 1_900), (2_299, 2_800)]);
        assert_eq!(features[0].transcripts, 0);
        assert!(features[1].exons.is_empty() && features[1].coding.is_empty());
    }

    /// Ensembl writes a gene, its transcripts, and each transcript's exons,
    /// CDS and untranslated regions. A gene is drawn once with every exon any
    /// transcript uses, or each transcript is, naming its gene.
    const ENSEMBL: &str = "\
##gff-version 3
7\tensembl\tgene\t1001\t9000\t.\t+\t.\tID=gene:G1;Name=GENE1
7\tensembl\tmRNA\t1001\t9000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;Name=GENE1-201
7\tensembl\tfive_prime_UTR\t1001\t1200\t.\t+\t.\tParent=transcript:T1
7\tensembl\texon\t1001\t1500\t.\t+\t.\tParent=transcript:T1
7\tensembl\tCDS\t1201\t1500\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
7\tensembl\texon\t4001\t4500\t.\t+\t.\tParent=transcript:T1
7\tensembl\tCDS\t4001\t4500\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
7\tensembl\texon\t8001\t9000\t.\t+\t.\tParent=transcript:T1
7\tensembl\tCDS\t8001\t8600\t.\t+\t1\tID=CDS:P1;Parent=transcript:T1
7\tensembl\tmRNA\t1001\t9000\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;Name=GENE1-202
7\tensembl\texon\t1001\t1500\t.\t+\t.\tParent=transcript:T2
7\tensembl\texon\t6001\t6300\t.\t+\t.\tParent=transcript:T2
7\tensembl\texon\t8001\t9000\t.\t+\t.\tParent=transcript:T2
7\tensembl\tncRNA_gene\t12001\t15000\t.\t-\t.\tID=gene:G2;Name=LINC1
7\tensembl\tlnc_RNA\t12001\t15000\t.\t-\t.\tID=transcript:T3;Parent=gene:G2;Name=LINC1-201
7\tensembl\texon\t12001\t12500\t.\t-\t.\tParent=transcript:T3
7\tensembl\texon\t14001\t15000\t.\t-\t.\tParent=transcript:T3
";

    #[test]
    fn a_gene_is_drawn_once_with_every_exon_its_transcripts_use() {
        let genes = read(ENSEMBL, "7:1-20000", None);
        assert_eq!(names(&genes), ["GENE1", "LINC1"]);
        let gene = &genes[0];
        assert_eq!((gene.start, gene.end), (1_000, 9_000));
        assert_eq!(
            gene.exons,
            [
                (1_000, 1_500),
                (4_000, 4_500),
                (6_000, 6_300),
                (8_000, 9_000)
            ]
        );
        assert_eq!(
            gene.coding,
            [(1_200, 1_500), (4_000, 4_500), (8_000, 8_600)]
        );
        assert_eq!(gene.transcripts, 2);
        assert_eq!(gene.strand, Strand::Forward);
        // Nothing codes in a long non-coding RNA, so its exons are drawn at
        // full height rather than as one long untranslated region.
        assert_eq!(genes[1].exons, [(12_000, 12_500), (14_000, 15_000)]);
        assert!(genes[1].coding.is_empty());
        assert_eq!(genes[1].transcripts, 0);
    }

    #[test]
    fn each_transcript_is_a_feature_naming_its_gene_when_asked_for() {
        let each = transcripts(ENSEMBL, &region("7:1-20000"), None).unwrap();
        assert_eq!(names(&each), ["GENE1-201", "GENE1-202", "LINC1-201"]);
        assert_eq!(
            each[0].exons,
            [(1_000, 1_500), (4_000, 4_500), (8_000, 9_000)]
        );
        assert_eq!(
            each[0].coding,
            [(1_200, 1_500), (4_000, 4_500), (8_000, 8_600)]
        );
        assert_eq!(
            each[1].exons,
            [(1_000, 1_500), (6_000, 6_300), (8_000, 9_000)]
        );
        assert!(each[1].coding.is_empty());
        assert_eq!(each[0].gene.as_deref(), Some("GENE1"));
        assert_eq!(each[2].gene.as_deref(), Some("LINC1"));
        assert_eq!(each[0].transcripts, 0);
    }

    /// A bacterial gene is written over its CDS, and is the arrow it always
    /// was: this is what keeps every figure without introns as it was.
    #[test]
    fn a_gene_over_its_own_cds_is_one_piece() {
        let text = "\
##gff-version 3
NC_000962.3\tRefSeq\tgene\t759807\t763325\t.\t+\t.\tID=gene-Rv0667;Name=rpoB
NC_000962.3\tRefSeq\tCDS\t759807\t763325\t.\t+\t0\tID=cds-1;Parent=gene-Rv0667;Name=NP_1
NC_000962.3\tRefSeq\tgene\t763370\t767320\t.\t+\t.\tID=gene-Rv0668;Name=rpoC
NC_000962.3\tRefSeq\tCDS\t763370\t767000\t.\t+\t0\tID=cds-2;Parent=gene-Rv0668
";
        let genes = read(text, "NC_000962.3:759,000-768,000", None);
        assert_eq!(
            genes[0],
            Feature::new(759_806, 763_325)
                .name("rpoB")
                .strand(Strand::Forward)
        );
        // A gene longer than its CDS has an end that does not code, and not an
        // intron: one piece, drawn thinner past the CDS.
        assert!(genes[1].exons.is_empty());
        assert_eq!(genes[1].coding, [(763_369, 767_000)]);
    }

    #[test]
    fn exons_whose_transcript_is_not_in_the_file_are_still_one_transcript() {
        let text = "\
##gff-version 3
chr2\t.\texon\t100\t300\t.\t-\t.\tParent=tx9
chr2\t.\texon\t701\t900\t.\t-\t.\tParent=tx9
";
        let genes = read(text, "chr2:1-1000", None);
        assert_eq!(names(&genes), ["tx9"]);
        assert_eq!(genes[0].exons, [(99, 300), (700, 900)]);
        assert_eq!(genes[0].strand, Strand::Reverse);
    }

    #[test]
    fn a_gene_under_an_operon_leaves_the_operon_drawn_as_it_was() {
        let text = "\
##gff-version 3
chr1\t.\toperon\t100\t2000\t.\t+\t.\tID=op1;Name=opA
chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1;Parent=op1;Name=a
chr1\t.\tCDS\t100\t900\t.\t+\t0\tParent=g1
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=g2;Parent=op1;Name=b
";
        let genes = read(text, "chr1:1-3000", None);
        assert_eq!(names(&genes), ["opA"]);
        assert!(genes[0].exons.is_empty());
    }

    #[test]
    fn a_gtf_codes_through_its_stop_codon() {
        // GTF writes the stop codon outside the CDS, and it codes.
        let text = "\
chr1\tensembl\ttranscript\t101\t1000\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\"; gene_name \"abc\";
chr1\tensembl\texon\t101\t400\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";
chr1\tensembl\tCDS\t201\t400\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";
chr1\tensembl\texon\t601\t1000\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";
chr1\tensembl\tCDS\t601\t797\t.\t+\t2\tgene_id \"G\"; transcript_id \"T\";
chr1\tensembl\tstop_codon\t798\t800\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";
";
        let genes = read(text, "chr1:1-2000", None);
        assert_eq!(names(&genes), ["abc"]);
        assert_eq!(genes[0].exons, [(100, 400), (600, 1_000)]);
        assert_eq!(genes[0].coding, [(200, 400), (600, 800)]);
    }

    /// UCSC writes a transcript per BED12 row: exons as blocks counted from
    /// the row's start, and the coding stretch as its thick span.
    #[test]
    fn a_bed12_row_is_read_with_its_blocks_and_its_thick_span() {
        let bed12 = "chr1\t1000\t9000\tNM_1\t0\t-\t1200\t8600\t0\t3\t500,500,1000,\t0,3000,7000,\n";
        let feature = &read(bed12, "chr1:1-10000", None)[0];
        assert_eq!(
            feature.exons,
            [(1_000, 1_500), (4_000, 4_500), (8_000, 9_000)]
        );
        assert_eq!(feature.coding, [(1_200, 8_600)]);
        assert_eq!(feature.strand, Strand::Reverse);

        // UCSC's default thick span is the whole row, which is one piece.
        let bed9 = "chr7\t100\t900\tPos1\t0\t+\t100\t900\t255,0,0\n";
        assert_eq!(
            read(bed9, "chr7:1-1000", None)[0],
            Feature::new(100, 900).name("Pos1").strand(Strand::Forward)
        );
        // Equal ends are how BED says nothing codes.
        let noncoding = "chr7\t100\t900\tx\t0\t+\t900\t900\n";
        assert!(read(noncoding, "chr7:1-1000", None)[0].coding.is_empty());
        // A narrowPeak spends the same columns on a signal and a p-value,
        // whole numbers as often as not, and they are no stretch of the row.
        for peak in [
            "chr7\t100\t900\tpeak1\t0\t.\t5.38\t12\t-1\t350\n",
            "chr7\t100\t900\tpeak2\t0\t.\t12\t30\t-1\t350\n",
        ] {
            let feature = &read(peak, "chr7:1-1000", None)[0];
            assert!(
                feature.coding.is_empty() && feature.exons.is_empty(),
                "{peak}"
            );
        }
        // Block lists one short, or a block past the row, are not read.
        let short = "chr1\t1000\t9000\tx\t0\t+\t1000\t9000\t0\t3\t500,500,\t0,3000,7000,\n";
        assert!(read(short, "chr1:1-10000", None)[0].exons.is_empty());
        let past = "chr1\t1000\t9000\tx\t0\t+\t1000\t9000\t0\t2\t500,5000,\t0,7000,\n";
        assert!(read(past, "chr1:1-10000", None)[0].exons.is_empty());
    }

    #[test]
    fn a_row_that_marks_part_of_a_sequence_is_still_drawn() {
        // Only a row from the first base describes the sequence: Ensembl's
        // chromosome line is left out, and a region further along is a region.
        let text = "\
##gff-version 3
2\tGRCh38\tchromosome\t1\t242193529\t.\t.\t.\tID=chromosome:2
2\tmanual\tregion\t100\t200\t.\t+\t.\tID=roi;Name=target
";
        assert_eq!(names(&read(text, "2:1-1000", None)), ["target"]);
    }

    #[test]
    fn a_part_whose_whole_is_not_in_the_file_is_drawn() {
        // Cut down to CDS rows, the genes they name are elsewhere, and a
        // file of nothing but parts is still a file of features.
        let text = "\
##gff-version 3
chr1\t.\tCDS\t100\t400\t.\t+\t0\tID=cds1;Parent=gene1;Name=dnaA
chr1\t.\tCDS\t500\t900\t.\t+\t0\tID=cds2;Parent=gene2;Name=dnaN
";
        assert_eq!(names(&read(text, "chr1:1-1000", None)), ["dnaA", "dnaN"]);
    }

    #[test]
    fn a_polypeptide_is_part_of_the_transcript_it_derives_from() {
        let text = "\
##gff-version 3
chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1;Name=abc
chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tpolypeptide\t150\t850\t.\t+\t.\tID=p1;Derives_from=t1;Name=ABC
";
        assert_eq!(names(&read(text, "chr1:1-1000", None)), ["abc"]);
    }

    /// GTF writes its attributes as `key "value";`, which no GFF3 key reads,
    /// so every gene came out nameless, and the gene beside its own CDS.
    #[test]
    fn a_gtf_is_drawn_once_a_gene_and_by_its_name() {
        let gencode = "\
chr1\tHAVANA\tgene\t1500\t2900\t.\t+\t.\tgene_id \"ENSG1\"; gene_name \"geneA\";
chr1\tHAVANA\ttranscript\t1500\t2900\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\"; gene_name \"geneA\";
chr1\tHAVANA\texon\t1500\t1900\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\"; exon_number 1;
chr1\tHAVANA\tCDS\t1550\t1900\t.\t+\t0\tgene_id \"ENSG1\"; transcript_id \"ENST1\";
chr1\tHAVANA\tgene\t3200\t4400\t.\t-\t.\tgene_id \"ENSG2\";
";
        assert_eq!(
            names(&read(gencode, "chr1:1-5000", None)),
            ["geneA", "ENSG2"]
        );

        // StringTie writes no gene rows: each transcript stands for its exons,
        // and the isoforms of one gene are told apart by their own names.
        let stringtie = "\
chr1\tStringTie\ttranscript\t100\t900\t.\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\";
chr1\tStringTie\texon\t100\t300\t.\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.1\";
chr1\tStringTie\ttranscript\t100\t700\t.\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.2\";
chr1\tStringTie\texon\t500\t700\t.\t+\t.\tgene_id \"STRG.1\"; transcript_id \"STRG.1.2\";
";
        // Its two transcripts are one gene, which says it merged them, and
        // each is a feature of its own when they are asked for one at a time.
        let genes = read(stringtie, "chr1:1-1000", None);
        assert_eq!(names(&genes), ["STRG.1"]);
        assert_eq!(genes[0].transcripts, 2);
        assert_eq!(genes[0].exons, [(99, 300), (499, 700)]);
        let each = transcripts(stringtie, &region("chr1:1-1000"), None).unwrap();
        assert_eq!(names(&each), ["STRG.1.1", "STRG.1.2"]);
        assert_eq!(each[0].gene.as_deref(), Some("STRG.1"));

        // A table browser GTF has exons and nothing above them, and every row
        // names its transcript and its gene: one gene of two exons, drawn
        // once, where each exon was drawn as a gene of its own.
        let exons = "\
chr1\thg38\texon\t100\t300\t.\t+\t.\tgene_id \"NM_1\"; transcript_id \"NM_1\"; gene_name \"abc\";
chr1\thg38\texon\t500\t700\t.\t+\t.\tgene_id \"NM_1\"; transcript_id \"NM_1\"; gene_name \"abc\";
";
        let genes = read(exons, "chr1:1-1000", None);
        assert_eq!(names(&genes), ["abc"]);
        assert_eq!(genes[0].exons, [(99, 300), (499, 700)]);
        assert_eq!((genes[0].start, genes[0].end), (99, 700));
    }

    #[test]
    fn a_broken_row_stops_the_file_even_where_it_would_have_been_left_out() {
        let text = "\
##gff-version 3
chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1;Name=abc
chr1\t.\tCDS\tten\t900\t.\t+\t0\tID=c1;Parent=g1
";
        let error = features(text, &region("chr1:1-1000"), None).unwrap_err();
        assert_eq!(error.line, 3);
    }

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

    /// Several names looked up in one pass find what each looked up alone
    /// finds, the loose match in any case included, in the order asked.
    #[test]
    fn several_names_in_one_pass_are_each_what_one_name_finds() {
        let each = named_each(GFF3, &["rpoC", "RPOB", "katG"]);
        assert_eq!(each.spans.len(), 3);
        for (index, name) in ["rpoC", "RPOB", "katG"].into_iter().enumerate() {
            let one = named(GFF3, name);
            assert_eq!(each.spans[index], one.spans, "{name}");
            assert_eq!(each.spelled[index], one.spelled, "{name}");
            assert_eq!(each.names, one.names, "{name}");
        }
        assert_eq!(
            each.spans[1],
            [("NC_000962.3".to_string(), 759_806, 763_325)]
        );
        assert_eq!(each.spelled[1].as_deref(), Some("rpoB"));
        assert!(each.spans[2].is_empty());
        assert_eq!(named_each(BED, &[]).spans.len(), 0);
    }

    /// What a ruler of codons reads: a gene written over its CDS codes from
    /// end to end, which the drawing reader folds into one piece and so into
    /// the same answer as a gene that codes nowhere.
    #[test]
    fn a_gene_over_its_cds_codes_from_end_to_end_where_one_that_has_none_does_not() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tgene\t101\t400\t.\t+\t.\tID=g1;Name=abcD\n\
                   chr1\t.\tCDS\t101\t400\t.\t+\t0\tParent=g1\n\
                   chr1\t.\tgene\t501\t800\t.\t-\t.\tID=g2;Name=efgH\n";
        let window = Region::new("chr1", 0, 1_000).unwrap();
        let drawn = transcripts(gff, &window, None).unwrap();
        assert!(drawn.iter().all(|feature| feature.coding.is_empty()));
        let read = coding(gff, &window, None).unwrap();
        assert_eq!(read.len(), 2);
        assert_eq!(read[0].feature.name.as_deref(), Some("abcD"));
        assert_eq!(read[0].feature.coding, [(100, 400)]);
        assert_eq!(read[0].feature.strand, Strand::Forward);
        assert_eq!(
            read[0].rows,
            [CdsRow {
                start: 100,
                end: 400,
                phase: Some(0),
                partial: false
            }]
        );
        assert!(
            read[1].feature.coding.is_empty() && read[1].rows.is_empty(),
            "a gene row alone codes nowhere"
        );
    }

    /// A CDS row under no gene codes over itself, and the gene beside it,
    /// under nothing either, codes nowhere: two features, one coding.
    #[test]
    fn a_cds_under_nothing_codes_over_itself() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tgene\t101\t400\t.\t+\t.\tID=g1_gene;Name=abcD\n\
                   chr1\t.\tCDS\t101\t400\t.\t+\t0\tID=g1;Name=abcD\n";
        let read = coding(gff, &Region::new("chr1", 0, 1_000).unwrap(), None).unwrap();
        let coded: Vec<&Cds> = read
            .iter()
            .filter(|cds| !cds.feature.coding.is_empty())
            .collect();
        assert_eq!(coded.len(), 1, "{read:?}");
        assert_eq!(coded[0].feature.coding, [(100, 400)]);
        assert_eq!(coded[0].rows.len(), 1);
    }

    /// GFF3 writes one feature in pieces as rows under one `ID=`, and a CDS
    /// under nothing written so is every piece of it, where the last row
    /// alone used to stand for the whole.
    #[test]
    fn a_cds_under_nothing_in_rows_under_one_id_is_every_row() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tCDS\t1001\t1100\t.\t+\t0\tID=c1;Name=abcA\n\
                   chr1\t.\tCDS\t1201\t1300\t.\t+\t2\tID=c1;Name=abcA\n";
        let read = coding(gff, &Region::new("chr1", 0, 2_000).unwrap(), None).unwrap();
        assert_eq!(read.len(), 1, "{read:?}");
        let cds = &read[0];
        assert_eq!((cds.feature.start, cds.feature.end), (1_000, 1_300));
        assert_eq!(cds.feature.coding, [(1_000, 1_100), (1_200, 1_300)]);
        assert_eq!(cds.feature.exons, [(1_000, 1_100), (1_200, 1_300)]);
        let phases: Vec<Option<u8>> = cds.rows.iter().map(|row| row.phase).collect();
        assert_eq!(phases, [Some(0), Some(2)]);
    }

    /// Two CDS rows sharing a base, as NCBI writes a ribosomal slippage, are
    /// one stretch to [`Feature::coding`] and two rows here, each with its
    /// phase; a GTF's start codon, inside the CDS, is no row of its own, and
    /// neither is a stop codon inside a CDS row, while one past it is.
    #[test]
    fn the_rows_of_a_cds_come_back_apart_with_their_phases() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tgene\t1001\t1600\t.\t+\t.\tID=g1;Name=gag-pol\n\
                   chr1\t.\tCDS\t1001\t1150\t.\t+\t0\tID=cds-1;Parent=g1\n\
                   chr1\t.\tCDS\t1150\t1600\t.\t+\t0\tID=cds-1;Parent=g1\n";
        let read = coding(gff, &Region::new("chr1", 0, 2_000).unwrap(), None).unwrap();
        assert_eq!(read[0].feature.coding, [(1_000, 1_600)]);
        let spans: Vec<(u64, u64)> = read[0]
            .rows
            .iter()
            .map(|row| (row.start, row.end))
            .collect();
        assert_eq!(spans, [(1_000, 1_150), (1_149, 1_600)]);
        let gtf = "chr1\t.\tCDS\t101\t397\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";\n\
                   chr1\t.\tstart_codon\t101\t103\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";\n\
                   chr1\t.\tstop_codon\t398\t400\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";\n\
                   chr1\t.\tCDS\t501\t900\t.\t-\t.\tgene_id \"h\"; transcript_id \"u\";\n\
                   chr1\t.\tstop_codon\t501\t503\t.\t-\t0\tgene_id \"h\"; transcript_id \"u\";\n";
        let read = coding(gtf, &Region::new("chr1", 0, 1_000).unwrap(), None).unwrap();
        let spans: Vec<(u64, u64)> = read[0]
            .rows
            .iter()
            .map(|row| (row.start, row.end))
            .collect();
        assert_eq!(spans, [(100, 397), (397, 400)]);
        assert_eq!(read[1].rows.len(), 1, "{read:?}");
        assert_eq!(read[1].rows[0].phase, None, "a phase of . is none");
    }

    /// Only a CDS under nothing with an ID of its own can have pieces a
    /// window through an index leaves out: one under a gene is spanned by
    /// the gene's row, one with no ID is one row, and a GTF row names its
    /// transcript.
    #[test]
    fn a_cds_under_nothing_with_an_id_is_one_whose_pieces_may_be_elsewhere() {
        let row = |attributes: &str| format!("chr1\t.\tCDS\t101\t400\t.\t+\t0\t{attributes}\n");
        assert!(pieced(&row("ID=c1;Name=abcA"), "chr1"));
        assert!(!pieced(&row("ID=c1;Name=abcA"), "chr2"));
        assert!(!pieced(&row("ID=c1;Parent=g1"), "chr1"));
        assert!(!pieced(&row("Name=abcA"), "chr1"));
        assert!(!pieced(
            &row("gene_id \"g1\"; transcript_id \"t1\";"),
            "chr1"
        ));
        let gene = "chr1\t.\tgene\t101\t400\t.\t+\t.\tID=g1\n";
        assert!(!pieced(gene, "chr1"));
    }

    /// A CDS cut short at its 5' end, as NCBI writes one at a contig's edge,
    /// says so at the start on the forward strand and at the end on the
    /// reverse; a 3' end cut short is no 5' end.
    #[test]
    fn a_cds_row_says_when_it_goes_on_past_its_5_prime_end() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tCDS\t1000\t1300\t.\t+\t1\tID=a;partial=true;start_range=.,1000\n\
                   chr1\t.\tCDS\t2001\t2301\t.\t-\t0\tID=b;partial=true;end_range=2301,.\n\
                   chr1\t.\tCDS\t3001\t3300\t.\t-\t0\tID=c;partial=true;start_range=.,3001\n";
        let read = coding(gff, &Region::new("chr1", 0, 4_000).unwrap(), None).unwrap();
        let said: Vec<(Option<u8>, bool)> = read
            .iter()
            .map(|cds| (cds.rows[0].phase, cds.rows[0].partial))
            .collect();
        assert_eq!(said, [(Some(1), true), (Some(0), true), (Some(0), false)]);
    }

    /// Each transcript comes back with its own pieces, and a CDS split by an
    /// intron is two of them, not the span between its ends.
    #[test]
    fn a_spliced_cds_comes_back_in_its_pieces_per_transcript() {
        let gff = "##gff-version 3\n\
                   chr1\t.\tgene\t1001\t9000\t.\t+\t.\tID=g1;Name=alpha\n\
                   chr1\t.\tmRNA\t1001\t9000\t.\t+\t.\tID=t1;Parent=g1\n\
                   chr1\t.\texon\t1001\t1500\t.\t+\t.\tParent=t1\n\
                   chr1\t.\texon\t8001\t9000\t.\t+\t.\tParent=t1\n\
                   chr1\t.\tCDS\t1201\t1500\t.\t+\t0\tParent=t1\n\
                   chr1\t.\tCDS\t8001\t8600\t.\t+\t0\tParent=t1\n";
        let read = coding(gff, &Region::new("chr1", 0, 10_000).unwrap(), None).unwrap();
        assert_eq!(read.len(), 1);
        assert_eq!(read[0].feature.coding, [(1_200, 1_500), (8_000, 8_600)]);
        assert_eq!(read[0].feature.exons, [(1_000, 1_500), (8_000, 9_000)]);
        assert_eq!(read[0].feature.gene.as_deref(), Some("alpha"));
        assert_eq!(read[0].rows.len(), 2);
    }

    /// A BED row's thick span is what codes, even from end to end, and a BED
    /// of six columns codes nowhere.
    #[test]
    fn a_bed_row_codes_over_its_thick_span_even_when_it_is_the_whole_row() {
        let bed = "chr1\t100\t400\tabcD\t0\t-\t100\t400\n\
                   chr1\t500\t800\tefgH\t0\t+\n";
        let window = Region::new("chr1", 0, 1_000).unwrap();
        let read = coding(bed, &window, None).unwrap();
        assert_eq!(read[0].feature.coding, [(100, 400)]);
        assert_eq!(read[0].feature.strand, Strand::Reverse);
        assert_eq!(
            read[0].rows,
            [CdsRow {
                start: 100,
                end: 400,
                phase: None,
                partial: false
            }]
        );
        assert!(read[1].feature.coding.is_empty() && read[1].rows.is_empty());
        // Drawn, the same row is one piece, as it was.
        assert!(features(bed, &window, None).unwrap()[0].coding.is_empty());
    }

    /// The table a CDS names is the one over the span asked about, on the
    /// sequence asked about, and a row that names none is passed over.
    #[test]
    fn the_translation_table_is_the_one_the_cds_over_the_span_names() {
        let gff = "##gff-version 3\n\
                   chrM\t.\tgene\t3307\t4262\t.\t+\t.\tID=g1;Name=ND1\n\
                   chrM\t.\tCDS\t3307\t4262\t.\t+\t0\tParent=g1;transl_table=2\n\
                   chr1\t.\tCDS\t3307\t4262\t.\t+\t0\tID=c2;transl_table=1\n\
                   chrM\t.\tCDS\t5000\t5100\t.\t+\t0\tID=c3\n";
        assert_eq!(translation_table(gff, "chrM", 3_306, 3_309), Some(2));
        assert_eq!(translation_table(gff, "chr1", 3_306, 3_309), Some(1));
        assert_eq!(translation_table(gff, "chrM", 4_999, 5_100), None);
        assert_eq!(
            translation_table(gff, "chrM", 4_262, 4_999),
            None,
            "half-open"
        );
        assert_eq!(
            translation_table("chrM\t3306\t4262\n", "chrM", 0, 9_000),
            None
        );
    }
}
