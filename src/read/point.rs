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
//! alternate allele go past without a word. The last of those are a gVCF's
//! reference blocks, a statement that nothing happened rather than a call,
//! which bcftools writes with an `ALT` of `<*>` and GATK with `<NON_REF>`. A
//! variant row of a gVCF names the placeholder too, after its own allele, as
//! `T,<NON_REF>`, and is a call of `T` and nothing else: the placeholder is
//! dropped allele by allele, so the indices a genotype names keep pointing at
//! the alleles they were written against. A row that does not parse stops the
//! read on its line, because a VCF that cannot be read is not a VCF and a
//! figure short of the calls it should have shows nothing wrong on its face.
//!
//! # A cohort's genotypes
//!
//! After the eight columns of a site, a VCF of several samples has a column
//! of keys, `FORMAT`, and one column per sample, named on the `#CHROM` line in
//! the same order. [`genotypes`] reads the `GT` of each, found by its name
//! among the keys rather than taken to be first, and [`samples`] reads the
//! names alone, from the header, without reading a row.
//!
//! An association table of two columns names no sequence, so every row in one
//! is taken to be about the sequence on display; three columns puts the name in
//! front and it is matched. A word where a position belongs is a header, but
//! only on the first line worth looking at, so the same word further down the
//! file is an error rather than a row that disappears.

use std::cmp::Ordering;

use crate::track::genotype::{placeholder, Unread};
use crate::{Association, Genotype, GenotypeSite, Region, Variant};

use super::{columns, lines, number, ReadError};

/// Reads calls from VCF text.
///
/// VCF `POS` is 1-based, so the variant lands at `POS - 1`. The category is
/// taken from the `ANN` or `BCSQ` consequence when one is there, and otherwise
/// from the shape of the call: a substitution, an insertion or a deletion,
/// which is `REF` against `ALT` and needs no annotation.
///
/// The value is the allele fraction, from `AF` in `INFO` when present, and
/// none when not: a call with no fraction is still a call, and stands full
/// height, but it is not a fraction of 1, which is what the axis would say.
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
        // A row with no alternate allele, or with only the placeholder for
        // one, is a reference block, which is most of what a gVCF holds, and
        // it is not a call.
        if reference_block(fields[4]) {
            continue;
        }

        let (reference, info) = (fields[3], fields[7]);
        // A whole genome VCF is a normal thing to hand over when only a window
        // is being drawn, so the rest of it stops here rather than being
        // carried into a track that would not draw it. What is kept is what
        // REF spells, not the anchor base alone: a deletion is written one base
        // to the left of the bases it removes, so a call anchored just outside
        // the window can still be a call about the window. The end saturates,
        // since a POS at the top of the number line leaves no room after it for
        // what REF spells, and the add panicked there instead of reading on.
        let spelled = pos.saturating_add(reference.len().max(1) as u64);
        if pos >= region.end() || spelled <= region.start() {
            continue;
        }
        let alternates: Vec<&str> = fields[4].split(',').collect();
        let fractions = allele_fractions(info, alternates.len(), line)?;
        for (index, alt) in alternates.iter().enumerate() {
            // A gVCF's variant rows name the placeholder after their own
            // alleles, and it is no call: drawn, `T,<NON_REF>` was a
            // substitution and a second lollipop called `non_ref` at one base.
            // Skipped here and not taken out of the list, so a fraction given
            // per allele stays with its allele.
            if placeholder(alt) {
                continue;
            }
            let category = consequence(info, alt).unwrap_or_else(|| shape(reference, alt));
            let mut call = Variant::new(pos).category(category);
            if let Some(fraction) = fractions.as_ref().and_then(|all| all[index]) {
                call = call.value(fraction);
            }
            calls.push(call);
        }
    }
    Ok(calls)
}

/// Whether an `ALT` column names no allele at all: `.`, or only the
/// placeholder a gVCF's reference blocks carry.
fn reference_block(alt: &str) -> bool {
    alt == "." || alt.split(',').all(placeholder)
}

/// The fields of a row, split as [`columns`] splits them, one at a time, so a
/// row of a thousand samples outside the window is not split past its `POS`.
fn fields_of(row: &str) -> impl Iterator<Item = &str> {
    let tabbed = row.contains('\t');
    let tabs = tabbed.then(|| row.split('\t'));
    let spaces = (!tabbed).then(|| row.split_whitespace());
    tabs.into_iter()
        .flatten()
        .chain(spaces.into_iter().flatten())
}

/// The `#CHROM` line of a VCF, numbered from one, split into its fields.
///
/// The last one in the header, which ends at the first line that is not a
/// comment: the header comes first, and a row is not looked at.
fn chrom_line(text: &str) -> Option<(usize, Vec<&str>)> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut found = None;
    for (index, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            found = Some((index + 1, columns(line)));
        }
    }
    found
}

/// The samples a VCF names, in the order of its columns: the names on the
/// `#CHROM` line after `FORMAT`.
///
/// Read off the header alone, so a caller can count the samples, or check a
/// name it was given, before any row is read. Empty for a VCF of sites only,
/// and for text with no `#CHROM` line.
///
/// ```
/// use karyon::read::point::samples;
///
/// let vcf = "##fileformat=VCFv4.2\n\
///            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
///            chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\t1\n";
/// assert_eq!(samples(vcf), ["S1", "S2"]);
/// ```
pub fn samples(text: &str) -> Vec<String> {
    chrom_line(text)
        .and_then(|(_, fields)| fields.get(9..).map(<[&str]>::to_vec))
        .unwrap_or_default()
        .into_iter()
        .map(str::to_string)
        .collect()
}

/// What a VCF's genotype columns hold over a window.
#[derive(Debug, Clone, PartialEq)]
pub struct Genotypes {
    /// The samples, in the order their rows are drawn: the order of the
    /// `#CHROM` line, or the order they were asked for in.
    pub samples: Vec<String>,
    /// The records in the window, in the order of the file, each with a call
    /// per sample in the order of `samples`.
    pub sites: Vec<GenotypeSite>,
    /// How many rows the file holds, on any sequence.
    pub records: usize,
}

/// Reads the genotype of every sample at every record of a VCF in a window.
///
/// Column 9 is `FORMAT`, the keys of the sample columns joined by `:`, and
/// `GT` is found among them by name rather than taken to be first. Columns 10
/// on are one sample each, named on the `#CHROM` line in the same order,
/// which has to be there: `bcftools view -H` leaves it out. A sample's column
/// is its values in the order of the keys, and a value left off its end, as a
/// caller writes `.` for a sample it has nothing for, is not there.
///
/// `GT` is allele indices joined by `/`, or by `|` where the call is phased,
/// with the leading `/` or `|` VCF 4.4 allows: 0 is `REF` and `i` the `i`th
/// allele of `ALT`, and `.` is a copy nobody called. One index is a haploid
/// call, two a diploid one, and more a polyploid one. A call with any copy
/// unknown, as `./1` is, is no call: a whole call is what the track draws,
/// and half of one read as a whole would draw a heterozygote as homozygous.
/// A row whose `FORMAT` has no `GT` leaves every sample of it without a call.
///
/// A record is kept by the rule [`variants`] keeps one by, when what `REF`
/// spells touches the window, and `POS` lands at `POS - 1`. A row whose `ALT`
/// is `.`, or only a placeholder, `<NON_REF>` or `<*>`, is a gVCF's reference
/// block and is skipped. `wanted` keeps only the samples it names, in the
/// order it names them.
///
/// # Errors
///
/// A VCF with no `#CHROM` line or no sample on it, and one naming a sample
/// twice; a sample `wanted` names that the file does not, or names twice. The
/// line that does not read: fewer than 8 columns, a POS that is not one, a row
/// with more or fewer samples than the `#CHROM` line names, a `GT` that is not
/// indices and dots, and an index past the alleles of its row. And a window
/// holding records of which none carries `GT`, since nothing is known there
/// of any sample.
///
/// ```
/// use karyon::read::point::genotypes;
/// use karyon::{GenotypeState, Region};
///
/// let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
///            chr1\t100\trs1\tC\tT,G\t.\t.\t.\tGT:DP\t0/1:20\t1/2:18\t./.:0\n";
/// let read = genotypes(vcf, &Region::parse("chr1:1-200")?, None)?;
/// let site = &read.sites[0];
/// assert_eq!(site.position, 99);
/// assert_eq!(site.call(0).state(), GenotypeState::Heterozygous);
/// assert_eq!(site.call(1).state(), GenotypeState::Alternate);
/// assert_eq!(site.call(2).state(), GenotypeState::NotCalled);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn genotypes(
    text: &str,
    region: &Region,
    wanted: Option<&[String]>,
) -> Result<Genotypes, ReadError> {
    let Some((header, names)) = chrom_line(text) else {
        return Err(ReadError::whole(
            "the VCF has no #CHROM line, which is where it names its samples; keep the \
             header, which bcftools view -H leaves out",
        ));
    };
    let held: &[&str] = names.get(9..).unwrap_or_default();
    if held.is_empty() {
        return Err(ReadError::at(
            header,
            "the VCF names no samples on its #CHROM line, so it holds no genotypes; \
             --variants draws its calls",
        ));
    }
    let mut seen = std::collections::HashSet::with_capacity(held.len());
    for name in held {
        if !seen.insert(*name) {
            return Err(ReadError::at(
                header,
                format!("the #CHROM line names {name} twice"),
            ));
        }
    }
    let columns: Vec<usize> = match wanted {
        None => (0..held.len()).collect(),
        Some(wanted) => {
            let mut picked: Vec<usize> = Vec::with_capacity(wanted.len());
            for name in wanted {
                let Some(at) = held.iter().position(|held| held == name) else {
                    return Err(ReadError::whole(format!(
                        "the #CHROM line names no sample called {name}"
                    )));
                };
                if picked.contains(&at) {
                    return Err(ReadError::whole(format!("{name} is asked for twice")));
                }
                picked.push(at);
            }
            picked
        }
    };
    let samples: Vec<String> = columns.iter().map(|at| held[*at].to_string()).collect();

    let mut sites = Vec::new();
    let mut records = 0;
    let mut typed = false;
    for (line, row) in lines(text) {
        records += 1;
        let mut fields = fields_of(row);
        let (Some(chrom), Some(pos), Some(id), Some(reference), Some(alt)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            return Err(ReadError::at(
                line,
                format!(
                    "a VCF line has at least 8 columns, this one has {}",
                    fields_of(row).count()
                ),
            ));
        };
        if chrom != region.seq() {
            continue;
        }
        let pos = position(pos, "POS", line)?;
        if reference_block(alt) {
            continue;
        }
        // Kept by the rule `variants` keeps a row by, saturating as it does.
        let spelled = pos.saturating_add(reference.len().max(1) as u64);
        if pos >= region.end() || spelled <= region.start() {
            continue;
        }

        // QUAL, FILTER, INFO, FORMAT and the samples, split only now that the
        // record is one the figure draws.
        let rest: Vec<&str> = fields.collect();
        if rest.len() + 5 < 8 {
            return Err(ReadError::at(
                line,
                format!(
                    "a VCF line has at least 8 columns, this one has {}",
                    rest.len() + 5
                ),
            ));
        }
        let given = rest.len().saturating_sub(4);
        if given != held.len() {
            return Err(ReadError::at(
                line,
                format!(
                    "this line has {given} sample{} and the #CHROM line names {}",
                    if given == 1 { "" } else { "s" },
                    held.len()
                ),
            ));
        }
        let gt = rest[3].split(':').position(|key| key == "GT");
        typed |= gt.is_some();
        let alternates: Vec<String> = alt.split(',').map(str::to_string).collect();

        let mut calls = Vec::with_capacity(columns.len());
        for at in &columns {
            let value = gt
                .and_then(|key| rest[4 + at].split(':').nth(key))
                .unwrap_or("");
            if value.is_empty() {
                calls.push(Genotype::not_called());
                continue;
            }
            let name = held[*at];
            let (call, largest) = Genotype::read(value).map_err(|unread| {
                ReadError::at(
                    line,
                    match unread {
                        Unread::Malformed => format!(
                            "sample {name}'s GT is {value:?}, and an allele is a number or ."
                        ),
                        Unread::TooMany(copies) => format!(
                            "sample {name}'s GT has {copies} copies, and a call is read with \
                             at most {}",
                            Genotype::MAX_COPIES
                        ),
                    },
                )
            })?;
            if largest as usize > alternates.len() {
                return Err(ReadError::at(
                    line,
                    format!(
                        "sample {name}'s GT {value} names allele {largest}, and the line has \
                         {} alternate allele{}",
                        alternates.len(),
                        if alternates.len() == 1 { "" } else { "s" }
                    ),
                ));
            }
            calls.push(call);
        }
        sites.push(GenotypeSite {
            position: pos,
            id: (id != ".").then(|| id.to_string()),
            reference: reference.to_string(),
            alternates,
            calls,
        });
    }
    if !sites.is_empty() && !typed {
        return Err(ReadError::whole(
            "no row in the window carries GT, so no sample's genotype is known there",
        ));
    }
    Ok(Genotypes {
        samples,
        sites,
        records,
    })
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
    /// The name of each variant in the window, as `(position, name)`, where
    /// the table has a column of them: `SNP` as PLINK and BOLT write it, `ID`
    /// as PLINK 2 and REGENIE do. Empty for a table with none.
    pub names: Vec<(u64, String)>,
}

impl Associations {
    /// The name of the variant at a 0-based position, where the table gives
    /// one and it is not a placeholder.
    pub fn name_at(&self, pos: u64) -> Option<&str> {
        self.names
            .iter()
            .find(|(at, _)| *at == pos)
            .map(|(_, name)| name.as_str())
            .filter(|name| !matches!(*name, "" | "." | "NA"))
    }
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
    let read = association_rows(text, Rows::In(region))?;
    Ok(Associations {
        points: read.points.into_iter().map(|(_, point)| point).collect(),
        p_values: read.p_values,
        names: read.names,
    })
}

/// Every row of an association table, on every sequence it names, as a scan
/// drawn across a whole genome reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomeAssociations {
    /// Each sequence the table names, in the order it first names them, with
    /// its points, 0-based on that sequence.
    pub sequences: Vec<(String, Vec<Association>)>,
    /// Whether the value column held p-values, converted on the way in.
    pub p_values: bool,
}

/// Reads every row of an association table, on every sequence, for a scan
/// drawn across a whole genome rather than over one window of it.
///
/// The table is read as [`association_table`] reads it, header, p-values and
/// all; only no row is left out for being somewhere else.
///
/// # Errors
///
/// What [`association_table`] refuses, and a table that names no sequence,
/// since a position on no sequence has no place on a genome.
///
/// ```
/// use karyon::read::point::genome_associations;
///
/// let table = "CHR BP P\n2 150 0.5\n1 300 1e-9\n1 900 0.01\n";
/// let read = genome_associations(table)?;
/// let names: Vec<&str> = read.sequences.iter().map(|(name, _)| name.as_str()).collect();
/// assert_eq!(names, ["2", "1"]);
/// assert_eq!(read.sequences[1].1.len(), 2);
/// assert!(read.p_values);
/// # Ok::<(), karyon::read::ReadError>(())
/// ```
pub fn genome_associations(text: &str) -> Result<GenomeAssociations, ReadError> {
    let read = association_rows(text, Rows::All)?;
    let mut sequences: Vec<(String, Vec<Association>)> = Vec::new();
    for (sequence, point) in read.points {
        let Some(sequence) = sequence else {
            return Err(ReadError::whole(
                "the table names no sequence, and a scan is drawn across a genome by the \
                 sequence each row is on; give it a column of them, or name the sequence \
                 as the place",
            ));
        };
        match sequences.iter_mut().find(|(named, _)| named == sequence) {
            Some((_, points)) => points.push(point),
            None => sequences.push((sequence.to_string(), vec![point])),
        }
    }
    Ok(GenomeAssociations {
        sequences,
        p_values: read.p_values,
    })
}

/// Which rows of a table are read: the ones in one window of one sequence, or
/// every one, on whatever sequence it names.
#[derive(Debug, Clone, Copy)]
enum Rows<'a> {
    In(&'a Region),
    All,
}

impl Rows<'_> {
    /// Whether a row on `sequence`, where it names one, is read at all.
    fn on(self, sequence: Option<&str>) -> bool {
        match self {
            Rows::In(region) => sequence.map_or(true, |name| name == region.seq()),
            Rows::All => true,
        }
    }

    /// Whether a row at `pos` on a sequence it is read on is kept.
    fn at(self, pos: u64) -> bool {
        match self {
            Rows::In(region) => region.contains(pos),
            Rows::All => true,
        }
    }
}

/// What a table held, each point with the sequence its row names.
struct Read<'t> {
    points: Vec<(Option<&'t str>, Association)>,
    p_values: bool,
    names: Vec<(u64, String)>,
}

/// Reads the rows `rows` asks for out of an association table.
fn association_rows<'t>(text: &'t str, rows: Rows<'_>) -> Result<Read<'t>, ReadError> {
    // A table of more columns, as association tools write them, is read by
    // the names its header gives the columns.
    if let Some((line, wide)) = hashed_header(text) {
        return wide.read(text, rows, line);
    }
    if let Some((line, head)) = lines(text).next() {
        let names = columns(head);
        if names.len() > 3 {
            return match Wide::of(&names) {
                Some(wide) => wide.read(text, rows, line),
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
    let mut points: Vec<(Option<&str>, u64, f64, usize)> = Vec::new();
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

        let sequence = sequence.map(str::trim);
        if !rows.on(sequence) {
            continue;
        }
        let pos = position(at, "position", line)?;
        if !rows.at(pos) {
            continue;
        }
        points.push((sequence, pos, number::<f64>(value, "value", line)?, line));
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
    for (sequence, pos, value, line) in points {
        if !p_values {
            converted.push((sequence, Association::new(pos, value)));
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
        converted.push((sequence, Association::from_p_value(pos, value)));
    }
    Ok(Read {
        points: converted,
        p_values,
        names: Vec::new(),
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
    id: Option<usize>,
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
        let id = find(&[
            "snp",
            "id",
            "rsid",
            "rs",
            "snpid",
            "marker",
            "markerid",
            "markername",
            "variant_id",
        ]);
        Some(Wide {
            sequence,
            position,
            value,
            p_values,
            id,
        })
    }

    /// The rows after the header, which is on line `header`.
    fn read<'t>(
        &self,
        text: &'t str,
        rows: Rows<'_>,
        header: usize,
    ) -> Result<Read<'t>, ReadError> {
        let widest = self
            .position
            .max(self.value)
            .max(self.sequence.unwrap_or(0));
        let mut points = Vec::new();
        let mut names = Vec::new();
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
            let sequence = self.sequence.map(|at| fields[at].trim());
            if !rows.on(sequence) {
                continue;
            }
            let pos = position(fields[self.position].trim(), "position", line)?;
            if !rows.at(pos) {
                continue;
            }
            // A test the tool could not run is written as NA, and has nothing
            // to draw.
            let value = fields[self.value].trim();
            if matches!(value, "" | "." | "-" | "NA" | "na" | "nan" | "NaN") {
                continue;
            }
            let value: f64 = number(value, "value", line)?;
            if let Some(name) = self.id.and_then(|at| fields.get(at)) {
                names.push((pos, name.trim().to_string()));
            }
            if !self.p_values {
                points.push((sequence, Association::new(pos, value)));
                continue;
            }
            if !(0.0..=1.0).contains(&value) {
                return Err(ReadError::at(
                    line,
                    format!("a p-value lies between 0 and 1, not {value}"),
                ));
            }
            points.push((sequence, Association::from_p_value(pos, value)));
        }
        Ok(Read {
            points,
            p_values: self.p_values,
            names,
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

/// The fraction each alternate allele gets, `None` for one that has none.
///
/// `AF` carries one number per alternate allele, so a multi-allelic row hands
/// each alternate its own. A single number is spread over all of them, since
/// that is what a caller writing one `AF` for the row means. Any other count
/// is a row that does not add up, and saying so beats guessing which allele
/// the numbers belong to.
///
/// A `.` is VCF's missing value, which bcftools writes for a fraction it could
/// not work out, for the whole row as `AF=.` or for one allele as `AF=0.5,.`.
/// That allele has no fraction, as a row with no `AF` has none, and is drawn
/// full height. Read as a number, it refused the whole track over one row.
fn allele_fractions(
    info: &str,
    alternates: usize,
    line: usize,
) -> Result<Option<Vec<Option<f64>>>, ReadError> {
    let Some(field) = info_field(info, "AF") else {
        // A call with no fraction is a call, and it gets a full height stem.
        return Ok(None);
    };
    let mut values = Vec::with_capacity(alternates);
    for value in field.split(',') {
        values.push(match value {
            "." => None,
            value => Some(number::<f64>(value, "AF", line)?),
        });
    }
    match values.len() {
        1 => Ok(Some(vec![values[0]; alternates])),
        given if given == alternates => Ok(Some(values)),
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

    /// A call with no `AF` is still a call, drawn full height, and it has no
    /// value: a value of one put an axis beside it, and then a title on the
    /// axis, saying the call was a fraction of one.
    #[test]
    fn a_call_with_no_fraction_is_still_a_call() {
        let text = "chrIV\t900\t.\tC\tT\t.\t.\tDP=30\n";
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants(text, &region).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].value, None);
    }

    /// `.` is the missing value bcftools writes for a fraction, for the row or
    /// for one allele: that allele has no fraction, as a row with no `AF` has
    /// none, and the others keep theirs. A word that is not the missing value
    /// is still refused.
    #[test]
    fn a_missing_fraction_is_no_fraction() {
        let region = Region::parse("chrIV:800-1000").unwrap();
        let calls = variants("chrIV\t900\t.\tC\tT\t.\t.\tAF=.\n", &region).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].value, None);
        let calls = variants("chrIV\t900\t.\tC\tT,G\t.\t.\tAF=.\n", &region).unwrap();
        assert_eq!(calls.len(), 2);
        assert!(calls.iter().all(|call| call.value.is_none()));
        let calls = variants("chrIV\t900\t.\tC\tT,G\t.\t.\tAF=0.5,.\n", &region).unwrap();
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].value, Some(0.5));
        assert_eq!(calls[1].value, None);
        let calls = variants("chrIV\t900\t.\tC\tT,G\t.\t.\tAF=.,0.25\n", &region).unwrap();
        assert_eq!(calls[0].value, None);
        assert_eq!(calls[1].value, Some(0.25));
        for refused in ["AF=NA", "AF=", "AF=.,.,.", "AF=..", "AF=0.5,"] {
            let text = format!("chrIV\t900\t.\tC\tT,G\t.\t.\t{refused}\n");
            assert!(variants(&text, &region).is_err(), "{refused}");
        }
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
    /// The tools that write a variant's name write it under their own names
    /// for the column, and a table with none names nothing.
    #[test]
    fn a_variant_is_named_by_the_column_its_tool_writes_names_in() {
        let region = Region::parse("1:1-1,000").unwrap();
        for table in [
            "CHR SNP BP A1 P\n1 rs1 150 A 0.5\n1 rs2 480 A 1e-9\n",
            "#CHROM POS ID P\n1 150 rs1 0.5\n1 480 rs2 1e-9\n",
            "CHROM GENPOS ID LOG10P\n1 150 rs1 0.3\n1 480 rs2 9\n",
        ] {
            let read = association_table(table, &region).unwrap();
            assert_eq!(read.name_at(479), Some("rs2"), "{table}");
            assert_eq!(read.name_at(100), None);
        }
        let unnamed = association_table("CHR BP A1 P\n1 150 A 0.5\n", &region).unwrap();
        assert!(unnamed.names.is_empty());
        let placeholder = association_table("CHR SNP BP A1 P\n1 . 150 A 0.5\n", &region).unwrap();
        assert_eq!(placeholder.name_at(149), None, "a dot names nothing");
    }

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

    /// A cohort of three, the header every genotype test reads under.
    const COHORT: &str = "##fileformat=VCFv4.2\n\
                          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n";

    fn cohort(rows: &str) -> String {
        format!("{COHORT}{rows}")
    }

    fn window() -> Region {
        Region::parse("chr1:1-1,000").unwrap()
    }

    /// The states of the first site's calls, in row order.
    fn states(rows: &str) -> Vec<crate::GenotypeState> {
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        read.sites[0]
            .calls
            .iter()
            .map(|call| call.state())
            .collect()
    }

    fn refusal(rows: &str) -> ReadError {
        genotypes(&cohort(rows), &window(), None).unwrap_err()
    }

    #[test]
    fn the_samples_are_the_columns_after_format_on_the_chrom_line() {
        assert_eq!(samples(COHORT), ["S1", "S2", "S3"]);
        // A VCF of sites alone names none, and neither does text with no
        // header at all.
        assert!(samples(RPOB).is_empty());
        assert!(samples("chr1\t100\t.\tC\tT\t.\t.\t.\n").is_empty());
        // The header is read and the rows are not: a row past it that looks
        // like a header changes nothing.
        let later = format!("{COHORT}chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\t1\t0\n#CHROM\tX\n");
        assert_eq!(samples(&later), ["S1", "S2", "S3"]);
    }

    #[test]
    fn a_haploid_one_and_a_diploid_one_one_are_both_all_alternate() {
        let read = genotypes(
            &cohort("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t1\t1/1\t0\n"),
            &window(),
            None,
        )
        .unwrap();
        let shares: Vec<Option<f64>> = read.sites[0].calls.iter().map(|c| c.share()).collect();
        assert_eq!(shares, [Some(1.0), Some(1.0), Some(0.0)]);
        assert_eq!(read.sites[0].call(0).copies(), 1);
        assert_eq!(read.sites[0].call(1).copies(), 2);
    }

    #[test]
    fn zero_one_is_half_and_zero_zero_zero_one_a_quarter() {
        let read = genotypes(
            &cohort("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0/0/1\t0/0\n"),
            &window(),
            None,
        )
        .unwrap();
        let shares: Vec<Option<f64>> = read.sites[0].calls.iter().map(|c| c.share()).collect();
        assert_eq!(shares, [Some(0.5), Some(0.25), Some(0.0)]);
    }

    /// Half a genotype is no call. Counted over the copies that were called,
    /// `./1` would be all alternate and drawn as a homozygote.
    #[test]
    fn a_partly_called_genotype_is_not_called() {
        use crate::GenotypeState::NotCalled;
        assert_eq!(
            states("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t./1\t0/.\t.|1\n"),
            [NotCalled, NotCalled, NotCalled]
        );
        let half = crate::Genotype::parse("./1").unwrap();
        assert_eq!((half.called(), half.alternate()), (1, 1));
        assert_eq!(half.share(), None);
    }

    /// A caller writes `.` for a sample it has nothing for, `./.` for a
    /// diploid one, or leaves a field empty, and keys past GT may be dropped.
    #[test]
    fn a_missing_sample_field_is_not_called() {
        use crate::GenotypeState::{NotCalled, Reference};
        assert_eq!(
            states("chr1\t100\t.\tC\tT\t.\t.\t.\tGT:DP\t.\t./.\t\n"),
            [NotCalled, NotCalled, NotCalled]
        );
        // GT second, and the whole field written `.`: there is no GT in it.
        assert_eq!(
            states("chr1\t100\t.\tC\tT\t.\t.\t.\tDP:GT\t.\t12:0/0\t3\n"),
            [NotCalled, Reference, NotCalled]
        );
    }

    #[test]
    fn phased_and_unphased_read_alike_and_phase_is_kept() {
        let read = genotypes(
            &cohort("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0|1\t0/1\t|0|1\n"),
            &window(),
            None,
        )
        .unwrap();
        let calls = &read.sites[0].calls;
        assert_eq!(calls[0].share(), calls[1].share());
        assert!(calls[0].phased());
        assert!(!calls[1].phased());
        // VCF 4.4's leading mark is the phase of the first copy, not a copy.
        assert_eq!(calls[2], calls[0]);
    }

    #[test]
    fn gt_need_not_be_the_first_format_key() {
        use crate::GenotypeState::{Alternate, Heterozygous, Reference};
        assert_eq!(
            states("chr1\t100\t.\tC\tT\t.\t.\t.\tDP:AD:GT\t9:4,5:0/1\t9:0,9:1/1\t9:9,0:0/0\n"),
            [Heterozygous, Alternate, Reference]
        );
    }

    #[test]
    fn a_row_without_gt_leaves_every_sample_not_called() {
        use crate::GenotypeState::{NotCalled, Reference};
        let rows = "chr1\t100\t.\tC\tT\t.\t.\t.\tDP\t9\t9\t9\n\
                    chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0\t0\t0\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        assert_eq!(read.sites.len(), 2);
        assert!(read.sites[0].calls.iter().all(|c| c.state() == NotCalled));
        assert!(read.sites[1].calls.iter().all(|c| c.state() == Reference));
    }

    /// A window of records that none of carries GT says nothing of anyone,
    /// and drawn it was a band of no calls that looked like a failed sample.
    #[test]
    fn a_window_with_no_gt_anywhere_is_refused() {
        let error = refusal("chr1\t100\t.\tC\tT\t.\t.\t.\tDP\t9\t9\t9\n");
        assert!(
            error
                .to_string()
                .contains("no row in the window carries GT"),
            "{error}"
        );
        // Outside the window, a row without GT is nobody's business.
        let elsewhere = cohort("chr2\t100\t.\tC\tT\t.\t.\t.\tDP\t9\t9\t9\n");
        assert!(genotypes(&elsewhere, &window(), None)
            .unwrap()
            .sites
            .is_empty());
    }

    #[test]
    fn an_allele_past_the_alternates_is_refused_on_its_line() {
        let error = refusal("chr1\t100\t.\tC\tT,G\t.\t.\t.\tGT\t0/1\t0/3\t0/0\n");
        assert_eq!(error.line, 3);
        assert!(
            error.to_string().contains(
                "sample S2's GT 0/3 names allele 3, and the line has 2 alternate alleles"
            ),
            "{error}"
        );
    }

    #[test]
    fn a_genotype_that_is_not_numbers_is_refused_on_its_line() {
        let error = refusal("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/x\t0/0\n");
        assert_eq!(error.line, 3);
        assert!(
            error
                .to_string()
                .contains("sample S2's GT is \"0/x\", and an allele is a number or ."),
            "{error}"
        );
        // A copy with nothing in it is no copy.
        assert!(refusal("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0/\t0\t0\n")
            .to_string()
            .contains("S1's GT"));
    }

    #[test]
    fn a_line_short_of_samples_is_refused_naming_both_counts() {
        let error = refusal("chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\n");
        assert_eq!(error.line, 3);
        assert!(
            error
                .to_string()
                .contains("this line has 2 samples and the #CHROM line names 3"),
            "{error}"
        );
    }

    #[test]
    fn a_sites_only_vcf_is_refused_and_names_variants() {
        let error = genotypes(RPOB, &rpob_region(), None).unwrap_err();
        assert!(
            error.to_string().contains("--variants draws its calls"),
            "{error}"
        );
        // And text with no header at all says to keep it.
        let bare = "chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\n";
        let error = genotypes(bare, &window(), None).unwrap_err();
        assert!(error.to_string().contains("no #CHROM line"), "{error}");
    }

    /// A gVCF's blocks name no allele but the placeholder, `<NON_REF>` from
    /// GATK and `<*>` from bcftools, and are not sites.
    #[test]
    fn a_reference_block_is_not_a_site() {
        let rows = "chr1\t100\t.\tC\t.\t.\t.\tEND=150\tGT\t0\t0\t0\n\
                    chr1\t200\t.\tC\t<NON_REF>\t.\t.\tEND=250\tGT\t0/0\t0/0\t0/0\n\
                    chr1\t300\t.\tC\t<*>\t.\t.\tEND=350\tGT\t0/0\t0/0\t0/0\n\
                    chr1\t400\t.\tC\tT,<NON_REF>\t.\t.\t.\tGT\t0/1\t0/0\t1/2\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        assert_eq!(read.sites.len(), 1);
        let site = &read.sites[0];
        assert_eq!(site.position, 399);
        // The placeholder stays in the list, so `2` still names it.
        assert_eq!(site.alternates, ["T", "<NON_REF>"]);
        assert_eq!(site.call(2).allele(1), Some(2));
        assert_eq!(read.records, 4);
    }

    /// A GATK variant row names the placeholder after its own allele, and
    /// is the call of that allele alone; a block of only the placeholder is
    /// no call. Dropping the row's placeholder, and not the allele after it,
    /// keeps a fraction given per allele with its own allele.
    #[test]
    fn variants_drawn_from_a_gvcf_do_not_include_its_reference_blocks() {
        let text = "chr1\t200\t.\tC\t<NON_REF>\t.\t.\tEND=250\n\
                    chr1\t280\t.\tC\t<*>\t.\t.\tEND=290\n\
                    chr1\t300\t.\tC\tT,<NON_REF>\t.\t.\tAF=0.4,0.0\n\
                    chr1\t320\t.\tC\t<NON_REF>,G\t.\t.\tAF=0.0,0.7\n";
        let calls = variants(text, &window()).unwrap();
        let drawn: Vec<(u64, Option<&str>, Option<f64>)> = calls
            .iter()
            .map(|call| (call.pos, call.category.as_deref(), call.value))
            .collect();
        assert_eq!(
            drawn,
            [
                (299, Some("substitution"), Some(0.4)),
                (319, Some("substitution"), Some(0.7))
            ]
        );
    }

    #[test]
    fn a_multi_allelic_row_is_one_site_and_one_two_has_no_reference_copy() {
        use crate::GenotypeState::{Alternate, Heterozygous};
        let rows = "chr1\t100\t.\tC\tT,G\t.\t.\t.\tGT\t1/2\t0/2\t2/2\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        assert_eq!(read.sites.len(), 1);
        let calls: Vec<_> = read.sites[0].calls.iter().map(|c| c.state()).collect();
        assert_eq!(calls, [Alternate, Heterozygous, Alternate]);
        assert_eq!(read.sites[0].call(0).allele(0), Some(1));
        assert_eq!(read.sites[0].call(0).allele(1), Some(2));
    }

    /// A joint caller writes `*` for the base a deletion upstream took away,
    /// and a gVCF's variant row keeps its placeholder: a copy naming either
    /// is not the reference, so it is an alternate copy like any other, and
    /// the indices still name the alleles they were written against.
    #[test]
    fn a_call_of_a_star_or_a_placeholder_is_an_alternate_copy() {
        use crate::GenotypeState::{Alternate, Heterozygous, Reference};
        let rows = "chr1\t100\t.\tA\tG,*\t.\t.\t.\tGT\t2/2\t0/2\t0/0\n\
                    chr1\t200\t.\tC\tT,<NON_REF>\t.\t.\t.\tGT\t1/2\t0/2\t2\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        assert_eq!(read.sites.len(), 2);
        let states = |site: usize| -> Vec<crate::GenotypeState> {
            read.sites[site].calls.iter().map(|c| c.state()).collect()
        };
        assert_eq!(states(0), [Alternate, Heterozygous, Reference]);
        assert_eq!(states(1), [Alternate, Heterozygous, Alternate]);
        assert_eq!(read.sites[0].alternates[1], "*");
        assert_eq!(read.sites[1].call(2).allele(0), Some(2));
    }

    #[test]
    fn wanted_samples_come_back_in_the_order_asked() {
        let rows = "chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\t1\t.\n";
        let wanted = ["S3".to_string(), "S1".to_string()];
        let read = genotypes(&cohort(rows), &window(), Some(&wanted)).unwrap();
        assert_eq!(read.samples, ["S3", "S1"]);
        let states: Vec<_> = read.sites[0].calls.iter().map(|c| c.state()).collect();
        assert_eq!(
            states,
            [
                crate::GenotypeState::NotCalled,
                crate::GenotypeState::Reference
            ]
        );
        let unknown = ["S9".to_string()];
        let error = genotypes(&cohort(rows), &window(), Some(&unknown)).unwrap_err();
        assert!(error.to_string().contains("no sample called S9"), "{error}");
        let twice = ["S1".to_string(), "S1".to_string()];
        let error = genotypes(&cohort(rows), &window(), Some(&twice)).unwrap_err();
        assert!(
            error.to_string().contains("S1 is asked for twice"),
            "{error}"
        );
    }

    #[test]
    fn a_sample_named_twice_on_the_chrom_line_is_refused() {
        let text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS1\n\
                    chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\t1\t0\n";
        let error = genotypes(text, &window(), None).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.to_string().contains("names S1 twice"), "{error}");
    }

    /// The rule `variants` keeps a row by: what REF spells touching the
    /// window, so a deletion anchored a base before it is a site of it.
    #[test]
    fn a_deletion_whose_ref_reaches_into_the_window_is_kept() {
        let rows = "chr1\t98\t.\tCTTT\tC\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n\
                    chr1\t90\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n";
        let region = Region::parse("chr1:100-200").unwrap();
        let read = genotypes(&cohort(rows), &region, None).unwrap();
        let positions: Vec<u64> = read.sites.iter().map(|site| site.position).collect();
        assert_eq!(positions, [97]);
        assert_eq!(
            variants(&cohort(rows), &region).unwrap().len(),
            1,
            "the same rule keeps the same rows"
        );
    }

    /// The largest POS a file can write is one base from the top of the
    /// number line, and what REF spells from there is past it. Both readers
    /// panicked on the add; the row is outside any window that ends before
    /// it, and is passed over as one.
    #[test]
    fn a_position_at_the_top_of_the_number_line_is_outside_the_window_and_not_a_panic() {
        let rows = "chr1\t18446744073709551615\t.\tCC\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n\
                    chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        let positions: Vec<u64> = read.sites.iter().map(|site| site.position).collect();
        assert_eq!(positions, [99]);
        let calls = variants(&cohort(rows), &window()).unwrap();
        let positions: Vec<u64> = calls.iter().map(|call| call.pos).collect();
        assert_eq!(positions, [99]);
    }

    /// A row outside the window is not split past its position, and a row on
    /// another sequence is not checked at all, as `variants` does not.
    #[test]
    fn rows_outside_the_window_are_not_read_past_where_they_are() {
        let rows = "chr2\t100\t.\tC\tT\t.\t.\t.\tGT\t0/x\n\
                    chr1\t5000\t.\tC\tT\t.\t.\t.\tGT\t0/x\n\
                    chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t0\t1\t0\n";
        let read = genotypes(&cohort(rows), &window(), None).unwrap();
        assert_eq!(read.sites.len(), 1);
        assert_eq!(read.records, 3);
    }
}
