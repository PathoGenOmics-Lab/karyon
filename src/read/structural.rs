//! Structural variants, from the VCF a caller writes them in.
//!
//! ```text
//! chrA  321682  .  T  <DEL>  6  PASS  SVTYPE=DEL;SVLEN=205;END=321887
//! chrA  321687  .  T  t[chrB:123457[  6  PASS  SVTYPE=BND;MATEID=bnd_W
//! ```
//!
//! # POS is the base before the event, so nothing is taken off it
//!
//! This is the one place in the directory where the 1-based to 0-based
//! subtraction the other readers do is wrong. A symbolic allele such as `<DEL>`
//! cannot be written against no reference base, so the specification puts POS on
//! the base *before* the event and REF on that base. The event therefore starts
//! at 1-based `POS + 1`, which counted from nought is `POS`, and END is the last
//! base of it, which counted from nought and made exclusive is `END`.
//!
//! Both conversions are the identity, and they are the identity for two
//! different reasons. Copying [`point`](super::point), which does take one off a
//! `POS`, would move every call two bases left.
//!
//! # A length that is absent is not a length of one
//!
//! A span comes from `SVLEN` where there is one, from `END` where there is not,
//! and from the length of REF for a variant spelled out in full. Where none of
//! the three says anything the row is refused, because the remaining answer is
//! `POS..POS + 1`, and that is a one-base call drawn at full confidence out of a
//! file that stated no length at all.
//!
//! `SVLEN` is taken as its absolute value. VCF 4.3 writes a deletion's length
//! negative and 4.4 writes it positive, and a reader that believes the sign on a
//! 4.3 file computes an end before its own start.
//!
//! # Two breakends are one rearrangement
//!
//! A `BND` record is one end of a join, and the other end is a second record.
//! The mate is in the ALT, as `t[chrB:123457[` and its three siblings, so the
//! pair is read off one row rather than chased through `MATEID`, and only the
//! row whose own position is the lower of the two becomes an arc. Otherwise the
//! reciprocal record draws the same arc again.
//!
//! A breakend whose mate is on another sequence cannot be drawn on an axis of
//! this one, and neither can a single breakend, whose other end is unknown by
//! definition. Both are counted rather than turned into a variant that starts
//! and ends in the same place, which is a glyph this crate already spends on a
//! confidently located insertion.

use crate::{Region, StructuralVariant, SvKind};

use super::{columns, lines, number, ReadError};

/// The calls a VCF holds, and what it held that is not among them.
#[derive(Debug, Clone, PartialEq)]
pub struct Calls {
    /// The variants, in the order the file listed them.
    pub variants: Vec<StructuralVariant>,
    /// Records in the file, before any filter.
    pub records: usize,
    /// Records naming another sequence.
    pub passed_over: usize,
    /// Records on this sequence that touch no base of the window.
    pub off_region: usize,
    /// Records that are not structural at all, having no symbolic allele and
    /// no `SVTYPE`.
    pub not_structural: usize,
    /// Breakends with nowhere to draw the other end.
    ///
    /// A mate on another sequence, or a single breakend, whose other end is
    /// unknown by definition. Not an error, and not a variant either: a
    /// rearrangement drawn with both ends in one place is a claim about a
    /// locus rather than about a join.
    pub half_a_join: usize,
    /// Breakend records folded onto the reciprocal record of the same join.
    pub reciprocal: usize,
    /// Records whose class this crate has no glyph for, `<CNV>` above all.
    ///
    /// A copy-number-variable region is not a gain, and drawing it as a
    /// duplication states which way it went when the file did not.
    pub unclassified: usize,
}

/// Reads structural calls out of a VCF.
///
/// `POS` is the base before a symbolic event, so it becomes the 0-based start
/// unchanged, and `END` becomes the exclusive end unchanged. The class comes
/// from the ALT allele where it is symbolic and from `SVTYPE` where it is not.
///
/// Rows naming another sequence than `region.seq()` are skipped.
///
/// # Errors
///
/// Returns the first row that is not a VCF row, whose `END` is before its `POS`,
/// whose length is nought, or which names a symbolic allele and no length at
/// all. A file that parses and holds no call inside the window is not an error
/// here: [`Calls`] says why, and the caller decides.
pub fn variants(text: &str, region: &Region) -> Result<Calls, ReadError> {
    let mut found = Calls {
        variants: Vec::new(),
        records: 0,
        passed_over: 0,
        off_region: 0,
        not_structural: 0,
        half_a_join: 0,
        reciprocal: 0,
        unclassified: 0,
    };

    for (at, line) in lines(text) {
        let cols = columns(line);
        if cols.len() < 8 {
            return Err(ReadError::at(
                at,
                format!(
                    "a VCF row has at least 8 columns, this one has {}",
                    cols.len()
                ),
            ));
        }
        found.records += 1;
        if cols[0] != region.seq() {
            found.passed_over += 1;
            continue;
        }

        let pos: u64 = number(cols[1], "POS", at)?;
        let info = cols[7];
        let alt = cols[4];

        let Some(kind) = kind(alt, info) else {
            // Either an ordinary small variant, which belongs in a variant
            // track, or a class this crate draws no glyph for.
            if symbolic(alt).is_some() || key(info, "SVTYPE").is_some() {
                found.unclassified += 1;
            } else {
                found.not_structural += 1;
            }
            continue;
        };

        let (start, end) = if kind == SvKind::Translocation {
            let Some((mate_seq, mate_pos)) = breakend(alt) else {
                found.half_a_join += 1;
                continue;
            };
            if mate_seq != region.seq() {
                found.half_a_join += 1;
                continue;
            }
            // One arc per join. The reciprocal record describes the same two
            // places from the other side, and drawn as well it doubles every
            // stroke on the figure.
            if mate_pos < pos {
                found.reciprocal += 1;
                continue;
            }
            (pos, mate_pos)
        } else {
            (pos, end(pos, cols[3], alt, info, kind, at)?)
        };

        if end < start {
            return Err(ReadError::at(
                at,
                format!("the call ends at {end} and starts at {start}"),
            ));
        }
        let mut variant = StructuralVariant::new(start, end, kind);
        // A support count that was never stated is not a count of nought.
        // Nought reads is the thinnest arc the track draws and a tooltip that
        // says no read supported the call, which is a finding.
        if let Some(reads) = support(info) {
            variant = variant.support(reads);
        }
        if let Some(name) = Some(cols[2]).filter(|id| *id != "." && !id.is_empty()) {
            variant = variant.name(name);
        }

        if !variant.touches(region.start(), region.end()) {
            found.off_region += 1;
            continue;
        }
        found.variants.push(variant);
    }

    Ok(found)
}

/// The class of a call, from its ALT where it can be and its `SVTYPE` after.
///
/// The ALT first because VCF 4.4 deprecated `SVTYPE`, and a bracketed ALT is
/// a breakend whatever the INFO column says.
fn kind(alt: &str, info: &str) -> Option<SvKind> {
    if alt.contains('[') || alt.contains(']') {
        return Some(SvKind::Translocation);
    }
    let word = symbolic(alt)
        .map(|word| word.split(':').next().unwrap_or(word).to_string())
        .or_else(|| key(info, "SVTYPE"))?;
    Some(match word.as_str() {
        "DEL" => SvKind::Deletion,
        "DUP" => SvKind::Duplication,
        "INV" => SvKind::Inversion,
        "INS" => SvKind::Insertion,
        "BND" | "TRA" => SvKind::Translocation,
        // `<CNV>` says the copy number varies and not which way, and this
        // crate has a glyph for a gain and none for either.
        _ => return None,
    })
}

/// The inside of a symbolic allele, or `None` for a spelled-out one.
fn symbolic(alt: &str) -> Option<&str> {
    alt.strip_prefix('<')?.strip_suffix('>')
}

/// The exclusive end of a call that covers reference.
fn end(
    pos: u64,
    reference: &str,
    alt: &str,
    info: &str,
    kind: SvKind,
    at: usize,
) -> Result<u64, ReadError> {
    // An insertion is sequence that is not in the reference, so it has one
    // breakpoint and no footprint. Its SVLEN is how much was inserted, and
    // spending that on the reference would draw over bases the file said
    // nothing about.
    if kind == SvKind::Insertion {
        return Ok(pos);
    }

    let stated = key(info, "SVLEN")
        .as_deref()
        .map(|word| {
            word.parse::<f64>()
                .map(f64::abs)
                .map_err(|_| ReadError::at(at, format!("SVLEN is not a number: {word:?}")))
        })
        .transpose()?;
    let ended = key(info, "END")
        .as_deref()
        .map(|word| number::<u64>(word, "END", at))
        .transpose()?;

    let end = match (stated, ended) {
        (Some(length), Some(end)) => {
            // A record disagreeing with itself is a record two callers wrote
            // half of, and picking the one that suits is how a figure ends up
            // stating a length no file holds.
            if (pos + length.round() as u64) != end {
                return Err(ReadError::at(
                    at,
                    format!(
                        "SVLEN puts the end at {}, and END says {end}",
                        pos + length.round() as u64
                    ),
                ));
            }
            end
        }
        (Some(length), None) => pos + length.round() as u64,
        (None, Some(end)) => end,
        (None, None) => {
            // A symbolic allele carries no sequence, so nothing but SVLEN or
            // END can say how far the event reached.
            if symbolic(alt).is_some() || reference == "." || reference.is_empty() {
                return Err(ReadError::at(
                    at,
                    "a symbolic call carries neither SVLEN nor END, so it states no length",
                ));
            }
            // Spelled out in full: the reference allele is the footprint, and
            // POS is its first base rather than the base before it.
            pos + reference.len() as u64 - 1
        }
    };
    if end == pos {
        return Err(ReadError::at(at, "the call covers no reference bases"));
    }
    Ok(end)
}

/// The mate of a breakend, out of its ALT.
///
/// The four spellings are `t[seq:pos[`, `t]seq:pos]`, `[seq:pos[t` and
/// `]seq:pos]t`, which differ in which side of the join is kept and not in
/// where the mate is. A single breakend is `t.` or `.t`, which names no mate.
fn breakend(alt: &str) -> Option<(String, u64)> {
    let inner = alt.split(['[', ']']).nth(1)?;
    let (seq, pos) = inner.rsplit_once(':')?;
    // Some callers write a contig in angle brackets even inside a breakend.
    let seq = seq.trim_start_matches('<').trim_end_matches('>');
    Some((seq.to_string(), pos.parse().ok()?))
}

/// How many reads a caller said supported the call, if it said.
fn support(info: &str) -> Option<u32> {
    ["SUPPORT", "PE", "SR", "RE", "DV"]
        .iter()
        .find_map(|name| key(info, name))
        .and_then(|word| word.parse().ok())
}

/// One `KEY=value` out of the INFO column, or a flag's own name.
fn key(info: &str, name: &str) -> Option<String> {
    info.split(';').find_map(|pair| {
        let (found, value) = pair.trim().split_once('=')?;
        (found.trim() == name).then(|| value.trim().to_string())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn window() -> Region {
        Region::new("chrA", 0, 500_000).unwrap()
    }

    #[test]
    fn pos_is_the_base_before_the_event_so_nothing_is_taken_off_it() {
        // The specification's own worked example: POS 2, SVLEN 2, END 4 is the
        // deletion of the two bases 1-based 3 and 4. Counted from nought that
        // is 2..4, and both conversions are the identity for different reasons.
        let text = "chrA\t2\t.\tT\t<DEL>\t6\tPASS\tSVTYPE=DEL;SVLEN=2;END=4\n";
        let found = variants(text, &Region::new("chrA", 0, 100).unwrap()).unwrap();
        assert_eq!(found.variants[0].start, 2);
        assert_eq!(found.variants[0].end, 4);
        assert_eq!(found.variants[0].span(), 2, "the span is not the SVLEN");
    }

    #[test]
    fn a_deletion_length_is_read_whichever_sign_its_spec_version_gives_it() {
        // VCF 4.3 writes it negative and 4.4 positive, and believing the sign
        // computes an end before the start on every 4.3 file in existence.
        let old = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVTYPE=DEL;SVLEN=-205\n";
        let new = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205\n";
        let a = variants(old, &window()).unwrap();
        let b = variants(new, &window()).unwrap();
        assert_eq!(a.variants[0].end, 1205);
        assert_eq!(a.variants, b.variants);
    }

    #[test]
    fn a_call_that_states_no_length_is_refused_rather_than_made_one_base() {
        // The remaining answer is POS..POS+1, which is a one-base deletion
        // drawn at full confidence out of a file that gave no length.
        let text = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVTYPE=DEL\n";
        let error = variants(text, &window()).unwrap_err();
        assert!(error.reason.contains("states no length"), "{error}");
    }

    #[test]
    fn a_record_that_disagrees_with_itself_stops_the_read() {
        let text = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205;END=9999\n";
        let error = variants(text, &window()).unwrap_err();
        assert!(error.reason.contains("END says 9999"), "{error}");

        let backwards = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tEND=100\n";
        let error = variants(backwards, &window()).unwrap_err();
        assert!(error.reason.contains("starts at 1000"), "{error}");
    }

    #[test]
    fn an_insertion_has_one_breakpoint_and_no_footprint() {
        // Its SVLEN is how much was inserted. Spent on the reference it draws a
        // hundred-base footprint over bases the file said nothing about.
        let text = "chrA\t1000\t.\tT\t<INS>\t6\tPASS\tSVTYPE=INS;SVLEN=100\n";
        let found = variants(text, &window()).unwrap();
        assert_eq!(found.variants[0].start, 1000);
        assert_eq!(found.variants[0].end, 1000);
        assert_eq!(found.variants[0].kind, SvKind::Insertion);
        assert!(!SvKind::Insertion.has_footprint());
    }

    #[test]
    fn a_class_this_crate_has_no_glyph_for_is_counted_and_not_guessed_at() {
        // A copy-number-variable region is not a gain. Drawn as a duplication
        // it states which way the copy number went, and the file did not.
        let text = "chrA\t1000\t.\tT\t<CNV>\t6\tPASS\tSVTYPE=CNV;END=2000\n";
        let found = variants(text, &window()).unwrap();
        assert!(found.variants.is_empty());
        assert_eq!(found.unclassified, 1);
    }

    #[test]
    fn a_breakend_finds_its_mate_in_its_own_alt_in_all_four_spellings() {
        for alt in [
            "t[chrA:4000[",
            "t]chrA:4000]",
            "[chrA:4000[t",
            "]chrA:4000]t",
        ] {
            let text = format!("chrA\t1000\t.\tT\t{alt}\t6\tPASS\tSVTYPE=BND\n");
            let found = variants(&text, &window()).unwrap();
            assert_eq!(found.variants.len(), 1, "{alt}");
            assert_eq!(found.variants[0].start, 1000, "{alt}");
            assert_eq!(found.variants[0].end, 4000, "{alt}");
            assert_eq!(found.variants[0].kind, SvKind::Translocation);
        }
    }

    #[test]
    fn a_join_is_one_arc_however_many_records_describe_it() {
        // Both ends of one rearrangement are in the file, as the format
        // requires. Drawn twice they double every stroke.
        let text = "\
chrA\t1000\t bnd_U\tT\tt[chrA:4000[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_V
chrA\t4000\tbnd_V\tG\t]chrA:1000]G\t6\tPASS\tSVTYPE=BND;MATEID=bnd_U
";
        let found = variants(text, &window()).unwrap();
        assert_eq!(found.variants.len(), 1);
        assert_eq!(found.reciprocal, 1);
        assert_eq!(
            (found.variants[0].start, found.variants[0].end),
            (1000, 4000)
        );
    }

    #[test]
    fn half_a_join_is_counted_rather_than_drawn_in_one_place() {
        // A mate elsewhere and a single breakend both leave one coordinate.
        // Built anyway they become start == end, which is the glyph this crate
        // already spends on a confidently located insertion.
        let text = "\
chrA\t1000\t.\tT\tt[chrZ:4000[\t6\tPASS\tSVTYPE=BND
chrA\t2000\t.\tT\tt.\t6\tPASS\tSVTYPE=BND
chrA\t3000\t.\tT\t.t\t6\tPASS\tSVTYPE=BND
";
        let found = variants(text, &window()).unwrap();
        assert!(found.variants.is_empty());
        assert_eq!(found.half_a_join, 3);
    }

    #[test]
    fn a_support_count_nobody_gave_is_absent_and_not_nought() {
        // Nought is the thinnest arc the track draws and a tooltip saying no
        // read supported the call, which is a finding rather than a silence.
        let bare = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205\n";
        assert_eq!(variants(bare, &window()).unwrap().variants[0].support, None);

        let counted = "chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205;PE=17\n";
        assert_eq!(
            variants(counted, &window()).unwrap().variants[0].support,
            Some(17)
        );
    }

    #[test]
    fn rows_elsewhere_are_counted_rather_than_dropped_in_silence() {
        let text = "\
##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chrA\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205
chrB\t1000\t.\tT\t<DEL>\t6\tPASS\tSVLEN=205
chrA\t2000\t.\tT\tTA\t6\tPASS\t.
";
        let found = variants(text, &window()).unwrap();
        assert_eq!(found.variants.len(), 1);
        assert_eq!(found.passed_over, 1);
        assert_eq!(found.not_structural, 1, "a small variant is not an SV");
        assert_eq!(found.records, 3);

        let far = variants(text, &Region::new("chrA", 300_000, 400_000).unwrap()).unwrap();
        assert!(far.variants.is_empty());
        assert_eq!(far.off_region, 1);
    }

    #[test]
    fn a_line_that_is_not_a_vcf_row_stops_the_read_and_says_which() {
        let error = variants("chrA\t1000\t.\tT\n", &window()).unwrap_err();
        assert_eq!(error.line, 1);
        assert!(error.reason.contains("at least 8 columns"), "{error}");
    }
}
