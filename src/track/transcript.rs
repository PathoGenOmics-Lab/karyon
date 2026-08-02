//! Where transcription starts, how far the leader runs, and where it stops.
//!
//! Every other interval in the crate is a piece of the reference. A
//! [`TranscriptionUnit`] is a molecule the reference does not contain, so its
//! glyph carries three statements no arrangement of gene models can make: where
//! the polymerase started, how much of the RNA runs before the first codon, and
//! how it stopped. [`Terminator::Unknown`] is the default of the third, so a 3'
//! end nobody has explained is a plain tick rather than a mechanism nobody
//! claimed.
//!
//! # The start site is the anchor, not the lower coordinate
//!
//! Several of the types nearby sort a reversed pair of coordinates on the way
//! in and carry on. This one cannot, because the two ends of a transcript are
//! not interchangeable: swapping them turns the molecule around. Everything
//! else is measured from the start site, which is why a coordinate on the wrong
//! side of it reads as a contradiction and not as a negative distance.
//!
//! # Rows are packed in pixels, not in coordinates
//!
//! A unit takes up more of the figure than its span does: the bent arrow rises
//! out of the transcript at the start site, the terminator stands over the 3'
//! end, and both are a fixed size in pixels whatever the zoom. Packing on the
//! coordinates alone would seat two transcripts a kilobase apart on one row at
//! every scale, and across a view megabases wide their glyphs would be drawn
//! through each other. The packer reserves a margin in pixels either side of
//! each unit instead, so the number of rows, and with it the height of the
//! track, changes with the width of the view.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::mix;
use crate::track::{strand_color, DrawContext, Strand, Track};

/// How a transcript is brought to an end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Terminator {
    /// A hairpin the RNA folds into by itself, drawn as a stem and a loop.
    Intrinsic,
    /// Rho has to catch up with the polymerase, drawn as a plain bar.
    RhoDependent,
    /// The 3' end is known but how it got there is not.
    #[default]
    Unknown,
}

/// One RNA molecule: where it starts, where it ends, and what it codes for.
#[derive(Debug, Clone, PartialEq)]
pub struct TranscriptionUnit {
    tss: u64,
    end: u64,
    strand: Strand,
    cds_start: Option<u64>,
    terminator: Terminator,
    name: Option<String>,
}

impl TranscriptionUnit {
    /// A transcript from its start site to its 3' end, both 0-based.
    ///
    /// On the reverse strand `tss` is the higher coordinate, because that is
    /// where transcription begins.
    pub fn new(tss: u64, end: u64, strand: Strand) -> Self {
        TranscriptionUnit {
            tss,
            end,
            strand,
            cds_start: None,
            terminator: Terminator::default(),
            name: None,
        }
    }

    /// The first base of the first codon of the first gene it carries.
    ///
    /// Setting it draws the 5' leader hollow, so its length is readable off the
    /// figure. A transcript whose `cds_start` equals its `tss` is leaderless, and
    /// how much of a collection is leaderless is the observation the figure exists
    /// to make. The fraction varies by organism from a handful to most of them.
    pub fn cds_start(mut self, at: u64) -> Self {
        self.cds_start = Some(at);
        self
    }

    /// How the transcript is terminated.
    pub fn terminator(mut self, terminator: Terminator) -> Self {
        self.terminator = terminator;
        self
    }

    /// Sets the name printed over the transcript.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The transcription start site.
    pub fn tss(&self) -> u64 {
        self.tss
    }

    /// The 3' end.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Which strand it is transcribed from.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// The reference span it covers, low coordinate first.
    pub fn extent(&self) -> (u64, u64) {
        (self.tss.min(self.end), self.tss.max(self.end))
    }

    /// How long the transcript is.
    pub fn len(&self) -> u64 {
        let (low, high) = self.extent();
        high - low
    }

    /// Whether the transcript covers no sequence at all.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bases between the start site and the start codon.
    ///
    /// `None` when no start codon was given, `Some(0)` for a leaderless
    /// transcript. A start codon on the wrong side of the start site is a
    /// contradiction rather than a negative leader, so it reads as `None`.
    pub fn leader(&self) -> Option<u64> {
        let cds = self.cds_start?;
        if self.strand == Strand::Reverse {
            self.tss.checked_sub(cds)
        } else {
            cds.checked_sub(self.tss)
        }
    }

    /// Whether translation begins at the very first transcribed base.
    pub fn is_leaderless(&self) -> bool {
        self.leader() == Some(0)
    }
}

/// Transcription units, one glyph set per RNA molecule.
///
/// Nothing else in the crate draws transcription, and the reason to draw it as a
/// track rather than as another interval is that the span itself is the claim:
/// from the bent arrow to the hairpin is **one molecule**, so the genes under it
/// are co-transcribed and a promoter mutation upstream of the arrow changes all
/// of them at once. A [`FeatureTrack`](crate::FeatureTrack) stacked below draws
/// the genes; this draws only what a gene model cannot say.
///
/// The 5' leader is drawn hollow, which makes its length a distance a reader can
/// measure against the ruler rather than a number in a caption. A leaderless
/// transcript has no hollow segment at all and its arrowhead lands flush on the
/// start codon, so it is a different picture rather than a different label.
///
/// ```
/// use karyon::{Figure, Region, Strand, Terminator, TranscriptionUnit, TranscriptionUnitTrack};
///
/// let leaderless = TranscriptionUnit::new(1_000, 2_400, Strand::Forward)
///     .cds_start(1_000)
///     .terminator(Terminator::Intrinsic)
///     .name("esxB");
/// assert!(leaderless.is_leaderless());
///
/// let led = TranscriptionUnit::new(3_000, 4_200, Strand::Forward).cds_start(3_064);
/// assert_eq!(led.leader(), Some(64));
///
/// let svg = Figure::new(Region::new("NC_000962.3", 0, 5_000).unwrap())
///     .push(TranscriptionUnitTrack::new(vec![leaderless, led]).label("transcripts"))
///     .to_svg();
/// assert!(svg.contains("esxB"));
/// ```
#[derive(Debug, Clone)]
pub struct TranscriptionUnitTrack {
    units: Vec<TranscriptionUnit>,
    row_height: f64,
    row_gap: f64,
    label: Option<String>,
    show_names: bool,
    color: Option<String>,
}

impl TranscriptionUnitTrack {
    /// A track over these transcripts.
    pub fn new(units: Vec<TranscriptionUnit>) -> Self {
        TranscriptionUnitTrack {
            units,
            row_height: 30.0,
            row_gap: 6.0,
            label: None,
            show_names: true,
            color: None,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one row, which has to hold the bent arrow as well as
    /// the transcript.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(12.0);
        self
    }

    /// Sets the gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Whether to print transcript names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Overrides the colour, which is otherwise the crate's strand colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// The transcripts.
    pub fn units(&self) -> &[TranscriptionUnit] {
        &self.units
    }

    /// How many of them begin at their own start codon.
    pub fn leaderless(&self) -> usize {
        self.units
            .iter()
            .filter(|unit| unit.is_leaderless())
            .count()
    }

    /// Which row each transcript goes on, packing so none overlaps another.
    ///
    /// A transcript needs room for its bent arrow before the start site and for
    /// its terminator after the 3' end, so the packer reserves a margin in
    /// pixels rather than packing on the coordinates alone.
    fn rows(&self, scale: &Scale) -> Vec<usize> {
        let margin = self.row_height * 0.75;
        let mut ends: Vec<f64> = Vec::new();
        let mut rows = Vec::with_capacity(self.units.len());
        for unit in &self.units {
            let (low, high) = unit.extent();
            let left = scale.x(low) - margin;
            let right = scale.x(high) + margin;
            let row = ends.iter().position(|end| *end <= left);
            match row {
                Some(row) => {
                    ends[row] = right;
                    rows.push(row);
                }
                None => {
                    ends.push(right);
                    rows.push(ends.len() - 1);
                }
            }
        }
        rows
    }

    /// How many rows the packing needs.
    fn row_count(&self, scale: &Scale) -> usize {
        self.rows(scale)
            .into_iter()
            .max()
            .map(|row| row + 1)
            .unwrap_or(1)
    }
}

impl Track for TranscriptionUnitTrack {
    fn height(&self, scale: &Scale) -> f64 {
        let rows = self.row_count(scale) as f64;
        rows * self.row_height + (rows - 1.0) * self.row_gap
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        if self.units.is_empty() {
            return;
        }
        let rows = self.rows(ctx.scale);
        for (unit, row) in self.units.iter().zip(&rows) {
            let top = ctx.band.y + *row as f64 * (self.row_height + self.row_gap);
            self.draw_unit(ctx, unit, top);
        }
    }
}

impl TranscriptionUnitTrack {
    fn draw_unit(&self, ctx: &mut DrawContext<'_>, unit: &TranscriptionUnit, top: f64) {
        let band = ctx.band;
        let theme = ctx.theme;
        let ink = self
            .color
            .clone()
            .unwrap_or_else(|| strand_color(unit.strand, theme).to_string());

        // The transcript runs along the lower third and the bent arrow rises
        // out of it, so a row reads as a molecule with a start rather than as a
        // bar with a decoration.
        let wire_y = top + self.row_height * 0.72;
        let arm_y = top + self.row_height * 0.26;
        let reverse = unit.strand == Strand::Reverse;
        let direction = if reverse { -1.0 } else { 1.0 };

        let tss_x = ctx.scale.x_center(unit.tss);
        let end_x = ctx.scale.x_center(unit.end);
        if tss_x.min(end_x) > band.right() + 20.0 || tss_x.max(end_x) < band.x - 20.0 {
            return;
        }

        // The leader is hollow and the coding part solid, so the distance from
        // the start site to the start codon can be measured off the ruler.
        // Only when the start codon is downstream of the start site: a leader
        // that came back `None` is a contradiction, not a distance to draw.
        let cds_x = unit
            .leader()
            .and(unit.cds_start)
            .map(|at| ctx.scale.x_center(at));
        match cds_x {
            Some(cds_x) if (cds_x - tss_x).abs() > 0.6 => {
                ctx.svg.line(
                    tss_x,
                    wire_y,
                    cds_x,
                    wire_y,
                    &mix(&ink, theme.surface(), 0.55),
                    2.6,
                );
                ctx.svg.line(cds_x, wire_y, end_x, wire_y, &ink, 2.6);
                // A tick where translation starts, so the hollow run has an end
                // and not just a colour change.
                ctx.svg
                    .line(cds_x, wire_y - 4.0, cds_x, wire_y + 4.0, &ink, 1.2);
            }
            _ => ctx.svg.line(tss_x, wire_y, end_x, wire_y, &ink, 2.6),
        }

        // The bent arrow: a stem standing on the start site base, then an arm
        // downstream with a head on it.
        let arm = (self.row_height * 0.55).min((end_x - tss_x).abs().max(6.0));
        let tip = tss_x + direction * arm;
        ctx.svg.path_stroked(
            &format!(
                "M {} {} L {} {} L {} {}",
                num(tss_x),
                num(wire_y),
                num(tss_x),
                num(arm_y),
                num(tip),
                num(arm_y)
            ),
            &ink,
            1.6,
        );
        let head = 4.5f64.min(arm.max(1.0));
        ctx.svg.polygon(
            &[
                (tip + direction * head, arm_y),
                (tip - direction * head * 0.2, arm_y - head * 0.62),
                (tip - direction * head * 0.2, arm_y + head * 0.62),
            ],
            &ink,
        );

        self.draw_terminator(ctx, unit, end_x, wire_y, direction, &ink);

        // The name goes above the wire and clear of the arrow's arm, in the
        // empty stretch the transcript leaves between them. Below the wire it
        // would sit on the next row, and on a single-row band the clip would eat
        // it outright.
        if self.show_names {
            if let Some(name) = &unit.name {
                let size = theme.font_size;
                let width = text_width(name, size);
                let after_arm = tss_x + direction * (arm + 6.0);
                let low = after_arm.min(end_x).max(band.x);
                let high = after_arm.max(end_x).min(band.right());
                if high - low >= width + 4.0 {
                    ctx.svg.text(
                        (low + high) / 2.0,
                        wire_y - 4.0,
                        name,
                        &theme.muted,
                        size,
                        Anchor::Middle,
                    );
                }
            }
        }
    }

    /// A hairpin, a bar, or nothing at all when how it ends is not known.
    fn draw_terminator(
        &self,
        ctx: &mut DrawContext<'_>,
        unit: &TranscriptionUnit,
        x: f64,
        wire_y: f64,
        direction: f64,
        ink: &str,
    ) {
        let stem = (self.row_height * 0.3).max(5.0);
        match unit.terminator {
            Terminator::Intrinsic => {
                // Stem and loop: the RNA folding back on itself, which is what
                // makes this terminator need no protein.
                let radius = (stem * 0.42).max(2.0);
                ctx.svg.line(x, wire_y, x, wire_y - stem + radius, ink, 1.6);
                ctx.svg
                    .circle_ringed(x, wire_y - stem, radius, ctx.theme.surface(), ink, 1.6);
            }
            Terminator::RhoDependent => {
                // A bar across the transcript: something else stopped it.
                ctx.svg
                    .line(x, wire_y - stem * 0.6, x, wire_y + stem * 0.6, ink, 2.0);
                ctx.svg.line(
                    x - direction * 3.0,
                    wire_y - stem * 0.6,
                    x + direction * 3.0,
                    wire_y - stem * 0.6,
                    ink,
                    1.2,
                );
            }
            Terminator::Unknown => {
                ctx.svg
                    .line(x, wire_y - 3.0, x, wire_y + 3.0, &ctx.theme.muted, 1.0);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn scale(start: u64, end: u64) -> Scale {
        Scale::new(&Region::new("chr", start, end).unwrap(), 0.0, 900.0)
    }

    #[test]
    fn a_leader_is_the_distance_from_the_start_site_to_the_start_codon() {
        let unit = TranscriptionUnit::new(1_000, 2_000, Strand::Forward).cds_start(1_064);
        assert_eq!(unit.leader(), Some(64));
        assert!(!unit.is_leaderless());
    }

    #[test]
    fn a_leaderless_transcript_starts_at_its_own_start_codon() {
        let unit = TranscriptionUnit::new(1_000, 2_000, Strand::Forward).cds_start(1_000);
        assert_eq!(unit.leader(), Some(0));
        assert!(unit.is_leaderless());
    }

    #[test]
    fn a_reverse_leader_is_measured_the_other_way() {
        // Transcription starts at 2,000 and runs down to 1,000, so the start
        // codon is at a lower coordinate than the start site.
        let unit = TranscriptionUnit::new(2_000, 1_000, Strand::Reverse).cds_start(1_950);
        assert_eq!(unit.leader(), Some(50));
        assert_eq!(unit.extent(), (1_000, 2_000));
        assert_eq!(unit.len(), 1_000);
    }

    #[test]
    fn a_start_codon_upstream_of_the_start_site_is_a_contradiction() {
        // Not a negative leader: a start codon before transcription begins
        // cannot be translated from this transcript at all.
        let forward = TranscriptionUnit::new(1_000, 2_000, Strand::Forward).cds_start(900);
        assert_eq!(forward.leader(), None);
        let reverse = TranscriptionUnit::new(2_000, 1_000, Strand::Reverse).cds_start(2_100);
        assert_eq!(reverse.leader(), None);
    }

    #[test]
    fn no_start_codon_means_no_leader_rather_than_a_leader_of_zero() {
        let unit = TranscriptionUnit::new(1_000, 2_000, Strand::Forward);
        assert_eq!(unit.leader(), None);
        assert!(!unit.is_leaderless(), "unknown is not leaderless");
    }

    #[test]
    fn leaderless_transcripts_are_counted() {
        let track = TranscriptionUnitTrack::new(vec![
            TranscriptionUnit::new(0, 500, Strand::Forward).cds_start(0),
            TranscriptionUnit::new(1_000, 1_500, Strand::Forward).cds_start(1_040),
            TranscriptionUnit::new(2_000, 2_500, Strand::Forward).cds_start(2_000),
            TranscriptionUnit::new(3_000, 3_500, Strand::Forward),
        ]);
        assert_eq!(track.leaderless(), 2);
    }

    #[test]
    fn overlapping_transcripts_go_on_separate_rows() {
        let track = TranscriptionUnitTrack::new(vec![
            TranscriptionUnit::new(1_000, 5_000, Strand::Forward),
            TranscriptionUnit::new(2_000, 6_000, Strand::Forward),
            TranscriptionUnit::new(20_000, 24_000, Strand::Forward),
        ]);
        let rows = track.rows(&scale(0, 30_000));
        assert_eq!(rows[0], 0);
        assert_eq!(rows[1], 1, "it overlaps the first");
        assert_eq!(rows[2], 0, "and this one is clear of both");
        assert_eq!(track.row_count(&scale(0, 30_000)), 2);
    }

    #[test]
    fn packing_is_in_pixels_so_the_glyphs_do_not_collide_when_zoomed_out() {
        // Two transcripts a kilobase apart share a row at 5 kb across and not
        // at 3 Mb, where their arrows and hairpins would sit on top of one
        // another even though the coordinates never overlap.
        let track = TranscriptionUnitTrack::new(vec![
            TranscriptionUnit::new(1_000, 2_000, Strand::Forward),
            TranscriptionUnit::new(3_000, 4_000, Strand::Forward),
        ]);
        assert_eq!(track.row_count(&scale(0, 5_000)), 1);
        assert_eq!(track.row_count(&scale(0, 3_000_000)), 2);
    }

    #[test]
    fn height_follows_the_packing() {
        let track = TranscriptionUnitTrack::new(vec![
            TranscriptionUnit::new(1_000, 5_000, Strand::Forward),
            TranscriptionUnit::new(2_000, 6_000, Strand::Forward),
        ]);
        let one = TranscriptionUnitTrack::new(vec![TranscriptionUnit::new(
            1_000,
            5_000,
            Strand::Forward,
        )]);
        assert!(track.height(&scale(0, 30_000)) > one.height(&scale(0, 30_000)));
    }

    #[test]
    fn a_hairpin_and_a_rho_terminator_are_different_pictures() {
        let render = |terminator: Terminator| {
            Figure::new(Region::new("chr", 0, 5_000).unwrap())
                .push(TranscriptionUnitTrack::new(vec![TranscriptionUnit::new(
                    1_000,
                    4_000,
                    Strand::Forward,
                )
                .terminator(terminator)]))
                .to_svg()
        };
        let hairpin = render(Terminator::Intrinsic);
        assert!(hairpin.contains("<circle"), "the loop is a circle");
        assert_ne!(hairpin, render(Terminator::RhoDependent));
        assert_ne!(
            render(Terminator::RhoDependent),
            render(Terminator::Unknown)
        );
    }

    #[test]
    fn a_leaderless_transcript_looks_different_from_one_with_a_leader() {
        let render = |cds: u64| {
            Figure::new(Region::new("chr", 0, 5_000).unwrap())
                .push(TranscriptionUnitTrack::new(vec![TranscriptionUnit::new(
                    1_000,
                    4_000,
                    Strand::Forward,
                )
                .cds_start(cds)]))
                .to_svg()
        };
        let leaderless = render(1_000);
        let led = render(1_400);
        assert_ne!(leaderless, led);
        assert!(
            led.matches("<line").count() > leaderless.matches("<line").count(),
            "the hollow leader and its tick are extra strokes"
        );
    }

    #[test]
    fn a_reverse_transcript_points_the_other_way() {
        let render = |strand: Strand, tss: u64, end: u64| {
            Figure::new(Region::new("chr", 0, 5_000).unwrap())
                .push(TranscriptionUnitTrack::new(vec![TranscriptionUnit::new(
                    tss, end, strand,
                )]))
                .to_svg()
        };
        assert_ne!(
            render(Strand::Forward, 1_000, 4_000),
            render(Strand::Reverse, 4_000, 1_000)
        );
    }

    #[test]
    fn a_transcript_off_screen_draws_nothing_of_itself() {
        let away = Figure::new(Region::new("chr", 0, 5_000).unwrap())
            .push(TranscriptionUnitTrack::new(vec![TranscriptionUnit::new(
                900_000,
                904_000,
                Strand::Forward,
            )]))
            .to_svg();
        assert!(!away.contains("<polygon"), "no arrowhead is in view");
    }

    #[test]
    fn an_empty_track_still_has_a_height() {
        let track = TranscriptionUnitTrack::new(vec![]);
        assert!(track.height(&scale(0, 100)) > 0.0);
        assert_eq!(track.leaderless(), 0);
    }
}
