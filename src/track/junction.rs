//! Splice junctions as arcs, each one weighted by the reads that crossed it.
//!
//! An RNA-seq aligner reports, for every intron it saw a read step over, how
//! many reads made that step. The figure the field reads them in is the sashimi
//! plot: an arc per intron over the coverage profile, thicker where more reads
//! crossed, with the count printed on it. What a reader is comparing is two
//! arcs over the same exon, so the ratio between them is the whole point.
//!
//! # Why this is not a structural variant track
//!
//! [`StructuralTrack`](crate::StructuralTrack) also draws arcs, also weights
//! them by support, and its arcs use the same quadratic, which this module
//! calls rather than copies. Three things separate them and none is a setting.
//!
//! What the arc joins. A structural variant joins two breakpoints, which are
//! bases. An intron is not at a base, it is the boundary between two, so the
//! feet here are at the left edge of a base rather than at its middle. At
//! twenty pixels a base that is ten pixels of drift, and the arc stops meeting
//! the step in the coverage profile underneath it.
//!
//! How high it goes. A structural variant arcs by the distance between its
//! ends, because reaching further is what it did. An intron reaching further is
//! not a bigger event, so height here carries nothing at all: arcs are put in
//! lanes so they miss each other, and [`Track::y_axis_width`] answers with
//! nought so that nothing invites a reader to measure them.
//!
//! What is printed. A structural variant keeps its support in the tooltip and
//! says in its own module doc that stroke weight is to be read as an ordering
//! and never as a length. That is the right answer there and not enough here,
//! because the ratio between two junctions is the finding, so the count is
//! printed over the apex.
//!
//! # Why this is not a split read track either
//!
//! A spliced alignment is one primary record whose CIGAR steps over the intron,
//! and it carries no `SA` tag.
//! [`read::split::reads`](crate::read::split::reads) counts such a record as
//! not split and emits nothing: measured, three spliced records in, nought
//! molecules out. The two are different events. A supplementary alignment is
//! one molecule in several places; a skipped region is one alignment that
//! stepped over the reference.
//!
//! And a split read track is one row per molecule. Four hundred reads over one
//! exon are four hundred rows there, which is the raw observation it exists to
//! show. Here they are one arc labelled 400, and the collapsing is the point.
//!
//! # A junction no read crossed is not an observation
//!
//! [`Junction::new`] is the only way to make one and the fields are private, so
//! a junction cannot be assembled with no support behind it. One that arrives
//! with nought reads anyway is counted and not drawn, because an arc on the
//! page says a junction was seen.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::feature::span_label;
use crate::track::{arc_path, DrawContext, Track};
use crate::{Region, Strand};

/// Which dinucleotides the intron began and ended with.
///
/// The six codes an aligner writes fold to four, because the pairs differ only
/// in which strand the intron is on and the strand is its own field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Motif {
    /// The common one, about ninety-nine intron in a hundred.
    GtAg,
    /// The minor U2 variant.
    GcAg,
    /// The U12 spliceosome's.
    AtAc,
    /// None of the above, which is usually an alignment artefact and
    /// occasionally the interesting one.
    Noncanonical,
}

impl Motif {
    /// How it is spelled in a tooltip.
    pub fn name(self) -> &'static str {
        match self {
            Motif::GtAg => "GT/AG",
            Motif::GcAg => "GC/AG",
            Motif::AtAc => "AT/AC",
            Motif::Noncanonical => "non-canonical",
        }
    }

    /// Whether a spliceosome is known to use it.
    pub fn is_canonical(self) -> bool {
        !matches!(self, Motif::Noncanonical)
    }
}

/// One intron, and the reads that crossed it.
///
/// The fields are private and [`Junction::new`] is the only way in, so a
/// junction cannot be assembled without stating how many reads crossed it.
#[derive(Debug, Clone, PartialEq)]
pub struct Junction {
    start: u64,
    end: u64,
    reads: u32,
    multi: Option<u32>,
    overhang: Option<u32>,
    strand: Strand,
    motif: Option<Motif>,
    annotated: Option<bool>,
    name: Option<String>,
}

impl Junction {
    /// An intron from `start` to `end`, crossed by `reads` reads.
    ///
    /// Half-open: `start` is the first base of the intron and `end` is one past
    /// its last, so the arc's feet are the two exon edges.
    pub fn new(start: u64, end: u64, reads: u32) -> Self {
        Junction {
            start: start.min(end),
            end: start.max(end),
            reads,
            multi: None,
            overhang: None,
            strand: Strand::Unknown,
            motif: None,
            annotated: None,
            name: None,
        }
    }

    /// Sets the multi-mapping reads, which are never added to the unique ones.
    pub fn multi(mut self, reads: u32) -> Self {
        self.multi = Some(reads);
        self
    }

    /// Sets the longest spliced overhang among the supporting reads.
    pub fn overhang(mut self, bases: u32) -> Self {
        self.overhang = Some(bases);
        self
    }

    /// Sets which strand the intron is on.
    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Sets the dinucleotides the intron began and ended with.
    pub fn motif(mut self, motif: Motif) -> Self {
        self.motif = Some(motif);
        self
    }

    /// Says whether an annotation already held this junction.
    ///
    /// Left alone it stays unstated, and unstated is not false: a file with no
    /// annotation column says nothing about novelty, and drawing that as novel
    /// would be a discovery nobody made.
    pub fn annotated(mut self, annotated: bool) -> Self {
        self.annotated = Some(annotated);
        self
    }

    /// Names the junction.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// First base of the intron, 0-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// One past the last base of the intron.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// How many uniquely mapping reads crossed it.
    pub fn reads(&self) -> u32 {
        self.reads
    }

    /// How many multi-mapping reads crossed it, where the file counted them.
    pub fn multi_reads(&self) -> Option<u32> {
        self.multi
    }

    /// The longest spliced overhang, where the file reported one.
    pub fn overhang_bases(&self) -> Option<u32> {
        self.overhang
    }

    /// Which strand the intron is on.
    pub fn on_strand(&self) -> Strand {
        self.strand
    }

    /// The dinucleotides it began and ended with, where the file said.
    pub fn intron_motif(&self) -> Option<Motif> {
        self.motif
    }

    /// Whether an annotation held it, or `None` where nothing was stated.
    pub fn in_annotation(&self) -> Option<bool> {
        self.annotated
    }

    /// The name it was given.
    pub fn label(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// How many bases the intron covers.
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Whether it covers none.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether a read was actually seen to cross it.
    ///
    /// An arc on the page says a junction was observed, so one nobody crossed
    /// is counted rather than drawn.
    pub fn is_observed(&self) -> bool {
        self.reads > 0 && !self.is_empty()
    }
}

/// Splice junctions drawn as arcs over a shared coordinate axis.
///
/// ```
/// use karyon::{plot, Junction, JunctionTrack};
///
/// let junctions = vec![
///     Junction::new(1_000, 3_000, 412).motif(karyon::Motif::GtAg),
///     // The minor isoform: the ratio between the two is the finding.
///     Junction::new(1_000, 5_000, 9),
/// ];
///
/// let svg = plot("chr1:1-6,000")
///     .expect("a region")
///     .add_junctions(junctions)
///     .label("junctions")
///     .to_svg();
///
/// assert!(svg.contains("412"));
/// assert!(svg.contains("9 reads"));
/// ```
#[derive(Debug, Clone)]
pub struct JunctionTrack {
    junctions: Vec<Junction>,
    label: Option<String>,
    height: f64,
    min_stroke: f64,
    max_stroke: f64,
    saturating_reads: u32,
    min_reads: u32,
    color: Option<String>,
    show_counts: bool,
}

impl JunctionTrack {
    /// A track over `junctions`.
    pub fn new(junctions: impl Into<Vec<Junction>>) -> Self {
        JunctionTrack {
            junctions: junctions.into(),
            label: None,
            height: 64.0,
            min_stroke: 0.9,
            max_stroke: 5.0,
            saturating_reads: 200,
            min_reads: 1,
            color: None,
            show_counts: true,
        }
    }

    /// Sets the name in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, pixels: f64) -> Self {
        self.height = crate::svg::finite_within(pixels, 16.0, 400.0, 64.0);
        self
    }

    /// Sets the thinnest and thickest an arc is drawn.
    pub fn strokes(mut self, min: f64, max: f64) -> Self {
        if min.is_finite() && max.is_finite() && min > 0.0 && max >= min {
            self.min_stroke = min;
            self.max_stroke = max;
        }
        self
    }

    /// How many reads make an arc as thick as it gets.
    pub fn saturating_reads(mut self, reads: u32) -> Self {
        self.saturating_reads = reads.max(1);
        self
    }

    /// Drops junctions crossed by fewer than `reads` reads.
    ///
    /// A filter, not a recommendation, and how many it dropped is printed on
    /// the figure: a threshold nobody can see is worse than no threshold. One
    /// by default, which drops only the junctions nobody crossed.
    pub fn min_reads(mut self, reads: u32) -> Self {
        self.min_reads = reads.max(1);
        self
    }

    /// Sets the ink the arcs are drawn in.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Prints or hides the read count over each arc.
    pub fn show_counts(mut self, show: bool) -> Self {
        self.show_counts = show;
        self
    }

    /// The junctions, in the order they were given.
    pub fn junctions(&self) -> &[Junction] {
        &self.junctions
    }

    /// The junctions that clear the threshold and were actually crossed.
    pub fn drawn(&self) -> impl Iterator<Item = &Junction> {
        self.junctions
            .iter()
            .filter(|j| j.is_observed() && j.reads >= self.min_reads)
    }

    /// How many were held back, by the threshold or for never being crossed.
    pub fn discarded(&self) -> usize {
        self.junctions.len() - self.drawn().count()
    }

    /// The most reads any drawn junction carries.
    pub fn busiest(&self) -> Option<u32> {
        self.drawn().map(Junction::reads).max()
    }

    /// How thick an arc carrying `reads` reads is drawn.
    ///
    /// Logarithmic, because junction counts inside one gene span three or four
    /// orders of magnitude. On a linear ramp every minor isoform sits on the
    /// floor and becomes indistinguishable from every other minor isoform,
    /// which is the comparison the figure is made for.
    fn weight(&self, reads: u32) -> f64 {
        let top = (1.0 + self.saturating_reads as f64).ln();
        let here = (1.0 + reads as f64).ln();
        let ratio = if top > 0.0 { here / top } else { 1.0 };
        self.min_stroke + (self.max_stroke - self.min_stroke) * ratio.clamp(0.0, 1.0)
    }

    /// Which lane each drawn junction goes in, so arcs miss each other.
    ///
    /// Packed on the arcs alone and never on the width of their labels. The
    /// figure asks a track how tall it is before it hands it a theme, so a
    /// packing that measured text would put the arcs in one set of lanes and
    /// reserve the room for another. A label that does not fit its chord is
    /// left out instead, which is what a structural variant already does.
    ///
    /// Busiest first, ties by position, so the lane a junction lands in comes
    /// from its own numbers rather than from the order a file listed it.
    fn lanes(&self, scale: &Scale, region: &Region) -> Vec<(usize, &Junction)> {
        // Only what this window draws. Packing over the whole track let a
        // junction a megabase away take a lane, and the arcs that are on screen
        // were flattened to make room for it: measured, one visible arc went
        // from forty-nine pixels tall to two and a half when twenty off-screen
        // junctions were added to the same track.
        let mut order: Vec<&Junction> = self
            .drawn()
            .filter(|j| j.end > region.start() && j.start < region.end())
            .collect();
        order.sort_by(|a, b| {
            b.reads
                .cmp(&a.reads)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });

        let mut ends: Vec<f64> = Vec::new();
        let mut placed = Vec::with_capacity(order.len());
        for junction in order {
            let x0 = scale.x(junction.start);
            let x1 = scale.x(junction.end);
            let lane = ends
                .iter()
                .position(|taken| *taken <= x0)
                .unwrap_or(ends.len());
            if lane == ends.len() {
                ends.push(x1);
            } else {
                ends[lane] = x1;
            }
            placed.push((lane, junction));
        }
        placed
    }

    /// What a pointer hovering one arc is told.
    fn tooltip(&self, junction: &Junction) -> String {
        let mut said = format!(
            "{}, {} read{}",
            span_label(junction.start, junction.end),
            crate::track::axis::group_thousands(junction.reads as u64),
            if junction.reads == 1 { "" } else { "s" }
        );
        if let Some(multi) = junction.multi {
            said.push_str(&format!(
                ", {} multi-mapping",
                crate::track::axis::group_thousands(multi as u64)
            ));
        }
        if let Some(motif) = junction.motif {
            said.push_str(&format!(", {}", motif.name()));
        }
        match junction.annotated {
            Some(true) => said.push_str(", annotated"),
            Some(false) => said.push_str(", novel"),
            None => {}
        }
        if let Some(name) = &junction.name {
            said.push_str(&format!(", {name}"));
        }
        said
    }
}

impl Track for JunctionTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    /// Nought, and on purpose. The height of an arc carries no measurement, so
    /// a strip of ticks beside this band would invite a reader to measure one.
    fn y_axis_width(&self, _theme: &Theme) -> f64 {
        0.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let baseline = band.bottom();
        let ink = self
            .color
            .clone()
            .unwrap_or_else(|| mix(&ctx.theme.accent, &ctx.theme.foreground, 0.15));

        ctx.svg.line(
            band.x,
            baseline,
            band.right(),
            baseline,
            &ctx.theme.rule,
            ctx.theme.tokens.hairline,
        );

        let placed = self.lanes(ctx.scale, ctx.region);
        let lanes = placed.iter().map(|(lane, _)| *lane).max().unwrap_or(0) + 1;
        let size = ctx.theme.font_size - 1.0;
        let room = (band.h - size - 4.0).max(6.0);
        let step = room / lanes as f64;

        for (lane, junction) in &placed {
            // The left edge of a base, not its middle. A junction is not at a
            // base, it is the boundary between two, and the arc has to meet the
            // step in whatever is drawn under it.
            let x0 = ctx.scale.x(junction.start);
            let x1 = ctx.scale.x(junction.end);
            let apex = baseline - (*lane as f64 + 1.0) * step;
            let weight = self.weight(junction.reads);

            ctx.svg.begin_titled(&self.tooltip(junction));
            if (x1 - x0).abs() < 1.0 {
                // Nothing to arch over at this zoom, so a spike rather than a
                // curve that would be drawn as a vertical line anyway.
                ctx.svg.line(x0, baseline, x0, apex, &ink, weight);
            } else {
                ctx.svg
                    .path_stroked(&arc_path(x0, x1, baseline, apex), &ink, weight);
            }

            if self.show_counts {
                let text = crate::track::axis::group_thousands(junction.reads as u64);
                // Elided where the chord is too narrow, which is what keeps the
                // packing free of the theme: a label that does not fit is left
                // out rather than given a lane of its own.
                if (x1 - x0).abs() > text_width(&text, size) + 6.0 {
                    ctx.svg.text(
                        (x0 + x1) / 2.0,
                        apex - 3.0,
                        &text,
                        &ctx.theme.foreground,
                        size,
                        Anchor::Middle,
                    );
                }
            }
            ctx.svg.end_group();
        }

        let held = self.discarded();
        if held > 0 {
            let text = format!("+{held} not drawn");
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 3.0,
                &text,
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 0, 6_000).unwrap()
    }

    fn drawn(track: JunctionTrack) -> String {
        Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    #[test]
    fn a_junction_no_read_crossed_is_counted_and_never_drawn() {
        // An arc on the page says a junction was observed.
        let track = JunctionTrack::new(vec![
            Junction::new(1_000, 2_000, 40),
            Junction::new(1_000, 3_000, 0),
        ]);
        assert_eq!(track.drawn().count(), 1);
        assert_eq!(track.discarded(), 1);
        let svg = drawn(track);
        assert!(svg.contains("+1 not drawn"), "the reader was not told");
    }

    #[test]
    fn the_feet_are_the_edges_of_the_intron_and_not_the_middles_of_its_bases() {
        // A junction is not at a base, it is the boundary between two, and the
        // arc has to meet the step in whatever is drawn under it. An intron
        // starting at the first base of the window puts that to the test
        // without the test having to know the figure's geometry: the left edge
        // of base nought is exactly where the band starts, and the middle of it
        // is half a base further in.
        let svg = drawn(JunctionTrack::new(vec![Junction::new(0, 2_000, 40)]));
        let band_x = svg
            .split("<line x1=\"")
            .nth(1)
            .and_then(|piece| piece.split('"').next())
            .expect("a baseline");
        let foot = svg
            .split("d=\"M")
            .nth(1)
            .and_then(|piece| piece.split(' ').next())
            .expect("an arc");
        assert_eq!(foot, band_x, "the arc's foot is not the edge of its base");
    }

    #[test]
    fn thickness_is_logarithmic_so_the_quiet_isoforms_stay_apart() {
        // Counts inside one gene span three or four orders of magnitude. On a
        // linear ramp every minor isoform sits on the floor together.
        let track = JunctionTrack::new(Vec::new());
        let one = track.weight(1);
        let ten = track.weight(10);
        let hundred = track.weight(100);
        assert!(ten - one > 0.3, "1 and 10 are indistinguishable");
        // Each tenfold step is about as big as the last, which is the property
        // a linear ramp has not got.
        let first = ten - one;
        let second = hundred - ten;
        assert!(
            (first - second).abs() < first * 0.5,
            "the steps are not even: {first} then {second}"
        );
    }

    #[test]
    fn arcs_that_overlap_go_in_different_lanes_and_the_busiest_goes_lowest() {
        let scale = crate::scale::Scale::new(&region(), 0.0, 600.0);
        let track = JunctionTrack::new(vec![
            Junction::new(1_000, 5_000, 3),
            Junction::new(1_000, 2_000, 400),
            Junction::new(3_000, 4_000, 200),
        ]);
        let lanes = track.lanes(&scale, &region());
        let lane_of = |reads: u32| {
            lanes
                .iter()
                .find(|(_, j)| j.reads() == reads)
                .map(|(lane, _)| *lane)
                .expect("placed")
        };
        // 400 and 200 do not overlap, so they share the lowest lane.
        assert_eq!(lane_of(400), 0);
        assert_eq!(lane_of(200), 0);
        // The one spanning both is pushed up.
        assert_eq!(lane_of(3), 1);
    }

    #[test]
    fn the_lane_comes_from_the_numbers_and_not_from_the_order_of_the_file() {
        let scale = crate::scale::Scale::new(&region(), 0.0, 600.0);
        let a = Junction::new(1_000, 5_000, 3);
        let b = Junction::new(1_000, 2_000, 400);
        let one = JunctionTrack::new(vec![a.clone(), b.clone()]);
        let other = JunctionTrack::new(vec![b, a]);
        let lanes = |t: &JunctionTrack| {
            let mut got: Vec<(usize, u32)> = t
                .lanes(&scale, &region())
                .into_iter()
                .map(|(lane, j)| (lane, j.reads()))
                .collect();
            got.sort();
            got
        };
        assert_eq!(lanes(&one), lanes(&other));
    }

    #[test]
    fn a_junction_off_the_screen_does_not_flatten_the_ones_on_it() {
        // Packing over the whole track let data a megabase away take a lane,
        // and the arcs actually drawn were squashed to make room for it.
        let alone = drawn(JunctionTrack::new(vec![Junction::new(1_000, 3_000, 400)]));
        let mut crowded = vec![Junction::new(1_000, 3_000, 400)];
        for i in 0..20u64 {
            crowded.push(Junction::new(1_000_000, 1_010_000 + i, 500));
        }
        let with = drawn(JunctionTrack::new(crowded));

        let apex = |svg: &str| -> f64 {
            svg.split("Q")
                .nth(1)
                .and_then(|piece| piece.split_whitespace().nth(1))
                .and_then(|value| value.parse::<f64>().ok())
                .expect("an arc")
        };
        assert_eq!(
            apex(&alone),
            apex(&with),
            "the visible arc changed shape because of data a megabase away"
        );
    }

    #[test]
    fn novelty_left_unstated_is_not_novelty_stated() {
        let quiet = Junction::new(100, 200, 5);
        assert_eq!(quiet.in_annotation(), None);
        let svg = drawn(JunctionTrack::new(vec![quiet]));
        assert!(
            !svg.contains("novel"),
            "an unstated junction was called novel"
        );

        let stated = Junction::new(100, 200, 5).annotated(false);
        assert!(drawn(JunctionTrack::new(vec![stated])).contains("novel"));
    }

    #[test]
    fn multi_mapping_reads_are_never_folded_into_the_thickness() {
        let track = JunctionTrack::new(vec![Junction::new(100, 200, 4).multi(9_000)]);
        let bare = JunctionTrack::new(vec![Junction::new(100, 200, 4)]);
        assert_eq!(track.weight(4), bare.weight(4));
        assert!(drawn(track).contains("9,000 multi-mapping"));
    }

    #[test]
    fn no_value_axis_is_offered_because_the_height_measures_nothing() {
        let track = JunctionTrack::new(vec![Junction::new(100, 200, 5)]);
        assert_eq!(track.y_axis_width(&crate::Theme::light()), 0.0);
    }

    #[test]
    fn a_threshold_is_visible_or_it_is_not_a_threshold() {
        let track = JunctionTrack::new(vec![
            Junction::new(100, 200, 40),
            Junction::new(300, 400, 2),
            Junction::new(500, 600, 1),
        ])
        .min_reads(5);
        assert_eq!(track.drawn().count(), 1);
        assert!(drawn(track).contains("+2 not drawn"));
    }
}
