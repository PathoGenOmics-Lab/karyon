//! Structural variants as arcs between their two breakpoints.
//!
//! A [`StructuralVariant`] is two 0-based positions and a [`SvKind`] saying what
//! happened between them. Everything the band draws follows from that pair:
//! where the arc springs from, how far it reaches, and whether there is any
//! reference sequence underneath it to mark at all.
//!
//! # Why an arc
//!
//! A structural variant is not a mark at a position, it is a statement that two
//! positions belong together: the two ends of a deletion, the source and the
//! destination of a duplication, the pair of junctions an inversion creates. A
//! variant track draws one point per call and so cannot say that, and a bar
//! spanning the event says only that something happened in the middle, which is
//! usually the one place nothing happened.
//!
//! So each call is drawn as an arc from one breakpoint to the other, springing
//! from the axis at both ends. How high it arches follows how far apart the two
//! ends are, so a local event sits low and a rearrangement across a chromosome
//! reaches over everything between, which is what it actually did.
//!
//! # What the height of an arc means
//!
//! Where the call sits among the calls it was given, and nothing more. The
//! tallest arc is always the widest event in the track, on screen or not, and
//! every other height is a fraction of that, so two panels drawn from two sets
//! of calls do not share a scale. Height is not proportional to the span even
//! within one track, so read it as an ordering and never as a length. Weight is
//! a separate reading: an arc is drawn heavier the more reads support the call.
//!
//! # What to put under it
//!
//! A [`CoverageTrack`](crate::CoverageTrack). Half of reading an SV call is
//! whether the depth agrees: a deletion with no drop under it and a duplication
//! with no rise are calls to argue with, and the arc and the profile have to be
//! on one axis for anyone to see that.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::feature::span_label;
use crate::track::{DrawContext, Track};

/// What a structural variant did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvKind {
    /// Sequence that is in the reference and not in the sample.
    Deletion,
    /// Sequence present more times in the sample than in the reference.
    Duplication,
    /// Sequence that is the other way round.
    Inversion,
    /// Two places joined that are not neighbours in the reference. The two
    /// breakpoints may be on different sequences, which on a
    /// [`Genome`](crate::Genome) axis is two points far apart.
    Translocation,
    /// Sequence in the sample that is not in the reference, which has one
    /// breakpoint rather than two.
    Insertion,
}

impl SvKind {
    /// The name this class goes by in a caption.
    pub fn name(self) -> &'static str {
        match self {
            SvKind::Deletion => "deletion",
            SvKind::Duplication => "duplication",
            SvKind::Inversion => "inversion",
            SvKind::Translocation => "translocation",
            SvKind::Insertion => "insertion",
        }
    }

    /// Whether the call covers a stretch of the reference.
    ///
    /// An insertion does not: it is sequence that is not in the reference at
    /// all, so it has one breakpoint and no footprint. A translocation joins
    /// two places rather than covering what is between them.
    pub fn has_footprint(self) -> bool {
        matches!(
            self,
            SvKind::Deletion | SvKind::Duplication | SvKind::Inversion
        )
    }
}

/// One structural variant call.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralVariant {
    /// First breakpoint, 0-based.
    pub start: u64,
    /// Second breakpoint. Equal to the start for an insertion.
    pub end: u64,
    /// What the call says happened.
    pub kind: SvKind,
    /// How many reads support it, which sets how heavy the arc is drawn.
    pub support: Option<u32>,
    /// A label, drawn over the apex when there is room.
    pub name: Option<String>,
}

impl StructuralVariant {
    /// A call joining `start` and `end`.
    pub fn new(start: u64, end: u64, kind: SvKind) -> Self {
        StructuralVariant {
            start: start.min(end),
            end: start.max(end),
            kind,
            support: None,
            name: None,
        }
    }

    /// Sets how many reads support the call.
    pub fn support(mut self, reads: u32) -> Self {
        self.support = Some(reads);
        self
    }

    /// Sets the label.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// How far apart the two breakpoints are.
    pub fn span(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the call reaches into a span.
    ///
    /// Either breakpoint inside it, or the two of them either side: an arc that
    /// leaves the view at one end still says something, and one that leaves it
    /// at both is the widest event there is, so asking only whether a
    /// breakpoint is on screen would hide exactly the calls that matter most.
    /// The end is inclusive, so a call with one breakpoint (an insertion, where
    /// the start and the end are the same) on the first base of the span is
    /// kept.
    pub fn touches(&self, start: u64, end: u64) -> bool {
        self.start < end && self.end >= start
    }
}

/// Structural variant calls, drawn as arcs.
///
/// ```
/// use karyon::{Figure, Region, StructuralTrack, StructuralVariant, SvKind};
///
/// let calls = vec![
///     StructuralVariant::new(1_200, 3_400, SvKind::Deletion).support(31),
///     StructuralVariant::new(5_000, 5_900, SvKind::Duplication).support(12),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 0, 8_000).unwrap())
///     .push(StructuralTrack::new(calls).label("SV"))
///     .to_svg();
/// assert!(svg.contains("<path"));
/// ```
#[derive(Debug, Clone)]
pub struct StructuralTrack {
    variants: Vec<StructuralVariant>,
    label: Option<String>,
    height: f64,
    colors: Vec<(SvKind, String)>,
    min_stroke: f64,
    max_stroke: f64,
    saturating_support: u32,
    show_footprints: bool,
    show_names: bool,
}

impl StructuralTrack {
    /// A track over `variants`.
    pub fn new(variants: impl Into<Vec<StructuralVariant>>) -> Self {
        StructuralTrack {
            variants: variants.into(),
            label: None,
            height: 90.0,
            colors: Vec::new(),
            min_stroke: 1.0,
            max_stroke: 3.4,
            saturating_support: 30,
            show_footprints: true,
            show_names: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    ///
    /// The tallest arc uses all of it, so this is also how much room the
    /// longest event gets to arch over the shorter ones.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(20.0);
        self
    }

    /// Sets the colour of one kind of call.
    pub fn color(mut self, kind: SvKind, color: impl Into<String>) -> Self {
        self.colors.retain(|(existing, _)| *existing != kind);
        self.colors.push((kind, color.into()));
        self
    }

    /// Sets how heavy an arc is at no support and at full support.
    pub fn stroke(mut self, min: f64, max: f64) -> Self {
        self.min_stroke = min.max(0.2);
        self.max_stroke = max.max(self.min_stroke);
        self
    }

    /// Sets the read count at which an arc is drawn at full weight.
    pub fn saturating_support(mut self, reads: u32) -> Self {
        self.saturating_support = reads.max(1);
        self
    }

    /// Draws or hides the bar along the axis under events that cover one.
    pub fn show_footprints(mut self, show: bool) -> Self {
        self.show_footprints = show;
        self
    }

    /// Draws or hides names over the arcs.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// The calls.
    pub fn variants(&self) -> &[StructuralVariant] {
        &self.variants
    }

    /// Calls of one kind.
    pub fn of_kind(&self, kind: SvKind) -> impl Iterator<Item = &StructuralVariant> {
        self.variants.iter().filter(move |call| call.kind == kind)
    }

    /// The widest call, which is the one whose arc reaches the top.
    pub fn widest(&self) -> u64 {
        self.variants
            .iter()
            .map(|call| call.span())
            .max()
            .unwrap_or(1)
            .max(1)
    }

    /// The colour of one kind of call.
    pub fn color_of(&self, kind: SvKind, theme: &Theme) -> String {
        if let Some((_, color)) = self.colors.iter().find(|(existing, _)| *existing == kind) {
            return color.clone();
        }
        let index = match kind {
            SvKind::Deletion => 0,
            SvKind::Duplication => 1,
            SvKind::Inversion => 2,
            SvKind::Translocation => 3,
            SvKind::Insertion => 4,
        };
        theme.color(index).to_string()
    }

    /// How heavy the arc of a call is drawn.
    fn stroke_of(&self, call: &StructuralVariant) -> f64 {
        let ratio = call
            .support
            .map(|reads| reads as f64 / self.saturating_support as f64)
            .unwrap_or(0.5)
            .clamp(0.0, 1.0);
        self.min_stroke + (self.max_stroke - self.min_stroke) * ratio
    }

    /// How high the arc of a call reaches, as a fraction of the band.
    ///
    /// Follows the square root of the span rather than the span itself. A
    /// collection of calls is mostly small ones with a few enormous ones, and
    /// on a linear scale the small ones all lie flat against the axis in a
    /// single indistinguishable line.
    pub fn arch(&self, call: &StructuralVariant) -> f64 {
        arch_against(call, self.widest() as f64)
    }
}

/// The height of one call's arc against a given widest span.
///
/// Split out from [`StructuralTrack::arch`] because the widest span is a
/// property of the whole track and the drawing needs it once, not once per
/// call: `widest` walks every call, so working it out inside the loop made
/// drawing quadratic. Measured before this, at a fixed window and width:
/// 25,000 calls took 0.34 seconds and 50,000 took 1.36, which is four times the
/// work for twice the data.
fn arch_against(call: &StructuralVariant, widest: f64) -> f64 {
    ((call.span().max(1) as f64) / widest)
        .sqrt()
        .clamp(0.08, 1.0)
}

impl Track for StructuralTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let baseline = band.bottom();
        ctx.svg.line(
            band.x,
            baseline - 0.5,
            band.right(),
            baseline - 0.5,
            &ctx.theme.rule,
            1.0,
        );

        // Widest first, so a small event inside a large one is drawn over it
        // rather than under. The small one is the one being looked at.
        let mut order: Vec<usize> = (0..self.variants.len()).collect();
        order.sort_by_key(|index| std::cmp::Reverse(self.variants[*index].span()));

        // Once, rather than once per call: this walks every call in the track.
        let widest = self.widest() as f64;

        for index in order {
            let call = &self.variants[index];
            if !call.touches(ctx.region.start(), ctx.region.end()) {
                continue;
            }
            let color = self.color_of(call.kind, ctx.theme);
            let x0 = ctx.scale.x_center(call.start);
            let x1 = ctx.scale.x_center(call.end);
            let apex = baseline - arch_against(call, widest) * band.h;

            // A call is a footprint, an arc and sometimes a name, and the whole
            // of it is one event. The group holds them together so a pointer
            // anywhere on the call answers for the call.
            ctx.svg.begin_titled(&tooltip(call));

            if self.show_footprints && call.kind.has_footprint() {
                // What the call covers, along the axis. It is a footprint and
                // not the mark itself: the interesting part of a deletion is
                // its two ends, and a bar says only that something happened in
                // the middle.
                ctx.svg.rect(
                    x0,
                    baseline - 3.5,
                    (x1 - x0).max(0.8),
                    3.0,
                    &mix(ctx.theme.surface(), &color, 0.35),
                );
            }

            if call.kind == SvKind::Insertion || (x1 - x0).abs() < 1.0 {
                // One breakpoint and nothing to arch over: a spike.
                ctx.svg
                    .line(x0, baseline, x0, apex, &color, self.stroke_of(call));
                ctx.svg.circle(x0, apex, self.stroke_of(call) * 0.9, &color);
            } else {
                // Springing from the axis at both ends. A control point twice
                // the height puts the apex of the curve at the height asked for.
                let d = crate::track::arc_path(x0, x1, baseline, apex);
                ctx.svg.path_stroked(&d, &color, self.stroke_of(call));

                if self.show_names {
                    if let Some(name) = &call.name {
                        let size = ctx.theme.font_size - 2.0;
                        if (x1 - x0).abs() > text_width(name, size) + 6.0 {
                            ctx.svg.text(
                                (x0 + x1) / 2.0,
                                apex - 3.0,
                                name,
                                &ctx.theme.muted,
                                size,
                                Anchor::Middle,
                            );
                        }
                    }
                }
            }

            ctx.svg.end_group();
        }
    }
}

/// What a reader hovering one call is told.
///
/// What happened, where it happened, and how much read evidence there was for
/// saying so. Support is the one number on the glyph that cannot be read off
/// it: stroke weight orders the calls but does not measure them, and the
/// height of an arc is an ordering too, so the count belongs here.
///
/// A call whose breakpoints coincide gets one coordinate rather than two. An
/// insertion is sequence the reference does not have, so `4,362,001 to
/// 4,362,001` would be a span where the data has a point.
fn tooltip(call: &StructuralVariant) -> String {
    let mut text = call.kind.name().to_string();
    text.push_str(", ");
    if call.start == call.end {
        text.push_str(&group_thousands(call.start.saturating_add(1)));
    } else {
        text.push_str(&span_label(call.start, call.end));
    }
    if let Some(reads) = call.support {
        text.push_str(", ");
        text.push_str(&group_thousands(reads as u64));
        text.push_str(if reads == 1 {
            " supporting read"
        } else {
            " supporting reads"
        });
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 0, 10_000).unwrap()
    }

    fn calls() -> Vec<StructuralVariant> {
        vec![
            StructuralVariant::new(1_000, 2_000, SvKind::Deletion).support(30),
            StructuralVariant::new(3_000, 8_000, SvKind::Duplication).support(6),
            StructuralVariant::new(4_000, 4_500, SvKind::Inversion),
            StructuralVariant::new(9_000, 9_000, SvKind::Insertion).support(14),
        ]
    }

    #[test]
    fn breakpoints_come_back_the_right_way_round() {
        let backwards = StructuralVariant::new(900, 100, SvKind::Deletion);
        assert_eq!((backwards.start, backwards.end), (100, 900));
        assert_eq!(backwards.span(), 800);
    }

    #[test]
    fn a_call_with_one_end_on_screen_is_still_drawn() {
        // An arc reaching out of the view says something, and dropping it would
        // quietly hide every event that leaves the window.
        let reaching = StructuralVariant::new(5_000, 900_000, SvKind::Translocation);
        assert!(reaching.touches(0, 10_000));
        let elsewhere = StructuralVariant::new(500_000, 900_000, SvKind::Deletion);
        assert!(!elsewhere.touches(0, 10_000));
    }

    #[test]
    fn a_call_that_engulfs_the_view_is_the_one_least_safe_to_drop() {
        // Neither breakpoint of a 499 kb deletion at 1,000..500,000 is inside
        // chr1:100,000-200,000, and the window is entirely deleted, so an empty
        // band would be the most misleading answer there is.
        let engulfing = StructuralVariant::new(1_000, 500_000, SvKind::Deletion).support(40);
        assert!(engulfing.touches(100_000, 200_000));
        let svg = Figure::new(Region::new("chr1", 100_000, 200_000).unwrap())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![engulfing]).label("SV"))
            .to_svg();
        assert_eq!(svg.matches("<path").count(), 1, "the arc crosses the band");
    }

    #[test]
    fn only_some_kinds_cover_a_stretch_of_the_reference() {
        assert!(SvKind::Deletion.has_footprint());
        assert!(SvKind::Duplication.has_footprint());
        assert!(SvKind::Inversion.has_footprint());
        // An insertion is sequence the reference does not have, and a
        // translocation joins two places rather than covering between them.
        assert!(!SvKind::Insertion.has_footprint());
        assert!(!SvKind::Translocation.has_footprint());
    }

    #[test]
    fn the_arch_follows_the_span_but_not_linearly() {
        let track = StructuralTrack::new(calls());
        assert_eq!(track.widest(), 5_000);
        let small = StructuralVariant::new(0, 500, SvKind::Deletion);
        let large = StructuralVariant::new(0, 5_000, SvKind::Deletion);
        // The drawing works the widest span out once and passes it in, so the
        // two ways of asking have to agree or the split silently changed the
        // picture.
        let widest = track.widest() as f64;
        for call in &track.variants {
            assert!(
                (track.arch(call) - arch_against(call, widest)).abs() < f64::EPSILON,
                "the hoisted arch disagrees with the method"
            );
        }
        assert!(track.arch(&large) > track.arch(&small));
        assert!(
            (track.arch(&large) - 1.0).abs() < 1e-9,
            "the widest reaches the top"
        );
        // A tenth of the span is a third of the height, not a tenth: on a
        // linear scale every small call lies flat against the axis in one
        // indistinguishable line.
        assert!(track.arch(&small) > 0.25);
    }

    #[test]
    fn a_call_of_no_width_still_arches_off_the_axis() {
        let track = StructuralTrack::new(calls());
        let point = StructuralVariant::new(500, 500, SvKind::Insertion);
        assert!(
            track.arch(&point) >= 0.08,
            "otherwise it is the axis itself"
        );
    }

    #[test]
    fn support_decides_how_heavy_an_arc_is() {
        let track = StructuralTrack::new(calls()).saturating_support(30);
        let strong = StructuralVariant::new(0, 100, SvKind::Deletion).support(30);
        let weak = StructuralVariant::new(0, 100, SvKind::Deletion).support(3);
        assert!(track.stroke_of(&strong) > track.stroke_of(&weak));
        assert_eq!(track.stroke_of(&strong), 3.4, "saturated");
        // A call with no support figure is drawn in the middle rather than at
        // either end, since neither would be honest.
        let unknown = StructuralVariant::new(0, 100, SvKind::Deletion);
        assert!(track.stroke_of(&unknown) > track.stroke_of(&weak));
        assert!(track.stroke_of(&unknown) < track.stroke_of(&strong));
    }

    #[test]
    fn every_kind_gets_a_colour_of_its_own() {
        let theme = Theme::light();
        let track = StructuralTrack::new(calls());
        let mut seen = Vec::new();
        for kind in [
            SvKind::Deletion,
            SvKind::Duplication,
            SvKind::Inversion,
            SvKind::Translocation,
            SvKind::Insertion,
        ] {
            let color = track.color_of(kind, &theme);
            assert!(!seen.contains(&color), "{} repeats a colour", kind.name());
            seen.push(color);
        }
        // And one can be overridden without disturbing the rest.
        let custom = StructuralTrack::new(calls()).color(SvKind::Deletion, "#123456");
        assert_eq!(custom.color_of(SvKind::Deletion, &theme), "#123456");
        assert_eq!(
            custom.color_of(SvKind::Inversion, &theme),
            track.color_of(SvKind::Inversion, &theme)
        );
    }

    #[test]
    fn an_arc_springs_from_the_axis_at_both_ends() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(calls()).show_footprints(false))
            .to_svg();
        // Three arcs; the insertion is a spike rather than a curve.
        assert_eq!(svg.matches("<path").count(), 3);
        assert_eq!(svg.matches('Q').count(), 3, "one control point each");
    }

    #[test]
    fn an_insertion_is_a_spike_rather_than_an_arc() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                5_000,
                5_000,
                SvKind::Insertion,
            )]))
            .to_svg();
        assert!(!svg.contains("<path"), "nothing to arch over");
        assert!(svg.contains("<line"));
        assert!(svg.contains("<circle"));
    }

    #[test]
    fn a_footprint_is_drawn_only_where_there_is_one() {
        let with = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(calls()))
            .to_svg();
        let without = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(calls()).show_footprints(false))
            .to_svg();
        // Three of the four calls cover a stretch of reference.
        assert_eq!(
            with.matches("<rect").count() - without.matches("<rect").count(),
            3
        );
    }

    #[test]
    fn a_small_call_inside_a_large_one_is_drawn_over_it() {
        let nested = vec![
            StructuralVariant::new(0, 9_000, SvKind::Duplication),
            StructuralVariant::new(4_000, 4_100, SvKind::Deletion),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(nested).show_footprints(false))
            .to_svg();
        let theme = Theme::light();
        let track = StructuralTrack::new(Vec::new());
        let big = svg
            .find(&track.color_of(SvKind::Duplication, &theme))
            .unwrap();
        let small = svg.find(&track.color_of(SvKind::Deletion, &theme)).unwrap();
        assert!(big < small, "the small call is the one being looked at");
    }

    #[test]
    fn a_name_is_written_only_where_the_arc_is_wide_enough() {
        let roomy = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                1_000,
                9_000,
                SvKind::Deletion,
            )
            .name("RD1")]))
            .to_svg();
        assert!(roomy.contains(">RD1</text>"));

        let cramped = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                1_000,
                1_020,
                SvKind::Deletion,
            )
            .name("a_long_event_name")]))
            .to_svg();
        assert!(!cramped.contains(">a_long_event_name</text>"));
    }

    /// Every group a track opens has to be closed by exactly one `end_group`.
    /// This track is the one at risk of getting it wrong, because a spike used
    /// to leave the loop early.
    fn groups_balance(svg: &str) -> bool {
        svg.matches("<g ").count() + svg.matches("<g>").count() == svg.matches("</g>").count()
    }

    #[test]
    fn a_call_says_what_happened_where_and_on_what_evidence() {
        let svg = Figure::new(Region::new("chr1", 4_340_000, 4_400_000).unwrap())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                4_350_000,
                4_359_600,
                SvKind::Deletion,
            )
            .support(62)
            .name("RD1")]))
            .to_svg();
        assert!(
            svg.contains("<title>deletion, 4,350,001 to 4,359,600, 62 supporting reads</title>"),
            "{svg}"
        );
        assert!(groups_balance(&svg));
    }

    #[test]
    fn a_call_with_one_breakpoint_gets_one_coordinate() {
        // An insertion is sequence the reference does not have, so a span
        // would be a claim about a stretch that is not there.
        let svg = Figure::new(Region::new("chr1", 4_340_000, 4_400_000).unwrap())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                4_362_000,
                4_362_000,
                SvKind::Insertion,
            )
            .support(17)]))
            .to_svg();
        assert!(
            svg.contains("<title>insertion, 4,362,001, 17 supporting reads</title>"),
            "{svg}"
        );
        // A spike draws a line and a head and leaves the loop by another path
        // than an arc does, and it still closes its group.
        assert!(groups_balance(&svg));
    }

    #[test]
    fn a_call_with_no_support_figure_claims_none() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                4_000,
                4_500,
                SvKind::Inversion,
            )]))
            .to_svg();
        assert!(
            svg.contains("<title>inversion, 4,001 to 4,500</title>"),
            "{svg}"
        );
    }

    #[test]
    fn one_supporting_read_is_a_read_and_not_reads() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                1_000,
                2_000,
                SvKind::Duplication,
            )
            .support(1)]))
            .to_svg();
        assert!(
            svg.contains("<title>duplication, 1,001 to 2,000, 1 supporting read</title>"),
            "{svg}"
        );
    }

    #[test]
    fn every_call_on_screen_is_named_exactly_once() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(calls()))
            .to_svg();
        // Four calls, four tooltips: the footprint, the arc and the name of one
        // call are one event and share one group.
        assert_eq!(svg.matches("<title>").count(), 4, "{svg}");
        assert!(groups_balance(&svg));
    }

    #[test]
    fn an_empty_track_still_draws_its_axis() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(Vec::new()).label("SV"))
            .to_svg();
        assert!(svg.contains("<line"));
        assert!(!svg.contains("NaN"));
    }
}
