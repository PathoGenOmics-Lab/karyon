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
//! Where the call sits among the calls in view, and nothing more. The widest
//! call lying wholly inside the view reaches the top of the band, less a
//! little room for its stroke and for a name written over it, and every other
//! height is a fraction of that; a call wider than it, which can only be one
//! that leaves the view, reaches the top as well. When no call lies wholly
//! inside, because every call in view crosses an edge or the view is in the
//! middle of one, the widest call reaching into the view sets the scale,
//! counted as no wider than the view.
//!
//! It used to be the widest call in the track, on screen or not, and one
//! translocation off to the side was enough to press every arc in view flat
//! against the axis: a 1.2 Mb event left a 10 kb deletion at a tenth of the
//! band, which was most of the band spent on nothing. The price of using the
//! height is that panning or zooming rescales the arcs, so two views do not
//! share a scale any more than two panels do. Height is not proportional to the
//! span even within one view, so read it as an ordering and never as a length.
//! Weight is a separate reading: an arc is drawn heavier the more reads support
//! the call.
//!
//! # What to put under it
//!
//! A [`CoverageTrack`](crate::CoverageTrack). Half of reading an SV call is
//! whether the depth agrees: a deletion with no drop under it and a duplication
//! with no rise are calls to argue with, and the arc and the profile have to be
//! on one axis for anyone to see that.

use crate::region::Region;
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
    /// The tallest arc in view uses all of it but the line its name needs, so
    /// this is also how much room the widest event gets to arch over the
    /// shorter ones.
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

    /// The widest call in the track, which is the one whose arc reaches the
    /// top in a view that holds every call.
    pub fn widest(&self) -> u64 {
        self.variants
            .iter()
            .map(|call| call.span())
            .max()
            .unwrap_or(1)
            .max(1)
    }

    /// The span the arcs in `region` are measured against.
    ///
    /// The widest call with both breakpoints inside the region and some width
    /// between them. Both ends are exclusive, the call's as the reader makes it
    /// and the region's, so a call zoomed to exactly its own extent is inside.
    /// An insertion has no width, and letting one set the scale would send
    /// every arc beside it to the top.
    ///
    /// When nothing qualifies, because every call in view crosses an edge of it
    /// or the view sits inside one, it is the widest call reaching into the
    /// region, counted as no wider than the region. Uncapped, a translocation
    /// running off one side would set the scale for a deletion running off the
    /// other, and press it flat again; capped, a call that fills the view
    /// reaches the top and the rest are measured against the view.
    pub fn widest_in(&self, region: &Region) -> u64 {
        let (start, end) = (region.start(), region.end());
        let inside = self
            .variants
            .iter()
            .filter(|call| call.start >= start && call.end <= end && call.span() > 0)
            .map(StructuralVariant::span)
            .max();
        inside
            .or_else(|| {
                self.variants
                    .iter()
                    .filter(|call| call.touches(start, end))
                    .map(|call| call.span().min(end - start))
                    .max()
            })
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

    /// How high the arc of a call reaches in a view that holds every call, as
    /// a fraction of the room the band has for arcs: the band less what
    /// [`StructuralTrack::height`] keeps free at its top.
    ///
    /// Follows the square root of the span rather than the span itself. A
    /// collection of calls is mostly small ones with a few enormous ones, and
    /// on a linear scale the small ones all lie flat against the axis in a
    /// single indistinguishable line.
    pub fn arch(&self, call: &StructuralVariant) -> f64 {
        arch_against(call, self.widest() as f64)
    }

    /// How high the arc of a call reaches in `region`, as a fraction of the
    /// room the band has for arcs. This is the height it is drawn at.
    pub fn arch_in(&self, call: &StructuralVariant, region: &Region) -> f64 {
        arch_against(call, self.widest_in(region) as f64)
    }

    /// Room kept free at the top of the band, so the tallest arc is not cut.
    ///
    /// Always the heaviest stroke and its apex marker. A line of type as well,
    /// but only for a name that will be drawn and would otherwise run into the
    /// top: a name on an insertion, on an arc too narrow for it, centred off
    /// the plot or over an arc that peaks well below the top takes no room, and
    /// reserving it anyway lowered every arc in view by a tenth of the band for
    /// a line nobody saw.
    fn headroom(&self, ctx: &DrawContext<'_>, widest: f64) -> f64 {
        let stroke = self.max_stroke + 1.0;
        if !self.show_names {
            return stroke;
        }
        let band = ctx.band;
        let size = ctx.theme.font_size - 2.0;
        let line = size + 4.0;
        let reach = (band.h - stroke).max(band.h * 0.5);
        let crowded = self.variants.iter().any(|call| {
            let Some(name) = &call.name else {
                return false;
            };
            if call.kind == SvKind::Insertion || !call.touches(ctx.region.start(), ctx.region.end())
            {
                return false;
            }
            let (x0, x1) = (ctx.scale.x_center(call.start), ctx.scale.x_center(call.end));
            let middle = (x0 + x1) / 2.0;
            let drawn = (x1 - x0).abs() > text_width(name, size) + 6.0
                && middle >= band.x
                && middle <= band.right();
            let apex = band.bottom() - arch_against(call, widest) * reach;
            drawn && apex - line < band.y
        });
        // Crowded means the stroke's room left the name short of the top, so
        // the line is the larger of the two whenever it is taken.
        if crowded {
            line
        } else {
            stroke
        }
    }
}

/// The height of one call's arc against a given widest span.
///
/// Shared by [`StructuralTrack::arch`] and [`StructuralTrack::arch_in`]. The
/// widest span depends on the track and the view but not on the call, so the
/// drawing works it out once per draw and not once per call: `widest_in` walks
/// every call, so working it out inside the loop made drawing quadratic. Measured before this, at a fixed window and width:
/// 25,000 calls took 0.34 seconds and 50,000 took 1.36, which is four times the
/// work for twice the data.
fn arch_against(call: &StructuralVariant, widest: f64) -> f64 {
    ((call.span().max(1) as f64) / widest)
        .sqrt()
        .clamp(0.08, 1.0)
}

impl Track for StructuralTrack {
    fn noun(&self) -> &str {
        "structural variants"
    }

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
        let widest = self.widest_in(ctx.region) as f64;
        // What is left for arcs once the tallest one's name has room, and
        // never less than half the band, however small the band is set.
        let reach = (band.h - self.headroom(ctx, widest)).max(band.h * 0.5);

        for index in order {
            let call = &self.variants[index];
            if !call.touches(ctx.region.start(), ctx.region.end()) {
                continue;
            }
            let color = self.color_of(call.kind, ctx.theme);
            let x0 = ctx.scale.x_center(call.start);
            let x1 = ctx.scale.x_center(call.end);
            let apex = baseline - arch_against(call, widest) * reach;

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
                        // Written only where the whole of it fits: over an
                        // arc wide enough, and under the top of the band,
                        // which a band set near its minimum height does not
                        // leave room for over its tallest arc.
                        let size = ctx.theme.font_size - 2.0;
                        let fits_above = apex - 3.0 - size * 0.75 >= band.y;
                        if (x1 - x0).abs() > text_width(name, size) + 6.0 && fits_above {
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
        let widest = track.widest_in(&region()) as f64;
        for call in &track.variants {
            assert!(
                (track.arch_in(call, &region()) - arch_against(call, widest)).abs() < f64::EPSILON,
                "the hoisted arch disagrees with the method"
            );
        }
        // A view holding every call is the whole track.
        assert_eq!(track.widest_in(&region()), track.widest());
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

    /// Where the arc stroked in `color` peaks, from its control point: a
    /// quadratic reaches half way from its ends to its control.
    fn apex_of(svg: &str, color: &str) -> f64 {
        let path = svg
            .split("<path d=\"M")
            .skip(1)
            .find(|path| path.contains(&format!(r#"stroke="{color}""#)))
            .expect("an arc in that colour");
        let d = &path[..path.find('"').unwrap()];
        let (start, curve) = d.split_once('Q').unwrap();
        let baseline: f64 = start.split(' ').nth(1).unwrap().parse().unwrap();
        let control: f64 = curve.split(' ').nth(1).unwrap().parse().unwrap();
        (baseline + control) / 2.0
    }

    #[test]
    fn a_call_far_off_to_the_side_does_not_flatten_the_arcs_in_view() {
        // A translocation reaching 900 kb away is 450 times the width of the
        // deletion in view. Measured against it, the deletion arched to a
        // twentieth of the band and lay on the axis with the rest.
        let deletion = StructuralVariant::new(1_000, 3_000, SvKind::Deletion);
        let far = StructuralVariant::new(5_000, 900_000, SvKind::Translocation);
        let track = StructuralTrack::new(vec![deletion.clone(), far.clone()]);
        assert!(track.arch(&deletion) < 0.1, "against the whole track");
        assert_eq!(track.widest_in(&region()), 2_000);
        assert!((track.arch_in(&deletion, &region()) - 1.0).abs() < 1e-9);
        assert!(
            (track.arch_in(&far, &region()) - 1.0).abs() < 1e-9,
            "wider than the view's widest, so it reaches the top too"
        );

        // Drawn beside the translocation, the deletion still peaks within its
        // headroom of the top of the band.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![deletion, far]).show_footprints(false))
            .to_svg();
        let band_top = 14.0;
        let deletion_color = track.color_of(SvKind::Deletion, &Theme::light());
        let apex = apex_of(&svg, &deletion_color);
        assert!(apex >= band_top, "{apex} leaves the band");
        assert!(
            apex < band_top + 90.0 * 0.1,
            "{apex} does not use the height"
        );
    }

    #[test]
    fn a_view_inside_one_call_takes_its_scale_from_that_call() {
        let engulfing = StructuralVariant::new(1_000, 500_000, SvKind::Deletion);
        let track = StructuralTrack::new(vec![engulfing.clone()]);
        let inside = Region::new("chr1", 100_000, 200_000).unwrap();
        // Counted as no wider than the view it fills.
        assert_eq!(track.widest_in(&inside), 100_000);
        assert!((track.arch_in(&engulfing, &inside) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn an_insertion_does_not_set_the_scale_for_the_arcs_beside_it() {
        // Of no width, it would put every arc next to it at the top.
        let point = StructuralVariant::new(2_000, 2_000, SvKind::Insertion);
        let small = StructuralVariant::new(4_000, 4_100, SvKind::Deletion);
        let large = StructuralVariant::new(5_000, 9_000, SvKind::Deletion);
        let track = StructuralTrack::new(vec![point.clone(), small.clone(), large]);
        assert_eq!(track.widest_in(&region()), 4_000);
        assert!(track.arch_in(&small, &region()) < 0.5);
        assert!(track.arch_in(&point, &region()) >= 0.08);
        // Nor does it when the only call with width leaves the view, which is
        // where the fallback would otherwise pick it up.
        let leaving = StructuralVariant::new(8_000, 18_000, SvKind::Deletion);
        let beside = StructuralTrack::new(vec![point.clone(), leaving]);
        assert_eq!(beside.widest_in(&region()), 10_000);
        assert!(beside.arch_in(&point, &region()) < 0.5);
        // Insertions alone still stand up off the axis.
        let alone = StructuralTrack::new(vec![point.clone()]);
        assert!(alone.arch_in(&point, &region()) > 0.08);
    }

    #[test]
    fn a_call_zoomed_to_its_own_extent_is_inside_the_view() {
        // The call's end is exclusive, as the reader makes it, and so is the
        // region's, so the two meet rather than miss by one.
        let deletion = StructuralVariant::new(1_000, 9_000, SvKind::Deletion);
        let duplication = StructuralVariant::new(3_000, 3_500, SvKind::Duplication);
        let track = StructuralTrack::new(vec![deletion, duplication.clone()]);
        let exactly = Region::new("chr1", 1_000, 9_000).unwrap();
        assert_eq!(track.widest_in(&exactly), 8_000);
        assert!((track.arch_in(&duplication, &exactly) - 0.25).abs() < 1e-9);
    }

    #[test]
    fn calls_crossing_both_edges_are_measured_against_the_view() {
        // Nothing lies wholly inside: a deletion runs off the right, and a
        // translocation 1.2 Mb long runs off the left. Uncapped, the
        // translocation set the scale and the deletion arched to a tenth.
        let deletion = StructuralVariant::new(4_350_000, 4_359_600, SvKind::Deletion);
        let far = StructuralVariant::new(3_120_000, 4_345_000, SvKind::Translocation);
        let track = StructuralTrack::new(vec![deletion.clone(), far.clone()]);
        let view = Region::new("chr1", 4_340_000, 4_355_000).unwrap();
        assert_eq!(track.widest_in(&view), 15_000);
        assert!((track.arch_in(&deletion, &view) - 0.8).abs() < 1e-9);
        assert!((track.arch_in(&far, &view) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_name_that_is_never_drawn_takes_no_room() {
        let color = StructuralTrack::new(Vec::new()).color_of(SvKind::Deletion, &Theme::light());
        let apex = |view: Region, calls: Vec<StructuralVariant>| {
            let svg = Figure::new(view)
                .show_region_label(false)
                .push(StructuralTrack::new(calls).show_footprints(false))
                .to_svg();
            apex_of(&svg, &color)
        };
        let deletion = StructuralVariant::new(2_000, 8_000, SvKind::Deletion);
        let alone = apex(region(), vec![deletion.clone()]);
        // A name on an insertion is never written.
        let insertion = StructuralVariant::new(1_000, 1_000, SvKind::Insertion).name("INS1");
        assert_eq!(apex(region(), vec![deletion.clone(), insertion]), alone);
        // Nor is one centred off the plot, on a call that only touches the view.
        let view = Region::new("chr1", 10_000, 20_000).unwrap();
        let inside = StructuralVariant::new(12_000, 18_000, SvKind::Deletion);
        let offside = StructuralVariant::new(0, 10_000, SvKind::Duplication).name("OFFSCREEN");
        assert!(offside.touches(view.start(), view.end()));
        assert_eq!(
            apex(view.clone(), vec![inside.clone(), offside]),
            apex(view, vec![inside])
        );
        // A name that is drawn over the tallest arc does take its line.
        let named = apex(region(), vec![deletion.name("DEL6")]);
        assert!(named > alone + 5.0, "{named} against {alone}");
    }

    #[test]
    fn a_stroke_heavier_than_a_line_of_type_keeps_its_own_room_under_a_name() {
        // A named arc at the top with a thirty pixel stroke: the name's line
        // is fourteen, and half the stroke is fifteen.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                StructuralTrack::new(vec![StructuralVariant::new(1_000, 9_000, SvKind::Deletion)
                    .support(100)
                    .name("RD1")])
                .stroke(1.0, 30.0)
                .show_footprints(false),
            )
            .to_svg();
        let color = StructuralTrack::new(Vec::new()).color_of(SvKind::Deletion, &Theme::light());
        let band_top = 14.0;
        let apex = apex_of(&svg, &color);
        assert!(
            apex - 15.0 >= band_top,
            "the stroke reaches {}",
            apex - 15.0
        );
    }

    #[test]
    fn a_heavy_apex_marker_stays_inside_the_band_with_names_on() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                StructuralTrack::new(vec![
                    StructuralVariant::new(5_000, 5_000, SvKind::Insertion).support(100),
                    StructuralVariant::new(7_000, 7_000, SvKind::Insertion).name("INS2"),
                ])
                .stroke(1.0, 20.0),
            )
            .to_svg();
        let band_top = 14.0;
        for circle in svg.split("<circle").skip(1) {
            let read = |key: &str| -> f64 {
                let at = circle.find(&format!(" {key}=\"")).unwrap() + key.len() + 3;
                circle[at..at + circle[at..].find('"').unwrap()]
                    .parse()
                    .unwrap()
            };
            let top = read("cy") - read("r");
            assert!(top >= band_top, "a marker reaches {top}, above {band_top}");
        }
    }

    #[test]
    fn a_band_at_its_smallest_drops_a_name_rather_than_cut_it() {
        use crate::style::RenderProfile;
        for profile in [
            RenderProfile::Manuscript,
            RenderProfile::Compact,
            RenderProfile::Web,
        ] {
            let figure = Figure::new(region())
                .show_region_label(false)
                .profile(profile)
                .push(
                    StructuralTrack::new(vec![StructuralVariant::new(
                        1_000,
                        9_000,
                        SvKind::Deletion,
                    )
                    .name("RD1")])
                    .height(20.0),
                );
            let svg = figure.to_svg();
            let clip_top: f64 = svg
                .split("<clipPath")
                .nth(1)
                .and_then(|clip| clip.split(" y=\"").nth(1))
                .and_then(|y| y.split('"').next())
                .and_then(|y| y.parse().ok())
                .unwrap();
            let Some(text) = svg
                .split(">RD1</text>")
                .next()
                .filter(|_| svg.contains(">RD1<"))
            else {
                continue;
            };
            let at = text.rfind("<text").unwrap();
            let tag = &text[at..];
            let number = |key: &str| -> f64 {
                let from = tag.find(&format!(" {key}=\"")).unwrap() + key.len() + 3;
                tag[from..from + tag[from..].find('"').unwrap()]
                    .parse()
                    .unwrap()
            };
            let top = number("y") - number("font-size") * 0.75;
            assert!(
                top >= clip_top,
                "{profile:?}: the name starts at {top}, cut at {clip_top}"
            );
        }
    }

    #[test]
    fn the_name_over_the_tallest_arc_is_not_cut_by_the_top_of_the_band() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(StructuralTrack::new(vec![StructuralVariant::new(
                1_000,
                9_000,
                SvKind::Deletion,
            )
            .name("RD1")]))
            .to_svg();
        let text = svg.split(">RD1</text>").next().unwrap();
        let y: f64 = text
            .rsplit("y=\"")
            .next()
            .and_then(|rest| rest.split('"').next())
            .and_then(|y| y.parse().ok())
            .expect("the name has a baseline");
        let size = Theme::light().font_size - 2.0;
        let band_top = 14.0;
        assert!(
            y - size * 0.75 >= band_top,
            "the name's letters start at {} above a band starting at {band_top}",
            y - size * 0.75
        );
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
