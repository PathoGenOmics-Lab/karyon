//! Annotated intervals from a BED or GFF file.
//!
//! A [`Feature`] is a span with an optional name, strand and colour: a gene, an
//! exon, a repeat, a primer. Coordinates are BED's, 0-based and half-open, so a
//! BED record goes in as it stands while a GFF or GenBank one has to be
//! converted first. The track's work is then to get a lot of them onto a few
//! rows without any two of them touching.
//!
//! # Rows follow the zoom, not the data
//!
//! Packing is first fit, leftmost first, and it is done in pixels rather than
//! in bases, so the same features take one row in a wide view and four in a
//! narrow one. A name too long to sit inside its feature is drawn to the right
//! of it on the same row, which is why the room a name needs is reserved during
//! the packing and not after it.
//!
//! Only the features in view are packed, which is also what sets the height, so
//! a cluster off the left edge cannot push the one gene on screen down a row.
//!
//! # Three places a colour can come from
//!
//! In order: the feature's own colour if it has one, then the track colour set
//! with [`FeatureTrack::color`], then [`strand_color`]. Falling through to the
//! last of those rather than to a track accent is what makes the default
//! useful, since five other tracks and the circular plots import that same
//! function from this module: a reverse feature drawn in the accent would wear
//! the colour that means forward in the pileup two bands down.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, Theme};
use crate::track::axis::group_thousands;
use crate::track::{DrawContext, Track};

/// A span as a reader reads it: 1-based, inclusive, thousands separated.
///
/// Coordinates are 0-based and half-open everywhere the crate computes, and a
/// tooltip is one of the two places that reaches a reader, the ruler being the
/// other. So `100..900` comes back as `101 to 900`, the coordinates the tick
/// labels print and the ones that go into a browser search box. It shares
/// [`group_thousands`] with the ruler for the same reason the ruler has one
/// unit across its whole width: two conventions on one figure have to be
/// decoded rather than read.
///
/// **This is the only place in the crate that writes a span.** Four tracks
/// once spelled `start + 1 .. end` out by hand and every one of them lost the
/// degenerate case with it: a zero-length interval has no last base to count
/// from one, so `100..100` came out `101 to 100`, a span running backwards on
/// a figure whose whole subject is direction. The floor below is what stops
/// that, and it only stops it for callers who come through here.
pub(crate) fn span_label(start: u64, end: u64) -> String {
    // Half-open in, inclusive out: the last base of `start..end` is `end - 1`,
    // which is `end` again once it is counted from one. The adds saturate
    // because a coordinate is a caller's number: counting the last base of the
    // range from one has nowhere to go, and a tooltip is not worth a panic.
    let last = end.max(start.saturating_add(1)) - 1;
    format!(
        "{} to {}",
        group_thousands(start.saturating_add(1)),
        group_thousands(last.saturating_add(1))
    )
}

/// How a strand is named in a tooltip, and nothing at all when it is unknown.
///
/// [`Strand::Unknown`] is drawn as a plain box precisely because there is
/// nothing to say about it, and a tooltip reading `unknown` would be a claim
/// where the glyph makes none.
pub(crate) fn strand_label(strand: Strand) -> &'static str {
    match strand {
        Strand::Forward => "forward",
        Strand::Reverse => "reverse",
        Strand::Unknown => "",
    }
}

/// What a reader hovering one feature is told: its name, its span, its strand.
///
/// The name leads because it is what was looked for. A feature without one
/// still gets its span, since where it is is the other half of the question and
/// the only half a nameless interval can answer.
///
/// A nameless one is given the noun `feature` in front of that span rather
/// than opening on a bare coordinate. Every tooltip in the crate is
/// `what it is, where it is`, and a name is what fills the first slot when
/// there is one; the fallback has to fill it too, or one glyph in a figure of
/// thirty answers a pointer in a different grammar from the rest.
pub(crate) fn feature_title(feature: &Feature) -> String {
    let mut title = String::new();
    match feature.name.as_deref().filter(|name| !name.is_empty()) {
        Some(name) => title.push_str(name),
        None => title.push_str("feature"),
    }
    title.push_str(", ");
    title.push_str(&span_label(feature.start, feature.end));
    let strand = strand_label(feature.strand);
    if !strand.is_empty() {
        title.push_str(", ");
        title.push_str(strand);
    }
    title
}

/// The colour a strand is drawn in, wherever a track colours by strand.
///
/// One convention for the whole crate, because it has to be. A figure holding
/// a read pileup with a methylation track under it, both coloured by strand,
/// would otherwise have blue meaning forward in one band and reverse in the
/// next, and nothing on the page would say so. An unknown strand is drawn as
/// forward, since a track that has to pick one may as well pick the common one.
///
/// A track that wants its own pair still has one: this is only the default.
pub fn strand_color(strand: Strand, theme: &Theme) -> &str {
    match strand {
        Strand::Reverse => theme.color(1),
        _ => theme.color(0),
    }
}

/// Which strand a feature sits on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Strand {
    /// Plus strand, drawn pointing right.
    Forward,
    /// Minus strand, drawn pointing left.
    Reverse,
    /// Unknown or not applicable, drawn as a plain box.
    #[default]
    Unknown,
}

impl Strand {
    /// Reads the `+`, `-` and `.` of BED and GFF, and anything else as
    /// [`Strand::Unknown`].
    pub fn from_symbol(symbol: char) -> Self {
        match symbol {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }
}

/// One annotated interval, in 0-based half-open coordinates.
///
/// A GFF file counts from 1 and includes its end, so a GFF line `start..end`
/// becomes `Feature::new(start - 1, end)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Feature {
    /// First base, 0-based.
    pub start: u64,
    /// One past the last base.
    pub end: u64,
    /// Label drawn on or beside the feature.
    pub name: Option<String>,
    /// Strand, which decides which way the arrow points.
    pub strand: Strand,
    /// Colour override, otherwise the track colour is used.
    pub color: Option<String>,
}

impl Feature {
    /// A feature spanning `start..end`.
    ///
    /// An end at or before the start is widened to a single base, so a
    /// zero-length record from a converter still shows up rather than silently
    /// vanishing. A start at the very top of the coordinate range has no next
    /// base to widen into, so it keeps the end it was given.
    pub fn new(start: u64, end: u64) -> Self {
        Feature {
            start,
            end: end.max(start.saturating_add(1)),
            name: None,
            strand: Strand::Unknown,
            color: None,
        }
    }

    /// Sets the label.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Sets the strand.
    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Sets a colour for this feature alone.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Length in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Always `false`: [`Feature::new`] guarantees at least one base.
    pub fn is_empty(&self) -> bool {
        false
    }
}

/// A row of features, packed so that nothing overlaps on screen.
///
/// Features that would collide are pushed onto extra rows, and the track grows
/// taller to fit them. Because collisions are measured in pixels and include
/// the space taken by labels, the number of rows changes with the zoom level:
/// this is why [`Track::height`] takes a [`Scale`].
///
/// ```
/// use karyon::{Feature, FeatureTrack, Figure, Region, Strand};
///
/// let genes = vec![
///     Feature::new(100, 900).name("rpoB").strand(Strand::Forward),
///     Feature::new(1200, 2000).name("katG").strand(Strand::Reverse),
/// ];
/// let svg = Figure::new(Region::parse("chr1:1-3000").unwrap())
///     .push(FeatureTrack::new(genes).label("genes"))
///     .to_svg();
/// assert!(svg.contains("rpoB"));
/// ```
#[derive(Debug, Clone)]
pub struct FeatureTrack {
    features: Vec<Feature>,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    color: Option<String>,
    show_names: bool,
}

impl FeatureTrack {
    /// A track holding `features`.
    pub fn new(features: impl Into<Vec<Feature>>) -> Self {
        FeatureTrack {
            features: features.into(),
            label: None,
            row_height: 14.0,
            row_gap: 3.0,
            color: None,
            show_names: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of a single row of features.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(2.0);
        self
    }

    /// Sets the vertical gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Sets the default colour for features without one of their own.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Draws or hides feature names.
    ///
    /// Hiding them also makes the track shorter, since names no longer take
    /// part in collision detection.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// The features in the track.
    pub fn features(&self) -> &[Feature] {
        &self.features
    }

    /// Assigns each feature to a row, first fit, left to right.
    ///
    /// Returns the row per feature in input order, plus the number of rows.
    ///
    /// Only what is on screen takes part, by the same test [`FeatureTrack::draw`]
    /// uses to skip a feature. Packing the rest decided the band height and the
    /// row of every visible feature from data the reader cannot see: three short
    /// features ending ten bases before the window opened reserved label room in
    /// pixels that reached into it, took the first three rows, and left the one
    /// gene in view floating on the fourth under three empty ones. A feature
    /// that merely overlaps an edge is kept, so its full pixel extent and its
    /// label room still count.
    fn layout(&self, scale: &Scale, theme: &Theme) -> (Vec<usize>, usize) {
        let mut rows = vec![0usize; self.features.len()];
        if self.features.is_empty() {
            return (rows, 1);
        }

        // The edges of the view as fractional positions rather than through
        // `Scale::bounds`, whose sum of the two overflows on a region running
        // to the top of the coordinate range. Every coordinate a genome uses is
        // far inside the range an f64 counts exactly.
        let view_start = scale.pos_at_x(scale.x0());
        let view_end = scale.pos_at_x(scale.x0() + scale.width());
        let mut order: Vec<usize> = (0..self.features.len())
            .filter(|&i| {
                self.features[i].end as f64 > view_start
                    && (self.features[i].start as f64) < view_end
            })
            .collect();
        order.sort_by_key(|&i| (self.features[i].start, self.features[i].end));

        // Horizontal breathing room between two features on the same row.
        let padding = 4.0;
        let mut row_ends: Vec<f64> = Vec::new();

        for &i in &order {
            let feature = &self.features[i];
            let left = scale.x(feature.start);
            let mut right = scale.x(feature.end).max(left + 2.0);
            if self.show_names {
                if let Some(name) = &feature.name {
                    let width = text_width(name, theme.font_size);
                    // A name that does not fit inside is drawn to the right,
                    // so it has to be reserved here or the next feature will
                    // sit on top of it.
                    if width + 6.0 > right - left {
                        right += width + 6.0;
                    }
                }
            }

            let mut placed = None;
            for (row, end) in row_ends.iter_mut().enumerate() {
                if *end + padding <= left {
                    *end = right;
                    placed = Some(row);
                    break;
                }
            }
            match placed {
                Some(row) => rows[i] = row,
                None => {
                    rows[i] = row_ends.len();
                    row_ends.push(right);
                }
            }
        }

        (rows, row_ends.len().max(1))
    }
}

impl Track for FeatureTrack {
    fn height(&self, scale: &Scale) -> f64 {
        // The theme only affects the height through the width of the labels,
        // and the default font size is what the figure will use unless the
        // caller changed it. Slight label crowding is a better failure mode
        // than threading the theme through every height computation.
        let (_, rows) = self.layout(scale, &Theme::default());
        rows as f64 * self.row_height + (rows.saturating_sub(1)) as f64 * self.row_gap
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let (rows, _) = self.layout(ctx.scale, ctx.theme);
        // Without an override, a feature takes the colour its strand takes
        // everywhere else on the sheet. Drawing every feature in one accent left
        // a reverse one wearing the colour that means forward in the pileup two
        // panels down, which is a quieter kind of wrong than a mislabelled axis
        // and just as misleading.
        let override_color = self.color.clone();
        let font = ctx.theme.font_size;
        let row_height = self.row_height * ctx.visual_scale;
        let row_gap = self.row_gap * ctx.visual_scale;

        for (i, feature) in self.features.iter().enumerate() {
            if feature.end <= ctx.region.start() || feature.start >= ctx.region.end() {
                continue;
            }

            let top = band.y + rows[i] as f64 * (row_height + row_gap);
            let bottom = top + row_height;
            let middle = (top + bottom) / 2.0;
            let left = ctx.scale.x(feature.start);
            let right = ctx.scale.x(feature.end).max(left + ctx.px(1.5));
            let color = feature.color.clone().unwrap_or_else(|| {
                override_color
                    .clone()
                    .unwrap_or_else(|| strand_color(feature.strand, ctx.theme).to_string())
            });

            // A gene with a name is the thing a reader most wants to point at,
            // so every feature carries one: the arrow and its label are one
            // glyph and answer together.
            //
            // Except when the feature is thinner than a pixel. `right` is
            // floored so a short feature stays visible, but a floor is not a
            // width: across a whole genome several thousand genes are a smear
            // a pointer cannot resolve, and naming each one would put a title
            // on a mark nobody can hit.
            let pointable = ctx.scale.x(feature.end) - left >= 1.0;
            if pointable {
                ctx.svg.begin_titled(&feature_title(feature));
            }

            // The arrowhead eats a third of a short feature but never more
            // than 8 pixels of a long one, so an interval stays a bar with a
            // point rather than becoming a triangle.
            let head = ((right - left) * 0.35).min(ctx.theme.tokens.arrow_size);
            match feature.strand {
                Strand::Forward if head > 1.0 => ctx.svg.polygon(
                    &[
                        (left, top),
                        (right - head, top),
                        (right, middle),
                        (right - head, bottom),
                        (left, bottom),
                    ],
                    &color,
                ),
                Strand::Reverse if head > 1.0 => ctx.svg.polygon(
                    &[
                        (right, top),
                        (left + head, top),
                        (left, middle),
                        (left + head, bottom),
                        (right, bottom),
                    ],
                    &color,
                ),
                _ => ctx.svg.rect_rounded(
                    left,
                    top,
                    right - left,
                    row_height,
                    ctx.theme.corner_radius,
                    &color,
                ),
            }

            // One exit from here on, so the group opened above is closed
            // exactly once however the name turns out.
            if let (true, Some(name)) = (self.show_names, &feature.name) {
                let width = text_width(name, font);
                let baseline = middle + font * 0.35;
                if width + 6.0 <= right - left {
                    ctx.svg.text(
                        (left + right) / 2.0,
                        baseline,
                        name,
                        contrast_ink(&color),
                        font,
                        Anchor::Middle,
                    );
                } else {
                    ctx.svg.text(
                        right + 3.0,
                        baseline,
                        name,
                        &ctx.theme.foreground,
                        font,
                        Anchor::Start,
                    );
                }
            }

            if pointable {
                ctx.svg.end_group();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn scale(region: &Region) -> Scale {
        Scale::new(region, 0.0, 1000.0)
    }

    #[test]
    fn zero_length_features_are_widened_to_one_base() {
        assert_eq!(Feature::new(100, 100).end, 101);
        assert_eq!(Feature::new(100, 50).end, 101);
        assert_eq!(Feature::new(100, 100).len(), 1);
    }

    #[test]
    fn strand_reads_bed_and_gff_symbols() {
        assert_eq!(Strand::from_symbol('+'), Strand::Forward);
        assert_eq!(Strand::from_symbol('-'), Strand::Reverse);
        assert_eq!(Strand::from_symbol('.'), Strand::Unknown);
        assert_eq!(Strand::from_symbol('?'), Strand::Unknown);
    }

    #[test]
    fn disjoint_features_share_one_row() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let track = FeatureTrack::new(vec![
            Feature::new(0, 100),
            Feature::new(300, 400),
            Feature::new(700, 800),
        ])
        .show_names(false);
        let (rows, count) = track.layout(&scale(&region), &Theme::default());
        assert_eq!(count, 1);
        assert_eq!(rows, vec![0, 0, 0]);
    }

    #[test]
    fn overlapping_features_stack_onto_extra_rows() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let track = FeatureTrack::new(vec![
            Feature::new(0, 500),
            Feature::new(100, 600),
            Feature::new(200, 700),
        ])
        .show_names(false);
        let (rows, count) = track.layout(&scale(&region), &Theme::default());
        assert_eq!(count, 3);
        assert_eq!(rows, vec![0, 1, 2]);
    }

    #[test]
    fn rows_are_assigned_in_input_order_not_sorted_order() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let track =
            FeatureTrack::new(vec![Feature::new(200, 700), Feature::new(0, 500)]).show_names(false);
        let (rows, count) = track.layout(&scale(&region), &Theme::default());
        assert_eq!(count, 2);
        // The leftmost feature takes row 0 even though it is listed second.
        assert_eq!(rows, vec![1, 0]);
    }

    #[test]
    fn labels_take_part_in_collision_detection() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let features = vec![
            Feature::new(0, 20).name("a_very_long_gene_name_here"),
            Feature::new(30, 50).name("another_long_gene_name"),
        ];
        let unnamed = FeatureTrack::new(features.clone()).show_names(false);
        let named = FeatureTrack::new(features);
        assert_eq!(unnamed.layout(&scale(&region), &Theme::default()).1, 1);
        assert_eq!(named.layout(&scale(&region), &Theme::default()).1, 2);
    }

    #[test]
    fn height_grows_with_the_number_of_rows() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let one = FeatureTrack::new(vec![Feature::new(0, 100)]).show_names(false);
        let two =
            FeatureTrack::new(vec![Feature::new(0, 500), Feature::new(100, 600)]).show_names(false);
        let s = scale(&region);
        assert_eq!(one.height(&s), 14.0);
        assert_eq!(two.height(&s), 14.0 * 2.0 + 3.0);
    }

    #[test]
    fn a_feature_off_the_left_edge_does_not_push_the_one_in_view_down_a_row() {
        // Three 10 bp features ending ten bases before the window opens are
        // never drawn, but their names reserved pixels that reached into it,
        // took rows 0 to 2 and left rpoB alone on row 3 under three empty ones.
        let region = Region::new("NC_000962.3", 761_000, 763_000).unwrap();
        let s = scale(&region);
        let visible = Feature::new(761_200, 762_400).name("rpoB");
        let mut all = vec![visible.clone()];
        all.extend((0..3u64).map(|i| {
            Feature::new(760_900 + i, 760_910 + i)
                .name("a_gene_with_a_very_long_name_indeed_xxxxxxxx")
        }));

        let alone = FeatureTrack::new(vec![visible]);
        let together = FeatureTrack::new(all);
        assert_eq!(alone.height(&s), 14.0);
        assert_eq!(together.height(&s), 14.0, "65.0 before the view filter");
        assert_eq!(together.layout(&s, &Theme::default()).0[0], 0);
    }

    #[test]
    fn a_feature_hanging_over_the_edge_still_takes_part_in_the_packing() {
        // Overlapping the view is enough to be packed, since the part of it on
        // screen is drawn and can be collided with.
        let region = Region::new("chr1", 1_000, 2_000).unwrap();
        let s = scale(&region);
        let track = FeatureTrack::new(vec![Feature::new(900, 1_500), Feature::new(1_200, 1_800)])
            .show_names(false);
        let (rows, count) = track.layout(&s, &Theme::default());
        assert_eq!(count, 2);
        assert_eq!(rows, vec![0, 1]);
    }

    #[test]
    fn a_feature_at_the_top_of_the_coordinate_range_keeps_the_end_it_was_given() {
        // `end.max(start + 1)` overflowed, so a record at u64::MAX aborted the
        // render instead of being drawn off screen.
        let feature = Feature::new(u64::MAX, u64::MAX);
        assert_eq!(feature.start, u64::MAX);
        assert_eq!(feature.end, u64::MAX);
        assert_eq!(
            span_label(u64::MAX, u64::MAX),
            "18,446,744,073,709,551,615 to 18,446,744,073,709,551,615"
        );
        let svg = Figure::new(Region::new("chr1", 0, 1000).unwrap())
            .push(FeatureTrack::new(vec![feature]))
            .to_svg();
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn a_region_running_to_the_top_of_the_coordinate_range_is_drawn_rather_than_refused() {
        // A guard rather than a closed defect: the new view filter reads the
        // edges off the scale as fractional positions, and `Scale::bounds`
        // adds its start to its span, which overflows at the top of the range.
        // Packing is not worth a panic, so the arithmetic has to be saturating.
        for (start, end) in [(0u64, u64::MAX), (1, u64::MAX), (u64::MAX - 10, u64::MAX)] {
            let region = Region::new("chr1", start, end).unwrap();
            let svg = Figure::new(region)
                .push(FeatureTrack::new(vec![Feature::new(start, start + 5)]))
                .to_svg();
            assert!(svg.contains("</svg>"), "{start}..{end}");
        }
    }

    #[test]
    fn an_empty_track_still_has_one_row_of_height() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let track = FeatureTrack::new(Vec::new());
        assert_eq!(track.height(&scale(&region)), 14.0);
    }

    #[test]
    fn a_span_is_quoted_the_way_the_ruler_prints_it() {
        // 0-based half-open in, 1-based inclusive out, grouped like the ticks.
        assert_eq!(span_label(759_806, 763_325), "759,807 to 763,325");
        assert_eq!(span_label(0, 1), "1 to 1");
        assert_eq!(span_label(99, 100), "100 to 100");
    }

    #[test]
    fn a_zero_length_span_still_reads_forwards() {
        // Written by hand as `start + 1 to end` this is `101 to 100`, a span
        // running backwards on a figure whose subject is direction. Four
        // tracks wrote it by hand and all four had the bug, which is why this
        // is the only span formatter left.
        assert_eq!(span_label(100, 100), "101 to 101");
        assert_eq!(span_label(0, 0), "1 to 1");
        // An end before the start collapses to the start rather than inverting.
        assert_eq!(span_label(100, 50), "101 to 101");
    }

    #[test]
    fn a_gene_is_named_by_its_name_its_span_and_its_strand() {
        let svg = Figure::new(Region::parse("NC_000962.3:759001-764000").unwrap())
            .show_region_label(false)
            .push(FeatureTrack::new(vec![Feature::new(759_806, 763_325)
                .name("rpoB")
                .strand(Strand::Forward)]))
            .to_svg();
        assert!(
            svg.contains("<title>rpoB, 759,807 to 763,325, forward</title>"),
            "{svg}"
        );
    }

    #[test]
    fn a_reverse_gene_says_so_and_a_strandless_one_says_nothing() {
        let render = |feature: Feature| {
            Figure::new(Region::new("chr1", 0, 2_000).unwrap())
                .show_region_label(false)
                .push(FeatureTrack::new(vec![feature]))
                .to_svg()
        };
        assert!(
            render(Feature::new(100, 900).name("katG").strand(Strand::Reverse))
                .contains("<title>katG, 101 to 900, reverse</title>")
        );
        // Unknown is drawn as a plain box because there is nothing to say, and
        // a tooltip reading `unknown` would be a claim the glyph does not make.
        assert!(
            render(Feature::new(100, 900).name("katG")).contains("<title>katG, 101 to 900</title>")
        );
    }

    #[test]
    fn a_feature_with_no_name_still_gets_its_span() {
        let svg = Figure::new(Region::new("chr1", 0, 2_000).unwrap())
            .show_region_label(false)
            .push(FeatureTrack::new(vec![
                Feature::new(100, 900).strand(Strand::Forward)
            ]))
            .to_svg();
        // The noun stands in for the missing name, so the tooltip still opens
        // on what it is rather than on a bare coordinate.
        assert!(
            svg.contains("<title>feature, 101 to 900, forward</title>"),
            "{svg}"
        );
    }

    #[test]
    fn a_feature_thinner_than_a_pixel_is_not_named() {
        // Two thousand genes across a bacterial genome are a smear a pointer
        // cannot resolve, and a title on each would be a title on a mark
        // nobody can hit as well as a quarter of the file.
        let genes: Vec<Feature> = (0..2_000)
            .map(|i| Feature::new(i * 2_000, i * 2_000 + 300).name(format!("g{i}")))
            .collect();
        let wide = Figure::new(Region::new("chr1", 0, 4_000_000).unwrap())
            .show_region_label(false)
            .push(FeatureTrack::new(genes.clone()).show_names(false))
            .to_svg();
        assert!(!wide.contains("<title>"), "a sub-pixel gene was named");

        // Zoomed in far enough to point at one, it is named again.
        let close = Figure::new(Region::new("chr1", 0, 20_000).unwrap())
            .show_region_label(false)
            .push(FeatureTrack::new(genes).show_names(false))
            .to_svg();
        assert!(close.contains("<title>g0, 1 to 300</title>"), "{close}");
    }

    #[test]
    fn every_group_a_feature_opens_is_closed_again() {
        let svg = Figure::new(Region::new("chr1", 0, 3_000).unwrap())
            .show_region_label(false)
            .push(FeatureTrack::new(vec![
                // Named and wide, named and narrow, and nameless.
                Feature::new(0, 900).name("wide").strand(Strand::Forward),
                Feature::new(1_000, 1_020)
                    .name("a_name_far_wider_than_its_feature")
                    .strand(Strand::Reverse),
                Feature::new(2_000, 2_800),
            ]))
            .to_svg();
        assert_eq!(
            svg.matches("<g").count(),
            svg.matches("</g>").count(),
            "{svg}"
        );
        assert_eq!(svg.matches("<title>").count(), 3);
    }
}
