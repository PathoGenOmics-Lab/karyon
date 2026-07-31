//! Two sequences compared: as a dotplot, or as ribbons between two bars.
//!
//! Both tracks draw the same thing, a set of alignment blocks, and differ only
//! in where they put the second sequence. A dotplot gives it the vertical axis,
//! which shows the shape of a rearrangement at a glance and costs a tall panel.
//! Ribbons give it a second bar under the first, which is compact and follows
//! one block at a time. Neither is a summary of the other, so both are here.
//!
//! The figure's region is always the query. The target is whatever the blocks
//! say it is, on its own scale.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Track};

/// One alignment between a stretch of the query and a stretch of the target.
///
/// Coordinates are 0-based half-open on both sequences, and both spans run
/// forwards: a reversed alignment is `target` ascending with `reversed` set,
/// which is how PAF records it and saves a conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentBlock {
    /// Start on the query, the sequence the figure's region is on.
    pub query_start: u64,
    /// End on the query.
    pub query_end: u64,
    /// Start on the target, the other sequence.
    pub target_start: u64,
    /// End on the target.
    pub target_end: u64,
    /// Whether the target runs backwards against the query.
    pub reversed: bool,
    /// Alignment identity from 0 to 1, if known.
    pub identity: Option<f64>,
}

impl AlignmentBlock {
    /// A forward block.
    pub fn new(query_start: u64, query_end: u64, target_start: u64, target_end: u64) -> Self {
        AlignmentBlock {
            query_start,
            query_end: query_end.max(query_start),
            target_start,
            target_end: target_end.max(target_start),
            reversed: false,
            identity: None,
        }
    }

    /// Marks the block as running backwards on the target.
    pub fn reversed(mut self, reversed: bool) -> Self {
        self.reversed = reversed;
        self
    }

    /// Sets the identity, from 0 to 1.
    pub fn identity(mut self, identity: f64) -> Self {
        self.identity = Some(identity.clamp(0.0, 1.0));
        self
    }

    /// Length of the block on the query.
    pub fn query_len(&self) -> u64 {
        self.query_end - self.query_start
    }

    /// Length of the block on the target.
    pub fn target_len(&self) -> u64 {
        self.target_end - self.target_start
    }

    /// Whether the block touches a query span.
    pub fn overlaps_query(&self, start: u64, end: u64) -> bool {
        self.query_start < end && self.query_end > start
    }
}

/// Where the second sequence is measured from and to.
#[derive(Debug, Clone, Copy, PartialEq)]
struct TargetRange {
    start: u64,
    end: u64,
}

impl TargetRange {
    fn fraction(&self, position: u64) -> f64 {
        let span = (self.end - self.start).max(1) as f64;
        ((position.max(self.start) - self.start) as f64 / span).clamp(0.0, 1.0)
    }
}

/// Picks the target span from what the blocks actually use, when none is given.
fn spanning_range(blocks: &[AlignmentBlock]) -> TargetRange {
    let start = blocks.iter().map(|b| b.target_start).min().unwrap_or(0);
    let end = blocks.iter().map(|b| b.target_end).max().unwrap_or(1);
    TargetRange {
        start,
        end: end.max(start + 1),
    }
}

/// Colour for a block, by orientation, which is the thing a reader is looking
/// for in a comparison.
fn block_color(
    block: &AlignmentBlock,
    theme: &Theme,
    forward: &Option<String>,
    reverse: &Option<String>,
) -> String {
    if block.reversed {
        reverse
            .clone()
            .unwrap_or_else(|| theme.color(1).to_string())
    } else {
        forward
            .clone()
            .unwrap_or_else(|| theme.color(0).to_string())
    }
}

/// Two sequences on two axes, alignments as diagonals.
///
/// A forward block runs bottom left to top right, a reversed one the other way,
/// and a rearrangement is whatever shape those make. The query is the figure's
/// region; the target gets the height of the band.
///
/// ```
/// use karyon::{AlignmentBlock, DotplotTrack, Figure, Region};
///
/// let blocks = vec![
///     AlignmentBlock::new(0, 4_000, 0, 4_000),
///     AlignmentBlock::new(4_000, 6_000, 6_000, 8_000).reversed(true),
/// ];
///
/// let svg = Figure::new(Region::new("query", 0, 10_000).unwrap())
///     .push(DotplotTrack::new(blocks).target_length(10_000).label("target"))
///     .to_svg();
/// assert!(svg.contains("<line"));
/// ```
#[derive(Debug, Clone)]
pub struct DotplotTrack {
    blocks: Vec<AlignmentBlock>,
    label: Option<String>,
    height: f64,
    target: Option<TargetRange>,
    forward_color: Option<String>,
    reverse_color: Option<String>,
    line_width: f64,
    show_frame: bool,
    show_scale: bool,
}

impl DotplotTrack {
    /// A dotplot of `blocks`.
    pub fn new(blocks: impl Into<Vec<AlignmentBlock>>) -> Self {
        DotplotTrack {
            blocks: blocks.into(),
            label: None,
            height: 190.0,
            target: None,
            forward_color: None,
            reverse_color: None,
            line_width: 1.6,
            show_frame: true,
            show_scale: true,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    ///
    /// A dotplot wants to be tall. The diagonal of a colinear alignment only
    /// looks like a diagonal when the panel is not a letterbox.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(20.0);
        self
    }

    /// Puts the whole target sequence on the vertical axis.
    pub fn target_length(mut self, length: u64) -> Self {
        self.target = Some(TargetRange {
            start: 0,
            end: length.max(1),
        });
        self
    }

    /// Puts one span of the target on the vertical axis.
    pub fn target_range(mut self, start: u64, end: u64) -> Self {
        self.target = Some(TargetRange {
            start,
            end: end.max(start + 1),
        });
        self
    }

    /// Sets the colour of forward blocks.
    pub fn forward_color(mut self, color: impl Into<String>) -> Self {
        self.forward_color = Some(color.into());
        self
    }

    /// Sets the colour of reversed blocks.
    pub fn reverse_color(mut self, color: impl Into<String>) -> Self {
        self.reverse_color = Some(color.into());
        self
    }

    /// Sets the stroke width of a diagonal.
    pub fn line_width(mut self, width: f64) -> Self {
        self.line_width = width.max(0.2);
        self
    }

    /// Draws or hides the box around the panel.
    pub fn show_frame(mut self, show: bool) -> Self {
        self.show_frame = show;
        self
    }

    /// Draws or hides the target coordinate labels.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// The blocks in the track.
    pub fn blocks(&self) -> &[AlignmentBlock] {
        &self.blocks
    }

    /// The target span on the vertical axis, given or inferred from the blocks.
    fn range(&self) -> TargetRange {
        self.target.unwrap_or_else(|| spanning_range(&self.blocks))
    }
}

impl Track for DotplotTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale || self.blocks.is_empty() {
            return 0.0;
        }
        let range = self.range();
        text_width(&format_position(range.end), theme.font_size - 1.0) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let range = self.range();

        if self.show_frame {
            ctx.svg
                .rect_outline(band.x, band.y, band.w, band.h, &ctx.theme.rule, 1.0);
        }

        // The target grows upwards, so the largest coordinate is at the top of
        // the band and the origin is at the bottom, as on any other y axis.
        let y_of = |position: u64| band.bottom() - range.fraction(position) * band.h;

        for block in &self.blocks {
            if !block.overlaps_query(ctx.region.start(), ctx.region.end()) {
                continue;
            }
            let color = block_color(block, ctx.theme, &self.forward_color, &self.reverse_color);
            let x0 = ctx.scale.x(block.query_start);
            let x1 = ctx.scale.x(block.query_end);
            // A reversed block descends as the query advances, which is the
            // anti-diagonal that makes an inversion recognisable.
            let (y0, y1) = if block.reversed {
                (y_of(block.target_end), y_of(block.target_start))
            } else {
                (y_of(block.target_start), y_of(block.target_end))
            };
            ctx.svg.line(x0, y0, x1, y1, &color, self.line_width);
        }

        if self.show_scale && ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            ctx.svg.text(
                right,
                band.y + size,
                &format_position(range.end),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
            ctx.svg.text(
                right,
                band.bottom(),
                &format_position(range.start),
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }
    }
}

/// Two bars joined by ribbons, one per alignment block.
///
/// The compact form of the same comparison. A reversed block crosses its own
/// ribbon, so an inversion is a twist rather than a shape you have to read off
/// two axes.
///
/// ```
/// use karyon::{AlignmentBlock, Figure, Region, SyntenyTrack};
///
/// let blocks = vec![AlignmentBlock::new(0, 4_000, 1_000, 5_000)];
///
/// let svg = Figure::new(Region::new("query", 0, 10_000).unwrap())
///     .push(
///         SyntenyTrack::new(blocks)
///             .target_length(10_000)
///             .names("H37Rv", "CDC1551"),
///     )
///     .to_svg();
/// assert!(svg.contains("H37Rv"));
/// ```
#[derive(Debug, Clone)]
pub struct SyntenyTrack {
    blocks: Vec<AlignmentBlock>,
    label: Option<String>,
    height: f64,
    bar_height: f64,
    target: Option<TargetRange>,
    forward_color: Option<String>,
    reverse_color: Option<String>,
    ribbon_opacity: f64,
    query_name: Option<String>,
    target_name: Option<String>,
    min_ribbon_width: f64,
}

impl SyntenyTrack {
    /// A ribbon plot of `blocks`.
    pub fn new(blocks: impl Into<Vec<AlignmentBlock>>) -> Self {
        SyntenyTrack {
            blocks: blocks.into(),
            label: None,
            height: 110.0,
            bar_height: 9.0,
            target: None,
            forward_color: None,
            reverse_color: None,
            ribbon_opacity: 0.45,
            query_name: None,
            target_name: None,
            min_ribbon_width: 1.5,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(24.0);
        self
    }

    /// Sets the thickness of the two sequence bars.
    pub fn bar_height(mut self, height: f64) -> Self {
        self.bar_height = height.max(2.0);
        self
    }

    /// Puts the whole target sequence on the lower bar.
    pub fn target_length(mut self, length: u64) -> Self {
        self.target = Some(TargetRange {
            start: 0,
            end: length.max(1),
        });
        self
    }

    /// Puts one span of the target on the lower bar.
    pub fn target_range(mut self, start: u64, end: u64) -> Self {
        self.target = Some(TargetRange {
            start,
            end: end.max(start + 1),
        });
        self
    }

    /// Names the two sequences, drawn beside their bars.
    pub fn names(mut self, query: impl Into<String>, target: impl Into<String>) -> Self {
        self.query_name = Some(query.into());
        self.target_name = Some(target.into());
        self
    }

    /// Sets the colour of forward ribbons.
    pub fn forward_color(mut self, color: impl Into<String>) -> Self {
        self.forward_color = Some(color.into());
        self
    }

    /// Sets the colour of reversed ribbons.
    pub fn reverse_color(mut self, color: impl Into<String>) -> Self {
        self.reverse_color = Some(color.into());
        self
    }

    /// Sets how solid a ribbon is, from 0 to 1.
    ///
    /// Ribbons overlap wherever blocks do, so they are drawn translucent: two
    /// crossing ribbons have to read as two, and a stack of them has to show
    /// that it is a stack.
    pub fn ribbon_opacity(mut self, opacity: f64) -> Self {
        self.ribbon_opacity = opacity.clamp(0.0, 1.0);
        self
    }

    /// The blocks in the track.
    pub fn blocks(&self) -> &[AlignmentBlock] {
        &self.blocks
    }

    /// The target span on the lower bar, given or inferred from the blocks.
    fn range(&self) -> TargetRange {
        self.target.unwrap_or_else(|| spanning_range(&self.blocks))
    }
}

impl Track for SyntenyTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let widest = [&self.query_name, &self.target_name]
            .iter()
            .filter_map(|name| name.as_ref())
            .map(|name| text_width(name, theme.font_size - 1.0))
            .fold(0.0f64, f64::max);
        if widest <= 0.0 {
            0.0
        } else {
            widest + 8.0
        }
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let range = self.range();
        let query_top = band.y;
        let target_top = band.bottom() - self.bar_height;
        let query_bottom = query_top + self.bar_height;

        let bar = ctx.theme.rule.clone();
        ctx.svg.rect_rounded(
            band.x,
            query_top,
            band.w,
            self.bar_height,
            ctx.theme.corner_radius,
            &bar,
        );
        ctx.svg.rect_rounded(
            band.x,
            target_top,
            band.w,
            self.bar_height,
            ctx.theme.corner_radius,
            &bar,
        );

        let x_target = |position: u64| band.x + range.fraction(position) * band.w;
        let middle = (query_bottom + target_top) / 2.0;

        for block in &self.blocks {
            if !block.overlaps_query(ctx.region.start(), ctx.region.end()) {
                continue;
            }
            let color = block_color(block, ctx.theme, &self.forward_color, &self.reverse_color);

            let qx0 = ctx.scale.x(block.query_start);
            let qx1 = ctx
                .scale
                .x(block.query_end)
                .max(qx0 + self.min_ribbon_width);
            // A reversed block joins the start of one span to the end of the
            // other, which is what puts the twist in the ribbon.
            let (tx0, tx1) = if block.reversed {
                (x_target(block.target_end), x_target(block.target_start))
            } else {
                (x_target(block.target_start), x_target(block.target_end))
            };
            let (tx0, tx1) = if (tx1 - tx0).abs() < self.min_ribbon_width {
                (tx0, tx0 + self.min_ribbon_width.copysign(tx1 - tx0 + 1e-9))
            } else {
                (tx0, tx1)
            };

            // Both edges leave their bar vertically and arrive vertically, so
            // the ribbon reads as a connection rather than a wedge.
            let d = format!(
                "M{} {} C{} {} {} {} {} {} L{} {} C{} {} {} {} {} {} Z",
                num(qx0),
                num(query_bottom),
                num(qx0),
                num(middle),
                num(tx0),
                num(middle),
                num(tx0),
                num(target_top),
                num(tx1),
                num(target_top),
                num(tx1),
                num(middle),
                num(qx1),
                num(middle),
                num(qx1),
                num(query_bottom)
            );
            ctx.svg.path(&d, &color, self.ribbon_opacity);

            // The blocks themselves, solid on their bars, so a thin ribbon
            // still shows exactly what it connects.
            ctx.svg
                .rect(qx0, query_top, qx1 - qx0, self.bar_height, &color);
            let (left, width) = if tx1 >= tx0 {
                (tx0, tx1 - tx0)
            } else {
                (tx1, tx0 - tx1)
            };
            ctx.svg
                .rect(left, target_top, width, self.bar_height, &color);
        }

        if ctx.axis.w > 0.0 {
            let size = ctx.theme.font_size - 1.0;
            let right = ctx.axis.right() - 4.0;
            if let Some(name) = &self.query_name {
                ctx.svg.text(
                    right,
                    query_top + self.bar_height / 2.0 + size * 0.35,
                    name,
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
            if let Some(name) = &self.target_name {
                ctx.svg.text(
                    right,
                    target_top + self.bar_height / 2.0 + size * 0.35,
                    name,
                    &ctx.theme.muted,
                    size,
                    Anchor::End,
                );
            }
        }
    }
}

/// A coordinate shortened for an axis label.
fn format_position(position: u64) -> String {
    if position >= 1_000_000 {
        let value = position as f64 / 1e6;
        if value == value.trunc() {
            format!("{} Mb", value as i64)
        } else {
            format!("{value:.1} Mb")
        }
    } else if position >= 1_000 {
        let value = position as f64 / 1e3;
        if value == value.trunc() {
            format!("{} kb", value as i64)
        } else {
            format!("{value:.1} kb")
        }
    } else {
        format!("{position}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("query", 0, 10_000).unwrap()
    }

    fn blocks() -> Vec<AlignmentBlock> {
        vec![
            AlignmentBlock::new(0, 4_000, 0, 4_000),
            AlignmentBlock::new(4_000, 6_000, 6_000, 8_000).reversed(true),
            AlignmentBlock::new(6_000, 9_000, 4_000, 7_000),
        ]
    }

    #[test]
    fn a_block_keeps_both_spans_ascending() {
        let block = AlignmentBlock::new(100, 50, 900, 400);
        assert_eq!(
            block.query_end, 100,
            "an inverted span collapses, not wraps"
        );
        assert_eq!(block.target_end, 900);
        assert_eq!(AlignmentBlock::new(0, 400, 0, 300).query_len(), 400);
        assert_eq!(AlignmentBlock::new(0, 400, 0, 300).target_len(), 300);
    }

    #[test]
    fn identity_is_clamped_to_a_fraction() {
        assert_eq!(
            AlignmentBlock::new(0, 1, 0, 1).identity(1.5).identity,
            Some(1.0)
        );
        assert_eq!(
            AlignmentBlock::new(0, 1, 0, 1).identity(-1.0).identity,
            Some(0.0)
        );
        assert_eq!(AlignmentBlock::new(0, 1, 0, 1).identity, None);
    }

    #[test]
    fn overlap_is_half_open_on_the_query() {
        let block = AlignmentBlock::new(100, 200, 0, 100);
        assert!(block.overlaps_query(150, 250));
        assert!(block.overlaps_query(0, 101));
        assert!(!block.overlaps_query(200, 300), "the end is exclusive");
        assert!(!block.overlaps_query(0, 100));
    }

    #[test]
    fn an_unspecified_target_spans_what_the_blocks_use() {
        let track = DotplotTrack::new(blocks());
        let range = track.range();
        assert_eq!(range.start, 0);
        assert_eq!(range.end, 8_000);
    }

    #[test]
    fn a_given_target_wins_over_the_inferred_one() {
        let whole = DotplotTrack::new(blocks()).target_length(50_000);
        assert_eq!(whole.range().end, 50_000);
        let window = DotplotTrack::new(blocks()).target_range(1_000, 2_000);
        assert_eq!(
            window.range(),
            TargetRange {
                start: 1_000,
                end: 2_000
            }
        );
    }

    #[test]
    fn an_empty_dotplot_infers_a_range_without_dividing_by_zero() {
        let track = DotplotTrack::new(Vec::new());
        let range = track.range();
        assert!(range.end > range.start);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn the_target_axis_puts_its_origin_at_the_bottom() {
        let range = TargetRange { start: 0, end: 100 };
        assert_eq!(range.fraction(0), 0.0);
        assert_eq!(range.fraction(50), 0.5);
        assert_eq!(range.fraction(100), 1.0);
        // Outside the range is clamped rather than drawn off the panel.
        assert_eq!(range.fraction(500), 1.0);
    }

    #[test]
    fn forward_and_reversed_blocks_take_different_colours() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                DotplotTrack::new(blocks())
                    .target_length(10_000)
                    .forward_color("#111111")
                    .reverse_color("#ff0000")
                    .show_scale(false),
            )
            .to_svg();
        assert!(svg.contains("#111111"));
        assert!(svg.contains("#ff0000"));
    }

    #[test]
    fn a_reversed_diagonal_descends_where_a_forward_one_climbs() {
        // Read the two endpoints of each diagonal out of the SVG and compare.
        let endpoints = |reversed: bool| {
            let block = AlignmentBlock::new(0, 5_000, 0, 5_000).reversed(reversed);
            let svg = Figure::new(region())
                .show_region_label(false)
                .push(
                    DotplotTrack::new(vec![block])
                        .target_length(10_000)
                        .show_frame(false)
                        .show_scale(false),
                )
                .to_svg();
            let start = svg.find("<line ").unwrap();
            let element = &svg[start..start + svg[start..].find("/>").unwrap()];
            let read = |key: &str| -> f64 {
                let at = element.find(key).unwrap() + key.len();
                let rest = &element[at..];
                rest[..rest.find('"').unwrap()].parse().unwrap()
            };
            (read("y1=\""), read("y2=\""))
        };

        let (fy1, fy2) = endpoints(false);
        let (ry1, ry2) = endpoints(true);
        assert!(fy2 < fy1, "a forward block climbs: {fy1} to {fy2}");
        assert!(ry2 > ry1, "a reversed block descends: {ry1} to {ry2}");
    }

    #[test]
    fn blocks_outside_the_query_region_are_not_drawn() {
        let far = vec![AlignmentBlock::new(90_000, 95_000, 0, 5_000)];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(
                DotplotTrack::new(far)
                    .target_length(10_000)
                    .show_frame(false),
            )
            .to_svg();
        assert!(!svg.contains("<line"));
    }

    #[test]
    fn the_dotplot_labels_both_ends_of_the_target() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(DotplotTrack::new(blocks()).target_range(0, 2_000_000))
            .to_svg();
        assert!(svg.contains("2 Mb"));
        assert!(svg.contains(">0</text>"));
    }

    #[test]
    fn synteny_names_size_the_axis_strip_and_get_drawn() {
        let theme = Theme::light();
        let anonymous = SyntenyTrack::new(blocks());
        let named = SyntenyTrack::new(blocks()).names("H37Rv", "CDC1551");
        assert_eq!(anonymous.y_axis_width(&theme), 0.0);
        assert!(named.y_axis_width(&theme) > 0.0);

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(named.target_length(10_000))
            .to_svg();
        assert!(svg.contains("H37Rv"));
        assert!(svg.contains("CDC1551"));
    }

    #[test]
    fn every_visible_block_gets_a_ribbon() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(SyntenyTrack::new(blocks()).target_length(10_000))
            .to_svg();
        assert_eq!(svg.matches("<path").count(), 3);
    }

    #[test]
    fn ribbons_are_translucent_so_a_crossing_reads_as_two() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(SyntenyTrack::new(blocks()).target_length(10_000))
            .to_svg();
        assert!(svg.contains("fill-opacity"));
    }

    #[test]
    fn an_empty_synteny_track_still_draws_its_two_bars() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(SyntenyTrack::new(Vec::new()).target_length(1_000))
            .to_svg();
        assert!(!svg.contains("<path"));
        // The page background, the two bars, and the clip path's own rect.
        assert_eq!(svg.matches("<rect").count(), 4);
    }

    #[test]
    fn positions_shorten_for_an_axis_label() {
        assert_eq!(format_position(0), "0");
        assert_eq!(format_position(999), "999");
        assert_eq!(format_position(4_000), "4 kb");
        assert_eq!(format_position(4_500), "4.5 kb");
        assert_eq!(format_position(2_000_000), "2 Mb");
        assert_eq!(format_position(4_411_532), "4.4 Mb");
    }
}
