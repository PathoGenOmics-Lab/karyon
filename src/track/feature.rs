//! Annotated intervals: genes, exons, repeats, primers, anything from a BED or
//! GFF file.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, Theme};
use crate::track::{DrawContext, Track};

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
    /// vanishing.
    pub fn new(start: u64, end: u64) -> Self {
        Feature {
            start,
            end: end.max(start + 1),
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
    fn layout(&self, scale: &Scale, theme: &Theme) -> (Vec<usize>, usize) {
        let mut rows = vec![0usize; self.features.len()];
        if self.features.is_empty() {
            return (rows, 1);
        }

        let mut order: Vec<usize> = (0..self.features.len()).collect();
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
        // everywhere else on the sheet. Drawing every gene in one accent left a
        // reverse gene wearing the colour that means forward in the pileup two
        // panels down, which is a quieter kind of wrong than a mislabelled axis
        // and just as misleading.
        let override_color = self.color.clone();
        let font = ctx.theme.font_size;

        for (i, feature) in self.features.iter().enumerate() {
            if feature.end <= ctx.region.start() || feature.start >= ctx.region.end() {
                continue;
            }

            let top = band.y + rows[i] as f64 * (self.row_height + self.row_gap);
            let bottom = top + self.row_height;
            let middle = (top + bottom) / 2.0;
            let left = ctx.scale.x(feature.start);
            let right = ctx.scale.x(feature.end).max(left + 1.5);
            let color = feature.color.clone().unwrap_or_else(|| {
                override_color
                    .clone()
                    .unwrap_or_else(|| strand_color(feature.strand, ctx.theme).to_string())
            });

            // The arrowhead eats a third of a short feature but never more
            // than 8 pixels of a long one, so a gene stays a bar with a point
            // rather than becoming a triangle.
            let head = ((right - left) * 0.35).min(8.0);
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
                    self.row_height,
                    ctx.theme.corner_radius,
                    &color,
                ),
            }

            if !self.show_names {
                continue;
            }
            let Some(name) = &feature.name else {
                continue;
            };
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn an_empty_track_still_has_one_row_of_height() {
        let region = Region::new("chr1", 0, 1000).unwrap();
        let track = FeatureTrack::new(Vec::new());
        assert_eq!(track.height(&scale(&region)), 14.0);
    }
}
