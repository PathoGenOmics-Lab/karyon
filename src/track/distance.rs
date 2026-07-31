//! A symmetric distance matrix, ordered by the tree beside it.
//!
//! # The idea
//!
//! A pairwise SNP distance matrix is the raw material of a transmission study,
//! and on its own it is a wall of numbers. Two things turn it into a figure:
//! colour, so the near pairs separate from the far ones without anyone reading
//! a digit, and an order, so the near pairs land next to each other instead of
//! being scattered across the wall by an alphabet.
//!
//! The order is what the tree is for. Sorted by sample name a cluster is
//! confetti; sorted by descent it is a block on the diagonal, and a block on
//! the diagonal is a claim you can argue with.
//!
//! # What the x axis is
//!
//! Nothing genomic. Both axes of this track are the sample list, so it does not
//! use the figure's shared coordinate system and there is no sense in putting
//! an [`AxisTrack`](crate::AxisTrack) under it. Cells are kept square, because
//! the two axes are the same axis and drawing them at different scales would
//! say otherwise. A matrix small enough that square cells would be enormous
//! stops at [`DistanceTrack::cell_size`] and leaves the rest of the width
//! empty rather than stretching.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::tree::{draw_tree, leaf_order, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// How far off the page the palest cell sits, so that a zero is still a cell.
const ZERO_TINT: f64 = 0.08;

/// A pairwise distance matrix.
///
/// ```
/// use karyon::{DistanceTrack, Figure, Region};
///
/// // Three isolates, SNP distances, row major.
/// let names = vec!["ERR001".to_string(), "ERR002".to_string(), "ERR003".to_string()];
/// let distances = vec![
///     0.0, 3.0, 44.0, //
///     3.0, 0.0, 41.0, //
///     44.0, 41.0, 0.0,
/// ];
///
/// let svg = Figure::new(Region::new("samples", 0, 3).unwrap())
///     .push(DistanceTrack::new(names, distances).label("SNP distance"))
///     .to_svg();
/// assert!(svg.contains("ERR001"));
/// ```
#[derive(Debug, Clone)]
pub struct DistanceTrack {
    names: Vec<String>,
    values: Vec<f64>,
    label: Option<String>,
    max: Option<f64>,
    hue: Option<String>,
    missing_color: Option<String>,
    cell_size: f64,
    show_row_names: bool,
    show_column_names: bool,
    show_values: bool,
    cluster_threshold: Option<f64>,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
}

impl DistanceTrack {
    /// A matrix over `names`, with `values` holding the rows one after another.
    ///
    /// `values[i * n + j]` is the distance from sample `i` to sample `j`, so a
    /// square matrix written out row by row goes straight in. Nothing is
    /// symmetrised: if the two halves disagree the figure shows that they do,
    /// which is more use than quietly averaging them.
    pub fn new(names: impl Into<Vec<String>>, values: impl Into<Vec<f64>>) -> Self {
        DistanceTrack {
            names: names.into(),
            values: values.into(),
            label: None,
            max: None,
            hue: None,
            missing_color: None,
            cell_size: 34.0,
            show_row_names: true,
            show_column_names: true,
            show_values: true,
            cluster_threshold: None,
            tree: None,
            tree_width: 90.0,
            tree_shape: TreeShape::Phylogram,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Pins the distance that reaches full colour.
    ///
    /// Pin it when comparing two matrices, or the same shade will mean two
    /// different distances in two figures next to each other.
    pub fn max(mut self, max: f64) -> Self {
        self.max = Some(max);
        self
    }

    /// Sets the hue the ramp runs to, which defaults to the theme accent.
    pub fn hue(mut self, hue: impl Into<String>) -> Self {
        self.hue = Some(hue.into());
        self
    }

    /// Sets the colour of a cell with no distance in it.
    pub fn missing_color(mut self, color: impl Into<String>) -> Self {
        self.missing_color = Some(color.into());
        self
    }

    /// Sets the largest a cell is allowed to get, in pixels.
    ///
    /// Cells are square and as wide as the space divided by the sample count,
    /// up to this. Past it the matrix simply stops rather than stretching into
    /// a grid of letterboxes, since a distance matrix drawn wider than it is
    /// tall implies its two axes are different things.
    pub fn cell_size(mut self, pixels: f64) -> Self {
        self.cell_size = pixels.max(2.0);
        self
    }

    /// Draws or hides the sample names down the left.
    pub fn show_row_names(mut self, show: bool) -> Self {
        self.show_row_names = show;
        self
    }

    /// Draws or hides the sample names along the top.
    pub fn show_column_names(mut self, show: bool) -> Self {
        self.show_column_names = show;
        self
    }

    /// Draws or hides the distance written inside each cell.
    ///
    /// Only ever drawn when the cell is big enough to hold it. SNP distances
    /// are small integers and worth reading exactly; a matrix of substitutions
    /// per site is not, and is better left to the colour.
    pub fn show_values(mut self, show: bool) -> Self {
        self.show_values = show;
        self
    }

    /// Outlines every pair at or under `threshold`.
    ///
    /// The usual question asked of a SNP distance matrix is not "how far apart
    /// is everything" but "which pairs are close enough to be recent
    /// transmission", and that question comes with a number attached. Drawing
    /// the number's answer is more honest than leaving the reader to judge it
    /// off a colour ramp.
    ///
    /// The threshold is a convention, not a fact: it belongs to a pathogen, a
    /// mutation rate and a sampling window, so put the one you are using in the
    /// caption.
    pub fn cluster_threshold(mut self, threshold: f64) -> Self {
        self.cluster_threshold = Some(threshold);
        self
    }

    /// Draws a phylogeny beside the matrix and sorts both axes to match it.
    ///
    /// Rows and columns move together, because they are the same list. A sample
    /// the tree does not name keeps its place at the end rather than being
    /// dropped.
    pub fn tree(mut self, tree: Tree) -> Self {
        let order = leaf_order(&tree, &self.names);
        let size = self.names.len();
        let mut values = vec![f64::NAN; size * size];
        for (row, from) in order.iter().enumerate() {
            for (column, to) in order.iter().enumerate() {
                if let Some(value) = self.distance(*from, *to) {
                    values[row * size + column] = value;
                }
            }
        }
        self.values = values;
        self.names = order
            .iter()
            .map(|index| self.names[*index].clone())
            .collect();
        self.tree = Some(tree);
        self
    }

    /// Sets how much of the strip the tree gets, in pixels.
    pub fn tree_width(mut self, width: f64) -> Self {
        self.tree_width = width.max(0.0);
        self
    }

    /// Chooses a phylogram or a cladogram for the tree beside the matrix.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// The tree beside the matrix, if there is one.
    pub fn attached_tree(&self) -> Option<&Tree> {
        self.tree.as_ref()
    }

    /// The sample names, in the order the matrix is drawn.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// How many samples the matrix covers.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Whether the matrix has no samples in it.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// The distance between two samples, or `None` when the matrix is short.
    pub fn distance(&self, from: usize, to: usize) -> Option<f64> {
        let size = self.names.len();
        if from >= size || to >= size {
            return None;
        }
        self.values
            .get(from * size + to)
            .copied()
            .filter(|value| value.is_finite())
    }

    /// The largest distance in the matrix.
    pub fn ceiling(&self) -> f64 {
        if let Some(max) = self.max {
            return max;
        }
        let largest = self
            .values
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .fold(0.0f64, f64::max);
        if largest > 0.0 {
            largest
        } else {
            1.0
        }
    }

    /// Every pair at or under the clustering threshold, each once.
    ///
    /// Empty without a threshold, since there is then nothing to be inside of.
    /// A sample is not paired with itself.
    pub fn clustered_pairs(&self) -> Vec<(usize, usize)> {
        let Some(threshold) = self.cluster_threshold else {
            return Vec::new();
        };
        let mut pairs = Vec::new();
        for from in 0..self.names.len() {
            for to in (from + 1)..self.names.len() {
                if self.distance(from, to).is_some_and(|d| d <= threshold) {
                    pairs.push((from, to));
                }
            }
        }
        pairs
    }

    /// Side of one cell, in pixels.
    fn cell(&self, scale: &Scale) -> f64 {
        let count = self.names.len().max(1) as f64;
        (scale.width() / count).min(self.cell_size).max(1.0)
    }

    /// Height of the strip above the matrix holding the column names.
    fn header(&self, theme: &Theme) -> f64 {
        if !self.show_column_names || self.names.is_empty() {
            return 0.0;
        }
        let size = theme.font_size - 2.0;
        self.names
            .iter()
            .map(|name| text_width(name, size))
            .fold(0.0f64, f64::max)
            + 6.0
    }
}

impl Track for DistanceTrack {
    fn height(&self, scale: &Scale) -> f64 {
        // The theme reaches the height only through the width of the column
        // names, and the default font size is the one the figure will use
        // unless the caller changed it, as in the feature track.
        let rows = self.names.len().max(1) as f64;
        rows * self.cell(scale) + self.header(&Theme::default())
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let tree = if self.tree.is_some() {
            self.tree_width
        } else {
            0.0
        };
        if !self.show_row_names || self.names.is_empty() {
            return tree;
        }
        let size = theme.font_size - 2.0;
        let widest = self
            .names
            .iter()
            .map(|name| text_width(name, size))
            .fold(0.0f64, f64::max);
        widest + 8.0 + tree
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        if self.names.is_empty() {
            return;
        }
        let cell = self.cell(ctx.scale);
        let header = self.header(ctx.theme);
        let top = band.y + header;
        let ceiling = self.ceiling();
        let hue = self.hue.clone().unwrap_or_else(|| ctx.theme.accent.clone());
        let missing = self
            .missing_color
            .clone()
            .unwrap_or_else(|| mix(ctx.theme.surface(), &ctx.theme.rule, 0.45));
        let name_size = (ctx.theme.font_size - 2.0).min(cell);

        if let Some(tree) = &self.tree {
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: top,
                w: (self.tree_width - 6.0).max(1.0),
                h: cell * self.names.len() as f64,
            };
            draw_tree(
                ctx.svg,
                tree,
                area,
                cell,
                top + cell / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &ctx.theme.foreground,
                    width: 1.1,
                },
            );
        }

        for row in 0..self.names.len() {
            let y = top + row as f64 * cell;
            for column in 0..self.names.len() {
                let x = band.x + column as f64 * cell;
                let color = match self.distance(row, column) {
                    Some(value) => {
                        let fraction = (value / ceiling).clamp(0.0, 1.0);
                        // The ramp starts a shade off the page, so that a zero
                        // distance is still visibly a cell and not a hole.
                        mix(
                            ctx.theme.surface(),
                            &hue,
                            ZERO_TINT + (1.0 - ZERO_TINT) * fraction,
                        )
                    }
                    None => missing.clone(),
                };
                ctx.svg.rect(x, y, cell, cell, &color);

                if self.show_values {
                    if let Some(value) = self.distance(row, column) {
                        let text = trim(value);
                        let size = (cell * 0.42).min(ctx.theme.font_size);
                        if text_width(&text, size) < cell - 2.0 && size >= 5.0 {
                            ctx.svg.text(
                                x + cell / 2.0,
                                y + cell / 2.0 + size * 0.35,
                                &text,
                                contrast_ink(&color),
                                size,
                                Anchor::Middle,
                            );
                        }
                    }
                }
            }

            if self.show_row_names && ctx.axis.w > 0.0 {
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    y + cell / 2.0 + name_size * 0.35,
                    &self.names[row],
                    &ctx.theme.muted,
                    name_size,
                    Anchor::End,
                );
            }
        }

        // Outlines go on top of every cell, or a neighbour would paint over
        // half of each one.
        for (from, to) in self.clustered_pairs() {
            let ink = ctx.theme.foreground.clone();
            for (row, column) in [(from, to), (to, from)] {
                ctx.svg.rect_outline(
                    band.x + column as f64 * cell + 0.5,
                    top + row as f64 * cell + 0.5,
                    cell - 1.0,
                    cell - 1.0,
                    &ink,
                    1.2,
                );
            }
        }

        if self.show_column_names {
            let size = ctx.theme.font_size - 2.0;
            for (column, name) in self.names.iter().enumerate() {
                // Standing the name on end is what lets a sample accession sit
                // over a column narrower than the accession is long.
                ctx.svg.text_rotated(
                    (
                        band.x + column as f64 * cell + cell / 2.0 + size * 0.35,
                        top - 3.0,
                    ),
                    -90.0,
                    name,
                    &ctx.theme.muted,
                    size,
                    Anchor::Start,
                );
            }
        }
    }
}

/// Two decimals at most, and none when they would be zeros.
fn trim(value: f64) -> String {
    let rounded = (value * 100.0).round() / 100.0;
    if rounded == rounded.trunc() {
        format!("{}", rounded as i64)
    } else {
        format!("{rounded}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn names() -> Vec<String> {
        ["ERR001", "ERR002", "ERR003", "ERR004"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// Two pairs: 001 with 002 at three SNPs, 003 with 004 at four.
    fn matrix() -> DistanceTrack {
        DistanceTrack::new(
            names(),
            vec![
                0.0, 3.0, 40.0, 42.0, //
                3.0, 0.0, 41.0, 43.0, //
                40.0, 41.0, 0.0, 4.0, //
                42.0, 43.0, 4.0, 0.0,
            ],
        )
    }

    fn region() -> Region {
        Region::new("samples", 0, 4).unwrap()
    }

    fn scale() -> Scale {
        Scale::new(&region(), 0.0, 400.0)
    }

    #[test]
    fn a_distance_is_looked_up_by_its_two_samples() {
        let track = matrix();
        assert_eq!(track.distance(0, 1), Some(3.0));
        assert_eq!(track.distance(1, 0), Some(3.0));
        assert_eq!(track.distance(0, 0), Some(0.0));
        assert_eq!(track.distance(0, 9), None, "past the end of the matrix");
    }

    #[test]
    fn a_short_matrix_reads_as_missing_rather_than_as_zero() {
        let track = DistanceTrack::new(names(), vec![0.0, 1.0]);
        assert_eq!(track.distance(0, 1), Some(1.0));
        assert_eq!(track.distance(2, 2), None);
        assert!(!track
            .clustered_pairs()
            .iter()
            .any(|(a, b)| *a == 2 || *b == 2));
    }

    #[test]
    fn the_ceiling_is_the_largest_distance_unless_pinned() {
        assert_eq!(matrix().ceiling(), 43.0);
        assert_eq!(matrix().max(100.0).ceiling(), 100.0);
        // An all-zero matrix would otherwise divide by zero.
        let flat = DistanceTrack::new(vec!["a".to_string()], vec![0.0]);
        assert_eq!(flat.ceiling(), 1.0);
    }

    #[test]
    fn a_threshold_finds_each_close_pair_once() {
        let pairs = matrix().cluster_threshold(12.0).clustered_pairs();
        assert_eq!(pairs, vec![(0, 1), (2, 3)]);
        // And without one there is nothing to be inside of.
        assert!(matrix().clustered_pairs().is_empty());
    }

    #[test]
    fn a_sample_is_not_in_a_cluster_with_itself() {
        let pairs = matrix().cluster_threshold(0.0).clustered_pairs();
        assert!(pairs.is_empty(), "the diagonal is not transmission");
    }

    #[test]
    fn the_tree_moves_rows_and_columns_together() {
        // The tree pairs 001 with 003, against the order they came in.
        let tree =
            Tree::parse_newick("((ERR001:0.01,ERR003:0.01):0.02,(ERR002:0.01,ERR004:0.01):0.02);")
                .unwrap();
        let track = matrix().tree(tree);
        assert_eq!(track.names(), ["ERR001", "ERR003", "ERR002", "ERR004"]);
        // Whatever the order, the distance between two samples is unchanged.
        let of = |a: &str, b: &str| {
            let index = |name: &str| track.names().iter().position(|n| n == name).unwrap();
            track.distance(index(a), index(b))
        };
        assert_eq!(of("ERR001", "ERR002"), Some(3.0));
        assert_eq!(of("ERR003", "ERR004"), Some(4.0));
        assert_eq!(of("ERR001", "ERR004"), Some(42.0));
        // And the diagonal is still the diagonal.
        for row in 0..4 {
            assert_eq!(track.distance(row, row), Some(0.0));
        }
    }

    #[test]
    fn a_sample_the_tree_does_not_name_keeps_its_place() {
        let tree = Tree::parse_newick("(ERR003,ERR001);").unwrap();
        let track = matrix().tree(tree);
        assert_eq!(track.names(), ["ERR003", "ERR001", "ERR002", "ERR004"]);
        assert_eq!(track.len(), 4, "nothing was dropped");
    }

    #[test]
    fn cells_are_square_and_stop_growing() {
        let track = matrix();
        // Four samples across four hundred pixels is a hundred each, which is
        // more than a cell is allowed to be.
        assert_eq!(track.cell(&scale()), 34.0);
        assert_eq!(
            track.height(&scale()),
            4.0 * 34.0 + track.header(&Theme::default())
        );

        // Enough samples and the width is what limits them.
        let many: Vec<String> = (0..100).map(|i| format!("s{i}")).collect();
        let wide = DistanceTrack::new(many, vec![1.0; 100 * 100]);
        assert_eq!(wide.cell(&scale()), 4.0);
    }

    #[test]
    fn the_matrix_is_drawn_and_named() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(matrix().label("SNP distance"))
            .to_svg();
        assert!(svg.contains("ERR001"));
        assert!(svg.contains("SNP distance"));
        assert!(svg.contains("rotate(-90)"), "columns are named on end");
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn close_pairs_are_outlined_on_both_sides_of_the_diagonal() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(matrix().cluster_threshold(12.0).show_values(false))
            .to_svg();
        // Two pairs, each outlined twice because the matrix is symmetric.
        assert_eq!(svg.matches(r#"fill="none" stroke="#).count(), 4);
    }

    #[test]
    fn values_are_only_written_when_they_fit() {
        let roomy = Figure::new(region())
            .show_region_label(false)
            .push(matrix().show_column_names(false))
            .to_svg();
        assert!(roomy.contains(">40</text>"));

        let cramped = Figure::new(region())
            .show_region_label(false)
            .push(matrix().show_column_names(false).cell_size(4.0))
            .to_svg();
        assert!(
            !cramped.contains(">40</text>"),
            "a digit wider than its cell"
        );
    }

    #[test]
    fn an_empty_matrix_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(DistanceTrack::new(Vec::new(), Vec::new()).label("none"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
