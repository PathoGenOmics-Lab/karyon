//! A sample by site matrix laid along the sequence.
//!
//! This is the plot that answers "which samples carry it": one row per sample,
//! one column per site, and a cell saying what that sample had there. A
//! genotype matrix out of a VCF is exactly this shape, and so is a pangenome
//! presence and absence matrix once its genes have coordinates.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::tree::{draw_tree, leaf_order, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// How far off the page the bottom of a sequential ramp sits.
///
/// Enough to see a cell, not enough to read it as a value.
const ZERO_TINT: f64 = 0.1;

/// How a cell value becomes a colour.
#[derive(Debug, Clone, PartialEq)]
pub enum CellScale {
    /// One hue, pale for zero and saturated for the maximum.
    ///
    /// The right choice for a quantity, and for presence and absence, where the
    /// two ends of the ramp are the whole story. A sequential ramp is one hue
    /// light to dark: two hues would imply a meaningful middle, and zero to one
    /// has no middle to mean anything.
    Sequential {
        /// Value that reaches full saturation. `None` takes the largest in the
        /// matrix.
        max: Option<f64>,
        /// The hue, defaulting to the theme accent.
        hue: Option<String>,
    },
    /// The value is an index into the categorical palette, rounded.
    ///
    /// For genotypes and other labels, where the numbers name things rather
    /// than measure them and a ramp would imply an order they do not have.
    Categorical,
}

impl Default for CellScale {
    fn default() -> Self {
        CellScale::Sequential {
            max: None,
            hue: None,
        }
    }
}

/// One sample, and what it had at every site.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixRow {
    /// Sample name, drawn beside the row.
    pub name: String,
    /// One value per site, in the same order as the sites. A non-finite value
    /// is missing data and is drawn as such.
    pub values: Vec<f64>,
}

impl MatrixRow {
    /// A row named `name`.
    pub fn new(name: impl Into<String>, values: impl Into<Vec<f64>>) -> Self {
        MatrixRow {
            name: name.into(),
            values: values.into(),
        }
    }

    /// The value at column `index`, if the row is long enough to have one.
    pub fn value(&self, index: usize) -> Option<f64> {
        self.values.get(index).copied().filter(|v| v.is_finite())
    }
}

/// A matrix of samples against sites.
///
/// ```
/// use karyon::{Figure, MatrixRow, MatrixTrack, Region};
///
/// // Three isolates typed at four variant sites.
/// let sites = vec![1_100, 1_400, 2_300, 2_800];
/// let rows = vec![
///     MatrixRow::new("ERR001", vec![1.0, 0.0, 1.0, 0.0]),
///     MatrixRow::new("ERR002", vec![1.0, 1.0, 1.0, 0.0]),
///     MatrixRow::new("ERR003", vec![0.0, 0.0, f64::NAN, 1.0]),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 1_000, 3_000).unwrap())
///     .push(MatrixTrack::new(sites, rows).label("isolates"))
///     .to_svg();
/// assert!(svg.contains("ERR001"));
/// ```
#[derive(Debug, Clone)]
pub struct MatrixTrack {
    sites: Vec<u64>,
    rows: Vec<MatrixRow>,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    scale: CellScale,
    missing_color: Option<String>,
    min_cell_width: f64,
    show_row_names: bool,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
}

impl MatrixTrack {
    /// A matrix over `sites`, one row per sample.
    ///
    /// `sites[j]` is the 0-based position of column `j`, and `rows[i].values[j]`
    /// is what sample `i` had there. Rows shorter than the site list simply
    /// stop; nothing is invented to fill them.
    pub fn new(sites: impl Into<Vec<u64>>, rows: impl Into<Vec<MatrixRow>>) -> Self {
        MatrixTrack {
            sites: sites.into(),
            rows: rows.into(),
            label: None,
            row_height: 11.0,
            row_gap: 1.0,
            scale: CellScale::default(),
            missing_color: None,
            min_cell_width: 3.0,
            show_row_names: true,
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

    /// Sets the height of one row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(1.0);
        self
    }

    /// Sets the gap between rows.
    ///
    /// A gap in the page colour is what separates two cells of the same shade,
    /// rather than a border drawn round each one, which would add ink that is
    /// not data.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Chooses how a value becomes a colour.
    pub fn scale(mut self, scale: CellScale) -> Self {
        self.scale = scale;
        self
    }

    /// Sets the colour of a missing cell.
    pub fn missing_color(mut self, color: impl Into<String>) -> Self {
        self.missing_color = Some(color.into());
        self
    }

    /// Sets how wide a cell is drawn at its narrowest.
    ///
    /// Variant sites are points, so at any useful zoom a cell would be a
    /// fraction of a pixel. The floor is what makes the matrix a matrix rather
    /// than a row of invisible slivers, and it means a cell's width says
    /// nothing about how much sequence it covers.
    pub fn min_cell_width(mut self, pixels: f64) -> Self {
        self.min_cell_width = pixels.max(0.5);
        self
    }

    /// Draws or hides the sample names.
    pub fn show_row_names(mut self, show: bool) -> Self {
        self.show_row_names = show;
        self
    }

    /// Draws a phylogeny beside the matrix and sorts the rows to match it.
    ///
    /// For a pangenome this is the difference between a figure and a shrug.
    /// Sorted by sample name, an accessory region present in one clade and
    /// absent from another is a speckle; sorted by descent it is a rectangle,
    /// and a rectangle is a claim about the biology.
    ///
    /// Rows match leaves by name. A sample the tree does not name keeps its
    /// place at the bottom rather than vanishing, since a row silently dropped
    /// from a figure is worse than a row out of order.
    pub fn tree(mut self, tree: Tree) -> Self {
        let names: Vec<String> = self.rows.iter().map(|row| row.name.clone()).collect();
        let order = leaf_order(&tree, &names);
        self.rows = order
            .iter()
            .map(|index| self.rows[*index].clone())
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

    /// The site positions.
    pub fn sites(&self) -> &[u64] {
        &self.sites
    }

    /// The rows.
    pub fn rows(&self) -> &[MatrixRow] {
        &self.rows
    }

    /// Largest finite value anywhere in the matrix.
    pub fn value_ceiling(&self) -> Option<f64> {
        self.rows
            .iter()
            .flat_map(|row| row.values.iter())
            .copied()
            .filter(|v| v.is_finite())
            .fold(None, |acc: Option<f64>, v| {
                Some(acc.map_or(v, |a| a.max(v)))
            })
    }

    /// Colour of one cell.
    fn cell_color(&self, value: f64, ceiling: f64, theme: &Theme) -> String {
        match &self.scale {
            CellScale::Categorical => theme.color(value.max(0.0).round() as usize).to_string(),
            CellScale::Sequential { hue, .. } => {
                let hue = hue.clone().unwrap_or_else(|| theme.accent.clone());
                let fraction = if ceiling > 0.0 {
                    (value / ceiling).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                // The ramp starts a step off the page rather than on it. A pure
                // surface at the bottom would make an absent genotype
                // indistinguishable from no cell at all, and "this sample does
                // not have it" and "this sample was not typed" are different
                // statements that must not look the same.
                let t = ZERO_TINT + (1.0 - ZERO_TINT) * fraction;
                mix(theme.surface(), &hue, t)
            }
        }
    }
}

impl Track for MatrixTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.rows.len().max(1) as f64;
        rows * self.row_height + (rows - 1.0).max(0.0) * self.row_gap
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
        if !self.show_row_names || self.rows.is_empty() {
            return tree;
        }
        let size = (theme.font_size - 2.0).min(self.row_height);
        let widest = self
            .rows
            .iter()
            .map(|row| text_width(&row.name, size))
            .fold(0.0f64, f64::max);
        widest + 8.0 + tree
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let ceiling = match &self.scale {
            CellScale::Sequential { max: Some(max), .. } => *max,
            _ => self.value_ceiling().unwrap_or(1.0),
        };
        let missing = self
            .missing_color
            .clone()
            .unwrap_or_else(|| mix(ctx.theme.surface(), &ctx.theme.rule, 0.45));
        let name_size = (ctx.theme.font_size - 2.0).min(self.row_height);

        // The tree takes the left of the strip and the names the right of it,
        // so a leaf, its name and its row of cells all sit on one line.
        if let Some(tree) = &self.tree {
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: band.y,
                w: (self.tree_width - 6.0).max(1.0),
                h: band.h,
            };
            draw_tree(
                ctx.svg,
                tree,
                area,
                self.row_height + self.row_gap,
                band.y + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &ctx.theme.foreground,
                    width: 1.1,
                    mirror: false,
                },
            );
        }

        for (index, row) in self.rows.iter().enumerate() {
            let top = band.y + index as f64 * (self.row_height + self.row_gap);
            // The gap between rows is left in the page colour rather than drawn.
            let cell_height = self.row_height;

            for (column, site) in self.sites.iter().enumerate() {
                if !ctx.region.contains(*site) {
                    continue;
                }
                let x = ctx.scale.x(*site);
                let width = (ctx.scale.x(site + 1) - x).max(self.min_cell_width);
                let color = match row.value(column) {
                    Some(value) => self.cell_color(value, ceiling, ctx.theme),
                    None => missing.clone(),
                };
                ctx.svg.rect(x, top, width, cell_height, &color);
            }

            if self.show_row_names && ctx.axis.w > 0.0 {
                ctx.svg.text(
                    ctx.axis.right() - 4.0,
                    top + cell_height / 2.0 + name_size * 0.35,
                    &row.name,
                    &ctx.theme.muted,
                    name_size,
                    Anchor::End,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn region() -> Region {
        Region::new("chr1", 1_000, 3_000).unwrap()
    }

    fn matrix() -> MatrixTrack {
        MatrixTrack::new(
            vec![1_100u64, 1_400, 2_300, 2_800],
            vec![
                MatrixRow::new("ERR001", vec![1.0, 0.0, 1.0, 0.0]),
                MatrixRow::new("ERR002", vec![1.0, 1.0, 1.0, 0.0]),
                MatrixRow::new("ERR003", vec![0.0, 0.0, f64::NAN, 1.0]),
            ],
        )
    }

    #[test]
    fn a_missing_value_is_missing_rather_than_zero() {
        let row = MatrixRow::new("s", vec![1.0, f64::NAN, 0.0]);
        assert_eq!(row.value(0), Some(1.0));
        assert_eq!(row.value(1), None, "not there is not the same as zero");
        assert_eq!(row.value(2), Some(0.0));
        assert_eq!(row.value(9), None, "past the end of the row");
    }

    #[test]
    fn the_ceiling_ignores_missing_values() {
        let track = MatrixTrack::new(
            vec![1u64, 2],
            vec![MatrixRow::new("s", vec![f64::INFINITY, 3.0])],
        );
        assert_eq!(track.value_ceiling(), Some(3.0));
    }

    #[test]
    fn an_empty_matrix_has_no_ceiling_and_one_row_of_height() {
        let track = MatrixTrack::new(Vec::new(), Vec::new());
        assert_eq!(track.value_ceiling(), None);
        let scale = Scale::new(&region(), 0.0, 100.0);
        assert_eq!(track.height(&scale), 11.0);
    }

    #[test]
    fn height_follows_the_row_count() {
        let scale = Scale::new(&region(), 0.0, 100.0);
        assert_eq!(matrix().height(&scale), 3.0 * 11.0 + 2.0 * 1.0);
    }

    #[test]
    fn a_sequential_ramp_reaches_the_hue_but_starts_off_the_page() {
        let theme = Theme::light();
        let track = matrix().scale(CellScale::Sequential {
            max: Some(1.0),
            hue: Some("#000000".to_string()),
        });
        assert_eq!(track.cell_color(1.0, 1.0, &theme), "#000000");
        // A zero cell is visible as a cell: "does not have it" must not look
        // like "was not typed", which in turn must not look like blank page.
        let zero = track.cell_color(0.0, 1.0, &theme);
        assert_ne!(zero, theme.background);
        // But it is still nearly the page, not a value in its own right.
        assert_eq!(zero, "#e6e6e6");
        // Halfway is halfway, not one end or the other.
        let middle = track.cell_color(0.5, 1.0, &theme);
        assert_ne!(middle, zero);
        assert_ne!(middle, "#000000");
    }

    #[test]
    fn a_categorical_scale_indexes_the_palette() {
        let theme = Theme::light();
        let track = matrix().scale(CellScale::Categorical);
        assert_eq!(track.cell_color(0.0, 1.0, &theme), theme.color(0));
        assert_eq!(track.cell_color(2.0, 1.0, &theme), theme.color(2));
        // A value between two categories rounds to one of them.
        assert_eq!(track.cell_color(1.9, 1.0, &theme), theme.color(2));
        // And a negative one cannot index out of the palette.
        assert_eq!(track.cell_color(-5.0, 1.0, &theme), theme.color(0));
    }

    #[test]
    fn row_names_go_in_the_axis_strip_and_size_it() {
        let theme = Theme::light();
        let named = matrix();
        let anonymous = matrix().show_row_names(false);
        assert!(named.y_axis_width(&theme) > 0.0);
        assert_eq!(anonymous.y_axis_width(&theme), 0.0);

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(named)
            .to_svg();
        assert!(svg.contains("ERR001"));
        assert!(svg.contains("ERR003"));
    }

    #[test]
    fn a_wider_name_asks_for_a_wider_strip() {
        let theme = Theme::light();
        let short = MatrixTrack::new(vec![1u64], vec![MatrixRow::new("a", vec![1.0])]);
        let long = MatrixTrack::new(
            vec![1u64],
            vec![MatrixRow::new("a_very_long_sample_name", vec![1.0])],
        );
        assert!(long.y_axis_width(&theme) > short.y_axis_width(&theme));
    }

    #[test]
    fn sites_outside_the_region_are_not_drawn() {
        let track = MatrixTrack::new(
            vec![1_100u64, 90_000],
            vec![MatrixRow::new("s", vec![1.0, 1.0])],
        )
        .show_row_names(false);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        // One cell, the page background, and the clip path's own rect.
        assert_eq!(svg.matches("<rect").count(), 3);
    }

    #[test]
    fn a_cell_is_never_narrower_than_the_floor() {
        // Two kilobases across the figure puts a single base well under a pixel.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(matrix().show_row_names(false).min_cell_width(6.0))
            .to_svg();
        let widths: Vec<f64> = svg
            .match_indices(r#"<rect x="#)
            .filter_map(|(index, _)| {
                let element = &svg[index..index + svg[index..].find("/>")?];
                let start = element.find(r#"width=""#)? + 7;
                let rest = &element[start..];
                let width: f64 = rest[..rest.find('"')?].parse().ok()?;
                // Skip the page background, which is the figure's own width.
                (width < 500.0).then_some(width)
            })
            .collect();
        assert!(!widths.is_empty());
        assert!(widths.iter().all(|w| *w >= 6.0), "{widths:?}");
    }

    #[test]
    fn a_short_row_simply_stops() {
        let track = MatrixTrack::new(
            vec![1_100u64, 1_400, 2_300],
            vec![MatrixRow::new("s", vec![1.0])],
        )
        .show_row_names(false);
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(!svg.contains("NaN"));
        // Three cells are drawn, one with data and two missing, plus the page
        // background and the clip path's own rect.
        assert_eq!(svg.matches("<rect").count(), 5);
    }
}
