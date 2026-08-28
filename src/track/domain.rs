//! Protein domains and other interval architectures, one sequence per row.
//!
//! A [`DomainTrack`] keeps the coordinate axis honest: domain boundaries use
//! the same half-open positions as every genomic interval in the crate. Attach
//! a phylogeny with [`DomainTrack::tree`] to sort architectures by descent and
//! make gains, losses and rearrangements read as clade-level patterns.

use std::collections::BTreeMap;

use crate::scale::Scale;
use crate::svg::{fit_text, text_width, Anchor};
use crate::theme::{contrast_ink, mix, wash};
use crate::track::traits::Traits;
use crate::track::tree::{draw_tree, leaf_order, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// One annotated interval in a protein or transcript architecture.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainFeature {
    /// Zero-based half-open start.
    pub start: u64,
    /// Zero-based half-open end.
    pub end: u64,
    /// Visible domain, motif, exon or repeat name.
    pub label: Option<String>,
    /// Explicit fill colour, or `None` for the categorical theme palette.
    pub color: Option<String>,
}

impl DomainFeature {
    /// Creates an interval from `start` up to but excluding `end`.
    pub fn new(start: u64, end: u64) -> Self {
        DomainFeature {
            start,
            end,
            label: None,
            color: None,
        }
    }

    /// Names the interval.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets an explicit fill colour.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }
}

/// Domain architecture for one named sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainArchitecture {
    /// Sequence name, matched exactly to an attached tree leaf.
    pub name: String,
    /// Sequence length in the track's coordinate units.
    pub length: u64,
    /// Domains, motifs, exons or other intervals in input order.
    pub features: Vec<DomainFeature>,
}

impl DomainArchitecture {
    /// Creates an empty architecture for a named sequence of `length` units.
    pub fn new(name: impl Into<String>, length: u64) -> Self {
        DomainArchitecture {
            name: name.into(),
            length,
            features: Vec::new(),
        }
    }

    /// Appends one interval feature.
    pub fn feature(mut self, feature: DomainFeature) -> Self {
        self.features.push(feature);
        self
    }
}

/// Aligned sequence architectures, optionally ordered by a phylogeny.
///
/// ```
/// use karyon::{DomainArchitecture, DomainFeature, DomainTrack, Figure, Region};
///
/// let rows = vec![
///     DomainArchitecture::new("A", 300)
///         .feature(DomainFeature::new(20, 110).label("kinase")),
///     DomainArchitecture::new("B", 300)
///         .feature(DomainFeature::new(160, 240).label("DNA-binding")),
/// ];
/// let svg = Figure::new(Region::new("protein", 0, 300).unwrap())
///     .push(DomainTrack::new(rows).label("domains"))
///     .to_svg();
/// assert!(svg.contains("kinase"));
/// ```
#[derive(Debug, Clone)]
pub struct DomainTrack {
    rows: Vec<DomainArchitecture>,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    show_names: bool,
    show_labels: bool,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
    traits: Traits,
}

impl DomainTrack {
    /// Creates a row-aligned architecture track.
    pub fn new(rows: impl Into<Vec<DomainArchitecture>>) -> Self {
        DomainTrack {
            rows: rows.into(),
            label: None,
            row_height: 14.0,
            row_gap: 3.0,
            show_names: true,
            show_labels: true,
            tree: None,
            tree_width: 100.0,
            tree_shape: TreeShape::Phylogram,
            traits: Traits::default(),
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one architecture row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = if height.is_finite() {
            height.max(3.0)
        } else {
            14.0
        };
        self
    }

    /// Sets the page-coloured gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = if gap.is_finite() { gap.max(0.0) } else { 3.0 };
        self
    }

    /// Draws or hides sequence names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Draws or hides domain names when an interval is wide enough.
    pub fn show_labels(mut self, show: bool) -> Self {
        self.show_labels = show;
        self
    }

    /// Draws a phylogeny beside the architectures and sorts rows by descent.
    ///
    /// Rows absent from the tree remain at the bottom instead of disappearing.
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

    /// Sets how much of the row-name strip the attached tree receives.
    pub fn tree_width(mut self, width: f64) -> Self {
        self.tree_width = if width.is_finite() {
            width.max(0.0)
        } else {
            100.0
        };
        self
    }

    /// Chooses a phylogram or cladogram for the attached tree.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// Attaches metadata columns drawn between the row names and the rows.
    ///
    /// Joined by name, so the columns follow whatever order the rows are
    /// already in.
    pub fn traits(mut self, traits: Traits) -> Self {
        self.traits = traits;
        self
    }

    /// The metadata columns attached to this track.
    pub fn attached_traits(&self) -> &Traits {
        &self.traits
    }

    /// Architectures in their current visual order.
    pub fn rows(&self) -> &[DomainArchitecture] {
        &self.rows
    }

    /// The tree attached to the architectures, when present.
    pub fn attached_tree(&self) -> Option<&Tree> {
        self.tree.as_ref()
    }

    fn categories(&self) -> BTreeMap<String, usize> {
        let mut categories = BTreeMap::new();
        for feature in self.rows.iter().flat_map(|row| &row.features) {
            let key = feature.label.clone().unwrap_or_default();
            let next = categories.len();
            categories.entry(key).or_insert(next);
        }
        categories
    }
}

impl Track for DomainTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.rows.len().max(1) as f64;
        rows * self.row_height + (rows - 1.0).max(0.0) * self.row_gap + self.traits.heading_height()
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &crate::Theme) -> f64 {
        let tree = self.tree.as_ref().map_or(0.0, |_| self.tree_width);
        let strip = self.traits.strip_width();
        if !self.show_names || self.rows.is_empty() {
            return tree + strip;
        }
        let size = (theme.font_size - 2.0).min(self.row_height);
        tree + strip
            + self
                .rows
                .iter()
                .map(|row| text_width(&row.name, size))
                .fold(0.0f64, f64::max)
            + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        // The band the rows are drawn in starts below the room the headings
        // of the metadata columns took, which `height` already asked for. The
        // bottom is unmoved, so everything measured from it still lands.
        let head = ctx.px(self.traits.heading_height());
        let strip = self.traits.strip_width();
        let band = Rect {
            y: ctx.band.y + head,
            h: (ctx.band.h - head).max(1.0),
            ..ctx.band
        };
        let categories = self.categories();
        let name_size = (ctx.theme.font_size - 2.0).min(self.row_height);
        let tree_color = mix(&ctx.theme.foreground, &ctx.theme.muted, 0.45);
        let backbone = mix(&ctx.theme.rule, &ctx.theme.foreground, 0.12);
        let zebra = mix(ctx.theme.surface(), &ctx.theme.rule, 0.16);

        if let Some(tree) = &self.tree {
            draw_tree(
                ctx.svg,
                tree,
                Rect {
                    x: ctx.axis.x + 2.0,
                    y: band.y,
                    w: (self.tree_width - 6.0).max(1.0),
                    h: band.h,
                },
                self.row_height + self.row_gap,
                band.y + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &tree_color,
                    width: ctx.theme.tokens.stroke.max(1.1),
                    mirror: false,
                },
            );
            ctx.svg.line(
                band.x - 2.0,
                band.y + 1.0,
                band.x - 2.0,
                band.bottom() - 1.0,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline.max(0.8),
            );
        }

        for (row_index, row) in self.rows.iter().enumerate() {
            let top = band.y + row_index as f64 * (self.row_height + self.row_gap);
            let centre = top + self.row_height / 2.0;
            if row_index % 2 == 1 {
                ctx.svg.rect(band.x, top, band.w, self.row_height, &zebra);
            }
            let visible_end = row.length.min(ctx.region.end());
            if visible_end > ctx.region.start() {
                ctx.svg.line(
                    ctx.scale.x(ctx.region.start()),
                    centre,
                    ctx.scale.x(visible_end),
                    centre,
                    &backbone,
                    ctx.theme.tokens.stroke.max(1.1),
                );
            }

            for feature in &row.features {
                let start = feature.start.max(ctx.region.start());
                let end = feature.end.min(ctx.region.end()).min(row.length);
                if end <= start {
                    continue;
                }
                let x = ctx.scale.x(start);
                let width = (ctx.scale.x(end) - x).max(0.8);
                let key = feature.label.clone().unwrap_or_default();
                let palette = *categories.get(&key).unwrap_or(&0);
                let color = feature
                    .color
                    .clone()
                    .unwrap_or_else(|| ctx.theme.color(palette).to_string());
                let face = if width >= 6.0 {
                    wash(&color, ctx.theme)
                } else {
                    color.clone()
                };
                let label = feature.label.as_deref().unwrap_or("domain");
                ctx.svg.begin_titled(&format!(
                    "{}; {label}; {} to {}",
                    row.name,
                    feature.start + 1,
                    feature.end
                ));
                ctx.svg.rect_rounded(
                    x,
                    top + 1.0,
                    width,
                    (self.row_height - 2.0).max(1.0),
                    ctx.theme.corner_radius.min(3.0),
                    &face,
                );
                if width >= 6.0 && self.row_height >= 7.0 {
                    ctx.svg.rect_rounded(
                        x + 1.5,
                        top + 3.0,
                        2.5,
                        (self.row_height - 6.0).max(1.0),
                        1.25,
                        &color,
                    );
                }
                if self.show_labels {
                    let size = (name_size - 1.0).max(6.0);
                    let visible = fit_text(label, width - 4.0, size);
                    if !visible.is_empty() {
                        ctx.svg.text(
                            x + width / 2.0,
                            centre + size * 0.34,
                            &visible,
                            contrast_ink(&face),
                            size,
                            Anchor::Middle,
                        );
                    }
                }
                ctx.svg.end_group();
            }

            if self.show_names && ctx.axis.w > 0.0 {
                ctx.svg.text(
                    ctx.axis.right() - strip - 4.0,
                    centre + name_size * 0.35,
                    &row.name,
                    &ctx.theme.foreground,
                    name_size,
                    Anchor::End,
                );
            }
        }

        let placed: Vec<(String, f64, f64)> = self
            .rows
            .iter()
            .enumerate()
            .map(|(index, row)| {
                (
                    row.name.clone(),
                    band.y + index as f64 * (self.row_height + self.row_gap),
                    self.row_height,
                )
            })
            .collect();
        self.traits.draw(
            ctx,
            Rect {
                x: ctx.axis.right() - strip,
                y: ctx.band.y,
                w: strip,
                h: ctx.band.h,
            },
            &placed,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Figure, Region};

    fn rows() -> Vec<DomainArchitecture> {
        vec![
            DomainArchitecture::new("B", 200).feature(DomainFeature::new(20, 80).label("ATPase")),
            DomainArchitecture::new("A", 200).feature(DomainFeature::new(90, 170).label("sensor")),
        ]
    }

    #[test]
    fn an_attached_tree_reorders_architectures() {
        let tree = Tree::parse_newick("(A:1,B:1);").unwrap();
        let track = DomainTrack::new(rows()).tree(tree);
        assert_eq!(track.rows()[0].name, "A");
        assert_eq!(track.rows()[1].name, "B");
    }

    #[test]
    fn labels_and_exact_boundaries_survive_in_tooltips() {
        let svg = Figure::new(Region::new("protein", 0, 200).unwrap())
            .show_region_label(false)
            .push(DomainTrack::new(rows()))
            .to_svg();
        assert!(svg.contains("ATPase"));
        assert!(svg.contains("21 to 80"));
    }
}
