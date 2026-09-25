//! Metadata columns: what is known about a row, drawn beside the row.
//!
//! A track drawn as rows answers "which ones". The question a reader asks next
//! is almost always "and what were they": which lineage, which host, which
//! treatment arm, which year. That answer is in a sample sheet rather than in
//! the file the track was drawn from, and it is what a [`TraitColumn`] draws:
//! one narrow strip per attribute, one cell per row, beside the rows it
//! describes.
//!
//! # Why this is not a track
//!
//! Every track in the crate is drawn on the shared coordinate axis. It is
//! placed at positions, and it moves when the region moves. An attribute has no
//! position at all: a sample's lineage is not at a base, and there is no zoom
//! level at which more of it comes into view. Given a track of its own it would
//! need an x it does not have, and the first pan would slide a sample's lineage
//! off the end of that sample's own row.
//!
//! So a column is not a track. It is drawn in the strip a track already
//! reserves to the left of the plotting area, the one the row names and the
//! dendrogram share, and it survives every pan and zoom untouched because
//! nothing in it was ever placed at a coordinate.
//!
//! # One vocabulary, wherever the rows came from
//!
//! These columns began beside a phylogeny, whose values came out of an
//! annotated Newick, and the tree still reads them from there. A matrix, an
//! alignment or a set of loci have rows too, and what is known about those rows
//! arrives as a table instead. Both end up as [`Annotations`] against a name,
//! so both are drawn by the code below rather than by two implementations that
//! would drift: the same lineage gets the same colour in a tree and in the
//! matrix beneath it, which is the whole reason to put them in one figure.
//!
//! # Colour is assigned by first appearance
//!
//! A column numbers each distinct value as it meets it, so a figure redrawn
//! from the same file colours the same way. Sorting the values instead
//! would recolour half a figure when a sample whose name sorts early is added,
//! and a figure that recolours itself cannot go in a paper.
//!
//! The palette has six colours. A column with more levels than that reuses one,
//! and two levels sharing a swatch is a figure that states something false, so
//! [`Traits::spread`] gives such a column [`TraitStyle::Symbol`], which carries
//! the level in a shape as well as a hue and separates twenty-four.
//!
//! # A missing value is drawn as missing
//!
//! Not as a colour, and not as a zero. A cell whose row says nothing about this
//! column is an empty outline, which is the one mark here that cannot be
//! mistaken for a level, and its tooltip says the word.

use std::collections::{BTreeMap, BTreeSet};

use crate::read::sheet::Sheet;
use crate::svg::{fit_text, fit_text_shrinking, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::legend::Legend;
use crate::track::{DrawContext, Rect};
use crate::tree::{AnnotationValue, Annotations};

/// How a trait column maps values to colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitScale {
    /// Each distinct value receives a categorical palette colour.
    Categorical,
    /// Numeric values form one continuous muted-to-accent ramp.
    Continuous,
}

/// Mark used for one metadata dataset beside or around the rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TraitStyle {
    /// One filled cell or annular sector per row.
    #[default]
    Strip,
    /// Numeric value encoded by bar length or radial height.
    Bar,
    /// Boolean or zero/non-zero value encoded by presence of a marker.
    Binary,
    /// Category encoded redundantly by both colour and marker shape.
    Symbol,
}

/// One metadata column drawn beside the rows of a track.
#[derive(Debug, Clone, PartialEq)]
pub struct TraitColumn {
    pub(crate) key: String,
    pub(crate) label: String,
    pub(crate) scale: TraitScale,
    pub(crate) style: TraitStyle,
    pub(crate) width: f64,
    pub(crate) ring_width: f64,
    pub(crate) show_values: bool,
    pub(crate) levels: Vec<String>,
}

impl TraitColumn {
    /// Builds a categorical column from annotation `key`.
    pub fn categorical(key: impl Into<String>) -> Self {
        let key = key.into();
        TraitColumn {
            label: key.clone(),
            key,
            scale: TraitScale::Categorical,
            style: TraitStyle::Strip,
            width: 56.0,
            ring_width: 10.0,
            show_values: true,
            levels: Vec::new(),
        }
    }

    /// Sets the order the levels are dealt their colours in: the first level
    /// named takes the palette's first colour, and a level met that is not
    /// named here takes the next one free.
    ///
    /// A column made by [`Traits::spread`] or added with [`Traits::column`]
    /// already carries the order its sheet lists the levels in. Handing that
    /// column to a phylogeny as well, with
    /// [`TreeTrack::trait_column`](crate::track::tree::TreeTrack::trait_column),
    /// is what makes a level one colour in both strips: the tree meets its
    /// tips in its own order, and without this each track dealt the palette
    /// in the order it met the levels.
    pub fn levels(mut self, order: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.levels = order.into_iter().map(Into::into).collect();
        self
    }

    /// The order the levels are dealt their colours in, if one was set.
    pub fn level_order(&self) -> &[String] {
        &self.levels
    }

    /// Builds a continuous column from numeric annotation `key`.
    pub fn continuous(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.scale = TraitScale::Continuous;
        column
    }

    /// Builds a numeric bar column or radial bar ring.
    pub fn bar(key: impl Into<String>) -> Self {
        let mut column = Self::continuous(key);
        column.style = TraitStyle::Bar;
        column.show_values = false;
        column
    }

    /// Builds a boolean presence/absence marker dataset.
    ///
    /// Boolean values and finite numbers are accepted; zero is absent and a
    /// non-zero number is present. Text is left missing rather than guessed.
    pub fn binary(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.style = TraitStyle::Binary;
        column.width = 28.0;
        column.show_values = false;
        column
    }

    /// Builds a categorical dataset encoded by colour and marker shape.
    pub fn symbol(key: impl Into<String>) -> Self {
        let mut column = Self::categorical(key);
        column.style = TraitStyle::Symbol;
        column.width = 32.0;
        column.show_values = false;
        column
    }

    /// Replaces the visible column heading without changing its metadata key.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Sets the cell width in pixels.
    pub fn width(mut self, width: f64) -> Self {
        self.width = if width.is_finite() {
            width.max(12.0)
        } else {
            56.0
        };
        self
    }

    /// Sets the thickness of this trait when drawn as a circular ring.
    pub fn ring_width(mut self, width: f64) -> Self {
        self.ring_width = if width.is_finite() {
            width.clamp(2.0, 24.0)
        } else {
            10.0
        };
        self
    }

    /// Draws or hides the value text inside each cell.
    pub fn show_values(mut self, show: bool) -> Self {
        self.show_values = show;
        self
    }

    /// Replaces the visual mark while retaining the column's value mapping.
    pub fn style(mut self, style: TraitStyle) -> Self {
        self.style = style;
        self
    }

    /// The annotation key read from each row.
    pub fn key(&self) -> &str {
        &self.key
    }

    /// The heading drawn over the column.
    pub fn heading(&self) -> &str {
        &self.label
    }

    /// The colour mapping used by this column.
    pub fn scale(&self) -> TraitScale {
        self.scale
    }

    /// The mark used in rectangular, radial and unrooted projections.
    pub fn trait_style(&self) -> TraitStyle {
        self.style
    }

    /// The cell width in pixels.
    pub fn cell_width(&self) -> f64 {
        self.width
    }
}

/// The levels and the range one column's values cover.
///
/// Built once per column from every value in it, because a colour is a
/// statement about a value's place among the others and cannot be worked out
/// from the value alone.
pub(crate) struct TraitDomain {
    pub(crate) categories: BTreeMap<String, usize>,
    /// The levels some value actually held, which a key names; the ones an
    /// order reserved and nobody here holds keep their colour and stay out of
    /// the key.
    met: BTreeSet<String>,
    pub(crate) minimum: f64,
    pub(crate) maximum: f64,
}

impl TraitDomain {
    pub(crate) fn new<'a>(values: impl IntoIterator<Item = &'a AnnotationValue>) -> Self {
        Self::ordered(&[], values)
    }

    /// The same, with `levels` dealt the palette first, in that order.
    pub(crate) fn ordered<'a>(
        levels: &[String],
        values: impl IntoIterator<Item = &'a AnnotationValue>,
    ) -> Self {
        let values: Vec<&AnnotationValue> = values.into_iter().collect();
        let mut categories = BTreeMap::new();
        for level in levels {
            let next = categories.len();
            categories.entry(level.clone()).or_insert(next);
        }
        let mut met = BTreeSet::new();
        for value in &values {
            let value = value.to_string();
            let next = categories.len();
            categories.entry(value.clone()).or_insert(next);
            met.insert(value);
        }
        let numeric: Vec<f64> = values
            .iter()
            .filter_map(|value| value.as_number())
            .filter(|value| value.is_finite())
            .collect();
        TraitDomain {
            categories,
            met,
            minimum: numeric.iter().copied().fold(f64::MAX, f64::min),
            maximum: numeric.iter().copied().fold(f64::MIN, f64::max),
        }
    }

    pub(crate) fn fraction(&self, value: Option<&AnnotationValue>) -> Option<f64> {
        let value = value?.as_number()?;
        if !value.is_finite() {
            return None;
        }
        Some(if self.maximum <= self.minimum {
            1.0
        } else {
            ((value - self.minimum) / (self.maximum - self.minimum)).clamp(0.0, 1.0)
        })
    }

    pub(crate) fn category(&self, value: Option<&AnnotationValue>) -> Option<usize> {
        self.categories.get(&value?.to_string()).copied()
    }

    pub(crate) fn color(
        &self,
        column: &TraitColumn,
        value: Option<&AnnotationValue>,
        theme: &Theme,
    ) -> Option<String> {
        match column.scale {
            TraitScale::Categorical => self
                .category(value)
                .map(|index| theme.color(index).to_string()),
            TraitScale::Continuous => self
                .fraction(value)
                .map(|fraction| mix(&theme.muted, &theme.accent, fraction)),
        }
    }

    /// The levels in the order their colours were assigned.
    ///
    /// The map is keyed by the text so that a lookup is a lookup, which puts
    /// its entries in the order the words sort. A legend has to name them in
    /// the order the palette went round instead, or the key and the strips
    /// disagree about which blue is which.
    pub(crate) fn levels(&self) -> Vec<(&str, usize)> {
        let mut levels: Vec<(&str, usize)> = self
            .categories
            .iter()
            .map(|(name, index)| (name.as_str(), *index))
            .collect();
        levels.sort_by_key(|(_, index)| *index);
        levels
    }

    /// The levels a key names: the ones some value held, in palette order.
    pub(crate) fn keyed(&self) -> Vec<(&str, usize)> {
        let mut levels = self.levels();
        levels.retain(|(level, _)| self.met.contains(*level));
        levels
    }
}

pub(crate) fn binary_state(value: Option<&AnnotationValue>) -> Option<bool> {
    match value? {
        AnnotationValue::Boolean(value) => Some(*value),
        AnnotationValue::Number(value) if value.is_finite() => Some(*value != 0.0),
        _ => None,
    }
}

/// One row a column is drawn against.
pub(crate) struct TraitRow<'a> {
    /// The row's own name, which is what the tooltip leads with.
    pub(crate) name: &'a str,
    /// Top of the cell in figure coordinates.
    pub(crate) top: f64,
    /// Height of the cell.
    pub(crate) height: f64,
    /// What this row holds for this column, or `None` where it holds nothing.
    pub(crate) value: Option<&'a AnnotationValue>,
}

/// Draws one column's heading and cells.
///
/// The one place a trait cell is drawn. A tree hands it terminal taxa and a
/// row-based track hands it rows; everything about how a level looks is
/// decided here, so the two cannot come out different.
pub(crate) fn draw_column(
    ctx: &mut DrawContext<'_>,
    column: &TraitColumn,
    domain: &TraitDomain,
    x: f64,
    heading_y: Option<f64>,
    rows: &[TraitRow<'_>],
) {
    let size = (ctx.theme.font_size - 2.0).max(6.0);

    if let Some(y) = heading_y {
        let heading = fit_text(&column.label, column.width, size);
        ctx.svg.text(
            x + column.width / 2.0,
            y,
            &heading,
            &ctx.theme.muted,
            size,
            Anchor::Middle,
        );
    }

    for row in rows {
        let TraitRow {
            name,
            top: y,
            height,
            value,
        } = *row;
        let fill = domain.color(column, value, ctx.theme);
        let displayed = value.map(ToString::to_string);
        let title = match &displayed {
            Some(value) => format!("{name}; {} {value}", column.key),
            None => format!("{name}; {} missing", column.key),
        };
        ctx.svg.begin_titled(&title);
        match column.style {
            TraitStyle::Strip => {
                if let Some(fill) = &fill {
                    ctx.svg.rect_rounded(
                        x,
                        y,
                        column.width,
                        height,
                        ctx.theme.corner_radius.min(2.0),
                        fill,
                    );
                } else {
                    ctx.svg.rect_outline(
                        x,
                        y,
                        column.width,
                        height,
                        &ctx.theme.rule,
                        ctx.theme.tokens.hairline,
                    );
                }
            }
            TraitStyle::Bar => {
                ctx.svg.rect_outline(
                    x,
                    y,
                    column.width,
                    height,
                    &ctx.theme.rule,
                    ctx.theme.tokens.hairline,
                );
                if let Some(fraction) = domain.fraction(value) {
                    ctx.svg.rect_rounded(
                        x,
                        y,
                        column.width * fraction,
                        height,
                        ctx.theme.corner_radius.min(2.0),
                        fill.as_deref().unwrap_or(&ctx.theme.accent),
                    );
                }
            }
            TraitStyle::Binary => match binary_state(value) {
                Some(true) => ctx.svg.circle_ringed(
                    x + column.width / 2.0,
                    y + height / 2.0,
                    (height * 0.28).clamp(1.4, 5.0),
                    &ctx.theme.accent,
                    &ctx.theme.background,
                    ctx.theme.tokens.hairline,
                ),
                Some(false) => ctx.svg.circle_ringed(
                    x + column.width / 2.0,
                    y + height / 2.0,
                    (height * 0.12).clamp(0.8, 2.0),
                    &ctx.theme.rule,
                    &ctx.theme.background,
                    ctx.theme.tokens.hairline,
                ),
                None => ctx.svg.rect_outline(
                    x,
                    y,
                    column.width,
                    height,
                    &ctx.theme.rule,
                    ctx.theme.tokens.hairline,
                ),
            },
            TraitStyle::Symbol => {
                if let Some(index) = domain.category(value) {
                    ctx.svg.symbol_ringed(
                        x + column.width / 2.0,
                        y + height / 2.0,
                        (height * 0.28).clamp(1.4, 5.0),
                        ctx.theme.symbol(index),
                        fill.as_deref().unwrap_or(&ctx.theme.accent),
                        &ctx.theme.background,
                        ctx.theme.tokens.hairline,
                    );
                } else {
                    ctx.svg.rect_outline(
                        x,
                        y,
                        column.width,
                        height,
                        &ctx.theme.rule,
                        ctx.theme.tokens.hairline,
                    );
                }
            }
        }
        if column.show_values && matches!(column.style, TraitStyle::Strip | TraitStyle::Bar) {
            let text = displayed.as_deref().unwrap_or(crate::tree::ABSENT);
            let (visible, size) = fit_text_shrinking(text, column.width - 4.0, size, size * 0.8);
            let ink = fill
                .as_deref()
                .filter(|_| column.style == TraitStyle::Strip)
                .map(contrast_ink)
                .unwrap_or(ctx.theme.muted.as_str());
            ctx.svg.text(
                x + column.width / 2.0,
                y + height / 2.0 + size * 0.35,
                &visible,
                ink,
                size,
                Anchor::Middle,
            );
        }
        ctx.svg.end_group();
    }
}

/// The number of levels beyond which a filled strip stops separating them.
///
/// The shipped palette has six colours. A theme may carry more, and a column
/// of seven levels then has seven distinct swatches, but the choice of mark is
/// made once when the columns are built and a theme arrives later, so it is
/// made against the palette the crate ships rather than against one it might
/// be handed.
const STRIP_LEVELS: usize = 6;

/// What is known about a track's rows, and the columns drawn from it.
///
/// The rows are keyed by name because that is the only thing a sheet and a
/// track share. A phylogeny attached to the track will already have put its
/// rows in the tree's order, and a value found by position would then be drawn
/// against the wrong sample and look exactly as convincing as the right one.
///
/// ```
/// use karyon::read;
/// use karyon::track::traits::Traits;
/// use karyon::{plot, MatrixRow, MatrixTrack};
///
/// let sheet = read::sheet::sheet(
///     "sample\tlineage\thost\tdepth\n\
///      S1\tL4\thuman\t72.5\n\
///      S2\tL2\tbovine\t61\n\
///      S3\tL4\t\t48.2\n",
/// )?;
/// let columns = sheet.columns.clone();
/// let traits = Traits::new(sheet.rows).spread(columns);
///
/// let rows = vec![
///     MatrixRow::new("S1", vec![1.0, 0.0]),
///     MatrixRow::new("S2", vec![0.0, 1.0]),
///     MatrixRow::new("S3", vec![1.0, 1.0]),
/// ];
/// let svg = plot("chr1:1-1,000")?
///     .add_track(MatrixTrack::new(vec![120, 340], rows).traits(traits))
///     .to_svg();
///
/// assert!(svg.contains("S1; lineage L4"));
/// // S3 has no host, and the figure says so rather than colouring it.
/// assert!(svg.contains("S3; host missing"));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Traits {
    rows: BTreeMap<String, Annotations>,
    /// The rows in the order their levels are met, which deals the palette.
    order: Vec<String>,
    columns: Vec<TraitColumn>,
    heading_room: f64,
    gap: f64,
}

impl Traits {
    /// Starts from what is known about each named row, with no columns yet.
    ///
    /// The levels of a column are dealt their colours in the order the names
    /// sort, which is the only order a map of rows has.
    /// [`Traits::from_sheet`] keeps the order of the file instead.
    pub fn new(rows: BTreeMap<String, Annotations>) -> Self {
        let order = rows.keys().cloned().collect();
        Traits {
            rows,
            order,
            columns: Vec::new(),
            heading_room: 52.0,
            gap: 2.0,
        }
    }

    /// Starts from a sample sheet, dealing each column's levels their colours
    /// in the order the file first gives them.
    ///
    /// So a sample appended to the end of the file never repaints the ones
    /// above it, whatever it is called, which sorting the names cannot
    /// promise: a new sample called `AAA` would have been met first.
    pub fn from_sheet(sheet: &Sheet) -> Self {
        let mut traits = Traits::new(sheet.rows.clone());
        traits.order = sheet.order.clone();
        traits
    }

    /// Adds one column drawn exactly as it was built.
    ///
    /// A categorical column given no [`TraitColumn::levels`] of its own is
    /// given the order these rows meet its levels in, so the column can be
    /// handed to a phylogeny too and colour each level the same there.
    pub fn column(mut self, column: TraitColumn) -> Self {
        let column = self.ordered(column);
        self.columns.push(column);
        self
    }

    /// A column carrying the order its levels are met in here, unless it
    /// already carries one.
    fn ordered(&self, column: TraitColumn) -> TraitColumn {
        if column.scale != TraitScale::Categorical || !column.levels.is_empty() {
            return column;
        }
        let met = self
            .domain(&column)
            .levels()
            .into_iter()
            .map(|(level, _)| level.to_string())
            .collect::<Vec<_>>();
        column.levels(met)
    }

    /// Every stated value of `key`, in the order the rows are met.
    fn stated<'a>(&'a self, key: &'a str) -> impl Iterator<Item = &'a AnnotationValue> + 'a {
        self.order
            .iter()
            .filter_map(|name| self.rows.get(name))
            .filter_map(move |held| held.get(key))
    }

    /// Adds a strip for each key, taking the mark from what the values are.
    ///
    /// A key whose every stated value is a number gets a ramp, because numbers
    /// with a ramp read as an order and numbers with a palette do not. Anything
    /// else gets a palette, and a palette of more than
    #[doc = concat!(stringify!(6), " levels")]
    /// gets [`TraitStyle::Symbol`] instead of a filled cell, so that the shape
    /// keeps two levels apart where the hue has run out and come round again.
    ///
    /// Keys are drawn in the order given, and a key no row mentions still gets
    /// a column: an attribute nobody in this figure has is a fact about the
    /// figure, and a column of empty outlines states it.
    pub fn spread(mut self, keys: impl IntoIterator<Item = impl Into<String>>) -> Self {
        for key in keys {
            let key = key.into();
            let stated: Vec<&AnnotationValue> = self.stated(&key).collect();
            let numeric =
                !stated.is_empty() && stated.iter().all(|value| value.as_number().is_some());

            let column = if numeric {
                TraitColumn::continuous(key)
            } else {
                let levels = TraitDomain::new(stated).categories.len();
                let column = TraitColumn::categorical(key);
                if levels > STRIP_LEVELS {
                    column.style(TraitStyle::Symbol)
                } else {
                    column
                }
            };
            let column = self.ordered(column.width(14.0).show_values(false));
            self.columns.push(column);
        }
        self
    }

    /// Sets the air between one column and the next, in pixels.
    ///
    /// Small on purpose. Columns further apart than they are wide stop being a
    /// strip and become three separate figures, and the thing a reader is
    /// looking for here is a block of one colour running down several rows.
    pub fn gap(mut self, gap: f64) -> Self {
        self.gap = if gap.is_finite() {
            gap.clamp(0.0, 24.0)
        } else {
            2.0
        };
        self
    }

    /// Sets the room above the strip its headings are turned on end in.
    ///
    /// Zero draws none, which is right where a legend names them instead or
    /// where there is only one column and the caption says what it is.
    pub fn heading_room(mut self, room: f64) -> Self {
        self.heading_room = if room.is_finite() { room.max(0.0) } else { 0.0 };
        self
    }

    /// The columns, in the order they are drawn.
    pub fn columns(&self) -> &[TraitColumn] {
        &self.columns
    }

    /// Whether there is anything to draw.
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// What the sheet holds about one row.
    pub fn values(&self, row: &str) -> Option<&Annotations> {
        self.rows.get(row)
    }

    /// How many of `rows` are named here.
    ///
    /// A join that names none of them draws every cell as an outline, which is
    /// a figure that looks finished and says nothing, so a caller checks this
    /// before drawing rather than after looking.
    pub fn covers<'a>(&self, rows: impl IntoIterator<Item = &'a str>) -> usize {
        rows.into_iter()
            .filter(|name| self.rows.contains_key(*name))
            .count()
    }

    /// The room the whole strip needs beside the rows, in pixels.
    pub fn strip_width(&self) -> f64 {
        if self.columns.is_empty() {
            return 0.0;
        }
        let cells: f64 = self.columns.iter().map(|column| column.width).sum();
        let gaps = (self.columns.len() as f64 - 1.0) * self.gap;
        cells + gaps + 8.0
    }

    /// The room the headings need above the rows, in pixels.
    pub fn heading_height(&self) -> f64 {
        if self.columns.is_empty() || self.heading_room <= 0.0 {
            return 0.0;
        }
        self.heading_room
    }

    /// A key naming every level and every ramp the columns drew.
    ///
    /// Nothing calls this on its own. A legend is a judgement about a figure
    /// rather than about a column, so the caller decides whether the figure
    /// needs one and where it goes, and this only spares them writing the
    /// colours down a second time and getting them wrong.
    ///
    /// It keys the strips this sheet draws beside the rows of a track, which
    /// number their levels in the order the sheet sorts its rows. A phylogeny
    /// does not draw from a sheet: it numbers the levels along its own walk,
    /// so its key is [`TreeTrack::legend`](crate::track::tree::TreeTrack::legend)
    /// and not this, which would name each colour beside another level.
    pub fn legend(&self, theme: &Theme) -> Legend {
        let mut legend = Legend::new();
        for column in &self.columns {
            let domain = self.domain(column);
            match column.scale {
                TraitScale::Continuous => {
                    // A column with no number in it has no range to show, and
                    // its ends would be the placeholders the count starts at.
                    if domain.minimum > domain.maximum {
                        continue;
                    }
                    let (low, high) = self.ramp_ends(column);
                    // The colours first and the values after, which is the
                    // order `Legend::ramp` takes them in. They were the other
                    // way round, so the ramp was painted with `fill="2015.17"`
                    // and labelled with two colour codes.
                    legend = legend.ramp(
                        column.label.clone(),
                        theme.muted.clone(),
                        theme.accent.clone(),
                        low,
                        high,
                    );
                }
                TraitScale::Categorical => {
                    for (level, index) in domain.keyed() {
                        legend = legend.key(
                            format!("{}: {level}", column.label),
                            theme.color(index).to_string(),
                        );
                    }
                }
            }
        }
        legend
    }

    /// The two numbers a continuous column's ramp runs between, as written.
    ///
    /// Taken from the column rather than written down beside it, so a legend
    /// cannot go on saying what the ramp used to be.
    pub fn ramp_ends(&self, column: &TraitColumn) -> (String, String) {
        let domain = self.domain(column);
        if domain.maximum <= domain.minimum {
            let one = crate::svg::text_rounded(domain.minimum, 3);
            return (one.clone(), one);
        }
        (
            crate::svg::text_rounded(domain.minimum, 3),
            crate::svg::text_rounded(domain.maximum, 3),
        )
    }

    /// The levels and the range one column covers over every row named here,
    /// dealt the palette in the column's own order where it has one.
    fn domain(&self, column: &TraitColumn) -> TraitDomain {
        TraitDomain::ordered(&column.levels, self.stated(&column.key))
    }

    /// Draws the strip beside rows that have already been laid out.
    ///
    /// `rows` is the drawn order with each row's top and height, which is the
    /// host track's arithmetic and not this module's: a row here lines up with
    /// a row there because it was given the same number, not because both
    /// worked it out.
    pub(crate) fn draw(&self, ctx: &mut DrawContext<'_>, area: Rect, rows: &[(String, f64, f64)]) {
        if self.columns.is_empty() {
            return;
        }
        let size = (ctx.theme.font_size - 2.0).max(6.0);
        let room = self.heading_height();
        let mut x = area.x + 4.0;

        for column in &self.columns {
            let domain = self.domain(column);
            if room > 0.0 {
                // Turned on end because a column is narrower than its name and
                // will stay that way: a strip wide enough to caption flat is a
                // strip wide enough to be mistaken for the data.
                let heading = fit_text(&column.label, ctx.px(room) - 6.0, size);
                ctx.svg.text_rotated(
                    (
                        x + column.width / 2.0 + size * 0.35,
                        area.y + ctx.px(room) - 6.0,
                    ),
                    -90.0,
                    &heading,
                    &ctx.theme.muted,
                    size,
                    Anchor::Start,
                );
            }

            let cells: Vec<TraitRow<'_>> = rows
                .iter()
                .map(|(name, top, height)| TraitRow {
                    name,
                    top: *top,
                    height: *height,
                    value: self.rows.get(name).and_then(|held| held.get(&column.key)),
                })
                .collect();
            draw_column(ctx, column, &domain, x, None, &cells);
            x += column.width + self.gap;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::read::sheet::sheet;
    use crate::region::Region;
    use crate::track::Track;
    use crate::{MatrixRow, MatrixTrack};

    const SHEET: &str = "\
sample\tlineage\thost\tdepth
A\tL4\thuman\t72.5
B\tL2\tbovine\t61
C\tL4\t\t48.2
D\tL1\thuman\t95
";

    fn traits() -> Traits {
        let held = sheet(SHEET).expect("a sheet");
        let columns = held.columns.clone();
        Traits::new(held.rows).spread(columns)
    }

    fn matrix() -> MatrixTrack {
        MatrixTrack::new(
            vec![10, 20],
            vec![
                MatrixRow::new("A", vec![1.0, 0.0]),
                MatrixRow::new("B", vec![0.0, 1.0]),
                MatrixRow::new("C", vec![1.0, 1.0]),
                MatrixRow::new("D", vec![0.0, 0.0]),
            ],
        )
    }

    fn drawn(track: MatrixTrack) -> String {
        Figure::new(Region::new("chr1", 0, 40).unwrap())
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    #[test]
    fn a_column_is_a_ramp_when_every_stated_value_is_a_number() {
        let traits = traits();
        let scales: Vec<TraitScale> = traits.columns().iter().map(|c| c.scale()).collect();
        assert_eq!(
            scales,
            [
                TraitScale::Categorical,
                TraitScale::Categorical,
                TraitScale::Continuous
            ]
        );
    }

    #[test]
    fn a_level_keeps_its_colour_when_a_row_is_added() {
        // The whole reason the domain numbers levels as it meets them. A file
        // with one more sample in it must not repaint the samples that were
        // already there, or two runs of the same figure disagree.
        let first = traits();
        let more = sheet(&format!("{SHEET}E\tL2\tbovine\t50\n")).expect("a sheet");
        let columns = more.columns.clone();
        let second = Traits::new(more.rows).spread(columns);

        let theme = Theme::light();
        for name in ["A", "B", "C", "D"] {
            let one = first.domain(&first.columns()[0]);
            let two = second.domain(&second.columns()[0]);
            let value = first.values(name).and_then(|held| held.get("lineage"));
            assert_eq!(
                one.color(&first.columns()[0], value, &theme),
                two.color(&second.columns()[0], value, &theme),
                "{name} changed colour"
            );
        }
    }

    #[test]
    fn a_ramp_in_the_key_runs_between_the_theme_colours_and_is_labelled_with_the_values() {
        // `Legend::ramp` takes the two colours and then the two values, and
        // the key passed them the other way round: the ramp was painted with
        // `fill="48.2"` and labelled with two colour codes.
        let theme = Theme::light();
        let legend = traits().legend(&theme);
        let ramps: Vec<_> = legend
            .items()
            .iter()
            .filter_map(|item| match item {
                crate::track::legend::LegendItem::Ramp {
                    label,
                    from,
                    to,
                    low,
                    high,
                } => Some((label, from, to, low, high)),
                _ => None,
            })
            .collect();
        assert_eq!(ramps.len(), 1, "one continuous column, one ramp");
        let (label, from, to, low, high) = ramps[0];
        assert_eq!(label, "depth");
        assert_eq!((from, to), (&theme.muted, &theme.accent));
        assert_eq!((low.as_str(), high.as_str()), ("48.2", "95"));

        let svg = Figure::new(Region::new("chr1", 0, 40).unwrap())
            .show_region_label(false)
            .push(crate::track::legend::LegendTrack::new(legend))
            .to_svg();
        assert!(
            !svg.contains("fill=\"48.2\""),
            "a value painted as a colour"
        );
        assert!(svg.contains(">48.2</text>"), "the low end is not written");
    }

    #[test]
    fn a_column_with_no_number_in_it_is_not_keyed_as_a_ramp() {
        // A column asked for as continuous and holding no number has no
        // range, and its ends were the placeholders an empty count starts at,
        // printed as a range three hundred digits long.
        let held = sheet("sample\tdepth\nA\t\nB\t\n").expect("a sheet");
        let traits = Traits::new(held.rows).column(TraitColumn::continuous("depth"));
        let legend = traits.legend(&Theme::light());
        assert!(legend.items().is_empty(), "{:?}", legend.items());
    }

    #[test]
    fn a_sample_appended_to_the_sheet_never_repaints_the_ones_above_it() {
        // Dealt in the order the names sort, a new sample called AA with a
        // new lineage was met first and took the first colour from L4, and
        // every sample of every lineage changed. Dealt in the order the file
        // lists them, an appended row only ever adds a colour.
        let dealt = |text: &str| {
            let held = sheet(text).expect("a sheet");
            let traits = Traits::from_sheet(&held).spread(["lineage"]);
            let theme = Theme::light();
            let column = &traits.columns()[0];
            let domain = traits.domain(column);
            ["A", "B", "D"].map(|name| {
                let value = traits.values(name).and_then(|held| held.get("lineage"));
                domain.color(column, value, &theme)
            })
        };
        assert_eq!(dealt(&format!("{SHEET}AA\tL7\thuman\t50\n")), dealt(SHEET));
        assert_eq!(dealt(SHEET)[0].as_deref(), Some(Theme::light().color(0)));
    }

    #[test]
    fn a_column_carries_the_order_its_sheet_lists_the_levels_in() {
        // What a phylogeny is handed, so it deals the palette the same way.
        let held = sheet(SHEET).expect("a sheet");
        let traits = Traits::from_sheet(&held).spread(["lineage", "depth"]);
        assert_eq!(traits.columns()[0].level_order(), ["L4", "L2", "L1"]);
        assert!(
            traits.columns()[1].level_order().is_empty(),
            "a ramp has no levels"
        );
        // An order given by hand is kept rather than replaced.
        let chosen = Traits::from_sheet(&held)
            .column(TraitColumn::categorical("lineage").levels(["L1", "L2", "L4"]));
        assert_eq!(chosen.columns()[0].level_order(), ["L1", "L2", "L4"]);
        let legend = chosen.legend(&Theme::light());
        let first = &legend.items()[0];
        let crate::track::legend::LegendItem::Key { label, color, .. } = first else {
            panic!("a categorical key: {first:?}");
        };
        assert_eq!(
            (label.as_str(), color.as_str()),
            ("lineage: L1", Theme::light().color(0))
        );
    }

    #[test]
    fn a_missing_value_is_an_outline_and_never_a_colour() {
        // C has no host. Drawn in any fill it would read as a level, and the
        // level it would read as is whichever one shares the colour.
        let svg = drawn(matrix().traits(traits()));
        assert!(svg.contains("C; host missing"), "no tooltip saying so");
        assert!(svg.contains("A; host human"));
        let missing = svg
            .split("C; host missing")
            .nth(1)
            .expect("a cell after the title");
        assert!(
            missing[..120].contains("stroke="),
            "the missing cell is filled rather than outlined: {}",
            &missing[..120]
        );
    }

    #[test]
    fn more_levels_than_the_palette_carries_get_a_shape_as_well_as_a_hue() {
        // Six colours go round, and two countries sharing a swatch is a figure
        // stating something untrue. A symbol separates twenty-four.
        let many = (0..9).map(|i| format!("S{i}\tC{i}\n")).collect::<String>();
        let held = sheet(&format!("sample\tcountry\n{many}")).expect("a sheet");
        let columns = held.columns.clone();
        let traits = Traits::new(held.rows).spread(columns);
        assert_eq!(traits.columns()[0].trait_style(), TraitStyle::Symbol);

        let held = sheet("sample\tcountry\nS0\tES\nS1\tFR\n").expect("a sheet");
        let columns = held.columns.clone();
        let few = Traits::new(held.rows).spread(columns);
        assert_eq!(few.columns()[0].trait_style(), TraitStyle::Strip);
    }

    #[test]
    fn the_strip_takes_room_from_the_track_rather_than_from_the_figure() {
        let bare = matrix();
        let with = matrix().traits(traits());
        let theme = Theme::light();
        assert!(
            with.y_axis_width(&theme) > bare.y_axis_width(&theme),
            "the strip was not reserved"
        );
        assert_eq!(
            with.y_axis_width(&theme) - bare.y_axis_width(&theme),
            traits().strip_width()
        );
    }

    #[test]
    fn headings_are_reserved_in_the_height_and_pushed_the_rows_down() {
        let scale = crate::scale::Scale::new(&Region::new("chr1", 0, 40).unwrap(), 0.0, 100.0);
        let bare = matrix().height(&scale);
        let with = matrix().traits(traits()).height(&scale);
        assert_eq!(with - bare, traits().heading_height());

        let none = matrix().traits(traits().heading_room(0.0));
        assert_eq!(none.height(&scale), bare);
    }

    #[test]
    fn a_sheet_that_names_nobody_is_a_strip_of_nothing_and_says_so() {
        // Not an error here: refusing belongs to whoever joined the two files.
        // What matters is that every cell reads as absent rather than as a
        // level, so the caller who did not check can still see it.
        let held = sheet("sample\tx\nZZ\ta\n").expect("a sheet");
        let columns = held.columns.clone();
        let traits = Traits::new(held.rows).spread(columns);
        assert_eq!(traits.covers(["A", "B", "C", "D"]), 0);
        let svg = drawn(matrix().traits(traits));
        assert_eq!(svg.matches("; x missing").count(), 4);
    }

    #[test]
    fn a_legend_names_every_level_and_the_two_ends_of_every_ramp() {
        let traits = traits();
        let legend = traits.legend(&Theme::light());
        let labels: Vec<String> = legend
            .items()
            .iter()
            .map(|item| format!("{item:?}"))
            .collect();
        let text = labels.join(" ");
        assert!(text.contains("lineage: L4"), "{text}");
        assert!(text.contains("lineage: L2"), "{text}");
        assert!(text.contains("lineage: L1"), "{text}");
        assert!(text.contains("host: human"), "{text}");
        // The ramp ends come off the column rather than being written down.
        assert_eq!(
            traits.ramp_ends(&traits.columns()[2]),
            ("48.2".into(), "95".into())
        );
    }

    #[test]
    fn a_column_of_one_number_has_no_gradient_and_says_so_at_both_ends() {
        let held = sheet("sample\tdepth\nA\t30\nB\t30\n").expect("a sheet");
        let columns = held.columns.clone();
        let traits = Traits::new(held.rows).spread(columns);
        assert_eq!(
            traits.ramp_ends(&traits.columns()[0]),
            ("30".into(), "30".into())
        );
    }

    #[test]
    fn the_columns_asked_for_are_the_columns_drawn_in_the_order_asked() {
        let held = sheet(SHEET).expect("a sheet");
        let traits = Traits::new(held.rows).spread(["host", "lineage"]);
        let names: Vec<&str> = traits.columns().iter().map(|c| c.key()).collect();
        assert_eq!(names, ["host", "lineage"]);
    }

    #[test]
    fn a_key_no_row_mentions_is_still_a_column() {
        // An attribute nobody in this figure has is a fact about the figure,
        // and a column of outlines states it. Dropping it would leave the
        // command that asked for it looking as though it had worked.
        let held = sheet(SHEET).expect("a sheet");
        let traits = Traits::new(held.rows).spread(["ward"]);
        assert_eq!(traits.columns().len(), 1);
        let svg = drawn(matrix().traits(traits));
        assert_eq!(svg.matches("; ward missing").count(), 4);
    }

    #[test]
    fn a_strip_of_no_columns_takes_no_room_at_all() {
        let none = Traits::default();
        assert_eq!(none.strip_width(), 0.0);
        assert_eq!(none.heading_height(), 0.0);
        assert!(none.is_empty());
    }
}
