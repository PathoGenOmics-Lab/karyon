//! A multiple sequence alignment, row by row.
//!
//! An [`MsaSequence`] is a name and the residues under it, and [`MsaTrack`] is
//! the stack of them. Two decisions about that stack are worth settling before
//! any of the styling: what the figure's coordinate actually counts, and how
//! much of the alignment is worth painting at all.
//!
//! # A column is not a position
//!
//! An alignment has its own coordinate system: column indices, gaps included.
//! Those are not genomic positions, and the two must not be confused, so the
//! figure's region is the column space. An alignment 900 columns wide is
//! `Region::new("alignment", 0, 900)`, and the ruler under it counts columns.
//! Ungapping a row back to reference coordinates is a real operation with real
//! decisions in it, and this crate does not do it silently.
//!
//! # Most of an alignment is not worth ink
//!
//! A wall of coloured residues is pretty and says very little, because in a
//! real alignment most cells agree and the agreement is the noise. The default
//! is therefore [`MsaDisplay::Differences`]: rows are drawn as a quiet bar and
//! only what disagrees with the comparison row is painted. Which row that is
//! matters as much as the display, and it is the consensus until
//! [`MsaTrack::compare_to`] names one, since comparing against a sequence no
//! more privileged than the rest turns the figure into a story about that
//! sequence. Conservation itself belongs above the alignment rather than inside
//! it, and [`LogoTrack::from_sequences`](crate::LogoTrack::from_sequences)
//! takes the same sequences this track does.
//!
//! # A deep alignment is not shown whole
//!
//! Rows are capped, at forty by default, and what the cap hid is counted in the
//! corner rather than left to be noticed, because a figure that quietly drops
//! half its samples is worse than one that is visibly truncated.
//! [`MsaTrack::max_rows`] moves the cap or lifts it. Within a row, runs of
//! equally coloured cells are merged into one rectangle on the way out, which
//! is what keeps a wide alignment from becoming a file no viewer will open.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{contrast_ink, mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::traits::Traits;
use crate::track::tree::{draw_tree, leaf_order, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// One aligned sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct MsaSequence {
    /// Name, drawn beside the row.
    pub name: String,
    /// Residues, one per alignment column, gaps included.
    pub residues: Vec<u8>,
}

impl MsaSequence {
    /// A sequence named `name`.
    pub fn new(name: impl Into<String>, residues: impl Into<Vec<u8>>) -> Self {
        MsaSequence {
            name: name.into(),
            residues: residues.into(),
        }
    }

    /// The residue in a column, if the row reaches that far.
    pub fn residue(&self, column: usize) -> Option<u8> {
        self.residues.get(column).copied()
    }

    /// How many columns the row occupies.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Whether the row has no residues at all.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
}

/// Whether a gap character is a gap.
pub fn is_gap(residue: u8) -> bool {
    matches!(residue, b'-' | b'.' | b'~' | b' ')
}

/// A count with its noun, in the number the count calls for.
///
/// The plural is given rather than derived, since the nouns a tooltip needs
/// are not all one letter apart from their singular. The number is grouped,
/// because every number in a tooltip is: an alignment thirty kilobases wide
/// reads `30,000 columns`, and a rule that grouped only the numbers a reader
/// found long would make the punctuation a fact about the data.
fn count(n: usize, one: &str, many: &str) -> String {
    let n = group_thousands(n as u64);
    if n == "1" {
        format!("{n} {one}")
    } else {
        format!("{n} {many}")
    }
}

/// What gets painted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MsaDisplay {
    /// Only cells that disagree with the comparison row, which is the default
    /// and is what makes a large alignment readable at all.
    Differences,
    /// Every residue, agreements included.
    Bases,
}

/// What a residue's colour means.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MsaColoring {
    /// Nucleotides, on the theme's base colours.
    Nucleotide,
    /// Amino acids, grouped into six physicochemical classes.
    Residue,
    /// One colour for every difference, when the identity of the residue is
    /// not the point.
    Uniform,
}

/// The six classes amino acids are grouped into here.
///
/// Six because that is how many hues the validated categorical palette has, and
/// the groupings are the conventional ones: cysteine sits with the hydrophobics,
/// histidine with the positives, tyrosine with the polars. Glycine and proline
/// keep their own classes, since a figure of a protein alignment is usually
/// looking for exactly those two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueClass {
    /// A, V, L, I, M, F, W, C.
    Hydrophobic,
    /// K, R, H.
    Positive,
    /// D, E.
    Negative,
    /// S, T, N, Q, Y.
    Polar,
    /// G.
    Glycine,
    /// P.
    Proline,
}

impl ResidueClass {
    /// The class of one residue, or `None` for a gap or an unknown symbol.
    pub fn of(residue: u8) -> Option<Self> {
        match residue.to_ascii_uppercase() {
            b'A' | b'V' | b'L' | b'I' | b'M' | b'F' | b'W' | b'C' => {
                Some(ResidueClass::Hydrophobic)
            }
            b'K' | b'R' | b'H' => Some(ResidueClass::Positive),
            b'D' | b'E' => Some(ResidueClass::Negative),
            b'S' | b'T' | b'N' | b'Q' | b'Y' => Some(ResidueClass::Polar),
            b'G' => Some(ResidueClass::Glycine),
            b'P' => Some(ResidueClass::Proline),
            _ => None,
        }
    }

    /// Which slot of the categorical palette this class takes.
    pub fn palette_index(self) -> usize {
        match self {
            ResidueClass::Hydrophobic => 0,
            ResidueClass::Positive => 1,
            ResidueClass::Negative => 2,
            ResidueClass::Polar => 3,
            ResidueClass::Glycine => 4,
            ResidueClass::Proline => 5,
        }
    }
}

/// A multiple sequence alignment.
///
/// ```
/// use karyon::{Figure, MsaSequence, MsaTrack, Region};
///
/// let rows = vec![
///     MsaSequence::new("H37Rv", b"ACGTACGTAC".to_vec()),
///     MsaSequence::new("CDC1551", b"ACGTTCGTAC".to_vec()),
///     MsaSequence::new("Beijing", b"ACGT-CGTAC".to_vec()),
/// ];
///
/// let svg = Figure::new(Region::new("alignment", 0, 10).unwrap())
///     .push(MsaTrack::new(rows).label("isolates"))
///     .to_svg();
/// assert!(svg.contains("H37Rv"));
/// ```
#[derive(Debug, Clone)]
pub struct MsaTrack {
    sequences: Vec<MsaSequence>,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    display: MsaDisplay,
    coloring: MsaColoring,
    compare_to: Option<usize>,
    max_rows: Option<usize>,
    show_names: bool,
    show_letters: bool,
    letter_threshold: f64,
    match_color: Option<String>,
    gap_color: Option<String>,
    uniform_color: Option<String>,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
    traits: Traits,
}

impl MsaTrack {
    /// A track holding `sequences`.
    pub fn new(sequences: impl Into<Vec<MsaSequence>>) -> Self {
        MsaTrack {
            sequences: sequences.into(),
            label: None,
            row_height: 12.0,
            row_gap: 1.0,
            display: MsaDisplay::Differences,
            coloring: MsaColoring::Nucleotide,
            compare_to: None,
            max_rows: Some(40),
            show_names: true,
            show_letters: true,
            letter_threshold: 7.0,
            match_color: None,
            gap_color: None,
            uniform_color: None,
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

    /// Sets the height of one row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(1.0);
        self
    }

    /// Sets the gap between rows.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Chooses between painting differences and painting everything.
    pub fn display(mut self, display: MsaDisplay) -> Self {
        self.display = display;
        self
    }

    /// Chooses what a colour means.
    pub fn coloring(mut self, coloring: MsaColoring) -> Self {
        self.coloring = coloring;
        self
    }

    /// Compares every row against row `index` rather than against the
    /// consensus.
    ///
    /// Pick a row when one of the sequences is the reference and the figure is
    /// about departures from it. Leave it alone when no row is privileged and
    /// the question is where the alignment disagrees with itself.
    pub fn compare_to(mut self, index: usize) -> Self {
        self.compare_to = Some(index);
        self
    }

    /// Caps how many rows are drawn, or removes the cap with `None`.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|r| r.max(1));
        self
    }

    /// Draws or hides the sequence names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Draws or hides residue letters once a column is wide enough.
    pub fn show_letters(mut self, show: bool) -> Self {
        self.show_letters = show;
        self
    }

    /// Sets the colour of a cell that agrees with the comparison row.
    pub fn match_color(mut self, color: impl Into<String>) -> Self {
        self.match_color = Some(color.into());
        self
    }

    /// Sets the colour of a gap.
    pub fn gap_color(mut self, color: impl Into<String>) -> Self {
        self.gap_color = Some(color.into());
        self
    }

    /// Sets the colour used under [`MsaColoring::Uniform`].
    pub fn uniform_color(mut self, color: impl Into<String>) -> Self {
        self.uniform_color = Some(color.into());
        self
    }

    /// Draws a phylogeny beside the alignment and sorts sequences by descent.
    ///
    /// Rows match leaves by exact name. Sequences absent from the tree remain
    /// at the bottom, and a selected comparison sequence follows its row when
    /// the alignment is reordered.
    pub fn tree(mut self, tree: Tree) -> Self {
        let names: Vec<String> = self
            .sequences
            .iter()
            .map(|sequence| sequence.name.clone())
            .collect();
        let order = leaf_order(&tree, &names);
        if let Some(comparison) = self.compare_to {
            self.compare_to = order.iter().position(|index| *index == comparison);
        }
        self.sequences = order
            .iter()
            .map(|index| self.sequences[*index].clone())
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

    /// Chooses a phylogram or cladogram for the tree beside the alignment.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// Attaches metadata columns drawn between the sequence names and the rows.
    ///
    /// Joined by name, so the columns follow whatever order the sequences are
    /// already in, including the one a phylogeny put them in.
    pub fn traits(mut self, traits: Traits) -> Self {
        self.traits = traits;
        self
    }

    /// The metadata columns attached to this track.
    pub fn attached_traits(&self) -> &Traits {
        &self.traits
    }

    /// The tree attached to the alignment, when present.
    pub fn attached_tree(&self) -> Option<&Tree> {
        self.tree.as_ref()
    }

    /// The sequences in the track.
    pub fn sequences(&self) -> &[MsaSequence] {
        &self.sequences
    }

    /// How many columns the alignment has, which is its longest row.
    pub fn columns(&self) -> usize {
        self.sequences
            .iter()
            .map(|row| row.residues.len())
            .max()
            .unwrap_or(0)
    }

    /// The most common non-gap residue in each column.
    ///
    /// Ties go to whichever residue was seen first, scanning rows in order, so
    /// the consensus of a given alignment never changes between runs. A column
    /// of nothing but gaps has no consensus and comes back as a gap.
    pub fn consensus(&self) -> Vec<u8> {
        (0..self.columns())
            .map(|column| {
                let mut tally: Vec<(u8, usize)> = Vec::new();
                for row in &self.sequences {
                    let Some(residue) = row.residue(column) else {
                        continue;
                    };
                    if is_gap(residue) {
                        continue;
                    }
                    let residue = residue.to_ascii_uppercase();
                    match tally.iter_mut().find(|(r, _)| *r == residue) {
                        Some((_, count)) => *count += 1,
                        None => tally.push((residue, 1)),
                    }
                }
                // Not max_by_key: that returns the last maximum, and the
                // promise above is the first, which is what keeps the same
                // alignment giving the same consensus every time.
                tally
                    .into_iter()
                    .reduce(|best, next| if next.1 > best.1 { next } else { best })
                    .map_or(b'-', |(residue, _)| residue)
            })
            .collect()
    }

    /// The row every other row is compared against.
    ///
    /// Either the one [`MsaTrack::compare_to`] named, or the consensus.
    pub fn comparison(&self) -> Vec<u8> {
        match self.compare_to.and_then(|index| self.sequences.get(index)) {
            Some(row) => row.residues.clone(),
            None => self.consensus(),
        }
    }

    /// How many rows are drawn and how many the cap hid.
    pub fn visible_rows(&self) -> (usize, usize) {
        match self.max_rows {
            Some(cap) if self.sequences.len() > cap => (cap, self.sequences.len() - cap),
            _ => (self.sequences.len(), 0),
        }
    }

    /// Colour of one cell, or `None` when nothing should be drawn there.
    fn cell_color(
        &self,
        residue: Option<u8>,
        against: Option<u8>,
        theme: &Theme,
    ) -> Option<String> {
        let Some(residue) = residue else {
            // Past the end of a short row: not a gap, just absent.
            return None;
        };
        if is_gap(residue) {
            return Some(
                self.gap_color
                    .clone()
                    .unwrap_or_else(|| mix(theme.surface(), &theme.rule, 0.5)),
            );
        }

        let agrees =
            against.is_some_and(|other| !is_gap(other) && other.eq_ignore_ascii_case(&residue));
        if self.display == MsaDisplay::Differences && agrees {
            return Some(
                self.match_color
                    .clone()
                    .unwrap_or_else(|| mix(theme.surface(), &theme.rule, 0.28)),
            );
        }

        Some(match self.coloring {
            MsaColoring::Nucleotide => theme.bases.of(residue).to_string(),
            MsaColoring::Residue => match ResidueClass::of(residue) {
                Some(class) => theme.color(class.palette_index()).to_string(),
                None => theme.bases.other.clone(),
            },
            MsaColoring::Uniform => self
                .uniform_color
                .clone()
                .unwrap_or_else(|| theme.accent.clone()),
        })
    }
}

impl Track for MsaTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        let (rows, _) = self.visible_rows();
        let rows = rows.max(1) as f64;
        rows * self.row_height + (rows - 1.0).max(0.0) * self.row_gap + self.traits.heading_height()
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let tree = self.tree.as_ref().map_or(0.0, |_| self.tree_width);
        let strip = self.traits.strip_width();
        if !self.show_names || self.sequences.is_empty() {
            return tree + strip;
        }
        let size = (theme.font_size - 2.0).min(self.row_height);
        let (rows, _) = self.visible_rows();
        self.sequences
            .iter()
            .take(rows)
            .map(|row| text_width(&row.name, size))
            .fold(0.0f64, f64::max)
            + 8.0
            + tree
            + strip
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let comparison = self.comparison();
        let (rows, hidden) = self.visible_rows();
        let head = ctx.px(self.traits.heading_height());
        let strip = self.traits.strip_width();
        let name_size = (ctx.theme.font_size - 2.0).min(self.row_height);
        let px_per_column = ctx.scale.px_per_bp();
        let draw_letters = self.show_letters && px_per_column >= self.letter_threshold;
        let tree_color = mix(&ctx.theme.foreground, &ctx.theme.muted, 0.45);

        let first = ctx.region.start() as usize;
        let last = (ctx.region.end() as usize).min(self.columns());

        if let Some(tree) = &self.tree {
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: band.y + head,
                w: (self.tree_width - 6.0).max(1.0),
                h: (band.h - head).max(1.0),
            };
            draw_tree(
                ctx.svg,
                tree,
                area,
                self.row_height + self.row_gap,
                band.y + head + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &tree_color,
                    width: ctx.theme.tokens.stroke.max(1.1),
                    mirror: false,
                },
            );
            ctx.svg.line(
                band.x - 2.0,
                band.y + head + 1.0,
                band.x - 2.0,
                band.bottom() - 1.0,
                &ctx.theme.rule,
                ctx.theme.tokens.hairline.max(0.8),
            );
        }

        for (index, row) in self.sequences.iter().take(rows).enumerate() {
            let top = band.y + head + index as f64 * (self.row_height + self.row_gap);

            // The row is named, and never the cell. A tooltip per residue is
            // rows times columns elements, which is tens of thousands of them
            // on an alignment worth drawing, and the name of a cell is its
            // row's name anyway.
            ctx.svg
                .begin_titled(&self.row_title(row, &comparison, first, last));

            // Neighbouring cells of the same colour are merged into one
            // rectangle. Most of an alignment agrees with itself, so this is
            // the difference between a figure and a file no viewer will open.
            let mut run: Option<(usize, String)> = None;
            for column in first..last {
                let color = self.cell_color(
                    row.residue(column),
                    comparison.get(column).copied(),
                    ctx.theme,
                );
                match (&mut run, color) {
                    (Some((_, current)), Some(next)) if *current == next => {}
                    (slot, next) => {
                        if let Some((start, current)) = slot.take() {
                            self.paint_run(ctx, start, column, top, &current);
                        }
                        *slot = next.map(|color| (column, color));
                    }
                }
            }
            if let Some((start, color)) = run {
                self.paint_run(ctx, start, last, top, &color);
            }

            if draw_letters {
                self.paint_letters(ctx, row, &comparison, first, last, top);
            }

            if self.show_names && ctx.axis.w > strip {
                ctx.svg.text(
                    ctx.axis.right() - strip - 4.0,
                    top + self.row_height / 2.0 + name_size * 0.35,
                    &row.name,
                    &ctx.theme.foreground,
                    name_size,
                    Anchor::End,
                );
            }

            ctx.svg.end_group();
        }

        if self.display == MsaDisplay::Bases && px_per_column >= 5.0 {
            for column in first.saturating_add(1)..last {
                let x = ctx.scale.x(column as u64);
                ctx.svg.line(
                    x,
                    band.y + head,
                    x,
                    band.bottom(),
                    ctx.theme.surface(),
                    0.65,
                );
            }
        }

        let placed: Vec<(String, f64, f64)> = self
            .sequences
            .iter()
            .take(rows)
            .enumerate()
            .map(|(index, row)| {
                (
                    row.name.clone(),
                    band.y + head + index as f64 * (self.row_height + self.row_gap),
                    self.row_height,
                )
            })
            .collect();
        self.traits.draw(
            ctx,
            Rect {
                x: ctx.axis.right() - strip,
                y: band.y,
                w: strip,
                h: band.h,
            },
            &placed,
        );

        if hidden > 0 {
            let text = format!("+{hidden} sequences not shown");
            let size = ctx.theme.font_size - 1.0;
            let width = text_width(&text, size) + 6.0;
            ctx.svg.rect_opacity(
                band.right() - width,
                band.bottom() - size - 3.0,
                width,
                size + 3.0,
                &ctx.theme.background,
                0.8,
            );
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 2.0,
                &text,
                &ctx.theme.muted,
                size,
                Anchor::End,
            );
        }
    }
}

impl MsaTrack {
    /// What a row is and how much of it departs from the comparison row.
    ///
    /// This is the claim the track makes about a row, so it is what a pointer
    /// resting on the row is told. Counting is over the columns in view rather
    /// than the whole alignment, since a tooltip that described columns nobody
    /// can see would disagree with the figure it sits on.
    ///
    /// Gaps are counted apart from mismatches because they are a third thing:
    /// a gap disagrees with the comparison row without carrying a residue that
    /// could have agreed, and it is painted in its own colour to say so.
    fn row_title(&self, row: &MsaSequence, comparison: &[u8], first: usize, last: usize) -> String {
        let mut columns = 0usize;
        let mut mismatches = 0usize;
        let mut gaps = 0usize;
        for column in first..last {
            let Some(residue) = row.residue(column) else {
                continue;
            };
            columns += 1;
            if is_gap(residue) {
                gaps += 1;
                continue;
            }
            let agrees = comparison
                .get(column)
                .copied()
                .is_some_and(|other| !is_gap(other) && other.eq_ignore_ascii_case(&residue));
            if !agrees {
                mismatches += 1;
            }
        }

        let mut title = String::new();
        if !row.name.is_empty() {
            title.push_str(&row.name);
            title.push_str(", ");
        }
        title.push_str(&count(columns, "column", "columns"));
        title.push_str(", ");
        title.push_str(&count(mismatches, "mismatch", "mismatches"));
        if gaps > 0 {
            title.push_str(", ");
            title.push_str(&count(gaps, "gap", "gaps"));
        }
        title
    }

    /// Paints one merged run of equally coloured columns.
    fn paint_run(&self, ctx: &mut DrawContext<'_>, from: usize, to: usize, top: f64, color: &str) {
        let x = ctx.scale.x(from as u64);
        let width = (ctx.scale.x(to as u64) - x).max(0.6);
        ctx.svg.rect(x, top, width, self.row_height, color);
    }

    /// Paints the residue letters of one row.
    fn paint_letters(
        &self,
        ctx: &mut DrawContext<'_>,
        row: &MsaSequence,
        comparison: &[u8],
        first: usize,
        last: usize,
        top: f64,
    ) {
        let size = (ctx.scale.px_per_bp() * 0.72).min(self.row_height * 0.85);
        for column in first..last {
            let Some(residue) = row.residue(column) else {
                continue;
            };
            if is_gap(residue) {
                continue;
            }
            let x = ctx.scale.x(column as u64);
            let width = ctx.scale.x(column as u64 + 1) - x;
            // The ink is chosen against the cell the letter sits on. Fixed to
            // the foreground it fell to two to one on the darker class
            // colours, which is under every contrast floor there is.
            let ink = self
                .cell_color(Some(residue), comparison.get(column).copied(), ctx.theme)
                .map_or_else(
                    || ctx.theme.foreground.clone(),
                    |cell| contrast_ink(&cell).to_string(),
                );
            // Drawn rather than stretched to the column. A logo letter has to
            // fill an exact box because its height is the datum; here it is
            // just a letter, and forcing an I to the width of a W turns it
            // into a block.
            ctx.svg.text(
                x + width / 2.0,
                top + self.row_height - (self.row_height - size) / 2.0,
                &(residue as char).to_ascii_uppercase().to_string(),
                &ink,
                size / ctx.theme.cap_height_ratio.max(0.1),
                Anchor::Middle,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn alignment() -> Vec<MsaSequence> {
        vec![
            MsaSequence::new("H37Rv", b"ACGTACGTAC".to_vec()),
            MsaSequence::new("CDC1551", b"ACGTTCGTAC".to_vec()),
            MsaSequence::new("Beijing", b"ACGT-CGTAC".to_vec()),
            MsaSequence::new("Erdman", b"ACGTTCGTAC".to_vec()),
        ]
    }

    fn region() -> Region {
        Region::new("alignment", 0, 10).unwrap()
    }

    #[test]
    fn gap_characters_are_all_the_usual_ones() {
        for gap in *b"-.~ " {
            assert!(is_gap(gap), "{} should be a gap", gap as char);
        }
        for residue in *b"ACXN" {
            assert!(!is_gap(residue));
        }
    }

    #[test]
    fn the_column_count_is_the_longest_row() {
        let ragged = vec![
            MsaSequence::new("a", b"ACGT".to_vec()),
            MsaSequence::new("b", b"AC".to_vec()),
        ];
        assert_eq!(MsaTrack::new(ragged).columns(), 4);
        assert_eq!(MsaTrack::new(Vec::new()).columns(), 0);
    }

    #[test]
    fn the_consensus_takes_the_commonest_residue_per_column() {
        // Column 4 is A, T, gap, T, so the consensus is T.
        let consensus = MsaTrack::new(alignment()).consensus();
        assert_eq!(consensus.len(), 10);
        assert_eq!(consensus[0], b'A');
        assert_eq!(consensus[4], b'T');
    }

    #[test]
    fn the_consensus_ignores_gaps_but_survives_a_column_of_them() {
        let rows = vec![
            MsaSequence::new("a", b"-A".to_vec()),
            MsaSequence::new("b", b"-A".to_vec()),
        ];
        let consensus = MsaTrack::new(rows).consensus();
        assert_eq!(consensus[0], b'-', "nothing but gaps has no consensus");
        assert_eq!(consensus[1], b'A');
    }

    #[test]
    fn a_tie_is_broken_by_which_row_came_first() {
        let rows = vec![
            MsaSequence::new("a", b"G".to_vec()),
            MsaSequence::new("b", b"C".to_vec()),
        ];
        let track = MsaTrack::new(rows);
        assert_eq!(track.consensus()[0], b'G');
        // And it does not wander between calls.
        assert_eq!(track.consensus(), track.consensus());
    }

    #[test]
    fn the_comparison_row_is_the_consensus_unless_one_is_named() {
        let track = MsaTrack::new(alignment());
        assert_eq!(track.comparison(), track.consensus());

        let pinned = MsaTrack::new(alignment()).compare_to(0);
        assert_eq!(pinned.comparison(), b"ACGTACGTAC".to_vec());

        // An index past the end falls back rather than panicking.
        let silly = MsaTrack::new(alignment()).compare_to(99);
        assert_eq!(silly.comparison(), silly.consensus());
    }

    #[test]
    fn agreement_is_case_insensitive() {
        let theme = Theme::light();
        let track = MsaTrack::new(alignment());
        let upper = track.cell_color(Some(b'A'), Some(b'A'), &theme);
        let lower = track.cell_color(Some(b'a'), Some(b'A'), &theme);
        assert_eq!(upper, lower, "a soft masked base still agrees");
    }

    #[test]
    fn a_gap_never_agrees_with_anything() {
        let theme = Theme::light();
        let track = MsaTrack::new(alignment());
        let gap = track.cell_color(Some(b'-'), Some(b'-'), &theme);
        let matched = track.cell_color(Some(b'A'), Some(b'A'), &theme);
        assert_ne!(gap, matched, "two gaps are not an agreement");
    }

    #[test]
    fn past_the_end_of_a_row_nothing_is_drawn() {
        let theme = Theme::light();
        let track = MsaTrack::new(alignment());
        assert_eq!(track.cell_color(None, Some(b'A'), &theme), None);
    }

    #[test]
    fn showing_bases_paints_what_showing_differences_leaves_quiet() {
        let theme = Theme::light();
        let differences = MsaTrack::new(alignment());
        let everything = MsaTrack::new(alignment()).display(MsaDisplay::Bases);
        let agreeing = Some(b'A');
        assert_ne!(
            differences.cell_color(agreeing, agreeing, &theme),
            everything.cell_color(agreeing, agreeing, &theme)
        );
        // A disagreement is painted either way.
        assert_eq!(
            differences.cell_color(Some(b'C'), agreeing, &theme),
            everything.cell_color(Some(b'C'), agreeing, &theme)
        );
    }

    #[test]
    fn residues_fall_into_six_classes() {
        assert_eq!(ResidueClass::of(b'A'), Some(ResidueClass::Hydrophobic));
        assert_eq!(ResidueClass::of(b'c'), Some(ResidueClass::Hydrophobic));
        assert_eq!(ResidueClass::of(b'K'), Some(ResidueClass::Positive));
        assert_eq!(ResidueClass::of(b'H'), Some(ResidueClass::Positive));
        assert_eq!(ResidueClass::of(b'E'), Some(ResidueClass::Negative));
        assert_eq!(ResidueClass::of(b'Y'), Some(ResidueClass::Polar));
        assert_eq!(ResidueClass::of(b'G'), Some(ResidueClass::Glycine));
        assert_eq!(ResidueClass::of(b'P'), Some(ResidueClass::Proline));
        assert_eq!(ResidueClass::of(b'-'), None);
        assert_eq!(ResidueClass::of(b'X'), None);
    }

    #[test]
    fn every_class_takes_a_different_palette_slot() {
        let classes = [
            ResidueClass::Hydrophobic,
            ResidueClass::Positive,
            ResidueClass::Negative,
            ResidueClass::Polar,
            ResidueClass::Glycine,
            ResidueClass::Proline,
        ];
        let mut slots: Vec<usize> = classes.iter().map(|c| c.palette_index()).collect();
        slots.sort_unstable();
        slots.dedup();
        assert_eq!(slots.len(), classes.len());
        // And they all fit the validated palette without wrapping onto a
        // colour some other class already has.
        assert!(slots
            .iter()
            .all(|slot| *slot < Theme::light().palette.len()));
    }

    #[test]
    fn a_deep_alignment_is_capped_and_says_so() {
        let rows: Vec<MsaSequence> = (0..60)
            .map(|i| MsaSequence::new(format!("s{i}"), b"ACGT".to_vec()))
            .collect();
        let track = MsaTrack::new(rows).max_rows(Some(10));
        assert_eq!(track.visible_rows(), (10, 50));

        let svg = Figure::new(Region::new("alignment", 0, 4).unwrap())
            .show_region_label(false)
            .push(track)
            .to_svg();
        assert!(svg.contains("sequences not shown"));
    }

    #[test]
    fn the_cap_can_be_lifted() {
        let rows: Vec<MsaSequence> = (0..60)
            .map(|i| MsaSequence::new(format!("s{i}"), b"ACGT".to_vec()))
            .collect();
        assert_eq!(MsaTrack::new(rows).max_rows(None).visible_rows(), (60, 0));
    }

    #[test]
    fn identical_neighbours_are_merged_into_one_rectangle() {
        // Every row identical to the consensus: each should collapse to a
        // single rectangle rather than one per column.
        let rows: Vec<MsaSequence> = (0..3)
            .map(|i| MsaSequence::new(format!("s{i}"), vec![b'A'; 400]))
            .collect();
        let svg = Figure::new(Region::new("alignment", 0, 400).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(rows).show_names(false))
            .to_svg();
        // Three rows, the page background, and the rect inside the clip path.
        assert_eq!(svg.matches("<rect").count(), 5);
    }

    #[test]
    fn a_run_breaks_where_the_alignment_disagrees() {
        let rows = vec![
            MsaSequence::new("a", vec![b'A'; 100]),
            MsaSequence::new("b", {
                let mut row = vec![b'A'; 100];
                row[50] = b'C';
                row
            }),
        ];
        let svg = Figure::new(Region::new("alignment", 0, 100).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(rows).show_names(false).show_letters(false))
            .to_svg();
        // Row one is a single run; row two is three, split by the mismatch.
        // Plus the page background and the clip path's own rect.
        assert_eq!(svg.matches("<rect").count(), 6);
    }

    #[test]
    fn names_go_in_the_axis_strip() {
        let theme = Theme::light();
        assert!(MsaTrack::new(alignment()).y_axis_width(&theme) > 0.0);
        assert_eq!(
            MsaTrack::new(alignment())
                .show_names(false)
                .y_axis_width(&theme),
            0.0
        );

        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MsaTrack::new(alignment()))
            .to_svg();
        assert!(svg.contains("CDC1551"));
    }

    #[test]
    fn an_empty_alignment_draws_nothing_and_does_not_panic() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MsaTrack::new(Vec::new()).label("none"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
        assert_eq!(MsaTrack::new(Vec::new()).visible_rows(), (0, 0));
    }

    #[test]
    fn a_row_is_named_once_and_a_cell_is_never_named() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MsaTrack::new(alignment()))
            .to_svg();
        // The consensus is ACGTTCGTAC, so H37Rv is the row that disagrees and
        // Beijing is the row with the gap.
        assert!(
            svg.contains("<title>H37Rv, 10 columns, 1 mismatch</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>CDC1551, 10 columns, 0 mismatches</title>"));
        assert!(svg.contains("<title>Beijing, 10 columns, 0 mismatches, 1 gap</title>"));
        assert!(svg.contains("<title>Erdman, 10 columns, 0 mismatches</title>"));
        // Four rows, four tooltips. Forty cells, no tooltips.
        assert_eq!(svg.matches("<title>").count(), 4);
    }

    #[test]
    fn a_count_in_a_tooltip_is_grouped_like_every_other_number() {
        // An alignment wide enough to reach four figures reads `1,200 columns`,
        // because a rule that grouped only the numbers a reader found long
        // would make the punctuation a fact about the data.
        let wide: Vec<u8> = b"ACGT".iter().cycle().take(1_200).copied().collect();
        let rows = vec![
            MsaSequence::new("ref", wide.clone()),
            MsaSequence::new("s1", wide),
        ];
        let svg = Figure::new(Region::new("aln", 0, 1_200).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(rows))
            .to_svg();
        assert!(
            svg.contains("<title>ref, 1,200 columns, 0 mismatches</title>"),
            "{svg}"
        );
        assert_eq!(svg.matches("<g").count(), svg.matches("</g>").count());
    }

    #[test]
    fn a_row_is_counted_over_the_columns_in_view() {
        let over = |start, end| {
            Figure::new(Region::new("alignment", start, end).unwrap())
                .show_region_label(false)
                .push(MsaTrack::new(alignment()))
                .to_svg()
        };
        // The one disagreement sits in column 4, so it is in the first half of
        // the alignment and not in the second.
        assert!(over(0, 5).contains("<title>H37Rv, 5 columns, 1 mismatch</title>"));
        assert!(over(5, 10).contains("<title>H37Rv, 5 columns, 0 mismatches</title>"));
    }

    #[test]
    fn a_row_that_stops_early_counts_only_the_columns_it_has() {
        let ragged = vec![
            MsaSequence::new("full", b"ACGTACGT".to_vec()),
            MsaSequence::new("short", b"ACGT".to_vec()),
        ];
        let svg = Figure::new(Region::new("alignment", 0, 8).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(ragged))
            .to_svg();
        assert!(
            svg.contains("<title>full, 8 columns, 0 mismatches</title>"),
            "{svg}"
        );
        assert!(svg.contains("<title>short, 4 columns, 0 mismatches</title>"));
    }

    #[test]
    fn only_the_rows_that_are_drawn_are_named() {
        let rows: Vec<MsaSequence> = (0..60)
            .map(|i| MsaSequence::new(format!("s{i}"), b"ACGT".to_vec()))
            .collect();
        let svg = Figure::new(Region::new("alignment", 0, 4).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(rows).max_rows(Some(10)))
            .to_svg();
        assert_eq!(svg.matches("<title>").count(), 10);
        assert!(svg.contains("<title>s9, 4 columns, 0 mismatches</title>"));
        assert!(!svg.contains("<title>s10,"), "a hidden row is not named");
    }

    #[test]
    fn letters_appear_only_when_a_column_is_wide_enough() {
        let wide = Figure::new(Region::new("alignment", 0, 10).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(alignment()).show_names(false))
            .to_svg();
        let narrow = Figure::new(Region::new("alignment", 0, 4_000).unwrap())
            .show_region_label(false)
            .push(MsaTrack::new(alignment()).show_names(false))
            .to_svg();
        // Letters are drawn rather than stretched to their column, so what
        // says a letter is there is the letter itself.
        assert!(wide.contains(">A</text>"), "no residues at ten columns");
        assert!(!narrow.contains(">A</text>"), "letters at four thousand");
        assert!(
            !wide.contains("textLength"),
            "an alignment letter keeps its own shape"
        );
    }

    #[test]
    fn an_attached_tree_reorders_rows_and_keeps_the_comparison_sequence() {
        let tree = Tree::parse_newick("(Beijing:1,(H37Rv:1,CDC1551:1):1,Erdman:1);").unwrap();
        let track = MsaTrack::new(alignment()).compare_to(1).tree(tree);
        let names: Vec<&str> = track
            .sequences()
            .iter()
            .map(|row| row.name.as_str())
            .collect();
        assert_eq!(names, ["Beijing", "H37Rv", "CDC1551", "Erdman"]);
        assert_eq!(track.compare_to, Some(2));
        assert_eq!(track.comparison(), b"ACGTTCGTAC");
        assert!(track.attached_tree().is_some());
    }

    #[test]
    fn the_tree_and_alignment_share_row_centres() {
        let tree = Tree::parse_newick("(Beijing:1,(H37Rv:1,CDC1551:1):1,Erdman:1);").unwrap();
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(MsaTrack::new(alignment()).tree(tree).tree_width(110.0))
            .to_svg();
        assert!(svg.matches("<line").count() >= 6, "{svg}");
        assert!(svg.contains("<title>Beijing, 10 columns"), "{svg}");
        assert!(!svg.contains("NaN"));
    }
}
