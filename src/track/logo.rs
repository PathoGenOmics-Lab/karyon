//! Sequence logos, including the enrichment and depletion variant.
//!
//! A logo turns a column of symbol frequencies into stacked letters whose
//! heights carry the numbers. Three scalings ship here, and they answer
//! different questions: how often, how conserved, and how different from a
//! background.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::track::{DrawContext, Rect, Track};

/// How symbol heights are derived from the weights in a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogoScaling {
    /// Every column is one unit tall and each symbol takes its share.
    ///
    /// Shows composition and nothing else: a column of four equal bases is as
    /// tall as a perfectly conserved one.
    Probability,
    /// The classic sequence logo. Column height is the information content,
    /// `log2(K) - H`, in bits, and each symbol takes its share of that.
    ///
    /// A perfectly conserved DNA column is 2 bits tall, a uniform one is flat.
    /// This is what `seqLogo` and WebLogo draw.
    InformationContent,
    /// Enrichment above the line, depletion below it.
    ///
    /// Each symbol scores `log2(p / q)` against a background `q`, and the
    /// column is then centred on its own median so that the line sits where a
    /// typical symbol would. Symbols scoring above the median stack upwards,
    /// those below stack downwards.
    ///
    /// This is the plot that `Logolas` calls an EDLogo, and its point is that a
    /// classic logo can only say "this symbol is common here". A column where
    /// one base is absent carries as much biology as one where a base is fixed,
    /// and only this scaling can draw it.
    ///
    /// Uses [`LogoTrack::smoothing`], because a symbol observed zero times has
    /// an infinitely negative score until something is added to it.
    EnrichmentDepletion,
}

/// Where the tallest symbol of a stack sits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StackOrder {
    /// Tallest symbol furthest from the baseline. The `seqLogo` and WebLogo
    /// convention, and the default.
    LargestOutside,
    /// Tallest symbol against the baseline.
    LargestInside,
}

/// The symbol weights at one position.
///
/// Weights are whatever you have: raw counts from an alignment, probabilities
/// from a position weight matrix, or expression values. Every scaling starts by
/// normalising the column, so only the ratios between symbols matter.
#[derive(Debug, Clone, PartialEq)]
pub struct LogoColumn {
    entries: Vec<(String, f64)>,
}

impl LogoColumn {
    /// A column from symbol and weight pairs.
    ///
    /// Symbols are arbitrary strings, not just single letters: amino acid
    /// codes, codons, k-mers or words all work. Negative and non-finite
    /// weights are treated as zero, and a repeated symbol accumulates.
    pub fn new<S, I>(entries: I) -> Self
    where
        S: Into<String>,
        I: IntoIterator<Item = (S, f64)>,
    {
        let mut column = LogoColumn {
            entries: Vec::new(),
        };
        for (symbol, weight) in entries {
            column.add(symbol, weight);
        }
        column
    }

    /// A DNA column, in the conventional A, C, G, T order.
    pub fn acgt(a: f64, c: f64, g: f64, t: f64) -> Self {
        LogoColumn::new([("A", a), ("C", c), ("G", g), ("T", t)])
    }

    /// Adds weight to a symbol, creating it if this is the first time.
    pub fn add(&mut self, symbol: impl Into<String>, weight: f64) {
        let weight = if weight.is_finite() && weight > 0.0 {
            weight
        } else {
            0.0
        };
        let symbol = symbol.into();
        match self.entries.iter_mut().find(|(s, _)| *s == symbol) {
            Some((_, existing)) => *existing += weight,
            None => self.entries.push((symbol, weight)),
        }
    }

    /// Weight of a symbol, zero if it is not in this column.
    pub fn weight(&self, symbol: &str) -> f64 {
        self.entries
            .iter()
            .find(|(s, _)| s == symbol)
            .map_or(0.0, |(_, w)| *w)
    }

    /// The symbols in this column, in insertion order.
    pub fn symbols(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(s, _)| s.as_str())
    }

    /// Sum of every weight.
    pub fn total(&self) -> f64 {
        self.entries.iter().map(|(_, w)| w).sum()
    }
}

/// One drawn column: the symbols going up and the symbols going down.
///
/// Heights are in the units of the scaling that produced them, bits for
/// [`LogoScaling::InformationContent`] and log2 fold change for
/// [`LogoScaling::EnrichmentDepletion`]. Both lists run from the baseline
/// outwards, so drawing them in order stacks them correctly.
#[derive(Debug, Clone, PartialEq)]
pub struct LogoStack {
    /// Position of this column, 0-based.
    pub pos: u64,
    /// Symbols above the baseline, and how tall each one is.
    pub up: Vec<(String, f64)>,
    /// Symbols below the baseline, and how far each one drops.
    pub down: Vec<(String, f64)>,
}

impl LogoStack {
    /// Total height above the baseline.
    pub fn up_total(&self) -> f64 {
        self.up.iter().map(|(_, h)| h).sum()
    }

    /// Total depth below the baseline.
    pub fn down_total(&self) -> f64 {
        self.down.iter().map(|(_, h)| h).sum()
    }
}

/// A sequence logo over consecutive positions.
///
/// ```
/// use karyon::{Figure, LogoColumn, LogoScaling, LogoTrack, Region};
///
/// // A motif read off an alignment, plotted as enrichment and depletion.
/// let logo = LogoTrack::from_sequences(0, &["ACGTA", "ACGTC", "ACGAA", "ACGTA"])
///     .scaling(LogoScaling::EnrichmentDepletion)
///     .label("motif");
///
/// let svg = Figure::new(Region::new("motif", 0, 5).unwrap())
///     .push(logo)
///     .to_svg();
/// assert!(svg.contains("textLength"));
///
/// // Or built by hand from a position weight matrix.
/// let by_hand = LogoTrack::new(0, vec![LogoColumn::acgt(0.7, 0.1, 0.1, 0.1)]);
/// assert_eq!(by_hand.stacks().len(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct LogoTrack {
    start: u64,
    columns: Vec<LogoColumn>,
    scaling: LogoScaling,
    order: StackOrder,
    height: f64,
    label: Option<String>,
    smoothing: f64,
    alphabet_size: Option<usize>,
    background: Vec<(String, f64)>,
    colors: Vec<(String, String)>,
    min_letter_width: f64,
    show_scale: bool,
}

impl LogoTrack {
    /// A logo whose `columns[i]` describes position `start + i`, 0-based.
    pub fn new(start: u64, columns: impl Into<Vec<LogoColumn>>) -> Self {
        LogoTrack {
            start,
            columns: columns.into(),
            scaling: LogoScaling::InformationContent,
            order: StackOrder::LargestOutside,
            height: 80.0,
            label: None,
            smoothing: 0.01,
            alphabet_size: None,
            background: Vec::new(),
            colors: Vec::new(),
            min_letter_width: 4.0,
            show_scale: true,
        }
    }

    /// A logo counted from aligned sequences, one column per position.
    ///
    /// Characters are upper-cased, so a soft-masked `acgt` counts as `ACGT`.
    /// Gap characters `-` and `.` are skipped rather than counted as a symbol,
    /// which means a column with gaps is normalised over the sequences that
    /// actually have a base there. Sequences may differ in length; each one
    /// contributes to as many columns as it has characters.
    pub fn from_sequences<S: AsRef<str>>(start: u64, sequences: &[S]) -> Self {
        let width = sequences
            .iter()
            .map(|s| s.as_ref().chars().count())
            .max()
            .unwrap_or(0);
        let mut columns = vec![
            LogoColumn {
                entries: Vec::new()
            };
            width
        ];
        for sequence in sequences {
            for (index, symbol) in sequence.as_ref().chars().enumerate() {
                if symbol == '-' || symbol == '.' {
                    continue;
                }
                columns[index].add(symbol.to_ascii_uppercase().to_string(), 1.0);
            }
        }
        LogoTrack::new(start, columns)
    }

    /// A logo from a weight matrix and the alphabet its rows are in.
    ///
    /// `columns[i][j]` is the weight of `alphabet[j]` at position `i`. Columns
    /// shorter than the alphabet are padded with zeros, and anything past the
    /// end of the alphabet is dropped.
    pub fn from_matrix(start: u64, alphabet: &[&str], columns: &[Vec<f64>]) -> Self {
        let built = columns
            .iter()
            .map(|weights| {
                LogoColumn::new(
                    alphabet
                        .iter()
                        .enumerate()
                        .map(|(j, symbol)| (*symbol, weights.get(j).copied().unwrap_or(0.0))),
                )
            })
            .collect::<Vec<_>>();
        let mut track = LogoTrack::new(start, built);
        track.alphabet_size = Some(alphabet.len());
        track
    }

    /// Chooses probability, information content or enrichment and depletion.
    pub fn scaling(mut self, scaling: LogoScaling) -> Self {
        self.scaling = scaling;
        self
    }

    /// Chooses which end of the stack the tallest symbol goes.
    pub fn order(mut self, order: StackOrder) -> Self {
        self.order = order;
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(4.0);
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets how much of a column's mass is spread evenly before scoring an
    /// EDLogo.
    ///
    /// This is a **fraction of the column total**, not a count, and that is
    /// deliberate: the same motif has to plot identically whether you hand it
    /// over as counts out of 500 or as a probability matrix summing to one. An
    /// absolute pseudocount cannot do that.
    ///
    /// It is also the knob that decides how loud "never observed" is. A symbol
    /// with zero weight has no finite log odds, so the smoothing stands in for
    /// the observation that has not happened yet, and it bounds how far that
    /// symbol can fall at roughly `log2(smoothing / K / q)`. With the default
    /// of 0.01 over a uniform four letter alphabet, an absent base bottoms out
    /// near -6.6, which keeps it on the page next to scores of one or two.
    /// Lower it to let absence speak louder, raise it to compress the plot.
    ///
    /// Ignored by the other two scalings, which need no such repair.
    pub fn smoothing(mut self, fraction: f64) -> Self {
        self.smoothing = fraction.max(1e-9);
        self
    }

    /// Sets how many symbols the alphabet has, including unobserved ones.
    ///
    /// It matters twice: information content is measured against `log2(K)`, and
    /// a uniform background is `1/K`. Left alone, `K` is the number of distinct
    /// symbols that actually appear, which is wrong for a DNA motif where one
    /// base never shows up. Set it to 4 for DNA and 20 for protein.
    pub fn alphabet_size(mut self, size: usize) -> Self {
        self.alphabet_size = Some(size);
        self
    }

    /// Sets the background frequencies an EDLogo scores against.
    ///
    /// Defaults to uniform. Give it the composition of your genome when that
    /// is the comparison you mean: in a GC-poor organism, a G is a bigger
    /// surprise than a uniform background admits. Values are normalised, and
    /// symbols left out fall back to uniform.
    pub fn background<S, I>(mut self, background: I) -> Self
    where
        S: Into<String>,
        I: IntoIterator<Item = (S, f64)>,
    {
        self.background = background
            .into_iter()
            .map(|(symbol, weight)| (symbol.into(), weight))
            .collect();
        self
    }

    /// Overrides the colour of one symbol.
    pub fn symbol_color(mut self, symbol: impl Into<String>, color: impl Into<String>) -> Self {
        self.colors.push((symbol.into(), color.into()));
        self
    }

    /// Sets how many pixels a column needs before letters are drawn.
    ///
    /// Below this the column becomes a stacked bar, which still carries the
    /// heights when a letter would be an illegible sliver.
    pub fn min_letter_width(mut self, pixels: f64) -> Self {
        self.min_letter_width = pixels;
        self
    }

    /// Shows or hides the scale annotation in the corner.
    pub fn show_scale(mut self, show: bool) -> Self {
        self.show_scale = show;
        self
    }

    /// Every symbol in the logo, in order of first appearance.
    ///
    /// The order is what fixes the colours, so adding a column of new symbols
    /// never recolours the ones already there.
    pub fn symbols(&self) -> Vec<String> {
        let mut seen: Vec<String> = Vec::new();
        for column in &self.columns {
            for symbol in column.symbols() {
                if !seen.iter().any(|s| s == symbol) {
                    seen.push(symbol.to_string());
                }
            }
        }
        for (symbol, _) in &self.background {
            if !seen.iter().any(|s| s == symbol) {
                seen.push(symbol.clone());
            }
        }
        seen
    }

    /// Size of the alphabet used for information content and uniform
    /// backgrounds.
    pub fn effective_alphabet_size(&self) -> usize {
        let observed = self.symbols().len().max(1);
        self.alphabet_size.unwrap_or(observed).max(observed)
    }

    /// The computed stacks, one per column.
    ///
    /// Exposed because the numbers are worth checking, and because a caller may
    /// want them for a table rather than a picture.
    pub fn stacks(&self) -> Vec<LogoStack> {
        let symbols = self.symbols();
        let alphabet = self.effective_alphabet_size();
        let background = self.background_probabilities(&symbols, alphabet);

        self.columns
            .iter()
            .enumerate()
            .map(|(index, column)| {
                let pos = self.start + index as u64;
                let heights = match self.scaling {
                    LogoScaling::Probability => probabilities(column, &symbols, 0.0, alphabet),
                    LogoScaling::InformationContent => {
                        let p = probabilities(column, &symbols, 0.0, alphabet);
                        let bits = information_content(&p, alphabet);
                        p.iter().map(|value| value * bits).collect()
                    }
                    LogoScaling::EnrichmentDepletion => {
                        let p = probabilities(column, &symbols, self.smoothing, alphabet);
                        median_centred_log_odds(&p, &background)
                    }
                };

                let mut up: Vec<(String, f64)> = Vec::new();
                let mut down: Vec<(String, f64)> = Vec::new();
                for (symbol, height) in symbols.iter().zip(&heights) {
                    if !height.is_finite() || height.abs() < 1e-9 {
                        continue;
                    }
                    if *height > 0.0 {
                        up.push((symbol.clone(), *height));
                    } else {
                        down.push((symbol.clone(), -*height));
                    }
                }
                self.sort_stack(&mut up);
                self.sort_stack(&mut down);
                LogoStack { pos, up, down }
            })
            .collect()
    }

    /// Orders a stack so that drawing it from the baseline outwards is correct.
    fn sort_stack(&self, stack: &mut [(String, f64)]) {
        stack.sort_by(|a, b| match self.order {
            // Ties broken on the symbol so the output stays deterministic.
            StackOrder::LargestOutside => {
                a.1.partial_cmp(&b.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.0.cmp(&b.0))
            }
            StackOrder::LargestInside => {
                b.1.partial_cmp(&a.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.0.cmp(&b.0))
            }
        });
    }

    /// Background probability per symbol, uniform unless one was given.
    fn background_probabilities(&self, symbols: &[String], alphabet: usize) -> Vec<f64> {
        let uniform = 1.0 / alphabet as f64;
        if self.background.is_empty() {
            return vec![uniform; symbols.len()];
        }
        let total: f64 = self
            .background
            .iter()
            .map(|(_, w)| if w.is_finite() && *w > 0.0 { *w } else { 0.0 })
            .sum();
        if total <= 0.0 {
            return vec![uniform; symbols.len()];
        }
        symbols
            .iter()
            .map(|symbol| {
                self.background
                    .iter()
                    .find(|(s, _)| s == symbol)
                    .map_or(uniform, |(_, w)| (w / total).max(1e-12))
            })
            .collect()
    }

    /// The largest extent above and below the baseline across every column.
    fn extents(&self, stacks: &[LogoStack]) -> (f64, f64) {
        let up = stacks
            .iter()
            .map(LogoStack::up_total)
            .fold(0.0f64, f64::max);
        let down = stacks
            .iter()
            .map(LogoStack::down_total)
            .fold(0.0f64, f64::max);
        (up, down)
    }

    fn color_of<'a>(&'a self, symbol: &str, index: usize, ctx: &'a DrawContext<'_>) -> String {
        if let Some((_, color)) = self.colors.iter().find(|(s, _)| s == symbol) {
            return color.clone();
        }
        let bytes = symbol.as_bytes();
        if bytes.len() == 1
            && matches!(
                bytes[0].to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'U'
            )
        {
            return ctx.theme.bases.of(bytes[0]).to_string();
        }
        ctx.theme.color(index).to_string()
    }
}

impl Track for LogoTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let stacks = self.stacks();
        let (max_up, max_down) = self.extents(&stacks);
        let span = max_up + max_down;
        if span <= 0.0 {
            ctx.svg.line(
                band.x,
                band.bottom(),
                band.right(),
                band.bottom(),
                &ctx.theme.rule,
                1.0,
            );
            return;
        }

        // One scale for both halves, with the baseline placed so that the
        // deepest depletion and the tallest enrichment both just fit. When
        // nothing is depleted this puts the baseline on the floor, which is
        // exactly a classic logo.
        let units_per_px = band.h / span;
        let baseline_y = band.y + max_up * units_per_px;
        ctx.svg.line(
            band.x,
            baseline_y,
            band.right(),
            baseline_y,
            &ctx.theme.rule,
            1.0,
        );

        let symbols = self.symbols();
        let index_of = |symbol: &str| symbols.iter().position(|s| s == symbol).unwrap_or(0);

        for stack in &stacks {
            if !ctx.region.contains(stack.pos) {
                continue;
            }
            let left = ctx.scale.x(stack.pos);
            let right = ctx.scale.x(stack.pos + 1);
            let padding = ((right - left) * 0.06).min(1.0);
            let x = left + padding;
            let width = (right - left) - 2.0 * padding;
            if width <= 0.0 {
                continue;
            }

            let mut cursor = baseline_y;
            for (symbol, height) in &stack.up {
                let pixels = height * units_per_px;
                let color = self.color_of(symbol, index_of(symbol), ctx);
                let cell = Rect {
                    x,
                    y: cursor - pixels,
                    w: width,
                    h: pixels,
                };
                self.draw_symbol(ctx, cell, symbol, &color);
                cursor -= pixels;
            }

            let mut cursor = baseline_y;
            for (symbol, height) in &stack.down {
                let pixels = height * units_per_px;
                let color = self.color_of(symbol, index_of(symbol), ctx);
                let cell = Rect {
                    x,
                    y: cursor,
                    w: width,
                    h: pixels,
                };
                self.draw_symbol(ctx, cell, symbol, &color);
                cursor += pixels;
            }
        }

        if self.show_scale {
            let text = self.scale_label(max_up.max(max_down));
            let size = ctx.theme.font_size - 1.0;
            ctx.svg.rect_opacity(
                band.x + 1.0,
                band.y,
                text_width(&text, size) + 5.0,
                size + 4.0,
                &ctx.theme.background,
                0.72,
            );
            ctx.svg.text(
                band.x + 3.0,
                band.y + size + 1.0,
                &text,
                &ctx.theme.muted,
                size,
                Anchor::Start,
            );
        }
    }
}

impl LogoTrack {
    /// Draws one symbol stretched to fill `cell`.
    fn draw_symbol(&self, ctx: &mut DrawContext<'_>, cell: Rect, symbol: &str, color: &str) {
        if cell.h <= 0.2 {
            return;
        }
        if cell.w < self.min_letter_width {
            // Too narrow to read a letter, so keep the height and drop the
            // glyph rather than drawing an illegible sliver.
            ctx.svg.rect(cell.x, cell.y, cell.w, cell.h, color);
            return;
        }
        let font_size = cell.h / ctx.theme.cap_height_ratio.max(0.1);
        ctx.svg
            .glyph(cell.x, cell.bottom(), cell.w, font_size, symbol, color);
    }

    fn scale_label(&self, extent: f64) -> String {
        match self.scaling {
            LogoScaling::Probability => "1".to_string(),
            LogoScaling::InformationContent => {
                format!("{} bits", trim(self.max_bits()))
            }
            LogoScaling::EnrichmentDepletion => format!("{} log2", trim(extent)),
        }
    }

    /// Information content of a perfectly conserved column, in bits.
    pub fn max_bits(&self) -> f64 {
        (self.effective_alphabet_size() as f64).log2()
    }
}

/// Normalises a column to probabilities over `symbols`.
///
/// `smoothing` is a fraction of the column total, spread evenly across the
/// whole alphabet including the symbols this column never saw. Expressing it
/// as a fraction is what makes the result identical for counts and for a
/// probability matrix of the same shape.
fn probabilities(
    column: &LogoColumn,
    symbols: &[String],
    smoothing: f64,
    alphabet: usize,
) -> Vec<f64> {
    let observed = column.total();
    if observed <= 0.0 {
        return vec![0.0; symbols.len()];
    }
    let per_symbol = smoothing * observed / alphabet as f64;
    let total = observed * (1.0 + smoothing);
    symbols
        .iter()
        .map(|symbol| (column.weight(symbol) + per_symbol) / total)
        .collect()
}

/// `log2(K) - H`, the bits a column carries over a uniform one.
fn information_content(probabilities: &[f64], alphabet: usize) -> f64 {
    let entropy: f64 = probabilities
        .iter()
        .filter(|p| **p > 0.0)
        .map(|p| -p * p.log2())
        .sum();
    ((alphabet as f64).log2() - entropy).max(0.0)
}

/// Log odds against a background, recentred on the column median.
fn median_centred_log_odds(probabilities: &[f64], background: &[f64]) -> Vec<f64> {
    let scores: Vec<f64> = probabilities
        .iter()
        .zip(background)
        .map(|(p, q)| (p.max(1e-12) / q.max(1e-12)).log2())
        .collect();
    let centre = median(&scores);
    scores.iter().map(|score| score - centre).collect()
}

/// Median of a slice, averaging the middle pair when the count is even.
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if sorted.is_empty() {
        return 0.0;
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let middle = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    }
}

/// Two decimals at most, and none at all when they would be zeros.
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

    fn dna(columns: Vec<LogoColumn>) -> LogoTrack {
        LogoTrack::new(0, columns).alphabet_size(4)
    }

    #[test]
    fn a_column_accumulates_repeated_symbols() {
        let column = LogoColumn::new([("A", 1.0), ("C", 2.0), ("A", 3.0)]);
        assert_eq!(column.weight("A"), 4.0);
        assert_eq!(column.weight("C"), 2.0);
        assert_eq!(column.weight("G"), 0.0);
        assert_eq!(column.total(), 6.0);
    }

    #[test]
    fn negative_and_non_finite_weights_become_zero() {
        let column = LogoColumn::new([("A", -5.0), ("C", f64::NAN), ("G", f64::INFINITY)]);
        assert_eq!(column.total(), 0.0);
    }

    #[test]
    fn probability_scaling_makes_every_column_one_unit_tall() {
        let track = dna(vec![
            LogoColumn::acgt(1.0, 1.0, 1.0, 1.0),
            LogoColumn::acgt(10.0, 0.0, 0.0, 0.0),
        ])
        .scaling(LogoScaling::Probability);
        for stack in track.stacks() {
            assert!((stack.up_total() - 1.0).abs() < 1e-9);
            assert_eq!(stack.down_total(), 0.0);
        }
    }

    #[test]
    fn a_conserved_dna_column_is_two_bits_and_a_uniform_one_is_flat() {
        let track = dna(vec![
            LogoColumn::acgt(10.0, 0.0, 0.0, 0.0),
            LogoColumn::acgt(5.0, 5.0, 5.0, 5.0),
        ]);
        let stacks = track.stacks();
        assert!((stacks[0].up_total() - 2.0).abs() < 1e-9);
        assert!(stacks[1].up_total() < 1e-9);
    }

    #[test]
    fn information_content_uses_the_declared_alphabet_not_the_observed_one() {
        // Only A and C ever appear. Over a 4 letter alphabet a pure A column
        // carries 2 bits; over a 2 letter one it would carry 1.
        let columns = vec![LogoColumn::new([("A", 10.0), ("C", 0.0)])];
        let declared = LogoTrack::new(0, columns.clone()).alphabet_size(4);
        let inferred = LogoTrack::new(0, columns);
        assert!((declared.stacks()[0].up_total() - 2.0).abs() < 1e-9);
        assert!((inferred.stacks()[0].up_total() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_uniform_column_has_nothing_to_enrich_or_deplete() {
        let track = dna(vec![LogoColumn::acgt(5.0, 5.0, 5.0, 5.0)])
            .scaling(LogoScaling::EnrichmentDepletion);
        let stack = &track.stacks()[0];
        assert!(stack.up.is_empty(), "{:?}", stack.up);
        assert!(stack.down.is_empty(), "{:?}", stack.down);
    }

    #[test]
    fn enrichment_goes_up_and_depletion_goes_down() {
        // A is over represented, T is absent, C and G sit at the median.
        let track = dna(vec![LogoColumn::acgt(70.0, 15.0, 15.0, 0.0)])
            .scaling(LogoScaling::EnrichmentDepletion);
        let stack = &track.stacks()[0];
        let up: Vec<&str> = stack.up.iter().map(|(s, _)| s.as_str()).collect();
        let down: Vec<&str> = stack.down.iter().map(|(s, _)| s.as_str()).collect();
        assert_eq!(up, vec!["A"]);
        assert_eq!(down, vec!["T"]);
    }

    #[test]
    fn an_absent_symbol_is_the_thing_a_classic_logo_cannot_show() {
        let columns = vec![LogoColumn::acgt(34.0, 33.0, 33.0, 0.0)];
        let classic = dna(columns.clone());
        let edlogo = dna(columns).scaling(LogoScaling::EnrichmentDepletion);
        // Near uniform over three bases: almost no information content.
        assert!(classic.stacks()[0].up_total() < 0.45);
        // But the missing T is unmistakable once it can go below the line.
        assert!(edlogo.stacks()[0].down_total() > 3.0);
    }

    #[test]
    fn counts_and_probabilities_of_the_same_shape_plot_identically() {
        // The reason smoothing is a fraction rather than a pseudocount: the
        // same motif handed over as counts out of 500 and as a probability
        // matrix has to produce the same figure.
        let counts = dna(vec![LogoColumn::acgt(350.0, 100.0, 50.0, 0.0)])
            .scaling(LogoScaling::EnrichmentDepletion);
        let probs = dna(vec![LogoColumn::acgt(0.7, 0.2, 0.1, 0.0)])
            .scaling(LogoScaling::EnrichmentDepletion);

        let a = &counts.stacks()[0];
        let b = &probs.stacks()[0];
        assert_eq!(a.up.len(), b.up.len());
        assert_eq!(a.down.len(), b.down.len());
        for ((sa, ha), (sb, hb)) in a.up.iter().chain(&a.down).zip(b.up.iter().chain(&b.down)) {
            assert_eq!(sa, sb);
            assert!((ha - hb).abs() < 1e-9, "{sa}: {ha} vs {hb}");
        }
    }

    #[test]
    fn less_smoothing_pushes_an_absent_symbol_further_down() {
        let columns = vec![LogoColumn::acgt(30.0, 30.0, 30.0, 0.0)];
        let loud = dna(columns.clone())
            .scaling(LogoScaling::EnrichmentDepletion)
            .smoothing(0.001);
        let quiet = dna(columns)
            .scaling(LogoScaling::EnrichmentDepletion)
            .smoothing(1.0);
        assert!(loud.stacks()[0].down_total() > quiet.stacks()[0].down_total());
    }

    #[test]
    fn a_background_can_make_a_common_symbol_unremarkable() {
        // 40 percent G in a column looks enriched against uniform, but not
        // against a genome that is 40 percent G.
        let columns = vec![LogoColumn::acgt(20.0, 20.0, 40.0, 20.0)];
        let uniform = dna(columns.clone()).scaling(LogoScaling::EnrichmentDepletion);
        let matched = dna(columns)
            .scaling(LogoScaling::EnrichmentDepletion)
            .background([("A", 0.2), ("C", 0.2), ("G", 0.4), ("T", 0.2)]);

        let g_uniform = uniform.stacks()[0]
            .up
            .iter()
            .find(|(s, _)| s == "G")
            .map_or(0.0, |(_, h)| *h);
        let g_matched = matched.stacks()[0]
            .up
            .iter()
            .find(|(s, _)| s == "G")
            .map_or(0.0, |(_, h)| *h);
        assert!(g_uniform > 0.9, "expected enrichment, got {g_uniform}");
        assert!(g_matched < 1e-6, "expected no enrichment, got {g_matched}");
    }

    #[test]
    fn stacks_run_from_the_baseline_outwards() {
        let track =
            dna(vec![LogoColumn::acgt(50.0, 30.0, 15.0, 5.0)]).scaling(LogoScaling::Probability);
        let heights: Vec<f64> = track.stacks()[0].up.iter().map(|(_, h)| *h).collect();
        assert!(
            heights.windows(2).all(|w| w[0] <= w[1]),
            "not ascending: {heights:?}"
        );

        let inside = dna(vec![LogoColumn::acgt(50.0, 30.0, 15.0, 5.0)])
            .scaling(LogoScaling::Probability)
            .order(StackOrder::LargestInside);
        let heights: Vec<f64> = inside.stacks()[0].up.iter().map(|(_, h)| *h).collect();
        assert!(
            heights.windows(2).all(|w| w[0] >= w[1]),
            "not descending: {heights:?}"
        );
    }

    #[test]
    fn counting_an_alignment_matches_counting_by_hand() {
        let track = LogoTrack::from_sequences(0, &["ACGT", "ACGA", "ACGA"]);
        let stacks = track.stacks();
        assert_eq!(stacks.len(), 4);
        assert_eq!(track.columns[3].weight("A"), 2.0);
        assert_eq!(track.columns[3].weight("T"), 1.0);
        assert_eq!(track.columns[0].weight("A"), 3.0);
    }

    #[test]
    fn counting_folds_case_and_skips_gaps() {
        let track = LogoTrack::from_sequences(0, &["acgt", "ACG-", "ACG."]);
        assert_eq!(track.columns[0].weight("A"), 3.0);
        assert_eq!(track.columns[3].weight("T"), 1.0);
        assert_eq!(track.columns[3].total(), 1.0);
        assert!(!track.symbols().iter().any(|s| s == "-"));
    }

    #[test]
    fn ragged_alignments_do_not_panic() {
        let track = LogoTrack::from_sequences(0, &["ACGTACGT", "AC", ""]);
        assert_eq!(track.stacks().len(), 8);
        assert_eq!(track.columns[0].total(), 2.0);
        assert_eq!(track.columns[7].total(), 1.0);
    }

    #[test]
    fn an_empty_alignment_produces_an_empty_logo() {
        let track = LogoTrack::from_sequences(0, &[] as &[&str]);
        assert!(track.stacks().is_empty());
        assert_eq!(track.effective_alphabet_size(), 1);
    }

    #[test]
    fn a_matrix_keeps_the_alphabet_it_was_given() {
        let track = LogoTrack::from_matrix(
            10,
            &["A", "C", "G", "T"],
            &[vec![0.7, 0.1, 0.1, 0.1], vec![0.25; 4]],
        );
        assert_eq!(track.effective_alphabet_size(), 4);
        assert_eq!(track.stacks()[0].pos, 10);
        assert_eq!(track.stacks()[1].pos, 11);
        assert!(track.stacks()[1].up_total() < 1e-9);
    }

    #[test]
    fn arbitrary_symbols_are_not_limited_to_single_letters() {
        let track = LogoTrack::new(
            0,
            vec![LogoColumn::new([("Ala", 5.0), ("Gly", 3.0), ("Trp", 1.0)])],
        )
        .scaling(LogoScaling::Probability);
        let symbols: Vec<String> = track.stacks()[0]
            .up
            .iter()
            .map(|(s, _)| s.clone())
            .collect();
        assert!(symbols.contains(&"Ala".to_string()));
        assert_eq!(symbols.len(), 3);
    }

    #[test]
    fn median_handles_both_parities() {
        assert_eq!(median(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[]), 0.0);
        assert_eq!(median(&[f64::NAN]), 0.0);
    }

    #[test]
    fn letters_are_drawn_as_stretched_glyphs() {
        let svg = Figure::new(Region::new("motif", 0, 6).unwrap())
            .show_region_label(false)
            .push(dna(vec![LogoColumn::acgt(9.0, 1.0, 0.0, 0.0); 6]))
            .to_svg();
        assert!(svg.contains("textLength"));
        assert!(svg.contains("lengthAdjust=\"spacingAndGlyphs\""));
        assert!(svg.contains(">A</text>"));
    }

    #[test]
    fn narrow_columns_fall_back_to_bars() {
        // Six thousand columns across 900 pixels leaves no room for a letter.
        let columns = vec![LogoColumn::acgt(9.0, 1.0, 0.0, 0.0); 6_000];
        let svg = Figure::new(Region::new("motif", 0, 6_000).unwrap())
            .show_region_label(false)
            .push(dna(columns))
            .to_svg();
        assert!(!svg.contains("textLength"));
        assert!(svg.contains("<rect"));
    }

    #[test]
    fn an_empty_logo_still_draws_a_baseline() {
        let svg = Figure::new(Region::new("motif", 0, 4).unwrap())
            .show_region_label(false)
            .push(LogoTrack::new(0, Vec::new()))
            .to_svg();
        assert!(svg.contains("<line"));
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn columns_outside_the_region_are_not_drawn() {
        let svg = Figure::new(Region::new("motif", 0, 4).unwrap())
            .show_region_label(false)
            .push(dna(vec![LogoColumn::acgt(9.0, 1.0, 0.0, 0.0)]).label("logo"))
            .to_svg();
        let drawn = svg.matches("</text>").count();

        let shifted = Figure::new(Region::new("motif", 0, 4).unwrap())
            .show_region_label(false)
            .push(
                LogoTrack::new(500, vec![LogoColumn::acgt(9.0, 1.0, 0.0, 0.0)])
                    .alphabet_size(4)
                    .label("logo"),
            )
            .to_svg();
        assert!(shifted.matches("</text>").count() < drawn);
    }
}
