//! Sequence logos, including the enrichment and depletion variant.
//!
//! A logo turns a column of symbol frequencies into stacked letters whose
//! heights carry the numbers. Which numbers is [`LogoScore`], and there are
//! seven: what is here, how conserved it is, and five ways of asking how it
//! differs from a background. Where the baseline sits is the separate question
//! [`Centering`] answers.

use crate::dash::{Dash, DashFit};
use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::Theme;
use crate::track::{DrawContext, Rect, Track};

/// What quantity a symbol's height stands for.
///
/// `p` is the symbol's share of its column and `q` its share of the background.
/// The first two scores ignore the background entirely and can only go upwards.
/// The rest compare against it, which is what lets a symbol hang below the
/// baseline, and they are the ones [`Centering`] applies to.
///
/// | Score | Height | Answers |
/// |:------|:-------|:--------|
/// | [`Probability`](LogoScore::Probability) | `p` | what is here |
/// | [`InformationContent`](LogoScore::InformationContent) | `p (log2 K - H)` | how conserved is this position |
/// | [`LogOdds`](LogoScore::LogOdds) | `log2(p / q)` | how many doublings from the background |
/// | [`KullbackLeibler`](LogoScore::KullbackLeibler) | `p log2(p / q)` | how much this symbol contributes to the divergence |
/// | [`Difference`](LogoScore::Difference) | `p - q` | how many percentage points of surplus |
/// | [`Ratio`](LogoScore::Ratio) | `p / q` | how many times the background |
/// | [`OddsRatio`](LogoScore::OddsRatio) | `log2(p/(1-p)) - log2(q/(1-q))` | how the odds shift |
///
/// The set follows the scoring schemes of `Logolas`, whose `log`, `probKL`,
/// `diff`, `ratio` and `log-odds` are the five background-relative rows above.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogoScore {
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
    /// `log2(p / q)`, the doublings between a symbol and its background.
    ///
    /// With the default median [`Centering`] this is the EDLogo, and it is the
    /// loudest of the background-relative scores: a symbol that never appears
    /// is bounded only by [`LogoTrack::smoothing`], so it can dwarf everything
    /// else on the plot. That is either the point or a nuisance depending on
    /// the question, and [`LogoScore::KullbackLeibler`] is the quieter answer.
    LogOdds,
    /// `p log2(p / q)`, the symbol's own contribution to the Kullback-Leibler
    /// divergence between the column and the background, in bits.
    ///
    /// The probability weight is what makes this well behaved: a symbol that is
    /// absent has an enormous log odds but contributes nothing to the
    /// divergence, because it is not there to contribute. Rare symbols
    /// therefore stay in proportion and the column total is the divergence
    /// itself. Use it when a handful of absent symbols would otherwise flatten
    /// the rest of the figure.
    KullbackLeibler,
    /// `p - q`, the surplus in plain probability.
    ///
    /// The one score whose units need no explanation in a figure legend, and
    /// the one that says least about a rare symbol: going from 0.001 to 0.01 is
    /// a tenfold enrichment and a difference of 0.009.
    Difference,
    /// `p / q`, how many times the background a symbol runs.
    ///
    /// Asymmetric by construction: enrichment is unbounded above while
    /// depletion is squeezed into the interval below one, so a doubling looks
    /// larger than a halving. Centre it, or prefer [`LogoScore::LogOdds`],
    /// which is this score's logarithm and does not have that problem.
    Ratio,
    /// `log2(p/(1-p)) - log2(q/(1-q))`, the shift in log odds.
    ///
    /// Compares a symbol against the rest of its own column rather than against
    /// its own share, which is the right question when the symbols compete for
    /// one slot.
    OddsRatio,
}

impl LogoScore {
    /// Whether the score is measured against a background.
    ///
    /// The two that are not, [`LogoScore::Probability`] and
    /// [`LogoScore::InformationContent`], cannot go below the baseline, take no
    /// [`Centering`], and need no smoothing.
    pub fn uses_background(self) -> bool {
        !matches!(self, LogoScore::Probability | LogoScore::InformationContent)
    }

    /// Unit suffix for the scale annotation.
    pub fn unit(self) -> &'static str {
        match self {
            LogoScore::Probability | LogoScore::Difference => "",
            LogoScore::InformationContent | LogoScore::KullbackLeibler => " bits",
            LogoScore::LogOdds | LogoScore::OddsRatio => " log2",
            LogoScore::Ratio => "x",
        }
    }
}

/// Where the baseline of a background-relative logo sits.
///
/// Scores are absolute until something decides what "no change" looks like on
/// the page. Centring subtracts a quantile of the column from every symbol in
/// it, which puts the baseline at a typical symbol and lets the rest sort
/// themselves above and below.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Centering {
    /// Leave the scores alone. The baseline then means "exactly the
    /// background", which is the honest reading of an absolute score, at the
    /// cost of columns that sit entirely on one side of it.
    None,
    /// Subtract a quantile of the column, between 0 and 1.
    ///
    /// `Quantile(0.5)` is the median and the `Logolas` default. Lower values
    /// push the baseline down and more symbols above it, higher values do the
    /// reverse. Interpolates the way R's default `quantile` does.
    Quantile(f64),
}

impl Centering {
    /// The median, which is what an EDLogo uses.
    pub fn median() -> Self {
        Centering::Quantile(0.5)
    }
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
/// Heights are in the units of the [`LogoScore`] that produced them: bits for
/// [`LogoScore::InformationContent`] and [`LogoScore::KullbackLeibler`],
/// doublings for [`LogoScore::LogOdds`]. Both lists run from the baseline
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
/// use karyon::{Figure, LogoColumn, LogoScore, LogoTrack, Region};
///
/// // A motif read off an alignment, plotted as enrichment and depletion.
/// let logo = LogoTrack::from_sequences(0, &["ACGTA", "ACGTC", "ACGAA", "ACGTA"])
///     .edlogo()
///     .label("motif");
///
/// let svg = Figure::new(Region::new("motif", 0, 5).unwrap())
///     .push(logo)
///     .to_svg();
/// assert!(svg.contains("textLength"));
///
/// // Or built by hand from a position weight matrix, scored any other way.
/// let by_hand = LogoTrack::new(0, vec![LogoColumn::acgt(0.7, 0.1, 0.1, 0.1)])
///     .score(LogoScore::KullbackLeibler);
/// assert_eq!(by_hand.stacks().len(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct LogoTrack {
    start: u64,
    columns: Vec<LogoColumn>,
    score: LogoScore,
    centering: Centering,
    order: StackOrder,
    height: f64,
    label: Option<String>,
    smoothing: f64,
    alphabet_size: Option<usize>,
    background: Vec<(String, f64)>,
    colors: Vec<(String, String)>,
    min_letter_width: f64,
    show_scale: bool,
    max_extent: Option<f64>,
    stabilization: Option<Dash>,
    sample_size: Option<f64>,
}

impl LogoTrack {
    /// A logo whose `columns[i]` describes position `start + i`, 0-based.
    pub fn new(start: u64, columns: impl Into<Vec<LogoColumn>>) -> Self {
        LogoTrack {
            start,
            columns: columns.into(),
            score: LogoScore::InformationContent,
            centering: Centering::median(),
            order: StackOrder::LargestOutside,
            height: 80.0,
            label: None,
            smoothing: 0.01,
            alphabet_size: None,
            background: Vec::new(),
            colors: Vec::new(),
            min_letter_width: 4.0,
            show_scale: true,
            max_extent: None,
            stabilization: None,
            sample_size: None,
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

    /// Chooses what a symbol's height stands for.
    pub fn score(mut self, score: LogoScore) -> Self {
        self.score = score;
        self
    }

    /// Chooses where the baseline of a background-relative logo sits.
    ///
    /// Ignored by [`LogoScore::Probability`] and
    /// [`LogoScore::InformationContent`], which have no background to be
    /// centred against.
    pub fn centering(mut self, centering: Centering) -> Self {
        self.centering = centering;
        self
    }

    /// The EDLogo: log odds against the background, centred on the median.
    ///
    /// Shorthand for [`LogoScore::LogOdds`] with [`Centering::median`], which
    /// is the combination `Logolas` popularised and the reason this track has
    /// two sides at all.
    pub fn edlogo(self) -> Self {
        self.score(LogoScore::LogOdds)
            .centering(Centering::median())
    }

    /// Shrinks the columns towards the background before scoring them, by
    /// [`Dash`], with the default settings.
    ///
    /// A logo drawn from four sequences and a logo drawn from four thousand
    /// look identical without this, and one of them is mostly noise. The
    /// shrinkage asks how much of each column the sample size can support,
    /// fitting that question across every column at once, so a well sampled
    /// position barely moves while a thin one is pulled towards the background.
    ///
    /// Two things follow from turning it on. It needs **counts**, so a track
    /// built from a probability matrix must also be told its
    /// [`sample_size`](LogoTrack::sample_size) or the shrinkage will read it as
    /// a single observation and flatten it. And [`smoothing`](LogoTrack::smoothing)
    /// stops being applied, because the shrunk composition already has no
    /// zeros in it and adding a second repair on top would double count.
    ///
    /// ```
    /// use karyon::LogoTrack;
    ///
    /// // Four sequences that happen to agree look perfectly conserved.
    /// let thin = ["ACGT", "ACGT", "ACGT", "ACGT"];
    /// let raw = LogoTrack::from_sequences(0, &thin).alphabet_size(4);
    /// let stabilized = LogoTrack::from_sequences(0, &thin)
    ///     .alphabet_size(4)
    ///     .stabilize();
    ///
    /// assert!((raw.stacks()[0].up_total() - 2.0).abs() < 1e-9);
    /// assert!(stabilized.stacks()[0].up_total() < 1.0);
    /// ```
    pub fn stabilize(self) -> Self {
        self.stabilize_with(Dash::new())
    }

    /// Shrinks the columns with a [`Dash`] you configured yourself.
    pub fn stabilize_with(mut self, dash: Dash) -> Self {
        self.stabilization = Some(dash);
        self
    }

    /// Tells the shrinkage how many observations each column stands on.
    ///
    /// Only needed when the weights are not counts. A position weight matrix of
    /// probabilities carries no record of the alignment it came from, and
    /// shrinkage is entirely a question about sample size, so it has to be
    /// supplied from outside. Each column is rescaled to this total before
    /// being shrunk; the scores themselves are unaffected, since they only ever
    /// saw proportions.
    pub fn sample_size(mut self, observations: f64) -> Self {
        self.sample_size = Some(observations.max(0.0));
        self
    }

    /// The shrinkage that was applied, or `None` when it is switched off.
    ///
    /// Worth reporting next to the figure: it carries the fitted mixture
    /// weights and, through [`DashFit::shrinkage`], how far each column
    /// actually moved.
    pub fn dash_fit(&self) -> Option<DashFit> {
        let dash = self.stabilization.as_ref()?;
        let symbols = self.symbols();
        if symbols.is_empty() {
            return None;
        }
        let background = self.background_probabilities(&symbols, self.effective_alphabet_size());
        let counts: Vec<Vec<f64>> = self
            .columns
            .iter()
            .map(|column| {
                let raw: Vec<f64> = symbols.iter().map(|s| column.weight(s)).collect();
                match self.sample_size {
                    Some(size) => {
                        let total: f64 = raw.iter().sum();
                        if total > 0.0 {
                            raw.iter().map(|w| w / total * size).collect()
                        } else {
                            raw
                        }
                    }
                    None => raw,
                }
            })
            .collect();
        Some(dash.fit(&counts, &background))
    }

    /// Pins how far the tallest stack reaches, in score units.
    ///
    /// Without it each track scales to what is on screen, which is the right
    /// default alone and the wrong one in a stack of panels: two logos of
    /// different heights read as comparable when they are not. Pin both to the
    /// same number and the comparison becomes real.
    pub fn max_extent(mut self, extent: f64) -> Self {
        self.max_extent = Some(extent.abs());
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
        // Shrinkage is fitted across every column at once, so it happens here
        // rather than inside the per column closure below.
        let stabilized = self.dash_fit();

        self.columns
            .iter()
            .enumerate()
            .map(|(index, column)| {
                let pos = self.start + index as u64;
                // Only a background-relative score needs the column repaired
                // against zeros, so the two absolute scores stay exact. A
                // shrunk column needs no repair at all: its posterior is
                // already strictly positive.
                let smoothing = if self.score.uses_background() && stabilized.is_none() {
                    self.smoothing
                } else {
                    0.0
                };
                let p = match &stabilized {
                    Some(fit) => fit.posterior[index].clone(),
                    None => probabilities(column, &symbols, smoothing, alphabet),
                };
                let raw = score_column(self.score, &p, &background, alphabet);
                let heights = if self.score.uses_background() {
                    centre(&raw, self.centering)
                } else {
                    raw
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

    /// How far the axis reaches above and below the baseline.
    ///
    /// The two absolute scores get the axis their convention asks for, fixed
    /// rather than fitted: probabilities run to one, and information content
    /// runs to `log2(K)`, the way WebLogo always draws a DNA logo from zero to
    /// two bits. Fitting those to the data would make two figures of the same
    /// motif look different for no reason. The background-relative scores have
    /// no such natural ceiling, so they fit what is on screen unless
    /// [`LogoTrack::max_extent`] says otherwise.
    fn extents(&self, stacks: &[LogoStack], region: &Region) -> (f64, f64) {
        if let Some(extent) = self.max_extent {
            let down = if self.score.uses_background() {
                extent
            } else {
                0.0
            };
            return (extent, down);
        }
        match self.score {
            LogoScore::Probability => (1.0, 0.0),
            LogoScore::InformationContent => (self.max_bits(), 0.0),
            _ => {
                let visible = || stacks.iter().filter(|s| region.contains(s.pos));
                let up = visible().map(LogoStack::up_total).fold(0.0f64, f64::max);
                let down = visible().map(LogoStack::down_total).fold(0.0f64, f64::max);
                (up, down)
            }
        }
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

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        if !self.show_scale {
            return 0.0;
        }
        // The widest thing the axis can print is a negative value with a unit
        // on it, so measure that rather than guessing.
        let sample = format!("-99.9{}", self.score.unit());
        text_width(&sample, theme.font_size - 1.0) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let stacks = self.stacks();
        let (max_up, max_down) = self.extents(&stacks, ctx.region);
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
            let unit = self.score.unit();
            let size = ctx.theme.font_size - 1.0;
            self.chip(
                ctx,
                band.y + size + 1.0,
                &format!("{}{}", trim(max_up), unit),
            );
            if max_down > 0.0 {
                let text = format!("-{}{}", trim(max_down), unit);
                self.chip(ctx, band.bottom() - 2.0, &text);
            }
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

    /// An axis value, in the strip the figure reserved to the left of the plot.
    ///
    /// It used to sit inside the band on a translucent patch, which is what you
    /// do when there is nowhere else to put it. There is somewhere else now.
    fn chip(&self, ctx: &mut DrawContext<'_>, baseline: f64, text: &str) {
        if ctx.axis.w <= 0.0 {
            return;
        }
        let size = ctx.theme.font_size - 1.0;
        ctx.svg.text(
            ctx.axis.right() - 4.0,
            baseline,
            text,
            &ctx.theme.muted,
            size,
            Anchor::End,
        );
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

/// Applies one score to a column of probabilities.
fn score_column(
    score: LogoScore,
    probabilities: &[f64],
    background: &[f64],
    alphabet: usize,
) -> Vec<f64> {
    // A probability can still be exactly zero here when the score takes no
    // smoothing, so every division and logarithm is floored.
    let floor = 1e-12;
    match score {
        LogoScore::Probability => probabilities.to_vec(),
        LogoScore::InformationContent => {
            let bits = information_content(probabilities, alphabet);
            probabilities.iter().map(|p| p * bits).collect()
        }
        LogoScore::LogOdds => probabilities
            .iter()
            .zip(background)
            .map(|(p, q)| (p.max(floor) / q.max(floor)).log2())
            .collect(),
        LogoScore::KullbackLeibler => probabilities
            .iter()
            .zip(background)
            .map(|(p, q)| {
                // The limit of p log p as p goes to zero, which is what makes
                // this score ignore a symbol that is not there rather than
                // being overwhelmed by it.
                if *p <= 0.0 {
                    0.0
                } else {
                    p * (p / q.max(floor)).log2()
                }
            })
            .collect(),
        LogoScore::Difference => probabilities
            .iter()
            .zip(background)
            .map(|(p, q)| p - q)
            .collect(),
        LogoScore::Ratio => probabilities
            .iter()
            .zip(background)
            .map(|(p, q)| p / q.max(floor))
            .collect(),
        LogoScore::OddsRatio => probabilities
            .iter()
            .zip(background)
            .map(|(p, q)| log_odds(*p, floor) - log_odds(*q, floor))
            .collect(),
    }
}

/// `log2(x / (1 - x))`, clamped so that a certainty stays finite.
fn log_odds(value: f64, floor: f64) -> f64 {
    let value = value.clamp(floor, 1.0 - floor);
    (value / (1.0 - value)).log2()
}

/// Moves the baseline of a column of scores.
fn centre(scores: &[f64], centering: Centering) -> Vec<f64> {
    match centering {
        Centering::None => scores.to_vec(),
        Centering::Quantile(q) => {
            let offset = quantile(scores, q);
            scores.iter().map(|score| score - offset).collect()
        }
    }
}

/// Quantile of a slice, interpolating the way R's default `quantile` does.
///
/// At `q = 0.5` this is the median, averaging the middle pair when the count is
/// even. Non-finite values are dropped rather than allowed to poison the sort.
fn quantile(values: &[f64], q: f64) -> f64 {
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if sorted.is_empty() {
        return 0.0;
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if sorted.len() == 1 {
        return sorted[0];
    }
    let q = q.clamp(0.0, 1.0);
    let position = (sorted.len() - 1) as f64 * q;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        return sorted[lower];
    }
    let weight = position - lower as f64;
    sorted[lower] + weight * (sorted[upper] - sorted[lower])
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
        .score(LogoScore::Probability);
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
        let track = dna(vec![LogoColumn::acgt(5.0, 5.0, 5.0, 5.0)]).edlogo();
        let stack = &track.stacks()[0];
        assert!(stack.up.is_empty(), "{:?}", stack.up);
        assert!(stack.down.is_empty(), "{:?}", stack.down);
    }

    #[test]
    fn enrichment_goes_up_and_depletion_goes_down() {
        // A is over represented, T is absent, C and G sit at the median.
        let track = dna(vec![LogoColumn::acgt(70.0, 15.0, 15.0, 0.0)]).edlogo();
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
        let edlogo = dna(columns).edlogo();
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
        let counts = dna(vec![LogoColumn::acgt(350.0, 100.0, 50.0, 0.0)]).edlogo();
        let probs = dna(vec![LogoColumn::acgt(0.7, 0.2, 0.1, 0.0)]).edlogo();

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
        let loud = dna(columns.clone()).edlogo().smoothing(0.001);
        let quiet = dna(columns).edlogo().smoothing(1.0);
        assert!(loud.stacks()[0].down_total() > quiet.stacks()[0].down_total());
    }

    #[test]
    fn a_background_can_make_a_common_symbol_unremarkable() {
        // 40 percent G in a column looks enriched against uniform, but not
        // against a genome that is 40 percent G.
        let columns = vec![LogoColumn::acgt(20.0, 20.0, 40.0, 20.0)];
        let uniform = dna(columns.clone()).edlogo();
        let matched =
            dna(columns)
                .edlogo()
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
            dna(vec![LogoColumn::acgt(50.0, 30.0, 15.0, 5.0)]).score(LogoScore::Probability);
        let heights: Vec<f64> = track.stacks()[0].up.iter().map(|(_, h)| *h).collect();
        assert!(
            heights.windows(2).all(|w| w[0] <= w[1]),
            "not ascending: {heights:?}"
        );

        let inside = dna(vec![LogoColumn::acgt(50.0, 30.0, 15.0, 5.0)])
            .score(LogoScore::Probability)
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
        .score(LogoScore::Probability);
        let symbols: Vec<String> = track.stacks()[0]
            .up
            .iter()
            .map(|(s, _)| s.clone())
            .collect();
        assert!(symbols.contains(&"Ala".to_string()));
        assert_eq!(symbols.len(), 3);
    }

    #[test]
    fn kullback_leibler_totals_the_divergence_of_the_column() {
        // Uncentred, the column sum of p log2(p/q) is the KL divergence. For
        // this column against a uniform background that is
        // 0.5*1 + 0.25*0 + 0.125*(-1) + 0.125*(-1) = 0.25 bits.
        let track = dna(vec![LogoColumn::acgt(0.5, 0.25, 0.125, 0.125)])
            .score(LogoScore::KullbackLeibler)
            .centering(Centering::None)
            .smoothing(1e-9);
        let stack = &track.stacks()[0];
        let divergence = stack.up_total() - stack.down_total();
        assert!((divergence - 0.25).abs() < 1e-6, "got {divergence}");
    }

    #[test]
    fn kullback_leibler_stops_an_absent_symbol_flattening_the_rest() {
        // The claim the documentation makes, tested as a reader would see it:
        // one column is near uniform but missing a T, the other carries a real
        // signal. Under log odds the missing T dwarfs the signal; weighting by
        // probability puts them back in order.
        let columns = vec![
            LogoColumn::acgt(34.0, 33.0, 33.0, 0.0), // absent T, nothing else
            LogoColumn::acgt(50.0, 25.0, 15.0, 10.0), // an actual gradient
        ];
        let extent = |track: &LogoTrack| {
            let stacks = track.stacks();
            let of = |i: usize| stacks[i].up_total().max(stacks[i].down_total());
            (of(0), of(1))
        };

        let (absent, signal) = extent(&dna(columns.clone()).edlogo());
        assert!(
            absent > 3.0 * signal,
            "log odds should be dominated by the absent base: {absent} vs {signal}"
        );

        let (absent, signal) = extent(
            &dna(columns)
                .score(LogoScore::KullbackLeibler)
                .centering(Centering::median()),
        );
        assert!(
            signal > 3.0 * absent,
            "the divergence should favour the real signal: {signal} vs {absent}"
        );
    }

    #[test]
    fn a_uniform_column_is_flat_under_every_background_relative_score() {
        for score in [
            LogoScore::LogOdds,
            LogoScore::KullbackLeibler,
            LogoScore::Difference,
            LogoScore::Ratio,
            LogoScore::OddsRatio,
        ] {
            let track = dna(vec![LogoColumn::acgt(25.0, 25.0, 25.0, 25.0)]).score(score);
            let stack = &track.stacks()[0];
            assert!(
                stack.up.is_empty() && stack.down.is_empty(),
                "{score:?} left {stack:?}"
            );
        }
    }

    #[test]
    fn every_background_relative_score_agrees_on_the_direction() {
        // A is over represented and G under it. Whatever the units, the signs
        // have to match or the plot is lying.
        for score in [
            LogoScore::LogOdds,
            LogoScore::KullbackLeibler,
            LogoScore::Difference,
            LogoScore::Ratio,
            LogoScore::OddsRatio,
        ] {
            let track = dna(vec![LogoColumn::acgt(70.0, 15.0, 5.0, 10.0)]).score(score);
            let stack = &track.stacks()[0];
            assert!(
                stack.up.iter().any(|(s, _)| s == "A"),
                "{score:?} did not enrich A"
            );
            assert!(
                stack.down.iter().any(|(s, _)| s == "G"),
                "{score:?} did not deplete G"
            );
        }
    }

    #[test]
    fn difference_is_a_surplus_in_plain_probability() {
        let track = dna(vec![LogoColumn::acgt(0.5, 0.25, 0.25, 0.0)])
            .score(LogoScore::Difference)
            .centering(Centering::None)
            .smoothing(1e-9);
        let stack = &track.stacks()[0];
        let a = stack.up.iter().find(|(s, _)| s == "A").unwrap().1;
        assert!((a - 0.25).abs() < 1e-6, "0.5 over a 0.25 background: {a}");
        // The surpluses and the deficits have to cancel.
        assert!((stack.up_total() - stack.down_total()).abs() < 1e-6);
    }

    #[test]
    fn an_uncentred_score_puts_the_baseline_at_the_background() {
        // Every base is at or above a uniform background except T, so with no
        // centring only T hangs below the line.
        let track = dna(vec![LogoColumn::acgt(30.0, 25.0, 30.0, 15.0)])
            .score(LogoScore::LogOdds)
            .centering(Centering::None);
        let stack = &track.stacks()[0];
        assert_eq!(
            stack
                .down
                .iter()
                .map(|(s, _)| s.as_str())
                .collect::<Vec<_>>(),
            vec!["T"]
        );
        assert_eq!(stack.up.len(), 2, "C sits on the line: {:?}", stack.up);
    }

    #[test]
    fn the_centring_quantile_moves_the_baseline() {
        let column = vec![LogoColumn::acgt(40.0, 30.0, 20.0, 10.0)];
        let low = dna(column.clone())
            .score(LogoScore::LogOdds)
            .centering(Centering::Quantile(0.1));
        let high = dna(column)
            .score(LogoScore::LogOdds)
            .centering(Centering::Quantile(0.9));
        // A low quantile sits under almost everything, so almost everything is
        // enriched relative to it, and the reverse at a high one.
        assert!(low.stacks()[0].up.len() > high.stacks()[0].up.len());
        assert!(low.stacks()[0].down.len() < high.stacks()[0].down.len());
    }

    #[test]
    fn the_information_content_axis_is_fixed_at_log2_of_the_alphabet() {
        // WebLogo always draws DNA from zero to two bits, so a figure of one
        // motif is comparable with a figure of another.
        let region = Region::new("motif", 0, 1).unwrap();
        let track = dna(vec![LogoColumn::acgt(30.0, 25.0, 25.0, 20.0)]);
        let (up, down) = track.extents(&track.stacks(), &region);
        assert_eq!(up, 2.0);
        assert_eq!(down, 0.0);
        assert!(track.stacks()[0].up_total() < 0.1, "a near flat column");
    }

    #[test]
    fn the_probability_axis_is_fixed_at_one() {
        let region = Region::new("motif", 0, 1).unwrap();
        let track =
            dna(vec![LogoColumn::acgt(30.0, 25.0, 25.0, 20.0)]).score(LogoScore::Probability);
        assert_eq!(track.extents(&track.stacks(), &region), (1.0, 0.0));
    }

    #[test]
    fn a_pinned_extent_overrides_the_fitted_one() {
        let region = Region::new("motif", 0, 1).unwrap();
        let track = dna(vec![LogoColumn::acgt(70.0, 15.0, 5.0, 10.0)])
            .edlogo()
            .max_extent(4.0);
        assert_eq!(track.extents(&track.stacks(), &region), (4.0, 4.0));

        let one_sided = dna(vec![LogoColumn::acgt(70.0, 15.0, 5.0, 10.0)])
            .score(LogoScore::Probability)
            .max_extent(0.5);
        assert_eq!(one_sided.extents(&one_sided.stacks(), &region), (0.5, 0.0));
    }

    #[test]
    fn a_fitted_extent_only_looks_at_what_is_on_screen() {
        let track = dna(vec![
            LogoColumn::acgt(26.0, 25.0, 25.0, 24.0), // barely anything
            LogoColumn::acgt(97.0, 1.0, 1.0, 1.0),    // very loud
        ])
        .edlogo();
        let stacks = track.stacks();
        let both = Region::new("motif", 0, 2).unwrap();
        let first_only = Region::new("motif", 0, 1).unwrap();
        assert!(track.extents(&stacks, &both).0 > track.extents(&stacks, &first_only).0);
    }

    #[test]
    fn scores_declare_whether_they_use_a_background() {
        assert!(!LogoScore::Probability.uses_background());
        assert!(!LogoScore::InformationContent.uses_background());
        for score in [
            LogoScore::LogOdds,
            LogoScore::KullbackLeibler,
            LogoScore::Difference,
            LogoScore::Ratio,
            LogoScore::OddsRatio,
        ] {
            assert!(score.uses_background(), "{score:?}");
        }
    }

    #[test]
    fn centring_never_touches_an_absolute_score() {
        let column = vec![LogoColumn::acgt(40.0, 30.0, 20.0, 10.0)];
        for score in [LogoScore::Probability, LogoScore::InformationContent] {
            let centred = dna(column.clone())
                .score(score)
                .centering(Centering::Quantile(0.9));
            let plain = dna(column.clone()).score(score).centering(Centering::None);
            assert_eq!(centred.stacks(), plain.stacks(), "{score:?}");
            assert!(centred.stacks()[0].down.is_empty());
        }
    }

    #[test]
    fn the_edlogo_shorthand_is_log_odds_centred_on_the_median() {
        let column = vec![LogoColumn::acgt(40.0, 30.0, 20.0, 10.0)];
        let shorthand = dna(column.clone()).edlogo();
        let spelled_out = dna(column)
            .score(LogoScore::LogOdds)
            .centering(Centering::median());
        assert_eq!(shorthand.stacks(), spelled_out.stacks());
    }

    #[test]
    fn stabilization_tells_a_thin_alignment_from_a_deep_one() {
        // The same perfectly conserved column, seen four times and four
        // thousand times. Without shrinkage both are two bits of certainty.
        let thin: Vec<String> = std::iter::repeat("ACGT".to_string()).take(4).collect();
        let deep: Vec<String> = std::iter::repeat("ACGT".to_string()).take(4000).collect();

        let raw = |seqs: &[String]| {
            LogoTrack::from_sequences(0, seqs).alphabet_size(4).stacks()[0].up_total()
        };
        let shrunk = |seqs: &[String]| {
            LogoTrack::from_sequences(0, seqs)
                .alphabet_size(4)
                .stabilize()
                .stacks()[0]
                .up_total()
        };

        assert!(
            (raw(&thin) - raw(&deep)).abs() < 1e-9,
            "raw cannot tell them apart"
        );
        assert!(shrunk(&deep) > 1.9, "deep evidence should survive");
        assert!(shrunk(&thin) < 1.0, "thin evidence should not");
    }

    #[test]
    fn stabilization_leaves_a_well_sampled_column_alone() {
        let columns = vec![LogoColumn::acgt(700.0, 150.0, 100.0, 50.0)];
        let raw = dna(columns.clone()).edlogo();
        let shrunk = dna(columns).edlogo().stabilize();
        let difference = (raw.stacks()[0].up_total() - shrunk.stacks()[0].up_total()).abs();
        assert!(difference < 0.15, "moved by {difference}");
    }

    #[test]
    fn a_probability_matrix_needs_its_sample_size_to_be_shrunk_properly() {
        // The same shape as probabilities, told it came from 5 sequences and
        // from 5000. Only the sample size can distinguish them.
        let columns = vec![LogoColumn::acgt(0.7, 0.15, 0.1, 0.05)];
        let few = dna(columns.clone()).stabilize().sample_size(5.0);
        let many = dna(columns).stabilize().sample_size(5000.0);
        assert!(
            few.stacks()[0].up_total() < many.stacks()[0].up_total(),
            "five observations should carry less than five thousand"
        );
        assert!(few.dash_fit().unwrap().shrinkage(0) > many.dash_fit().unwrap().shrinkage(0));
    }

    #[test]
    fn the_fit_is_reported_rather_than_hidden() {
        let track = dna(vec![LogoColumn::acgt(7.0, 1.0, 1.0, 1.0)]).stabilize();
        let fit = track
            .dash_fit()
            .expect("stabilised track should report a fit");
        assert_eq!(fit.weights.len(), fit.concentrations.len());
        assert!((fit.weights.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert_eq!(fit.posterior.len(), 1);
        assert!(fit.shrinkage(0) > 0.0);

        let off = dna(vec![LogoColumn::acgt(7.0, 1.0, 1.0, 1.0)]);
        assert!(off.dash_fit().is_none());
    }

    #[test]
    fn a_shrunk_column_needs_no_smoothing() {
        // Smoothing exists to keep an absent symbol finite. The posterior is
        // already strictly positive, so changing the smoothing must not change
        // a stabilised figure at all.
        let columns = vec![LogoColumn::acgt(50.0, 25.0, 25.0, 0.0)];
        let a = dna(columns.clone()).edlogo().stabilize().smoothing(1e-9);
        let b = dna(columns).edlogo().stabilize().smoothing(0.5);
        assert_eq!(a.stacks(), b.stacks());
    }

    #[test]
    fn stabilization_keeps_an_absent_symbol_finite() {
        let track = dna(vec![LogoColumn::acgt(30.0, 30.0, 30.0, 0.0)])
            .edlogo()
            .stabilize();
        let stack = &track.stacks()[0];
        for (_, height) in stack.up.iter().chain(&stack.down) {
            assert!(height.is_finite() && *height < 100.0, "{height}");
        }
        assert!(stack.down.iter().any(|(s, _)| s == "T"));
    }

    #[test]
    fn stabilization_is_deterministic() {
        let columns = vec![
            LogoColumn::acgt(7.0, 1.0, 1.0, 1.0),
            LogoColumn::acgt(2.0, 2.0, 3.0, 3.0),
        ];
        let once = dna(columns.clone()).edlogo().stabilize().stacks();
        let twice = dna(columns).edlogo().stabilize().stacks();
        assert_eq!(once, twice);
    }

    #[test]
    fn stabilization_survives_an_empty_track() {
        let track = LogoTrack::new(0, Vec::new()).stabilize();
        assert!(track.dash_fit().is_none());
        assert!(track.stacks().is_empty());
    }

    #[test]
    fn the_median_quantile_handles_both_parities() {
        assert_eq!(quantile(&[1.0, 2.0, 3.0], 0.5), 2.0);
        assert_eq!(quantile(&[1.0, 2.0, 3.0, 4.0], 0.5), 2.5);
        assert_eq!(quantile(&[], 0.5), 0.0);
        assert_eq!(quantile(&[f64::NAN], 0.5), 0.0);
        assert_eq!(quantile(&[7.0], 0.5), 7.0);
    }

    #[test]
    fn quantiles_interpolate_the_way_r_does() {
        let values = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(quantile(&values, 0.0), 1.0);
        assert_eq!(quantile(&values, 1.0), 5.0);
        assert_eq!(quantile(&values, 0.25), 2.0);
        // Position 0.375 * 4 = 1.5, halfway between the second and third.
        assert_eq!(quantile(&values, 0.375), 2.5);
    }

    #[test]
    fn a_quantile_out_of_range_is_clamped_rather_than_panicking() {
        let values = [1.0, 2.0, 3.0];
        assert_eq!(quantile(&values, -5.0), 1.0);
        assert_eq!(quantile(&values, 5.0), 3.0);
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
