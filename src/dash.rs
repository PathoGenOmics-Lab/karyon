//! Dirichlet adaptive shrinkage: how much of a column to believe.
//!
//! A logo drawn from four aligned sequences and a logo drawn from four thousand
//! look identical, and one of them is mostly noise. Four sequences that happen
//! to agree at a position produce a perfectly conserved column, two full bits,
//! with no hint that the evidence is thin. That is not a plotting problem, it
//! is an estimation problem: the column proportions are estimates, and small
//! samples make bad ones.
//!
//! This module shrinks those estimates towards a background before anything is
//! drawn. Each column of counts is modelled as multinomial with a prior that is
//! a mixture of Dirichlet distributions, all centred on the background and
//! differing only in how tightly: a very concentrated component says the column
//! is background noise, a diffuse one says it is whatever the data says. A
//! column with thousands of observations overrules the prior and barely moves;
//! a column with five gets pulled most of the way home.
//!
//! The method is the one `Logolas` calls dash, from Dey, Xie and Stephens,
//! *A new sequence logo plot to highlight enrichment and depletion*,
//! BMC Bioinformatics 19:473 (2018).
//!
//! ```
//! use karyon::dash::Dash;
//!
//! // Three columns of a four letter alphabet. The first is well sampled and
//! // lopsided, the second is the same shape from a tenth of the data, the
//! // third is flat.
//! let counts = vec![
//!     vec![900.0, 40.0, 30.0, 30.0],
//!     vec![9.0, 1.0, 0.0, 0.0],
//!     vec![250.0, 250.0, 250.0, 250.0],
//! ];
//! let fit = Dash::new().fit(&counts, &[0.25; 4]);
//!
//! // The thin column is pulled towards the background harder than the thick
//! // one, even though they have the same shape.
//! assert!(fit.shrinkage(1) > fit.shrinkage(0));
//! ```
//!
//! # Counts, not proportions
//!
//! [`Dash::fit`] takes counts and it means it, and nothing in the signature
//! catches a caller who normalised first. What comes back is not an error and
//! does not look like one: a figure with a faint pattern in it, drawn from data
//! that had a strong one. That is the first thing to check when a logo comes
//! out flatter than the alignment behind it.
//!
//! # A column is judged next to its neighbours
//!
//! The mixture weights are fitted across every column at once, which is what
//! makes this empirical Bayes rather than a guess, and it means a column
//! borrows credibility from the others. The same thin column shrinks less among
//! forty that say the same kind of thing than it does on its own, and the flat
//! third column above is part of what the other two are judged against. The
//! price is a fit that a small enough set of columns could sway, and
//! [`Dash::null_weight`] is the ballast against it.

/// The fitted shrinkage and everything needed to report it.
///
/// A figure that quietly moved the numbers should be able to say by how much,
/// so the fit is returned rather than hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct DashFit {
    /// Concentration of each mixture component, largest first. The first is
    /// infinite: the component that says a column is pure background.
    pub concentrations: Vec<f64>,
    /// Fitted weight of each component, summing to one. Weight on the leading
    /// components is weight on "this data is mostly noise".
    pub weights: Vec<f64>,
    /// Shrunk composition of each column, in the order the columns arrived.
    pub posterior: Vec<Vec<f64>>,
    /// Proportions as observed, for comparison with [`DashFit::posterior`].
    pub observed: Vec<Vec<f64>>,
    /// Background every component is centred on.
    pub background: Vec<f64>,
    /// EM iterations actually run, at most [`Dash::max_iterations`].
    pub iterations: usize,
}

impl DashFit {
    /// How far column `index` moved towards the background, from 0 to 1.
    ///
    /// Zero means the data was left alone, one means it was replaced by the
    /// background entirely. Measured as total variation distance travelled
    /// over the distance there was to travel, and zero when the column already
    /// sat on the background.
    pub fn shrinkage(&self, index: usize) -> f64 {
        let (Some(observed), Some(posterior)) =
            (self.observed.get(index), self.posterior.get(index))
        else {
            return 0.0;
        };
        let moved = total_variation(observed, posterior);
        let available = total_variation(observed, &self.background);
        if available <= 1e-12 {
            return 0.0;
        }
        (moved / available).clamp(0.0, 1.0)
    }
}

fn total_variation(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum::<f64>() / 2.0
}

/// The shrinkage fitter.
///
/// The defaults are the ones `Logolas` uses and are a reasonable starting
/// point: a concentration grid spanning "identical to the background" down to
/// "flat over the simplex", and a prior on the mixture weights that leans
/// towards the first of those.
#[derive(Debug, Clone, PartialEq)]
pub struct Dash {
    concentrations: Vec<f64>,
    null_weight: f64,
    max_iterations: usize,
    tolerance: f64,
}

impl Dash {
    /// A fitter with the default concentration grid.
    pub fn new() -> Self {
        Dash {
            concentrations: vec![f64::INFINITY, 100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0],
            null_weight: 10.0,
            max_iterations: 200,
            tolerance: 1e-8,
        }
    }

    /// Sets the concentration grid.
    ///
    /// Each entry is the total pseudo-count of one Dirichlet component. Large
    /// values pin a column to the background, one leaves it flat over the
    /// simplex, and values below one push it out towards the corners, further
    /// from the background than the data alone would go. Infinity is the point
    /// mass at the background and is added if missing. Entries are sorted
    /// largest first and non-positive ones are dropped.
    pub fn concentrations(mut self, concentrations: impl Into<Vec<f64>>) -> Self {
        let mut grid: Vec<f64> = concentrations
            .into()
            .into_iter()
            .filter(|c| *c > 0.0 && !c.is_nan())
            .collect();
        if !grid.iter().any(|c| c.is_infinite()) {
            grid.push(f64::INFINITY);
        }
        grid.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        grid.dedup();
        self.concentrations = grid;
        self
    }

    /// Sets how strongly the mixture weights are pulled towards the component
    /// that calls everything background.
    ///
    /// This is the concentration of a Dirichlet prior on the weights, one on
    /// every component and this on the first. It is what stops a handful of
    /// columns from talking the fit into believing them. Set it to 1 for a
    /// flat prior and the maximum likelihood weights.
    pub fn null_weight(mut self, weight: f64) -> Self {
        self.null_weight = weight.max(1.0);
        self
    }

    /// Sets the EM iteration cap.
    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = iterations.max(1);
        self
    }

    /// Sets the EM convergence threshold on the weights.
    pub fn tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Shrinks every column of `counts` towards `background`.
    ///
    /// `counts[n][i]` is how often symbol `i` was seen in column `n`, and they
    /// really do have to be counts: the whole method turns on how much data
    /// stands behind a column, so handing it proportions summing to one is
    /// telling it there was one observation. `background` is normalised, and a
    /// zero entry is floored so that a symbol seen against a background that
    /// forbids it does not produce an infinite likelihood.
    pub fn fit(&self, counts: &[Vec<f64>], background: &[f64]) -> DashFit {
        let symbols = background.len();
        let background = normalise_background(background);
        if symbols == 0 || counts.is_empty() {
            return DashFit {
                concentrations: self.concentrations.clone(),
                weights: vec![1.0; self.concentrations.len().max(1)],
                posterior: Vec::new(),
                observed: Vec::new(),
                background,
                iterations: 0,
            };
        }

        let components = self.concentrations.len();
        let log_likelihood: Vec<Vec<f64>> = counts
            .iter()
            .map(|column| {
                self.concentrations
                    .iter()
                    .map(|c| log_marginal(column, &background, *c))
                    .collect()
            })
            .collect();

        let (weights, iterations) = self.expectation_maximisation(&log_likelihood, components);

        let posterior = counts
            .iter()
            .zip(&log_likelihood)
            .map(|(column, logs)| {
                let responsibility = responsibilities(logs, &weights);
                let mut mean = vec![0.0; symbols];
                for (k, share) in responsibility.iter().enumerate() {
                    if *share <= 0.0 {
                        continue;
                    }
                    let component = component_mean(column, &background, self.concentrations[k]);
                    for (slot, value) in mean.iter_mut().zip(&component) {
                        *slot += share * value;
                    }
                }
                renormalise(mean)
            })
            .collect();

        let observed = counts
            .iter()
            .map(|column| {
                let total: f64 = column.iter().filter(|c| c.is_finite()).sum();
                if total <= 0.0 {
                    background.clone()
                } else {
                    column.iter().map(|c| c.max(0.0) / total).collect()
                }
            })
            .collect();

        DashFit {
            concentrations: self.concentrations.clone(),
            weights,
            posterior,
            observed,
            background,
            iterations,
        }
    }

    /// Fits the mixture weights shared by every column.
    fn expectation_maximisation(
        &self,
        log_likelihood: &[Vec<f64>],
        components: usize,
    ) -> (Vec<f64>, usize) {
        // A Dirichlet prior on the weights, leaning on the null component.
        let prior: Vec<f64> = (0..components)
            .map(|k| if k == 0 { self.null_weight } else { 1.0 })
            .collect();
        let prior_excess: f64 = prior.iter().map(|p| p - 1.0).sum();

        let mut weights = vec![1.0 / components as f64; components];
        let mut iterations = 0;

        for _ in 0..self.max_iterations {
            iterations += 1;
            let mut totals = vec![0.0; components];
            for logs in log_likelihood {
                for (k, share) in responsibilities(logs, &weights).iter().enumerate() {
                    totals[k] += share;
                }
            }

            let denominator = log_likelihood.len() as f64 + prior_excess;
            let updated: Vec<f64> = totals
                .iter()
                .zip(&prior)
                .map(|(total, p)| ((total + p - 1.0) / denominator).max(0.0))
                .collect();
            let updated = renormalise(updated);

            let change: f64 = updated
                .iter()
                .zip(&weights)
                .map(|(a, b)| (a - b).abs())
                .sum();
            weights = updated;
            if change < self.tolerance {
                break;
            }
        }

        (weights, iterations)
    }
}

impl Default for Dash {
    fn default() -> Self {
        Dash::new()
    }
}

/// Posterior weight of each component for one column, given the mixture.
fn responsibilities(log_likelihood: &[f64], weights: &[f64]) -> Vec<f64> {
    let logs: Vec<f64> = log_likelihood
        .iter()
        .zip(weights)
        .map(|(log_lik, weight)| {
            if *weight <= 0.0 || !log_lik.is_finite() {
                f64::NEG_INFINITY
            } else {
                weight.ln() + log_lik
            }
        })
        .collect();

    let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        // Nothing is possible, which happens only for a degenerate background.
        // Fall back to the prior rather than producing zeros over zeros.
        return renormalise(weights.to_vec());
    }
    let unnormalised: Vec<f64> = logs.iter().map(|l| (l - max).exp()).collect();
    renormalise(unnormalised)
}

/// Log marginal likelihood of a column under one Dirichlet component.
///
/// For a finite concentration this is the Dirichlet-multinomial, dropping the
/// multinomial coefficient because it is the same for every component:
///
/// ```text
/// lnB(c) - lnB(N + c) + sum_i [lnG(x_i + c q_i) - lnG(c q_i)]
/// ```
///
/// An infinite concentration is the point mass at the background, whose
/// marginal is the plain multinomial likelihood `sum_i x_i ln q_i`.
fn log_marginal(counts: &[f64], background: &[f64], concentration: f64) -> f64 {
    let total: f64 = counts
        .iter()
        .filter(|c| c.is_finite())
        .map(|c| c.max(0.0))
        .sum();

    if concentration.is_infinite() {
        return counts
            .iter()
            .zip(background)
            .map(|(x, q)| {
                let x = if x.is_finite() { x.max(0.0) } else { 0.0 };
                if x == 0.0 {
                    0.0
                } else {
                    x * q.ln()
                }
            })
            .sum();
    }

    let mut value = ln_gamma(concentration) - ln_gamma(total + concentration);
    for (x, q) in counts.iter().zip(background) {
        let x = if x.is_finite() { x.max(0.0) } else { 0.0 };
        let alpha = concentration * q;
        value += ln_gamma(x + alpha) - ln_gamma(alpha);
    }
    value
}

/// Posterior mean composition of a column under one component.
fn component_mean(counts: &[f64], background: &[f64], concentration: f64) -> Vec<f64> {
    if concentration.is_infinite() {
        return background.to_vec();
    }
    let total: f64 = counts
        .iter()
        .filter(|c| c.is_finite())
        .map(|c| c.max(0.0))
        .sum();
    let denominator = total + concentration;
    counts
        .iter()
        .zip(background)
        .map(|(x, q)| {
            let x = if x.is_finite() { x.max(0.0) } else { 0.0 };
            (x + concentration * q) / denominator
        })
        .collect()
}

fn normalise_background(background: &[f64]) -> Vec<f64> {
    let floor = 1e-9;
    let cleaned: Vec<f64> = background
        .iter()
        .map(|q| if q.is_finite() { q.max(floor) } else { floor })
        .collect();
    renormalise(cleaned)
}

fn renormalise(values: Vec<f64>) -> Vec<f64> {
    let total: f64 = values.iter().filter(|v| v.is_finite()).sum();
    if total <= 0.0 || !total.is_finite() {
        let share = 1.0 / values.len().max(1) as f64;
        return vec![share; values.len()];
    }
    values.iter().map(|v| v.max(0.0) / total).collect()
}

/// Natural log of the gamma function, by the Lanczos approximation.
///
/// Rust has no `lgamma` in the standard library and this crate has no
/// dependencies, so here it is. Accurate to roughly fifteen digits over the
/// positive reals, which is far more than a shrinkage weight needs.
pub fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    if !x.is_finite() || x <= 0.0 {
        return f64::INFINITY;
    }
    if x < 0.5 {
        // Reflection, for the small alphas a diffuse component produces.
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }

    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() < tolerance
    }

    #[test]
    fn ln_gamma_matches_the_factorials() {
        assert!(close(ln_gamma(1.0), 0.0, 1e-12));
        assert!(close(ln_gamma(2.0), 0.0, 1e-12));
        assert!(close(ln_gamma(3.0), 2.0f64.ln(), 1e-12));
        assert!(close(ln_gamma(5.0), 24.0f64.ln(), 1e-12));
        assert!(close(ln_gamma(11.0), (3_628_800.0f64).ln(), 1e-10));
    }

    #[test]
    fn ln_gamma_handles_the_half_integers_and_small_arguments() {
        // Gamma(1/2) is the square root of pi.
        assert!(close(
            ln_gamma(0.5),
            std::f64::consts::PI.sqrt().ln(),
            1e-12
        ));
        // Gamma(1/4), the smallest alpha a uniform four letter background makes.
        assert!(close(ln_gamma(0.25), 3.625_609_908_221_908f64.ln(), 1e-10));
        assert!(ln_gamma(0.0).is_infinite());
        assert!(ln_gamma(-1.0).is_infinite());
    }

    #[test]
    fn the_recurrence_holds() {
        // ln G(x + 1) = ln x + ln G(x), across the reflection boundary too.
        for x in [0.3, 0.5, 0.9, 1.5, 4.2, 37.0, 500.0] {
            assert!(
                close(ln_gamma(x + 1.0), x.ln() + ln_gamma(x), 1e-9),
                "failed at {x}"
            );
        }
    }

    #[test]
    fn the_infinite_component_is_the_multinomial_likelihood() {
        let counts = vec![3.0, 1.0, 0.0, 0.0];
        let background = vec![0.25; 4];
        let expected = 4.0 * 0.25f64.ln();
        assert!(close(
            log_marginal(&counts, &background, f64::INFINITY),
            expected,
            1e-12
        ));
    }

    #[test]
    fn the_infinite_component_returns_the_background_unchanged() {
        let mean = component_mean(&[9.0, 1.0, 0.0, 0.0], &[0.25; 4], f64::INFINITY);
        assert_eq!(mean, vec![0.25; 4]);
    }

    #[test]
    fn a_finite_component_is_the_dirichlet_posterior_mean() {
        // (x + c q) / (N + c), by hand: (9 + 10*0.25) / (10 + 10).
        let mean = component_mean(&[9.0, 1.0, 0.0, 0.0], &[0.25; 4], 10.0);
        assert!(close(mean[0], 11.5 / 20.0, 1e-12));
        assert!(close(mean[2], 2.5 / 20.0, 1e-12));
        assert!(close(mean.iter().sum::<f64>(), 1.0, 1e-12));
    }

    #[test]
    fn a_thin_column_shrinks_further_than_a_thick_one_of_the_same_shape() {
        let counts = vec![vec![900.0, 40.0, 30.0, 30.0], vec![9.0, 0.4, 0.3, 0.3]];
        let fit = Dash::new().fit(&counts, &[0.25; 4]);
        assert!(
            fit.shrinkage(1) > fit.shrinkage(0),
            "thin {} should exceed thick {}",
            fit.shrinkage(1),
            fit.shrinkage(0)
        );
        // And the well sampled column keeps most of its shape.
        assert!(fit.posterior[0][0] > 0.8, "{:?}", fit.posterior[0]);
    }

    #[test]
    fn more_of_the_same_evidence_means_less_shrinkage() {
        let shape = [0.6, 0.2, 0.1, 0.1];
        let mut previous = f64::INFINITY;
        for scale in [5.0, 20.0, 100.0, 1000.0] {
            let counts = vec![shape.iter().map(|p| p * scale).collect::<Vec<f64>>()];
            let fit = Dash::new().fit(&counts, &[0.25; 4]);
            let shrinkage = fit.shrinkage(0);
            assert!(
                shrinkage <= previous + 1e-9,
                "shrinkage rose from {previous} to {shrinkage} at n = {scale}"
            );
            previous = shrinkage;
        }
        assert!(previous < 0.05, "a thousand observations still moved a lot");
    }

    #[test]
    fn a_posterior_is_always_a_composition() {
        let counts = vec![
            vec![10.0, 0.0, 0.0, 0.0],
            vec![0.0; 4],
            vec![1.0, 1.0, 1.0, 1.0],
            vec![f64::NAN, 5.0, -3.0, 2.0],
        ];
        let fit = Dash::new().fit(&counts, &[0.25; 4]);
        for posterior in &fit.posterior {
            assert!(close(posterior.iter().sum::<f64>(), 1.0, 1e-9));
            assert!(
                posterior.iter().all(|p| *p > 0.0 && p.is_finite()),
                "{posterior:?}"
            );
        }
    }

    #[test]
    fn an_empty_column_becomes_the_background() {
        let fit = Dash::new().fit(&[vec![0.0; 4]], &[0.4, 0.3, 0.2, 0.1]);
        for (posterior, background) in fit.posterior[0].iter().zip(&fit.background) {
            assert!(close(*posterior, *background, 1e-9));
        }
    }

    #[test]
    fn the_weights_are_a_distribution() {
        let counts = vec![vec![50.0, 30.0, 15.0, 5.0], vec![25.0; 4]];
        let fit = Dash::new().fit(&counts, &[0.25; 4]);
        assert_eq!(fit.weights.len(), fit.concentrations.len());
        assert!(close(fit.weights.iter().sum::<f64>(), 1.0, 1e-9));
        assert!(fit.weights.iter().all(|w| *w >= 0.0));
        assert!(fit.iterations >= 1);
    }

    #[test]
    fn flat_data_puts_its_weight_on_the_null_component() {
        // Twenty columns that are all background. The fit should conclude that
        // the differences it sees are noise.
        let counts = vec![vec![25.0, 25.0, 25.0, 25.0]; 20];
        let fit = Dash::new().fit(&counts, &[0.25; 4]);
        assert!(
            fit.weights[0] > 0.5,
            "null weight was only {}",
            fit.weights[0]
        );
    }

    #[test]
    fn strong_data_pulls_weight_off_the_null_component() {
        let counts = vec![vec![970.0, 10.0, 10.0, 10.0]; 20];
        let fit = Dash::new().fit(&counts, &[0.25; 4]);
        let diffuse: f64 = fit.weights[1..].iter().sum();
        assert!(diffuse > fit.weights[0], "{:?}", fit.weights);
        assert!(fit.posterior[0][0] > 0.9);
    }

    #[test]
    fn strength_is_borrowed_across_columns() {
        // The same thin column, once alone and once among columns that all say
        // the same thing. Agreement should buy it some credibility.
        let lonely = Dash::new().fit(&[vec![7.0, 1.0, 1.0, 1.0]], &[0.25; 4]);
        let crowd_counts = vec![vec![7.0, 1.0, 1.0, 1.0]; 40];
        let crowd = Dash::new().fit(&crowd_counts, &[0.25; 4]);
        assert!(
            crowd.shrinkage(0) < lonely.shrinkage(0),
            "crowd {} should be under lonely {}",
            crowd.shrinkage(0),
            lonely.shrinkage(0)
        );
    }

    #[test]
    fn a_flat_prior_shrinks_less_than_a_null_biased_one() {
        let counts = vec![vec![7.0, 1.0, 1.0, 1.0]; 5];
        let biased = Dash::new().fit(&counts, &[0.25; 4]);
        let flat = Dash::new().null_weight(1.0).fit(&counts, &[0.25; 4]);
        assert!(flat.shrinkage(0) <= biased.shrinkage(0));
    }

    #[test]
    fn the_grid_always_contains_the_null_component_and_is_ordered() {
        let dash = Dash::new().concentrations(vec![5.0, 50.0, -1.0, 5.0]);
        assert!(dash.concentrations[0].is_infinite());
        assert_eq!(dash.concentrations, vec![f64::INFINITY, 50.0, 5.0]);
    }

    #[test]
    fn a_background_can_be_something_other_than_uniform() {
        // A column matching a skewed background should barely move.
        let background = [0.4, 0.4, 0.1, 0.1];
        let counts = vec![vec![40.0, 40.0, 10.0, 10.0]; 10];
        let fit = Dash::new().fit(&counts, &background);
        assert!(fit.shrinkage(0) < 0.2);
        for (posterior, expected) in fit.posterior[0].iter().zip(&background) {
            assert!(close(*posterior, *expected, 0.02), "{:?}", fit.posterior[0]);
        }
    }

    #[test]
    fn fitting_nothing_does_not_panic() {
        let fit = Dash::new().fit(&[], &[0.25; 4]);
        assert!(fit.posterior.is_empty());
        assert_eq!(fit.iterations, 0);
        let fit = Dash::new().fit(&[vec![]], &[]);
        assert!(fit.posterior.is_empty());
    }

    #[test]
    fn the_fit_is_deterministic() {
        let counts = vec![vec![50.0, 30.0, 15.0, 5.0], vec![9.0, 1.0, 0.0, 0.0]];
        let first = Dash::new().fit(&counts, &[0.25; 4]);
        let second = Dash::new().fit(&counts, &[0.25; 4]);
        assert_eq!(first, second);
    }
}
