//! The coordinate ruler.
//!
//! Every other track draws data. This one draws the coordinate that the data is
//! read against, and takes no data of its own: everything it needs is the
//! region already on display. What it has to get right is turning that region
//! into a handful of numbers a reader can use without translating them first.
//!
//! # The labels are 1-based on purpose
//!
//! Positions are 0-based and half-open everywhere else in the crate, and the
//! ruler is the one place where that choice would reach a reader. A tick
//! reading `761,100` has to be the coordinate that goes into a browser search
//! box or a samtools region string, so tick positions are chosen in 1-based
//! space and converted back to `pos - 1` on the way to the [`Scale`].
//!
//! Whether a tick belongs on the boundary between two bases or in the middle of
//! one is the caller's to settle, since only the caller knows how far the
//! figure is zoomed in: [`AxisTrack::center_on_bases`] moves it.
//!
//! # One step and one unit for the whole ruler
//!
//! [`AxisTrack::tick_spacing`] asks for a density and not for a step, and the
//! step it settles on is rounded until the labels come out round. The same
//! rounding chooses which codon numbers [`CodonTrack`](crate::CodonTrack)
//! prints, which is why it is defined in this module and not in that one.
//!
//! The unit is picked from the largest label and the step together, because
//! either one alone gets it wrong. Whichever unit wins is then used for every
//! label on that ruler, since an axis that changes unit half way across has to
//! be decoded rather than read.

use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::track::{DrawContext, Track};

/// A ruler of genomic coordinates, with tick spacing chosen for the zoom level.
///
/// Ticks land on round 1-based coordinates, the numbers a reader would type
/// into a genome browser, and are labelled in bp, kb or Mb depending on how
/// much sequence is on screen. Put one at the bottom of a figure, or at the top
/// and the bottom of a tall one.
///
/// ```
/// use karyon::{AxisTrack, Figure, Region};
///
/// let svg = Figure::new(Region::parse("chr1:1-10000").unwrap())
///     .push(AxisTrack::new())
///     .to_svg();
/// assert!(svg.contains("kb"));
/// ```
#[derive(Debug, Clone)]
pub struct AxisTrack {
    height: f64,
    tick_length: f64,
    target_spacing: f64,
    label: Option<String>,
    center_on_bases: bool,
}

impl AxisTrack {
    /// A ruler with the default height of 26 pixels.
    pub fn new() -> Self {
        AxisTrack {
            height: 26.0,
            tick_length: 5.0,
            target_spacing: 110.0,
            label: None,
            center_on_bases: false,
        }
    }

    /// Puts each tick in the middle of its base rather than on its left edge.
    ///
    /// A ruler marks boundaries, which is right when a base is a fraction of a
    /// pixel wide. Once a base is a column you can see, as in a sequence logo
    /// or a short motif, the number belongs under the column it counts rather
    /// than on the line between two of them.
    pub fn center_on_bases(mut self, center: bool) -> Self {
        self.center_on_bases = center;
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(1.0);
        self
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the target distance between ticks in pixels.
    ///
    /// This is a target, not a rule: the tick step is rounded to 1, 2 or 5
    /// times a power of ten so the labels stay round numbers.
    pub fn tick_spacing(mut self, pixels: f64) -> Self {
        self.target_spacing = pixels.max(20.0);
        self
    }
}

impl Default for AxisTrack {
    fn default() -> Self {
        AxisTrack::new()
    }
}

impl Track for AxisTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let rule_y = band.y + 0.5;
        let font = ctx.theme.font_size;
        let tick_length = self.tick_length * ctx.visual_scale;
        let label_y = rule_y + tick_length + font;

        ctx.svg.line(
            band.x,
            rule_y,
            band.right(),
            rule_y,
            &ctx.theme.foreground,
            ctx.theme.tokens.stroke,
        );

        let ticks = self.ticks(ctx.region.display_start(), ctx.region.display_end(), band.w);
        // One unit for the whole ruler: an axis that switches from kb to Mb
        // half way across is unreadable.
        let unit = tick_unit(ticks.positions.last().copied().unwrap_or(0), ticks.step);
        for pos in ticks.positions {
            // Ticks are computed in 1-based space so their labels are round
            // numbers; the scale works in 0-based space.
            let x = if self.center_on_bases {
                ctx.scale.x_center(pos - 1)
            } else {
                ctx.scale.x(pos - 1)
            };
            ctx.svg.line(
                x,
                rule_y,
                x,
                rule_y + tick_length,
                &ctx.theme.foreground,
                ctx.theme.tokens.stroke,
            );

            let text = unit.format(pos);
            // Keep the first and last labels from hanging off the figure.
            let half = text_width(&text, font) / 2.0;
            let (anchor, tx) = if x - half < band.x {
                (Anchor::Start, band.x)
            } else if x + half > band.right() {
                (Anchor::End, band.right())
            } else {
                (Anchor::Middle, x)
            };
            ctx.svg
                .text(tx, label_y, &text, &ctx.theme.muted, font, anchor);
        }
    }
}

struct Ticks {
    positions: Vec<u64>,
    step: u64,
}

impl AxisTrack {
    /// Tick positions as 1-based inclusive coordinates.
    fn ticks(&self, first: u64, last: u64, width: f64) -> Ticks {
        let span = (last - first + 1) as f64;
        let target = (width / self.target_spacing).max(2.0);
        let step = nice_step(span, target);

        let mut positions = Vec::new();
        // First multiple of step at or after `first`, never 0.
        let mut pos = first.div_ceil(step) * step;
        if pos == 0 {
            pos = step;
        }
        while pos <= last && positions.len() < 512 {
            positions.push(pos);
            // A region running to the top of the coordinate range has a last
            // tick with no next one. Unchecked, the step wrapped and the ruler
            // carried on with hundreds of labels counting backwards.
            pos = match pos.checked_add(step) {
                Some(next) => next,
                None => break,
            };
        }
        // A region narrower than one step gets its own edges labelled instead
        // of nothing at all.
        if positions.len() < 2 {
            positions = vec![first, last];
            positions.dedup();
        }
        Ticks { positions, step }
    }
}

/// Rounds a raw tick interval to the nearest 1, 2 or 5 times a power of ten.
///
/// The thresholds sit at 1.5, 3 and 7 rather than at 1, 2 and 5, so an interval
/// of 2.06 becomes 2 and not 5. Rounding up at the boundary halves the number
/// of ticks for a rounding error's worth of input.
pub(crate) fn nice_step(span: f64, target_ticks: f64) -> u64 {
    if !span.is_finite() || span <= 0.0 || target_ticks <= 0.0 {
        return 1;
    }
    let raw = span / target_ticks;
    if raw <= 1.0 {
        return 1;
    }
    let magnitude = 10f64.powf(raw.log10().floor());
    let normalised = raw / magnitude;
    let multiplier = if normalised <= 1.5 {
        1.0
    } else if normalised <= 3.0 {
        2.0
    } else if normalised <= 7.0 {
        5.0
    } else {
        10.0
    };
    ((multiplier * magnitude).round() as u64).max(1)
}

/// The unit every label on one ruler is printed in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TickUnit {
    divisor: u64,
    suffix: &'static str,
    decimals: usize,
}

impl TickUnit {
    fn format(&self, pos: u64) -> String {
        if self.divisor <= 1 {
            group_thousands(pos)
        } else {
            format!(
                "{:.*}{}",
                self.decimals,
                pos as f64 / self.divisor as f64,
                self.suffix
            )
        }
    }
}

/// Picks bp, kb or Mb from how far along the sequence the ruler sits and how
/// fine its step is.
///
/// Both matter. Position alone would print `5.000 kb` for a 20 bp window at
/// position 5000, and step alone would print `1500 kb` where `1.5 Mb` belongs.
fn tick_unit(max_pos: u64, step: u64) -> TickUnit {
    // Three decimals of a megabase is kilobase resolution, still readable.
    // Three decimals of a kilobase is just a base with a decimal point in it,
    // so kilobases are held to two.
    if max_pos >= 1_000_000 {
        if let Some(decimals) = decimals_for(step, 1_000_000, 3) {
            return TickUnit {
                divisor: 1_000_000,
                suffix: " Mb",
                decimals,
            };
        }
    }
    if max_pos >= 10_000 {
        if let Some(decimals) = decimals_for(step, 1_000, 2) {
            return TickUnit {
                divisor: 1_000,
                suffix: " kb",
                decimals,
            };
        }
    }
    TickUnit {
        divisor: 1,
        suffix: "",
        decimals: 0,
    }
}

/// Decimal places needed to write `step` in units of `unit`.
///
/// `None` means the unit is too coarse for this step even at `max_decimals`, so
/// the caller should try a finer one.
fn decimals_for(step: u64, unit: u64, max_decimals: usize) -> Option<usize> {
    (0..=max_decimals).find(|decimals| {
        10u64
            .checked_pow(*decimals as u32)
            .and_then(|factor| step.checked_mul(factor))
            .is_some_and(|scaled| scaled % unit == 0)
    })
}

/// `1234567` as `1,234,567`.
///
/// Shared with the tracks rather than copied into them: a tooltip quotes the
/// coordinate a reader would type, and a coordinate that is grouped on the
/// ruler and ungrouped in a tooltip is two conventions on one figure.
///
/// # The rule the whole crate follows
///
/// **Every whole number a tooltip writes goes through here**, coordinates and
/// counts alike, whatever its magnitude. `26 of 26 sites called` and
/// `12,480 of 12,500 sites called` are the same sentence, and a rule that
/// grouped only the numbers a reader happened to find long would make the
/// tooltip's punctuation a fact about the data rather than a convention.
///
/// There is one exception, and it is written down where it is taken: a
/// [`SnpTrack`](crate::SnpTrack) column prints its own position under itself
/// in whatever the caller was counting in, and the tooltip echoes that label
/// character for character rather than disagreeing with the figure it sits on.
pub(crate) fn group_thousands(value: u64) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steps_are_round_numbers() {
        assert_eq!(nice_step(1000.0, 5.0), 200);
        assert_eq!(nice_step(10_000.0, 5.0), 2000);
        assert_eq!(nice_step(30.0, 6.0), 5);
        assert_eq!(nice_step(80.0, 10.0), 10);
    }

    #[test]
    fn a_step_just_over_a_boundary_is_not_rounded_all_the_way_up() {
        // 1_234_567 / 6 is 205_761, which is 2.06 hundred-thousands. Rounding
        // up to 5 would leave the ruler with two ticks.
        assert_eq!(nice_step(1_234_567.0, 6.0), 200_000);
    }

    #[test]
    fn steps_never_drop_below_one_base() {
        assert_eq!(nice_step(5.0, 10.0), 1);
        assert_eq!(nice_step(0.0, 10.0), 1);
        assert_eq!(nice_step(f64::NAN, 10.0), 1);
    }

    #[test]
    fn ticks_land_on_multiples_of_the_step() {
        let axis = AxisTrack::new();
        let ticks = axis.ticks(1, 10_000, 800.0);
        assert!(ticks.positions.len() >= 2);
        for pos in &ticks.positions {
            assert_eq!(pos % ticks.step, 0, "{pos} is not a multiple of the step");
            assert!((1..=10_000).contains(pos));
        }
    }

    #[test]
    fn a_three_base_region_gets_one_tick_per_base() {
        let axis = AxisTrack::new();
        let ticks = axis.ticks(1_000_001, 1_000_003, 800.0);
        assert_eq!(ticks.step, 1);
        assert_eq!(ticks.positions, vec![1_000_001, 1_000_002, 1_000_003]);
    }

    #[test]
    fn a_ruler_running_to_the_top_of_the_coordinate_range_stops_at_its_last_tick() {
        // Region::new("chr1", 0, u64::MAX) is a legal region, so the ruler over
        // it has to end rather than wrap: the step used to overflow, which
        // panicked in a debug build and drew 503 ticks counting backwards in a
        // release one.
        let axis = AxisTrack::new();
        let ticks = axis.ticks(1, u64::MAX, 872.0);
        assert_eq!(ticks.step, 2_000_000_000_000_000_000);
        assert_eq!(ticks.positions.len(), 9);
        assert_eq!(ticks.positions[0], 2_000_000_000_000_000_000);
        assert_eq!(
            ticks.positions[8], 18_000_000_000_000_000_000,
            "the last tick is the last multiple of the step that fits"
        );
        assert!(ticks.positions.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn a_single_base_region_falls_back_to_labelling_its_own_edge() {
        let axis = AxisTrack::new();
        let ticks = axis.ticks(42, 42, 800.0);
        assert_eq!(ticks.positions, vec![42]);
    }

    #[test]
    fn the_unit_follows_both_the_position_and_the_step() {
        assert_eq!(tick_unit(1_500_000, 500_000).format(1_500_000), "1.5 Mb");
        assert_eq!(tick_unit(4_000_000, 1_000_000).format(2_000_000), "2 Mb");
        assert_eq!(tick_unit(1_234_000, 1_000).format(1_234_000), "1.234 Mb");
        assert_eq!(tick_unit(20_000, 500).format(2_500), "2.5 kb");
        assert_eq!(tick_unit(900, 100).format(750), "750");
    }

    #[test]
    fn a_zoomed_in_window_far_along_a_chromosome_keeps_base_labels() {
        // A 20 bp window at position 5,000,000: the position says Mb, but a
        // 5 bp step cannot be written in Mb, so the ruler stays in bases.
        let unit = tick_unit(5_000_020, 5);
        assert_eq!(unit.divisor, 1);
        assert_eq!(unit.format(5_000_010), "5,000,010");
    }

    #[test]
    fn a_unit_too_coarse_for_the_step_is_rejected() {
        assert_eq!(decimals_for(500_000, 1_000_000, 3), Some(1));
        assert_eq!(decimals_for(1_000, 1_000_000, 3), Some(3));
        assert_eq!(decimals_for(100, 1_000_000, 3), None);
        assert_eq!(decimals_for(5, 1_000, 2), None);
        assert_eq!(decimals_for(1, 1, 0), Some(0));
    }

    #[test]
    fn centring_shifts_every_tick_by_half_a_base() {
        use crate::figure::Figure;
        use crate::region::Region;

        let region = Region::new("motif", 0, 10).unwrap();
        let edges = Figure::new(region.clone())
            .show_region_label(false)
            .push(AxisTrack::new())
            .to_svg();
        let centres = Figure::new(region)
            .show_region_label(false)
            .push(AxisTrack::new().center_on_bases(true))
            .to_svg();

        let first_tick = |svg: &str| -> f64 {
            let start = svg.find("<line x1=\"").unwrap() + "<line x1=\"".len();
            let rest = &svg[start..];
            // Skip the baseline, which spans the plot, and take the next line.
            let next = rest.find("<line x1=\"").unwrap() + "<line x1=\"".len();
            let rest = &rest[next..];
            rest[..rest.find('"').unwrap()].parse().unwrap()
        };

        let half_base = (900.0 - 16.0 - 18.0) / 10.0 / 2.0;
        assert!((first_tick(&centres) - first_tick(&edges) - half_base).abs() < 1e-6);
    }

    #[test]
    fn thousands_are_grouped() {
        assert_eq!(group_thousands(0), "0");
        assert_eq!(group_thousands(999), "999");
        assert_eq!(group_thousands(1_000), "1,000");
        assert_eq!(group_thousands(1_234_567), "1,234,567");
    }
}
