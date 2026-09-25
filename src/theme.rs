//! Colours, fonts and spacing shared by every track in a figure.
//!
//! A figure builds one [`Theme`] and hands every track a `&Theme` when its turn
//! to draw comes, so appearance is settled once for the whole picture rather
//! than argued track by track. Nothing here draws: it is the values a drawing
//! reads from, and the few functions for deriving a value the theme does not
//! already carry.
//!
//! # A track never writes a colour down
//!
//! A shade the theme has no field for is mixed from one it does, which is what
//! [`mix`], [`wash`] and [`contrast_ink`] are for and why they sit beside the
//! theme rather than in the tracks that call them. A hardcoded grey is the same
//! grey on a white page and on a black one, and it cannot be right on both.
//!
//! # Nothing turns a light figure dark
//!
//! [`Theme::dark`] is a second set of measured values, not [`Theme::light`] put
//! through a transform, and there is no call that converts one into the other.
//! The page a figure is going on is therefore a decision taken before it is
//! built: a track takes its ink from whatever it is handed and never asks which
//! way the page runs.
//!
//! # Where convention and legibility disagree
//!
//! The categorical palette is measured for colour vision deficiency rather than
//! chosen by eye, and [`Theme::color`] hands its entries out in a fixed order
//! with [`Theme::accent`] first, so a figure with one series and a figure with
//! six agree about what the first colour is. The two themes hand out the same
//! six hues in the same order, each stepped for its own page, so a figure keeps
//! its colours' names when it changes page. The nucleotide colours are the one
//! place where the measurement does not decide: convention is the default, and
//! [`BaseColors::colorblind_safe`] ships beside it for the figures where that
//! default costs too much.
//!
//! # Type
//!
//! Text is set in Inter, with the numbers a reader reads off a coordinate, the
//! ruler and the locus, in JetBrains Mono. Neither is required: both stacks
//! fall back to fonts that are on every system, and text is measured wide
//! enough for whichever of them draws it (see [`text_width`](crate::svg::text_width)).

use crate::style::{Emphasis, LinePattern, MarkStyle, RenderProfile, Symbol, VisualTokens};

/// Visual settings applied to a whole [`Figure`](crate::Figure).
///
/// Every field is public: start from [`Theme::light`] or [`Theme::dark`] and
/// overwrite what you need, rather than building one from scratch.
#[derive(Debug, Clone, PartialEq)]
pub struct Theme {
    /// Page background. Set it to `"none"` for a transparent SVG.
    pub background: String,
    /// Titles, branches and text drawn on marks: the strongest ink.
    pub foreground: String,
    /// Secondary text: track names, coordinates, tick labels, legends.
    pub muted: String,
    /// Quiet rules: baselines, the ruler and guides.
    pub rule: String,
    /// Default colour for a track that was not given one.
    pub accent: String,
    /// Categorical colours, cycled by [`Theme::color`].
    pub palette: Vec<String>,
    /// Per-nucleotide colours used by the sequence track.
    pub bases: BaseColors,
    /// Colour of an insertion mark in a read pileup, kept off the categorical
    /// palette so it never impersonates a series.
    pub insertion: String,
    /// Corner radius of a data mark, in pixels.
    ///
    /// Rounded ends are most of the difference between a figure that looks
    /// drawn and one that looks emitted. Set it to zero for square corners.
    pub corner_radius: f64,
    /// Font stack written on the root `<svg>` element.
    pub font_family: String,
    /// Font stack for coordinates: the ruler, the locus beside the title.
    ///
    /// Digits of one width line a column of positions up, and a position is
    /// read digit by digit against the one beside it.
    pub mono_family: String,
    /// Size of tick labels and in-plot annotations, in pixels.
    pub font_size: f64,
    /// Size of track labels in the left gutter, in pixels.
    pub label_font_size: f64,
    /// Size of the figure title, in pixels.
    pub title_font_size: f64,
    /// Height of a capital letter as a fraction of the font size.
    ///
    /// A sequence logo stretches each letter to an exact box, which means
    /// working back from the box height to a font size. The default suits the
    /// Arial-compatible metrics of [`Theme::font_family`]; change both together
    /// or logo letters will sit slightly proud of their boxes.
    pub cap_height_ratio: f64,
    /// Shared measurements for marks, guides, legends and annotations.
    pub tokens: VisualTokens,
}

impl Theme {
    /// Dark ink on a white page, for manuscripts and slides.
    ///
    /// The ink is a deep indigo rather than a neutral black, the violet of the
    /// dark theme's page taken down to where it reads as black on paper.
    pub fn light() -> Self {
        Theme {
            background: "#ffffff".into(),
            foreground: "#1a1233".into(),
            muted: "#5b5480".into(),
            rule: "#e4def5".into(),
            accent: "#1634c2".into(),
            palette: noche(),
            bases: BaseColors::default(),
            insertion: "#8e44ad".into(),
            corner_radius: 3.5,
            font_family: TEXT_FONTS.into(),
            mono_family: MONO_FONTS.into(),
            font_size: 11.5,
            label_font_size: 11.0,
            title_font_size: 17.0,
            cap_height_ratio: 0.72,
            tokens: VisualTokens::default(),
        }
    }

    /// Light ink on a deep indigo page, for dark-mode documents and screens.
    ///
    /// The page is the one the documentation site draws its figures on at
    /// night, so a figure made with this sits on it without a seam.
    pub fn dark() -> Self {
        Theme {
            background: "#120b2b".into(),
            foreground: "#efeaff".into(),
            muted: "#b3a9e0".into(),
            rule: "#2e2558".into(),
            accent: "#6275fc".into(),
            palette: noche_dark(),
            bases: BaseColors::default(),
            insertion: "#c59fe1".into(),
            corner_radius: 3.5,
            font_family: TEXT_FONTS.into(),
            mono_family: MONO_FONTS.into(),
            font_size: 11.5,
            label_font_size: 11.0,
            title_font_size: 17.0,
            cap_height_ratio: 0.72,
            tokens: VisualTokens::default(),
        }
    }

    /// Colour `index` of the categorical palette, wrapping around at the end.
    pub fn color(&self, index: usize) -> &str {
        if self.palette.is_empty() {
            return &self.accent;
        }
        &self.palette[index % self.palette.len()]
    }

    /// Starts a theme from one of the named output profiles.
    pub fn for_profile(profile: RenderProfile) -> Self {
        let theme = if profile.is_dark() {
            Theme::dark()
        } else {
            Theme::light()
        };
        theme.scaled(profile.visual_scale())
    }

    /// Resolves semantic prominence into measurable, colour-independent marks.
    pub fn mark_style(&self, emphasis: Emphasis) -> MarkStyle {
        match emphasis {
            Emphasis::Muted => MarkStyle {
                stroke_width: self.tokens.hairline,
                marker_radius: self.tokens.marker_radius * 0.8,
                opacity: 0.55,
                pattern: LinePattern::Dotted,
            },
            Emphasis::Normal => MarkStyle {
                stroke_width: self.tokens.stroke,
                marker_radius: self.tokens.marker_radius,
                opacity: 0.9,
                pattern: LinePattern::Solid,
            },
            Emphasis::Primary => MarkStyle {
                stroke_width: self.tokens.strong_stroke,
                marker_radius: self.tokens.marker_radius * 1.15,
                opacity: 1.0,
                pattern: LinePattern::Solid,
            },
            Emphasis::Alert => MarkStyle {
                stroke_width: self.tokens.strong_stroke * 1.15,
                marker_radius: self.tokens.marker_radius * 1.25,
                opacity: 1.0,
                pattern: LinePattern::Dashed,
            },
        }
    }

    /// Returns a categorical point shape, wrapping in the same way as colours.
    pub fn symbol(&self, index: usize) -> Symbol {
        const SYMBOLS: [Symbol; 4] = [
            Symbol::Circle,
            Symbol::Square,
            Symbol::Diamond,
            Symbol::Triangle,
        ];
        SYMBOLS[index % SYMBOLS.len()]
    }

    /// Returns the same visual system with all typographic chrome scaled by
    /// `factor`.
    ///
    /// Data coordinates and the canvas size do not change. This is intended
    /// for a larger presentation or a compact multi-panel figure without
    /// making callers adjust three font sizes and the corner radius by hand.
    /// Non-finite factors are ignored and the lower bound prevents invisible
    /// text and degenerate corners.
    pub fn scaled(mut self, factor: f64) -> Self {
        let factor = if factor.is_finite() {
            factor.max(0.25)
        } else {
            1.0
        };
        self.corner_radius *= factor;
        self.font_size *= factor;
        self.label_font_size *= factor;
        self.title_font_size *= factor;
        self.tokens = self.tokens.scaled(factor);
        self
    }

    /// A real colour standing in for the page, to blend a tint against.
    ///
    /// [`Theme::background`] is allowed to be `"none"`, which is a transparent
    /// page rather than a colour, and [`mix`] blending against something that
    /// is not a colour gives back something that is not a colour: the mark
    /// comes out `fill="none"` and disappears. A quiet cell in a matrix, the
    /// bar under a base that agrees with the reference and the zebra tint on a
    /// panel of variable sites are all tints towards the page, so they ask for
    /// this instead.
    ///
    /// A transparent page is whatever it ends up composited onto, which cannot
    /// be known from here. The theme's own ink is the one clue it carries, so
    /// dark ink implies a light page and light ink a dark one.
    pub fn surface(&self) -> &str {
        if parse_hex(&self.background).is_some() {
            &self.background
        } else {
            contrast_ink(&self.foreground)
        }
    }
}

impl Default for Theme {
    fn default() -> Self {
        Theme::light()
    }
}

/// Text stack: Inter, then fonts every system has.
const TEXT_FONTS: &str = "Inter, Liberation Sans, Arial, Helvetica, sans-serif";

/// Coordinate stack: JetBrains Mono, then monospaced fonts every system has.
const MONO_FONTS: &str = "JetBrains Mono, Liberation Mono, Menlo, Consolas, monospace";

/// Six hues for a white page: indigo, pink, teal, ochre, plum and orange.
///
/// Chosen by slot, one hue range each so the palette reads as six named
/// colours, and inside each range the lightness and chroma that keep the worst
/// pair furthest apart. Measured in OKLab hundredths after simulating
/// protanopia, deuteranopia and tritanopia (Machado, Oliveira and Fernandes,
/// 2009, at full severity), the closest pair under any of them is 10.7, the
/// ochre against the orange under deuteranopia. The Okabe-Ito set this
/// replaces came to 6.7. Every entry reaches a contrast of 3.1 against white,
/// which is the floor for a mark that has to be seen rather than read.
fn noche() -> Vec<String> {
    [
        "#1634c2", "#e63f9f", "#108169", "#b78a2c", "#69437c", "#b0540e",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

/// The same six hues stepped for the dark page, in the same order.
///
/// A dark page wants its marks lighter and a narrower band of lightness, so
/// these are measured against it rather than lightened from the light set. The
/// closest pair is 14.4, the yellow against the orange under deuteranopia, and
/// the lowest contrast against the page is 4.9. The set this replaces measured
/// 1.6 under deuteranopia, its green against its pink, which a reader with
/// that deficiency cannot tell apart.
fn noche_dark() -> Vec<String> {
    [
        "#6275fc", "#e4488c", "#41f5ea", "#fcdc67", "#b7a1f5", "#ee9a69",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

/// Colours for the four nucleotides plus a fallback for everything else.
///
/// Two sets ship: [`BaseColors::conventional`], which is the default and is
/// what a reader expects, and [`BaseColors::colorblind_safe`], which is what
/// survives being read by everyone. The difference is measured, not asserted;
/// both constructors carry the numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseColors {
    /// Colour of adenine.
    pub a: String,
    /// Colour of cytosine.
    pub c: String,
    /// Colour of guanine.
    pub g: String,
    /// Colour of thymine, also used for uracil.
    pub t: String,
    /// Colour of N and any other IUPAC or unexpected symbol.
    pub other: String,
}

impl BaseColors {
    /// Colour for one base, case insensitive, with `U` treated as `T`.
    pub fn of(&self, base: u8) -> &str {
        match base.to_ascii_uppercase() {
            b'A' => &self.a,
            b'C' => &self.c,
            b'G' => &self.g,
            b'T' | b'U' => &self.t,
            _ => &self.other,
        }
    }
}

impl BaseColors {
    /// The conventional colours, the ones a reader of genome figures already
    /// has in their head.
    ///
    /// They are not safe for every reader. Measured pairwise, adenine and
    /// guanine come out 1.7 apart under protanopia, on a scale where 8 is the
    /// floor for telling two colours apart. That is the transition pair, the
    /// commonest substitution there is, so the convention costs a real reader a
    /// real distinction. Convention is still the default because a figure that
    /// recolours the bases surprises everyone, but see
    /// [`BaseColors::colorblind_safe`] before sending anything to print.
    pub fn conventional() -> Self {
        BaseColors {
            a: "#33a02c".into(),
            c: "#1f78b4".into(),
            g: "#e08214".into(),
            t: "#e31a1c".into(),
            other: "#9e9e9e".into(),
        }
    }

    /// Four bases that stay apart for every reader.
    ///
    /// The closest pair is 11.0 in OKLab hundredths under deuteranopia, against
    /// a floor of 8, and every other pair is further apart than that. The cost
    /// is in the shades, not in which base is which: adenine is green in both
    /// sets, a bluish green here, guanine turns a yellower orange and cytosine
    /// barely moves. Thymine moves furthest, from red to vermillion, so a
    /// reader expecting a red thymine finds an orange-red one: label it or say
    /// so in the caption.
    pub fn colorblind_safe() -> Self {
        BaseColors {
            a: "#009e73".into(),
            c: "#0072b2".into(),
            g: "#e69f00".into(),
            t: "#d55e00".into(),
            other: "#9e9e9e".into(),
        }
    }
}

impl Default for BaseColors {
    fn default() -> Self {
        BaseColors::conventional()
    }
}

/// Blends two `#rrggbb` colours, `t` running from all of `a` to all of `b`.
///
/// What it is for: a greyscale that belongs to its theme. Cytogenetic bands are
/// specified as shades from white to black, and hardcoding those shades makes a
/// dark figure look like a light one someone forgot to invert. Mixing the
/// theme's own background and foreground gives the same ladder in whichever
/// direction the page runs.
pub fn mix(a: &str, b: &str, t: f64) -> String {
    let (Some(from), Some(to)) = (parse_hex(a), parse_hex(b)) else {
        return a.to_string();
    };
    // `clamp` propagates a NaN, and `NaN as u8` saturates to nought, so an
    // amount that is not a number used to leave the ramp entirely and come
    // back pure black: darker than `b`, on a scale where darker means more.
    // Whatever the caller meant by it, the answer has to be a colour on this
    // ramp, and of the two ends the pale one is the one that does not read as
    // the strongest evidence on the page.
    let t = if t.is_nan() { 0.0 } else { t.clamp(0.0, 1.0) };
    let channel = |from: u8, to: u8| (from as f64 + (to as f64 - from as f64) * t).round() as u8;
    format!(
        "#{:02x}{:02x}{:02x}",
        channel(from.0, to.0),
        channel(from.1, to.1),
        channel(from.2, to.2)
    )
}

/// Reads a `#rrggbb` string, and nothing else.
fn parse_hex(color: &str) -> Option<(u8, u8, u8)> {
    let hex = color.strip_prefix('#')?;
    if hex.len() != 6 {
        return None;
    }
    Some((
        u8::from_str_radix(&hex[0..2], 16).ok()?,
        u8::from_str_radix(&hex[2..4], 16).ok()?,
        u8::from_str_radix(&hex[4..6], 16).ok()?,
    ))
}

/// A pale wash of a colour against the page.
///
/// What a large filled shape wears. A gene arrow or a matrix cell at full
/// saturation shouts, and eight of them make a colour chart rather than a
/// figure; the hue belongs in the edge and a wash of it in the body. Small
/// marks are the other way round and want the colour undiluted, which is why
/// this is a function rather than a rule.
pub fn wash(color: &str, theme: &Theme) -> String {
    mix(color, theme.surface(), WASH)
}

/// How far a wash is pulled towards the page.
const WASH: f64 = 0.62;

/// Black or white, whichever stays readable on top of `color`.
///
/// Feature labels are drawn inside coloured boxes whose colour the caller
/// chooses, so the ink has to be picked per box rather than fixed by the theme.
/// Anything that is not a `#rrggbb` string is assumed dark.
///
/// The choice is the WCAG contrast ratio of each ink against `color`, computed
/// from sRGB relative luminance. Rec. 601 luma over the gamma-encoded channels
/// is a different quantity and is not monotone in relative luminance across
/// hues, so on a mid saturated orange or green it picks the less legible of the
/// two by a wide margin.
pub fn contrast_ink(color: &str) -> &'static str {
    let Some((r, g, b)) = parse_hex(color) else {
        return "#ffffff";
    };
    let luma = relative_luminance(r, g, b);
    let dark = relative_luminance(DARK_INK_RGB.0, DARK_INK_RGB.1, DARK_INK_RGB.2);
    // Contrast ratio is (lighter + 0.05) / (darker + 0.05), and `color` always
    // sits between the two inks, so each ratio is one division. White has a
    // relative luminance of exactly 1.
    let against_white = 1.05 / (luma + 0.05);
    let against_dark = (luma + 0.05) / (dark + 0.05);
    if against_dark > against_white {
        DARK_INK
    } else {
        "#ffffff"
    }
}

/// The dark ink [`contrast_ink`] hands out, which is also
/// [`Theme::light`]'s foreground, and its channels.
const DARK_INK: &str = "#1a1233";
const DARK_INK_RGB: (u8, u8, u8) = (0x1a, 0x12, 0x33);

/// sRGB relative luminance, as WCAG 2.x defines it: undo the transfer function
/// on each channel, then weight them by the response of the eye.
fn relative_luminance(r: u8, g: u8, b: u8) -> f64 {
    let linear = |c: u8| {
        let s = c as f64 / 255.0;
        if s <= 0.04045 {
            s / 12.92
        } else {
            ((s + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * linear(r) + 0.7152 * linear(g) + 0.0722 * linear(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every ramp in the crate goes through here, so an amount that is not a
    /// number used to be one line away from a black mark on any of them. Black
    /// is not on this ramp: it is past the dark end, which on a scale where
    /// darker means more reads as beyond the strongest value anyone measured.
    #[test]
    fn an_amount_that_is_not_a_number_stays_on_the_ramp() {
        let pale = "#ffffff";
        let dark = "#1b1f23";
        assert_eq!(mix(pale, dark, f64::NAN), pale, "NaN left the ramp");
        // The infinities already clamped, and still do.
        assert_eq!(mix(pale, dark, f64::INFINITY), dark);
        assert_eq!(mix(pale, dark, f64::NEG_INFINITY), pale);
        // And the ordinary range is untouched by the guard.
        assert_eq!(mix(pale, dark, 0.0), pale);
        assert_eq!(mix(pale, dark, 1.0), dark);
        assert_eq!(mix(pale, dark, 0.5), "#8d8f91");
    }

    #[test]
    fn a_transparent_page_still_offers_a_colour_to_blend_against() {
        // Otherwise every tint towards the page comes back as `none`, and the
        // marks that carry "present but quiet" vanish from the figure.
        let mut light = Theme::light();
        assert_eq!(light.surface(), "#ffffff");
        light.background = "none".into();
        assert_eq!(light.surface(), "#ffffff", "dark ink means a light page");
        assert_ne!(mix(light.surface(), &light.rule, 0.2), "none");

        let mut dark = Theme::dark();
        assert_eq!(dark.surface(), Theme::dark().background);
        dark.background = "none".into();
        assert_eq!(dark.surface(), DARK_INK, "light ink means a dark page");
    }

    #[test]
    fn ink_flips_with_the_brightness_of_the_box() {
        assert_eq!(contrast_ink("#000000"), "#ffffff");
        assert_eq!(contrast_ink("#0072b2"), "#ffffff");
        assert_eq!(contrast_ink("#ffffff"), DARK_INK);
        assert_eq!(contrast_ink("#f0e442"), DARK_INK);
    }

    #[test]
    fn ink_follows_the_contrast_ratio_and_not_gamma_encoded_luma() {
        // Rec. 601 luma over the encoded channels put white on all four of
        // these. The WCAG ratios say otherwise for the green and the orange:
        // #33a02c is 3.384 against white and 4.898 against the dark ink, and
        // #e08214 is 2.848 against white, under the 3:1 floor, against 5.820.
        let bases = BaseColors::conventional();
        assert_eq!(contrast_ink(&bases.a), DARK_INK, "A #33a02c");
        assert_eq!(contrast_ink(&bases.c), "#ffffff", "C #1f78b4");
        assert_eq!(contrast_ink(&bases.g), DARK_INK, "G #e08214");
        assert_eq!(contrast_ink(&bases.t), "#ffffff", "T #e31a1c");

        let safe = BaseColors::colorblind_safe();
        assert_eq!(contrast_ink(&safe.a), DARK_INK, "A #009e73");
        assert_eq!(contrast_ink(&safe.c), "#ffffff", "C #0072b2");
        assert_eq!(contrast_ink(&safe.g), DARK_INK, "G #e69f00");
        assert_eq!(contrast_ink(&safe.t), DARK_INK, "T #d55e00");
    }

    #[test]
    fn the_ink_it_picks_is_the_one_with_the_higher_contrast_ratio() {
        // The property, checked over the whole cube rather than a palette.
        let ratio = |x: f64, y: f64| {
            let (hi, lo) = if x > y { (x, y) } else { (y, x) };
            (hi + 0.05) / (lo + 0.05)
        };
        let dark = relative_luminance(DARK_INK_RGB.0, DARK_INK_RGB.1, DARK_INK_RGB.2);
        for r in (0..=255).step_by(17) {
            for g in (0..=255).step_by(17) {
                for b in (0..=255).step_by(17) {
                    let (r, g, b) = (r as u8, g as u8, b as u8);
                    let hex = format!("#{r:02x}{g:02x}{b:02x}");
                    let luma = relative_luminance(r, g, b);
                    let want = if ratio(luma, dark) > ratio(luma, 1.0) {
                        DARK_INK
                    } else {
                        "#ffffff"
                    };
                    assert_eq!(contrast_ink(&hex), want, "{hex}");
                }
            }
        }
    }

    #[test]
    fn ink_falls_back_to_white_for_colours_it_cannot_read() {
        assert_eq!(contrast_ink("red"), "#ffffff");
        assert_eq!(contrast_ink("#fff"), "#ffffff");
        assert_eq!(contrast_ink("#zzzzzz"), "#ffffff");
    }

    /// sRGB channel to linear light.
    fn linear(c: u8) -> f64 {
        let s = c as f64 / 255.0;
        if s <= 0.04045 {
            s / 12.92
        } else {
            ((s + 0.055) / 1.055).powf(2.4)
        }
    }

    /// Linear RGB to OKLab, in hundredths, which is the unit the palettes'
    /// comments quote.
    fn oklab([r, g, b]: [f64; 3]) -> [f64; 3] {
        let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
        let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
        let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
        [
            100.0 * (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s),
            100.0 * (1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s),
            100.0 * (0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s),
        ]
    }

    /// The worst pair of a palette in OKLab hundredths, over normal vision and
    /// protanopia, deuteranopia and tritanopia simulated at full severity
    /// (Machado, Oliveira and Fernandes, 2009).
    fn worst_pair(palette: &[String]) -> f64 {
        const VISIONS: [[[f64; 3]; 3]; 4] = [
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [
                [0.152286, 1.052583, -0.204868],
                [0.114503, 0.786281, 0.099216],
                [-0.003882, -0.048116, 1.051998],
            ],
            [
                [0.367322, 0.860646, -0.227968],
                [0.280085, 0.672501, 0.047413],
                [-0.011820, 0.042940, 0.968881],
            ],
            [
                [1.255528, -0.076749, -0.178779],
                [-0.078411, 0.930809, 0.147602],
                [0.004733, 0.691367, 0.303900],
            ],
        ];
        let seen = |colour: &str, vision: &[[f64; 3]; 3]| {
            let (r, g, b) = parse_hex(colour).unwrap();
            let rgb = [linear(r), linear(g), linear(b)];
            oklab(std::array::from_fn(|row| {
                (0..3)
                    .map(|col| vision[row][col] * rgb[col])
                    .sum::<f64>()
                    .clamp(0.0, 1.0)
            }))
        };
        let mut worst = f64::INFINITY;
        for (i, a) in palette.iter().enumerate() {
            for b in &palette[i + 1..] {
                for vision in &VISIONS {
                    let (p, q) = (seen(a, vision), seen(b, vision));
                    let distance =
                        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2))
                            .sqrt();
                    worst = worst.min(distance);
                }
            }
        }
        worst
    }

    #[test]
    fn every_pair_of_either_palette_stays_apart_for_every_reader() {
        // The numbers the palettes' comments quote, held to. The dark palette
        // these replaced said it had been validated and came to 1.6.
        let light = worst_pair(&Theme::light().palette);
        let dark = worst_pair(&Theme::dark().palette);
        assert!((light - 10.7).abs() < 0.05, "light {light:.2}");
        assert!((dark - 14.4).abs() < 0.05, "dark {dark:.2}");
        // And every mark can be seen on its own page, 3:1 being the floor for
        // something that has to be seen rather than read.
        for theme in [Theme::light(), Theme::dark()] {
            let (r, g, b) = parse_hex(&theme.background).unwrap();
            let page = relative_luminance(r, g, b);
            for colour in &theme.palette {
                let (r, g, b) = parse_hex(colour).unwrap();
                let mark = relative_luminance(r, g, b);
                let (hi, lo) = if mark > page {
                    (mark, page)
                } else {
                    (page, mark)
                };
                assert!(
                    (hi + 0.05) / (lo + 0.05) >= 3.0,
                    "{colour} on {}",
                    theme.background
                );
            }
        }
    }

    #[test]
    fn the_two_palettes_name_the_same_hues_in_the_same_order() {
        // A figure keeps the names of its colours when it changes page: the
        // first series is the indigo one on either. Checked on the hue angle
        // in OKLab, which is what "the same colour" means across a lightness
        // step.
        let hue = |colour: &str| {
            let (r, g, b) = parse_hex(colour).unwrap();
            let [_, a, b] = oklab([linear(r), linear(g), linear(b)]);
            b.atan2(a).to_degrees().rem_euclid(360.0)
        };
        for (light, dark) in Theme::light().palette.iter().zip(&Theme::dark().palette) {
            let apart = (hue(light) - hue(dark)).abs();
            let apart = apart.min(360.0 - apart);
            assert!(
                apart < 30.0,
                "{light} and {dark} are {apart:.0} degrees apart"
            );
        }
    }

    #[test]
    fn palette_wraps_instead_of_panicking() {
        let theme = Theme::light();
        let n = theme.palette.len();
        assert_eq!(theme.color(0), theme.color(n));
        assert_eq!(theme.color(n + 3), theme.color(3));
    }

    #[test]
    fn palette_falls_back_to_accent_when_emptied() {
        let mut theme = Theme::light();
        theme.palette.clear();
        assert_eq!(theme.color(2), theme.accent);
    }

    #[test]
    fn bases_are_case_insensitive_and_map_u_to_t() {
        let bases = BaseColors::default();
        assert_eq!(bases.of(b'a'), bases.of(b'A'));
        assert_eq!(bases.of(b'u'), bases.of(b'T'));
        assert_eq!(bases.of(b'N'), bases.other);
        assert_eq!(bases.of(b'-'), bases.other);
    }

    #[test]
    fn the_accent_is_the_first_slot_of_its_own_palette() {
        // The accent is what a single series gets. If it were not slot one, a
        // one series figure and a two series figure would disagree about what
        // the first colour is.
        for theme in [Theme::light(), Theme::dark()] {
            assert_eq!(theme.accent, theme.palette[0]);
        }
    }

    #[test]
    fn the_two_palettes_are_different_sets_not_a_flip() {
        // A dark background wants a narrower lightness band, so the dark steps
        // were chosen against it rather than derived from the light ones.
        assert_ne!(Theme::light().palette, Theme::dark().palette);
        assert_eq!(Theme::light().palette.len(), Theme::dark().palette.len());
    }

    #[test]
    fn the_colourblind_safe_bases_differ_from_the_conventional_ones() {
        let conventional = BaseColors::conventional();
        let safe = BaseColors::colorblind_safe();
        assert_ne!(conventional.a, safe.a);
        assert_eq!(BaseColors::default(), conventional);
        // Both still answer for every symbol.
        for base in *b"ACGTUN" {
            assert!(safe.of(base).starts_with('#'));
        }
    }

    #[test]
    fn every_palette_entry_is_a_hex_colour() {
        for theme in [Theme::light(), Theme::dark()] {
            for color in &theme.palette {
                assert_eq!(color.len(), 7, "{color} is not #rrggbb");
                assert!(color.starts_with('#'), "{color} is not #rrggbb");
                assert!(
                    color[1..].chars().all(|c| c.is_ascii_hexdigit()),
                    "{color} is not #rrggbb"
                );
            }
        }
    }

    #[test]
    fn scaling_keeps_the_palette_and_moves_the_visual_sizes_together() {
        let original = Theme::light();
        let scaled = original.clone().scaled(1.5);

        assert_eq!(scaled.palette, original.palette);
        assert_eq!(scaled.font_size, original.font_size * 1.5);
        assert_eq!(scaled.label_font_size, original.label_font_size * 1.5);
        assert_eq!(scaled.title_font_size, original.title_font_size * 1.5);
        assert_eq!(scaled.corner_radius, original.corner_radius * 1.5);
        assert_eq!(scaled.tokens.stroke, original.tokens.stroke * 1.5);
    }

    #[test]
    fn profiles_and_emphasis_are_resolved_in_one_place() {
        assert_eq!(
            Theme::for_profile(RenderProfile::Dark).background,
            Theme::dark().background
        );
        assert!(
            Theme::for_profile(RenderProfile::Presentation).font_size > Theme::light().font_size
        );
        assert!(
            Theme::light().mark_style(Emphasis::Primary).stroke_width
                > Theme::light().mark_style(Emphasis::Muted).stroke_width
        );
        assert_ne!(Theme::light().symbol(0), Theme::light().symbol(1));
        assert_eq!(Theme::light().symbol(0), Theme::light().symbol(4));
    }

    #[test]
    fn an_invalid_visual_scale_falls_back_to_the_original_sizes() {
        let original = Theme::light();
        assert_eq!(original.clone().scaled(f64::NAN), original);
        assert_eq!(original.clone().scaled(f64::INFINITY), original);
        assert_eq!(
            original.clone().scaled(0.0).font_size,
            original.font_size * 0.25
        );
    }
}
