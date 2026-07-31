//! Colours, fonts and spacing shared by every track in a figure.

/// Visual settings applied to a whole [`Figure`](crate::Figure).
///
/// Every field is public: start from [`Theme::light`] or [`Theme::dark`] and
/// overwrite what you need, rather than building one from scratch.
#[derive(Debug, Clone, PartialEq)]
pub struct Theme {
    /// Page background. Set it to `"none"` for a transparent SVG.
    pub background: String,
    /// Titles and tick labels.
    pub foreground: String,
    /// Secondary text: coordinates, track labels, legends.
    pub muted: String,
    /// Axis lines, baselines and tick marks.
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
    /// Helvetica and Arial metrics of [`Theme::font_family`]; change both
    /// together or logo letters will sit slightly proud of their boxes.
    pub cap_height_ratio: f64,
}

impl Theme {
    /// Dark ink on a white page, for manuscripts and slides.
    pub fn light() -> Self {
        Theme {
            background: "#ffffff".into(),
            foreground: "#1b1f23".into(),
            muted: "#6b7280".into(),
            rule: "#c8ccd1".into(),
            accent: "#0072b2".into(),
            palette: okabe_ito(),
            bases: BaseColors::default(),
            insertion: "#8e44ad".into(),
            corner_radius: 2.0,
            font_family: "Helvetica, Arial, sans-serif".into(),
            font_size: 11.0,
            label_font_size: 11.0,
            title_font_size: 14.0,
            cap_height_ratio: 0.72,
        }
    }

    /// Light ink on a dark page, for terminals and dark-mode documents.
    pub fn dark() -> Self {
        Theme {
            background: "#14181d".into(),
            foreground: "#e6edf3".into(),
            muted: "#9aa4b0".into(),
            rule: "#3a424c".into(),
            accent: "#3987e5".into(),
            palette: okabe_ito_dark(),
            bases: BaseColors::default(),
            insertion: "#8e44ad".into(),
            corner_radius: 2.0,
            font_family: "Helvetica, Arial, sans-serif".into(),
            font_size: 11.0,
            label_font_size: 11.0,
            title_font_size: 14.0,
            cap_height_ratio: 0.72,
        }
    }

    /// Colour `index` of the categorical palette, wrapping around at the end.
    pub fn color(&self, index: usize) -> &str {
        if self.palette.is_empty() {
            return &self.accent;
        }
        &self.palette[index % self.palette.len()]
    }
}

impl Default for Theme {
    fn default() -> Self {
        Theme::light()
    }
}

/// Six hues from the Okabe-Ito palette, ordered and trimmed so that every pair
/// stays apart under colour vision deficiency, not just neighbouring ones.
///
/// Six is where the measurement stopped, not where the eye got bored. Every
/// pair of these clears a CVD separation of 6.7 in OKLab hundredths, inside the
/// band that is sound when identity is also carried by something other than
/// colour, which here is always a legend or a letter. A seventh hue could not
/// be added without some pair collapsing: an olive against the vermillion came
/// out at 1.8 under protanopia, indistinguishable.
fn okabe_ito() -> Vec<String> {
    [
        "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#7b3294",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

/// The dark mode steps, chosen against the dark surface rather than lightened
/// from the light ones.
///
/// A dark background wants a narrower lightness band than a light one, so an
/// automatic flip of the light palette fails: half of it lands outside the band
/// and two of the entries stop being distinguishable. These were stepped for
/// the dark surface and validated there.
fn okabe_ito_dark() -> Vec<String> {
    [
        "#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9",
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
    /// is that green now means adenine rather than the thymine-red a reader
    /// expects, so label it or say so in the caption.
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
    let t = t.clamp(0.0, 1.0);
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

/// Black or white, whichever stays readable on top of `color`.
///
/// Feature labels are drawn inside coloured boxes whose colour the caller
/// chooses, so the ink has to be picked per box rather than fixed by the theme.
/// Anything that is not a `#rrggbb` string is assumed dark.
pub fn contrast_ink(color: &str) -> &'static str {
    let Some(hex) = color.strip_prefix('#') else {
        return "#ffffff";
    };
    if hex.len() != 6 {
        return "#ffffff";
    }
    let channel = |range: std::ops::Range<usize>| u8::from_str_radix(&hex[range], 16).ok();
    let (Some(r), Some(g), Some(b)) = (channel(0..2), channel(2..4), channel(4..6)) else {
        return "#ffffff";
    };
    // Rec. 601 luma, good enough to choose between two inks.
    let luma = 0.299 * r as f64 + 0.587 * g as f64 + 0.114 * b as f64;
    if luma > 150.0 {
        "#1b1f23"
    } else {
        "#ffffff"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ink_flips_with_the_brightness_of_the_box() {
        assert_eq!(contrast_ink("#000000"), "#ffffff");
        assert_eq!(contrast_ink("#0072b2"), "#ffffff");
        assert_eq!(contrast_ink("#ffffff"), "#1b1f23");
        assert_eq!(contrast_ink("#f0e442"), "#1b1f23");
    }

    #[test]
    fn ink_falls_back_to_white_for_colours_it_cannot_read() {
        assert_eq!(contrast_ink("red"), "#ffffff");
        assert_eq!(contrast_ink("#fff"), "#ffffff");
        assert_eq!(contrast_ink("#zzzzzz"), "#ffffff");
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
}
