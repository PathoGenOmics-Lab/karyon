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
            accent: "#56b4e9".into(),
            palette: okabe_ito_dark(),
            bases: BaseColors::default(),
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

/// The Okabe-Ito qualitative palette, ordered so the first colours stay
/// distinguishable under the common forms of colour vision deficiency.
fn okabe_ito() -> Vec<String> {
    [
        "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9", "#8c564b", "#606060",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

/// Same hues, lifted so they hold up against a dark background.
fn okabe_ito_dark() -> Vec<String> {
    [
        "#56b4e9", "#e8833a", "#33c295", "#e39ec1", "#f0bf6c", "#8fd3f4", "#c39a8c", "#b6bec9",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

/// Colours for the four nucleotides plus a fallback for everything else.
///
/// The defaults follow the IGV convention, which most readers of a genome
/// figure already have in their heads. Green and red are not a colourblind-safe
/// pair; swap the fields if that matters more than matching convention.
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

impl Default for BaseColors {
    fn default() -> Self {
        BaseColors {
            a: "#33a02c".into(),
            c: "#1f78b4".into(),
            g: "#e08214".into(),
            t: "#e31a1c".into(),
            other: "#9e9e9e".into(),
        }
    }
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
