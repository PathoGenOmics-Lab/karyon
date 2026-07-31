//! Several drawings in one document.
//!
//! A [`Figure`](crate::Figure) is one stack of tracks over one coordinate
//! system, and a [`Rings`](crate::Rings) plot is one circular sequence. A paper
//! figure is usually several of those with letters on them, and so is an
//! overview of everything a library can draw. Panels stack finished drawings
//! into a single SVG without any of them having to know about the others, and
//! it takes anything that implements [`Drawing`], so a linear
//! stack and a circle go on the same sheet.
//!
//! Each figure is embedded as a nested `<svg>` inside a translated group, so
//! nothing is reparsed and nothing is moved around afterwards. The one thing a
//! panel does not keep to itself is its ids: an id in SVG belongs to the whole
//! document, so every figure is rendered with a prefix of its own and its clips
//! go on meaning what they meant.

use std::fs;
use std::io;
use std::path::Path;

use crate::rings::Drawing;
use crate::svg::{escape, num, Anchor, SvgWriter};
use crate::theme::Theme;

/// One figure placed on a sheet.
struct Panel {
    svg: String,
    width: f64,
    height: f64,
    label: Option<String>,
    caption: Option<String>,
}

/// A column of figures rendered into one document.
///
/// ```
/// use karyon::{AxisTrack, CoverageTrack, Figure, Panels, Region};
///
/// let region = Region::new("chr1", 0, 1_000).unwrap();
/// let depth: Vec<f64> = (0..1_000).map(|i| (i % 40) as f64).collect();
///
/// let svg = Panels::new()
///     .push(&Figure::new(region.clone()).push(CoverageTrack::new(0, depth)), "A")
///     .push(&Figure::new(region).push(AxisTrack::new()), "B")
///     .to_svg();
///
/// // Two nested figures inside one document.
/// assert_eq!(svg.matches("<svg").count(), 3);
/// ```
pub struct Panels {
    panels: Vec<Panel>,
    gap: f64,
    margin: f64,
    theme: Theme,
    title: Option<String>,
}

impl Panels {
    /// An empty sheet.
    pub fn new() -> Self {
        Panels {
            panels: Vec::new(),
            gap: 18.0,
            margin: 14.0,
            theme: Theme::light(),
            title: None,
        }
    }

    /// Sets the vertical gap between panels.
    pub fn gap(mut self, gap: f64) -> Self {
        self.gap = gap.max(0.0);
        self
    }

    /// Sets the whitespace around the whole sheet.
    pub fn margin(mut self, margin: f64) -> Self {
        self.margin = margin.max(0.0);
        self
    }

    /// Replaces the theme, which is used for the sheet's own title and letters
    /// rather than for the panels, each of which keeps its own.
    pub fn theme(mut self, theme: Theme) -> Self {
        self.theme = theme;
        self
    }

    /// Sets a title across the top of the sheet.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Adds a figure with a panel letter.
    ///
    /// The figure is rendered now and stored, so later changes to it do not
    /// reach the sheet.
    pub fn push(self, figure: &impl Drawing, label: impl Into<String>) -> Self {
        self.add(figure, Some(label.into()), None)
    }

    /// Adds a figure with a letter and a caption under it.
    pub fn push_captioned(
        self,
        figure: &impl Drawing,
        label: impl Into<String>,
        caption: impl Into<String>,
    ) -> Self {
        self.add(figure, Some(label.into()), Some(caption.into()))
    }

    /// Adds a figure with no letter.
    pub fn push_bare(self, figure: &impl Drawing) -> Self {
        self.add(figure, None, None)
    }

    fn add(
        mut self,
        figure: &impl Drawing,
        label: Option<String>,
        caption: Option<String>,
    ) -> Self {
        let (width, height) = figure.dimensions();
        // Each panel gets an id space of its own. Without this the second
        // panel's `url(#karyon-clip-0)` would resolve to the first panel's
        // clipping rectangle, and its tracks would be cropped to a band
        // belonging to a different figure.
        let prefix = format!("p{}-", self.panels.len());
        self.panels.push(Panel {
            svg: figure.to_svg_with_id_prefix(&prefix),
            width,
            height,
            label,
            caption,
        });
        self
    }

    /// Width of the strip on the left holding the panel letters.
    ///
    /// A letter drawn on top of a panel is a letter drawn on top of data, and
    /// on top of an opaque page colour it is a letter that is not there at all,
    /// so it gets a column of its own.
    fn letter_gutter(&self) -> f64 {
        if self.panels.iter().any(|panel| panel.label.is_some()) {
            self.theme.title_font_size + 6.0
        } else {
            0.0
        }
    }

    /// How many panels the sheet holds.
    pub fn len(&self) -> usize {
        self.panels.len()
    }

    /// Whether the sheet is empty.
    pub fn is_empty(&self) -> bool {
        self.panels.is_empty()
    }

    /// Width and height of the finished sheet.
    pub fn dimensions(&self) -> (f64, f64) {
        let width = self
            .panels
            .iter()
            .map(|panel| panel.width)
            .fold(0.0f64, f64::max)
            + self.letter_gutter()
            + self.margin * 2.0;
        let caption_room = self.theme.font_size + 4.0;
        let content: f64 = self
            .panels
            .iter()
            .map(|panel| {
                panel.height
                    + if panel.caption.is_some() {
                        caption_room
                    } else {
                        0.0
                    }
            })
            .sum();
        let gaps = self.gap * (self.panels.len().saturating_sub(1)) as f64;
        let header = if self.title.is_some() {
            self.theme.title_font_size + 12.0
        } else {
            0.0
        };
        (
            width.max(1.0),
            (content + gaps + header + self.margin * 2.0).max(1.0),
        )
    }

    /// Renders the sheet.
    pub fn to_svg(&self) -> String {
        let (width, height) = self.dimensions();
        let gutter = self.letter_gutter();
        let left = self.margin + gutter;
        let mut svg = SvgWriter::new();

        let mut y = self.margin;
        if let Some(title) = &self.title {
            svg.text_bold(
                self.margin,
                y + self.theme.title_font_size,
                title,
                &self.theme.foreground,
                self.theme.title_font_size,
                Anchor::Start,
            );
            y += self.theme.title_font_size + 12.0;
        }

        let mut body = String::new();
        for panel in &self.panels {
            // The panel goes in untouched, inside a group that moves it, and it
            // goes in first: a figure paints its own page colour, so anything
            // written before it would end up underneath and invisible.
            body.push_str(&format!(
                r#"<g transform="translate({} {})">{}</g>"#,
                num(left),
                num(y),
                panel.svg
            ));

            if let Some(label) = &panel.label {
                svg.text_bold(
                    self.margin,
                    y + self.theme.title_font_size,
                    label,
                    &self.theme.foreground,
                    self.theme.title_font_size,
                    Anchor::Start,
                );
            }
            y += panel.height;

            if let Some(caption) = &panel.caption {
                svg.text(
                    left,
                    y + self.theme.font_size,
                    caption,
                    &self.theme.muted,
                    self.theme.font_size,
                    Anchor::Start,
                );
                y += self.theme.font_size + 4.0;
            }
            y += self.gap;
        }

        // Assembled by hand rather than through `finish`, because the sheet's
        // own lettering has to sit over the panels and the panels are not the
        // writer's to hold.
        let (defs, overlay) = svg.into_parts();
        let mut out = String::with_capacity(body.len() + overlay.len() + 512);
        out.push_str(&format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}" font-family="{}">"#,
            num(width),
            num(height),
            num(width),
            num(height),
            escape(&self.theme.font_family)
        ));
        if !defs.is_empty() {
            out.push_str(&format!("<defs>{defs}</defs>"));
        }
        if self.theme.background != "none" {
            out.push_str(&format!(
                r#"<rect x="0" y="0" width="{}" height="{}" fill="{}"/>"#,
                num(width),
                num(height),
                self.theme.background
            ));
        }
        out.push_str(&body);
        out.push_str(&overlay);
        out.push_str("</svg>");
        out
    }

    /// Renders the sheet and writes it to `path`.
    ///
    /// # Errors
    ///
    /// Returns whatever [`fs::write`] returns.
    pub fn save_svg(&self, path: impl AsRef<Path>) -> io::Result<()> {
        fs::write(path, self.to_svg())
    }
}

impl Default for Panels {
    fn default() -> Self {
        Panels::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;
    use crate::track::{AxisTrack, CoverageTrack};

    fn figure() -> Figure {
        let region = Region::new("chr1", 0, 1_000).unwrap();
        Figure::new(region)
            .show_region_label(false)
            .push(CoverageTrack::new(0, vec![5.0; 1_000]))
    }

    #[test]
    fn an_empty_sheet_is_a_valid_document() {
        let svg = Panels::new().to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(svg.ends_with("</svg>"));
        assert!(Panels::new().is_empty());
        assert_eq!(Panels::new().len(), 0);
    }

    #[test]
    fn every_panel_is_nested_whole() {
        let sheet = Panels::new().push(&figure(), "A").push(&figure(), "B");
        let svg = sheet.to_svg();
        // The sheet's own root, plus one nested root per panel.
        assert_eq!(svg.matches("<svg ").count(), 3);
        assert_eq!(svg.matches("</svg>").count(), 3);
        assert_eq!(sheet.len(), 2);
    }

    #[test]
    fn a_panel_is_the_same_picture_on_a_sheet_as_on_its_own() {
        // Its ids are namespaced, and that is the only difference.
        let alone = figure().to_svg_with_id_prefix("p0-");
        let sheet = Panels::new().push_bare(&figure()).to_svg();
        assert!(
            sheet.contains(&alone),
            "the panel was rewritten on its way onto the sheet"
        );
    }

    #[test]
    fn two_panels_never_claim_the_same_id() {
        // An id in SVG is document-wide, so a duplicate does not merely look
        // untidy: `url(#id)` resolves to the first match, and the second
        // panel's tracks would be clipped to a band from the first one.
        let svg = Panels::new()
            .push(&figure(), "A")
            .push(&figure(), "B")
            .to_svg();
        let ids: Vec<&str> = svg
            .match_indices(r#" id=""#)
            .map(|(index, prefix)| {
                let rest = &svg[index + prefix.len()..];
                &rest[..rest.find('"').unwrap()]
            })
            .collect();
        assert!(!ids.is_empty(), "the figures do clip, so there are ids");
        let mut unique = ids.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), ids.len(), "{ids:?}");
    }

    #[test]
    fn a_letter_is_drawn_over_the_sheet_rather_than_under_a_panel() {
        let svg = Panels::new().push(&figure(), "A").to_svg();
        let panel = svg.find("<g transform=\"translate(").unwrap();
        let letter = svg.find(">A</text>").unwrap();
        // A figure paints its own opaque page colour. A letter written before
        // the panel is a letter nobody ever sees.
        assert!(letter > panel, "the letter went in underneath the panel");
    }

    #[test]
    fn letters_get_a_column_of_their_own() {
        let bare = Panels::new().push_bare(&figure()).dimensions().0;
        let lettered = Panels::new().push(&figure(), "A").dimensions().0;
        assert!(
            lettered > bare,
            "a letter on top of the panel is a letter on top of data"
        );
    }

    #[test]
    fn panels_are_moved_rather_than_redrawn() {
        let svg = Panels::new()
            .push(&figure(), "A")
            .push(&figure(), "B")
            .to_svg();
        assert_eq!(svg.matches("<g transform=\"translate(").count(), 2);
    }

    #[test]
    fn the_sheet_is_as_tall_as_its_panels_and_gaps() {
        let (_, one) = Panels::new().push_bare(&figure()).dimensions();
        let (_, two) = Panels::new()
            .push_bare(&figure())
            .push_bare(&figure())
            .dimensions();
        let (_, panel) = figure().dimensions();
        assert!(
            (two - one - panel - 18.0).abs() < 1e-9,
            "one gap between two"
        );
    }

    #[test]
    fn the_sheet_is_as_wide_as_its_widest_panel() {
        let narrow = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .width(400.0)
            .push(AxisTrack::new());
        let wide = Figure::new(Region::new("chr1", 0, 100).unwrap())
            .width(900.0)
            .push(AxisTrack::new());
        let (width, _) = Panels::new()
            .push_bare(&narrow)
            .push_bare(&wide)
            .margin(10.0)
            .dimensions();
        assert_eq!(width, 900.0 + 20.0);
    }

    #[test]
    fn letters_and_captions_are_drawn() {
        let svg = Panels::new()
            .title("everything")
            .push_captioned(&figure(), "A", "a coverage profile")
            .to_svg();
        assert!(svg.contains(">A</text>"));
        assert!(svg.contains("a coverage profile"));
        assert!(svg.contains(">everything</text>"));
    }

    #[test]
    fn a_caption_makes_the_sheet_taller() {
        let bare = Panels::new().push_bare(&figure()).dimensions().1;
        let captioned = Panels::new()
            .push_captioned(&figure(), "A", "words")
            .dimensions()
            .1;
        assert!(captioned > bare);
    }

    #[test]
    fn a_later_change_to_a_figure_does_not_reach_the_sheet() {
        let first = figure();
        let sheet = Panels::new().push_bare(&first);
        let before = sheet.to_svg();
        // The figure is consumed by value elsewhere, but the point stands: the
        // sheet holds a rendering, not a reference.
        let after = sheet.to_svg();
        assert_eq!(before, after);
    }
}
