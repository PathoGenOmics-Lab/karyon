//! The whole chromosome, and a mark showing which part of it you are looking
//! at.
//!
//! Every other track in this crate maps its data through the shared
//! [`Scale`], which is what keeps their x axes aligned. This one deliberately
//! does not. An ideogram exists to answer "where am I", and a track that only
//! showed the region on display could not answer it: it would be a picture of
//! the window, drawn inside the window. So it draws the whole chromosome across
//! the plotting area and puts a marker where the region sits.
//!
//! That makes it the one track whose horizontal position means something
//! different from its neighbours', which is worth knowing before you put a
//! ruler under it.

use crate::scale::Scale;
use crate::svg::{num, text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// How far the region marker stands proud of the chromosome, in pixels.
///
/// The chromosome is inset by this much inside the band so the overhang has
/// somewhere to go: a track that draws outside its band is a track whose
/// figure clips the difference away.
const MARKER_PAD: f64 = 2.0;

/// How a cytogenetic band was stained, which is what decides its shade.
///
/// The names are the ones in the UCSC `cytoBand` table, so a row of that file
/// converts without a lookup table of your own.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stain {
    /// Unstained, the palest band.
    Gneg,
    /// Lightly stained.
    Gpos25,
    /// Half stained.
    Gpos50,
    /// Heavily stained.
    Gpos75,
    /// Fully stained, the darkest band.
    Gpos100,
    /// The centromere.
    Acen,
    /// A variable region, conventionally drawn lighter than its neighbours.
    Gvar,
    /// An acrocentric stalk.
    Stalk,
}

impl Stain {
    /// Reads the `gieStain` column of a UCSC `cytoBand` row.
    ///
    /// Unknown values become [`Stain::Gneg`], the palest shade, on the grounds
    /// that inventing a dark band is worse than drawing nothing.
    pub fn from_name(name: &str) -> Self {
        match name.trim().to_ascii_lowercase().as_str() {
            "gpos25" => Stain::Gpos25,
            "gpos50" => Stain::Gpos50,
            "gpos75" => Stain::Gpos75,
            "gpos" | "gpos100" => Stain::Gpos100,
            "acen" => Stain::Acen,
            "gvar" => Stain::Gvar,
            "stalk" => Stain::Stalk,
            _ => Stain::Gneg,
        }
    }

    /// How dark the band is, from 0 for the page colour to 1 for full ink.
    ///
    /// The centromere has no shade of its own because it is drawn in its own
    /// colour rather than on the grey ladder.
    pub fn shade(self) -> f64 {
        match self {
            Stain::Gneg => 0.0,
            Stain::Gpos25 => 0.25,
            Stain::Gpos50 => 0.5,
            Stain::Gpos75 => 0.72,
            Stain::Gpos100 => 0.9,
            Stain::Gvar => 0.18,
            Stain::Stalk => 0.35,
            Stain::Acen => 0.0,
        }
    }

    /// Whether this band is part of the centromere.
    pub fn is_centromere(self) -> bool {
        matches!(self, Stain::Acen)
    }
}

/// One cytogenetic band, in 0-based half-open coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct Band {
    /// First base, 0-based.
    pub start: u64,
    /// One past the last base.
    pub end: u64,
    /// Band name, such as `p31.3`.
    pub name: Option<String>,
    /// How it was stained.
    pub stain: Stain,
}

impl Band {
    /// A band spanning `start..end`.
    pub fn new(start: u64, end: u64, stain: Stain) -> Self {
        Band {
            start,
            end: end.max(start),
            name: None,
            stain,
        }
    }

    /// Sets the band name.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Length in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Whether the band covers no bases.
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

/// A chromosome drawn end to end, with the region on display marked on it.
///
/// ```
/// use karyon::{Band, Figure, IdeogramTrack, Region, Stain};
///
/// let bands = vec![
///     Band::new(0, 2_300_000, Stain::Gneg).name("p13"),
///     Band::new(2_300_000, 2_800_000, Stain::Acen),
///     Band::new(2_800_000, 5_000_000, Stain::Gpos75).name("q11"),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 3_000_000, 3_010_000).unwrap())
///     .push(IdeogramTrack::new(5_000_000, bands).label("chr1"))
///     .to_svg();
/// assert!(svg.contains("<path"));
/// ```
#[derive(Debug, Clone)]
pub struct IdeogramTrack {
    length: u64,
    bands: Vec<Band>,
    label: Option<String>,
    height: f64,
    centromere_color: Option<String>,
    outline: Option<String>,
    marker_color: Option<String>,
    show_marker: bool,
    highlight: Option<(u64, u64)>,
    show_band_names: bool,
}

impl IdeogramTrack {
    /// A chromosome `length` bases long, made of `bands`.
    ///
    /// Bands need not tile the chromosome: gaps simply draw as the page colour,
    /// and an empty list gives a bare outline, which is what a bacterial genome
    /// with no cytogenetics to speak of will want.
    pub fn new(length: u64, bands: impl Into<Vec<Band>>) -> Self {
        IdeogramTrack {
            length: length.max(1),
            bands: bands.into(),
            label: None,
            height: 22.0,
            centromere_color: None,
            outline: None,
            marker_color: None,
            show_marker: true,
            highlight: None,
            show_band_names: false,
        }
    }

    /// A bare chromosome with no banding.
    pub fn bare(length: u64) -> Self {
        IdeogramTrack::new(length, Vec::new())
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(4.0);
        self
    }

    /// Sets the colour of the centromere.
    pub fn centromere_color(mut self, color: impl Into<String>) -> Self {
        self.centromere_color = Some(color.into());
        self
    }

    /// Sets the colour of the chromosome outline.
    pub fn outline(mut self, color: impl Into<String>) -> Self {
        self.outline = Some(color.into());
        self
    }

    /// Sets the colour of the marker showing the region on display.
    pub fn marker_color(mut self, color: impl Into<String>) -> Self {
        self.marker_color = Some(color.into());
        self
    }

    /// Shows or hides that marker.
    pub fn show_marker(mut self, show: bool) -> Self {
        self.show_marker = show;
        self
    }

    /// Marks a span other than the region the figure is showing.
    ///
    /// Useful when the ideogram is context for something the figure does not
    /// otherwise draw: a whole gene, say, while the tracks below show one exon
    /// of it.
    pub fn highlight(mut self, start: u64, end: u64) -> Self {
        self.highlight = Some((start, end.max(start + 1)));
        self
    }

    /// Draws band names under the chromosome where they fit.
    pub fn show_band_names(mut self, show: bool) -> Self {
        self.show_band_names = show;
        self
    }

    /// Length of the chromosome in bases.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// The bands, in the order given.
    pub fn bands(&self) -> &[Band] {
        &self.bands
    }

    /// Where the centromere sits, as the span of the `acen` bands.
    ///
    /// `None` when no band is stained as a centromere, which is the normal
    /// answer for a bacterial chromosome and for an unbanded ideogram.
    pub fn centromere(&self) -> Option<(u64, u64)> {
        let mut span: Option<(u64, u64)> = None;
        for band in self.bands.iter().filter(|b| b.stain.is_centromere()) {
            span = Some(match span {
                Some((start, end)) => (start.min(band.start), end.max(band.end)),
                None => (band.start, band.end),
            });
        }
        span
    }

    /// Fraction of the chromosome at a position, clamped to the ends.
    fn fraction(&self, position: u64) -> f64 {
        (position as f64 / self.length as f64).clamp(0.0, 1.0)
    }

    /// Colour of a band under a theme.
    fn band_color(&self, band: &Band, theme: &Theme) -> String {
        if band.stain.is_centromere() {
            return self
                .centromere_color
                .clone()
                .unwrap_or_else(|| "#b03a2e".to_string());
        }
        // The grey ladder is mixed from the theme's own ink and page, so a dark
        // figure gets a dark ladder rather than a light one left inverted.
        mix(theme.surface(), &theme.foreground, band.stain.shade())
    }
}

/// The outline of one arm: a rectangle with either end optionally rounded off
/// into a half circle.
fn arm_path(x: f64, y: f64, w: f64, h: f64, round_left: bool, round_right: bool) -> String {
    if w <= 0.0 || h <= 0.0 {
        return String::new();
    }
    let radius = (h / 2.0).min(w / 2.0);
    let (left, right) = (x, x + w);
    let (top, bottom) = (y, y + h);
    let mut d = String::with_capacity(96);

    if round_left {
        d.push_str(&format!("M{} {}", num(left + radius), num(top)));
    } else {
        d.push_str(&format!("M{} {}", num(left), num(top)));
    }

    if round_right {
        d.push_str(&format!(" L{} {}", num(right - radius), num(top)));
        d.push_str(&format!(
            " A{} {} 0 0 1 {} {}",
            num(radius),
            num(radius),
            num(right - radius),
            num(bottom)
        ));
    } else {
        d.push_str(&format!(" L{} {}", num(right), num(top)));
        d.push_str(&format!(" L{} {}", num(right), num(bottom)));
    }

    if round_left {
        d.push_str(&format!(" L{} {}", num(left + radius), num(bottom)));
        d.push_str(&format!(
            " A{} {} 0 0 1 {} {}",
            num(radius),
            num(radius),
            num(left + radius),
            num(top)
        ));
    } else {
        d.push_str(&format!(" L{} {}", num(left), num(bottom)));
        d.push_str(&format!(" L{} {}", num(left), num(top)));
    }

    d.push('Z');
    d
}

impl Track for IdeogramTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band_area = ctx.band;
        // Room under the chromosome for band names, when they are wanted, and
        // room above and below it for the marker to stand proud of it. Without
        // the second, the marker is drawn outside the band and the figure's
        // clip takes the overhang off again, so the pointer reads as flush and
        // the two pixels of intent never reach the page.
        let pad = if self.show_marker { MARKER_PAD } else { 0.0 };
        let names_room = if self.show_band_names {
            ctx.theme.font_size + 2.0
        } else {
            0.0
        };
        let body_height = (self.height - names_room - pad * 2.0).max(6.0);
        let top = band_area.y + pad;
        let bottom = top + body_height;
        let left = band_area.x;
        let width = band_area.w;
        let outline = self
            .outline
            .clone()
            .unwrap_or_else(|| ctx.theme.muted.clone());

        // The chromosome spans the plotting area regardless of the region, which
        // is the whole point of the track.
        let at = |position: u64| left + self.fraction(position) * width;

        let arms: Vec<(f64, f64, bool, bool)> = match self.centromere() {
            Some((start, end)) => {
                let (p_end, q_start) = (at(start), at(end));
                vec![
                    (left, p_end - left, true, false),
                    (q_start, left + width - q_start, false, true),
                ]
            }
            None => vec![(left, width, true, true)],
        };

        for (arm_x, arm_w, round_left, round_right) in &arms {
            if *arm_w <= 0.0 {
                continue;
            }
            let path = arm_path(*arm_x, top, *arm_w, body_height, *round_left, *round_right);
            if path.is_empty() {
                continue;
            }

            // Bands are painted inside the arm's own silhouette, so a rounded
            // end stays rounded however the banding falls across it.
            ctx.svg.begin_clip_path(&path);
            for cytoband in &self.bands {
                if cytoband.stain.is_centromere() || cytoband.is_empty() {
                    continue;
                }
                let x0 = at(cytoband.start);
                let x1 = at(cytoband.end);
                if x1 <= *arm_x || x0 >= arm_x + arm_w {
                    continue;
                }
                let color = self.band_color(cytoband, ctx.theme);
                ctx.svg
                    .rect(x0, top, (x1 - x0).max(0.4), body_height, &color);
            }
            ctx.svg.end_group();
            ctx.svg.path_stroked(&path, &outline, 1.0);
        }

        // The centromere, drawn as its own shape rather than as a band, so the
        // waist of the chromosome reads as a waist.
        if let Some((start, end)) = self.centromere() {
            let (x0, x1) = (at(start), at(end));
            let middle = (top + bottom) / 2.0;
            let centre = (x0 + x1) / 2.0;
            let color = self.band_color(&Band::new(start, end, Stain::Acen), ctx.theme);
            // An hourglass with the pinch in the middle, which is what makes a
            // chromosome read as having a waist. Both edges converge on the
            // centre and open out again, so the shape is symmetric about it.
            let waist = body_height * 0.3;
            ctx.svg.polygon(
                &[
                    (x0, top),
                    (centre, middle - waist / 2.0),
                    (x1, top),
                    (x1, bottom),
                    (centre, middle + waist / 2.0),
                    (x0, bottom),
                ],
                &color,
            );
        }

        if self.show_band_names {
            let font = ctx.theme.font_size - 1.0;
            let baseline = bottom + font + 1.0;
            for cytoband in &self.bands {
                let Some(name) = &cytoband.name else {
                    continue;
                };
                let x0 = at(cytoband.start);
                let x1 = at(cytoband.end);
                if x1 - x0 < text_width(name, font) + 4.0 {
                    continue;
                }
                ctx.svg.text(
                    (x0 + x1) / 2.0,
                    baseline,
                    name,
                    &ctx.theme.muted,
                    font,
                    Anchor::Middle,
                );
            }
        }

        if self.show_marker {
            let (start, end) = self
                .highlight
                .unwrap_or((ctx.region.start(), ctx.region.end()));
            let color = self
                .marker_color
                .clone()
                .unwrap_or_else(|| "#d7263d".to_string());
            let x0 = at(start);
            // A ten kilobase window on a chromosome is a fraction of a pixel, so
            // the marker has a minimum width. It is a pointer, not a
            // measurement, and one too thin to see would be neither.
            let x1 = at(end).max(x0 + 3.0);
            ctx.svg.rect_opacity(
                x0,
                top - MARKER_PAD,
                x1 - x0,
                body_height + MARKER_PAD * 2.0,
                &color,
                0.22,
            );
            ctx.svg.rect_outline(
                x0,
                top - MARKER_PAD,
                x1 - x0,
                body_height + MARKER_PAD * 2.0,
                &color,
                1.2,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::region::Region;

    fn banded() -> Vec<Band> {
        vec![
            Band::new(0, 2_000_000, Stain::Gneg).name("p13"),
            Band::new(2_000_000, 2_400_000, Stain::Gpos75).name("p12"),
            Band::new(2_400_000, 2_600_000, Stain::Acen),
            Band::new(2_600_000, 4_000_000, Stain::Gpos50).name("q11"),
            Band::new(4_000_000, 5_000_000, Stain::Gneg).name("q12"),
        ]
    }

    fn region() -> Region {
        Region::new("chr1", 3_000_000, 3_010_000).unwrap()
    }

    #[test]
    fn stain_names_come_straight_from_the_ucsc_column() {
        assert_eq!(Stain::from_name("gneg"), Stain::Gneg);
        assert_eq!(Stain::from_name("gpos25"), Stain::Gpos25);
        assert_eq!(Stain::from_name("gpos100"), Stain::Gpos100);
        assert_eq!(Stain::from_name("gpos"), Stain::Gpos100);
        assert_eq!(Stain::from_name("ACEN"), Stain::Acen);
        assert_eq!(Stain::from_name(" stalk "), Stain::Stalk);
        assert_eq!(Stain::from_name("gvar"), Stain::Gvar);
    }

    #[test]
    fn an_unknown_stain_is_pale_rather_than_invented() {
        assert_eq!(Stain::from_name("whatever"), Stain::Gneg);
        assert_eq!(Stain::from_name(""), Stain::Gneg);
        assert_eq!(Stain::Gneg.shade(), 0.0);
    }

    #[test]
    fn the_grey_ladder_climbs() {
        let ladder = [
            Stain::Gneg.shade(),
            Stain::Gpos25.shade(),
            Stain::Gpos50.shade(),
            Stain::Gpos75.shade(),
            Stain::Gpos100.shade(),
        ];
        assert!(
            ladder.windows(2).all(|w| w[0] < w[1]),
            "not monotone: {ladder:?}"
        );
        assert!(ladder.iter().all(|s| (0.0..=1.0).contains(s)));
    }

    #[test]
    fn the_centromere_is_the_span_of_the_acen_bands() {
        let track = IdeogramTrack::new(5_000_000, banded());
        assert_eq!(track.centromere(), Some((2_400_000, 2_600_000)));
    }

    #[test]
    fn two_acen_bands_merge_into_one_span() {
        let bands = vec![
            Band::new(1_000, 2_000, Stain::Acen),
            Band::new(2_000, 3_500, Stain::Acen),
        ];
        let track = IdeogramTrack::new(10_000, bands);
        assert_eq!(track.centromere(), Some((1_000, 3_500)));
    }

    #[test]
    fn a_chromosome_without_a_centromere_has_none() {
        assert_eq!(IdeogramTrack::bare(4_400_000).centromere(), None);
        let unbanded = vec![Band::new(0, 100, Stain::Gneg)];
        assert_eq!(IdeogramTrack::new(1_000, unbanded).centromere(), None);
    }

    #[test]
    fn positions_map_to_fractions_of_the_whole_chromosome() {
        let track = IdeogramTrack::bare(1_000);
        assert_eq!(track.fraction(0), 0.0);
        assert_eq!(track.fraction(500), 0.5);
        assert_eq!(track.fraction(1_000), 1.0);
        // Past the end is clamped rather than drawn off the chromosome.
        assert_eq!(track.fraction(9_999), 1.0);
    }

    #[test]
    fn a_zero_length_chromosome_cannot_happen() {
        assert_eq!(IdeogramTrack::bare(0).length(), 1);
    }

    #[test]
    fn the_shade_ladder_follows_the_theme_it_is_drawn_on() {
        let track = IdeogramTrack::new(100, banded());
        let dark_band = Band::new(0, 10, Stain::Gpos100);
        let on_light = track.band_color(&dark_band, &Theme::light());
        let on_dark = track.band_color(&dark_band, &Theme::dark());
        // The darkest band is near the ink of its own page, which is dark on a
        // light theme and light on a dark one.
        assert_ne!(on_light, on_dark);
        assert!(on_light.starts_with('#') && on_dark.starts_with('#'));
    }

    #[test]
    fn the_centromere_keeps_its_own_colour_on_either_theme() {
        let track = IdeogramTrack::new(100, banded()).centromere_color("#123456");
        let acen = Band::new(0, 10, Stain::Acen);
        assert_eq!(track.band_color(&acen, &Theme::light()), "#123456");
        assert_eq!(track.band_color(&acen, &Theme::dark()), "#123456");
    }

    #[test]
    fn an_arm_can_be_rounded_at_one_end_only() {
        let both = arm_path(0.0, 0.0, 100.0, 20.0, true, true);
        let left = arm_path(0.0, 0.0, 100.0, 20.0, true, false);
        let neither = arm_path(0.0, 0.0, 100.0, 20.0, false, false);
        assert_eq!(both.matches('A').count(), 2);
        assert_eq!(left.matches('A').count(), 1);
        assert_eq!(neither.matches('A').count(), 0);
        for path in [both, left, neither] {
            assert!(path.starts_with('M') && path.ends_with('Z'));
        }
    }

    #[test]
    fn a_degenerate_arm_produces_no_path() {
        assert!(arm_path(0.0, 0.0, 0.0, 20.0, true, true).is_empty());
        assert!(arm_path(0.0, 0.0, 100.0, 0.0, true, true).is_empty());
    }

    #[test]
    fn the_ideogram_ignores_the_shared_scale_and_draws_the_whole_chromosome() {
        // Two figures on the same chromosome, showing windows at opposite ends.
        // The chromosome must come out identical; only the marker moves.
        let left_window = Region::new("chr1", 100_000, 110_000).unwrap();
        let right_window = Region::new("chr1", 4_800_000, 4_810_000).unwrap();
        let svg_of = |region| {
            Figure::new(region)
                .show_region_label(false)
                .push(IdeogramTrack::new(5_000_000, banded()))
                .to_svg()
        };
        let left = svg_of(left_window);
        let right = svg_of(right_window);

        assert_ne!(left, right, "the marker has to move");
        // The arms are the same paths in both, since they do not depend on the
        // region at all.
        let arms = |svg: &str| {
            svg.match_indices("<path")
                .map(|(i, _)| svg[i..i + 120].to_string())
                .collect::<Vec<_>>()
        };
        assert_eq!(arms(&left), arms(&right));
    }

    #[test]
    fn a_tiny_window_still_gets_a_visible_marker() {
        // Ten kilobases of a five megabase chromosome is a fiftieth of a pixel.
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()))
            .to_svg();
        let widths: Vec<f64> = svg
            .match_indices(r#"<rect x="#)
            .filter_map(|(index, _)| {
                let element = &svg[index..index + svg[index..].find("/>")?];
                if !element.contains("#d7263d") {
                    return None;
                }
                let start = element.find(r#"width=""#)? + 7;
                let rest = &element[start..];
                rest[..rest.find('"')?].parse().ok()
            })
            .collect();
        assert!(!widths.is_empty(), "no marker drawn");
        assert!(widths.iter().all(|w| *w >= 3.0), "marker widths {widths:?}");
    }

    #[test]
    fn the_marker_can_be_pointed_somewhere_else_entirely() {
        let here = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()))
            .to_svg();
        let there = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()).highlight(100_000, 900_000))
            .to_svg();
        assert_ne!(here, there);
    }

    #[test]
    fn nothing_the_ideogram_draws_leaves_its_own_band() {
        // The marker stands proud of the chromosome on purpose. Drawn outside
        // the band, the figure's clip takes the overhang straight back off and
        // the pointer reads as flush, so the chromosome is inset instead.
        let region = Region::new("chr1", 3_000_000, 3_010_000).unwrap();
        let svg = Figure::new(region)
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()).height(18.0))
            .to_svg();

        let clip: Vec<f64> = svg
            .split("<clipPath")
            .nth(1)
            .unwrap()
            .split('"')
            .filter_map(|piece| piece.parse::<f64>().ok())
            .collect();
        let (top, height) = (clip[1], clip[3]);

        for chunk in svg.split("<rect ").skip(1) {
            let head = &chunk[..chunk.find("/>").unwrap()];
            let value = |key: &str| -> Option<f64> {
                let at = head.find(key)? + key.len();
                head[at..].split('"').next()?.parse().ok()
            };
            let (Some(y), Some(h)) = (value("y=\""), value("height=\"")) else {
                continue;
            };
            if h > height + 1.0 {
                continue; // the page background, which is not in the band
            }
            assert!(y >= top - 0.01, "a rect starts {} above the band", top - y);
            assert!(
                y + h <= top + height + 0.01,
                "a rect ends {} below the band",
                y + h - top - height
            );
        }
    }

    #[test]
    fn the_marker_can_be_turned_off() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()).show_marker(false))
            .to_svg();
        assert!(!svg.contains("#d7263d"));
    }

    #[test]
    fn a_bare_chromosome_is_one_arm_and_still_draws() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::bare(4_411_532).label("H37Rv"))
            .to_svg();
        // One silhouette, so one clip path from the track plus the figure's own.
        assert_eq!(svg.matches("<clipPath").count(), 2);
        assert_eq!(
            svg.matches(r#"<path d="#).count(),
            2,
            "the silhouette, once as a clip and once as an outline"
        );
        assert!(!svg.contains("NaN"));
    }

    #[test]
    fn a_band_narrower_than_its_name_goes_unlabelled() {
        // One roomy band and one hairline band, both named.
        let bands = vec![
            Band::new(0, 4_000_000, Stain::Gneg).name("q11"),
            Band::new(4_000_000, 4_010_000, Stain::Gpos75).name("q12.33"),
        ];
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, bands).show_band_names(true))
            .to_svg();
        assert!(svg.contains("q11"), "the roomy band should be labelled");
        assert!(
            !svg.contains("q12.33"),
            "a name wider than its band is dropped, not squeezed"
        );
    }

    #[test]
    fn band_names_are_off_unless_asked_for() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()))
            .to_svg();
        assert!(!svg.contains("q11"));
    }

    #[test]
    fn bands_are_clipped_to_the_silhouette_of_their_arm() {
        let svg = Figure::new(region())
            .show_region_label(false)
            .push(IdeogramTrack::new(5_000_000, banded()))
            .to_svg();
        // One clip path per arm, plus the one the figure puts round the track.
        assert_eq!(svg.matches("<clipPath").count(), 3);
    }
}
