//! The sequences of a [`Genome`], drawn as the axis they are.
//!
//! A ruler of global coordinates under a concatenated genome would be a ruler
//! of a coordinate system nothing else uses: nobody has ever quoted a position
//! as "1,437,902 bases into the assembly". So this track labels the sequences
//! instead, and marks where each one ends, which is the thing a reader of a
//! genome-wide figure actually needs to know.

use crate::genome::Genome;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::{DrawContext, Track};

/// The sequence bar of a genome-wide figure.
///
/// ```
/// use karyon::{Figure, Genome, GenomeTrack};
///
/// let genome = Genome::new([("chr1", 900_000u64), ("chr2", 600_000), ("chr3", 500_000)]);
/// let svg = Figure::new(genome.region())
///     .push(GenomeTrack::new(genome).label("assembly"))
///     .to_svg();
/// assert!(svg.contains("chr1"));
/// ```
#[derive(Debug, Clone)]
pub struct GenomeTrack {
    genome: Genome,
    label: Option<String>,
    height: f64,
    color: Option<String>,
    show_names: bool,
    min_name_width: f64,
}

impl GenomeTrack {
    /// A bar of the sequences in `genome`.
    pub fn new(genome: Genome) -> Self {
        GenomeTrack {
            genome,
            label: None,
            height: 22.0,
            color: None,
            show_names: true,
            min_name_width: 4.0,
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the band height in pixels.
    pub fn height(mut self, height: f64) -> Self {
        self.height = height.max(6.0);
        self
    }

    /// Sets the colour the alternating blocks are shaded from.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Draws or hides the sequence names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Sets how much room a name needs before it is written.
    ///
    /// An assembly of two hundred contigs has two hundred names and room for
    /// perhaps twelve. The ones that do not fit are left out rather than
    /// overprinted, and [`GenomeTrack::named`] says how many were written.
    pub fn min_name_width(mut self, pixels: f64) -> Self {
        self.min_name_width = pixels.max(0.0);
        self
    }

    /// The genome.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// How many sequence names fit at this scale, and how many did not.
    pub fn named(&self, scale: &Scale, theme: &Theme) -> (usize, usize) {
        if !self.show_names {
            return (0, self.genome.len());
        }
        let size = theme.font_size - 1.0;
        let mut written = 0usize;
        for (name, start, end) in self.genome.spans() {
            let width = scale.x(end) - scale.x(start);
            if width >= text_width(name, size) + self.min_name_width {
                written += 1;
            }
        }
        (written, self.genome.len() - written)
    }
}

impl Track for GenomeTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let hue = self
            .color
            .clone()
            .unwrap_or_else(|| ctx.theme.muted.clone());
        // Alternating rather than one colour per sequence: an assembly has more
        // contigs than any palette has colours, and what a reader needs from
        // the bar is where one ends and the next starts, which two shades say
        // as well as two hundred would.
        let shades = [
            mix(ctx.theme.surface(), &hue, 0.30),
            mix(ctx.theme.surface(), &hue, 0.16),
        ];
        let size = ctx.theme.font_size - 1.0;
        let bar = self.height.min(band.h);
        let top = band.y + (band.h - bar) / 2.0;

        for (index, (name, start, end)) in self.genome.spans().into_iter().enumerate() {
            let x0 = ctx.scale.x(start);
            let x1 = ctx.scale.x(end).max(x0 + 0.5);
            ctx.svg
                .rect(x0, top, x1 - x0, bar, &shades[index % shades.len()]);

            if self.show_names && x1 - x0 >= text_width(name, size) + self.min_name_width {
                ctx.svg.text(
                    (x0 + x1) / 2.0,
                    top + bar / 2.0 + size * 0.35,
                    name,
                    &ctx.theme.foreground,
                    size,
                    Anchor::Middle,
                );
            }
        }

        let (_, unnamed) = self.named(ctx.scale, ctx.theme);
        if unnamed > 0 && self.show_names {
            // A bar of unlabelled blocks looks like a bar of labelled ones
            // with the labels forgotten. Saying how many is the difference.
            // Inside the bar, because the band is what the track is clipped to
            // and a line above it is a line nobody sees.
            ctx.svg.text(
                band.right() - 4.0,
                top + bar / 2.0 + (size - 1.0) * 0.35,
                &format!("+{unnamed} unlabelled"),
                &ctx.theme.muted,
                size - 1.0,
                Anchor::End,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;

    fn genome() -> Genome {
        Genome::new([("chr1", 900_000u64), ("chr2", 600_000), ("chr3", 500_000)])
    }

    #[test]
    fn every_sequence_gets_a_block_and_they_alternate() {
        let svg = Figure::new(genome().region())
            .show_region_label(false)
            .push(GenomeTrack::new(genome()))
            .to_svg();
        let theme = Theme::light();
        let dark = mix(theme.surface(), &theme.muted, 0.30);
        let pale = mix(theme.surface(), &theme.muted, 0.16);
        assert_eq!(svg.matches(&dark).count(), 2, "chr1 and chr3");
        assert_eq!(svg.matches(&pale).count(), 1, "chr2");
    }

    #[test]
    fn names_are_written_where_they_fit() {
        let svg = Figure::new(genome().region())
            .show_region_label(false)
            .push(GenomeTrack::new(genome()))
            .to_svg();
        for name in ["chr1", "chr2", "chr3"] {
            assert!(svg.contains(&format!(">{name}</text>")), "missing {name}");
        }
    }

    #[test]
    fn contigs_too_narrow_to_name_are_counted_rather_than_overprinted() {
        // Two hundred contigs and room for a handful of names.
        let many = Genome::new((0..200).map(|index| (format!("contig_{index}"), 20_000u64)));
        let theme = Theme::light();
        let scale = Scale::new(&many.region(), 0.0, 800.0);
        let track = GenomeTrack::new(many);
        let (written, unnamed) = track.named(&scale, &theme);
        assert_eq!(written, 0, "four pixels each");
        assert_eq!(unnamed, 200);

        let svg = Figure::new(track.genome().region())
            .show_region_label(false)
            .push(GenomeTrack::new(track.genome().clone()))
            .to_svg();
        assert!(svg.contains("+200 unlabelled"), "a silent drop is a lie");
        // And it is inside the band, not on the line above it that the clip
        // would take away.
        let at = svg.find("+200 unlabelled").unwrap();
        let head = &svg[..at];
        let y: f64 = head[head.rfind("y=\"").unwrap() + 3..]
            .split('"')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let clip_top: f64 = svg[svg.find("<clipPath").unwrap()..]
            .split("y=\"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert!(y > clip_top, "the count sits above the band at {y}");
    }

    #[test]
    fn turning_names_off_does_not_leave_the_count_behind() {
        let svg = Figure::new(genome().region())
            .show_region_label(false)
            .push(GenomeTrack::new(genome()).show_names(false))
            .to_svg();
        assert!(!svg.contains("unlabelled"));
        assert!(!svg.contains(">chr1</text>"));
    }

    #[test]
    fn an_empty_genome_draws_nothing_and_does_not_panic() {
        let empty = Genome::new(Vec::<(String, u64)>::new());
        let svg = Figure::new(empty.region())
            .show_region_label(false)
            .push(GenomeTrack::new(empty).label("assembly"))
            .to_svg();
        assert!(svg.starts_with("<svg "));
        assert!(!svg.contains("NaN"));
    }
}
