# Writing a track

Add a track type the crate does not have. The thirty-six that ship are
implementations of one small trait with no privileged access to the figure, so
yours is one more file in your project, not a fork of karyon.
{ .k-lead }

## The entry test

Before writing anything, one question decides whether the thing is a track at
all.

!!! warning "Does `draw` read `ctx.scale`?"
    A track lives on the figure's shared coordinate axis, so its `draw` places
    marks through `ctx.scale`, the one mapping from position to pixel that
    every band is drawn with. If it does not, its x axis is a sample list, a
    category or a count, and the plot is a bar chart, a line chart or a
    heatmap that was handed genomic data. A general plotting library draws
    those better.

The axis is usually genomic position, but any shared integer coordinate
counts: an alignment column, a raw-signal sample, or the time pivots that
[PhylodynamicTrack](../tracks/evolution-surveillance.md#phylodynamictrack) and
[SurveillanceTrack](../tracks/evolution-surveillance.md#surveillancetrack)
share. A plot whose x is a count of genomes or a list of sample names fails,
which is why karyon draws no accumulation curve and no clustered frequency
heatmap. A rarefaction over a presence matrix is a statistic to compute, not a
plot type, and it needs no `Track`.

Thirty-one of the thirty-six shipped tracks draw through `ctx.scale`. The five
that do not each answer for it in their own module documentation:

- [IdeogramTrack](../tracks/whole-genome.md#ideogramtrack) draws the whole
  sequence across the plotting area on purpose. A track that showed only the
  region on display could not say where the region is: it would be a picture of
  the window, drawn inside the window.
- [TreeTrack](../tracks/phylogeny.md#treetrack) and
  [TanglegramTrack](../tracks/phylogeny.md#tanglegramtrack) measure
  evolutionary distance across. What they share with their neighbours is the
  other axis, because a leaf is a row.
- [SnpTrack](../tracks/variation.md#snptrack) lays out its own columns. Its x
  is a site index, since throwing the invariant columns away is the point of
  the panel, and no shared ruler survives that.
- [LegendTrack](../tracks/scales-keys.md#legendtrack) carries no coordinates
  at all. It is a band that has to be stacked, sized and clipped like every
  other, which is exactly what a track is.

An exception costs a paragraph of explanation, which is about the right price.

### And a thing that is not a track

Metadata columns fail the entry test in the other direction. A sample's
lineage is not at a base, and no zoom level brings more of it into view. Drawn
as a track it would need an x nobody has, and the first pan would slide a
sample's lineage off the end of that sample's own row.

So `Traits` is not a track. It attaches through `traits` to the six tracks
drawn as a row per named thing,
[MatrixTrack](../tracks/variation.md#matrixtrack),
[MsaTrack](../tracks/comparison.md#msatrack),
[SnpTrack](../tracks/variation.md#snptrack),
[CladeTrack](../tracks/phylogeny.md#cladetrack),
[DomainTrack](../tracks/comparison.md#domaintrack) and
[LocusTrack](../tracks/comparison.md#locustrack), and is drawn in the strip
those tracks already reserve left of the plotting area, beside the row names
and the dendrogram. It survives every pan and zoom untouched, because nothing
in it was placed at a coordinate. A phylogeny draws the same columns
[beside its tips](../guide/phylogenetics.md#draw-time-branches-and-sample-traits-together).

## The trait

```rust
pub trait Track {
    fn height(&self, scale: &Scale) -> f64;
    fn label(&self) -> Option<&str> { None }
    fn y_axis_width(&self, _theme: &Theme) -> f64 { 0.0 }
    fn on_coordinates(&self) -> bool { true }
    fn draw(&self, ctx: &mut DrawContext<'_>);
}
```

Two methods are required and three have defaults that are right most of the
time. A track says how tall it wants to be and how to draw itself, and never
decides where it sits: the `Figure` stacks the bands in push order, reserves
the label gutter and the value axis strip, hands each track the rectangle it
may paint in, and clips it there.

| Method | Answers | Default |
|:--|:--|:--|
| `height(&Scale)` | how tall the band is, in pixels before the figure's scale factor | required |
| `label()` | the name in the left gutter, which the figure draws | `None` |
| `y_axis_width(&Theme)` | how much room this track wants for a [value axis](#asking-for-a-value-axis) | `0.0` |
| `on_coordinates()` | whether the track is laid on the shared coordinate axis | `true` |
| `draw(&mut DrawContext)` | the marks | required |

Return `false` from `on_coordinates` only when the track's x is not the
figure's coordinate at all, as a phylogeny's is not. `plot()` and the command
line leave the automatic ruler off a figure where nothing is on the
coordinates, since it would measure nothing.

## A whole track

A tick per site, for data whose only quantity is "here": restriction sites,
primer landing sites, integration sites, the positions a peak caller kept.
Here is all of it, as a program you can run:

```rust
use karyon::{plot, DrawContext, Scale, Track};

/// A tick per site, for data whose only quantity is "here".
pub struct RugTrack {
    sites: Vec<u64>,
    label: Option<String>,
}

impl RugTrack {
    /// A rug over 0-based positions.
    pub fn new(sites: impl Into<Vec<u64>>) -> Self {
        RugTrack { sites: sites.into(), label: None }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }
}

impl Track for RugTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        16.0
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let stroke = ctx.px(1.0);
        // Half a stroke inside the band, or the clip takes half the rule.
        let floor = band.bottom() - stroke / 2.0;
        ctx.svg.line(band.x, floor, band.right(), floor, &ctx.theme.rule, stroke);

        for &site in &self.sites {
            // Only what is on screen: the cost follows the width of the
            // figure, not the length of the list.
            if !ctx.region.contains(site) {
                continue;
            }
            let x = ctx.scale.x_center(site);
            // The tooltip names the site the way a reader counts, from 1.
            ctx.svg.begin_titled(&format!("site {}", site + 1));
            ctx.svg.line(x, band.y, x, floor, &ctx.theme.accent, stroke);
            ctx.svg.end_group();
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    plot("chr1:1-2000")?
        .add_coverage(vec![40.0; 2000])
        .label("depth")
        .add_track(RugTrack::new(vec![119, 480, 1_755]).label("cut sites"))
        .save("rug.svg")?;
    Ok(())
}
```

The whole of `draw` is a dozen lines. `plot()` takes the track through
`add_track`, `Figure::push` takes it as it is, and `Figure::push_boxed` or
`Plot::add_boxed` take a `Box<dyn Track>` when the stack is built at runtime
and its arms have different types. `Plot::label` and `Plot::adjust` reach only
the shipped tracks, so name yours before handing it over, as above.

## What `DrawContext` holds

```rust
pub struct DrawContext<'a> {
    pub svg: &'a mut SvgWriter,
    pub scale: &'a Scale,
    pub theme: &'a Theme,
    pub band: Rect,
    pub axis: Rect,
    pub region: &'a Region,
    pub visual_scale: f64,
}
```

`svg`
: Where to write. A track never builds SVG text itself, which keeps escaping
  and number formatting in one place.

`scale`
: The shared horizontal mapping, identical for every track in the figure:
  `x(pos)` is the left edge of a base, `x_center(pos)` its middle,
  `pos_at_x(x)` the way back, `bp_per_px()` and `px_per_bp()` the zoom, and
  `bounds()` the region on display. This is what makes a track a track.

`theme`
: Shared colours, fonts and measurements, already scaled for the figure.
  `theme.accent` for a single series, `theme.color(i)` for a categorical one,
  `theme.rule` for baselines, `theme.muted` for secondary text, and
  `theme.tokens` for stroke widths and marker sizes. Blend a tint towards
  `theme.surface()`, not `theme.background`: the background may be `"none"` on
  a transparent page, and nothing can be blended with that.

`band`
: The rectangle this track may paint in: `x`, `y`, `w`, `h` and the derived
  `right()`, `bottom()` and `mid_y()`. `y` grows downwards.

`axis`
: The strip immediately left of the band, as wide as this track asked for in
  `y_axis_width`, and zero wide unless it asked.

`region`
: The region on display, 0-based and half-open. `region.contains(pos)` is the
  cheap way to skip what is off screen.

`visual_scale`
: The figure's factor for track geometry, from its profile and its density.
  The band already includes it; `ctx.px(value)` applies it to a measurement of
  your own, such as a stroke or a gap. Genomic x still comes only from `scale`.

## The clip edge

The figure opens one clip per track, over its band and the axis strip it asked
for, and closes it when `draw` returns. A track that asked for no axis is
clipped to its band alone. Nothing a track draws can reach a neighbour, which
makes overhang free: `scale` does not clamp, so a gene that starts before the
window is drawn whole with its left edge off the page, and the clip decides how
much shows. The track's label is the one thing drawn outside, by the figure.

The clip has one consequence that bites the same way every time. A stroke is
centred on its path, so a one pixel rule drawn exactly on `band.bottom()` keeps
only the half inside the band and renders at half the weight of the same rule
anywhere else. Shipped tracks that rule a band's edge set the rule half a
stroke inside, as the rug above does with `band.bottom() - stroke / 2.0`.

Text has the same problem from the other side. A label whose baseline sits on
`band.bottom()` loses the tails of its g and y, so lift it by a descender.

## Height, and why it takes a scale

Most tracks return a constant. Some cannot, because the room they need depends
on what is on screen: `FeatureTrack` packs features into rows by collisions
measured in pixels, so features that share a row zoomed in need a row each
zoomed out, and `PileupTrack` packs only the reads in view. `height` is called
before `draw` and is not handed the region, so it reads the view off the
scale:

```rust
impl RugTrack {
    /// How deep the ticks stack where sorted sites fall within three pixels
    /// of each other, in the region on display.
    fn tiers(&self, scale: &Scale) -> usize {
        let (start, end) = scale.bounds();
        let (mut tiers, mut run, mut last) = (1, 0, f64::NEG_INFINITY);
        for &site in self.sites.iter().filter(|&&s| s >= start && s < end) {
            let x = scale.x_center(site);
            run = if x - last < 3.0 { run + 1 } else { 1 };
            tiers = tiers.max(run);
            last = x;
        }
        tiers
    }
}
```

`height` then returns `self.tiers(scale) as f64 * 8.0`, and `draw` asks
`self.tiers(ctx.scale)` the same question. Keep that arithmetic in one method:
`height` runs during layout and may run more than once, since asking a figure
for its dimensions lays it out again, so it has to be a pure function of the
track and the scale. A different answer the second time would hand `draw` a
band that no longer fits what it draws.

The number `height` returns is in pixels before the figure's scale factor,
which multiplies it by `ctx.visual_scale`. Inside `draw`, measure the band or
use `ctx.px(8.0)` rather than repeating the bare constant.

## Asking for a value axis

A track that returns more than zero from `y_axis_width` gets `ctx.axis`, a
strip between the labels and the plotting area, to draw its own ticks and
numbers in. The figure reserves the widest request across every track for all
of them, so the plotting areas still start at the same x.

```rust
use karyon::svg::text_width;
use karyon::{Anchor, DrawContext, Scale, Theme, Track};

/// A rug that prints how many of its sites are in view.
pub struct CountedRug {
    sites: Vec<u64>,
}

impl Track for CountedRug {
    fn height(&self, _scale: &Scale) -> f64 {
        16.0
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        text_width(&self.sites.len().to_string(), theme.font_size - 1.0) + 8.0
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        let here: Vec<u64> = self
            .sites
            .iter()
            .copied()
            .filter(|site| ctx.region.contains(*site))
            .collect();
        for &site in &here {
            let x = ctx.scale.x_center(site);
            ctx.svg.line(x, band.y, x, band.bottom(), &ctx.theme.accent, ctx.px(1.0));
        }
        if ctx.axis.w > 0.0 {
            ctx.svg.text(
                ctx.axis.right() - 4.0,
                band.mid_y() + ctx.theme.font_size * 0.35,
                &here.len().to_string(),
                &ctx.theme.muted,
                ctx.theme.font_size - 1.0,
                Anchor::End,
            );
        }
    }
}
```

Ask for the width of the widest label the track will print, measured with
`text_width` at the size it will be drawn, plus a little padding. Draw into the
strip only when `ctx.axis.w > 0.0`: a request that can come out at zero, from a
setting that hides the scale or nothing to label, gets a strip of no width.

Zero is right for most tracks. A pileup or a sequence has no value to put a
number on, and a strip asked for and left empty narrows every plotting area in
the figure.

## `SvgWriter`

The writer has exactly the elements the tracks need and no more, which is the
whole reason the crate has no dependencies and writes plain SVG 1.1 that opens
unchanged in a browser, in Inkscape and in Illustrator.

| Calls | Draw |
|:--|:--|
| `rect`, `rect_opacity`, `rect_rounded`, `rect_rounded_opacity`, `rect_rounded_edged`, `rect_outline` | rectangles: filled, translucent, rounded, edged or outlined |
| `circle`, `symbol` | a dot, or a categorical `Symbol`: circle, square, diamond or triangle |
| `circle_ringed`, `symbol_ringed` | the same inside a ring of page colour, which keeps overlapping marks apart |
| `line`, `polyline`, `polygon`, `polygon_edged` | strokes and filled shapes |
| `line_pattern`, `polyline_pattern`, `path_stroked_pattern` | dashed or dotted strokes, from a `LinePattern` |
| `path`, `path_stroked` | a ready-made `d` attribute, filled or stroked |
| `text`, `text_bold`, `text_haloed`, `text_rotated` | labels; `y` is the baseline, not the centre or the top |
| `glyph` | one symbol stretched to an exact box, which is what a sequence logo is made of |
| `begin_titled`, `begin_titled_inert` | a group with a tooltip; the inert one lets the pointer through to the marks beneath |
| `begin_clip`, `begin_clip_path`, `end_group` | a clipped group of your own inside the band; `end_group` closes any group |

Three behaviours are worth relying on. Degenerate shapes are skipped rather
than written, so a zero-width rectangle, a negative radius or a coordinate that
came out NaN produces no element instead of a broken one. Text is escaped, so a
feature name containing `&` or `<` cannot break the document. Numbers are
written with at most three decimals and no trailing zeros, the biggest single
lever on the size of a genome-wide figure.

A few helpers sit beside it. `karyon::svg::text_width` measures a string as
the wider of Inter and of the Helvetica widths Liberation Sans and Arial share,
character by character, so a label measured before it is drawn fits in
whichever face the reader's machine sets it; `text_width_strong` does the same
for semibold text and `mono_width` for the monospaced stack.
`karyon::svg::num` is the formatter the writer uses, for a track that builds a
`d` attribute of its own, and `karyon::fit_text` shortens text with an
ellipsis to fit a width. `Anchor` says whether the `x` you pass is the start,
the middle or the end of the text, and `TextStyle` sets a face, a weight or
extra spacing on one run of text with `text_styled`.

## Pitfalls

- **Draw only what is on screen.** The cost of a figure should follow its
  width, not the size of the input. `region.contains`, `scale.bounds()` and
  `scale.pos_at_x` are the tools, and [Scale](scale.md) is the longer
  argument.
- **Bin above one base per pixel.** When `scale.bp_per_px()` is above one, say
  what a pixel column means rather than overdrawing it. `Aggregate` is one
  answer and an envelope of the extremes another; both are honest in a way
  that "the last value drawn wins" is not.
- **Close every group you open.** The figure ends the track's clip with the
  first `end_group` after `draw` returns. A group left open takes that call
  instead, the clip stays open, and every track after it is drawn inside this
  one's clip and cut away.
- **Put the exact value in a tooltip.** A colour or a height only
  approximates a number; wrap each mark in `begin_titled` with the value
  itself, as the shipped tracks do.
- **Never write a colour down.** Take ink from `ctx.theme` and derive shades
  with `karyon::theme::mix`, `wash` and `contrast_ink`, so the track is right
  on a light page and on a dark one.
- **Scale your own measurements.** Size strokes, gaps and marks with `ctx.px`,
  so a presentation profile or a compact density moves them with everything
  else.
- **Be deterministic.** No hash map iteration order, no clock, no unstable sort
  without a tie-break. The same input has to give byte-identical output, or the
  figure in the paper and the figure in the repository stop being the same
  figure.
- **Coordinates in, pixels out.** Take 0-based, half-open positions in the
  constructors, like the rest of the crate, and let the scale do the
  arithmetic. [Coordinates](coordinates.md) is the whole convention.

## Where a track is not the answer

A circular sequence maps position to an angle, so it is not a stack of bands.
`Rings` draws one, and its rings implement the parallel trait `Ring`:
`thickness`, an optional `gap` of 5 pixels by default, and a `draw` handed a
`RingContext`, which carries a `Polar` and the ring's `inner` and `outer` radii
where a track gets a `Scale` and a `Rect`. A ring cannot be pushed onto a
`Figure`. What the two share is `Panels`, the sheet that holds figures, rings
and maps side by side through the small `Drawing` trait: state your dimensions,
and render yourself with the ids you are given.

## Where next

<div class="grid cards" markdown>

-   **[Track catalogue](../tracks/index.md)**

    The thirty-six tracks already written to this trait.

-   **[Scale](scale.md)**

    The binning and zoom thresholds a new track works out for itself.

-   **[Coordinates](coordinates.md)**

    0-based, half-open positions in; 1-based labels where a reader looks.

-   **[Contributing](../about/contributing.md)**

    What a track shipped in the crate needs on top of the trait.

</div>
