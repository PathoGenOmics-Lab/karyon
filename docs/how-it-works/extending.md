# Writing a track

The [thirty track types](../tracks.md) the crate ships are implementations
of one trait with no privileged access to the figure. `CoverageTrack` and
`PileupTrack` see exactly what a track written outside the crate sees, so a
track type the crate does not have is not a fork of it, it is one more file in
the project that needs it. This page is that trait: what the figure hands a
track, what it does on the track's behalf, and what it leaves to the track.

## The entry test

Before anything else, one question decides whether the thing belongs here at
all.

!!! warning "Does `draw` read `ctx.scale`?"
    If it does not, the x axis is a sample list, a category or a count, and the
    plot is a bar chart, a line chart or a heatmap that happened to be handed
    genomic data. A general plotting library already draws those better.

That is the rule the crate holds itself to, and three tracks were removed under
it rather than kept for the sake of a longer list.

Five of the thirty do not read `ctx.scale`, and each of them answers for it
in its own module doc: `IdeogramTrack` draws the whole sequence so that it can
say where the region is, `TreeTrack` and `TanglegramTrack` measure evolutionary
distance across and share the other axis instead, `SnpTrack` lays out its own
columns because throwing the invariant ones away is the point of the panel, and
`LegendTrack` carries no coordinates at all. The list, with the reason for each,
is on the [tracks page](../tracks.md#the-entry-test). An exception costs a
paragraph of explanation, which is about the right price.

## The trait

```rust
pub trait Track {
    fn height(&self, scale: &Scale) -> f64;
    fn label(&self) -> Option<&str> { None }
    fn y_axis_width(&self, theme: &Theme) -> f64 { 0.0 }
    fn draw(&self, ctx: &mut DrawContext<'_>);
}
```

Two methods are required and two have defaults that are right most of the time.
A track knows how tall it wants to be and how to draw itself, and it never
decides where it sits: the `Figure` stacks the bands, reserves the gutter and
the axis strip, hands each track the rectangle it may paint in, and clips it
there.

## A whole track

Here is a complete one. It draws a tick per position, which is what anything
whose only quantity is "here" needs: restriction sites, primer landing sites,
integration sites, the positions a peak caller kept.

```rust
use karyon::{DrawContext, Scale, Track};

/// A tick per site, for data whose only quantity is "here".
pub struct RugTrack {
    sites: Vec<u64>,
    label: Option<String>,
    height: f64,
}

impl RugTrack {
    /// A rug over 0-based positions.
    pub fn new(sites: impl Into<Vec<u64>>) -> Self {
        RugTrack { sites: sites.into(), label: None, height: 16.0 }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }
}

impl Track for RugTrack {
    fn height(&self, _scale: &Scale) -> f64 {
        self.height
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        // Half a stroke inside the band: a rule centred on the clip edge keeps
        // only half its ink and renders at half the weight it asked for.
        let floor = band.bottom() - 0.5;
        ctx.svg.line(band.x, floor, band.right(), floor, &ctx.theme.rule, 1.0);

        for site in &self.sites {
            // Only what is on screen, so the cost follows the width of the
            // figure and not the length of the list.
            if !ctx.region.contains(*site) {
                continue;
            }
            let x = ctx.scale.x_center(*site);
            ctx.svg.line(x, band.y, x, floor, &ctx.theme.accent, 1.0);
        }
    }
}
```

The whole of `draw` is a dozen lines. It goes into a figure the same way
anything else does, from either layer:

```rust
let svg = plot("chr1:1-2000")?
    .add_coverage(vec![40.0; 2000])
    .label("depth")
    .add_track(RugTrack::new(vec![119, 480, 1_755]).label("cut sites"))
    .to_svg();
```

`Figure::push` takes it too, and `Figure::push_boxed` takes a `Box<dyn Track>`
when the stack is built at runtime and the arms have different types.

## What `DrawContext` holds

```rust
pub struct DrawContext<'a> {
    pub svg: &'a mut SvgWriter,
    pub scale: &'a Scale,
    pub theme: &'a Theme,
    pub band: Rect,
    pub axis: Rect,
    pub region: &'a Region,
}
```

`svg`
: Where to write. A track never builds SVG strings itself, which is what keeps
  escaping and number formatting in one place.

`scale`
: The shared horizontal mapping, identical for every track in the figure. This
  is the thing that makes it a genomic track.

`theme`
: Shared colours, fonts and the corner radius. Reach for `theme.accent` for a
  single series, `theme.color(i)` for a categorical one, `theme.rule` for
  baselines, `theme.muted` for secondary text, and `theme.surface()` rather than
  `theme.background` when blending a tint towards the page, because the
  background is allowed to be `"none"` and blending against that gives back
  something that is not a colour.

`band`
: The rectangle this track may paint in. `Rect` carries `x`, `y`, `w`, `h` and
  the three derived edges `right()`, `bottom()` and `mid_y()`. `y` grows
  downwards.

`axis`
: The strip immediately to the left of the band, as wide as this track asked for
  in `y_axis_width`. Zero width unless it asked.

`region`
: The region on display, 0-based half-open. `region.contains(pos)` is the cheap
  way to skip data that is off screen.

## The clip edge

The figure opens one clip per track, covering the axis strip and the band
together, and closes it after `draw` returns. Nothing a track draws can reach a
neighbour, so a stem that overshoots or a label that overruns is a cosmetic
problem rather than a corrupted figure.

The clip has one consequence worth knowing, because it bites the same way every
time. A stroke is centred on its path, so a one pixel rule drawn exactly on
`band.bottom()` keeps only the half of its width that falls inside the band, and
renders at half the weight of an identical rule drawn anywhere else. The tracks
in the crate that put a rule on a band boundary all set it half a stroke inside:

```rust
let floor = band.bottom() - 0.5;
ctx.svg.line(band.x, floor, band.right(), floor, &ctx.theme.rule, 1.0);
```

Text has the same problem from the other side. A label whose baseline sits on
`band.bottom()` has the tail of its `g` and `y` cut off, so lift it by a
descender's worth.

## Height, and why it takes a scale

Most tracks return a constant. Some cannot, because how much room they need
depends on what is on screen: `FeatureTrack` packs features into rows and
collisions are measured in pixels, so zooming out merges two rows into one and
zooming in splits them again.

A track whose height follows the view reads the view off the scale, since
`height` is called before `draw` and is not handed the region:

```rust
fn height(&self, scale: &Scale) -> f64 {
    let (start, end) = scale.bounds();   // the region, 0-based half-open
    // ... work out how many rows the data in that window needs
}
```

`height` runs during layout, once per render, so it may do real work:
`FeatureTrack` runs its whole packing there and again in `draw`. It has to be a
pure function of the track and the scale, though. Returning a different number
the second time would mean the bands no longer match the heights the figure laid
the stack out with.

## Asking for a value axis

A track that returns more than zero from `y_axis_width` is given
`DrawContext::axis`, a strip between the track labels and the plotting area, and
may draw its own ticks and numbers there. The figure reserves the **widest**
request across every track and gives that width to all of them, so the plotting
areas still start at the same x and the tracks still line up.

Ask for exactly the width of the widest label the track will print, measured
with `text_width`, plus a little padding. Then check `ctx.axis.w > 0.0` before
drawing into it, because a caller may have turned the axis off:

```rust
use karyon::svg::text_width;
use karyon::{Anchor, DrawContext, Scale, Theme, Track};

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
        let here = self.sites.iter().filter(|s| ctx.region.contains(**s)).count();
        for site in self.sites.iter().filter(|s| ctx.region.contains(**s)) {
            let x = ctx.scale.x_center(*site);
            ctx.svg.line(x, band.y, x, band.bottom(), &ctx.theme.accent, 1.0);
        }
        if ctx.axis.w > 0.0 {
            ctx.svg.text(
                ctx.axis.right() - 4.0,
                band.mid_y() + ctx.theme.font_size * 0.35,
                &here.to_string(),
                &ctx.theme.muted,
                ctx.theme.font_size - 1.0,
                Anchor::End,
            );
        }
    }
}
```

The default of zero is the right answer for most tracks. A pileup or a sequence
has no value to put a number on, and a strip asked for and left empty narrows
every plotting area in the figure.

## `SvgWriter`

The writer has exactly the elements the tracks need and no more. That is the
whole reason the crate has no dependencies and the whole reason the output is
plain SVG 1.1 that opens unchanged in a browser, in Inkscape and in Illustrator.

| Call | Draws |
|:-----|:------|
| `rect`, `rect_opacity`, `rect_rounded`, `rect_outline` | rectangles, filled or stroked |
| `circle`, `circle_ringed` | a dot, and a dot with a page-coloured ring that keeps overlapping marks legible |
| `line`, `polyline`, `polygon` | strokes and filled shapes |
| `path`, `path_stroked` | a ready-made `d` attribute, filled or stroked |
| `text`, `text_bold`, `text_rotated` | labels. `y` is the text baseline, not its centre or its top |
| `glyph` | one symbol stretched to an exact box with `textLength`, which is what a sequence logo is made of |
| `begin_clip`, `begin_clip_path`, `end_group` | a clipped group, if a track needs one of its own inside its band |

Three behaviours are worth relying on. Degenerate shapes are skipped rather than
written, so a zero-width rectangle, a negative radius or a coordinate that came
out `NaN` produces no element instead of a broken one. Text content is escaped,
so a feature name containing `&` or `<` cannot break the document. And numbers
are written with at most three decimals and no trailing zeros, which is the
single biggest lever on the size of a genome-wide figure.

Three helpers live alongside it. `karyon::svg::text_width` gives the advance
width of a string using Helvetica's own metrics, which are exact for the default
font stack, so a label that is measured before it is drawn actually fits.
`karyon::svg::num` is the same coordinate formatter the writer uses, for a track
that assembles a `d` attribute of its own. And `Anchor` picks whether the `x`
you pass is the left edge, the middle or the right edge of the text.

## Things to hold to

- **Draw only what is on screen.** The cost of a figure should follow its width,
  not the size of the input. `region.contains`, `scale.bounds()` and
  `scale.pos_at_x` are the three tools for that, and
  [Scale](scale.md) is the longer version of the argument.
- **Bin above one base per pixel.** Data denser than the output has to say what
  a pixel column means rather than overdrawing it. `Aggregate` is one answer, an
  envelope of the extremes is another, and both are honest in a way that "the
  last value drawn wins" is not.
- **Be deterministic.** No hash iteration order, no clock, no unstable sort
  without a tie-break. The same input has to render byte-identical output, or
  the figure in the paper and the figure in the repository stop being the same
  figure.
- **Coordinates in, pixels out.** Take 0-based half-open positions in the
  constructors, like the rest of the crate, and let the scale do the arithmetic.
  [Coordinates](coordinates.md) is the whole of that convention.

## Where a track is not the answer

`Rings` draws a circular sequence, and its rings implement the parallel trait,
`Ring`: `thickness`, an optional `gap`, and a `draw` handed a `RingContext` that
carries a `Polar` and two radii instead of a `Scale` and a `Rect`. A ring maps
position to an angle, so it is not a `Track` and cannot be pushed onto a
`Figure`. What the two have in common is `Panels`, the sheet that holds figures
of either kind, and `Drawing` is exactly the small surface `Panels` needs from
both: state your dimensions, and render yourself with these ids.

## Next

- [Tracks](../tracks.md), for the thirty already written to this trait.
- [Scale](scale.md), for the binning a new track has to do for itself.
- [Contributing](../about/contributing.md), for what a track shipped by the
  crate needs on top of the trait.
