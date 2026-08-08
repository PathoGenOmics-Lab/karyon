# Theming

A `Theme` is the colours, fonts and spacing every track in a
[figure](figure.md) shares. It is a plain struct with every field public, so
there is no builder to learn and nothing hidden: start from `Theme::light` or
`Theme::dark` and overwrite what you need. This page is every field, why the two
palettes hold the colours they hold, and the two methods and three functions that
blend a colour against the page.

```rust
use karyon::{Figure, Region, Theme};

let figure = Figure::new(Region::new("chr1", 0, 1_000)?).theme(Theme::dark());
```

`Theme` derives `Debug`, `Clone` and `PartialEq`, and its `Default` is
`Theme::light`. `Plot::theme` and `Rings::theme` take one too; `Panels::theme`
takes one for the sheet's own title and letters, while each panel keeps the
theme it was rendered with.

## The two themes

=== "Light"

    ![A coverage profile with a dropout, the reference sequence, two gene models and variants coloured by consequence, all over one coordinate axis](../assets/figures/example.svg)

=== "Dark"

    ![The same locus figure on a dark page: the coverage profile, reference sequence, gene models and variants redrawn in the dark palette's steps against a near-black background](../assets/figures/example-dark.svg)

The dark theme is **a selected set of colours, not an inversion** of the light
one. Flipping a light palette is the obvious approach and it fails: a dark
background wants a narrower lightness band than a light one, so half of a
flipped palette lands outside that band and two of the entries stop being
distinguishable. The dark steps were chosen against the dark surface and
validated there. The two palettes are the same length and share no entries.

## Every field

| Field | Type | What it is |
|:------|:-----|:-----------|
| `background` | `String` | The page. Set it to `"none"` for a transparent SVG. |
| `foreground` | `String` | Titles, primary axes and tick marks. |
| `muted` | `String` | Secondary text: coordinates, track labels, legends. |
| `rule` | `String` | Quiet baselines, guides and secondary separators. |
| `accent` | `String` | Default colour for a track that was not given one. |
| `palette` | `Vec<String>` | Categorical colours, cycled by `Theme::color`. |
| `bases` | `BaseColors` | Per-nucleotide colours, used by every track that paints a base: the reference sequence, a pileup mismatch, an alignment column, a logo letter. |
| `insertion` | `String` | An insertion mark in a read pileup, deliberately off the categorical palette so it never impersonates a series. |
| `corner_radius` | `f64` | Corner radius of a data mark, in pixels. |
| `font_family` | `String` | Font stack written on the root `<svg>` element. |
| `font_size` | `f64` | Tick labels and in-plot annotations. |
| `label_font_size` | `f64` | Track labels in the left gutter. |
| `title_font_size` | `f64` | The figure title. |
| `cap_height_ratio` | `f64` | Height of a capital letter as a fraction of the font size. |
| `tokens` | `VisualTokens` | Shared stroke, marker, tick, feature, legend and spacing measurements. |

And the values the two constructors give them:

| Field | `Theme::light()` | `Theme::dark()` |
|:------|:-----------------|:----------------|
| `background` | `#ffffff` | `#14181d` |
| `foreground` | `#1b1f23` | `#e6edf3` |
| `muted` | `#4b5563` | `#aab4c0` |
| `rule` | `#d7dce2` | `#3a424c` |
| `accent` | `#0072b2` | `#3987e5` |
| `insertion` | `#8e44ad` | `#8e44ad` |
| `corner_radius` | `2.5` | `2.5` |
| `font_family` | `Liberation Sans, Arial, Helvetica, sans-serif` | `Liberation Sans, Arial, Helvetica, sans-serif` |
| `font_size` | `12.0` | `12.0` |
| `label_font_size` | `12.0` | `12.0` |
| `title_font_size` | `18.0` | `18.0` |
| `cap_height_ratio` | `0.72` | `0.72` |

`corner_radius` is the one field that is neither a colour nor a size of text.
Rounded ends are most of the difference between a figure that looks drawn and
one that looks emitted. Set it to zero for square corners. Features, reads in a
pileup, legend swatches, variable-site cells and open reading frames all use it,
and the radius is clamped to half the smaller side, so a short bar becomes a
lozenge rather than losing its geometry to an oversized corner.

`Theme::scaled(factor)` returns the same palette and font stack with all three
font sizes, the corner radius and every visual token scaled together.
`Figure`, `Plot`, `Rings` and `Panels` expose `visual_scale`; named
`RenderProfile` values provide a coordinated starting point. See the
[visual-system guide](visual-system.md).

## The categorical palette

The default is six hues from the Okabe-Ito palette, ordered and trimmed so that
**every pair** stays apart under colour vision deficiency, not just neighbouring
ones.

| Slot | Light | Dark |
|:-----|:------|:-----|
| 0 | `#0072b2` | `#3987e5` |
| 1 | `#d55e00` | `#d95926` |
| 2 | `#009e73` | `#199e70` |
| 3 | `#cc79a7` | `#c98500` |
| 4 | `#e69f00` | `#d55181` |
| 5 | `#7b3294` | `#9085e9` |

A palette that survives colour vision deficiency is the default because two
categories a reader cannot tell apart are one category. Both palettes here were
run through a colour vision validator rather than chosen by eye, and the numbers
decided the outcome. Every pair clears a separation of 6.7 in OKLab hundredths,
inside the band that is sound when identity is also carried by something other
than colour, and in this crate it always is: a legend, a letter, or a label.

Six is where the measurement stopped, not where the eye got bored. A seventh hue
could not be added without some pair collapsing: an olive against the vermillion
came out at 1.8 under protanopia, which is indistinguishable. The palette was
eight slots once, and two of those were a brown and a grey 1.8 apart under
protanopia, so any figure with seven or eight categories was drawing two of them
the same colour. Cutting it to six is what fixed that.

`Theme::color(index)` reads a slot and wraps at the end, so a track with more
categories than the palette has reuses colours rather than panicking. When the
palette is emptied it falls back to `accent`.

```rust
use karyon::Theme;

let theme = Theme::light();
assert_eq!(theme.color(0), theme.accent);
assert_eq!(theme.color(6), theme.color(0));
```

`accent` is slot zero of its own palette in both themes, and that is not a
coincidence. `accent` is what a single unlabelled series gets; if it were not
slot zero, a one-series figure and a two-series figure would disagree about what
the first colour is.

## Nucleotides

`BaseColors` holds five strings, one per nucleotide plus `other` for `N` and any
other IUPAC or unexpected symbol. `BaseColors::of(base)` reads one, case
insensitively, with `U` treated as `T`.

| Base | `conventional()` | `colorblind_safe()` |
|:-----|:-----------------|:--------------------|
| A | `#33a02c` | `#009e73` |
| C | `#1f78b4` | `#0072b2` |
| G | `#e08214` | `#e69f00` |
| T, U | `#e31a1c` | `#d55e00` |
| other | `#9e9e9e` | `#9e9e9e` |

Nucleotides are the honest exception to the rule above. `conventional()` is the
default, because a figure that recolours the bases surprises every reader who
has ever opened a genome browser. It is also not safe: measured pairwise,
**adenine and guanine sit 1.7 apart under protanopia**, on a scale where 8 is
the floor for telling two colours apart, and that is the transition pair, the
commonest substitution there is. So the convention costs a real reader a real
distinction.

`colorblind_safe()` fixes it, with a closest pair of 11.0 under deuteranopia.
The cost is that green now means adenine rather than the thymine red a reader
expects, so label it or say so in the caption.

```rust
use karyon::{BaseColors, Theme};

let mut theme = Theme::light();
theme.bases = BaseColors::colorblind_safe();
```

## Blending against the page

One method and three free functions do the colour arithmetic the tracks need.
`wash` is re-exported at the crate root; it, `mix` and `contrast_ink` all live
in `karyon::theme`.

**`Theme::surface()`** is a real colour standing in for the page. `background` is
allowed to be `"none"`, which is a transparent page rather than a colour, and
blending against something that is not a colour gives back something that is not
a colour: the mark comes out `fill="none"` and disappears. A quiet cell in a
matrix, the bar under a base that agrees with the reference and the zebra tint
on a panel of variable sites are all tints towards the page, so they ask for
`surface()` instead of `background`. A transparent page is whatever it ends up
composited onto, which cannot be known from inside the theme, so the theme's own
ink is the one clue it has: dark ink implies a light page, light ink a dark one.

**`mix(a, b, t)`** blends two `#rrggbb` colours, `t` running from all of `a` at
zero to all of `b` at one. Anything that is not a six-digit hex string comes
back as `a` unchanged. What it is for: a greyscale that belongs to its theme.
Cytogenetic bands are specified as shades from white to black, and hardcoding
those shades makes a dark figure look like a light one somebody forgot to
invert. Mixing the theme's own background and foreground gives the same ladder
in whichever direction the page runs.

**`wash(color, theme)`** is a pale version of a colour against the page,
`mix(color, theme.surface(), 0.62)`. It is what a large filled shape wears. A
gene arrow or a matrix cell at full saturation shouts, and eight of them make a
colour chart rather than a figure; the hue belongs in the edge and a wash of it
in the body. Small marks are the other way round and want the colour undiluted,
which is why this is a function and not a rule.

**`contrast_ink(color)`** returns `#1b1f23` or `#ffffff`, whichever stays
readable on top of `color`, using Rec. 601 luma. Feature labels are drawn inside
coloured boxes whose colour the caller chooses, so the ink has to be picked per
box rather than fixed by the theme. Anything that is not a `#rrggbb` string is
assumed dark and gets white.

```rust
use karyon::theme::{contrast_ink, mix};
use karyon::{wash, Theme};

let theme = Theme::dark();
let body = wash(theme.color(0), &theme);            // the fill of a large shape
let halfway = mix(theme.surface(), &theme.foreground, 0.5);
let ink = contrast_ink(theme.color(0));             // text on top of that colour
```

## Strand

`strand_color(strand, theme)` is where the strand convention lives: the reverse
strand gets slot one of the palette, and everything else, forward and unknown,
gets slot zero. Every track that paints by strand goes through it, so replacing
the palette moves all of them together.

```rust
use karyon::{strand_color, Strand, Theme};

let theme = Theme::light();
assert_eq!(strand_color(Strand::Forward, &theme), theme.color(0));
assert_eq!(strand_color(Strand::Reverse, &theme), theme.color(1));
assert_eq!(strand_color(Strand::Unknown, &theme), theme.color(0));
```

## Fonts and text measurement

`font_family` is written as an attribute on the root `<svg>` element and nothing
else. No font is embedded and no font is fetched, which is part of why the
output is a file that opens unchanged in a browser, in Inkscape and in
Illustrator. It also means the stack names fonts the viewer has to have; a
machine without them substitutes, and text will set slightly differently there.

That matters more than it sounds, because the crate has to know how wide a
string will be before it draws it. `karyon::svg::text_width(text, font_size)`
answers that, and the widths it uses are Arial-compatible advance widths.
Liberation Sans and Arial are the first two fonts in the default stack and are
metrically compatible. For
the default theme the measurement is exact rather than approximate. One flat
width per character would under-reserve for a run of capitals by about a fifth,
which is precisely what a column of sample accessions is, and a label that
overruns the space reserved for it gets clipped. Characters outside printable
ASCII fall back to a wide default, so an accented name reserves a little too
much rather than too little.

```rust
use karyon::svg::text_width;
use karyon::Theme;

let theme = Theme::light();
let room = text_width("ERR1234567", theme.label_font_size);
```

If you change `font_family` to a stack with different metrics, `text_width` will
disagree in the third significant figure. That is usually harmless, and it is
worth checking on the longest label in the figure.

`cap_height_ratio` is the one font setting that is not a size. A sequence logo
stretches each letter to an exact box, which means working back from a box
height to a font size, and that division needs the cap height. The default of
`0.72` suits Liberation Sans and Arial. **Change it and `font_family` together**, or
logo letters will sit slightly proud of their boxes.

## A transparent background

```rust
use karyon::Theme;

let mut theme = Theme::light();
theme.background = "none".into();
```

The page rectangle is skipped entirely when `background` is exactly the string
`"none"`, so the SVG composites onto whatever is behind it. That is what a
figure meant to sit on a coloured slide, inside another document, or on a page
that has a dark mode wants, since a white rectangle behind it would show.

```rust
use karyon::{Feature, FeatureTrack, Figure, Region, Theme};

let mut theme = Theme::light();
theme.background = "none".into();

let bands = vec![
    Feature::new(0, 96).color("#0072b2"),
    Feature::new(0, 62).color("#d55e00"),
    Feature::new(0, 80).color("#009e73"),
];

Figure::new(Region::new("mark", 0, 100)?)
    .width(120.0)
    .theme(theme)
    .label_width(0.0)
    .show_region_label(false)
    .push(FeatureTrack::new(bands).row_height(18.0).show_names(false))
    .save_svg("mark.svg")?;
```

Note that `foreground` and `muted` still have to be readable against whatever
the page turns out to be. `surface()` covers the tints, as described above, but
it can only guess from the ink, so a transparent figure meant for a dark page
wants `Theme::dark()` as its starting point with the background removed.

## Building a theme of your own

Start from a constructor and overwrite. The struct has fourteen fields and it
has gained fields between releases, so a struct literal is code that will break
on an upgrade for no reason.

```rust
use karyon::Theme;

let mut theme = Theme::light();
theme.font_family = "Inter, Helvetica, Arial, sans-serif".into();
theme.font_size = 12.0;
theme.label_font_size = 12.0;
theme.corner_radius = 0.0;
theme.palette = vec!["#4c72b0".into(), "#dd8452".into(), "#55a868".into()];
theme.accent = theme.palette[0].clone();
```

Four things to keep in mind when replacing a palette:

1. **Set `accent` to the first slot.** Everything that draws a single series
   without being told a colour uses `accent`, and everything that draws several
   uses `color(0)` upwards. If they disagree, a figure with one series and a
   figure with two disagree about their first colour.
2. **Every entry has to be `#rrggbb`.** `mix`, `wash` and `contrast_ink` all
   parse six-digit hex and nothing else, and they degrade quietly rather than
   failing: `mix` hands back its first argument and `contrast_ink` assumes a
   dark box. A named CSS colour will draw, but anything blended from it will
   not be what you meant.
3. **Check the pairs, not the neighbours.** A palette is used for categories
   that appear in any order, so the pair a reader has to separate is rarely the
   pair sitting next to each other in the list.
4. **A shorter palette is fine.** `color` wraps, so three entries cycle through
   three colours. Whether that is legible depends on whether a legend or a label
   is also carrying the identity.

!!! warning "Adding a field is a breaking change"
    `Theme` has gained fields before, `corner_radius`, `insertion` and
    `cap_height_ratio` among them. Code that builds one with
    `Theme { .. }` has to be updated each time; code that starts from
    `Theme::light()` never does.

## Next

- [Figure](figure.md), for where a theme is attached and what else the figure
  decides.
- [Writing a track](../how-it-works/extending.md), for how a track reads one.
- [Recipes](../recipes.md), for the same stack rendered light and dark.
