# The karyon logo

A chromosome, one region of it opened up, and the name drawn in that region
with its coverage underneath.

Everything is drawn on one grid of square cells, the way a figure is drawn on
one coordinate axis. At the top is a chromosome ideogram with its G-bands, its
centromere and a region boxed in pink. Dotted lines open that region out across
the name, whose letters are cut into horizontal bands like the chromosome's.
Under the name runs its depth, with a dropout under the r, a variant call in the
dropout, and a ruler. Read from top to bottom it is a karyon figure: one region,
tracks stacked on it, and every row lined up.

## Files

| File | Use |
|:-----|:----|
| `karyon.svg` | The logo on a light page: this README and the home page of the site. |
| `karyon-dark.svg` | The logo on a dark page. The chromosome, the zoom lines and the ruler are lighter; the name and the coverage do not change. |
| `docs/assets/karyon_logo.svg`, `karyon_logo_dark.svg` | The same two, where the site can serve them. |
| `docs/assets/karyon_mark.svg` | The k over its own coverage, with the dropout and the call: the site header and the drawer. The same file serves both schemes. |
| `docs/assets/favicon.svg` | The mark on a night tile, so it reads on a light or a dark tab bar. |

## Construction

| | |
|:--|:--|
| Letters | Inter ExtraBold (SIL Open Font License) at an x-height of 19 cells: a cell is inked when the letter covers at least half of it, after the letter is shifted by whichever fraction of a cell leaves the fewest cells half covered, so its stems come out whole cells wide |
| Bands | 3 cells of ink and 1 of gap, so five bands fill the x-height; the ascender of the k and the descender of the y are stretched to end on a band edge, 27 cells up and 8 down |
| Spacing | 2 cells between letters, 1 fewer from r to y and from y to o |
| Chromosome | 53 cells long and 7 high, with rounded ends, a pinch at the centromere and G-bands in four strengths |
| Coverage | one column under every column of the name, 1 cell deep under the r |
| Size | the logo is 117 by 79 cells with a margin of 2, the mark 26 by 35, the favicon 40 by 40 |

Each element is one outline: the name, the coverage, the call, the ruler, the
zoom, and the chromosome with its bands painted over it. Where the colour
changes from row to row, in the name and the coverage, the outline is filled
with a gradient that steps at the row edges instead of being cut into a shape
per row. Where two elements touch, as the zoom touches the box and, in the
mark, the call touches its floor, the one painted first runs on under the
other. Two shapes of different colours that meet edge to edge with nothing
under them show a hairline in a PDF viewer or on paper, and no two shapes here
do, so the files can go into a PDF, a poster or a vector editor as they are.
There is no font and no filter, so nothing has to be installed to see them.

Show them at a whole number of screen pixels per cell. Scaled by a fraction, a
pixel drawing loses some of the one-pixel gaps between its bands, which is why
the site draws the logo 351 or 234 pixels wide and the mark 35 high, and this
README 351 wide.

## Colours

| | |
|:--|:--|
| Name | every row takes its colour from its height, from `#2b3be8` at the top of the k through `#6b3fe0` and `#b04ed8` to `#ff4fa3` at the foot of the y |
| Coverage | the same ramp, from indigo at its base to pink at the top of its highest column |
| Region, centromere, call | `#ff4fa3` |
| Chromosome, light page | `#ddd6fb`, banded in `#b9a8f5`, `#8f6fe8`, `#6b3fe0` and `#3b1f8f` |
| Chromosome, dark page | `#d9d2ff`, banded in `#b9a8f5`, `#8f6fe8`, `#b04ed8` and `#6b3fe0` |
| Zoom lines and ruler | `#b04ed8` and `#3b1f8f` on a light page, `#d9d2ff` on a dark one |
| Tile (favicon) | `#120b2b` |

These colours belong to the logo and to nothing else. The figures and the site
keep the Okabe-Ito palette the crate draws with, which is chosen to stay legible
with a colour vision deficiency.

Keep the grid, the bands, and the colour running from the top of the k to the
foot of the y. Do not smooth the logo, stretch it, or show it at a fraction of a
pixel per cell.
