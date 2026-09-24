# The karyon logo

The name, set with its k drawn in bands.

The k is cut into seven horizontal bands of one height, five across the x-height
and two up the ascender. That is the idea karyon is built on, rows that line up
across everything drawn on them, and it is also how a chromosome ideogram is
drawn. The middle band is the centromere. It is the only band where the stem,
the arm and the leg meet, it narrows where they join the stem the way a
chromosome narrows at its centromere, and it is the one drawn in the accent.

## Files

| File | Use |
|:-----|:----|
| `karyon.svg` | The logotype on a light page: this README, the home page of the site. |
| `karyon-dark.svg` | The logotype on a dark page. |
| `docs/assets/karyon_logo.svg`, `karyon_logo_dark.svg` | The same two, where the site can serve them. |
| `docs/assets/karyon_mark.svg`, `karyon_mark_dark.svg` | The banded k alone: the site header, and anywhere the name is already written beside it. |
| `docs/assets/favicon.svg` | The dark mark on a tile of the dark page, so it reads on a light or a dark tab bar. |

## Construction

Measured in units of the x-height, which is 100.

| | |
|:--|:--|
| Bands | 17.28 high with a 3.4 gap, so five fill the x-height exactly; the ascender is 141.36 |
| Stem | 35 wide, a little over Inter ExtraBold's 31, for the ink the gaps take out |
| Arm and leg | 41 wide horizontally, running 0.8 across for every unit up or down, mirror images about the middle of the x-height |
| Centromere | the arm and the leg start 2.5 inside the stem, which leaves the waist |

The rest of the name is Inter ExtraBold (SIL Open Font License), tracked by -2
units, with the ink traps at its joins closed so they do not read as chips at a
logo's size. Everything is outlines; no font has to be installed to see it.

## Colours

| | Light page | Dark page |
|:--|:--|:--|
| k | `#0072b2` | `#3d9be9` |
| Centromere | `#d55e00` | `#ef7d2b` |
| Letters | `#1b1f23` | `#e6edf3` |
| Tile (favicon) | | `#14181d` |

The blue and the vermillion are the site's primary and accent, from the
Okabe-Ito palette the crate draws with; the ink and the tile are the crate's own
light foreground and dark background.

Keep the bands, and keep the centromere the only band in the accent. Do not
recolour the logotype band by band, stretch it, or set the name in another face
beside the banded k.
