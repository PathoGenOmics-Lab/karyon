# PDF

A figure is drawn once, as SVG, and a PDF is that SVG read back and written as
the PDF operators that draw the same thing. This page says what that keeps,
what it changes, and how to check a PDF before it goes to a journal.
{ .k-lead }

```bash
karyon rpoB reads.bam genes.gff3 calls.vcf.gz -o rpoB.pdf
```

```rust
use karyon::plot;

let figure = plot("NC_000962.3:761,000-763,000")?
    .add_coverage(depth)
    .label("depth")
    .into_figure();
figure.save_pdf("rpoB.pdf")?;
let pdf = figure.to_pdf(); // pdf.bytes, and pdf.notes
```

`Figure`, `Panels`, `Rings`, `Map` and `PhyloMap` each have `to_pdf()` and
`save_pdf(path)` beside `to_svg()` and `save_svg(path)`, and `Plot::save` writes
PDF when the name it is given ends in `.pdf`. `Pdf::from_svg` converts any SVG
written in the same terms, such as one a `Drawing` of your own returns.

## Read from the SVG, not drawn twice

Every track draws through one writer, `SvgWriter`, and the writer already
speaks SVG: a track hands it path data as a `d` string, and a fade comes back
from it as a `url(#id)` to fill with. A second backend behind a trait would
still have to read path data, colours and references, which is most of the
work, and it would break every track written outside the crate. Reading the SVG
keeps one drawing and one source of truth.

What holds the two together is a test, not care. Every method of the writer is
called once and the result converted, and an element or attribute the
converter does not read comes back as a note that fails the test by name. A
second test reads the writer's own source for its methods, so a new one fails
until that first test calls it.

## Size

A pixel is three quarters of a point, which is CSS's 96 pixels to the inch and
the size Inkscape and `rsvg-convert` give the same SVG. A figure 900 pixels wide
is a page 675 points, 9.4 inches, wide; ask for `--width` in pixels as you would
for the SVG. Inside the file every coordinate is the SVG's own, under one matrix
that scales and flips the page.

A PDF page can be at most 14,400 points, 200 inches, a side. A figure wider than
19,200 pixels is written smaller with a `/UserUnit` saying how much larger a
unit is, which every reader since Acrobat 7 honours; one that does not shows the
whole figure smaller rather than part of it.

## Text

The SVG asks for Inter and JetBrains Mono, and the PDF sets its text in three
faces it carries itself, each embedded in the file as a subset of the letters
the figure draws:

| The SVG's text | The PDF's |
|:--|:--|
| the default stack, and any other that is not monospaced | Liberation Sans |
| weight 600 and up, and `bold` | Liberation Sans Bold |
| the monospaced stack: coordinates, sequence, accessions | Liberation Mono |

Greek, the arrows and the relations karyon writes, ω, →, ≤, ≥ and ≈, come from
the same face as the letters around them. A serif stack, and a bold monospaced
one, ask for a face the PDF does not carry; karyon writes neither, and an SVG of
your own that does is set in the nearest of the three, with a note.

Each face is written as a composite font over one TrueType font holding only the
glyphs the page uses, with their widths and a map from each glyph back to its
character, so every reader draws the same letters, the text stays text, to be
searched, copied and read aloud, and a preflight check that wants every font
embedded finds them all, as [Checking a PDF](#checking-a-pdf) shows. The faces
are the Liberation fonts 2.00.1, which the crate carries in `fonts/` under the
SIL Open Font License; a PDF made with them is yours, under whatever terms you
choose. A subset is 4 to 30 kilobytes, so a figure's PDF is about that much
larger for each face it uses. A character none of the faces has, such as a
sample name in Chinese or Arabic, is drawn as a question mark, and the command
line says which characters on standard error, as `Pdf::notes` does in Rust. The
document's title and description keep every character either way.

Where a label sits was settled before the PDF existed, by the widths karyon
measures text with, and Liberation is drawn to the widths of Helvetica and
Arial, and Liberation Mono to Courier's. The measure is never narrower than
Liberation Sans, nor than Liberation Sans Bold for bold text, for printable
ASCII and for the characters beyond it that karyon writes itself: the ellipsis,
the en dash, the superscript two, the multiplication sign and the middle dot,
and ω, ≤, ≥, ≈ and →. For the monospaced stack it is Liberation Mono's 600
thousandths of an em exactly. A label that fitted in the SVG therefore fits in
the PDF, and every label is anchored with the widths of the face that draws it,
so one anchored at its end still ends where the SVG ended it. A character
karyon does not write, such as an accented letter in a sample name, is measured
at six tenths of an em, and Liberation Sans draws an accented capital wider than
that, an Ö at 778 thousandths, so a label of them can still run a little past
its room at the end it is not anchored by.

## Fades

A coverage profile is filled with a fade, strong under the line and faint at
the baseline. In the PDF a fade is a column of colour stretched over the box of
the shape it fills, with its opacity in a mask alongside, and clipped to the
shape: an image, which every reader tried draws in full. Opacity along a
shading needs a soft-mask group instead, and that was cut short by Quartz,
which macOS Preview draws with, and left out by Inkscape's import.

## What does not carry over

- **Hover titles.** A group's `<title>` is what a browser shows under the
  pointer; a page has no pointer, so they are left out. The figure's own title
  and description are kept, as the PDF's title and subject and as the
  alternative text of the one tagged figure on the page, which is where a
  screen reader looks.
- **The faces.** Liberation Sans and Liberation Mono, embedded, not Inter and
  JetBrains Mono, as above.
- **Size on disk.** The file is not compressed. Most figures come out about
  the size of their SVG, with the subsets of the faces on top, and a figure of
  many dots up to four times it, since each circle is four curves: the
  [association scan](../your-data/scan.md) is 3.9 times its SVG.

## No PNG

There is no PNG writer, and a name ending in `.png` is refused. A raster needs a
rasterizer, which turns glyph outlines into pixels and takes a few thousand
lines of anti-aliasing, stroking and compositing, for something every PDF and
SVG tool already does with hinted text:

```bash
pdftoppm -png -singlefile -r 300 rpoB.pdf rpoB
rsvg-convert -d 300 -p 300 -o rpoB.png rpoB.svg
```

## Checking a PDF

```bash
pdffonts rpoB.pdf
gs -o /dev/null -sDEVICE=nullpage rpoB.pdf
pdftotext rpoB.pdf -
```

`pdffonts` lists each face the figure uses, Liberation Sans, Liberation Sans
Bold and Liberation Mono, as a `CID TrueType` font with a subset tag before its
name, such as `AECRZT+LiberationSans`, and `yes` under `emb`, `sub` and `uni`:
embedded, a subset, and mapped back to the text. That is what a journal's
upload check, and the font check of the PDF/X and PDF/A standards that printers
and archives hold a file to, want of every font. Ghostscript marks anything it
had to repair with a line of `****`, and has no such line for a sound file;
leave out its `-q`, which hides those lines along with the rest. `pdftotext`
gives the labels back, Greek and arrows included.
