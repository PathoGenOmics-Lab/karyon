# Scale

Every track in a figure draws through one `Scale`, and asks it how much
sequence a pixel covers before it draws anything. This page shows that map, how
each track reduces what falls under one pixel column, and when that choice
changes what a figure says.
{ .k-lead }

<figure class="k-plate" markdown>
![An association scan across a draft assembly of twelve contigs laid end to end on one axis: association points banded by contig, with a peak on the fourth contig crossing the significance line, a depth profile that dips at every contig join, and the contig names underneath](../assets/figures/example-genomewide.svg){ width="940" height="288" loading="lazy" }
</figure>

A four-megabase genome drawn at the default width of 900 pixels leaves about
800 pixel columns for the tracks, so each column covers about 5,000 bases. A
track that drew one mark per base would write four million marks to fill 800
columns: a file no viewer opens, and a picture no reader can tell from the
right one. So every track draws at the resolution of the output rather than
the resolution of the input, and a per-base depth profile over four megabases
comes out as about 800 points.

## The one map

```rust
use karyon::{Region, Scale};

let region = Region::new("NC_000962.3", 0, 4_000_000)?;
let scale = Scale::new(&region, 12.0, 872.0);

assert_eq!(scale.bounds(), (0, 4_000_000));
assert!((scale.bp_per_px() - 4_587.2).abs() < 0.05);
```

`Scale` is a linear map from positions to x coordinates in the image. The
figure builds one from the [region](coordinates.md) and the plotting area, and
hands the same one to every track. That sharing is the whole reason the tracks
line up: no track computes its own x.

| Method | Answers |
|:--|:--|
| `x(pos)` | the left edge of the base at 0-based `pos` |
| `x_center(pos)` | the middle of that base |
| `x_at(pos)` | the left edge of a fractional position |
| `pos_at_x(x)` | the inverse: which fractional position lies under a pixel |
| `bp_per_px()` | positions covered by one pixel; above 1, the data has to be binned |
| `px_per_bp()` | pixels covered by one position; above roughly 8, there is room for a letter |
| `bounds()` | the region on display, 0-based and half-open |
| `x0()`, `width()` | where the plotting area starts, and how wide it is |

- **It does not clamp.** A position outside the region maps outside the
  plotting area, and that is the useful answer: a gene that begins before the
  window is drawn as the whole rectangle it is, and the clip over the track's
  band decides how much of it shows. Clamping would stop the rectangle at the
  border and say the gene ends there.
- **Thresholds are in pixels per position.** Every decision that changes with
  zoom is a comparison against `px_per_bp()` or `bp_per_px()`, so a threshold
  means the same thing at any image width, and a track finds out it has room
  for letters by asking rather than by being told the figure got wider.
- **It has no opinion about the unit.** Bases, alignment columns and samples of
  a raw signal are all position `n` to it, which is why the same track code
  works on all three (see [the region as a coordinate
  system](coordinates.md#the-region-is-a-coordinate-system)).

`bounds()` exists because a track is asked how tall it is before it is asked
to draw, and some tracks cannot answer without knowing what is on screen: a
pileup packs only the reads in view, so its height follows the view.

## One point per pixel column

[`CoverageTrack`](../tracks/signal-sequence.md#coveragetrack) holds one value
per base, which is the shape `samtools depth` output arrives in. It does not
walk that array. It walks the pixel columns of its band, asks `pos_at_x` which
bases fall under each one, reduces them to a single number and emits one
point:

```text
for column in 0 .. band width:
    x     = band.x + column
    first = floor(scale.pos_at_x(x))
    last  = ceil(scale.pos_at_x(x + 1))
    point = reduce(values[first .. last])    # by the track's Aggregate
```

The size of the output therefore follows the width of the figure, not the
length of the input. Four megabases of per-base depth, a four-million-base
sequence track and a ruler render to about 21 KB of SVG. A test fails if that
document reaches 100 KB, so a change that starts emitting a mark per base
breaks the build rather than quietly getting slower.

## Aggregate, and when it changes the answer

Reducing a column to one number is a choice, and `Aggregate` makes it explicit
because no default is right for every question.

| `Aggregate` | A column reads | Use it when |
|:--|:--|:--|
| `Max`, the default | the highest value under it | hunting for spikes and duplications |
| `Mean` | the mean of the values under it | reading the overall level |
| `Min` | the lowest value under it | hunting for dropouts and deletions |

The choice stops being cosmetic exactly when the feature you are looking for is
narrower than a pixel. Take a 100-base dropout to zero in an otherwise flat
profile at 60x, drawn across four megabases at the default width, where one
column holds about 5,000 bases:

| `Aggregate` | The column holding the dropout reads | The dropout |
|:--|:--|:--|
| `Max` | 60 | invisible |
| `Mean` | 58.8 | a 2% dip, inside the noise of a real profile |
| `Min` | 0 | the full depth of the hole |

=== "Rust"

    ```rust
    use karyon::{plot, Aggregate};

    let svg = plot("NC_000962.3:1-4,000,000")?
        .add_coverage(depth)
        .label("depth")
        .adjust(|track| track.aggregate(Aggregate::Min))
        .to_svg();
    ```

=== "Command line"

    ```bash
    karyon NC_000962.3:1-4,000,000 \
      --coverage depth.bedgraph --label depth --aggregate min -o depth.svg
    ```

`Max` is the default because a spike is the commoner thing to hunt for, and
because it never shows a value the data does not contain. It is also the one
that quietly hides a deletion. The locus figure below uses `Aggregate::Min` for
that reason: at two and a half bases per pixel, the dropout is the thing worth
not smoothing away.

<figure class="k-plate" markdown>
![The rpoB locus over two kilobases: a depth profile whose dropout reaches zero because each pixel column takes its minimum, a reference track too zoomed out for bases that says so, the rpoB gene and its resistance-determining region, and variants coloured by consequence](../assets/figures/example.svg){ width="900" height="304" loading="lazy" }
</figure>

Two smaller rules go with it:

- **A value that is not finite is ignored**, so one `NaN` does not erase its
  column.
- **A column with no finite value under it is missing, not zero.** No point is
  emitted for it, so missing data is never drawn as zero depth. In the `Bars`
  style the column stays empty; in `Area` and `Line` the outline runs straight
  on to the next column that has data.

A base that a bedGraph or a `samtools depth` file leaves out is a different
case. `CoverageTrack::from_spans`, like `--coverage`, starts from zero over the
whole region, because zero is what a depth file means by leaving a base out.

## Past single-base resolution

[`SequenceTrack`](../tracks/signal-sequence.md#sequencetrack) climbs three
levels of detail, the way a genome browser does:

| Pixels per base | What it draws |
|:--|:--|
| 7 or more (`letter_threshold`) | the letter, in the base's colour |
| 0.6 or more (`block_threshold`) | a coloured block per base, overdrawn by a quarter pixel so antialiasing leaves no hairline gaps |
| fewer | the words `zoom in to see bases`, and nothing else |

The third level is not decoration. A million bases drawn as blocks is a million
rectangles; with the hint instead, a figure holding a million-base sequence
track is well under 4 KB, which a test checks. Above the floor the cost follows
the figure as well, because only the bases inside the region are visited: a
whole reference sequence can be handed to a track that shows two hundred of
its bases.

Two other tracks climb a ladder of their own:

- **[`CodonTrack`](../tracks/scales-keys.md#codontrack)** draws translated
  letters while a codon is at least 0.95 of the font size wide (11.4 pixels in
  the default theme) and a sequence is attached. Narrower than that it draws
  alternating cells while a codon is at least 2.5 pixels wide, so the partition
  stays visible, and below that a plain bar over the coding sequence, so the
  codon numbers underneath still have something to belong to.
- **[`DynseqTrack`](../tracks/signal-sequence.md#dynseqtrack)** draws each base
  as its letter from 5 pixels per base (`letter_threshold`) and as a bar from 1
  pixel per base. Below that it draws an envelope of the extremes in one
  neutral ink, never in a base colour, since a column spanning forty bases has
  no base.

??? info "Why `DynseqTrack` takes no `Aggregate`"
    Every answer `Aggregate` offers is wrong for a signed score. A maximum
    hides a strong negative, a minimum hides a strong positive, and a mean
    cancels a +2 against a -2 into a zero that says the model ignored the
    place. So the track keeps both extremes of each column, or draws nothing.

## How the other tracks bin

The same discipline takes a different shape wherever the data are denser than
the output:

- **[`WindowTrack`](../tracks/signal-sequence.md#windowtrack)** keeps the
  lowest and the highest value in each pixel column, in the default
  `WindowStyle::Steps`, rather than one summary: a signed statistic that swings
  both ways inside a column has to show both swings. Runs of columns holding
  the same pair become one rectangle, since drawn a pixel at a time every join
  leaves a seam. `WindowStyle::Line` draws the mean of each column instead.
- **[`SquiggleTrack`](../tracks/reads-molecules.md#squiggletrack)** draws each
  column as the range of the samples under it while there are fewer than
  `point_threshold` pixels per sample (4 by default). That is the envelope an
  oscilloscope draws: the extremes are honest, and the shape between them is
  not there. From the threshold up the samples are drawn as a line, and from
  two and a half times the threshold each one gets a ringed dot.
- **[`MsaTrack`](../tracks/comparison.md#msatrack)** merges neighbouring cells
  of one colour in a row into one rectangle. Most of an alignment agrees with
  itself, so a row that agrees for eighty columns costs one rectangle rather
  than eighty.
- **[`PileupTrack`](../tracks/reads-molecules.md#pileuptrack)** looks for
  mismatches only once a base is at least a fifth of a pixel wide
  (`mismatch_threshold`, 0.2 by default). Below that, finding one would mean
  walking every base of every read to draw something nobody can see.
- **[`VariantTrack`](../tracks/variation.md#varianttrack)** draws at most one
  tick per pixel column for each category with `VariantStyle::Tick`
  (`--style tick`), and names none of them in a tooltip. Two hundred thousand
  calls over four megabases come to about 150 KB as ticks, and to about 49 MB
  as lollipops, which draw every call.
- **[`ManhattanTrack`](../tracks/variation.md#manhattantrack)** does not bin: it
  draws every point it is given, so a scan of a million tests is a document of
  a million marks.

## Height follows the zoom

`Track::height` is handed the `Scale` because some tracks cannot know how tall
they are without it. Collisions are measured in pixels, not in bases, so how
many rows a track needs is a property of the zoom.

- **[`FeatureTrack`](../tracks/annotation.md#featuretrack)** packs only the
  features in view: sorted by start, first fit from the left, with four pixels
  between neighbours on a row. A name too long to sit inside its feature is
  drawn to its right, so the width of the name is reserved while packing, or
  the next feature would sit on top of it. Zoom out, and two genes that had a
  row each collide and the track grows a row; zoom in, and it shrinks again.
- **`PileupTrack`** packs only the reads that overlap the view, so a window
  holding two reads gets a band two reads tall. It stops at forty rows by
  default and writes how many reads that hid on the figure, rather than
  dropping them quietly. `max_rows(None)`, or `--max-rows all`, lifts the cap.

<figure class="k-plate" markdown>
![A read pileup under a depth profile, a candidate call and the reference: reads coloured by strand, mismatches painted where reads disagree with the reference, two reads with a deletion, one with an insertion, and faded reads of low mapping quality](../assets/figures/example-pileup.svg){ width="920" height="472" loading="lazy" }
</figure>

## Where next

<div class="grid cards" markdown>

-   **[Coordinates](coordinates.md)**

    The convention every position arrives in, and what each file format's
    numbers become.

-   **[Writing a track](extending.md)**

    What the figure does on a track's behalf, and what it leaves to the track.

-   **[Track catalogue](../tracks/index.md)**

    Every track, with the options that move the thresholds on this page.

</div>
