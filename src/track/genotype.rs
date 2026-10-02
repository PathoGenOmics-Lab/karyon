//! The call of each sample at each site, a row per sample, on the coordinates.
//!
//! A cohort's VCF carries, after its sites, one column per sample saying what
//! that sample was called there. This track draws it the way a genome browser
//! does: one row per sample and one cell per record, each at its own position
//! on the shared axis, so a column of cells stands under the lollipop a
//! [`VariantTrack`](crate::VariantTrack) draws for the same record and under
//! the gene it falls in.
//!
//! # What a cell says
//!
//! A call is a list of copies, one per chromosome set, each naming an allele:
//! `1` for a haploid genome, `0/1` for a diploid one, `0/0/0/1` for a
//! tetraploid. What every ploidy has in common is the share of its copies that
//! are not the reference, so that is what the colour is. A haploid `1` and a
//! diploid `1/1` are both all alternate and drawn at full strength, `0/1` at
//! half, `0/0/0/1` at a quarter. A multi-allelic `1/2` carries no copy of the
//! reference, so it is all alternate too, and its tooltip says which two.
//!
//! Four things are kept apart, since each is a different statement: a call of
//! the reference, which is a short quiet bar as an agreement is in a
//! [`SnpTrack`](crate::SnpTrack); a call carrying an alternate allele, a full
//! cell in the hue; a sample with no call, a pale full cell; and a stretch
//! with no record in it, which is the page. A call with any copy unknown, as
//! `./1` is, is no call: half a genotype read as a whole one would draw a
//! heterozygote as homozygous.
//!
//! # Many sites under one pixel
//!
//! A cell is drawn at a floor width, so at the zoom a cohort is read at most
//! cells would sit on top of each other. Twenty thousand sites of two hundred
//! samples over a megabase, drawn a cell at a time by
//! [`MatrixTrack`](crate::MatrixTrack), were an SVG of 258 MB and four million
//! rectangles. So the track works out which cells would overlap before it
//! draws any. A site whose cell touches no other is drawn as a cell, centred
//! on its base. Sites whose cells overlap form a cluster, and a cluster is
//! painted a pixel at a time, as a dense [`SnpTrack`](crate::SnpTrack) is, but
//! on the coordinates: each pixel takes the share of alternate copies among the
//! calls under it, in eight steps, and any alternate copy at all is at least
//! the first step, so one heterozygote among forty references is not drawn as
//! a reference. The same twenty thousand sites are 0.94 MB this way, in
//! 14,509 rectangles, forty rows of them, since the rows stop at forty and
//! count the rest.

use std::ops::Range;

use crate::region::Region;
use crate::scale::Scale;
use crate::svg::{text_width, Anchor};
use crate::theme::{mix, Theme};
use crate::track::axis::group_thousands;
use crate::track::traits::Traits;
use crate::track::tree::{draw_tree, leaf_order, tree_beside_rows, TreeShape, TreeStyle};
use crate::track::{DrawContext, Rect, Track};
use crate::tree::Tree;

/// The steps a cell's share of alternate copies is drawn in, from the first
/// alternate copy to all of them.
const LEVELS: usize = 8;

/// How tall the bar of a reference call is, as a share of the row.
const REFERENCE_BAR: f64 = 0.36;

/// An allele index that names no allele: a copy written `.`.
const NO_ALLELE: u16 = u16::MAX;

/// The longest allele a tooltip spells out before it says how much more
/// there is: a structural call's REF can run to thousands of bases.
const SPELLED: usize = 12;

/// What one sample was called at one site.
///
/// Eight bytes and no heap, since a cohort is many of them: fifty thousand
/// sites of a thousand samples is fifty million calls. It keeps how many
/// copies the call has, how many of them are known, how many are not the
/// reference, whether it is phased, and the alleles of its first two copies,
/// which is every allele of a haploid or a diploid call. A call of more copies
/// keeps the count of its alternate copies whole, which is what its colour is
/// drawn from, and the alleles of its first two.
///
/// An allele is an index into the site's alleles: 0 is the reference and `i`
/// the `i`th alternate, as a VCF writes it. Any copy whose index is not 0 is
/// an alternate copy, whichever allele it names: `*`, the base a deletion
/// upstream took away, and a placeholder such as `<NON_REF>`, which says the
/// base is not the reference without saying what it is, are both not the
/// reference.
///
/// ```
/// use karyon::{Genotype, GenotypeState};
///
/// assert_eq!(Genotype::parse("0/1").unwrap().state(), GenotypeState::Heterozygous);
/// assert_eq!(Genotype::parse("1").unwrap().share(), Some(1.0));
/// assert_eq!(Genotype::parse("1/1").unwrap().share(), Some(1.0));
/// assert_eq!(Genotype::parse("0/0/0/1").unwrap().share(), Some(0.25));
/// // Half a genotype is not a genotype.
/// assert_eq!(Genotype::parse("./1").unwrap().state(), GenotypeState::NotCalled);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Genotype {
    copies: u8,
    called: u8,
    alternate: u8,
    phased: bool,
    alleles: [u16; 2],
}

/// Which of the four things a call is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenotypeState {
    /// Every copy is the reference.
    Reference,
    /// Some copies are the reference and some are not.
    Heterozygous,
    /// No copy is the reference, which `1/2` is as much as `1/1`.
    Alternate,
    /// At least one copy is unknown, or there is no call at all.
    NotCalled,
}

/// Why a `GT` value did not read, which the reader turns into a sentence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Unread {
    /// A copy that is neither a number nor `.`, or no copy at all.
    Malformed,
    /// More copies than a call keeps count of.
    TooMany(usize),
}

impl Genotype {
    /// The most copies a call is read with.
    pub const MAX_COPIES: usize = u8::MAX as usize;

    /// Reads a `GT` value as a VCF writes it: allele indices or `.`, joined by
    /// `/` where the call is unphased and `|` where it is phased, with the
    /// leading `/` or `|` VCF 4.4 allows. `None` for anything else, or for a
    /// call of more than [`Genotype::MAX_COPIES`] copies.
    pub fn parse(gt: &str) -> Option<Genotype> {
        Genotype::read(gt).ok().map(|(call, _)| call)
    }

    /// The call and the largest allele index it names, which the reader holds
    /// against the alleles its row has.
    pub(crate) fn read(gt: &str) -> Result<(Genotype, u32), Unread> {
        let joint = |c: char| c == '/' || c == '|';
        // VCF 4.4 may write the phase of the first copy in front of it, as
        // `|0|1`, which is the same call as `0|1`.
        let body = gt.strip_prefix(joint).unwrap_or(gt);
        let marks: Vec<char> = gt.chars().filter(|c| joint(*c)).collect();
        let copies = body.split(joint).count();
        if copies > Genotype::MAX_COPIES {
            return Err(Unread::TooMany(copies));
        }
        let mut call = Genotype::not_called();
        let mut largest = 0u32;
        for (copy, token) in body.split(joint).enumerate() {
            let allele = match token {
                "." => NO_ALLELE,
                digits if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) => {
                    let index: u32 = digits.parse().unwrap_or(u32::MAX);
                    largest = largest.max(index);
                    // An index this large is past the alleles of any row that
                    // can be written, and is refused there; kept below the
                    // value that means no allele so it is not read as one.
                    index.min(u32::from(NO_ALLELE) - 1) as u16
                }
                _ => return Err(Unread::Malformed),
            };
            if allele != NO_ALLELE {
                call.called += 1;
                call.alternate += u8::from(allele != 0);
            }
            if let Some(kept) = call.alleles.get_mut(copy) {
                *kept = allele;
            }
        }
        call.copies = copies as u8;
        call.phased = !marks.is_empty() && marks.iter().all(|mark| *mark == '|');
        Ok((call, largest))
    }

    /// No call at all: a sample with nothing written for it.
    pub fn not_called() -> Genotype {
        Genotype {
            copies: 0,
            called: 0,
            alternate: 0,
            phased: false,
            alleles: [NO_ALLELE; 2],
        }
    }

    /// A haploid call of `allele`, where `u16::MAX` is a copy nobody called.
    pub fn haploid(allele: u16) -> Genotype {
        Genotype::from_copies(&[allele])
    }

    /// An unphased diploid call of `a` and `b`, where `u16::MAX` is a copy
    /// nobody called.
    pub fn diploid(a: u16, b: u16) -> Genotype {
        Genotype::from_copies(&[a, b])
    }

    fn from_copies(alleles: &[u16]) -> Genotype {
        let mut call = Genotype::not_called();
        call.copies = alleles.len() as u8;
        for (copy, allele) in alleles.iter().enumerate() {
            if *allele != NO_ALLELE {
                call.called += 1;
                call.alternate += u8::from(*allele != 0);
            }
            call.alleles[copy] = *allele;
        }
        call
    }

    /// How many copies the call has: 1 for a haploid call, 2 for a diploid
    /// one, and 0 for no call at all.
    pub fn copies(&self) -> usize {
        usize::from(self.copies)
    }

    /// How many of its copies name an allele.
    pub fn called(&self) -> usize {
        usize::from(self.called)
    }

    /// How many of its copies name an allele other than the reference.
    pub fn alternate(&self) -> usize {
        usize::from(self.alternate)
    }

    /// Whether every copy names an allele.
    pub fn is_called(&self) -> bool {
        self.copies > 0 && self.called == self.copies
    }

    /// The share of its copies that are not the reference, which is what the
    /// track draws: the one reading that means the same thing for a haploid,
    /// a diploid and a polyploid call. `None` unless every copy is known.
    pub fn share(&self) -> Option<f64> {
        self.is_called()
            .then(|| f64::from(self.alternate) / f64::from(self.copies))
    }

    /// Which of the four things this call is.
    pub fn state(&self) -> GenotypeState {
        if !self.is_called() {
            GenotypeState::NotCalled
        } else if self.alternate == 0 {
            GenotypeState::Reference
        } else if self.alternate == self.copies {
            GenotypeState::Alternate
        } else {
            GenotypeState::Heterozygous
        }
    }

    /// Whether the call was written phased, with `|`.
    pub fn phased(&self) -> bool {
        self.phased
    }

    /// The allele index of copy 0 or copy 1, or `None` for a copy that is not
    /// there, is not known, or is past the two a call keeps.
    pub fn allele(&self, copy: usize) -> Option<u16> {
        if copy >= self.copies().min(2) {
            return None;
        }
        Some(self.alleles[copy]).filter(|allele| *allele != NO_ALLELE)
    }

    /// The step of the eight a call is drawn in: 0 for the reference, and for
    /// any alternate copy at least 1. `None` for no call.
    fn level(&self) -> Option<usize> {
        self.is_called()
            .then(|| level(self.alternate(), self.copies()))
    }
}

/// The step `alternate` copies of `copies` are drawn in, any alternate copy
/// at least the first.
fn level(alternate: usize, copies: usize) -> usize {
    if copies == 0 || alternate == 0 {
        return 0;
    }
    let step = (alternate as f64 / copies as f64 * LEVELS as f64).round() as usize;
    step.clamp(1, LEVELS)
}

/// One record of a VCF, and what each sample was called there.
///
/// The site is a column of the figure and its calls are indexed by row, as a
/// [`SnpSite`](crate::SnpSite)'s alleles are.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeSite {
    /// Where the record is anchored, 0-based.
    pub position: u64,
    /// The record's name, where the file gives one.
    pub id: Option<String>,
    /// The reference allele, as `REF` spells it.
    pub reference: String,
    /// The alternate alleles, in the order `ALT` lists them, so an allele
    /// index `i` names `alternates[i - 1]`.
    pub alternates: Vec<String>,
    /// What each sample was called, in row order.
    pub calls: Vec<Genotype>,
}

impl GenotypeSite {
    /// A site at a 0-based `position`.
    pub fn new(
        position: u64,
        reference: impl Into<String>,
        alternates: impl IntoIterator<Item = impl Into<String>>,
        calls: impl Into<Vec<Genotype>>,
    ) -> Self {
        GenotypeSite {
            position,
            id: None,
            reference: reference.into(),
            alternates: alternates.into_iter().map(Into::into).collect(),
            calls: calls.into(),
        }
    }

    /// Names the record.
    pub fn id(mut self, name: impl Into<String>) -> Self {
        self.id = Some(name.into());
        self
    }

    /// What sample `row` was called here, and no call past the end of the
    /// list: nothing is invented to fill it.
    pub fn call(&self, row: usize) -> Genotype {
        self.calls
            .get(row)
            .copied()
            .unwrap_or_else(Genotype::not_called)
    }

    /// The allele an index names, as the file spells it.
    fn allele(&self, index: u16) -> &str {
        match index {
            0 => &self.reference,
            i => self
                .alternates
                .get(usize::from(i) - 1)
                .map_or("?", String::as_str),
        }
    }
}

/// Whether an alternate allele is a placeholder, which stands for any allele
/// the record does not name: `<NON_REF>` as GATK writes it and `<*>` as
/// bcftools does.
pub(crate) fn placeholder(alt: &str) -> bool {
    matches!(alt, "<NON_REF>" | "<*>")
}

/// An allele as a tooltip writes it: in full up to twelve characters, and
/// past that its first twelve and how many more there are.
fn spelled(allele: &str) -> String {
    let length = allele.chars().count();
    if length <= SPELLED {
        return allele.to_string();
    }
    let head: String = allele.chars().take(SPELLED).collect();
    format!("{head}+{}", length - SPELLED)
}

/// What one stretch of the band is drawn as, worked out once for every row.
#[derive(Debug, Clone, PartialEq)]
enum Cell {
    /// A site whose cell touches no other: a cell at the floor width,
    /// centred on its base, `left` pixels from the start of the band.
    Single { site: usize, left: f64, width: f64 },
    /// Sites whose cells would overlap, painted a pixel at a time.
    Pooled(Vec<Stretch>),
}

/// Whole pixels of a cluster, counted from the start of the band, and the
/// sites under the pixel of them that holds any.
#[derive(Debug, Clone, PartialEq)]
struct Stretch {
    from: f64,
    to: f64,
    sites: Range<usize>,
}

/// How one stretch of a row is drawn: no call, or a step of the eight, 0
/// being the reference.
type Class = Option<usize>;

/// What a figure's rows drew, which the key names and nothing else.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct Drawn {
    reference: bool,
    not_called: bool,
    /// Steps 1 to 8, at index 1 to 8.
    levels: [bool; LEVELS + 1],
    pooled: bool,
}

/// The calls of many samples at many sites, a row per sample.
///
/// ```
/// use karyon::{Figure, Genotype, GenotypeSite, GenotypeTrack, Region};
///
/// let samples = vec!["S1".to_string(), "S2".to_string(), "S3".to_string()];
/// let sites = vec![
///     GenotypeSite::new(1_100, "C", ["T"], vec![
///         Genotype::diploid(0, 1),
///         Genotype::diploid(1, 1),
///         Genotype::diploid(0, 0),
///     ]),
///     GenotypeSite::new(2_300, "G", ["A"], vec![
///         Genotype::diploid(0, 0),
///         Genotype::parse("./.").unwrap(),
///         Genotype::diploid(0, 1),
///     ]),
/// ];
///
/// let svg = Figure::new(Region::new("chr1", 1_000, 3_000).unwrap())
///     .push(GenotypeTrack::new(samples, sites).label("cohort"))
///     .to_svg();
/// assert!(svg.contains("S2, 1,101, C&gt;T: T/T (alternate)"));
/// ```
#[derive(Debug, Clone)]
pub struct GenotypeTrack {
    samples: Vec<String>,
    sites: Vec<GenotypeSite>,
    label: Option<String>,
    row_height: f64,
    row_gap: f64,
    min_cell_width: f64,
    max_rows: Option<usize>,
    show_names: bool,
    color: Option<String>,
    tree: Option<Tree>,
    tree_width: f64,
    tree_shape: TreeShape,
    traits: Traits,
}

impl GenotypeTrack {
    /// A row per sample in `samples`, and a cell per site.
    ///
    /// `sites[j].calls[i]` is what sample `i` was called at site `j`. The
    /// sites are put in order of position, keeping the order of two at one
    /// position, since a VCF need not be sorted and a cluster is worked out
    /// from neighbours.
    pub fn new(samples: impl Into<Vec<String>>, sites: impl Into<Vec<GenotypeSite>>) -> Self {
        let mut sites = sites.into();
        sites.sort_by_key(|site| site.position);
        GenotypeTrack {
            samples: samples.into(),
            sites,
            label: None,
            row_height: 11.0,
            row_gap: 1.0,
            min_cell_width: 3.0,
            max_rows: Some(40),
            show_names: true,
            color: None,
            tree: None,
            tree_width: 90.0,
            tree_shape: TreeShape::Phylogram,
            traits: Traits::default(),
        }
    }

    /// Sets the text shown in the left gutter.
    pub fn label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Sets the height of one row.
    pub fn row_height(mut self, height: f64) -> Self {
        self.row_height = height.max(1.0);
        self
    }

    /// Sets the gap between rows, which is left in the page colour.
    pub fn row_gap(mut self, gap: f64) -> Self {
        self.row_gap = gap.max(0.0);
        self
    }

    /// Sets how wide a cell is drawn at its narrowest, which is also how
    /// close two sites can be before their cells are pooled.
    pub fn min_cell_width(mut self, pixels: f64) -> Self {
        self.min_cell_width = pixels.max(0.5);
        self
    }

    /// Caps how many sample rows are drawn; `None` draws every one.
    pub fn max_rows(mut self, rows: Option<usize>) -> Self {
        self.max_rows = rows.map(|rows| rows.max(1));
        self
    }

    /// Draws or hides the sample names.
    pub fn show_names(mut self, show: bool) -> Self {
        self.show_names = show;
        self
    }

    /// Sets the hue of an alternate call, the theme's ink by default.
    ///
    /// The ink rather than the accent, because the accent is the first colour
    /// of the palette a `--traits` strip deals from: beside a strip of
    /// lineages, the alternate calls were the colour of one lineage, and a
    /// reader matching colours read a genotype as a lineage. A call is a
    /// quantity, how many copies are not the reference, and is drawn in the
    /// one colour no category takes.
    pub fn color(mut self, color: impl Into<String>) -> Self {
        self.color = Some(color.into());
        self
    }

    /// Draws a phylogeny beside the rows and puts the rows in the order of
    /// its tips.
    ///
    /// Sorted as the header lists them, the alleles a clade shares are a
    /// speckle down the rows; sorted by descent they are a block. A sample
    /// the tree does not name keeps its row at the bottom rather than
    /// vanishing, since a row dropped from a figure without a word is worse
    /// than a row out of order.
    pub fn tree(mut self, tree: Tree) -> Self {
        let order = leaf_order(&tree, &self.samples);
        self.samples = order
            .iter()
            .map(|index| self.samples[*index].clone())
            .collect();
        for site in &mut self.sites {
            site.calls = order.iter().map(|index| site.call(*index)).collect();
        }
        self.tree = Some(tree);
        self
    }

    /// Sets how much of the strip the tree gets, in pixels.
    pub fn tree_width(mut self, width: f64) -> Self {
        self.tree_width = width.max(0.0);
        self
    }

    /// Chooses a phylogram or a cladogram for the tree beside the rows.
    pub fn tree_shape(mut self, shape: TreeShape) -> Self {
        self.tree_shape = shape;
        self
    }

    /// Attaches metadata columns drawn between the sample names and the
    /// cells, joined by name, so they follow whatever order the rows are in.
    pub fn traits(mut self, traits: Traits) -> Self {
        self.traits = traits;
        self
    }

    /// The samples, in the order their rows are drawn.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// The sites, in order of position.
    pub fn sites(&self) -> &[GenotypeSite] {
        &self.sites
    }

    /// How many rows are drawn and how many the cap hid.
    pub fn visible_rows(&self) -> (usize, usize) {
        match self.max_rows {
            Some(cap) if self.samples.len() > cap => (cap, self.samples.len() - cap),
            _ => (self.samples.len(), 0),
        }
    }

    /// The metadata columns attached to this track.
    pub fn attached_traits(&self) -> &Traits {
        &self.traits
    }

    /// The tree beside the rows, if there is one.
    pub fn attached_tree(&self) -> Option<&Tree> {
        self.tree.as_ref()
    }

    /// The sites on display, which are a run of the list since it is sorted.
    fn shown(&self, region: &Region) -> Range<usize> {
        let from = self
            .sites
            .partition_point(|site| site.position < region.start());
        let to = self
            .sites
            .partition_point(|site| site.position < region.end());
        from..to.max(from)
    }

    /// How many of the tree's tips have no row drawn here.
    fn tips_without_rows(&self) -> usize {
        self.tree.as_ref().map_or(0, |tree| {
            tree_beside_rows(tree, &self.samples, self.visible_rows().0).1
        })
    }

    /// The line under the rows that says what was left out, or nothing.
    fn note(&self) -> Option<String> {
        let mut notes = Vec::new();
        let (_, hidden) = self.visible_rows();
        if hidden > 0 {
            notes.push(format!("+{} more", group_thousands(hidden as u64)));
        }
        let without = self.tips_without_rows();
        if without > 0 {
            notes.push(crate::track::tips_without_rows(without));
        }
        (!notes.is_empty()).then(|| notes.join(", "))
    }

    /// Where every site of the window is drawn, worked out once for all the
    /// rows, and the same for the key as for the figure, since both ask it
    /// here.
    ///
    /// A site's cell is centred on its base at the floor width, the wider of
    /// `min_cell_width` and a base. Consecutive sites whose cells overlap
    /// chain into a cluster. Inside a cluster the sites are pooled by the
    /// pixel their centre falls in, and every pixel of the cluster, from the
    /// left edge of its first cell to the right edge of its last, belongs to
    /// the nearest pixel that holds a site, halfway between two going to the
    /// left one: so a cluster is painted from end to end with no gap, and no
    /// pixel is drawn twice. Pixels are counted from the start of the band,
    /// which is all the figure and its key agree on.
    fn cells(&self, region: &Region, px_per_bp: f64) -> Vec<Cell> {
        let shown = self.shown(region);
        if shown.is_empty() || !px_per_bp.is_finite() || px_per_bp <= 0.0 {
            return Vec::new();
        }
        let width = self.min_cell_width.max(px_per_bp);
        let start = region.start();
        let centre =
            |site: usize| (self.sites[site].position - start) as f64 * px_per_bp + px_per_bp / 2.0;

        let mut cells = Vec::new();
        let mut first = shown.start;
        while first < shown.end {
            let mut last = first + 1;
            while last < shown.end && centre(last) - centre(last - 1) < width {
                last += 1;
            }
            if last - first == 1 {
                cells.push(Cell::Single {
                    site: first,
                    left: centre(first) - width / 2.0,
                    width,
                });
            } else {
                cells.push(Cell::Pooled(self.pool(first..last, &centre, width)));
            }
            first = last;
        }
        cells
    }

    /// The stretches of one cluster: each pixel holding a site, with every
    /// pixel nearer it than any other.
    fn pool(
        &self,
        cluster: Range<usize>,
        centre: &dyn Fn(usize) -> f64,
        width: f64,
    ) -> Vec<Stretch> {
        let mut held: Vec<(f64, Range<usize>)> = Vec::new();
        for site in cluster.clone() {
            let pixel = centre(site).floor();
            match held.last_mut() {
                Some((at, sites)) if *at == pixel => sites.end = site + 1,
                _ => held.push((pixel, site..site + 1)),
            }
        }
        let from = (centre(cluster.start) - width / 2.0).floor();
        let to = (centre(cluster.end - 1) + width / 2.0).ceil();
        let mut stretches = Vec::with_capacity(held.len());
        let mut left = from;
        for (index, (pixel, sites)) in held.iter().enumerate() {
            let right = match held.get(index + 1) {
                Some((next, _)) => ((pixel + next) / 2.0).floor() + 1.0,
                None => to,
            };
            stretches.push(Stretch {
                from: left,
                to: right.max(left),
                sites: sites.clone(),
            });
            left = right;
        }
        stretches
    }

    /// How a stretch of a row is drawn: the share of alternate copies among
    /// the calls under it that are calls, in eight steps.
    fn class(&self, row: usize, sites: Range<usize>) -> Class {
        let (mut copies, mut alternate) = (0usize, 0usize);
        for site in &self.sites[sites] {
            let call = site.call(row);
            if call.is_called() {
                copies += call.copies();
                alternate += call.alternate();
            }
        }
        (copies > 0).then(|| level(alternate, copies))
    }

    /// What the rows on display drew, walked the way the drawing walks them.
    fn drawn(&self, cells: &[Cell]) -> Drawn {
        let mut drawn = Drawn::default();
        let (rows, _) = self.visible_rows();
        let mut mark = |class: Class| match class {
            None => drawn.not_called = true,
            Some(0) => drawn.reference = true,
            Some(step) => drawn.levels[step] = true,
        };
        for row in 0..rows {
            for cell in cells {
                match cell {
                    Cell::Single { site, .. } => mark(self.sites[*site].call(row).level()),
                    Cell::Pooled(stretches) => {
                        for stretch in stretches {
                            mark(self.class(row, stretch.sites.clone()));
                        }
                    }
                }
            }
        }
        drawn.pooled = rows > 0 && cells.iter().any(|cell| matches!(cell, Cell::Pooled(_)));
        drawn
    }

    /// The colour of a step of alternate copies, from a fifth of the way off
    /// the page at the first to the hue at all of them.
    fn step_color(step: usize, hue: &str, theme: &Theme) -> String {
        mix(
            theme.surface(),
            hue,
            0.2 + 0.8 * step as f64 / LEVELS as f64,
        )
    }

    /// The colour of a reference call's bar.
    ///
    /// Mixed from the muted ink rather than from the rule, as the quiet marks
    /// of a matrix or a panel of variable sites are, because those measure
    /// 1.06 and 1.13 to one against the light page and are close to not
    /// being there. Mixed from the muted ink the bar is 2.07 to one, and the
    /// pale cell of no call 1.31 to one, with 1.58 between the two of them as
    /// well as the difference of shape.
    fn reference_color(theme: &Theme) -> String {
        mix(theme.surface(), &theme.muted, 0.45)
    }

    /// The colour of a cell with no call.
    fn missing_color(theme: &Theme) -> String {
        mix(theme.surface(), &theme.muted, 0.18)
    }

    fn hue(&self, theme: &Theme) -> String {
        self.color
            .clone()
            .unwrap_or_else(|| theme.foreground.clone())
    }

    /// Paints one run of a row, `x` and `width` in the figure's pixels.
    #[allow(clippy::too_many_arguments)]
    fn paint(
        &self,
        ctx: &mut DrawContext<'_>,
        class: Class,
        x: f64,
        width: f64,
        top: f64,
        hue: &str,
    ) {
        let theme = ctx.theme;
        match class {
            // A reference call is a short bar: the quiet state, which is most
            // of a cohort and the one thing a reader can tell without asking.
            Some(0) => ctx.svg.rect(
                x,
                top + self.row_height * (1.0 - REFERENCE_BAR) / 2.0,
                width,
                self.row_height * REFERENCE_BAR,
                &Self::reference_color(theme),
            ),
            None => ctx
                .svg
                .rect(x, top, width, self.row_height, &Self::missing_color(theme)),
            Some(step) => {
                let color = Self::step_color(step, hue, theme);
                ctx.svg.rect(x, top, width, self.row_height, &color);
            }
        }
    }

    /// What a reader hovering a row is told: the sample, and how many of the
    /// sites on display it was called at and carries an alternate allele at.
    ///
    /// Counted over the sites the figure draws, with the same predicate, so
    /// the tooltip of a zoomed figure does not count calls that are off it.
    fn row_tooltip(&self, row: usize, shown: Range<usize>) -> String {
        let name = self.samples.get(row).map_or("", String::as_str);
        let total = shown.len();
        if total == 0 {
            return name.to_string();
        }
        let calls = self.sites[shown].iter().map(|site| site.call(row));
        let (called, carrying) = calls.fold((0usize, 0usize), |(called, carrying), call| {
            (
                called + usize::from(call.is_called()),
                carrying + usize::from(call.is_called() && call.alternate() > 0),
            )
        });
        format!(
            "{name}, {} of {} site{} called, {} carrying an alternate allele",
            group_thousands(called as u64),
            group_thousands(total as u64),
            if total == 1 { "" } else { "s" },
            group_thousands(carrying as u64),
        )
    }

    /// What a reader hovering one cell is told, for a call that carries an
    /// alternate allele: where the site is, what it is, and the alleles of
    /// the call spelled out, so `1/2` says which two.
    ///
    /// A call of more than two copies keeps the alleles of its first two, so
    /// its tooltip writes those, how many more copies there are, and how many
    /// of all of them are alternate, which is what its colour was drawn from.
    fn cell_tooltip(&self, site: &GenotypeSite, row: usize) -> String {
        let call = site.call(row);
        let name = self.samples.get(row).map_or("", String::as_str);
        let named: Vec<&str> = site
            .alternates
            .iter()
            .map(String::as_str)
            .filter(|alt| !placeholder(alt))
            .collect();
        let alternates = if named.is_empty() {
            site.alternates.join(",")
        } else {
            named.join(",")
        };
        let joint = if call.phased() { "|" } else { "/" };
        let mut bases: Vec<String> = Vec::new();
        for copy in 0..call.copies().min(2) {
            bases.push(match call.allele(copy) {
                Some(index) => spelled(site.allele(index)),
                None => ".".to_string(),
            });
        }
        let mut alleles = bases.join(joint);
        if call.copies() > 2 {
            alleles.push_str(&format!(
                " and {} more copies, {} of {} alternate",
                call.copies() - 2,
                call.alternate(),
                call.copies()
            ));
        }
        let state = match call.state() {
            GenotypeState::Alternate => "alternate",
            GenotypeState::Heterozygous => "heterozygous",
            GenotypeState::Reference => "reference",
            GenotypeState::NotCalled => "no call",
        };
        let deleted = (0..call.copies().min(2))
            .filter_map(|copy| call.allele(copy))
            .any(|index| site.allele(index) == "*");
        let id = site
            .id
            .as_deref()
            .map_or(String::new(), |id| format!(", {id}"));
        format!(
            "{name}, {}{id}, {}>{}: {alleles} ({state}{})",
            group_thousands(site.position.saturating_add(1)),
            spelled(&site.reference),
            alternates,
            if deleted {
                "; * is a deletion upstream that takes this base"
            } else {
                ""
            }
        )
    }
}

impl Track for GenotypeTrack {
    fn noun(&self) -> &str {
        "genotypes by sample"
    }

    fn height(&self, _scale: &Scale) -> f64 {
        let rows = self.visible_rows().0.max(1) as f64;
        // A line under the rows for what was left out, in the room the track
        // asked for, since anything drawn outside it is clipped away.
        let note = if self.note().is_some() {
            Theme::default().font_size + 2.0
        } else {
            0.0
        };
        rows * self.row_height
            + (rows - 1.0).max(0.0) * self.row_gap
            + self.traits.heading_height()
            + note
    }

    fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    /// The states the rows drew, each with the mark it is drawn as, while
    /// every cell is one call: a reference call, a heterozygous one, an
    /// alternate one and no call, and only those of them the figure holds.
    /// Once any pixel is an average of several calls, or a call is a share
    /// no diploid makes, the steps are named as a ramp instead, from the
    /// first alternate copy to all of them.
    fn key(
        &self,
        region: &Region,
        px_per_bp: f64,
        theme: &Theme,
    ) -> Option<crate::track::legend::Legend> {
        let cells = self.cells(region, px_per_bp);
        let drawn = self.drawn(&cells);
        let hue = self.hue(theme);
        let mut key = crate::track::legend::Legend::new();
        if drawn.reference {
            key = key.line("reference call", Self::reference_color(theme));
        }
        let half = LEVELS / 2;
        let ramp = drawn.pooled
            || (1..=LEVELS).any(|step| drawn.levels[step] && step != half && step != LEVELS);
        if ramp {
            if drawn.levels.iter().any(|drawn| *drawn) {
                key = key.ramp(
                    "alternate copies among the calls",
                    Self::step_color(1, &hue, theme),
                    Self::step_color(LEVELS, &hue, theme),
                    "some",
                    "all",
                );
            }
        } else {
            if drawn.levels[half] {
                key = key.key("heterozygous call", Self::step_color(half, &hue, theme));
            }
            if drawn.levels[LEVELS] {
                key = key.key("alternate call", Self::step_color(LEVELS, &hue, theme));
            }
        }
        if drawn.not_called {
            key = key.key("no call", Self::missing_color(theme));
        }
        (!key.is_empty()).then_some(key)
    }

    fn y_axis_width(&self, theme: &Theme) -> f64 {
        let tree = if self.tree.is_some() {
            self.tree_width
        } else {
            0.0
        };
        let strip = self.traits.strip_width();
        let (rows, _) = self.visible_rows();
        if !self.show_names || rows == 0 {
            return tree + strip;
        }
        let size = (theme.font_size - 2.0).min(self.row_height);
        let widest = self
            .samples
            .iter()
            .take(rows)
            .map(|name| text_width(name, size))
            .fold(0.0f64, f64::max);
        widest + 8.0 + tree + strip
    }

    fn draw(&self, ctx: &mut DrawContext<'_>) {
        let band = ctx.band;
        // Scaled, because the band was, as a matrix steps over its headings.
        let head = ctx.px(self.traits.heading_height());
        let strip = self.traits.strip_width();
        let (rows, _) = self.visible_rows();
        let pitch = self.row_height + self.row_gap;
        let name_size = (ctx.theme.font_size - 2.0).min(self.row_height);
        let hue = self.hue(ctx.theme);
        let x0 = ctx.scale.x0();
        let shown = self.shown(ctx.region);
        let cells = self.cells(ctx.region, ctx.scale.px_per_bp());

        // The tree takes the left of the strip and the names the right of it,
        // so a tip, its name and its row of cells sit on one line.
        if let Some(tree) = self
            .tree
            .as_ref()
            .and_then(|tree| tree_beside_rows(tree, &self.samples, rows).0)
        {
            let area = Rect {
                x: ctx.axis.x + 2.0,
                y: band.y + head,
                w: (self.tree_width - 6.0).max(1.0),
                h: (band.h - head).max(1.0),
            };
            draw_tree(
                ctx.svg,
                &tree,
                area,
                pitch,
                band.y + head + self.row_height / 2.0,
                TreeStyle {
                    shape: self.tree_shape,
                    color: &ctx.theme.foreground,
                    width: 1.1,
                    mirror: false,
                },
            );
        }

        for row in 0..rows {
            let top = band.y + head + row as f64 * pitch;
            ctx.svg.begin_titled(&self.row_tooltip(row, shown.clone()));
            for cell in &cells {
                match cell {
                    Cell::Single { site, left, width } => {
                        let site = &self.sites[*site];
                        let call = site.call(row);
                        // Only a call carrying an alternate allele is named.
                        // A reference call and no call are what a reader can
                        // tell from the mark, and they are most of a cohort:
                        // naming them would be most of the file for the least
                        // it could say.
                        let named = matches!(
                            call.state(),
                            GenotypeState::Heterozygous | GenotypeState::Alternate
                        );
                        if named {
                            ctx.svg.begin_titled(&self.cell_tooltip(site, row));
                        }
                        self.paint(ctx, call.level(), x0 + left, *width, top, &hue);
                        if named {
                            ctx.svg.end_group();
                        }
                    }
                    Cell::Pooled(stretches) => {
                        // Runs of one class are one rectangle, and a pooled
                        // pixel carries no tooltip of its own: the row's says
                        // what a reader hovering it can be told.
                        let mut at = 0;
                        while at < stretches.len() {
                            let class = self.class(row, stretches[at].sites.clone());
                            let mut end = at + 1;
                            while end < stretches.len()
                                && self.class(row, stretches[end].sites.clone()) == class
                            {
                                end += 1;
                            }
                            let from = stretches[at].from;
                            let to = stretches[end - 1].to;
                            self.paint(ctx, class, x0 + from, to - from, top, &hue);
                            at = end;
                        }
                    }
                }
            }
            if self.show_names && ctx.axis.w > strip {
                if let Some(name) = self.samples.get(row) {
                    ctx.svg.text(
                        ctx.axis.right() - strip - 4.0,
                        top + self.row_height / 2.0 + name_size * 0.35,
                        name,
                        &ctx.theme.muted,
                        name_size,
                        Anchor::End,
                    );
                }
            }
            ctx.svg.end_group();
        }

        // After the rows and outside their groups, so a cell of the strip
        // carries its own tooltip rather than its row's.
        let placed: Vec<(String, f64, f64)> = self
            .samples
            .iter()
            .take(rows)
            .enumerate()
            .map(|(row, name)| {
                (
                    name.clone(),
                    band.y + head + row as f64 * pitch,
                    self.row_height,
                )
            })
            .collect();
        self.traits.draw(
            ctx,
            Rect {
                x: ctx.axis.right() - strip,
                y: band.y,
                w: strip,
                h: band.h,
            },
            &placed,
        );

        if let Some(note) = self.note() {
            ctx.svg.text(
                band.right() - 3.0,
                band.bottom() - 2.0,
                &note,
                &ctx.theme.muted,
                ctx.theme.font_size - 1.0,
                Anchor::End,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::figure::Figure;
    use crate::track::legend::LegendItem;
    use crate::track::variant::{Variant, VariantTrack};

    /// One rectangle of a document.
    #[derive(Debug, Clone, PartialEq)]
    struct Mark {
        x: f64,
        y: f64,
        w: f64,
        h: f64,
        fill: String,
    }

    /// Every filled rectangle narrower than the page, in document order: the
    /// cells, and not the page behind them or the band's clip.
    fn marks(svg: &str) -> Vec<Mark> {
        svg.split("<rect ")
            .skip(1)
            .filter_map(|rest| {
                let element = format!(" {}", &rest[..rest.find("/>")?]);
                let attr = |name: &str| -> Option<String> {
                    let at = element.find(&format!(" {name}=\""))? + name.len() + 3;
                    Some(element[at..at + element[at..].find('"')?].to_string())
                };
                Some(Mark {
                    x: attr("x")?.parse().ok()?,
                    y: attr("y")?.parse().ok()?,
                    w: attr("width")?.parse().ok()?,
                    h: attr("height")?.parse().ok()?,
                    fill: attr("fill")?,
                })
            })
            .filter(|mark| mark.w < 500.0)
            .collect()
    }

    fn names(count: usize) -> Vec<String> {
        (1..=count).map(|n| format!("S{n}")).collect()
    }

    fn gt(text: &str) -> Genotype {
        Genotype::parse(text).unwrap()
    }

    fn site(position: u64, calls: &[&str]) -> GenotypeSite {
        GenotypeSite::new(
            position,
            "C",
            ["T"],
            calls.iter().map(|call| gt(call)).collect::<Vec<_>>(),
        )
    }

    fn drawn(region: Region, track: GenotypeTrack) -> String {
        Figure::new(region)
            .show_region_label(false)
            .push(track)
            .to_svg()
    }

    /// WCAG's contrast ratio between two `#rrggbb` colours.
    fn contrast(one: &str, other: &str) -> f64 {
        fn channel(value: f64) -> f64 {
            if value <= 0.04045 {
                value / 12.92
            } else {
                ((value + 0.055) / 1.055).powf(2.4)
            }
        }
        fn luminance(color: &str) -> f64 {
            let parse =
                |at: usize| u8::from_str_radix(&color[at..at + 2], 16).unwrap_or(0) as f64 / 255.0;
            0.2126 * channel(parse(1)) + 0.7152 * channel(parse(3)) + 0.0722 * channel(parse(5))
        }
        let (a, b) = (luminance(one), luminance(other));
        let (high, low) = if a > b { (a, b) } else { (b, a) };
        (high + 0.05) / (low + 0.05)
    }

    /// A cohort is tens of millions of these, so the size is the point.
    #[test]
    fn a_call_is_eight_bytes() {
        assert_eq!(std::mem::size_of::<Genotype>(), 8);
    }

    /// Twenty thousand sites of two hundred samples drawn through a matrix
    /// were 257,862,650 bytes. Drawn here, where the cells would overlap, the
    /// pixel is the unit, and forty rows of a pixel's worth of runs are under
    /// a megabyte.
    #[test]
    fn twenty_thousand_sites_of_two_hundred_samples_stay_small() {
        let mut state = 20_261_001u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        let sites: Vec<GenotypeSite> = (0..20_000u64)
            .map(|index| {
                let calls: Vec<Genotype> = (0..200)
                    .map(|_| match next() % 10 {
                        0 => Genotype::diploid(0, 1),
                        1 => Genotype::diploid(1, 1),
                        2 => Genotype::not_called(),
                        _ => Genotype::diploid(0, 0),
                    })
                    .collect();
                GenotypeSite::new(index * 50, "C", ["T"], calls)
            })
            .collect();
        let svg = drawn(
            Region::parse("chr1:1-1,000,000").unwrap(),
            GenotypeTrack::new(names(200), sites),
        );
        assert!(svg.len() < 2_000_000, "{} bytes", svg.len());
        let rects = svg.matches("<rect").count();
        assert!(rects < 40 * 1_000, "{rects} rectangles");
        // Forty rows, and the rest counted rather than drawn.
        assert!(svg.contains(">+160 more</text>"), "the cap went missing");
    }

    /// Alone, a call is a cell at the floor width centred on its base, which
    /// is where a lollipop for the same record stands: a column of calls is
    /// under its variant, and not half a cell to one side of it.
    #[test]
    fn an_isolated_site_is_a_cell_centred_on_its_base_at_the_floor_width() {
        let theme = Theme::light();
        let at = 2_000u64;
        let svg = Figure::new(Region::new("chr1", 1_000, 3_000).unwrap())
            .show_region_label(false)
            .push(VariantTrack::new(vec![Variant::new(at)]).show_scale(false))
            .push(GenotypeTrack::new(names(1), vec![site(at, &["1/1"])]).show_names(false))
            .to_svg();
        let circle = &svg[svg.find("<circle").unwrap()..];
        let start = circle.find("cx=\"").unwrap() + 4;
        let cx: f64 = circle[start..start + circle[start..].find('"').unwrap()]
            .parse()
            .unwrap();
        let cells: Vec<Mark> = marks(&svg)
            .into_iter()
            .filter(|mark| mark.fill == theme.foreground && mark.h == 11.0)
            .collect();
        assert_eq!(cells.len(), 1, "{cells:?}");
        let cell = &cells[0];
        // Two kilobases over seven hundred pixels is a third of a pixel a
        // base, so the floor is what binds.
        assert_eq!(cell.w, 3.0);
        assert!(
            (cell.x + cell.w / 2.0 - cx).abs() < 0.002,
            "the cell is centred at {}, the lollipop at {cx}",
            cell.x + cell.w / 2.0
        );
    }

    /// A reference call is a short bar, and a heterozygous call, an alternate
    /// call and no call are full cells in three colours, none of them the
    /// page and none of them each other, in either theme.
    #[test]
    fn reference_heterozygous_alternate_and_not_called_are_four_different_marks() {
        let one = site(2_000, &["0/0", "0/1", "1/1", "./."]);
        let svg = drawn(
            Region::new("chr1", 1_000, 3_000).unwrap(),
            GenotypeTrack::new(names(4), vec![one]).show_names(false),
        );
        let cells = marks(&svg);
        assert_eq!(cells.len(), 4, "{cells:?}");
        let (reference, rest) = cells.split_first().unwrap();
        assert!(
            (reference.h - 11.0 * REFERENCE_BAR).abs() < 0.01,
            "the reference is the short bar: {reference:?}"
        );
        assert!(rest.iter().all(|mark| mark.h == 11.0), "{rest:?}");

        for theme in [Theme::light(), Theme::dark()] {
            let hue = theme.foreground.clone();
            let colors = [
                GenotypeTrack::reference_color(&theme),
                GenotypeTrack::step_color(LEVELS / 2, &hue, &theme),
                GenotypeTrack::step_color(LEVELS, &hue, &theme),
                GenotypeTrack::missing_color(&theme),
            ];
            for (a, one) in colors.iter().enumerate() {
                assert_ne!(one, theme.surface(), "a mark the colour of the page");
                for other in &colors[a + 1..] {
                    assert_ne!(one, other, "two states drawn alike");
                }
            }
        }
        // The quiet states mixed from the rule were 1.06 and 1.13 to one
        // against the light page, close to not there at all.
        let light = Theme::light();
        let reference = GenotypeTrack::reference_color(&light);
        let missing = GenotypeTrack::missing_color(&light);
        assert!(contrast(&reference, light.surface()) > 2.0);
        assert!(contrast(&missing, light.surface()) > 1.3);
        assert!(contrast(&reference, &missing) > 1.5);
    }

    /// No record is the page, and no call is a cell: "nothing was called
    /// here" and "nothing is here" are different statements.
    #[test]
    fn no_site_draws_nothing_and_not_called_draws_a_cell() {
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let none = drawn(region.clone(), GenotypeTrack::new(names(2), Vec::new()));
        assert!(marks(&none).is_empty(), "{none}");
        let uncalled = drawn(
            region,
            GenotypeTrack::new(names(1), vec![site(2_000, &["."])]),
        );
        let cells = marks(&uncalled);
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].fill, GenotypeTrack::missing_color(&Theme::light()));
        assert_eq!(cells[0].h, 11.0);
    }

    /// Forty-one sites under one pixel, one sample heterozygous at one of
    /// them: the pixel is drawn at the first step of the hue, not as the
    /// reference that forty of the forty-one are.
    #[test]
    fn sites_closer_than_the_floor_are_pooled_and_any_alternate_is_at_least_the_first_step() {
        let theme = Theme::light();
        let region = Region::parse("chr1:1-1,000,000").unwrap();
        let pooled = |middle: &str| -> Vec<Mark> {
            let sites: Vec<GenotypeSite> = (0..41u64)
                .map(|at| site(500_000 + at, &[if at == 20 { middle } else { "0/0" }]))
                .collect();
            marks(&drawn(
                region.clone(),
                GenotypeTrack::new(names(1), sites).show_names(false),
            ))
        };
        let carrying = pooled("0/1");
        assert!(!carrying.is_empty());
        let first = GenotypeTrack::step_color(1, &theme.foreground, &theme);
        assert!(
            carrying
                .iter()
                .all(|mark| mark.fill == first && mark.h == 11.0),
            "{carrying:?}"
        );
        // With no alternate copy among them, the same pixel is the bar.
        let quiet = pooled("0/0");
        assert!(!quiet.is_empty());
        assert!(
            quiet
                .iter()
                .all(|mark| mark.fill == GenotypeTrack::reference_color(&theme)),
            "{quiet:?}"
        );
    }

    /// A cluster is painted from the left edge of its first cell to the
    /// right edge of its last, every pixel once: each run starts where the
    /// one before it ends.
    #[test]
    fn a_cluster_is_painted_from_its_first_cell_to_its_last_with_no_gap() {
        let region = Region::parse("chr1:1-10,000").unwrap();
        // Twenty bases apart is about a pixel and a half, under the floor.
        let sites: Vec<GenotypeSite> = [
            ("1/1", 2_000u64),
            ("0/0", 2_020),
            ("1/1", 2_040),
            ("./.", 2_060),
        ]
        .iter()
        .map(|(call, at)| site(*at, &[call]))
        .collect();
        let track = GenotypeTrack::new(names(1), sites).show_names(false);
        let cells = track.cells(&region, 700.0 / 10_000.0);
        assert!(
            matches!(cells.as_slice(), [Cell::Pooled(stretches)] if stretches.len() == 4),
            "{cells:?}"
        );
        let svg = drawn(region, track);
        let runs = marks(&svg);
        assert_eq!(runs.len(), 4, "{runs:?}");
        for pair in runs.windows(2) {
            assert!(
                (pair[0].x + pair[0].w - pair[1].x).abs() < 1e-6,
                "a gap or an overlap: {pair:?}"
            );
        }
        // Whole pixels from the start of the band, so each run is too.
        for run in &runs {
            assert!(
                (run.w - run.w.round()).abs() < 1e-6,
                "{run:?} is not whole pixels"
            );
        }
    }

    /// While every cell is one call the key names the states; once a pixel
    /// is an average of several, it is a ramp of the steps, and no call.
    #[test]
    fn the_key_names_states_while_every_cell_is_one_call_and_is_a_ramp_once_any_is_pooled() {
        let theme = Theme::light();
        let sites = vec![site(1_100, &["0/0", "0/1"]), site(1_600, &["1/1", "./."])];
        let track = GenotypeTrack::new(names(2), sites);
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let label = |item: &LegendItem| match item {
            LegendItem::Key { label, .. } | LegendItem::Ramp { label, .. } => label.clone(),
        };
        let key = track.key(&region, 0.35, &theme).unwrap();
        let labels: Vec<String> = key.items().iter().map(label).collect();
        assert_eq!(
            labels,
            [
                "reference call",
                "heterozygous call",
                "alternate call",
                "no call"
            ]
        );
        // Zoomed out until the two cells overlap, they are one pixel's worth.
        let wide = Region::parse("chr1:1-1,000,000").unwrap();
        let key = track.key(&wide, 0.0007, &theme).unwrap();
        assert!(
            key.items()
                .iter()
                .any(|item| matches!(item, LegendItem::Ramp { .. })),
            "{key:?}"
        );
        let labels: Vec<String> = key.items().iter().map(label).collect();
        assert!(
            !labels.contains(&"heterozygous call".to_string()),
            "{labels:?}"
        );
    }

    /// A haploid cohort has no heterozygous call to key, and a key naming one
    /// would be a colour drawn nowhere.
    #[test]
    fn the_key_names_only_the_states_drawn() {
        let theme = Theme::light();
        let track = GenotypeTrack::new(names(2), vec![site(1_100, &["0", "1"])]);
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let key = track.key(&region, 0.35, &theme).unwrap();
        let labels: Vec<String> = key
            .items()
            .iter()
            .filter_map(|item| match item {
                LegendItem::Key { label, .. } => Some(label.clone()),
                LegendItem::Ramp { .. } => None,
            })
            .collect();
        assert_eq!(labels, ["reference call", "alternate call"]);
        // A window with no record in it keys nothing.
        let empty = Region::new("chr1", 5_000, 6_000).unwrap();
        assert!(track.key(&empty, 0.7, &theme).is_none());
    }

    /// The key is of the rows on display: a state held only by rows past
    /// the cap is drawn nowhere, and keyed it would be a colour with no cell.
    #[test]
    fn a_state_only_the_rows_past_the_cap_hold_is_not_keyed() {
        let theme = Theme::light();
        let mut calls = vec!["0/1"; 40];
        calls.extend(["./."; 5]);
        let track = GenotypeTrack::new(names(45), vec![site(1_100, &calls)]);
        assert_eq!(track.visible_rows(), (40, 5));
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let labels = |track: &GenotypeTrack| -> Vec<String> {
            let key = track.key(&region, 0.35, &theme).unwrap();
            key.items()
                .iter()
                .map(|item| match item {
                    LegendItem::Key { label, .. } | LegendItem::Ramp { label, .. } => label.clone(),
                })
                .collect()
        };
        assert_eq!(labels(&track), ["heterozygous call"]);
        // Lifted, the rows that hold it are drawn, and it is keyed.
        assert_eq!(
            labels(&track.clone().max_rows(None)),
            ["heterozygous call", "no call"]
        );
    }

    /// A row counts what is on the page: three of six sites in the window
    /// are three, and a call off it is not counted as called or carried.
    #[test]
    fn a_row_tooltip_counts_only_the_sites_in_the_window() {
        let sites: Vec<GenotypeSite> = [100u64, 200, 300, 400, 500, 600]
            .iter()
            .map(|at| site(*at, &["0/1"]))
            .collect();
        let svg = drawn(
            Region::parse("chr1:1-350").unwrap(),
            GenotypeTrack::new(names(1), sites),
        );
        assert!(
            svg.contains("<title>S1, 3 of 3 sites called, 3 carrying an alternate allele</title>"),
            "{svg}"
        );
    }

    /// A call with a copy missing is no call, and the row does not count it
    /// as carrying an alternate allele either, though the copy it has is one:
    /// carrying is counted among the calls, as the cell's colour is.
    #[test]
    fn a_partly_called_genotype_is_not_counted_as_carrying() {
        assert_eq!(gt("./1").alternate(), 1, "the copy it has is alternate");
        let track = GenotypeTrack::new(names(1), vec![site(100, &["./1"]), site(200, &["0/1"])]);
        assert_eq!(
            track.row_tooltip(0, 0..2),
            "S1, 1 of 2 sites called, 1 carrying an alternate allele"
        );
    }

    /// Only a cell that is one call and carries an alternate allele is named:
    /// a reference call and no call are what the mark already says, and a
    /// pooled pixel is several calls, which its row's tooltip counts.
    #[test]
    fn only_a_single_cell_that_carries_an_alternate_is_named() {
        let sites = vec![site(1_100, &["0/0", "0/1", "./.", "1/1"])];
        let svg = drawn(
            Region::new("chr1", 1_000, 3_000).unwrap(),
            GenotypeTrack::new(names(4), sites),
        );
        assert_eq!(
            svg.matches("<title>").count(),
            4 + 2,
            "four rows, two calls"
        );
        assert!(svg.contains("<title>S2, 1,101, C&gt;T: C/T (heterozygous)</title>"));
        assert!(svg.contains("<title>S4, 1,101, C&gt;T: T/T (alternate)</title>"));

        let pooled: Vec<GenotypeSite> = (0..10u64).map(|at| site(500_000 + at, &["1/1"])).collect();
        let dense = drawn(
            Region::parse("chr1:1-1,000,000").unwrap(),
            GenotypeTrack::new(names(1), pooled),
        );
        assert_eq!(dense.matches("<title>").count(), 1, "{dense}");
    }

    /// A call of more copies keeps the alleles of two, and says so, with the
    /// count of its alternate copies that its colour was drawn from.
    #[test]
    fn a_polyploid_call_says_how_many_copies_it_did_not_spell_out() {
        let call = gt("0/1/1/2");
        assert_eq!((call.copies(), call.alternate()), (4, 3));
        assert_eq!(
            (call.allele(0), call.allele(1), call.allele(2)),
            (Some(0), Some(1), None)
        );
        let track = GenotypeTrack::new(
            names(1),
            vec![GenotypeSite::new(999, "C", ["T", "G"], vec![call])],
        );
        assert_eq!(
            track.cell_tooltip(&track.sites[0], 0),
            "S1, 1,000, C>T,G: C/T and 2 more copies, 3 of 4 alternate (heterozygous)"
        );
    }

    /// A phased call is written with the bar the file wrote it with, so a
    /// reader hovering it can tell which haplotype carries the allele.
    #[test]
    fn a_phased_call_keeps_its_bar_in_the_tooltip() {
        let track = GenotypeTrack::new(names(2), vec![site(999, &["0|1", "0/1"])]);
        assert_eq!(
            track.cell_tooltip(&track.sites[0], 0),
            "S1, 1,000, C>T: C|T (heterozygous)"
        );
        assert_eq!(
            track.cell_tooltip(&track.sites[0], 1),
            "S2, 1,000, C>T: C/T (heterozygous)"
        );
    }

    /// Any copy that is not the reference is alternate, whichever allele it
    /// names: the base a deletion upstream took away and the placeholder that
    /// names no base are both not the reference, and the tooltip says which.
    #[test]
    fn a_star_and_a_placeholder_are_alternate_copies_and_named_as_such() {
        let star = GenotypeSite::new(299, "A", ["G", "*"], vec![gt("2/2"), gt("0/2")]);
        assert_eq!(star.call(0).state(), GenotypeState::Alternate);
        assert_eq!(star.call(1).state(), GenotypeState::Heterozygous);
        let gvcf = GenotypeSite::new(699, "C", ["T", "<NON_REF>"], vec![gt("1/2")]);
        assert_eq!(gvcf.call(0).state(), GenotypeState::Alternate);
        let track = GenotypeTrack::new(names(2), vec![star, gvcf]);
        assert_eq!(
            track.cell_tooltip(&track.sites[0], 0),
            "S1, 300, A>G,*: */* (alternate; * is a deletion upstream that takes this base)"
        );
        // The placeholder is no allele the record names, so the record is
        // written without it, and the call that names it says so.
        assert_eq!(
            track.cell_tooltip(&track.sites[1], 0),
            "S1, 700, C>T: T/<NON_REF> (alternate)"
        );
    }

    /// Forty rows are drawn and the rest are counted, on a line the track
    /// asked for room for.
    #[test]
    fn rows_past_the_cap_are_counted_not_drawn() {
        let calls: Vec<&str> = vec!["0/1"; 45];
        let track = GenotypeTrack::new(names(45), vec![site(1_100, &calls)]);
        assert_eq!(track.visible_rows(), (40, 5));
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let scale = Scale::new(&region, 0.0, 700.0);
        let rows = 40.0 * 11.0 + 39.0;
        assert_eq!(
            track.height(&scale),
            rows + Theme::default().font_size + 2.0
        );
        let svg = drawn(region, track.clone());
        assert!(svg.contains(">+5 more</text>"), "{svg}");
        assert!(svg.contains(">S40</text>") && !svg.contains(">S41</text>"));
        // Lifted, every row is drawn, and there is nothing to say.
        let all = track.max_rows(None);
        assert_eq!(all.visible_rows(), (45, 0));
        assert_eq!(all.height(&scale), 45.0 * 11.0 + 44.0);
    }

    /// The tree orders the rows by its tips and takes each sample's calls
    /// with it; a sample the tree does not name keeps its row at the bottom.
    #[test]
    fn a_tree_orders_the_rows_and_a_sample_it_does_not_name_stays_at_the_bottom() {
        let samples: Vec<String> = ["a", "b", "z", "c"].iter().map(|s| s.to_string()).collect();
        let one = GenotypeSite::new(
            1_100,
            "C",
            ["T"],
            vec![gt("0"), gt("1"), gt("."), Genotype::haploid(0)],
        );
        let tree = Tree::parse_newick("((c:1,a:1):1,(b:1,d:1):1);").unwrap();
        let track = GenotypeTrack::new(samples, vec![one]).tree(tree);
        assert_eq!(track.samples(), ["c", "a", "b", "z"]);
        let states: Vec<GenotypeState> = (0..4)
            .map(|row| track.sites()[0].call(row).state())
            .collect();
        assert_eq!(
            states,
            [
                GenotypeState::Reference,
                GenotypeState::Reference,
                GenotypeState::Alternate,
                GenotypeState::NotCalled
            ]
        );
        let svg = drawn(Region::new("chr1", 1_000, 3_000).unwrap(), track);
        assert!(svg.contains("1 tip of the tree has no row"), "{svg}");
    }

    /// A tip of the tree with no row is said on the line under the rows, and
    /// the track asks room for that line with no row hidden by the cap as
    /// well: without it, the line is drawn over the last row.
    #[test]
    fn a_tip_with_no_row_is_given_its_line_under_the_rows() {
        let samples: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let tree = Tree::parse_newick("((c:1,a:1):1,(b:1,d:1):1);").unwrap();
        let bare = GenotypeTrack::new(samples, vec![site(1_100, &["0", "1", "0"])]);
        let beside = bare.clone().tree(tree);
        assert_eq!(beside.visible_rows(), (3, 0), "no row is hidden");
        let scale = Scale::new(&Region::new("chr1", 1_000, 3_000).unwrap(), 0.0, 700.0);
        assert_eq!(
            beside.height(&scale),
            bare.height(&scale) + Theme::default().font_size + 2.0
        );
    }

    /// A VCF need not be sorted, and a figure of one in any order is the
    /// figure of it sorted, byte for byte, every time.
    #[test]
    fn unsorted_sites_are_drawn_in_position_order_and_the_same_input_draws_the_same_document() {
        let unsorted = vec![
            site(2_500, &["1"]),
            site(1_100, &["0"]),
            site(1_900, &["."]),
        ];
        let sorted = vec![
            site(1_100, &["0"]),
            site(1_900, &["."]),
            site(2_500, &["1"]),
        ];
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let track = GenotypeTrack::new(names(1), unsorted);
        let positions: Vec<u64> = track.sites().iter().map(|site| site.position).collect();
        assert_eq!(positions, [1_100, 1_900, 2_500]);
        let first = drawn(region.clone(), track.clone());
        assert_eq!(first, drawn(region.clone(), track));
        assert_eq!(first, drawn(region, GenotypeTrack::new(names(1), sorted)));
    }

    /// A sheet beside the rows heads its columns above them, and the track
    /// asks for that room as well as the rows': without it, the rows are
    /// pushed down under the headings and the last of them past the foot of
    /// the band, where the clip takes it away.
    #[test]
    fn a_sheet_beside_the_rows_is_given_room_for_its_headings() {
        let held = crate::read::sheet::sheet("sample\tlineage\nS1\tL1\nS2\tL2\n").unwrap();
        let traits = Traits::from_sheet(&held).strips(["lineage"]);
        let head = traits.heading_height();
        assert!(head > 0.0, "a column with a heading");
        let region = Region::new("chr1", 1_000, 3_000).unwrap();
        let scale = Scale::new(&region, 0.0, 700.0);
        let bare = GenotypeTrack::new(names(2), vec![site(1_100, &["1", "1"])]);
        let beside = bare.clone().traits(traits);
        assert_eq!(beside.height(&scale), bare.height(&scale) + head);
    }

    #[test]
    fn an_empty_track_has_one_row_of_height() {
        let scale = Scale::new(&Region::new("chr1", 0, 100).unwrap(), 0.0, 100.0);
        assert_eq!(
            GenotypeTrack::new(Vec::new(), Vec::new()).height(&scale),
            11.0
        );
    }
}
