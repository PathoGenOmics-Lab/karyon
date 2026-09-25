//! The bridge the documentation site's playground runs the command line over.
//!
//! `karyon::cli::stack::build` takes a closure that answers with a source's
//! text rather than opening a path, which is what makes this possible at all: a
//! browser has no disk, and the same grammar a shell drives from one is driven
//! here from whatever the page is holding. Nothing in the library changes to
//! run here, and nothing here is in the library.
//!
//! # Why this is a crate of its own
//!
//! Two reasons, and they are the same reason twice. `karyon` has no
//! dependencies and forbids `unsafe`; passing a string between JavaScript and
//! wasm means raw pointers, and doing it with `wasm-bindgen` means a build
//! toolchain and a dependency tree. Both belong outside a library that has
//! neither, so the pointers are here, they are counted on one hand, and the
//! protocol they speak is written down below.
//!
//! # The protocol
//!
//! One buffer in, one buffer out, both length prefixed, every number a
//! little-endian `u32`, every string UTF-8. The caller allocates with
//! [`alloc`], writes into the memory wasm exports, calls [`render`], reads the
//! answer, and frees both with [`dealloc`].
//!
//! ```text
//! in   [u32 argc]  argc x ([u32 len][bytes])        the command line, one word each
//!      [u32 filec] filec x ([u32 len][name] [u32 len][body])
//!
//! out  [u8 ok]     1 for a figure, 0 for a message
//!      [u32 len][bytes]                              the SVG, or what went wrong
//! ```
//!
//! A framing rather than a text format because there is no parser here to
//! disagree with: a file may hold any byte, a path may hold a space, and a
//! command line word may hold both.
//!
//! # The committed figures
//!
//! The figures the site shows are files the examples write, and the functions
//! that build them are in `examples/figures/`, one file per example. Those
//! files are compiled in here as they stand, which is what lets a page draw a
//! committed figure itself instead of showing the file: in the page's own light
//! or dark, laid out at the width of the column it sits in, over whatever
//! window the reader has moved to. It is the code that wrote the file, so a
//! page that asks for nothing but the figure gets the file back byte for byte,
//! and a test below holds it to that.
//!
//! [`figure`] draws one. The buffer in names it by its file under `assets/`
//! without the `.svg`, and says what the page wants different; the buffer out
//! is the one [`render`] answers with.
//!
//! ```text
//! in   [u32 len][name]         the figure, as in example-genomewide
//!      [u32 len][theme]        light or dark
//!      [u32 len][background]   a colour for the page under it, or empty for the theme's
//!      [u32 width]             in pixels, or 0 for the figure's own
//!      [u32 len][region]       a locus as samtools writes one, or empty for its own
//!      [u32 len][prefix]       what every id in it starts with, or empty
//!
//! out  [u8 ok][u32 len][bytes] the SVG, or what went wrong
//! ```
//!
//! Two figures inlined into one page share one id space, so each is given a
//! prefix of its own, and a figure given one knows it is going inside another
//! document and leaves out the `<title>` and `<desc>` that would name it there.
//!
//! Only a stack of tracks has a width to change, and only one whose tracks lay
//! their marks along a region has a window to move. A sheet, a circle, a map
//! or a tree keeps its own and says nothing, so a page can ask every figure the
//! same question. The two figures that exist to show one theme, the dark
//! example and the sheet of the web profile, keep that theme the same way, and
//! its page colour with it, whatever the page is running in. Which figure
//! moves is what [`figures`] is for: it takes nothing, or an empty buffer, and
//! answers with every figure there is and the window each one is drawn over,
//! so a page knows which it can pan and zoom along the genome, and how far.
//!
//! ```text
//! out  [u8 ok][u32 len]        1, and the length of what follows
//!      [u32 count]
//!      count x ([u32 len][name] [u8 moves] [u32 len][region])
//! ```
//!
//! `moves` is 1 for a figure drawn over a region and 0 for one that is not,
//! whose region is then empty. Handed back to [`figure`] as it is, a region
//! draws the figure the file holds.
//!
//! Answering that for every figure is building every figure, about a fifth of
//! a second, and a page that shows three of them has no use for the other
//! forty. [`figure_region`] answers it for one: the buffer in is the name as
//! [`figure`] takes it, and the answer is that figure's entry of the list.
//!
//! ```text
//! in   [u32 len][name]
//! out  [u8 ok][u32 len]        1, and the length of what follows
//!      [u8 moves] [u32 len][region]
//! ```

use std::cell::RefCell;
use std::collections::BTreeMap;

use karyon::cli::{args, stack};
use karyon::{Region, Theme, Tree};

// The examples' own files, which the library's examples compile as well, and
// there they are held to the oldest compiler the library supports. This crate
// declares no such version, so without saying it here Clippy asks for
// functions newer than that, and a change the examples cannot make would fail
// the lint here and nowhere else. The number is `rust-version` in the root
// Cargo.toml.
#[clippy::msrv = "1.74"]
mod committed;

/// Hands the caller a buffer of `len` bytes to write into.
///
/// # Safety
///
/// The pointer is only valid until it is passed to [`dealloc`] with the same
/// length. Writing past `len` is the caller's own undefined behaviour.
#[no_mangle]
pub extern "C" fn alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::with_capacity(len);
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

/// Gives back a buffer [`alloc`] handed out, or an answer from any of the others.
///
/// # Safety
///
/// `ptr` must be one this module returned and `len` the length it was made
/// with, and neither may be used again afterwards.
#[no_mangle]
pub unsafe extern "C" fn dealloc(ptr: *mut u8, len: usize) {
    if !ptr.is_null() && len > 0 {
        drop(Vec::from_raw_parts(ptr, 0, len));
    }
}

/// Runs one command line over the files the caller supplied.
///
/// Returns a buffer in the shape the module documentation gives, whose first
/// byte says whether the rest is a figure or the reason there is not one. The
/// caller reads the length back out of the buffer and frees it with
/// [`dealloc`].
///
/// # Safety
///
/// `ptr` and `len` must describe a buffer written in the input shape above.
#[no_mangle]
pub unsafe extern "C" fn render(ptr: *const u8, len: usize) -> *mut u8 {
    let input = std::slice::from_raw_parts(ptr, len);
    match run(input) {
        Ok(svg) => answer(true, &svg),
        Err(message) => answer(false, &message),
    }
}

/// The whole of what the playground does, with the pointers already gone.
fn run(mut input: &[u8]) -> Result<String, String> {
    let argv = strings(&mut input).ok_or("the command line is not in the shape this expects")?;
    let count = number(&mut input).ok_or("the file list is not in the shape this expects")?;
    let mut files: BTreeMap<String, String> = BTreeMap::new();
    for _ in 0..count {
        let name = text(&mut input).ok_or("a file name is not in the shape this expects")?;
        let body = text(&mut input).ok_or("a file body is not in the shape this expects")?;
        files.insert(name, body);
    }

    let request = args::parse(&argv).map_err(|error| error.to_string())?;
    let invocation = match request {
        args::Request::Draw(invocation) => invocation,
        // A page has nowhere to print to and no exit code, so the two requests
        // that are not a figure are answered as text rather than performed.
        args::Request::Help => return Err("--help prints to a terminal".to_string()),
        // `karyon::VERSION` and not this crate's own, which is the shim's.
        args::Request::Version => return Err(format!("karyon {}", karyon::VERSION)),
    };

    stack::build_with(
        &invocation,
        |source| match source {
            args::Source::Path(path) => {
                let name = path.display().to_string();
                files.get(&name).cloned().ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        // The one error this front end phrases differently from
                        // a shell, because a page's files are a list a reader
                        // can see rather than a directory they have to go and
                        // look in.
                        format!("no such file here; this page holds {}", named(&files)),
                    )
                })
            }
            // There is no pipe into a browser tab. Saying so is better than
            // answering with the empty string, which reads as an empty file.
            args::Source::Stdin => Err(std::io::Error::other("nothing is piped into a page")),
        },
        remembered,
    )
    .map_err(|error| error.to_string())
}

/// Draws one of the committed figures the way the page asks for it.
///
/// Returns a buffer in the shape [`render`] returns, holding the SVG or the
/// reason there is not one, and freed the same way. The input is the shape the
/// module documentation gives under the committed figures.
///
/// # Safety
///
/// `ptr` and `len` must describe a buffer written in that input shape.
#[no_mangle]
pub unsafe extern "C" fn figure(ptr: *const u8, len: usize) -> *mut u8 {
    let input = std::slice::from_raw_parts(ptr, len);
    match drawn(input) {
        Ok(svg) => answer(true, &svg),
        Err(message) => answer(false, &message),
    }
}

/// The widest figure drawn, in pixels.
///
/// The number the command line holds `--width` to, `MAX_WIDTH` in
/// src/cli/args.rs, and for its reason: a width becomes a column of pixels per
/// band, and one wider than any page asks for an allocation that fails, which
/// in a module built to abort is a trap rather than a message. A test holds the
/// two numbers together.
const MAX_WIDTH: usize = 100_000;

/// The whole of what [`figure`] does, with the pointer already gone.
fn drawn(mut input: &[u8]) -> Result<String, String> {
    let name = text(&mut input).ok_or("the figure's name is not in the shape this expects")?;
    let scheme = text(&mut input).ok_or("the theme is not in the shape this expects")?;
    let background = text(&mut input).ok_or("the background is not in the shape this expects")?;
    let width = number(&mut input).ok_or("the width is not in the shape this expects")?;
    let window = text(&mut input).ok_or("the region is not in the shape this expects")?;
    let prefix = text(&mut input).ok_or("the id prefix is not in the shape this expects")?;

    let build =
        committed::builder(&name).ok_or_else(|| format!("no committed figure is called {name}"))?;
    let mut theme = match scheme.as_str() {
        "light" => Theme::light(),
        "dark" => Theme::dark(),
        _ => return Err(format!("a figure is drawn light or dark, not {scheme}")),
    };
    if !background.is_empty() {
        // Written into a `fill` attribute as it stands, the way the command
        // line writes `--color`, so it is refused on the same characters: any
        // of these would end the attribute early and leave a page that does
        // not parse, and no spelling of a colour needs one.
        if background.contains(['"', '\'', '<', '>', '&']) {
            return Err(format!("{background:?} is not a colour, as in '#1e2129'"));
        }
        theme.background = background;
    }
    let width = match width {
        0 => None,
        px if px > MAX_WIDTH => {
            return Err(format!(
                "a figure {px} pixels wide is wider than any page; the most is {MAX_WIDTH}"
            ))
        }
        px => Some(px as f64),
    };
    let region = if window.is_empty() {
        None
    } else {
        Some(Region::parse(&window).map_err(|error| error.to_string())?)
    };
    // An id is written as it stands and pointed at from `url(#...)`, so the
    // prefix keeps to the characters that can go in both without quoting.
    if !prefix
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
    {
        return Err(format!(
            "an id prefix is letters, digits, '-' and '_', and {prefix:?} is not"
        ));
    }

    Ok(build(&theme, width, region.as_ref()).to_svg_with_id_prefix(&prefix))
}

/// Lists the committed figures, and the window each one moves over.
///
/// Returns a buffer in the shape the module documentation gives under the
/// committed figures, freed like the one [`render`] returns. Nothing is read
/// from the input, which may be empty.
#[no_mangle]
pub extern "C" fn figures(_ptr: *const u8, _len: usize) -> *mut u8 {
    CATALOGUE.with(|kept| framed(true, kept.get_or_init(catalogue)))
}

/// Says whether one committed figure moves along the genome, and over what.
///
/// The answer is that figure's entry of what [`figures`] answers with, framed
/// the same way; a page asks this once for each figure it holds rather than
/// having every figure built to find out about a few.
///
/// # Safety
///
/// `ptr` and `len` must describe a buffer holding one length-prefixed name.
#[no_mangle]
pub unsafe extern "C" fn figure_region(ptr: *const u8, len: usize) -> *mut u8 {
    let mut input = std::slice::from_raw_parts(ptr, len);
    let Some(name) = text(&mut input) else {
        return answer(false, "the figure's name is not in the shape this expects");
    };
    let Some(build) = committed::builder(&name) else {
        return answer(false, &format!("no committed figure is called {name}"));
    };
    framed(true, &region_entry(build))
}

/// One figure's entry of the list: whether it moves, and over what.
fn region_entry(build: committed::Builder) -> Vec<u8> {
    let mut out = Vec::new();
    match build(&Theme::light(), None, None).region() {
        Some(region) => {
            out.push(1);
            put_text(&mut out, &region.to_string());
        }
        None => {
            out.push(0);
            put_text(&mut out, "");
        }
    }
    out
}

thread_local! {
    /// What [`figures`] answers with, worked out the first time it is asked.
    ///
    /// The answer cannot change while the module is loaded, and working it out
    /// is building every figure there is: about a fifth of a second in a
    /// browser, most of it the circular chromosome and the gallery that holds a
    /// second one, and a page that asks again should not pay for it again.
    static CATALOGUE: std::cell::OnceCell<Vec<u8>> = const { std::cell::OnceCell::new() };
}

/// The body of what [`figures`] answers with.
///
/// Each figure is built to be asked, rather than the answers being written
/// down beside the list, because a written answer is a second copy of the
/// region every builder already holds, and it would be the copy that went
/// stale. Asked, the figure answers with the region it was drawn over, and
/// that is the one thing a page cannot get wrong by handing it back.
fn catalogue() -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(&(committed::FIGURES.len() as u32).to_le_bytes());
    for (name, build) in committed::FIGURES {
        put_text(&mut out, name);
        out.extend_from_slice(&region_entry(*build));
    }
    out
}

/// Hands back the coordinates a phylogeny is drawn at, rather than a drawing.
///
/// # Why a second way out
///
/// The figures this crate draws are SVG, and an SVG of a million tips is not a
/// thing a browser will move under a hand: the elements alone run to millions.
/// A viewer that wants to fly over such a tree needs the positions and a camera
/// of its own, which is what every viewer of that size has.
///
/// What it does not need is a second layout. The numbers below are the ones
/// `Tree::layout` works out and the ones an SVG of the same tree is drawn from,
/// so a canvas and a figure from a shell are the same picture at different
/// resolutions, and there is still one place where a tree's shape is decided.
///
/// # The shape of the answer
///
/// One buffer, little-endian, every array `count` long:
///
/// ```text
/// [u8 1]                     ok, or 0 and a message as `render` gives
/// [u32 count]
/// [f32 x]      x count       across: depth from the root, or a position
/// [f32 y]      x count       down: the row, or a position
/// [u32 parent] x count       what this hangs from, 0xFFFFFFFF where nothing
/// [u32 start]  x count       where this node's name begins in the blob
/// [u32 len]    x count       how long it is, zero where it has none
/// [u32 bytes][u8 bytes]      the names, run together
/// [u32 tips][u32 tip x tips] the order to read the terminals in, or none
/// ```
///
/// The third number in is which projection to lay out for: 0 for the one with
/// a root and 1 for the one without. They are the same arrays read two ways.
/// Rooted, `x` is a depth, `y` is a row, and `parent` is the tree's own; the
/// tip list is empty, because the rows already say what order to read them in.
/// Unrooted, `x` and `y` are a position in the plane, `parent` is the
/// neighbour on the way back to the middle, and the tip list is the order the
/// terminals come round it, which is what stands in for rows there.
///
/// # Safety
///
/// As [`render`].
#[no_mangle]
pub unsafe extern "C" fn layout(ptr: *const u8, len: usize) -> *mut u8 {
    let input = std::slice::from_raw_parts(ptr, len);
    match positions(input) {
        Ok(buffer) => {
            let mut out = buffer;
            let at = out.as_mut_ptr();
            std::mem::forget(out);
            at
        }
        Err(message) => answer(false, &message),
    }
}

/// One column of traits, resolved to what a figure would draw: the key, the
/// label, the levels as (value, light colour, dark colour), and one level index
/// per node with `u32::MAX` where a node carries nothing.
struct Strip {
    key: String,
    label: String,
    levels: Vec<(String, String, String)>,
    of: Vec<u32>,
}

fn put_text(out: &mut Vec<u8>, text: &str) {
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(text.as_bytes());
}

/// The sheet resolved against the tree the way the command line resolves it.
///
/// A tree reads its strips out of its own annotations, so the sheet is copied
/// onto the tips it names first and the crate is then asked what it would draw.
/// Working the levels out here instead would be a second opinion about which
/// blue is which, and the whole point of this page is that it is not one.
fn resolve_strips(tree: &Tree, body: &str) -> Vec<Strip> {
    let Ok(held) = karyon::read::sheet::sheet(body) else {
        return Vec::new();
    };
    let mut tree = tree.clone();
    for name in tree.leaf_names() {
        let (Some(values), Some(node)) = (held.rows.get(&name), tree.node_named(&name)) else {
            continue;
        };
        let values: Vec<(String, karyon::AnnotationValue)> = values
            .iter()
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();
        if let Some(into) = tree.annotations_mut(node) {
            for (key, value) in values {
                into.insert(key, value);
            }
        }
    }
    let spread = karyon::track::traits::Traits::new(held.rows.clone()).spread(held.columns.clone());
    let mut track = karyon::track::tree::TreeTrack::new(tree);
    for column in spread.columns() {
        track = track.trait_column(column.clone());
    }
    let light = track.strips(&karyon::Theme::light());
    let dark = track.strips(&karyon::Theme::dark());
    light
        .into_iter()
        .zip(dark)
        .map(|(pale, deep)| Strip {
            key: pale.key,
            label: pale.label,
            levels: pale
                .levels
                .iter()
                .zip(deep.levels.iter())
                .map(|(one, other)| (one.value.clone(), one.color.clone(), other.color.clone()))
                .collect(),
            of: pale
                .of
                .iter()
                .map(|at| at.map_or(u32::MAX, |index| index as u32))
                .collect(),
        })
        .collect()
}

/// The whole of what [`layout`] does, with the pointers already gone.
fn positions(mut input: &[u8]) -> Result<Vec<u8>, String> {
    let name = text(&mut input).ok_or("the file name is not in the shape this expects")?;
    let body = text(&mut input).ok_or("the file body is not in the shape this expects")?;
    let cladogram = number(&mut input).unwrap_or(0) == 1;
    let rootless = number(&mut input).unwrap_or(0) == 1;
    // The sheet, if one was dropped. Empty means there is none, which is not
    // the same as one that turned out to hold nothing.
    let _sheet_name = text(&mut input).unwrap_or_default();
    let sheet_body = text(&mut input).unwrap_or_default();

    let tree = match remembered(&name, body.trim()) {
        Some(tree) => tree,
        None => {
            Tree::parse_newick(body.trim()).map_err(|cause| format!("--tree {name}: {cause}"))?
        }
    };
    let nodes = tree.nodes();
    let count = nodes.len();

    // Both walks come back in the order they were made and the page wants them
    // by node, so they go back where they belong.
    let mut x = vec![0f32; count];
    let mut y = vec![0f32; count];
    let mut hangs: Vec<u32> = nodes
        .iter()
        .map(|clade| clade.parent.map_or(u32::MAX, |at| at as u32))
        .collect();
    let mut order: Vec<u32> = Vec::new();
    if rootless {
        let laid = tree.unrooted(cladogram);
        // A node with no position is one the walk could not reach, which means
        // the file is in more than one piece. It keeps the parent it had.
        for spot in &laid.spots {
            if spot.node < count {
                x[spot.node] = spot.x as f32;
                y[spot.node] = spot.y as f32;
                hangs[spot.node] = spot.toward.map_or(u32::MAX, |at| at as u32);
            }
        }
        order = laid.terminals.iter().map(|tip| *tip as u32).collect();
    } else {
        for placement in &tree.layout(cladogram) {
            if placement.node < count {
                x[placement.node] = placement.depth as f32;
                y[placement.node] = placement.row as f32;
            }
        }
    }

    // The trait columns, resolved by the crate rather than by the page: which
    // level every node is at and which colour that level gets, in both schemes,
    // so the strips beside the names are the strips the figure would draw.
    let strips = if sheet_body.trim().is_empty() {
        Vec::new()
    } else {
        resolve_strips(&tree, &sheet_body)
    };

    let mut names: Vec<u8> = Vec::new();
    let mut start = Vec::with_capacity(count);
    let mut length = Vec::with_capacity(count);
    for clade in nodes {
        start.push(names.len() as u32);
        match &clade.name {
            Some(text) => {
                length.push(text.len() as u32);
                names.extend_from_slice(text.as_bytes());
            }
            None => length.push(0u32),
        }
    }

    // Exactly the room the answer takes, because the caller frees it by that
    // size: a Vec that had to grow would hand back a capacity the free does not
    // match, and freeing a block by the wrong size traps the whole module.
    let mut out = Vec::with_capacity(1 + 4 + count * 20 + 4 + names.len() + 4 + order.len() * 4);
    out.push(1u8);
    out.extend_from_slice(&(count as u32).to_le_bytes());
    for value in &x {
        out.extend_from_slice(&value.to_le_bytes());
    }
    for value in &y {
        out.extend_from_slice(&value.to_le_bytes());
    }
    for value in &hangs {
        out.extend_from_slice(&value.to_le_bytes());
    }
    for value in &start {
        out.extend_from_slice(&value.to_le_bytes());
    }
    for value in &length {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out.extend_from_slice(&(names.len() as u32).to_le_bytes());
    out.extend_from_slice(&names);
    out.extend_from_slice(&(order.len() as u32).to_le_bytes());
    for value in &order {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out.extend_from_slice(&(strips.len() as u32).to_le_bytes());
    for strip in &strips {
        put_text(&mut out, &strip.key);
        put_text(&mut out, &strip.label);
        out.extend_from_slice(&(strip.levels.len() as u32).to_le_bytes());
        for level in &strip.levels {
            put_text(&mut out, &level.0);
            put_text(&mut out, &level.1);
            put_text(&mut out, &level.2);
        }
        for at in &strip.of {
            out.extend_from_slice(&at.to_le_bytes());
        }
    }
    // The caller frees this by its length, so the block it frees has to be
    // exactly that long: a Vec with room to spare hands back a capacity the
    // free does not match, and freeing by the wrong size traps the module.
    out.shrink_to_fit();
    debug_assert_eq!(
        out.len(),
        out.capacity(),
        "the answer is freed by its length"
    );
    Ok(out)
}

/// The phylogenies read lately, kept for the next call.
///
/// A shell runs the program once and reads each file once. A page runs it again
/// on every move, and reading the file is most of the work: a million tip tree
/// takes 361 ms to read of a 578 ms figure in a browser, against 189 for
/// folding it and drawing sixty rows of it.
///
/// Two of them, because the figures that take a phylogeny at all take either
/// one or two, and a single slot would have a tanglegram evicting its own left
/// tree with its right one on every frame and never hitting.
///
/// The text is kept beside the tree and compared in full rather than hashed. A
/// hash of twenty four megabytes costs about what it saves, and the failure it
/// would leave behind is the worst kind: the tree you had before, drawn under
/// the name of the one you asked for. Comparing is a length check and a
/// memcmp, and it is exact.
const KEPT: usize = 2;

/// Below this a tree is read again rather than kept.
///
/// Keeping one costs a copy of the text and a copy of the tree, which is worth
/// paying when reading it is a third of a second and not when it is a tenth of
/// a millisecond. The playground's own examples are all far under this, so
/// nothing there pays for a viewer's benefit.
const WORTH_KEEPING: usize = 1 << 20;

thread_local! {
    static TREES: RefCell<Vec<(String, String, Tree)>> = const { RefCell::new(Vec::new()) };
}

/// Answers with a tree already read, when it is the same one.
fn remembered(name: &str, text: &str) -> Option<Tree> {
    if text.len() < WORTH_KEEPING {
        return None;
    }
    TREES.with(|kept| {
        let mut kept = kept.borrow_mut();
        if let Some(at) = kept
            .iter()
            .position(|(had, body, _)| had == name && body == text)
        {
            // Most recently used first, so two trees in turn both stay.
            let entry = kept.remove(at);
            let tree = entry.2.clone();
            kept.insert(0, entry);
            return Some(tree);
        }
        let tree = Tree::parse_newick(text).ok()?;
        let answer = tree.clone();
        kept.insert(0, (name.to_string(), text.to_string(), tree));
        kept.truncate(KEPT);
        Some(answer)
    })
}

/// The files a page is holding, for the error that says one is missing.
fn named(files: &BTreeMap<String, String>) -> String {
    if files.is_empty() {
        return "no files".to_string();
    }
    files.keys().cloned().collect::<Vec<_>>().join(", ")
}

/// Reads a little-endian `u32` off the front.
fn number(input: &mut &[u8]) -> Option<usize> {
    let (head, rest) = input.split_at_checked(4)?;
    *input = rest;
    Some(u32::from_le_bytes(head.try_into().ok()?) as usize)
}

/// Reads a length-prefixed UTF-8 string off the front.
fn text(input: &mut &[u8]) -> Option<String> {
    let len = number(input)?;
    let (head, rest) = input.split_at_checked(len)?;
    *input = rest;
    String::from_utf8(head.to_vec()).ok()
}

/// Reads a length-prefixed list of length-prefixed strings off the front.
fn strings(input: &mut &[u8]) -> Option<Vec<String>> {
    let count = number(input)?;
    (0..count).map(|_| text(input)).collect()
}

/// Packs an answer into the buffer shape the caller reads.
fn answer(ok: bool, body: &str) -> *mut u8 {
    framed(ok, body.as_bytes())
}

/// Packs any bytes into that shape: the flag, the length, and the bytes.
fn framed(ok: bool, bytes: &[u8]) -> *mut u8 {
    let mut out = Vec::with_capacity(5 + bytes.len());
    out.push(u8::from(ok));
    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
    out.extend_from_slice(bytes);
    let ptr = out.as_mut_ptr();
    std::mem::forget(out);
    ptr
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    /// Builds the input buffer the way the page does, so the protocol is
    /// exercised rather than described.
    fn packed(argv: &[&str], files: &[(&str, &str)]) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::new();
        fn push(out: &mut Vec<u8>, text: &str) {
            out.extend_from_slice(&(text.len() as u32).to_le_bytes());
            out.extend_from_slice(text.as_bytes());
        }
        out.extend_from_slice(&(argv.len() as u32).to_le_bytes());
        for word in argv {
            push(&mut out, word);
        }
        out.extend_from_slice(&(files.len() as u32).to_le_bytes());
        for (name, body) in files {
            push(&mut out, name);
            push(&mut out, body);
        }
        out
    }

    #[test]
    fn a_whole_command_line_runs_with_no_disk_under_it() {
        let input = packed(
            &["chr1:1-60", "--coverage", "depth.bg", "--label", "depth"],
            &[("depth.bg", "chr1\t0\t60\t7\n")],
        );
        let svg = run(&input).expect("a figure");
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("depth"));
    }

    #[test]
    fn a_file_the_page_is_not_holding_says_what_it_is_holding() {
        // A shell can be told to go and look. A page cannot, so the error names
        // the list instead of the directory.
        let input = packed(
            &["chr1:1-60", "--coverage", "missing.bg"],
            &[("depth.bg", "chr1\t0\t60\t7\n")],
        );
        let error = run(&input).unwrap_err();
        assert!(error.contains("this page holds depth.bg"), "{error}");
    }

    /// A page runs the walk a shell runs, so a genome handed to `--sequence`
    /// is read by the region here too: the record the region names is drawn,
    /// and a region naming none of them is refused in the shell's words.
    #[test]
    fn a_genome_is_read_by_the_region_here_as_it_is_at_a_shell() {
        let genome = format!(">chrA\n{}\n>chrB\n{}\n", "A".repeat(60), "G".repeat(60));
        let files = [("genome.fa", genome.as_str())];

        let svg =
            run(&packed(&["chrB:1-60", "--sequence", "genome.fa"], &files)).expect("a figure");
        assert!(svg.contains(">G</text>"), "chrB was not the record drawn");
        assert!(
            !svg.contains(">A</text>"),
            "chrA's bases reached a chrB figure"
        );

        let error = run(&packed(&["chrC:1-60", "--sequence", "genome.fa"], &files)).unwrap_err();
        assert_eq!(
            error,
            "--sequence genome.fa has no record called chrC; it has chrA, chrB"
        );
    }

    #[test]
    fn a_command_line_that_makes_no_sense_comes_back_as_a_message() {
        let error = run(&packed(&["--nonsense"], &[])).unwrap_err();
        assert!(error.contains("unknown flag"), "{error}");

        let error = run(&packed(&["chr1:1-60", "--coverage", "-"], &[])).unwrap_err();
        assert!(error.contains("piped into a page"), "{error}");
    }

    #[test]
    fn a_buffer_that_is_not_the_shape_this_expects_is_a_message_and_not_a_panic() {
        for cut in [0usize, 1, 3, 4, 7] {
            let whole = packed(&["chr1:1-60"], &[("a", "b")]);
            let error = run(&whole[..cut.min(whole.len())]).unwrap_err();
            assert!(error.contains("shape this expects"), "cut {cut}: {error}");
        }
    }

    #[test]
    fn the_answer_says_whether_it_is_a_figure_before_it_says_anything_else() {
        let packed = answer(true, "<svg/>");
        // Safety: the pointer came from `answer` two lines above and is freed
        // on the line after it is read, which is what the page does too.
        let (ok, body) = unsafe {
            let ok = *packed;
            let len = u32::from_le_bytes(
                std::slice::from_raw_parts(packed.add(1), 4)
                    .try_into()
                    .unwrap(),
            ) as usize;
            let body =
                String::from_utf8(std::slice::from_raw_parts(packed.add(5), len).to_vec()).unwrap();
            dealloc(packed, 5 + len);
            (ok, body)
        };
        assert_eq!(ok, 1);
        assert_eq!(body, "<svg/>");
    }

    /// The one way a cache of files can be wrong is by answering with the
    /// wrong file, and the way it happens is a name reused for new contents.
    /// So the text is compared and not only the name, and this checks it in
    /// both directions: the same text is answered from memory, and a different
    /// text under the same name is not.
    #[test]
    fn a_tree_is_only_reused_when_the_text_is_the_same() {
        let long = |tips: usize, prefix: &str| {
            let mut parts: Vec<String> = (0..tips).map(|i| format!("{prefix}{i}:0.1")).collect();
            while parts.len() > 1 {
                let mut up = Vec::new();
                let mut at = 0;
                while at + 1 < parts.len() {
                    up.push(format!("({},{}):0.1", parts[at], parts[at + 1]));
                    at += 2;
                }
                if parts.len() % 2 == 1 {
                    up.push(parts[parts.len() - 1].clone());
                }
                parts = up;
            }
            format!("{};", parts[0])
        };
        // Over the size worth keeping, or nothing is kept at all.
        let first = long(90_000, "a");
        let second = long(90_000, "b");
        assert!(
            first.len() > WORTH_KEEPING,
            "the fixture has to be big enough"
        );

        let names = |tree: &Tree| tree.leaf_names().join(",");
        let one = remembered("t.nwk", &first).expect("a tree comes back");
        let again = remembered("t.nwk", &first).expect("and comes back again");
        assert_eq!(names(&one), names(&again), "the same text is the same tree");

        let other = remembered("t.nwk", &second).expect("a different text is read afresh");
        assert!(
            names(&other).starts_with('b'),
            "the same name with new contents must not answer with the old tree"
        );

        // Two of them stay, which is what a tanglegram needs: it hands over a
        // left tree and a right tree on every frame, and one slot would have
        // each evicting the other and never hitting. Asked from outside, a
        // miss and a hit both answer with the right tree, so the store itself
        // is what says which happened.
        let held = TREES.with(|kept| kept.borrow().len());
        assert_eq!(held, 2, "both trees are kept, not just the last");
        let back = remembered("t.nwk", &first).expect("the first is still held");
        assert!(
            names(&back).starts_with('a'),
            "and it is still the right one"
        );

        // And a small one is not kept, so nothing pays for a copy it will not
        // use.
        assert!(remembered("small.nwk", "((a:0.1,b:0.1):0.1,c:0.1);").is_none());
    }
    /// A sheet dropped on the page has to reach the canvas with the colours the
    /// figure would give it, or the page draws one picture and exports another.
    #[test]
    fn a_sheet_comes_back_as_the_colours_the_figure_would_draw() {
        const TREE: &str = "((a:1,b:1):1,c:1);";
        const SHEET: &str = "name\tplace\na\tValencia\nb\tValencia\nc\tLisbon\n";
        let mut input: Vec<u8> = Vec::new();
        fn push(out: &mut Vec<u8>, text: &str) {
            out.extend_from_slice(&(text.len() as u32).to_le_bytes());
            out.extend_from_slice(text.as_bytes());
        }
        push(&mut input, "t.nwk");
        push(&mut input, TREE);
        input.extend_from_slice(&0u32.to_le_bytes());
        input.extend_from_slice(&0u32.to_le_bytes());
        push(&mut input, "s.tsv");
        push(&mut input, SHEET);

        let out = positions(&input).expect("a tree and a sheet");
        let mut at = 1usize;
        let take = |at: &mut usize, out: &[u8]| {
            let value = u32::from_le_bytes(out[*at..*at + 4].try_into().unwrap());
            *at += 4;
            value
        };
        let word = |at: &mut usize, out: &[u8]| {
            let len = u32::from_le_bytes(out[*at..*at + 4].try_into().unwrap()) as usize;
            *at += 4;
            let text = String::from_utf8(out[*at..*at + len].to_vec()).unwrap();
            *at += len;
            text
        };
        let count = take(&mut at, &out) as usize;
        at += count * 20;
        let blob = take(&mut at, &out) as usize;
        at += blob;
        let tips = take(&mut at, &out) as usize;
        at += tips * 4;

        let columns = take(&mut at, &out) as usize;
        assert_eq!(columns, 1, "one column in the sheet, one column back");
        assert_eq!(word(&mut at, &out), "place");
        assert_eq!(word(&mut at, &out), "place");
        let levels = take(&mut at, &out) as usize;
        assert_eq!(levels, 2, "two places");
        let mut named = Vec::new();
        for _ in 0..levels {
            let value = word(&mut at, &out);
            let light = word(&mut at, &out);
            let dark = word(&mut at, &out);
            assert!(
                light.starts_with('#'),
                "a light colour that is not one: {light}"
            );
            assert!(
                dark.starts_with('#'),
                "a dark colour that is not one: {dark}"
            );
            named.push((value, light, dark));
        }
        let mut of = Vec::with_capacity(count);
        for _ in 0..count {
            of.push(take(&mut at, &out));
        }
        assert_eq!(at, out.len(), "the buffer is exactly as long as it says");

        // And the colours are the crate's own palette, in the order it hands
        // them out, rather than any this file chose.
        let light = karyon::Theme::light();
        let dark = karyon::Theme::dark();
        for (at, (_, pale, deep)) in named.iter().enumerate() {
            assert_eq!(
                pale,
                light.color(at),
                "level {at} is not the light palette's"
            );
            assert_eq!(deep, dark.color(at), "level {at} is not the dark palette's");
        }

        let tree = Tree::parse_newick(TREE).unwrap();
        let a = tree.node_named("a").unwrap();
        let b = tree.node_named("b").unwrap();
        let c = tree.node_named("c").unwrap();
        assert_ne!(of[a], u32::MAX, "a tip the sheet names carries nothing");
        assert_eq!(
            of[a], of[b],
            "two tips in one place are at different levels"
        );
        assert_ne!(of[a], of[c], "two places share a level");
        assert_eq!(named[of[a] as usize].0, "Valencia");
        assert_eq!(named[of[c] as usize].0, "Lisbon");
    }

    /// The coordinates a page flies over have to be the coordinates the crate
    /// The same buffer read the other way, which is the one the page has no
    /// other source for: an unrooted tree has no rows to sort by, so if the
    /// order of the terminals does not come over the wire it does not exist.
    #[test]
    fn the_rootless_layout_carries_what_rows_would_have() {
        const TREE: &str = "[&U] (((a:1,b:1):1,c:1):1,(d:1,e:1):1,f:1);";
        let mut input: Vec<u8> = Vec::new();
        fn push(out: &mut Vec<u8>, text: &str) {
            out.extend_from_slice(&(text.len() as u32).to_le_bytes());
            out.extend_from_slice(text.as_bytes());
        }
        push(&mut input, "u.nwk");
        push(&mut input, TREE);
        input.extend_from_slice(&0u32.to_le_bytes());
        input.extend_from_slice(&1u32.to_le_bytes());

        let out = positions(&input).expect("a tree this small lays out");
        assert_eq!(out[0], 1);
        let mut at = 1usize;
        let take = |at: &mut usize, out: &[u8]| {
            let value = u32::from_le_bytes(out[*at..*at + 4].try_into().unwrap());
            *at += 4;
            value
        };
        let count = take(&mut at, &out) as usize;
        let tree = Tree::parse_annotated_newick(TREE).unwrap();
        assert_eq!(count, tree.nodes().len());

        let mut x = Vec::with_capacity(count);
        for _ in 0..count {
            x.push(f32::from_le_bytes(out[at..at + 4].try_into().unwrap()));
            at += 4;
        }
        let mut y = Vec::with_capacity(count);
        for _ in 0..count {
            y.push(f32::from_le_bytes(out[at..at + 4].try_into().unwrap()));
            at += 4;
        }
        let mut toward = Vec::with_capacity(count);
        for _ in 0..count {
            toward.push(take(&mut at, &out));
        }
        for _ in 0..count * 2 {
            take(&mut at, &out);
        }
        let blob = take(&mut at, &out) as usize;
        at += blob;
        let tips = take(&mut at, &out) as usize;
        let mut order = Vec::with_capacity(tips);
        for _ in 0..tips {
            order.push(take(&mut at, &out) as usize);
        }
        assert_eq!(
            take(&mut at, &out),
            0,
            "no sheet was sent, so no strips come back"
        );
        assert_eq!(at, out.len(), "the buffer is exactly as long as it says");

        let laid = tree.unrooted(false);
        assert_eq!(order, laid.terminals, "the tip order is not the crate's");
        assert_eq!(tips, laid.terminals.len());
        for spot in &laid.spots {
            assert!((x[spot.node] as f64 - spot.x).abs() < 1e-5, "x moved");
            assert!((y[spot.node] as f64 - spot.y).abs() < 1e-5, "y moved");
            assert_eq!(
                toward[spot.node],
                spot.toward.map_or(u32::MAX, |at| at as u32),
                "what the branch hangs from moved"
            );
        }
        // Exactly one node has nothing to go toward, and it is the middle.
        let loose: Vec<usize> = (0..count)
            .filter(|node| toward[*node] == u32::MAX)
            .collect();
        assert_eq!(
            loose,
            vec![laid.centre],
            "the middle is not the only loose end"
        );
        // And what comes over the wire is not the rooted layout wearing a hat.
        let rooted = tree.layout(false);
        let same = rooted
            .iter()
            .all(|placement| (y[placement.node] as f64 - placement.row).abs() < 1e-5);
        assert!(!same, "the rootless layout came back as rows");
    }

    /// draws at, or the canvas and the figure are two different pictures. This
    /// takes the buffer apart again and checks it against what `Tree::layout`
    /// says, which is the one place the shape is decided.
    #[test]
    fn the_layout_buffer_says_what_the_crate_says() {
        const TREE: &str = "((a:0.5,b:0.25):0.25,c:1.0);";
        let mut input: Vec<u8> = Vec::new();
        fn push(out: &mut Vec<u8>, text: &str) {
            out.extend_from_slice(&(text.len() as u32).to_le_bytes());
            out.extend_from_slice(text.as_bytes());
        }
        push(&mut input, "t.nwk");
        push(&mut input, TREE);
        input.extend_from_slice(&0u32.to_le_bytes());

        let out = positions(&input).expect("a tree this small lays out");
        assert_eq!(out[0], 1, "the first byte says it worked");
        let mut at = 1usize;
        let take = |at: &mut usize, out: &[u8]| {
            let value = u32::from_le_bytes(out[*at..*at + 4].try_into().unwrap());
            *at += 4;
            value
        };
        let count = take(&mut at, &out) as usize;

        let tree = Tree::parse_newick(TREE).unwrap();
        assert_eq!(count, tree.nodes().len());

        let mut x = Vec::with_capacity(count);
        for _ in 0..count {
            x.push(f32::from_le_bytes(out[at..at + 4].try_into().unwrap()));
            at += 4;
        }
        let mut y = Vec::with_capacity(count);
        for _ in 0..count {
            y.push(f32::from_le_bytes(out[at..at + 4].try_into().unwrap()));
            at += 4;
        }
        let mut parent = Vec::with_capacity(count);
        for _ in 0..count {
            parent.push(take(&mut at, &out));
        }
        let mut start = Vec::with_capacity(count);
        for _ in 0..count {
            start.push(take(&mut at, &out) as usize);
        }
        let mut length = Vec::with_capacity(count);
        for _ in 0..count {
            length.push(take(&mut at, &out) as usize);
        }
        let blob = take(&mut at, &out) as usize;
        let names = out[at..at + blob].to_vec();
        at += blob;
        let names = &names[..];
        let tips = take(&mut at, &out) as usize;
        assert_eq!(
            tips, 0,
            "a rooted layout sends no tip order: the rows are one"
        );
        let strips = take(&mut at, &out) as usize;
        assert_eq!(strips, 0, "no sheet was sent, so no strips come back");
        assert_eq!(at, out.len(), "the buffer is exactly as long as it says");

        // Every coordinate is the one the crate worked out.
        for placement in tree.layout(false) {
            assert!(
                (x[placement.node] as f64 - placement.depth).abs() < 1e-5,
                "node {} depth",
                placement.node
            );
            assert!(
                (y[placement.node] as f64 - placement.row).abs() < 1e-5,
                "node {} row",
                placement.node
            );
        }
        // And so is every parent and every name.
        for (node, clade) in tree.nodes().iter().enumerate() {
            assert_eq!(parent[node], clade.parent.map_or(u32::MAX, |at| at as u32));
            let held =
                std::str::from_utf8(&names[start[node]..start[node] + length[node]]).unwrap();
            assert_eq!(
                held,
                clade.name.as_deref().unwrap_or(""),
                "node {node} name"
            );
        }

        // A cladogram is a different question and gives a different answer.
        let mut asked: Vec<u8> = Vec::new();
        push(&mut asked, "t.nwk");
        push(&mut asked, TREE);
        asked.extend_from_slice(&1u32.to_le_bytes());
        let flat = positions(&asked).expect("a cladogram lays out too");
        assert_ne!(flat, out, "asking for a cladogram changes the depths");
    }

    /// A request for [`figure`], packed the way a page packs one.
    fn ask(
        name: &str,
        theme: &str,
        background: &str,
        width: u32,
        region: &str,
        prefix: &str,
    ) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::new();
        put_text(&mut out, name);
        put_text(&mut out, theme);
        put_text(&mut out, background);
        out.extend_from_slice(&width.to_le_bytes());
        put_text(&mut out, region);
        put_text(&mut out, prefix);
        out
    }

    /// A figure drawn with nothing asked of it but a theme.
    fn plain(name: &str, theme: &str) -> String {
        drawn(&ask(name, theme, "", 0, "", "")).unwrap_or_else(|error| panic!("{name}: {error}"))
    }

    /// What `assets/` holds for a figure, the file the examples wrote.
    fn committed_file(name: &str) -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../assets")
            .join(format!("{name}.svg"));
        std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
    }

    /// The text between one `key` and the `end` after it, for every `key`.
    fn after<'a>(svg: &'a str, key: &str, end: char) -> Vec<&'a str> {
        svg.match_indices(key)
            .map(|(at, found)| {
                let rest = &svg[at + found.len()..];
                &rest[..rest.find(end).unwrap_or(rest.len())]
            })
            .collect()
    }

    /// The words a figure writes, tick labels included, in the order written.
    fn words(svg: &str) -> Vec<&str> {
        svg.split("<text")
            .skip(1)
            .filter_map(|piece| {
                let open = piece.find('>')?;
                let close = piece.find("</text>")?;
                piece.get(open + 1..close)
            })
            .collect()
    }

    #[test]
    fn every_committed_figure_is_listed_and_nothing_else_is() {
        let folder = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../assets");
        let on_disk: BTreeSet<String> = std::fs::read_dir(&folder)
            .expect("assets/ is where the figures are committed")
            .filter_map(|entry| {
                let name = entry.ok()?.file_name().into_string().ok()?;
                name.strip_suffix(".svg").map(str::to_string)
            })
            .collect();
        let listed: BTreeSet<String> = committed::FIGURES
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();
        let unlisted: Vec<&String> = on_disk.difference(&listed).collect();
        let uncommitted: Vec<&String> = listed.difference(&on_disk).collect();
        assert!(
            unlisted.is_empty() && uncommitted.is_empty(),
            "committed and not listed: {unlisted:?}; listed and not committed: {uncommitted:?}"
        );
        assert_eq!(
            listed.len(),
            committed::FIGURES.len(),
            "a name listed twice"
        );
        assert!(
            committed::FIGURES
                .windows(2)
                .all(|pair| pair[0].0 < pair[1].0),
            "the list is not sorted by name"
        );
    }

    /// The claim the whole arrangement rests on: the page and the file are one
    /// figure. A figure asked for in light with nothing else changed is the
    /// file CI holds the examples to, byte for byte.
    #[test]
    fn asked_for_nothing_a_figure_is_the_file_it_stands_in_for() {
        for (name, _) in committed::FIGURES {
            let svg = plain(name, "light");
            let file = committed_file(name);
            if svg != file {
                let parted = svg
                    .bytes()
                    .zip(file.bytes())
                    .position(|(one, other)| one != other)
                    .unwrap_or(svg.len().min(file.len()));
                panic!(
                    "{name} is not assets/{name}.svg: they part at byte {parted} of {}",
                    file.len()
                );
            }
        }
    }

    #[test]
    fn every_figure_draws_in_both_themes_and_takes_the_one_it_is_given() {
        let light_palette = Theme::light().palette;
        for (name, _) in committed::FIGURES {
            let light = plain(name, "light");
            let dark = plain(name, "dark");
            for svg in [&light, &dark] {
                assert!(
                    svg.starts_with("<svg ") && svg.ends_with("</svg>"),
                    "{name}"
                );
            }
            // The two that exist to show one theme keep it, and every other
            // one is drawn in the page's.
            let keeps = matches!(*name, "example-dark" | "example-visual-system");
            assert_eq!(light == dark, keeps, "{name}: the theme asked for");
            if keeps {
                continue;
            }
            // Every page in it is dark, the sheet's own and each panel's: a
            // panel that kept the light theme paints a white rectangle in the
            // middle of a dark sheet.
            let pages = after(&dark, r#"<rect x="0" y="0" width=""#, '/');
            assert!(!pages.is_empty(), "{name} paints no page");
            for page in pages {
                assert!(
                    page.ends_with(&format!(r#"fill="{}""#, Theme::dark().background)),
                    "{name} in dark has a page that is not: {page}"
                );
            }
            // Not one of the light palette's colours is left in a dark figure,
            // which is what says that no colour in it was written down as a
            // literal the theme could not reach.
            for colour in &light_palette {
                assert!(
                    !dark.contains(colour.as_str()),
                    "{name} in dark still draws {colour}"
                );
            }
        }
    }

    #[test]
    fn a_region_half_as_wide_moves_the_figure_and_nothing_else_moves() {
        for (name, build) in committed::FIGURES {
            let whole = plain(name, "light");
            match build(&Theme::light(), None, None).region().cloned() {
                Some(own) => {
                    let half =
                        Region::new(own.seq(), own.start(), own.start() + (own.len() / 2).max(1))
                            .unwrap();
                    let moved = drawn(&ask(name, "light", "", 0, &half.to_string(), ""))
                        .unwrap_or_else(|error| panic!("{name} over {half}: {error}"));
                    assert!(
                        after(&moved, "viewBox=\"", '"')[0] != after(&whole, "viewBox=\"", '"')[0]
                            || words(&moved) != words(&whole),
                        "{name} over {half} is drawn as it is over {own}"
                    );
                }
                None => {
                    // A drawing with no window along a genome has none to
                    // move, so a region asked of it changes nothing at all.
                    let asked = drawn(&ask(name, "light", "", 0, "chr1:1-100", ""))
                        .unwrap_or_else(|error| panic!("{name}: {error}"));
                    assert_eq!(asked, whole, "{name} has no window and moved one anyway");
                }
            }
        }
    }

    /// Moving the window moves what is looked at, not what is there. The
    /// close-up is the overview's own data over sixty of its bases, so the
    /// overview moved onto those sixty has to spell out the same sequence: data
    /// laid from the left edge of whatever window is asked for would spell the
    /// first sixty bases instead, and look just as plausible.
    #[test]
    fn a_window_moved_over_the_overview_reads_the_bases_the_close_up_reads() {
        let bases = |svg: &str| -> String {
            words(svg)
                .into_iter()
                .filter(|word| matches!(*word, "A" | "C" | "G" | "T"))
                .collect()
        };
        let close = plain("example-zoom", "light");
        let window = committed::builder("example-zoom").unwrap()(&Theme::light(), None, None)
            .region()
            .unwrap()
            .to_string();
        let moved = drawn(&ask("example", "light", "", 0, &window, "")).unwrap();
        assert_eq!(
            bases(&close).len(),
            60,
            "the close-up spells its sixty bases"
        );
        assert_eq!(bases(&moved), bases(&close), "over {window}");
    }

    #[test]
    fn a_width_is_taken_by_a_stack_of_tracks_and_by_nothing_else() {
        for (name, build) in committed::FIGURES {
            let whole = plain(name, "light");
            let narrow = drawn(&ask(name, "light", "", 700, "", ""))
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            // A figure is the one drawing that has a plotting area to line up,
            // and the one that has a width of its own.
            if build(&Theme::light(), None, None)
                .content_anchor()
                .is_some()
            {
                assert_eq!(after(&narrow, " width=\"", '"')[0], "700", "{name}");
            } else {
                assert_eq!(narrow, whole, "{name} has no width and changed anyway");
            }
        }
    }

    #[test]
    fn a_background_is_the_page_under_the_figure_and_under_every_panel() {
        let dark = Theme::dark().background;
        for name in [
            "example-genomewide",
            "example-circular",
            "gallery",
            "example-maps",
        ] {
            let svg = drawn(&ask(name, "dark", "#1e2129", 0, "", ""))
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert!(svg.contains(r##"fill="#1e2129""##), "{name}");
            assert!(
                !svg.contains(&dark),
                "{name} still paints the theme's own page"
            );
        }
    }

    #[test]
    fn a_prefix_reaches_every_id_in_every_figure() {
        for (name, _) in committed::FIGURES {
            let svg = drawn(&ask(name, "light", "", 0, "", "f7-"))
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            let written = after(&svg, " id=\"", '"');
            for id in &written {
                assert!(id.starts_with("f7-"), "{name} wrote the id {id}");
            }
            for id in after(&svg, "url(#", ')') {
                assert!(id.starts_with("f7-"), "{name} points at {id}");
                assert!(
                    written.contains(&id),
                    "{name} points at {id}, which is not there"
                );
            }
            // Going inside a page, the figure itself does not name itself.
            let root = &svg[..svg.find('>').unwrap()];
            assert!(!root.contains("aria-labelledby"), "{name}: {root}");
        }
    }

    #[test]
    fn the_list_says_which_figures_move_and_over_what() {
        let body = catalogue();
        let mut at = &body[..];
        assert_eq!(number(&mut at), Some(committed::FIGURES.len()));
        for (name, build) in committed::FIGURES {
            assert_eq!(text(&mut at).as_deref(), Some(*name));
            let (moves, rest) = at.split_first().expect("a byte for whether it moves");
            at = rest;
            let region = text(&mut at).expect("a region, empty or not");
            let own = build(&Theme::light(), None, None)
                .region()
                .map(|region| region.to_string());
            assert_eq!(*moves == 1, own.is_some(), "{name}");
            assert_eq!(region, own.unwrap_or_default(), "{name}");
            // Handed back as it is, the region draws the file.
            if *moves == 1 {
                let svg = drawn(&ask(name, "light", "", 0, &region, "")).unwrap();
                assert!(
                    svg == committed_file(name),
                    "{name} over {region} is not the file"
                );
            }
        }
        assert!(at.is_empty(), "the list is exactly as long as it says");

        // And the export hands back that body, framed the way render's is.
        let packed = figures(std::ptr::null(), 0);
        // Safety: the pointer came from `figures` on the line above and is
        // freed by its length once read, which is what the page does too.
        let (ok, read) = unsafe {
            let ok = *packed;
            let len = u32::from_le_bytes(
                std::slice::from_raw_parts(packed.add(1), 4)
                    .try_into()
                    .unwrap(),
            ) as usize;
            let read = std::slice::from_raw_parts(packed.add(5), len).to_vec();
            dealloc(packed, 5 + len);
            (ok, read)
        };
        assert_eq!(ok, 1);
        assert_eq!(read, body);
    }

    #[test]
    fn one_figure_answers_with_its_own_entry_of_the_list() {
        // The same bytes as its entry in `figures`, for a figure that moves and
        // for one that does not, and a name nobody committed is a message.
        let body = catalogue();
        let mut at = &body[..];
        number(&mut at).unwrap();
        let mut entries = BTreeMap::new();
        for _ in committed::FIGURES {
            let name = text(&mut at).unwrap();
            let start = at;
            let (_, rest) = at.split_first().unwrap();
            at = rest;
            text(&mut at).unwrap();
            entries.insert(name, start[..start.len() - at.len()].to_vec());
        }
        let ask_region = |name: &str| {
            let mut input = Vec::new();
            put_text(&mut input, name);
            // Safety: the pointers are this module's own, used and freed the
            // way the page uses and frees them.
            unsafe {
                let into = alloc(input.len());
                std::ptr::copy_nonoverlapping(input.as_ptr(), into, input.len());
                let out = figure_region(into, input.len());
                dealloc(into, input.len());
                let ok = *out;
                let len = u32::from_le_bytes(
                    std::slice::from_raw_parts(out.add(1), 4)
                        .try_into()
                        .unwrap(),
                ) as usize;
                let read = std::slice::from_raw_parts(out.add(5), len).to_vec();
                dealloc(out, 5 + len);
                (ok, read)
            }
        };
        for name in ["example-genomewide", "example-circular"] {
            let (ok, read) = ask_region(name);
            assert_eq!(ok, 1, "{name}");
            assert_eq!(&read, &entries[name], "{name}");
        }
        let (ok, read) = ask_region("example-nothing");
        assert_eq!(ok, 0);
        assert!(String::from_utf8(read).unwrap().contains("example-nothing"));
    }

    #[test]
    fn a_figure_comes_back_through_the_pointers_a_page_uses() {
        let input = ask("example-zoom", "dark", "", 640, "", "zoom-");
        // Safety: every pointer here is one this module handed out, written
        // within its length and freed by it, in the order the page does it.
        let (ok, svg) = unsafe {
            let into = alloc(input.len());
            std::ptr::copy_nonoverlapping(input.as_ptr(), into, input.len());
            let out = figure(into, input.len());
            dealloc(into, input.len());
            let ok = *out;
            let len = u32::from_le_bytes(
                std::slice::from_raw_parts(out.add(1), 4)
                    .try_into()
                    .unwrap(),
            ) as usize;
            let svg =
                String::from_utf8(std::slice::from_raw_parts(out.add(5), len).to_vec()).unwrap();
            dealloc(out, 5 + len);
            (ok, svg)
        };
        assert_eq!(ok, 1, "{svg}");
        assert!(svg.starts_with(r#"<svg xmlns="http://www.w3.org/2000/svg" width="640" "#));
        assert!(svg.contains(&Theme::dark().background));
        assert!(svg.contains(r#"id="zoom-karyon-clip-0""#));
    }

    #[test]
    fn a_request_that_makes_no_sense_comes_back_as_a_message() {
        let refused = |input: Vec<u8>| drawn(&input).unwrap_err();
        assert!(
            refused(ask("example-nothing", "light", "", 0, "", "")).contains("no committed figure")
        );
        assert!(refused(ask("example", "sepia", "", 0, "", "")).contains("light or dark"));
        assert!(
            refused(ask("example", "light", "red\" onload=\"x", 0, "", ""))
                .contains("not a colour")
        );
        assert!(refused(ask("example", "light", "", 0, "chr1", "")).contains("chr1"));
        assert!(refused(ask("example", "light", "", 0, "", "a b")).contains("id prefix"));
        assert!(
            refused(ask("example", "light", "", 100_001, "", "")).contains("wider than any page")
        );
        // Cut anywhere, a request is a message and never a panic.
        let whole = ask("example", "light", "", 0, "", "p-");
        for cut in 0..whole.len() {
            let error = drawn(&whole[..cut]).unwrap_err();
            assert!(error.contains("shape this expects"), "cut {cut}: {error}");
        }
    }

    /// The version the examples are linted against here is a copy of the one
    /// the library promises, so the two are compared rather than trusted: a
    /// root Cargo.toml that moved on without this would have Clippy here asking
    /// the examples for functions the library cannot use yet, or not asking for
    /// ones it can.
    #[test]
    fn the_examples_are_linted_here_against_the_version_the_library_promises() {
        let quoted = |text: &str, key: &str| -> Option<String> {
            let at = text.find(key)? + key.len();
            let rest = &text[at..];
            let open = rest.find('"')? + 1;
            let close = rest[open..].find('"')? + open;
            Some(rest[open..close].to_string())
        };
        let manifest = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../Cargo.toml"),
        )
        .unwrap();
        let promised = quoted(&manifest, "\nrust-version").expect("the library states a version");
        let linted = quoted(include_str!("lib.rs"), "#[clippy::msrv").expect("the attribute");
        assert_eq!(linted, promised);
    }

    /// The width limit is a second copy of the command line's, so the two are
    /// checked against each other rather than trusted to agree.
    #[test]
    fn the_widest_figure_is_the_widest_the_command_line_draws() {
        let argv = |width: usize| -> Vec<String> {
            ["chr1:1-200", "--coverage", "d.bg", "--width"]
                .iter()
                .map(|word| word.to_string())
                .chain([width.to_string()])
                .collect()
        };
        assert!(args::parse(&argv(MAX_WIDTH)).is_ok());
        assert!(args::parse(&argv(MAX_WIDTH + 1)).is_err());
    }
}
