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

use std::cell::RefCell;
use std::collections::BTreeMap;

use karyon::cli::{args, stack};
use karyon::Tree;

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

/// Gives back a buffer [`alloc`] or [`render`] handed out.
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

/// The whole of what [`layout`] does, with the pointers already gone.
fn positions(mut input: &[u8]) -> Result<Vec<u8>, String> {
    let name = text(&mut input).ok_or("the file name is not in the shape this expects")?;
    let body = text(&mut input).ok_or("the file body is not in the shape this expects")?;
    let cladogram = number(&mut input).unwrap_or(0) == 1;
    let rootless = number(&mut input).unwrap_or(0) == 1;

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
    let bytes = body.as_bytes();
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
}
