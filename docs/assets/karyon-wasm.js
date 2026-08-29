// The bridge to the program, shared by every page that runs it.
//
// The protocol is written down in playground/src/lib.rs and repeated here only
// where the code would otherwise be a row of magic offsets: one buffer in and
// one buffer out, every number a little-endian u32, every string UTF-8.
//
// One copy, because there are two pages that run the program and two copies of
// a protocol are two things that drift. No framework and no build step: the
// crate has no dependencies, and a page that needed a bundler to demonstrate a
// program that needs nothing would be making the wrong point.

window.karyon = (function () {
  "use strict";

  // Worked out from where this script itself was loaded from, rather than
  // written down. The site is served from a sub-path on GitHub Pages and from
  // the root on a laptop, and a page under a folder resolves a bare relative
  // path against the folder rather than against the site.
  var here = (document.currentScript && document.currentScript.src) || "";
  var WASM = here
    ? here.replace(/[^/]*$/, "karyon_playground.wasm")
    : "karyon_playground.wasm";

  var wasm = null;
  var arriving = null;
  var encoder = new TextEncoder();
  var decoder = new TextDecoder();

  // Below sixty bases the ruler has nothing left to label.
  var MIN_SPAN = 60;
  // And above this the parser refuses the region, because a per-base track
  // keeps eight bytes a base and a `Vec` on a 32-bit target cannot address
  // more. It is the same number as `MAX_SPAN` in src/cli/args.rs and has to
  // stay so: a wheel that can reach a span the program will not draw is a
  // wheel that stops working partway round.
  var MAX_SPAN = 268435455;

  function load() {
    if (arriving) return arriving;
    arriving = fetch(WASM)
      .then(function (response) {
        if (!response.ok) throw new Error(response.status + " fetching the program");
        return response.arrayBuffer();
      })
      // `instantiate` rather than `instantiateStreaming`, which needs the file
      // to arrive as application/wasm and fails outright when a host serves it
      // as something else.
      .then(function (bytes) { return WebAssembly.instantiate(bytes, {}); })
      .then(function (built) { wasm = built.instance.exports; return wasm; });
    return arriving;
  }

  function ready() {
    return wasm !== null;
  }

  // ---------------------------------------------------------------- command

  // A command line, not a shell. Quoting is handled because a title has spaces
  // in it and every example in the documentation shows one; everything else a
  // shell does is a shell's business and is not reimplemented here.
  function words(text) {
    var out = [];
    var word = "";
    var quote = null;
    var started = false;
    for (var i = 0; i < text.length; i++) {
      var c = text[i];
      if (quote) {
        if (c === quote) quote = null;
        else { word += c; started = true; }
      } else if (c === '"' || c === "'") {
        quote = c;
        started = true;
      } else if (c === "\\" && text[i + 1] === "\n") {
        i++;
      } else if (/\s/.test(c)) {
        if (started) { out.push(word); word = ""; started = false; }
      } else {
        word += c;
        started = true;
      }
    }
    if (started) out.push(word);
    return out;
  }

  var LOCUS = /^(.+):([\d,]+)-([\d,]+)$/;

  // The four flags that stand on their own. Every other flag takes the word
  // after it, and that word is not the region however much it looks like one:
  // `--label 'chr1:5-9'` written before the locus was read as the locus, and a
  // drag then rewrote the label and left the figure where it was.
  var ALONE = ["--axis", "--fade-by-mapq", "--log", "--no-axis", "--no-counts",
               "--no-names", "--no-region-label"];

  // The flags that open a track. Everything between one of these and the next
  // describes that track, which is what makes a command a stack rather than a
  // set of settings, and it is the only way to tell one track's `--style` from
  // another's.
  //
  // Written out rather than worked out, and the reason is worth having. The
  // structural test is that a track flag is a `--word` followed by a file, and
  // this page knows its own file names, but `--traits`, `--with-tree`,
  // `--against`, `--links` and `--with-sequence` all take a file and describe
  // the track before them rather than opening one. So that reading needs a
  // list of its own exceptions, which is a second copy of a fact exactly as
  // this is, and a worse one: when it falls behind, a modifier is mistaken for
  // a track, a segment ends early, and a control quietly reads a flag off the
  // wrong track. When this list falls behind, a Rust test fails by name.
  //
  // That test is `the_playground_knows_every_track_the_parser_does` in
  // src/cli/args.rs, and it reads this file and compares both directions, the
  // way the help text and `Kind::ALL` are already held together.
  var TRACKS = [
    "--coverage", "--copy-number", "--dynseq", "--junctions", "--sequence",
    "--features", "--variants", "--windows", "--manhattan", "--tree",
    "--msa", "--snps", "--ideogram", "--matrix", "--pileup",
    "--synteny", "--dotplot", "--orfs", "--logo", "--tanglegram",
    "--clades", "--loci", "--methylation", "--structural", "--split-reads",
    "--bisulfite", "--domains", "--axis",
  ];

  /// Where the words describing one track begin and end.
  ///
  /// `after` names the track a modifier belongs to, and the modifier can only
  /// be the one written between that track's flag and the next track's. Two
  /// tracks in a command may carry the same modifier with different values,
  /// and without this every helper below found whichever came first.
  ///
  /// `null` where no track was named, which is what a figure option passes,
  /// and then the whole command is the segment.
  ///
  /// A track that was named and is not in the command is an empty segment
  /// rather than the whole of it. The reader can delete a track by hand while
  /// its control is still on the page, and answering with some other track's
  /// value would be the same mistake this function exists to stop.
  function segment(argv, after) {
    if (!after) return null;
    var from = argv.indexOf(after);
    if (from < 0) return [0, 0];
    var to = argv.length;
    for (var i = from + 1; i < argv.length; i++) {
      if (TRACKS.indexOf(argv[i]) >= 0) {
        to = i;
        break;
      }
    }
    return [from, to];
  }

  // The region a command names, which is its one positional word.
  function locus(text) {
    var argv = words(text);
    for (var i = 0; i < argv.length; i++) {
      var word = argv[i];
      if (word.charAt(0) === "-" && word !== "-") {
        if (ALONE.indexOf(word) < 0) i++;
        continue;
      }
      var m = LOCUS.exec(word);
      if (m) {
        return {
          word: word,
          seq: m[1],
          // 1-based inclusive, which is what a person reads and types.
          start: parseInt(m[2].replace(/,/g, ""), 10),
          end: parseInt(m[3].replace(/,/g, ""), 10),
        };
      }
    }
    return null;
  }

  function grouped(n) {
    return String(n).replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }

  // Rewrites the region a command names, leaving the rest of it alone.
  //
  // The span asked for is kept and the window is slid, rather than the start
  // being clamped where it is and the end left where it was. Clamping alone
  // made a drag towards the first base of a sequence shrink the window: from
  // `chr1:1-1,000` five drags gave 700, 490, 343, 240 and 168 bases, so the
  // figure zoomed itself in while the reader was only moving it sideways.
  function retarget(text, start, end) {
    var where = locus(text);
    if (!where) return text;
    start = Math.round(start);
    end = Math.round(end);
    var span = Math.min(MAX_SPAN, Math.max(MIN_SPAN, end - start + 1));
    if (start < 1) start = 1;
    end = start + span - 1;
    return text.replace(where.word, where.seq + ":" + grouped(start) + "-" + grouped(end));
  }

  // Moves the window by a fraction of its own span, which is what a drag of
  // that fraction of the figure's width comes to.
  function panned(text, fraction) {
    var where = locus(text);
    if (!where) return text;
    var bases = (where.end - where.start + 1) * fraction;
    return retarget(text, where.start - bases, where.end - bases);
  }

  // Multiplies the span, holding `at` (0 to 1 across the figure) still, so
  // zooming reads as magnifying what is under the pointer rather than as
  // jumping somewhere near it.
  function zoomed(text, factor, at) {
    var where = locus(text);
    if (!where) return text;
    var span = where.end - where.start + 1;
    var next = Math.min(MAX_SPAN, Math.max(MIN_SPAN, span * factor));
    var anchor = where.start + at * span;
    return retarget(text, anchor - at * next, anchor - at * next + next - 1);
  }

  // Clamps a region to the stretch a set of files actually covers.
  //
  // A window taken past the data comes back as the program's own refusal,
  // which is right when someone typed it and wrong when they were turning a
  // control: a control that can be put somewhere there is nothing to see is a
  // control that looks broken. So a control that knows its data says so, and
  // this keeps the window inside it.
  function within(text, bounds) {
    if (!bounds) return text;
    var where = locus(text);
    if (!where) return text;
    var room = bounds.to - bounds.from + 1;
    var span = Math.min(where.end - where.start + 1, room);
    // The floor is the example's own, measured against its own files: the
    // narrowest window that still draws wherever it is put. A track refuses a
    // window with nothing of its own in it, which is right, and a control that
    // can be turned to a window like that is a control that looks broken.
    if (bounds.min) span = Math.min(room, Math.max(span, bounds.min));
    var start = where.start;
    if (start < bounds.from) start = bounds.from;
    if (start + span - 1 > bounds.to) start = bounds.to - span + 1;
    return retarget(text, start, start + span - 1);
  }

  // ------------------------------------------------------------------ flags
  //
  // A control is a flag with a value, so setting one is rewriting one word of
  // the command and leaving every other word alone. The command stays the
  // thing that decides, and a reader watching it can see what a control did.

  /// The value a flag was given on the track `after` names, or null.
  function flagOf(text, flag, after) {
    var argv = words(text);
    var span = segment(argv, after) || [0, argv.length];
    for (var i = span[0]; i < span[1]; i++) {
      if (argv[i] === flag) return i + 1 < argv.length ? argv[i + 1] : "";
    }
    return null;
  }

  /// Whether the track `after` names carries this flag.
  function hasFlag(text, flag, after) {
    var argv = words(text);
    var span = segment(argv, after) || [0, argv.length];
    return argv.slice(span[0], span[1]).indexOf(flag) >= 0;
  }

  /// Sets, replaces or removes a flag and its value.
  ///
  /// `value` of `null` takes the flag out, `true` puts it in on its own, and a
  /// string puts it in with that after it. A flag that describes a track is
  /// written next to the one it describes rather than at the end, since where
  /// a word sits is what binds it in this grammar.
  function setFlag(text, flag, value, after) {
    var argv = words(text);
    var alone = ALONE.indexOf(flag) >= 0;

    // Taken out of the track it belongs to rather than out of the command.
    // Removing the first one in the command deleted a coverage track's
    // `--style area` on the way to setting a variants track's `--style tick`,
    // and left both controls reading the same word back.
    var span = segment(argv, after) || [0, argv.length];
    var at = argv.indexOf(flag, span[0]);
    if (at >= 0 && at < span[1]) argv.splice(at, alone ? 1 : 2);
    if (value === null || value === false) return join(argv);

    var piece = alone || value === true ? [flag] : [flag, String(value)];

    // A modifier describes the track before it, so where it goes is not a
    // detail: appended to the end, `--aggregate` landed on whichever track
    // happened to be last and the program refused it, correctly, as meaning
    // nothing to a variants track.
    var anchor = after ? argv.indexOf(after) : -1;
    if (anchor < 0) return join(argv.concat(piece));

    var into = anchor + (ALONE.indexOf(after) >= 0 ? 1 : 2);
    argv.splice.apply(argv, [into, 0].concat(piece));
    return join(argv);
  }

  /// Puts a command back together, quoting only what has to be quoted.
  function join(argv) {
    return argv
      .map(function (word) {
        return /\s/.test(word) ? "'" + word + "'" : word;
      })
      .join(" ");
  }

  // Material writes the palette onto `body`, and a figure drawn light on a dark
  // page is not a figure with a light theme, it is a hole in the page.
  function dark() {
    return document.body.getAttribute("data-md-color-scheme") === "slate";
  }

  // ----------------------------------------------------------------- buffers

  function pack(argv, list) {
    var parts = [];
    var total = 0;
    function u32(n) {
      var b = new Uint8Array(4);
      new DataView(b.buffer).setUint32(0, n, true);
      parts.push(b);
      total += 4;
    }
    function str(s) {
      var bytes = encoder.encode(s);
      u32(bytes.length);
      parts.push(bytes);
      total += bytes.length;
    }
    u32(argv.length);
    argv.forEach(str);
    u32(list.length);
    list.forEach(function (file) { str(file.name); str(file.body); });

    var out = new Uint8Array(total);
    var at = 0;
    parts.forEach(function (part) { out.set(part, at); at += part.length; });
    return out;
  }

  // Runs one command line over a list of `{name, body}` files.
  //
  // Answers `{ok, body, ms}`: a figure, or the sentence the command would print
  // at a shell. A refusal is an ordinary outcome and not an exception.
  function run(command, list, room) {
    if (!wasm) return { ok: false, body: "the program has not arrived", ms: 0 };
    var argv = words(command);
    // Two things a page knows and the command does not: how wide the figure has
    // room to be, and which way the page is running. Both are supplied only
    // where the command did not say, so writing either one overrules the page
    // rather than fighting it.
    // A floor of 320 was above every phone, so the one reader who most needed
    // the figure drawn to fit got the 900 pixel default squashed by CSS to
    // about a third of it instead. Nothing here has to guard the low end:
    // `Figure::width` raises a width that would leave no plotting area to the
    // smallest one that does, so the only value worth refusing is a pane that
    // has not been laid out yet.
    if (argv.indexOf("--width") < 0 && room > 0) {
      argv = argv.concat(["--width", String(Math.round(room))]);
    }
    if (argv.indexOf("--theme") < 0 && dark()) {
      argv = argv.concat(["--theme", "dark"]);
    }

    var input = pack(argv, list);
    var started = performance.now();
    var inPtr = wasm.alloc(input.length);
    new Uint8Array(wasm.memory.buffer, inPtr, input.length).set(input);

    // A trap is not an error the program returned, it is the program stopping,
    // and it escaped to nobody: both callers had already written the new
    // command down before drawing, so the page went on showing the previous
    // figure, the previous region and the previous timing underneath a region
    // it had never drawn. Caught here, it becomes an ordinary refusal, which
    // both callers already know how to show.
    var outPtr;
    try {
      outPtr = wasm.render(inPtr, input.length);
    } catch (error) {
      wasm.dealloc(inPtr, input.length);
      return {
        ok: false,
        body: "the program stopped on this command (" + error.message + ")",
        ms: performance.now() - started,
      };
    }
    wasm.dealloc(inPtr, input.length);

    // `memory.buffer` is detached and replaced whenever wasm grows its heap,
    // so it is read again here rather than reused from above.
    var view = new DataView(wasm.memory.buffer);
    var ok = view.getUint8(outPtr) === 1;
    var len = view.getUint32(outPtr + 1, true);
    var body = decoder.decode(new Uint8Array(wasm.memory.buffer, outPtr + 5, len));
    wasm.dealloc(outPtr, 5 + len);
    return { ok: ok, body: body, ms: performance.now() - started };
  }

  // Calls `then` whenever the reader changes the site's light or dark setting,
  // which Material does by rewriting an attribute rather than by reloading.
  function onScheme(then) {
    new MutationObserver(then).observe(document.body, {
      attributes: true,
      attributeFilter: ["data-md-color-scheme"],
    });
  }

  return {
    load: load,
    ready: ready,
    run: run,
    words: words,
    locus: locus,
    grouped: grouped,
    retarget: retarget,
    panned: panned,
    zoomed: zoomed,
    within: within,
    flagOf: flagOf,
    hasFlag: hasFlag,
    setFlag: setFlag,
    join: join,
    dark: dark,
    onScheme: onScheme,
    MIN_SPAN: MIN_SPAN,
    MAX_SPAN: MAX_SPAN,
  };
})();
