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
  var ALONE = ["--axis", "--log", "--no-axis", "--no-region-label"];

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
    if (argv.indexOf("--width") < 0 && room > 320) {
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
    dark: dark,
    onScheme: onScheme,
    MIN_SPAN: MIN_SPAN,
    MAX_SPAN: MAX_SPAN,
  };
})();
