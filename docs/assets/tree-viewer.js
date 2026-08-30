// A tree viewer: the crate works out where every branch goes, the page flies
// over it.
//
// The layout is computed once and never again. That is the whole idea, and it
// is what the earlier version of this page got wrong: it drew a figure, let a
// transform move it, and then threw the figure away and drew another one, so
// every gesture ended in a wait. Measured before it changed, one drawing froze
// the page for 77 ms at twenty thousand tips and 2,183 at a million.
//
// Now the program is asked once for the coordinates, in
// `karyon::cli`'s wasm bridge, and the page keeps a camera. Moving the camera
// costs a repaint of what is on screen, which is bounded by the screen and not
// by the file. Nothing is recomputed, so nothing is waited for.
//
// The figure is still the crate's. What a reader saves with Export is karyon's
// own SVG of the view they are looking at, drawn by the same program a shell
// runs, so the canvas is a way of getting about and not a second opinion about
// where the branches are.

(function () {
  "use strict";

  var K = window.karyon;
  var el = {};
  var painter = null;
  var tree = { name: "tree.nwk", body: "" };
  var sheet = null;
  var worker = null;
  // One counter for each kind of job. They shared one, so pressing Export
  // while a layout was still in flight bumped it and the layout's answer was
  // thrown away when it came back.
  // A tree that has been read but not yet drawn. Not `pending`, which this file
  // already spends on the resize timer.
  var candidate = null;
  // Whether the layout in flight is for a tree that is not on screen yet.
  var fresh = true;
  // What to put back once the tree that is being built has arrived.
  var building = null;
  // What the search box was last asked for, so asking again means the next one.
  var lastSought = null;
  var asked = 0;
  var drawn = 0;
  var cladogram = false;
  // Which of the three the reader has asked for. The first two are the same
  // layout read two ways; the third is a different walk, so asking for it
  // means asking the program again.
  var projection = "rows";
  var pending = null;

  // ------------------------------------------------------------- the program

  function here(name) {
    var src = (document.currentScript && document.currentScript.src) || "";
    if (!src) {
      var all = document.getElementsByTagName("script");
      for (var i = 0; i < all.length; i++) {
        if (/tree-viewer\.js/.test(all[i].src)) src = all[i].src;
      }
    }
    return src ? src.replace(/[^/]*$/, name) : name;
  }

  function startWorker() {
    worker = new Worker(here("tree-worker.js"));
    worker.addEventListener("message", function (event) {
      arrived(event.data);
    });
    // A worker that cannot be started at all raises this and nothing else, and
    // without it the page waited for an answer that was never coming.
    worker.addEventListener("error", function (event) {
      working(false);
      fail(
        "the program could not be started in this browser" +
          (event && event.message ? ": " + event.message : "")
      );
    });
  }

  function send(job, giving) {
    if (!worker) startWorker();
    worker.postMessage(job, giving || []);
  }

  function files() {
    var list = [{ name: tree.name, body: tree.body }];
    if (sheet) list.push({ name: sheet.name, body: sheet.body });
    return list;
  }

  function working(yes) {
    // On the plot once there is one, and on the panel before that. The sweep
    // used to live only on the plot, which is inside the panel that stays
    // hidden until a tree arrives, so the whole of the first wait showed
    // nothing at all.
    el.plot.classList.toggle("tv-working", yes);
    el.drop.classList.toggle("tv-working", yes);
  }

  function arrived(message) {
    if (message.kind === "layout" && message.id !== asked) return;
    if (message.kind === "drawn" && message.id !== drawn) return;
    working(false);

    if (message.kind === "layout") {
      if (!message.ok) {
        // What was on screen stays on screen and stays true, and the state
        // keeps holding the tree that worked.
        candidate = null;
        fail(message.body);
        return;
      }
      var arrivedFresh = fresh;
      lastSought = null;
      if (candidate) {
        tree = candidate;
        candidate = null;
        sheet = null;
        paintSheet();
      }
      // A layout that came back for the tree already on screen keeps the
      // window onto it.
      painter.load(message, !arrivedFresh);
      // The rootless walk arrives as its own layout, so which projection the
      // page is looking through is settled by which one came back.
      painter.shape(projection);
      el.drop.hidden = true;
      el.app.hidden = false;
      el.error.hidden = true;
      el.count.textContent = message.count.toLocaleString("en") + " nodes";
      repaint();
      if (building) {
        building();
        building = null;
        // The panel the button was on has just been hidden, so a keyboard is
        // otherwise left at the top of the document with no sign anything
        // happened.
        var landing = el.projection.querySelector('[aria-checked="true"]');
        if (landing) landing.focus();
      }
      // A tree with no branch lengths draws as one flat line, and the control
      // that fixes it is right there unremarked.
      var measured = false;
      for (var i = 0; i < message.count && !measured; i++) {
        if (message.x[i] !== 0) measured = true;
      }
      say();
      el.note.textContent = measured
        ? ""
        : "this file carries no branch lengths, so every tip sits at the root; Cladogram counts them instead";
      el.note.hidden = measured;
      return;
    }

    if (message.kind === "drawn") {
      if (message.ok) {
        save(message.body);
        el.error.hidden = true;
      } else {
        // Named, because the same words twice running look like nothing
        // happened at all.
        fail("that view could not be saved: " + message.body);
      }
    }
  }

  // What is under the hand, put beside it. Null clears it.
  function name(where) {
    if (!where) {
      el.tip.hidden = true;
      el.plot.classList.remove("tv-onnode");
      return;
    }
    var found = painter.at(where.x, where.y);
    if (!found) {
      el.tip.hidden = true;
      el.plot.classList.remove("tv-onnode");
      return;
    }
    var tips = found.tips === 1 ? "one tip" : found.tips.toLocaleString("en") + " tips";
    el.tip.innerHTML = "";
    var head = document.createElement("b");
    head.textContent = found.name || (found.tips === 1 ? "an unnamed tip" : "a branch");
    var under = document.createElement("small");
    under.textContent = tips + " beyond it, at depth " + found.depth.toPrecision(3);
    el.tip.appendChild(head);
    el.tip.appendChild(under);
    el.tip.hidden = false;
    el.plot.classList.add("tv-onnode");
    // Beside the pointer, and on whichever side has the room.
    var box = el.plot.getBoundingClientRect();
    var wide = el.tip.getBoundingClientRect().width;
    var left = found.x + 14;
    if (left + wide > box.width - 8) left = found.x - wide - 14;
    el.tip.style.left = Math.max(4, left) + "px";
    el.tip.style.top = Math.max(4, Math.min(box.height - 40, found.y + 12)) + "px";
  }

  function fail(text) {
    el.error.textContent = text;
    el.error.hidden = false;
  }

  // ---------------------------------------------------------------- painting

  var theme = { branch: "#1b1f23", muted: "#4b5563", font: "system-ui, sans-serif" };

  function readTheme() {
    var dark = K.dark();
    theme.branch = dark ? "#e6edf3" : "#1b1f23";
    theme.muted = dark ? "#aab4c0" : "#4b5563";
    // The rail stands behind the tree, so its ink is quieter than the tree's.
    // Quieter is not invisible: the first pair tried here measured 1.96 to 1 on
    // white, and a silhouette nobody can see is the whole point of the rail
    // thrown away. These are 3.75 and 4.92.
    // The rail and the dial are controls, so the line round them has to clear
    // the three to one a boundary needs. The pair here was 1.4, which is a line
    // nobody can find with a mouse; these measure 3.03 and 3.01.
    theme.frame = dark ? "#62666b" : "#91959a";
    theme.faint = dark ? "#79838f" : "#7b8591";
    // The mark is the one thing on the canvas that is not a branch, and it is
    // the page's own accent rather than a fourth hue: the same colour the
    // working bar uses to say which part of this is live.
    theme.window = dark ? "rgba(232, 131, 58, 0.24)" : "rgba(213, 94, 0, 0.18)";
    theme.edge = dark ? "#e8833a" : "#d55e00";
    // Something for the dial to sit on, since the disc behind it would show
    // through and the small tree would be read as part of the big one.
    // A step off the page rather than exactly it, so the dial reads as a plate
    // laid on the canvas and not as a hole cut in it.
    theme.plate = dark ? "#1e2327" : "#f4f5f7";
  }

  // Painted on the spot rather than on a frame callback: a browser does not
  // run frame callbacks in a tab it is not showing, and a viewer that comes up
  // blank until it is looked at is a viewer that looks broken.
  function repaint() {
    if (!painter || !painter.loaded()) return;
    readTheme();
    var report = painter.paint(theme);
    // A wheel notch that lands on nothing is a wheel notch that went too far,
    // and only the drawing knows: a window can sit squarely on the tree and
    // still hold the gap between two branches. One step back, once.
    if (!report.drawn && painter.stepBack()) report = painter.paint(theme);
    var at = painter.looking();
    // How many rows are on screen, and not how many nodes fall inside them: a
    // parent sits on a row of its own between its children, so counting nodes
    // said 39,999 rows of a tree with 20,000 tips. The page is written in
    // English, so the number is grouped in English wherever the browser is.
    var rows = at ? at.rows : 0;
    el.rowsOut.textContent = at
      ? rows.toLocaleString("en") + (rows === 1 ? " row" : " rows") + " in view"
      : "";
    // Below a certain width there is no room for the small picture, and the
    // hint underneath went on describing one and telling the reader to click
    // it. The paint says whether it drew one, so the page can stop claiming it.
    el.app.classList.toggle("tv-nomap", !report.rail);
    el.detail.textContent =
      report.stride > 1
        ? "1 row in " + report.stride + " and its ancestors"
        : report.labels + (report.labels === 1 ? " name" : " names");
    say();
  }

  // The command a shell would type for the view on screen, which is also what
  // Export asks the program for.
  function command() {
    var at = painter && painter.looking();
    var argv = ["tree:1-1", "--tree", tree.name];
    if (cladogram) argv = argv.concat(["--shape", "cladogram"]);
    if (projection === "disc") argv = argv.concat(["--projection", "circular"]);
    if (projection === "spread") argv = argv.concat(["--projection", "unrooted"]);
    // The rootless view holds the whole tree however far it is zoomed: there
    // are no rows to have some of. Asking for a subtree there would hand back
    // a figure of something else.
    if (projection !== "spread" && at && at.first && at.last && at.first !== at.last) {
      argv = argv.concat(["--focus", at.first + "," + at.last]);
    }
    argv = argv.concat(["--max-rows", String(Math.max(8, Math.min(2000, at ? at.rows : 60)))]);
    // What the reader searched for, marked on the figure they asked for. A
    // search that matched a hundred thousand tips is not a highlight, so only a
    // handful are carried and the rest stay a thing you look at on the page.
    var hits = painter && painter.found ? painter.found() : { count: 0, names: [] };
    if (hits.count && hits.names.length && hits.count <= 12) {
      argv = argv.concat(["--highlight", hits.names.join(",")]);
    }
    if (sheet) argv = argv.concat(["--traits", sheet.name]);
    argv.push("--no-region-label");
    return argv;
  }

  function say() {
    el.command.textContent = K.join(command());
  }

  function exportSvg() {
    drawn += 1;
    working(true);
    el.error.hidden = true;
    send({ kind: "draw", id: drawn, command: K.join(command()), room: 900 });
  }

  function save(svg) {
    var blob = new Blob([svg], { type: "image/svg+xml" });
    var url = URL.createObjectURL(blob);
    var link = document.createElement("a");
    var called = tree.name.replace(/\.[^.]*$/, "") + ".svg";
    link.href = url;
    link.download = called;
    el.note.textContent = "saved " + called;
    el.note.hidden = false;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  // --------------------------------------------------------------- the hand

  function hand() {
    var dragging = false;
    var scrubbing = false;
    // Which pointer owns the gesture. Without it a second finger overwrites the
    // first one's mode and the first release ends the gesture for both.
    var owner = null;
    var from = { x: 0, y: 0 };

    // Where a pointer is, in the canvas's own pixels. The canvas and not the
    // box around it: that box carries a one pixel border, and measuring from it
    // put every click on the rail a row of pixels low, which at two million
    // rows is three thousand of them.
    function at(event) {
      var box = el.canvas.getBoundingClientRect();
      return { x: event.clientX - box.left, y: event.clientY - box.top };
    }

    el.plot.addEventListener(
      "wheel",
      function (event) {
        if (!painter.loaded()) return;
        event.preventDefault();
        // A line-mode wheel reports a handful of lines where a pixel-mode one
        // reports tens of pixels, and treating them alike makes a mouse either
        // useless or violent next to a trackpad.
        var step = event.deltaMode === 1 ? event.deltaY * 16 : event.deltaY;
        var here = at(event);
        var box = { wide: el.canvas.clientWidth, tall: el.canvas.clientHeight };
        // On the rail a wheel does what a wheel does to a scrollbar: it moves
        // the window on from where it is rather than changing how much of the
        // tree is in it. A zoom there would be anchored on the wrong row, since
        // the two sides of the canvas are at different scales.
        if (painter.onMap(here.x, here.y)) {
          // On the rail a wheel does what it does to a scrollbar. On the dial
          // there is nothing to scroll, so it zooms the view it stands for,
          // about the middle of that view rather than about the dial.
          if (painter.shapeNow() !== "rows") {
            painter.zoomAt(box.wide / 2, box.tall / 2, Math.exp(-step * 0.002));
          }
          else painter.panBy(0, -step);
        } else painter.zoomAt(here.x, here.y, Math.exp(-step * 0.002));
        repaint();
      },
      { passive: false }
    );

    el.plot.addEventListener("pointerdown", function (event) {
      if (event.button !== 0 || !painter.loaded() || owner !== null) return;
      var here = at(event);
      owner = event.pointerId;
      // Which half of the canvas the gesture began on decides what it is for
      // the whole of its life, so a drag that starts on the rail and wanders
      // onto the tree goes on moving the window.
      if (painter.onMap(here.x, here.y)) {
        scrubbing = true;
        painter.jumpTo(here.x, here.y);
        repaint();
      } else {
        dragging = true;
        el.plot.classList.add("tv-dragging");
      }
      from = { x: event.clientX, y: event.clientY };
      el.plot.setPointerCapture(event.pointerId);
    });

    el.plot.addEventListener("pointermove", function (event) {
      if (owner !== null && event.pointerId !== owner) return;
      if (scrubbing) {
        var where = at(event);
        painter.jumpTo(where.x, where.y);
        repaint();
        return;
      }
      if (dragging) {
        painter.panBy(event.clientX - from.x, event.clientY - from.y);
        from = { x: event.clientX, y: event.clientY };
        repaint();
        return;
      }
      // Not a gesture, just a hand passing over. The cursor says which of the
      // things underneath it would answer, and a branch says what it is.
      if (painter.loaded()) {
        var over = at(event);
        var onMap = painter.onMap(over.x, over.y);
        var rows = painter.shapeNow() === "rows";
        el.plot.classList.toggle("tv-onrail", onMap && rows);
        el.plot.classList.toggle("tv-ondial", onMap && !rows);
        name(onMap ? null : over);
      }
    });

    // lostpointercapture as well as the two endings: capture can be taken away
    // without either of them firing, and that used to strand a pan. Now it
    // would strand a mode, which is worse.
    ["pointerup", "pointercancel", "lostpointercapture"].forEach(function (kind) {
      el.plot.addEventListener(kind, function (event) {
        if (owner !== null && event.pointerId !== owner) return;
        if (!dragging && !scrubbing) {
          owner = null;
          return;
        }
        dragging = false;
        scrubbing = false;
        owner = null;
        el.plot.classList.remove("tv-dragging");
        if (el.plot.hasPointerCapture(event.pointerId)) {
          el.plot.releasePointerCapture(event.pointerId);
        }
      });
    });

    // Off the canvas, so a name does not hang about over a picture the hand
    // has left.
    el.plot.addEventListener("pointerleave", function () { name(null); });

    el.plot.addEventListener("click", function (event) {
      if (!painter.loaded() || dragging || scrubbing) return;
      var where = at(event);
      if (painter.onMap(where.x, where.y)) return;
      var found = painter.at(where.x, where.y);
      if (!found) return;
      // Take the clade. On a tree of two million this is the only way to get
      // from the whole of it to one part without knowing a name to search for.
      painter.focusOn(found.node);
      repaint();
      say();
    });

    el.plot.addEventListener("dblclick", function (event) {
      var here = at(event);
      if (!painter.loaded()) return;
      if (painter.onMap(here.x, here.y)) return;
      painter.zoomAt(here.x, here.y, 2.4);
      repaint();
    });

    el.fit.addEventListener("click", function () {
      painter.home();
      repaint();
    });

    // Every gesture this viewer had was a pointer gesture, so a reader without
    // one could reach the canvas and then do nothing with it. These are the
    // same three moves the hand has: go somewhere, get closer, start again.
    el.canvas.addEventListener("keydown", function (event) {
      if (!painter.loaded()) return;
      var box = el.canvas.getBoundingClientRect();
      var step = event.shiftKey ? 0.4 : 0.12;
      // Sideways does nothing in the rectangular view, where the depth axis is
      // fitted to what is drawn on every paint, so those two keys are left for
      // the browser to scroll the page with rather than swallowed.
      var flat = painter.shapeNow() !== "rows";
      var by = { x: 0, y: 0 };
      if (event.key === "ArrowLeft") { if (!flat) return; by.x = box.width * step; }
      else if (event.key === "ArrowRight") { if (!flat) return; by.x = -box.width * step; }
      else if (event.key === "ArrowUp") by.y = box.height * step;
      else if (event.key === "ArrowDown") by.y = -box.height * step;
      else if (event.key === "PageUp") by.y = box.height * 0.9;
      else if (event.key === "PageDown") by.y = -box.height * 0.9;
      else if (event.key === "+" || event.key === "=") {
        painter.zoomAt(box.width / 2, box.height / 2, 1.6);
      } else if (event.key === "-" || event.key === "_") {
        painter.zoomAt(box.width / 2, box.height / 2, 1 / 1.6);
      } else if (event.key === "Home") {
        painter.home();
      } else {
        return;
      }
      event.preventDefault();
      if (by.x || by.y) painter.panBy(by.x, by.y);
      repaint();
    });
  }

  // ------------------------------------------------------------- the files

  // Which of the two things a dropped file is. It used to ask whether the first
  // four kilobytes held a bracket and a semicolon, so a tree it did not
  // recognise was filed as a table without a word, and the command line then
  // named the same file as both the tree and the traits. It asks the other way
  // round now: a table is a first line of fields with no bracket in it, and
  // anything else is offered to the parser, which can say what is wrong with it
  // far better than two characters can.
  function looksLikeATable(text) {
    var head = text.slice(0, 4096);
    var first = "";
    var lines = head.split("\n");
    for (var i = 0; i < lines.length && !first; i++) {
      if (lines[i].trim()) first = lines[i];
    }
    if (!first) return false;
    if (first.indexOf("(") >= 0) return false;
    return first.indexOf("\t") >= 0 || first.indexOf(",") >= 0;
  }

  function take(text, name) {
    if (!looksLikeATable(text)) {
      load(text, name);
      return;
    }
    if (!tree.body) {
      fail("that reads as a table rather than a tree; drop the phylogeny first and the sheet after it");
      return;
    }
    var called = name || "traits.tsv";
    if (called === tree.name) {
      fail("that file is already the tree; a sheet of traits has to be a second file");
      return;
    }
    sheet = { name: called, body: text };
    paintSheet();
    el.error.hidden = true;
    say();
  }

  function paintSheet() {
    if (!el.sheet) return;
    el.sheet.hidden = !sheet;
    if (sheet) el.sheetName.textContent = sheet.name;
  }

  function load(text, name) {
    // Held rather than committed. Writing the candidate into `tree` before
    // anything had read it meant a file that turned out not to be a tree still
    // counted as one: the guard that says "drop the phylogeny first" stopped
    // firing, and the next thing dropped was taken as a sheet for a tree that
    // did not exist.
    candidate = { name: name || "tree.nwk", body: text };
    send({ kind: "files", files: [{ name: candidate.name, body: candidate.body }] });
    // Through the one place that knows which layout to ask for. Asking here as
    // well left a new tree drawn with a root while the button and the exported
    // command both said it had none, and the command under the picture has to
    // be the picture.
    relayout();
  }

  function relayout() {
    var reading = candidate || tree;
    if (!reading.body) return;
    fresh = !!candidate;
    // The old refusal was about the old attempt. Leaving it up meant a second
    // failure with the same words looked like nothing had happened at all.
    el.error.hidden = true;
    asked += 1;
    working(true);
    send({
      kind: "layout",
      id: asked,
      name: reading.name,
      body: reading.body,
      cladogram: cladogram,
      rootless: projection === "spread",
    });
  }

  // ------------------------------------------------------------- examples

  function example(tips) {
    var seed = 20260829;
    function next() {
      seed = (seed * 1103515245 + 12345) % 2147483648;
      return seed / 2147483648;
    }
    var lineages = ["L1", "L2", "L3", "L4", "L4.9", "L2.2.1"];
    var parts = new Array(tips);
    for (var i = 0; i < tips; i++) {
      var lineage = lineages[Math.floor(next() * lineages.length)];
      parts[i] = lineage + "_" + String(i).padStart(5, "0") + ":" + (next() * 0.02).toFixed(4);
    }
    // Paired up a level at a time, which is linear: joining two at random out
    // of a list and putting the pair back is quadratic, and at a million tips
    // that is the difference between a moment and a locked page.
    while (parts.length > 1) {
      var up = [];
      for (var j = 0; j + 1 < parts.length; j += 2) {
        up.push("(" + parts[j] + "," + parts[j + 1] + ")" +
          (0.5 + next() * 0.5).toFixed(2) + ":" + (next() * 0.02).toFixed(4));
      }
      if (parts.length % 2) up.push(parts[parts.length - 1]);
      parts = up;
    }
    return parts[0] + ";\n";
  }

  // --------------------------------------------------------------- wiring

  function start() {
    ["plot", "canvas", "search", "rowsOut", "detail", "count", "command", "error",
      "drop", "app", "file", "paste", "usePaste", "fit", "export", "sheet",
      "sheetName", "dropSheet", "lengths", "projection", "note", "another", "tip"].forEach(function (name) {
      el[name] = document.getElementById("tv-" + name.toLowerCase());
    });
    if (!el.canvas || !el.lengths || !el.projection) return;

    painter = window.karyonCanvas.make(el.canvas);
    hand();

    el.search.addEventListener("input", function () {
      // The refusal was about what was typed before, and what is typed now is
      // different, so it stops being true the moment a key lands.
      el.error.hidden = true;
    });

    el.search.addEventListener("keydown", function (event) {
      if (event.key !== "Enter") return;
      event.preventDefault();
      var wanted = el.search.value.trim();
      if (!wanted) return;
      // The same words again means the next one of them, not the first again.
      if (wanted === lastSought) {
        var moved = painter.nextFound(event.shiftKey ? -1 : 1);
        if (moved) {
          el.note.textContent =
            moved.of.toLocaleString("en") + " tips match, at " + (moved.at + 1);
          el.note.hidden = false;
          repaint();
          return;
        }
      }
      lastSought = wanted;
      // Exactly, then ignoring case, then a prefix, then anywhere in the name.
      // A name typed from memory is rarely typed to the letter, and a reader
      // looking for a lineage or a country is not typing a name at all.
      var how = ["exact", "loose", "starts", "in"];
      var got = false;
      for (var i = 0; i < how.length && !got; i++) got = painter.goTo(wanted, how[i]);
      if (got) {
        var found = painter.found();
        el.error.hidden = true;
        el.note.hidden = false;
        el.note.textContent = found.count === 1
          ? "one tip matches; Enter again to search elsewhere"
          : found.count.toLocaleString("en") + " tips match, at 1; Enter again for the next";
        repaint();
      } else {
        fail("no tip here is called " + wanted + ", or has it in its name");
      }
    });

    window.karyonRadio.make(el.lengths, "data-lengths", function (value) {
      var wanted = value === "cladogram";
      if (wanted === cladogram) return;
      cladogram = wanted;
      relayout();
    });

    window.karyonRadio.make(el.projection, "data-projection", function (value) {
      if (value === projection) return;
      var wasRootless = projection === "spread";
      projection = value;
      // Rectangular and circular are the same rows and depths read two ways, so
      // one is a repaint. The rootless walk is a different walk and the program
      // has to do it, so going in or out of it costs a round trip.
      if (wasRootless || projection === "spread") relayout();
      else {
        painter.shape(projection);
        repaint();
        say();
      }
    });

    el.export.addEventListener("click", exportSvg);

    el.another.addEventListener("change", function () {
      var file = el.another.files[0];
      if (file) file.text().then(function (text) { take(text, file.name); });
      el.another.value = "";
    });

    el.file.addEventListener("change", function () {
      var file = el.file.files[0];
      if (file) file.text().then(function (text) { take(text, file.name); });
    });
    el.usePaste.addEventListener("click", function () {
      if (el.paste.value.trim()) take(el.paste.value.trim(), "pasted.nwk");
    });
    if (el.dropSheet) {
      el.dropSheet.addEventListener("click", function () {
        sheet = null;
        paintSheet();
        say();
      });
    }

    document.querySelectorAll("[data-example]").forEach(function (button) {
      button.addEventListener("click", function () {
        var tips = Number(button.getAttribute("data-example"));
        button.disabled = true;
        var was = button.textContent;
        button.textContent = "building " + was;
        // Put back when the tree is drawn, not when the string is built: the
        // string is the fast part, and restoring the label there made the
        // button look finished while most of the wait was still to come.
        building = function () {
          button.disabled = false;
          button.textContent = was;
        };
        setTimeout(function () { load(example(tips), tips + "-tips.nwk"); }, 0);
      });
    });

    ["dragover", "drop"].forEach(function (kind) {
      document.addEventListener(kind, function (event) { event.preventDefault(); });
    });
    document.addEventListener("drop", function (event) {
      var file = event.dataTransfer && event.dataTransfer.files[0];
      if (file) file.text().then(function (text) { take(text, file.name); });
    });

    K.onScheme(repaint);
    window.addEventListener("resize", function () {
      if (pending) clearTimeout(pending);
      pending = setTimeout(repaint, 80);
    });
    startWorker();
    el.drop.hidden = false;
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
