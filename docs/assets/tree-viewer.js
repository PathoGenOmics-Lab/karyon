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
  var asked = 0;
  var cladogram = false;
  // Which of the three the reader has asked for. The first two are the same
  // layout read two ways; the third is a different walk, so asking for it
  // means asking the program again.
  var projection = "rows";
  var PROJECTIONS = [
    { key: "rows", label: "Rectangular" },
    { key: "disc", label: "Circular" },
    { key: "spread", label: "Unrooted" },
  ];
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
    el.plot.classList.toggle("tv-working", yes);
  }

  function arrived(message) {
    if (message.id !== asked) return;
    working(false);

    if (message.kind === "layout") {
      if (!message.ok) {
        fail(message.body);
        return;
      }
      painter.load(message);
      // The rootless walk arrives as its own layout, so which projection the
      // page is looking through is settled by which one came back.
      painter.shape(projection);
      el.drop.hidden = true;
      el.app.hidden = false;
      el.error.hidden = true;
      el.count.textContent = message.count.toLocaleString() + " nodes";
      repaint();
      return;
    }

    if (message.kind === "drawn") {
      if (message.ok) {
        save(message.body);
        el.error.hidden = true;
      } else {
        fail(message.body);
      }
    }
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
    theme.frame = dark ? "#30363d" : "#d7dbe0";
    theme.faint = dark ? "#79838f" : "#7b8591";
    // The mark is the one thing on the canvas that is not a branch, and it is
    // the page's own accent rather than a fourth hue: the same colour the
    // working bar uses to say which part of this is live.
    theme.window = dark ? "rgba(232, 131, 58, 0.24)" : "rgba(213, 94, 0, 0.18)";
    theme.edge = dark ? "#e8833a" : "#d55e00";
    // Something for the dial to sit on, since the disc behind it would show
    // through and the small tree would be read as part of the big one.
    theme.plate = dark ? "#161a1d" : "#ffffff";
  }

  // Painted on the spot rather than on a frame callback: a browser does not
  // run frame callbacks in a tab it is not showing, and a viewer that comes up
  // blank until it is looked at is a viewer that looks broken.
  function repaint() {
    if (!painter || !painter.loaded()) return;
    readTheme();
    var report = painter.paint(theme);
    var at = painter.looking();
    el.rowsOut.textContent = at
      ? at.rows.toLocaleString() + (at.rows === 1 ? " row" : " rows") + " in view"
      : "";
    el.detail.textContent =
      report.stride > 1
        ? "1 row in " + report.stride + " and its ancestors"
        : report.labels + " names";
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
    if (sheet) argv = argv.concat(["--traits", sheet.name]);
    argv.push("--no-region-label");
    return argv;
  }

  function say() {
    el.command.textContent = K.join(command());
  }

  function exportSvg() {
    asked += 1;
    working(true);
    send({ kind: "draw", id: asked, command: K.join(command()), room: 900 });
  }

  function save(svg) {
    var blob = new Blob([svg], { type: "image/svg+xml" });
    var url = URL.createObjectURL(blob);
    var link = document.createElement("a");
    link.href = url;
    link.download = "tree.svg";
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
      // two things underneath it would answer.
      if (painter.loaded()) {
        var over = at(event);
        var onMap = painter.onMap(over.x, over.y);
        var rows = painter.shapeNow() === "rows";
        el.plot.classList.toggle("tv-onrail", onMap && rows);
        el.plot.classList.toggle("tv-ondial", onMap && !rows);
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
  }

  // ------------------------------------------------------------- the files

  function looksLikeATree(text) {
    var head = text.slice(0, 4096);
    return head.indexOf("(") >= 0 && head.indexOf(";") >= 0;
  }

  function take(text, name) {
    if (looksLikeATree(text)) {
      load(text, name);
      return;
    }
    if (!tree.body) {
      fail("that reads as a table rather than a tree; drop the phylogeny first and the sheet after it");
      return;
    }
    sheet = { name: name || "traits.tsv", body: text };
    paintSheet();
    say();
  }

  function paintSheet() {
    if (!el.sheet) return;
    el.sheet.hidden = !sheet;
    if (sheet) el.sheetName.textContent = sheet.name;
  }

  function load(text, name) {
    tree = { name: name || "tree.nwk", body: text };
    sheet = null;
    paintSheet();
    asked += 1;
    working(true);
    send({ kind: "files", files: files() });
    send({ kind: "layout", id: asked, name: tree.name, body: text, cladogram: cladogram });
  }

  function relayout() {
    if (!tree.body) return;
    asked += 1;
    working(true);
    send({
      kind: "layout",
      id: asked,
      name: tree.name,
      body: tree.body,
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
      "sheetName", "dropSheet", "shape", "round"].forEach(function (name) {
      el[name] = document.getElementById("tv-" + name.toLowerCase());
    });
    if (!el.canvas) return;

    painter = window.karyonCanvas.make(el.canvas);
    hand();

    el.search.addEventListener("keydown", function (event) {
      if (event.key !== "Enter") return;
      event.preventDefault();
      var wanted = el.search.value.trim();
      if (!wanted) return;
      if (painter.goTo(wanted)) {
        el.error.hidden = true;
        repaint();
      } else {
        fail("no tip here is called " + wanted);
      }
    });

    el.round.addEventListener("click", function () {
      var at = 0;
      for (var i = 0; i < PROJECTIONS.length; i++) {
        if (PROJECTIONS[i].key === projection) at = i;
      }
      var next = PROJECTIONS[(at + 1) % PROJECTIONS.length];
      var wasRootless = projection === "spread";
      projection = next.key;
      el.round.setAttribute("aria-pressed", String(projection !== "rows"));
      el.round.textContent = next.label;
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

    el.shape.addEventListener("click", function () {
      cladogram = !cladogram;
      el.shape.setAttribute("aria-pressed", String(cladogram));
      el.shape.textContent = cladogram ? "Cladogram" : "Phylogram";
      relayout();
    });

    el.export.addEventListener("click", exportSvg);

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
        setTimeout(function () {
          load(example(tips), tips + "-tips.nwk");
          button.disabled = false;
          button.textContent = was;
        }, 0);
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
