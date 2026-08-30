// A tree viewer: the same program the command line runs, run again on every
// move.
//
// A phylogeny of a million tips has no picture. What it has is a picture of
// the part you are looking at, and the way to look at another part is to draw
// that one instead. So nothing here draws: it works out the command line the
// figure would be drawn by and hands it to karyon, which folds the clades that
// do not fit and answers with an SVG of a few hundred rows whatever the tree
// holds. Panning is a new command. Zooming is a new command. Opening a clade
// is a new command with --focus on it.
//
// That is why it does not fall over on a big tree: the drawing is bounded by
// the rows asked for and not by the file. It is also why there is no second
// copy of the layout here to disagree with the first one.

(function () {
  "use strict";

  // The protocol and the wasm loading live in karyon-wasm.js, which the
  // playground and the home page run over too.
  var K = window.karyon;

  var el = {};
  // Where we are: a stack of clades opened, innermost last. Each entry is what
  // --focus takes and what to show in the trail.
  var trail = [];
  var tree = { name: "tree.nwk", body: "" };
  // A sheet of metadata, drawn as strips beside the tips or rings around them.
  var sheet = null;
  var rows = 400;
  var projection = "rectangular";
  var pending = null;
  // Where the hand has put the picture since it was last drawn: a scale about
  // the top left of the stage and an offset in page pixels. It is a CSS
  // transform and nothing else, so a wheel notch costs a composite and not a
  // render of the tree.
  var view = { k: 1, tx: 0, ty: 0 };
  var settling = null;
  // The camera's range. The floor is whatever scale makes the whole drawing
  // fit the window, worked out after each drawing rather than fixed, and the
  // ceiling is well past life size so a reader can get inside a cherry.
  var fit = 1;
  var ZOOM_MAX = 8;

  // The most rows worth putting in the page at once.
  //
  // Drawing is charged for the walk over the tree and barely at all for the
  // rows it ends up writing: on a two hundred thousand tip tree, sixty rows
  // took about 500 ms and two thousand took 502. What two thousand does cost is
  // the page itself, 16,005 elements and 87 ms to put them in, and five
  // thousand costs 323 ms, which is where it stops being worth it.
  //
  // This is the whole of why the viewer is not a slideshow. All of it is in the
  // page, the camera moves over rows that are already there, and the program is
  // asked for nothing until the reader has used that detail up.
  var DETAIL_CAP = 2000;
  // How hard the reader has pulled against the far end of the zoom, and when
  // that last took them up a level.
  var pull = 0;
  var lastPull = 0;
  var lastPop = 0;

  // ------------------------------------------------------------ the command

  // The command line this view would be drawn by, which is also the one to
  // copy into a terminal and get the same figure.
  // How many rows the window would show at the size they are drawn, which is
  // what the reader's own setting means, and how many are actually drawn.
  function drawnRows() {
    return Math.min(DETAIL_CAP, Math.max(8, Math.round(rows)));
  }

  function command(steps) {
    var path = steps || trail;
    var argv = ["tree:1-1", "--tree", tree.name, "--projection", projection];
    argv = argv.concat(["--max-rows", String(drawnRows())]);
    var at = path[path.length - 1];
    if (at) argv = argv.concat(["--focus", at.focus]);
    if (sheet) argv = argv.concat(["--traits", sheet.name]);
    argv.push("--no-region-label");
    return argv;
  }

  // What the program is given to read. A sheet is a second file and not a
  // second program: the same command line a shell would type names both.
  function files() {
    var list = [{ name: tree.name, body: tree.body }];
    if (sheet) list.push({ name: sheet.name, body: sheet.body });
    return list;
  }

  // ------------------------------------------------------------- the program

  // It runs in a worker. Drawing is not free: a million tip figure takes a
  // couple of seconds, and on this thread those are seconds in which the wheel
  // does not turn and the drag does not land. Measured before it moved: one
  // drawing froze the page for 77 ms at twenty thousand tips, 346 at two
  // hundred thousand and 2,183 at a million, and a zoom asks for up to three.
  //
  // The files go over once and stay there, because posting megabytes of Newick
  // with every request puts the copying back on this thread, which is the cost
  // being removed.
  var worker = null;
  var asked = 0;
  var showing = 0;

  function start_worker() {
    worker = new Worker(here("tree-worker.js"));
    worker.addEventListener("message", function (event) {
      arrived(event.data);
    });
  }

  // The worker sits beside this script, wherever the site is served from.
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

  function send(job) {
    if (!worker) start_worker();
    worker.postMessage(job);
  }

  function sendFiles() {
    send({ kind: "files", files: files() });
  }

  // An answer that is not the newest is thrown away. A drag can outrun a
  // drawing, and painting an old one over a new one is worse than waiting.
  function arrived(message) {
    if (message.id !== asked) return;
    working(false);

    if (message.kind === "drawn") {
      if (message.ok) {
        el.stage.innerHTML = message.body;
        reset();
        el.error.hidden = true;
        showing = message.id;
        el.timing.textContent = Math.round(message.ms) + " ms";
      } else {
        el.error.textContent = message.body;
        el.error.hidden = false;
        if (message.onFailure === "pop") {
          trail.pop();
          paintTrail();
        }
      }
      return;
    }

    if (message.kind === "chosen") {
      if (!message.focus) return;
      trail.push({ focus: message.focus, label: message.label });
      el.stage.innerHTML = message.body;
      reset();
      el.error.hidden = true;
      el.timing.textContent = Math.round(message.ms) + " ms";
      // The line under the figure is what a shell would type to get it, so it
      // has to be the line that drew what is on screen. A zoom that chose its
      // own clade used to leave the previous one written there.
      el.command.textContent = K.join(command());
      paintTrail();
    }
  }

  function working(yes) {
    el.plot.classList.toggle("tv-working", yes);
  }

  function draw() {
    if (!tree.body) return;
    var argv = command();
    var room = el.plot.clientWidth - 24;
    if (room < 40) room = 900;
    asked += 1;
    working(true);
    send({ kind: "draw", id: asked, command: K.join(argv), room: room });
    el.command.textContent = K.join(argv);
    paintTrail();
  }

  // Redrawing on every notch of a slider would queue commands faster than they
  // run on a big tree, so a burst of moves collapses into one drawing.
  //
  // A timer and not requestAnimationFrame, which a browser does not run at all
  // while its tab is in the background: the viewer came up blank in a tab that
  // had not been looked at yet, and stayed blank until it was.
  function later() {
    if (pending) return;
    pending = setTimeout(function () {
      pending = null;
      draw();
    }, 0);
  }

  // ------------------------------------------------------------- the trail

  function paintTrail() {
    el.trail.textContent = "";
    var whole = document.createElement("button");
    whole.type = "button";
    whole.className = "tv-crumb";
    whole.textContent = "whole tree";
    whole.addEventListener("click", function () {
      trail = [];
      later();
    });
    el.trail.appendChild(whole);
    trail.forEach(function (step, index) {
      var crumb = document.createElement("button");
      crumb.type = "button";
      crumb.className = "tv-crumb";
      crumb.textContent = step.label;
      crumb.addEventListener("click", function () {
        trail = trail.slice(0, index + 1);
        later();
      });
      el.trail.appendChild(crumb);
    });
  }

  // ----------------------------------------------------------- the viewport

  // A map, not a slideshow. The wheel and the drag move a CSS transform on the
  // stage, which costs a composite and never a render, so the picture keeps up
  // with the hand on a million tip tree exactly as it does on a hundred. What
  // the transform cannot do is add detail: past a certain magnification the
  // names are simply the same names drawn larger. So when the hand stops, the
  // view is read back, the clade under it is worked out, and the program is
  // asked for that clade at full detail. Zoom for the feel, redraw for the
  // detail, which is how a tiled map works and for the same reason.

  function apply() {
    el.stage.style.transform =
      "translate(" + view.tx + "px," + view.ty + "px) scale(" + view.k + ")";
    var times = view.k / fit;
    el.zoom.textContent = times > 1.02 ? "\u00d7" + (times < 10 ? times.toFixed(1) : Math.round(times)) : "";
    el.fit.disabled = times <= 1.001 && Math.abs(view.tx) < 1 && Math.abs(view.ty) < 1;
  }

  // Back to seeing all of it. The whole drawing is more than the window holds,
  // so this is a scale and not an identity.
  function reset() {
    var box = el.plot.getBoundingClientRect();
    var tall = el.stage.offsetHeight || 1;
    var wide = el.stage.offsetWidth || 1;
    fit = Math.min(1, box.height / tall, box.width / wide);
    if (!isFinite(fit) || fit <= 0) fit = 1;
    view = { k: fit, tx: 0, ty: 0 };
    contain();
    apply();
  }

  // Keeps the picture from being thrown off the edge of its own window. At a
  // scale of one it cannot move at all, which is what makes the wheel feel
  // anchored rather than slippery.
  function contain() {
    var box = el.plot.getBoundingClientRect();
    var w = el.stage.offsetWidth * view.k;
    var h = el.stage.offsetHeight * view.k;
    var slack = 40;
    view.tx = Math.min(slack, Math.max(box.width - w - slack, view.tx));
    view.ty = Math.min(slack, Math.max(box.height - h - slack, view.ty));
    if (w <= box.width) view.tx = 0;
    if (h <= box.height) view.ty = 0;
  }

  function zoomAt(px, py, factor) {
    var next = Math.min(ZOOM_MAX, Math.max(fit, view.k * factor));
    if (next === view.k) {
      // Already as far out as the picture goes, so the gesture means the level
      // above rather than nothing at all. It has to be pulled against, though:
      // one notch of a wheel is a twitch, and treating each notch as a level
      // took a view of five hundred tips back to all twenty thousand in a
      // single flick, skipping every level in between.
      if (factor < 1 && trail.length) {
        var now = Date.now();
        // Reset when the notches stop coming, not when a level was last left:
        // measuring the gap from the last pop meant the count restarted on
        // every notch until a pop happened, so it never reached one and the
        // view would not come back out at all.
        if (now - lastPull > 400) pull = 0;
        lastPull = now;
        pull += 1 - factor;
        if (pull > 0.6 && now - lastPop > 350) {
          pull = 0;
          lastPop = now;
          trail.pop();
          later();
        }
      }
      return;
    }
    pull = 0;
    var ratio = next / view.k;
    view.tx = px - (px - view.tx) * ratio;
    view.ty = py - (py - view.ty) * ratio;
    view.k = next;
    contain();
    apply();
    settle();
  }

  // --------------------------------------------------------- reading it back

  // Which names the view is looking at. Measured from the laid-out page rather
  // than from the SVG's own coordinates, so it costs nothing to be right about
  // the transform.
  //
  // Not simply "is the label on screen", which was the first try and was
  // wrong: the names are written down the right hand edge, so zooming into the
  // branches on the left put every one of them off screen and the view had
  // nothing to name. A row is what is being looked at, and a row reaches
  // across the whole figure. On a rectangle that means the height of the label
  // and nothing about its x. On a disc a row is a spoke, so the test is
  // whether the line from the middle of the drawing out to the name passes
  // through the window.
  function visibleTips() {
    var box = el.plot.getBoundingClientRect();
    var stage = el.stage.getBoundingClientRect();
    var mid = { x: (stage.left + stage.right) / 2, y: (stage.top + stage.bottom) / 2 };
    var radial = projection !== "rectangular";
    var seen = [];
    var all = el.stage.querySelectorAll("text");
    for (var i = 0; i < all.length; i++) {
      var body = all[i].textContent.trim();
      // A scale bar prints a number and nothing else; everything else written
      // beside a row is a name, whether the row is one tip or a folded clade
      // saying which tip it starts from.
      if (!body || /^[\d.,]+$/.test(body)) continue;
      var r = all[i].getBoundingClientRect();
      var on;
      if (!radial) {
        on = r.bottom > box.top && r.top < box.bottom;
      } else {
        on = false;
        var px = (r.left + r.right) / 2;
        var py = (r.top + r.bottom) / 2;
        for (var t = 0; t <= 1.0001 && !on; t += 0.1) {
          var x = mid.x + (px - mid.x) * t;
          var y = mid.y + (py - mid.y) * t;
          on = x > box.left && x < box.right && y > box.top && y < box.bottom;
        }
      }
      seen.push({ name: body.split(" +")[0].split(" (")[0], on: on });
    }
    return seen;
  }

  // The longest unbroken run of rows the view is looking at.
  //
  // It is not allowed to wrap, and that is deliberate rather than an
  // oversight. On a disc the names go round, so a wedge near where the ring
  // was started sees an arc whose two ends are the first and last tips in the
  // file. Those two are as far apart in the tree as two tips can be and the
  // smallest clade holding both of them is the root, so asking for it draws
  // the whole tree again: zooming a corner of a twenty thousand tip circle
  // asked for L1_17152,L3_00512 and got all twenty thousand back. Splitting
  // the arc at the seam and keeping the longer half asks for a clade that is
  // really there.
  function longestRun(tips) {
    var best = { at: -1, len: 0 };
    var start = -1;
    for (var i = 0; i <= tips.length; i++) {
      var on = i < tips.length && tips[i].on;
      if (on && start < 0) start = i;
      if (!on && start >= 0) {
        if (i - start > best.len) best = { at: start, len: i - start };
        start = -1;
      }
    }
    return best.len > 0 ? tips.slice(best.at, best.at + best.len) : [];
  }

  // After the hand stops: if the view has closed in on part of the picture,
  // ask the program for that part instead of magnifying what is already drawn.
  function settle() {
    if (settling) clearTimeout(settling);
    settling = setTimeout(refine, 220);
  }

  // The windows worth asking about, widest first.
  //
  // Two rows next to each other on the page can be as far apart in the tree as
  // two rows get: the clade holding a run that straddles its own deepest split
  // is the clade you are already looking at, and asking for it draws the same
  // figure again. That is what made a disc stop resolving after two gestures.
  // So when the whole of the view names nothing smaller, the middle of it is
  // tried, then the middle of that. The reader is looking at the middle.
  function windows(run) {
    var out = [];
    var seen = {};
    [run.length, Math.ceil(run.length / 2), Math.ceil(run.length / 4)].forEach(function (want) {
      var wide = Math.max(1, Math.min(run.length, want));
      var at = Math.floor((run.length - wide) / 2);
      var slice = run.slice(at, at + wide);
      var key = slice[0].name + "," + slice[slice.length - 1].name;
      if (seen[key]) return;
      seen[key] = 1;
      out.push(slice);
    });
    return out;
  }

  function refine() {
    settling = null;
    // Only once the detail already in the page has been used up. The point of
    // drawing eight times what the window holds is that most gestures need
    // nothing from the program: a reader who has zoomed twice is still looking
    // at rows that are already drawn, and asking for them again would replace
    // a picture with the same picture.
    if (view.k < fit * 4) return;
    var tips = visibleTips();
    var on = longestRun(tips);
    // Nothing to name when none of it is on screen, and nothing to gain until
    // the view has closed right in on the rows that are drawn.
    if (!on.length || on.length > tips.length * 0.12) return;

    var here_now = trail.length ? trail[trail.length - 1].focus : null;
    var candidates = [];
    windows(on).forEach(function (slice) {
      var first = slice[0].name;
      var last = slice[slice.length - 1].name;
      var focus = first === last ? first : first + "," + last;
      if (focus === here_now) return;
      candidates.push({
        focus: focus,
        label: first === last ? first : first + " to " + last,
        command: K.join(command(trail.concat([{ focus: focus, label: "" }]))),
      });
    });
    if (!candidates.length) {
      reset();
      return;
    }

    // The choosing happens in the worker. It is up to three drawings, and the
    // clade holding two rows can be the clade already on screen, so a step that
    // drew exactly what was there would waste the gesture.
    var room = el.plot.clientWidth - 24;
    if (room < 40) room = 900;
    asked += 1;
    working(true);
    send({
      kind: "choose",
      id: asked,
      candidates: candidates,
      was: K.tipsAccountedFor(el.stage.innerHTML),
      room: room,
    });
  }

  // --------------------------------------------------------------- the hand  // --------------------------------------------------------------- the hand

  function hand() {
    var dragging = false;
    var moved = 0;
    var from = { x: 0, y: 0 };

    el.plot.addEventListener(
      "wheel",
      function (event) {
        if (!el.stage.firstChild) return;
        event.preventDefault();
        var box = el.plot.getBoundingClientRect();
        // A line-mode wheel reports a handful of lines where a pixel-mode one
        // reports tens of pixels, and treating them the same makes a mouse
        // either useless or violent next to a trackpad.
        var step = event.deltaMode === 1 ? event.deltaY * 16 : event.deltaY;
        zoomAt(event.clientX - box.left, event.clientY - box.top, Math.exp(-step * 0.002));
      },
      { passive: false }
    );

    el.plot.addEventListener("pointerdown", function (event) {
      if (event.button !== 0 || !el.stage.firstChild) return;
      dragging = true;
      moved = 0;
      from = { x: event.clientX, y: event.clientY };
      el.plot.setPointerCapture(event.pointerId);
      el.plot.classList.add("tv-dragging");
    });

    el.plot.addEventListener("pointermove", function (event) {
      if (!dragging) return;
      var dx = event.clientX - from.x;
      var dy = event.clientY - from.y;
      moved += Math.abs(dx) + Math.abs(dy);
      from = { x: event.clientX, y: event.clientY };
      view.tx += dx;
      view.ty += dy;
      contain();
      apply();
    });

    ["pointerup", "pointercancel"].forEach(function (kind) {
      el.plot.addEventListener(kind, function (event) {
        if (!dragging) return;
        dragging = false;
        el.plot.classList.remove("tv-dragging");
        if (el.plot.hasPointerCapture(event.pointerId)) {
          el.plot.releasePointerCapture(event.pointerId);
        }
        if (moved > 4) settle();
      });
    });

    el.plot.addEventListener("dblclick", function (event) {
      var box = el.plot.getBoundingClientRect();
      zoomAt(event.clientX - box.left, event.clientY - box.top, 2.2);
    });

    el.fit.addEventListener("click", reset);
    if (el.dropSheet) {
      el.dropSheet.addEventListener("click", function () {
        sheet = null;
        paintSheet();
        sendFiles();
        later();
      });
    }

    // A drag that happened to end on a triangle is a drag, not a click on it.
    el.plot.addEventListener(
      "click",
      function (event) {
        if (moved > 4) {
          event.stopPropagation();
          moved = 0;
        }
      },
      true
    );
  }

  // --------------------------------------------------------------- opening

  // A folded triangle's tooltip reads "clade (46 tips), t123 to t168", and the
  // pair at the end is what names it: the tips under a node are a run, so the
  // first and the last of them pick out one clade and no other. That is the
  // whole of the hit testing. Nothing is stored in the SVG for this page's
  // benefit, and the figure a reader saves is the figure they were looking at.
  var SPAN = /^(.*) \((.+)\), (.+) to (.+)$/;

  function openAt(title) {
    var found = SPAN.exec(title);
    if (!found) return false;
    var label = found[1] === "clade" ? found[3] + " to " + found[4] : found[1];
    trail.push({ focus: found[3] + "," + found[4], label: label + ", " + found[2] });
    later();
    return true;
  }

  function titleOf(node) {
    while (node && node !== el.plot && node !== el.stage.parentNode) {
      if (node.tagName === "g") {
        var title = node.querySelector(":scope > title");
        if (title) return title.textContent;
      }
      node = node.parentNode;
    }
    return null;
  }

  // ---------------------------------------------------------------- wiring

  function find() {
    var wanted = el.search.value.trim();
    if (!wanted) return;
    // Asked for, and taken back if the program has never heard of it. The
    // program is the only thing that knows which names the file holds, and it
    // says what it does hold when it does not hold this one, so there is no
    // index here to fall out of step with the file.
    trail.push({ focus: wanted, label: wanted });
    var argv = command();
    var room = el.plot.clientWidth - 24;
    if (room < 40) room = 900;
    asked += 1;
    working(true);
    send({
      kind: "draw",
      id: asked,
      command: K.join(argv),
      room: room,
      onFailure: "pop",
    });
    el.command.textContent = K.join(argv);
    paintTrail();
  }

  // Which of the two kinds of file this is, decided by looking at it rather
  // than by the name on it: a phylogeny comes out of a hundred programs under a
  // hundred suffixes, and a sheet is a table.
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
      el.error.textContent =
        "that reads as a table rather than a tree; drop the phylogeny first and the sheet after it";
      el.error.hidden = false;
      return;
    }
    sheet = { name: name || "traits.tsv", body: text };
    paintSheet();
    sendFiles();
    later();
  }

  function paintSheet() {
    if (!el.sheet) return;
    el.sheet.hidden = !sheet;
    if (sheet) el.sheetName.textContent = sheet.name;
  }

  function load(text, name) {
    tree = { name: name || "tree.nwk", body: text };
    // Whether the program can read the file at all is answered by the first
    // drawing rather than by a check of its own: a refusal arrives in the same
    // place either way, and asking twice is asking a million tip tree to be
    // read twice.
    trail = [];
    sheet = null;
    paintSheet();
    el.drop.hidden = true;
    el.app.hidden = false;
    el.error.hidden = true;
    sendFiles();
    later();
  }

  // ------------------------------------------------------------- examples

  // Trees made here rather than fetched, so the page carries no data files and
  // works the same offline. The big one is the point of the viewer: it is the
  // size at which a phylogeny stops having a picture.
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

  function start() {
    ["plot", "stage", "trail", "search", "rows", "rowsOut", "command", "timing", "error", "drop", "app", "file", "paste", "usePaste", "fit", "zoom", "sheet", "sheetName", "dropSheet"].forEach(
      function (name) {
        el[name] = document.getElementById("tv-" + name.toLowerCase());
      }
    );
    if (!el.plot) return;

    document.querySelectorAll("[data-projection]").forEach(function (button) {
      button.addEventListener("click", function () {
        projection = button.getAttribute("data-projection");
        document.querySelectorAll("[data-projection]").forEach(function (other) {
          other.setAttribute("aria-pressed", String(other === button));
        });
        later();
      });
    });

    el.rows.addEventListener("input", function () {
      rows = Number(el.rows.value);
      el.rowsOut.textContent = rows + " rows";
      later();
    });

    el.plot.addEventListener("click", function (event) {
      var title = titleOf(event.target);
      if (title) openAt(title);
    });

    hand();

    el.search.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        find();
      }
    });

    el.file.addEventListener("change", function () {
      var file = el.file.files[0];
      if (!file) return;
      file.text().then(function (text) { take(text, file.name); });
    });
    document.querySelectorAll("[data-example]").forEach(function (button) {
      button.addEventListener("click", function () {
        var tips = Number(button.getAttribute("data-example"));
        button.disabled = true;
        button.textContent = "building " + button.textContent;
        // Yielded to first, so the button is seen to change before a million
        // tips are built and the page stops answering for a moment.
        setTimeout(function () {
          load(example(tips), tips + "-tips.nwk");
          button.disabled = false;
          button.textContent = button.getAttribute("data-label");
        }, 0);
      });
    });
    el.usePaste.addEventListener("click", function () {
      if (el.paste.value.trim()) take(el.paste.value.trim(), "pasted.nwk");
    });

    ["dragover", "drop"].forEach(function (kind) {
      document.addEventListener(kind, function (event) { event.preventDefault(); });
    });
    document.addEventListener("drop", function (event) {
      var file = event.dataTransfer && event.dataTransfer.files[0];
      if (file) file.text().then(function (text) { take(text, file.name); });
    });

    K.onScheme(later);
    // The page does not fetch the program any more: the worker does, once, and
    // this thread only ever builds command lines and moves a transform.
    start_worker();
    el.drop.hidden = false;
    window.addEventListener("resize", later);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
