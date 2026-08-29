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
  var rows = 60;
  var projection = "rectangular";
  var pending = null;

  // ------------------------------------------------------------ the command

  // The command line this view would be drawn by, which is also the one to
  // copy into a terminal and get the same figure.
  function command() {
    var argv = ["tree:1-1", "--tree", tree.name, "--projection", projection];
    argv = argv.concat(["--max-rows", String(rows)]);
    var at = trail[trail.length - 1];
    if (at) argv = argv.concat(["--focus", at.focus]);
    argv.push("--no-region-label");
    return argv;
  }

  function draw() {
    if (!K.ready() || !tree.body) return;
    var argv = command();
    var room = el.plot.clientWidth - 24;
    var answer = K.run(K.join(argv), [{ name: tree.name, body: tree.body }], room);
    if (answer.ok) {
      el.plot.innerHTML = answer.body;
      el.error.hidden = true;
    } else {
      el.error.textContent = answer.body;
      el.error.hidden = false;
    }
    el.command.textContent = K.join(argv);
    el.timing.textContent = Math.round(answer.ms) + " ms";
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
    while (node && node !== el.plot) {
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
    // The program answers whether the name is there, and says what is when it
    // is not, so there is no index here to fall out of step with the file.
    var argv = ["tree:1-1", "--tree", tree.name, "--focus", wanted, "--max-rows", "2"];
    var answer = K.run(K.join(argv), [{ name: tree.name, body: tree.body }], 320);
    if (!answer.ok) {
      el.error.textContent = answer.body;
      el.error.hidden = false;
      return;
    }
    trail.push({ focus: wanted, label: wanted });
    later();
  }

  function load(text, name) {
    tree = { name: name || "tree.nwk", body: text };
    // Drawn once at the smallest cap there is, purely to find out whether the
    // program can read the file at all: a refusal here is a refusal a reader
    // can act on, and one after the viewer has opened looks like the viewer.
    //
    // No tip count is taken from it. Counting them in this page means parsing
    // Newick twice and disagreeing with the program over a quoted label, and
    // the counts that matter are the ones each clade already carries in its
    // own tooltip.
    var answer = K.run(
      K.join(["tree:1-1", "--tree", tree.name, "--max-rows", "2", "--no-region-label"]),
      [{ name: tree.name, body: text }],
      320
    );
    if (!answer.ok) {
      el.error.textContent = answer.body;
      el.error.hidden = false;
      return;
    }
    trail = [];
    el.drop.hidden = true;
    el.app.hidden = false;
    el.error.hidden = true;
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
    ["plot", "trail", "search", "rows", "rowsOut", "command", "timing", "error", "drop", "app", "file", "paste", "usePaste"].forEach(
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

    el.search.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        find();
      }
    });

    el.file.addEventListener("change", function () {
      var file = el.file.files[0];
      if (!file) return;
      file.text().then(function (text) { load(text, file.name); });
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
      if (el.paste.value.trim()) load(el.paste.value.trim(), "pasted.nwk");
    });

    ["dragover", "drop"].forEach(function (kind) {
      document.addEventListener(kind, function (event) { event.preventDefault(); });
    });
    document.addEventListener("drop", function (event) {
      var file = event.dataTransfer && event.dataTransfer.files[0];
      if (file) file.text().then(function (text) { load(text, file.name); });
    });

    K.onScheme(later);
    K.load().then(function () {
      el.drop.hidden = false;
      window.addEventListener("resize", later);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
