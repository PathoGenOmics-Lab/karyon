// A figure on a page that the program draws, and that can be taken hold of.
//
// It carries no copy of anything. The command, the files, the mount and the
// controls are all read out of the page, so the text a reader sees and the
// figure they are looking at cannot drift apart.
//
// Until the program arrives the mount holds a picture drawn in advance, which
// is what the page is correct with JavaScript off, with the program missing,
// and for the moment the fetch takes.

(function () {
  "use strict";

  var K = window.karyon;
  var el = {};
  var files = [];
  var home = "";
  var live = false;
  var settle = null;

  function read(name) {
    return document.querySelector("[data-karyon-" + name + "]");
  }

  function command() {
    return el.command.textContent;
  }

  function setCommand(text) {
    el.command.textContent = text;
  }

  function say(region, rest) {
    // Debounced, so a drag writes one sentence into the live region when it
    // ends rather than sixty as it goes.
    clearTimeout(settle);
    settle = setTimeout(function () {
      el.status.textContent = region ? region + "  ·  " + rest : rest;
    }, 120);
  }

  function draw() {
    if (!K.ready()) return;
    var room = el.plot.clientWidth - 2;
    var answer = K.run(command(), files, room);

    if (answer.ok) {
      el.plot.innerHTML = answer.body;
      el.plot.classList.remove("is-refused");
      var where = K.locus(command());
      say(
        where
          ? where.seq + ":" + K.grouped(where.start) + "-" + K.grouped(where.end) +
            " (" + K.grouped(where.end - where.start + 1) + " bases)"
          : "",
        files.length + " files, drawn in " + answer.ms.toFixed(answer.ms < 10 ? 1 : 0) + " ms"
      );
    } else {
      // A refusal is an ordinary outcome here and the page says so a paragraph
      // in advance. It is the program's own sentence, printed as a sentence.
      el.plot.textContent = "karyon: " + answer.body;
      el.plot.classList.add("is-refused");
      say("", "the command did not draw. Reset puts the window back.");
    }
  }

  function move(next) {
    setCommand(next);
    draw();
  }

  // ------------------------------------------------------------ interaction

  var origin = null;

  function onDown(event) {
    if (!live || event.button !== 0) return;
    el.plot.focus();
    origin = event.clientX;
    el.plot.setPointerCapture(event.pointerId);
    el.plot.classList.add("is-dragging");
  }

  function onMove(event) {
    if (origin === null) return;
    var moved = event.clientX - origin;
    if (Math.abs(moved) < 1) return;
    origin = event.clientX;
    move(K.panned(command(), moved / Math.max(1, el.plot.clientWidth)));
  }

  function onUp(event) {
    if (origin === null) return;
    origin = null;
    el.plot.classList.remove("is-dragging");
    if (el.plot.hasPointerCapture(event.pointerId)) {
      el.plot.releasePointerCapture(event.pointerId);
    }
  }

  function onWheel(event) {
    // Only once the figure has focus, so a wheel over a figure nobody has
    // clicked scrolls the page. A figure that eats the scroll is a figure a
    // reader cannot get past.
    if (!live || document.activeElement !== el.plot) return;
    event.preventDefault();
    var box = el.plot.getBoundingClientRect();
    var at = Math.min(1, Math.max(0, (event.clientX - box.left) / box.width));
    move(K.zoomed(command(), event.deltaY > 0 ? 1.25 : 0.8, at));
  }

  function onKey(event) {
    if (!live) return;
    var step = 0.1;
    if (event.key === "ArrowLeft") move(K.panned(command(), step));
    else if (event.key === "ArrowRight") move(K.panned(command(), -step));
    else if (event.key === "+" || event.key === "=") move(K.zoomed(command(), 0.8, 0.5));
    else if (event.key === "-" || event.key === "_") move(K.zoomed(command(), 1.25, 0.5));
    else if (event.key === "Home") move(home);
    else return;
    event.preventDefault();
  }

  function start() {
    el.plot = read("plot");
    el.command = document.querySelector("[data-karyon-command] code");
    el.status = read("status");
    el.name = read("name");
    el.hint = read("hint");
    el.zoomIn = read("in");
    el.zoomOut = read("out");
    el.reset = read("reset");
    if (!el.plot || !el.command || !K) return;

    // The files are read out of the page, so what is printed and what is drawn
    // are one thing.
    files = [].slice.call(document.querySelectorAll("[data-karyon-file]")).map(function (node) {
      return { name: node.getAttribute("data-karyon-file"), body: node.textContent };
    });
    home = command();

    // The program is two hundred and seventy kilobytes over the wire, which is
    // more than everything else this page fetches put together, and the figure
    // it would redraw is already on the page as a real image. This section is
    // eighteen hundred pixels down, so a reader who never scrolls to it was
    // paying for it in full. It is fetched when the figure comes into view, or
    // the moment a reader reaches for it, whichever happens first.
    var asked = false;
    function arrive() {
      if (asked) return;
      asked = true;
      K.load()
        .then(function () {
          live = true;
          el.plot.tabIndex = 0;
          el.plot.classList.add("is-live");
          if (el.name) el.name.textContent = "karyon, running in this page";
          if (el.hint) el.hint.textContent = "drag, or use the arrow keys and + and -";
          [el.zoomIn, el.zoomOut, el.reset].forEach(function (button) {
            if (button) button.disabled = false;
          });
          draw();
        })
        .catch(function (error) {
          // The picture drawn in advance stays where it is, which is the whole
          // reason it is a real image rather than a fallback.
          say("", "the program did not arrive (" + error.message + "), so this is the picture drawn in advance");
        });
    }

    if (window.IntersectionObserver) {
      // A margin of one screen, so on an ordinary scroll the program is already
      // there by the time the figure is.
      new IntersectionObserver(function (entries, observer) {
        if (!entries.some(function (e) { return e.isIntersecting; })) return;
        observer.disconnect();
        arrive();
      }, { rootMargin: "100% 0px" }).observe(el.plot);
    } else {
      arrive();
    }

    // Reaching for it counts as asking for it, for a reader who tabs to a
    // control rather than scrolling the figure into the middle of the screen.
    ["pointerdown", "keydown", "focusin"].forEach(function (kind) {
      el.plot.addEventListener(kind, arrive, { once: false });
    });
    [el.zoomIn, el.zoomOut, el.reset].forEach(function (button) {
      if (button) button.addEventListener("focus", arrive);
    });

    el.plot.addEventListener("pointerdown", onDown);
    el.plot.addEventListener("pointermove", onMove);
    el.plot.addEventListener("pointerup", onUp);
    el.plot.addEventListener("pointercancel", onUp);
    el.plot.addEventListener("wheel", onWheel, { passive: false });
    el.plot.addEventListener("keydown", onKey);

    if (el.zoomIn) el.zoomIn.addEventListener("click", function () { move(K.zoomed(command(), 0.8, 0.5)); });
    if (el.zoomOut) el.zoomOut.addEventListener("click", function () { move(K.zoomed(command(), 1.25, 0.5)); });
    if (el.reset) el.reset.addEventListener("click", function () { move(home); });

    K.onScheme(draw);
    var wide = null;
    window.addEventListener("resize", function () {
      clearTimeout(wide);
      wide = setTimeout(draw, 150);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
