/* Every committed figure on a page, drawn by the program itself.

   The files under assets/figures are what the examples write, and the code
   that wrote each one is compiled into the page's copy of the program. So a
   figure is drawn here rather than shown: in the page's own light or dark,
   laid out at the width of the column it sits in, and redrawn under the
   reader's hand. A figure over a stretch of genome moves along it and zooms
   into it, down to the bases, because every frame is the program run again
   over a new window; a tree, a map, a circle or a sheet of panels zooms as a
   picture. Anything the figure names on a mark, it says under the pointer.

   A figure is drawn on exactly what is behind it, the page or the card it
   sits in, in that page's light or dark, so it is part of the page rather
   than a picture laid on it.

   A figure a page prints the command for, drawn in advance under it in the
   page's light and in its dark, is drawn the same way from that command, over
   the example files the site publishes under data/: the program's own parser
   says which files the command reads, and they are fetched when the figure
   comes near.

   The file stays in the page as it was written, and is what shows with
   JavaScript off or if the program never arrives. Nothing is fetched but the
   program and those example files, and nothing is sent anywhere.

   No dependencies, in a crate that has none. */

(function () {
  "use strict";

  var K = window.karyon;
  if (!K || !K.figure || !("IntersectionObserver" in window)) return;

  var FIGURE = /assets\/figures\/([A-Za-z0-9_-]+)\.svg(?:[?#].*)?$/;

  /* How far in a figure over the genome can go, in bases: below this the
     ruler has nothing left to label. And how far in a picture can go. */
  var MIN_SPAN = K.MIN_SPAN || 60;
  var MAX_PICTURE_ZOOM = 8;

  /* A figure starts drawing this far before it scrolls into view, so it is
     there by the time the reader is. */
  var AHEAD = "600px 0px";

  var serial = 0;
  var all = [];
  var tip = null;

  function scheme() {
    return K.dark() ? "dark" : "light";
  }

  /* ------------------------------------------------------------- ground */

  function channels(css) {
    var m = /rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)(?:\s*[,/]\s*([\d.]+%?))?\s*\)/.exec(css || "");
    if (!m) return null;
    var alpha = m[4] === undefined ? 1 : m[4].slice(-1) === "%" ? parseFloat(m[4]) / 100 : parseFloat(m[4]);
    return [parseFloat(m[1]), parseFloat(m[2]), parseFloat(m[3]), alpha];
  }

  function hex(c) {
    return "#" + c.slice(0, 3).map(function (v) {
      var h = Math.round(clamp(v, 0, 255)).toString(16);
      return h.length < 2 ? "0" + h : h;
    }).join("");
  }

  /* The colour a figure is drawn on: whatever is behind it, the page or the
     card it sits in. The first element up the tree that paints a background
     settles it, with any translucent ones between laid over it, so a figure
     inside a glassy panel is drawn on the colour the eye sees there. Asked
     each time the figure is drawn, since the page changes it between light
     and dark. */
  function groundOf(host) {
    var layers = [];
    for (var node = host.parentElement; node && node.nodeType === 1; node = node.parentElement) {
      var c = channels(getComputedStyle(node).backgroundColor);
      if (!c || c[3] === 0) continue;
      layers.push(c);
      if (c[3] >= 0.999) break;
    }
    var base = layers.length && layers[layers.length - 1][3] >= 0.999 ? layers.pop() : [255, 255, 255, 1];
    for (var i = layers.length - 1; i >= 0; i--) {
      var a = layers[i][3];
      base = [
        layers[i][0] * a + base[0] * (1 - a),
        layers[i][1] * a + base[1] * (1 - a),
        layers[i][2] * a + base[2] * (1 - a),
        1,
      ];
    }
    return hex(base);
  }

  function regionText(view) {
    return view.seq + ":" + view.start + "-" + view.end;
  }

  function clamp(value, low, high) {
    return Math.min(high, Math.max(low, value));
  }

  /* ------------------------------------------------------------- tooltip */

  /* One for the page, following the pointer, holding whatever the mark under
     it says about itself. The figure writes that as a `<title>`, which a
     browser shows as a native tooltip a second late and in a system font; it
     is moved into an attribute when the figure is drawn, so there is one
     tooltip and it is this one. */
  function tooltip() {
    if (tip) return tip;
    tip = document.createElement("div");
    tip.className = "k-fig-tip";
    tip.setAttribute("role", "status");
    tip.hidden = true;
    document.body.appendChild(tip);
    return tip;
  }

  function showTip(text, x, y) {
    var t = tooltip();
    if (t.textContent !== text) t.textContent = text;
    t.hidden = false;
    var pad = 14;
    var w = t.offsetWidth;
    var h = t.offsetHeight;
    var left = x + pad + w > window.innerWidth - 8 ? x - pad - w : x + pad;
    var top = y + pad + h > window.innerHeight - 8 ? y - pad - h : y + pad;
    t.style.transform = "translate(" + Math.round(left) + "px," + Math.round(top) + "px)";
  }

  function hideTip() {
    if (tip) tip.hidden = true;
  }

  /* -------------------------------------------------------------- sources */

  /* What a figure is drawn from. A committed figure is compiled into the
     program and named by its file; a command is the words a page prints,
     over the example files they name. Both answer the same three questions:
     is it ready, draw it this way, and does it run along a genome. */

  function Committed(stem) {
    this.stem = stem;
  }
  Committed.prototype.prepare = function () {
    return Promise.resolve();
  };
  Committed.prototype.draw = function (opts) {
    return K.figure(this.stem, opts);
  };
  Committed.prototype.region = function () {
    return K.figureRegion(this.stem);
  };

  var fetches = {};

  /* An example file's bytes, fetched once however many figures read it. */
  function fetched(name) {
    if (!fetches[name]) {
      fetches[name] = fetch(K.data(name)).then(function (response) {
        if (!response.ok) throw new Error(response.status + " fetching " + name);
        return response.arrayBuffer();
      }).then(function (buffer) {
        return { name: name, body: new Uint8Array(buffer) };
      });
    }
    return fetches[name];
  }

  function Command(argv) {
    this.argv = argv;
    this.files = null;
    this.ready = null;
  }
  Command.prototype.prepare = function () {
    if (this.ready) return this.ready;
    var self = this;
    var named = K.commandFiles(this.argv);
    this.ready = named.ok
      ? Promise.all(named.files.map(fetched)).then(function (files) { self.files = files; })
      : Promise.reject(new Error(named.body));
    return this.ready;
  };
  Command.prototype.draw = function (opts) {
    return K.command(this.argv, this.files || [], opts);
  };
  Command.prototype.region = function () {
    return K.commandRegion(this.argv, this.files || []);
  };

  /* ------------------------------------------------------------- a figure */

  function Figure(img, stem, options) {
    this.img = img;
    this.stem = stem;
    this.source = options.source || new Committed(stem);
    this.alt = img.getAttribute("alt") || "";
    this.thumb = !!options.thumb;
    this.wide = !!options.wide;
    this.prefix = "k" + ++serial + "-";
    this.home = null;
    this.view = null;
    this.moves = false;
    this.box = null;
    this.full = null;
    this.drawn = "";
    this.width = 0;
    this.svg = null;
    this.pending = false;
    this.pointers = {};
    this.drag = null;

    var host = document.createElement("div");
    host.className = "k-fig" + (this.thumb ? " k-fig--thumb" : "");
    if (img.classList.contains("k-wide")) host.classList.add("k-fig--wide");
    host.dataset.state = "waiting";
    var stage = document.createElement("div");
    stage.className = "k-fig__stage";
    img.parentNode.insertBefore(host, img);
    host.appendChild(stage);
    stage.appendChild(img);
    /* A figure drawn in advance twice, on the light page and on the dark, is
       both pictures until the program draws it once. */
    (options.pair || []).forEach(function (other) { stage.appendChild(other); });
    if (options.pair) host.classList.add("k-fig--pair");
    this.host = host;
    this.stage = stage;

    if (!this.thumb) this.controls();
  }

  Figure.prototype.controls = function () {
    var self = this;
    var host = this.host;
    host.tabIndex = 0;
    host.setAttribute("role", "group");
    host.setAttribute("aria-roledescription", "interactive figure");
    host.setAttribute("aria-label", this.alt);

    var bar = document.createElement("div");
    bar.className = "k-fig__bar";
    bar.innerHTML =
      '<span class="k-fig__where" aria-live="polite"></span>' +
      '<button type="button" class="k-fig__key" data-act="out" aria-label="Zoom out">&minus;</button>' +
      '<button type="button" class="k-fig__key" data-act="in" aria-label="Zoom in">+</button>' +
      '<button type="button" class="k-fig__key" data-act="reset">Reset</button>' +
      (this.wide
        ? '<button type="button" class="k-fig__key" data-act="close" aria-label="Close">&times;</button>'
        : '<button type="button" class="k-fig__key" data-act="open" aria-label="Open larger">Larger</button>') +
      '<button type="button" class="k-fig__key" data-act="save" aria-label="Save this view as SVG">SVG</button>';
    host.appendChild(bar);
    this.bar = bar;
    this.where = bar.querySelector(".k-fig__where");

    bar.addEventListener("click", function (event) {
      var key = event.target.closest("[data-act]");
      if (!key) return;
      var act = key.getAttribute("data-act");
      if (act === "in") self.zoomAt(0.5, null);
      else if (act === "out") self.zoomAt(2, null);
      else if (act === "reset") self.reset();
      else if (act === "open") openLarger(self);
      else if (act === "close") closeLarger();
      else if (act === "save") self.save();
    });

    host.addEventListener("keydown", function (event) {
      if (event.target !== host) return;
      var handled = true;
      if (event.key === "+" || event.key === "=") self.zoomAt(0.5, null);
      else if (event.key === "-" || event.key === "_") self.zoomAt(2, null);
      else if (event.key === "0") self.reset();
      else if (event.key === "ArrowLeft") self.nudge(-0.1);
      else if (event.key === "ArrowRight") self.nudge(0.1);
      else handled = false;
      if (handled) event.preventDefault();
    });

    var stage = this.stage;
    stage.addEventListener("pointerdown", function (event) { self.down(event); });
    stage.addEventListener("pointermove", function (event) { self.move(event); });
    stage.addEventListener("pointerup", function (event) { self.up(event); });
    stage.addEventListener("pointercancel", function (event) { self.up(event); });
    stage.addEventListener("pointerleave", function () { hideTip(); });
    stage.addEventListener("dblclick", function (event) {
      if (!self.svg) return;
      event.preventDefault();
      self.zoomAt(event.shiftKey ? 2 : 0.5, event);
    });
    /* The wheel scrolls the page, as it does everywhere else on it. Zooming
       takes a pinch, which a trackpad sends as a wheel with ctrl held, or the
       wheel with ctrl or cmd held on purpose. */
    stage.addEventListener(
      "wheel",
      function (event) {
        if (!self.svg || !(event.ctrlKey || event.metaKey)) return;
        event.preventDefault();
        self.zoomAt(Math.exp(clamp(event.deltaY, -60, 60) / 120), event);
      },
      { passive: false }
    );
  };

  /* Draws the figure the way it is asked for now, if that is not how it was
     last drawn. */
  Figure.prototype.draw = function (force) {
    if (!K.ready()) return;
    var width = this.thumb ? 0 : Math.round(this.stage.clientWidth);
    if (!this.thumb && width < 40) return;
    var ground = groundOf(this.host);
    var want = [scheme(), ground, width, this.view ? regionText(this.view) : ""].join("|");
    if (!force && want === this.drawn) return;
    var answer = this.source.draw({
      theme: scheme(),
      background: ground,
      width: width,
      region: this.view ? regionText(this.view) : "",
      prefix: this.prefix,
    });
    if (!answer.ok) {
      this.host.dataset.state = "static";
      this.host.title = "karyon: " + answer.body;
      return;
    }
    this.drawn = want;
    this.width = width;
    this.stage.innerHTML = answer.body;
    var svg = this.stage.firstElementChild;
    this.svg = svg;
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", this.alt);
    svg.setAttribute("focusable", "false");
    this.full = svg.getAttribute("viewBox").split(/[\s,]+/).map(Number);
    if (this.thumb) {
      svg.setAttribute(
        "preserveAspectRatio",
        this.host.classList.contains("k-fig--wide") ? "xMidYMin slice" : "xMinYMin slice"
      );
    } else {
      this.titles();
      this.plotArea();
      if (!this.moves && this.box) this.applyBox();
    }
    this.host.style.setProperty("--k-fig-ground", ground);
    this.host.dataset.scheme = scheme();
    this.host.dataset.state = "live";
    this.say();
  };

  /* Whatever a mark says about itself moves from its `<title>` to an
     attribute the page's own tooltip reads. The figure as a whole keeps its
     name on the element, for anything that reads the page aloud. */
  Figure.prototype.titles = function () {
    var titles = this.svg.querySelectorAll("title");
    for (var i = 0; i < titles.length; i++) {
      var title = titles[i];
      var owner = title.parentNode;
      if (!owner || owner === this.svg) continue;
      owner.setAttribute("data-k-tip", title.textContent);
      owner.removeChild(title);
    }
  };

  /* Where the plotting area is, in the figure's own units: every band is
     clipped to it, with its value axis to the left of it where it has one,
     and the ruler has none. So the rightmost left edge is where the genome
     starts, and every band ends at the same right edge. */
  Figure.prototype.plotArea = function () {
    var rects = this.svg.querySelectorAll("clipPath > rect");
    var left = -Infinity;
    var right = -Infinity;
    for (var i = 0; i < rects.length; i++) {
      var x = parseFloat(rects[i].getAttribute("x"));
      var w = parseFloat(rects[i].getAttribute("width"));
      if (!isFinite(x) || !isFinite(w)) continue;
      left = Math.max(left, x);
      right = Math.max(right, x + w);
    }
    this.plot = isFinite(left) && right > left ? { left: left, right: right } : null;
  };

  /* How many of the figure's own units one pixel of the page is. */
  Figure.prototype.scale = function () {
    var box = this.svg.getBoundingClientRect();
    var units = this.box ? this.box[2] : this.full[2];
    return box.width > 0 ? units / box.width : 1;
  };

  /* A point of the page, in the figure's own units. */
  Figure.prototype.unitAt = function (event) {
    var box = this.svg.getBoundingClientRect();
    var view = this.box || this.full;
    return {
      x: view[0] + (event.clientX - box.left) * (view[2] / box.width),
      y: view[1] + (event.clientY - box.top) * (view[3] / box.height),
    };
  };

  /* ------------------------------------------------------------ moving */

  Figure.prototype.zoomAt = function (factor, event) {
    if (!this.svg) return;
    if (this.moves && this.view) {
      var span = this.view.end - this.view.start + 1;
      var homeSpan = this.home.end - this.home.start + 1;
      var next = clamp(Math.round(span * factor), Math.min(MIN_SPAN, homeSpan), homeSpan);
      if (next === span) return;
      var at = 0.5;
      if (event && this.plot) {
        var unit = this.unitAt(event);
        at = clamp((unit.x - this.plot.left) / (this.plot.right - this.plot.left), 0, 1);
      }
      var anchor = this.view.start + at * (span - 1);
      this.setView(Math.round(anchor - at * (next - 1)), next);
      return;
    }
    var view = this.box || this.full.slice();
    var full = this.full;
    var w = clamp(view[2] * factor, full[2] / MAX_PICTURE_ZOOM, full[2]);
    var h = w * (full[3] / full[2]);
    var point = event ? this.unitAt(event) : { x: view[0] + view[2] / 2, y: view[1] + view[3] / 2 };
    var fx = (point.x - view[0]) / view[2];
    var fy = (point.y - view[1]) / view[3];
    this.box = [point.x - fx * w, point.y - fy * h, w, h];
    this.applyBox();
  };

  Figure.prototype.setView = function (start, span) {
    var first = this.home.start;
    var last = this.home.end;
    start = clamp(start, first, last - span + 1);
    this.view = { seq: this.home.seq, start: start, end: start + span - 1 };
    this.schedule();
  };

  Figure.prototype.nudge = function (fraction) {
    if (this.moves && this.view) {
      var span = this.view.end - this.view.start + 1;
      this.setView(this.view.start + Math.round(span * fraction), span);
    } else if (this.box) {
      this.box[0] += this.box[2] * fraction;
      this.applyBox();
    }
  };

  Figure.prototype.reset = function () {
    if (this.moves && this.home) {
      this.view = { seq: this.home.seq, start: this.home.start, end: this.home.end };
      this.schedule();
    } else if (this.box) {
      this.box = null;
      this.svg.setAttribute("viewBox", this.full.join(" "));
      this.say();
    }
  };

  /* A picture zoomed in keeps its place on the page: the element keeps the
     size it was laid out at and shows less of the drawing, and never more
     than the drawing there is. */
  Figure.prototype.applyBox = function () {
    var full = this.full;
    var b = this.box;
    if (!b) return;
    b[2] = clamp(b[2], full[2] / MAX_PICTURE_ZOOM, full[2]);
    b[3] = b[2] * (full[3] / full[2]);
    b[0] = clamp(b[0], full[0], full[0] + full[2] - b[2]);
    b[1] = clamp(b[1], full[1], full[1] + full[3] - b[3]);
    if (b[2] >= full[2] - 1e-6) {
      this.box = null;
      this.svg.setAttribute("viewBox", full.join(" "));
    } else {
      this.svg.setAttribute("viewBox", b.map(function (v) { return v.toFixed(2); }).join(" "));
    }
    this.say();
  };

  /* One redraw per frame, however many moves arrived in it. */
  Figure.prototype.schedule = function () {
    var self = this;
    if (this.pending) return;
    this.pending = true;
    requestAnimationFrame(function () {
      self.pending = false;
      self.draw();
    });
  };

  Figure.prototype.say = function () {
    if (!this.where) return;
    if (this.moves && this.view) {
      var span = this.view.end - this.view.start + 1;
      this.where.textContent =
        this.view.seq + ":" + K.grouped(this.view.start) + "-" + K.grouped(this.view.end) +
        "  ·  " + K.grouped(span) + " bases";
      this.host.classList.toggle("k-fig--moved", !(this.view.start === this.home.start && this.view.end === this.home.end));
    } else {
      var zoom = this.box ? this.full[2] / this.box[2] : 1;
      this.where.textContent = zoom > 1.001 ? "×" + zoom.toFixed(zoom < 10 ? 1 : 0) : "";
      this.host.classList.toggle("k-fig--moved", zoom > 1.001);
    }
  };

  /* ----------------------------------------------------------- pointers */

  Figure.prototype.down = function (event) {
    if (!this.svg || event.button > 0) return;
    this.pointers[event.pointerId] = { x: event.clientX, y: event.clientY };
    var ids = Object.keys(this.pointers);
    if (ids.length === 1) {
      this.drag = {
        x: event.clientX,
        y: event.clientY,
        view: this.view ? { start: this.view.start, end: this.view.end } : null,
        box: this.box ? this.box.slice() : null,
        moved: false,
      };
    } else if (ids.length === 2) {
      var a = this.pointers[ids[0]];
      var b = this.pointers[ids[1]];
      this.pinch = { distance: Math.hypot(a.x - b.x, a.y - b.y) };
      this.drag = null;
    }
    this.stage.setPointerCapture(event.pointerId);
  };

  Figure.prototype.move = function (event) {
    if (!this.svg) return;
    if (this.pointers[event.pointerId]) {
      this.pointers[event.pointerId] = { x: event.clientX, y: event.clientY };
    }
    var ids = Object.keys(this.pointers);
    if (ids.length === 2 && this.pinch) {
      var a = this.pointers[ids[0]];
      var b = this.pointers[ids[1]];
      var distance = Math.hypot(a.x - b.x, a.y - b.y);
      if (distance > 0 && this.pinch.distance > 0) {
        var middle = { clientX: (a.x + b.x) / 2, clientY: (a.y + b.y) / 2 };
        this.zoomAt(this.pinch.distance / distance, middle);
        this.pinch.distance = distance;
      }
      return;
    }
    var drag = this.drag;
    if (drag) {
      var dx = event.clientX - drag.x;
      var dy = event.clientY - drag.y;
      if (!drag.moved && Math.hypot(dx, dy) < 4) return;
      drag.moved = true;
      this.host.classList.add("k-fig--dragging");
      hideTip();
      if (this.moves && drag.view && this.plot) {
        var span = drag.view.end - drag.view.start + 1;
        var bases = (dx * this.scale() * span) / (this.plot.right - this.plot.left);
        this.setView(Math.round(drag.view.start - bases), span);
      } else if (drag.box) {
        var s = this.scale();
        this.box = [drag.box[0] - dx * s, drag.box[1] - dy * s, drag.box[2], drag.box[3]];
        this.applyBox();
      }
      return;
    }
    var owner = event.target.closest && event.target.closest("[data-k-tip]");
    if (owner && this.svg.contains(owner)) showTip(owner.getAttribute("data-k-tip"), event.clientX, event.clientY);
    else hideTip();
  };

  Figure.prototype.up = function (event) {
    delete this.pointers[event.pointerId];
    if (Object.keys(this.pointers).length < 2) this.pinch = null;
    if (!Object.keys(this.pointers).length) {
      this.drag = null;
      this.host.classList.remove("k-fig--dragging");
    }
  };

  /* ---------------------------------------------------------------- save */

  /* The view on screen as a file of its own: drawn again without a prefix,
     so it names itself the way a file the program wrote does. */
  Figure.prototype.save = function () {
    var answer = this.source.draw({
      theme: scheme(),
      width: this.moves ? this.width : 0,
      region: this.view ? regionText(this.view) : "",
      prefix: "",
    });
    if (!answer.ok) return;
    var blob = new Blob([answer.body], { type: "image/svg+xml" });
    var link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = this.stem + (this.view ? "-" + this.view.start + "-" + this.view.end : "") + ".svg";
    document.body.appendChild(link);
    link.click();
    setTimeout(function () {
      URL.revokeObjectURL(link.href);
      link.remove();
    }, 0);
  };

  /* Asks the program, once, whether this figure runs along the genome. */
  Figure.prototype.learn = function () {
    if (this.thumb || this.learnt) return;
    this.learnt = true;
    var answer = this.source.region();
    if (!answer.ok || !answer.moves) return;
    var home = K.locus(answer.region);
    if (!home) return;
    this.moves = true;
    this.home = { seq: home.seq, start: home.start, end: home.end };
    this.view = { seq: home.seq, start: home.start, end: home.end };
    this.host.classList.add("k-fig--moves");
  };

  /* ------------------------------------------------------ larger, in place */

  var larger = null;

  function openLarger(fig) {
    closeLarger();
    var dialog = document.createElement("dialog");
    dialog.className = "k-fig-dialog";
    dialog.setAttribute("aria-label", fig.alt);
    var img = fig.img.cloneNode(true);
    dialog.appendChild(img);
    document.body.appendChild(dialog);
    var big = new Figure(img, fig.stem, { wide: true, source: fig.source });
    big.host.kFigure = big;
    big.opener = fig;
    larger = { dialog: dialog, figure: big };
    dialog.addEventListener("close", function () {
      if (larger && larger.dialog === dialog) larger = null;
      dialog.remove();
      fig.host.focus();
    });
    dialog.addEventListener("click", function (event) {
      if (event.target === dialog) dialog.close();
    });
    dialog.showModal();
    big.learn();
    if (fig.view) big.view = { seq: fig.view.seq, start: fig.view.start, end: fig.view.end };
    big.draw(true);
    if (fig.box) {
      big.box = fig.box.slice();
      big.applyBox();
    }
    big.host.focus();
  }

  function closeLarger() {
    if (larger) larger.dialog.close();
  }

  /* --------------------------------------------------------------- page */

  function candidates() {
    var found = [];
    var imgs = document.querySelectorAll(".md-typeset img");
    for (var i = 0; i < imgs.length; i++) {
      var img = imgs[i];
      var match = FIGURE.exec(img.getAttribute("src") || "");
      if (!match) continue;
      if (img.closest(".k-live") || img.closest(".k-fig")) continue;
      found.push({ img: img, stem: match[1], thumb: !!img.closest("a") });
    }
    return found;
  }

  /* A figure drawn in advance under the command a page prints for it: a
     `figure.k-start` right after its code block, holding the picture for the
     light page and the one for the dark. The command is the block's words,
     less `karyon` and where it writes to, which names the file a view is
     saved as. */
  function commands() {
    var found = [];
    var figures = document.querySelectorAll(".md-typeset figure.k-start");
    for (var i = 0; i < figures.length; i++) {
      var figure = figures[i];
      if (figure.querySelector(".k-fig")) continue;
      var block = figure.previousElementSibling;
      var code = block && block.matches(".language-bash") ? block.querySelector("code") : null;
      if (!code) continue;
      var words = K.words(code.textContent);
      if (words[0] !== "karyon") continue;
      var argv = [];
      var name = "figure";
      for (var w = 1; w < words.length; w++) {
        if (words[w] === "-o" || words[w] === "--output") {
          name = (words[w + 1] || name).replace(/\.svg$/, "");
          w++;
        } else {
          argv.push(words[w]);
        }
      }
      var light = figure.querySelector("img.k-light");
      var dark = figure.querySelector("img.k-dark");
      if (!light) continue;
      found.push({
        img: light,
        stem: name,
        thumb: false,
        source: new Command(argv),
        pair: dark ? [dark] : [],
      });
    }
    return found;
  }

  function each(then) {
    for (var i = 0; i < all.length; i++) then(all[i]);
    if (larger) then(larger.figure);
  }

  function start() {
    var found = candidates().concat(commands());
    if (!found.length) return;
    found.forEach(function (item) {
      var fig = new Figure(item.img, item.stem, {
        thumb: item.thumb,
        source: item.source,
        pair: item.pair,
      });
      fig.host.kFigure = fig;
      all.push(fig);
    });

    /* Fetched now rather than when the first figure scrolls near, since a
       page that holds figures is going to want the program, and every moment
       it takes to arrive is a moment a figure shows as a blank. If it never
       arrives, every figure shows the file it was written as. */
    K.load().then(null, function () {
      all.forEach(function (fig) {
        if (fig.host.dataset.state === "waiting") fig.host.dataset.state = "static";
      });
    });

    var shown = new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          var fig = entry.target.kFigure;
          fig.visible = entry.isIntersecting;
          if (!entry.isIntersecting) return;
          K.load()
            .then(function () { return fig.source.prepare(); })
            .then(
              function () {
                fig.learn();
                fig.draw();
              },
              function (error) {
                fig.host.dataset.state = "static";
                if (error && error.message) fig.host.title = "karyon: " + error.message;
              }
            );
        });
      },
      { rootMargin: AHEAD }
    );
    all.forEach(function (fig) { shown.observe(fig.host); });

    /* The page changing between light and dark redraws what is in view now,
       and everything else when it comes into view. */
    K.onScheme(function () {
      each(function (fig) {
        if (fig.visible || fig === (larger && larger.figure)) fig.draw();
        else fig.drawn = "";
      });
    });

    /* A column that changes width redraws the figures that are laid out to
       it; a sheet, a circle or a map keeps its own width and only scales. */
    if ("ResizeObserver" in window) {
      var settle = null;
      var resized = new ResizeObserver(function () {
        clearTimeout(settle);
        settle = setTimeout(function () {
          each(function (fig) {
            if (!fig.svg || fig.thumb) return;
            var width = Math.round(fig.stage.clientWidth);
            if (Math.abs(width - fig.width) >= 4 && (fig.visible || fig.wide)) fig.draw();
          });
        }, 150);
      });
      all.forEach(function (fig) { resized.observe(fig.stage); });
    }
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
